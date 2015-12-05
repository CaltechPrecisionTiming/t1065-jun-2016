#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TPostScript.h>
#include <TAxis.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TF1.h>
#include <TLine.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TMapFile.h>
#include <TPaveStats.h>


#include <fstream>
#include <iomanip>
#include <string>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>  
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

//LOCAL INCLUDES
#include "Aux.hh"

int graphic_init();
int FindMin( int n, short *a);
float GausFit_MeanTime(TGraphErrors * pulse, const float index_first, const float index_last);

TStyle* style;


int 
main(int argc, char **argv){

  double off_mean[4][9][1024];

  FILE* fp1;
  char stitle[200];
  int dummy;

  std::string boardNumber = "1";
  if (argc > 3) boardNumber = argv[3];
  std::cout << "Using Calibration files for board number " << boardNumber << "\n";

  //**************************************
  //Load Voltage Calibration
  //**************************************
  for( int i = 0; i < 4; i++){
    sprintf( stitle, "v1740_bd%s_group_%d_offset.txt", boardNumber.c_str(), i);

    fp1 = fopen( stitle, "r");
    printf("offset data : %s\n", stitle);

    for( int k = 0; k < 1024; k++)      
      for( int j = 0; j < 9; j++){      
	dummy = fscanf( fp1, "%lf ", &off_mean[i][j][k]);       
	if( k < 2 && 0)
	  printf("%5d  %8.4f\n", j, off_mean[i][j][k]);
      }
  
    fclose(fp1);
  }


  //**************************************
  //Load Time Calibration
  //**************************************
  double fdummy;
  double tcal_dV[4][1024];
  for( int i = 0; i < 4; i++){
    sprintf( stitle, "v1740_bd%s_group_%d_dV.txt", boardNumber.c_str(), i);

    fp1 = fopen( stitle, "r");
    printf("dV data : %s\n", stitle);

    for( int k = 0; k < 1024; k++)      
	dummy = fscanf( fp1, "%lf %lf %lf %lf %lf ", 
		        &fdummy, &fdummy, &fdummy, &fdummy, &tcal_dV[i][k]);       
    fclose(fp1);
  }
  double dV_sum[4] = {0, 0, 0, 0};
  for( int i = 0; i < 4; i++)
    for( int j = 0; j < 1024; j++)
    dV_sum[i] += tcal_dV[i][j];


  double tcal[4][1024];
  for( int i = 0; i < 4; i++)
    for( int j = 0; j < 1024; j++)
      tcal[i][j] = tcal_dV[i][j]/dV_sum[i]*200.0;

  


  char   title[200];  
  sprintf( title, "%s-full.root", argv[1]);

  TFile* file = new TFile( title, "RECREATE", "CAEN V1742");
  TTree* tree = new TTree("pulse", "Wave Form");

  int event;
  short   b_c[4][9][1024], tc[4]; 
  float time[4][1024];
  short channel[36][1024];
  float amp[36];
  float gauspeak[36];
  int t[36864];

  tree->Branch("event", &event, "event/I");
  tree->Branch("tc",   tc, "tc[4]/s");
  tree->Branch("b_c",  b_c, "b_c[36864]/s"); //this is for 9 channels per group
  tree->Branch("t",  t, "t[36864]/I");
  tree->Branch("channel", channel, "channel[36][1024]/S");
  tree->Branch("time", time, "time[4][1024]/F");
  tree->Branch("amp", amp, "amp[36]/F");
  tree->Branch("gauspeak", gauspeak, "gauspeak[36]/F");
  //  tree->Branch("b_c",  b_c, "b_c[32768]/s"); //this is for 8 channels

  uint   event_header;
  uint   temp[3];
  ushort samples[9][1024];

  // loop over root files
  sprintf( title, "%s", argv[1]);

  FILE* fpin = fopen( title, "r");


  for( int i  = i; i < 36864; i++ ) t[i] = i;

  for( int eventn = 0; eventn < atoi(argv[2]); eventn++){  
    // printf("---- loop  %5d\n", loop);
    event = eventn;

    //Reading First Header Word
    dummy = fread( &event_header, sizeof(uint), 1, fpin);
    uint evtSize =  event_header & 0x0fffffff;
    //Reading Second Header Word
    dummy = fread( &event_header, sizeof(uint), 1, fpin);  
    uint grM     = event_header & 0x0f;
    uint pattern = (event_header >> 8) & 0x2fff;
    uint bID     = (event_header >> 27) & 0x1f;
    //Reading Third Header Word
    dummy = fread( &event_header, sizeof(uint), 1, fpin);  
    //Reading Fourth Header Word
    dummy = fread( &event_header, sizeof(uint), 1, fpin);  

    
    //--------------------------------
    //Parsing group Mask into channels
    //--------------------------------
    bool _isGR_On[4];
    _isGR_On[0] = (grM & 0x01);
    _isGR_On[1] = (grM & 0x02);
    _isGR_On[2] = (grM & 0x04);
    _isGR_On[3] = (grM & 0x08);
    
    int ActiveGroupsN = 0;
    for ( int l = 0; l < 4; l++ )
      {
	if ( _isGR_On[l] ) ActiveGroupsN++;
      }
    
    std::cout << "------------------------------" << std::endl;
    std::cout << "evtSize: " << evtSize << " grM: " << grM << " pattern: " << pattern 
	      << " bID: " << bID  << " number of Active groups: " << ActiveGroupsN << std::endl;
    std::cout << "------------------------------" << std::endl;
    
    for( int group = 0; group < ActiveGroupsN; group++){
      dummy = fread( &event_header, sizeof(uint), 1, fpin);  

      ushort tcn = (event_header >> 20) & 0xfff;
      tc[group] = tcn;

      int nsample = (event_header & 0xfff)/3;
      std::cout << "group #"<< group << "; nsample: " << nsample << std::endl;
      //Define Time coordinate
      time[group][0] = 0.0;
      for( int i = 1; i < 1024; i++){
	time[group][i] = float(i);
	time[group][i] = float(tcal[group][(i-1+tcn)%1024] + time[group][i-1]);
	//std::cout << i << " " << time[i] << "\n";
      }      

      for(int i = 0; i < nsample; i++){
	dummy = fread( &temp, sizeof(uint), 3, fpin);  
	samples[0][i] =  temp[0] & 0xfff;
	samples[1][i] = (temp[0] >> 12) & 0xfff;
	samples[2][i] = (temp[0] >> 24) | ((temp[1] & 0xf) << 8);
	samples[3][i] = (temp[1] >>  4) & 0xfff;
	samples[4][i] = (temp[1] >> 16) & 0xfff;
	samples[5][i] = (temp[1] >> 28) | ((temp[2] & 0xff) << 4);
	samples[6][i] = (temp[2] >>  8) & 0xfff;
	samples[7][i] =  temp[2] >> 20;	
      }

      for(int j = 0; j < nsample/8; j++){
	fread( &temp, sizeof(uint), 3, fpin);  
	samples[8][j*8+0] =  temp[0] & 0xfff;
	samples[8][j*8+1] = (temp[0] >> 12) & 0xfff;
	samples[8][j*8+2] = (temp[0] >> 24) | ((temp[1] & 0xf) << 8);
	samples[8][j*8+3] = (temp[1] >>  4) & 0xfff;
	samples[8][j*8+4] = (temp[1] >> 16) & 0xfff;
	samples[8][j*8+5] = (temp[1] >> 28) | ((temp[2] & 0xff) << 4);
	samples[8][j*8+6] = (temp[2] >>  8) & 0xfff;
	samples[8][j*8+7] =  temp[2] >> 20;
      }

      for(int i = 0; i < 9; i++) {
	for(int j = 0; j < 1024; j++) {
	  b_c[group][i][j] = (double)samples[i][j];
	  channel[group*9 + i][j] = (short)((double)(samples[i][j]) - (double)(off_mean[group][i][(j+tcn)%1024]));
	}
      }
      
      double amplitude[8][1024];
      for( int i = 0; i < 8; i++)       
	for( int j = 0; j < 1024; j++){
	  amplitude[i][j] = (double)samples[i][j] - off_mean[group][i][(j+tcn)%1024];  
	  // samples[i][j] -= off_mean[group][i][(j+tcn)%1024];  
	  //amplitude[i][j] += 235;	  	  
	}

      dummy = fread( &event_header, sizeof(uint), 1, fpin);

      // time stamping
      for(int i = 0; i < 9; i++) {
	
	int index_min = FindMin (1024, channel[group*9 + i]); // return index of the min	
	TGraphErrors* pulse = GetTGraph( channel[group*9 + i], time[group] );

	Double_t min =0.; Double_t low_edge =0.; Double_t high_edge =0.; Double_t y = 0.;

	pulse->GetPoint(index_min, min, y);
	amp[group*9 + i] = y; // get amplitude
	
	pulse->GetPoint(index_min-4, low_edge, y); // get the time of the low edge of the fit range
	pulse->GetPoint(index_min+5, high_edge, y);  // get the time of the upper edge of the fit range
	
	float timepeak =  GausFit_MeanTime(pulse, low_edge, high_edge); // get the time stamp
	gauspeak[group*9 + i] = timepeak;
      }
    }
    
    tree->Fill();
  }

  fclose(fpin);

  file->Write();
  file->Close();

  return 0;
  printf("dummy = %d\n", dummy);
}



int
graphic_init( void){

  style = new TStyle("style", "style");
  
  style->SetLabelFont(132,"X");
  style->SetLabelFont(132,"Y");
  style->SetTitleFont(132,"X");
  style->SetTitleFont(132,"Y");
  style->SetTitleFont(132,"");
  style->SetTitleFontSize( 0.07);
  style->SetStatFont(132);
  style->GetAttDate()->SetTextFont(132);

  style->SetStatW(0.20);
  style->SetStatH(0.23);

  style->SetFuncColor(2);
  style->SetFuncWidth(2);
  style->SetLineWidth(2);
  
  style->SetOptFile(0);
  style->SetOptTitle(1);

  style->SetFrameBorderMode(0);
  style->SetCanvasBorderMode(0);
  style->SetPadBorderMode(0);
  style->SetTitleStyle(4000);
  style->SetPadColor(0);
  style->SetCanvasColor(0);

  style->SetTitleFillColor(0);

  style->SetTitleBorderSize(0);
  //  style->SetTitleX(0.3);
  //  style->SetTitleY(0.06);
  style->SetStatColor(0);
  style->SetStatBorderSize(1);

  style->SetOptStat("emri");
  // style->SetOptStat(1);
  style->SetOptFit(1);
  style->SetTitleOffset( 1.0,"Y");

  style->SetMarkerStyle(20);
  style->SetMarkerSize( 0.3);
  style->SetMarkerColor(4);

  // style->SetOptDate(21);

  //  style->SetPadGridX(1);
  //  style->SetPadGridY(1);


  style->cd();

  return 0;
}


////////////////////////////////////////////
// find minimum of the pulse
// aa added protection against pulses with single high bin
////////////////////////////////////////////
int FindMin( int n, short *a) {
  
  if (n <= 0 || !a) return -1;
  float xmin = a[5];
  int loc = 0;
  for  (int i = 5; i < n-5; i++) {
    if (xmin > a[i] && a[i+1] < 0.5*a[i])  {
      xmin = a[i];
      loc = i;
    }
  }
  
  return loc;
}

// find the mean time from gaus fit
float GausFit_MeanTime(TGraphErrors* pulse, const float index_first, const float index_last)
{
  TF1* fpeak = new TF1("fpeak","gaus", index_first, index_last);
  pulse->Fit("fpeak","Q","", index_first, index_last);
  
  float timepeak = fpeak->GetParameter(1);
  delete fpeak;
  
  return timepeak;
}
