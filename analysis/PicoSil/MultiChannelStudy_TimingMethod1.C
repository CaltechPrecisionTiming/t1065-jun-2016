
#include <iostream>
#include <fstream> 
#include <sstream>
#include <map>
#include <utility>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include <math.h> 




//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
TH1F* NormalizeHist(TH1F *originalHist) {
  TH1F* hist = (TH1F*)originalHist->Clone((string(originalHist->GetName())+"_normalized").c_str());
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return hist;
}



struct Pixel
{
  int   index;
  float charge;
  float time;
};

void DoMultiChannelStudy( string filename , string outputFilename) {
  // TFile *inputfile = TFile::Open("t1065-jun-2016-65To80.root","READ"); //lead 6X0, 3cm away, 32GeV
  // TFile *inputfile = TFile::Open("t1065-jun-2016-84To89.root","READ"); //lead 6X0, 0cm away, 32GeV
  // TFile *inputfile = TFile::Open("t1065-jun-2016-91.dat-full.root","READ"); //lead 6X0, 0cm away, 32GeV, after turning off pixel telescope
   // TFile *inputfile = TFile::Open("t1065-jun-2016-90.dat-full.root","READ"); //lead 6X0, 3.4cm away, 32GeV, before turning off pixel telescope
  // TFile *inputfile = TFile::Open("t1065-jun-2016-55To59.root","READ"); //tungsten 6X0, 3.4cm away, 32 GeV
  // TFile *inputfile = TFile::Open("t1065-jun-2016-34To35.root","READ"); //tungsten 6X0, 0cm away, 8 GeV
  //TFile *inputfile = TFile::Open("t1065-jun-2016-32To33.root","READ"); //tungsten 2X0, 0cm away, 32 GeV
  //TFile *inputfile = TFile::Open("t1065-jun-2016-81.dat-full.root","READ"); //lead 6X0, 0cm away, 32GeV, before turning off pixel telescope
  //TFile *inputfile = TFile::Open("t1065-jun-2016-94.dat-full.root","READ"); //lead 6X0, 1cm away, 32GeV, before turning off pixel telescope
  // TFile *inputfile = TFile::Open("t1065-jun-2016-59.dat-full.root","READ"); //tungsten 6X0, 3.4cm away, 32 GeV
  TFile *inputfile = TFile::Open(filename.c_str(),"READ"); //tungsten 6X0, 3.4cm away, 32 GeV

  TTree *tree = (TTree*)inputfile->Get("pulse");

  // get the variables from the ntuple
  float amp[36];
  float integral[36];
  float gauspeak[36];
  float linearTime45[36];

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("gauspeak",1);
  tree->SetBranchStatus("linearTime45",1);
  tree->SetBranchStatus("amp",1);
  tree->SetBranchStatus("int",1);
  
  tree->SetBranchAddress("gauspeak",gauspeak);
  tree->SetBranchAddress("linearTime45",linearTime45);
  tree->SetBranchAddress("amp",amp);
  tree->SetBranchAddress("int",integral);

  TH1F *histTotalCharge;
  TH1F *histTOFCenter;
  TH1F *histTOFFlatAvg;
  TH1F *histTOFChargeWeightedAvg;
  TH1F *histTOFRingOneFlatAvg;
  TH1F *histTOFRingOneChargeWeightedAvg;
  TH1F *histNChannelRingOne;
  TH1F *histChargeRingOne;
  TH1F *histChargeCenterOverTotalCharge;
  TH1F *histChargeRingOneOverTotalCharge;
  TH1F *histDeltaTCombined;

  TH1F *histDeltaT[7];

  for(int j=0; j<7; j++) histDeltaT[j]= new TH1F(Form("histDeltaT_%d",j),"; Time [ns];Number of Events", 100, 4,5.5);
  
  histTotalCharge = new TH1F("histTotalCharge","; Time [ns];Number of Events", 200, 0,100);
  histTOFCenter = new TH1F("histTOFCenter","; Time [ns];Number of Events", 200, -6,-4);
  histTOFFlatAvg = new TH1F("histTOFFlatAvg","; Time [ns];Number of Events", 200, -6,-4);
  histTOFChargeWeightedAvg = new TH1F("histTOFChargeWeightedAvg","; Time [ns];Number of Events", 200, -6,-4);
  histTOFRingOneFlatAvg = new TH1F("histTOFRingOneFlatAvg","; Time [ns];Number of Events", 200, -10,10);
  histTOFRingOneChargeWeightedAvg = new TH1F("histTOFRingOneChargeWeightedAvg","; Time [ns];Number of Events", 200, -10,10);
  histNChannelRingOne = new TH1F("histNChannelRingOne","; Time [ns];Number of Events", 7, -0.5,6.5);
  histChargeRingOne = new TH1F("histChargeRingOne","; Time [ns];Number of Events", 200, 0,100);
  histChargeCenterOverTotalCharge = new TH1F("histChargeCenterOverTotalCharge","; Time [ns];Number of Events", 50, 0,1);
  histChargeRingOneOverTotalCharge = new TH1F("histChargeRingOneOverTotalCharge","; Time [ns];Number of Events", 50, 0,1);
  histDeltaTCombined= new TH1F("histDeltaTCombined","; Time [ns];Number of Events", 100, 4,5.5);

  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();

  std::cout<<"Number of events in Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss = gauspeak[0];
    float photekAmp = amp[0];
    float photekCharge = integral[0];

    float centerAmp = amp[1];
    float centerCharge = integral[1];
    float centerTime = linearTime45[1];

    float RingOneCharge = 0;
    int NumberOfChannelsInRingOne = 0;
    float RingOneTimeFlatAvg = 0.0;
    float RingOneTimeChargeWeightedAvg = 0.0;
    float TotalCharge = 0;
    int NumberOfChannels = 0;
    float TimeFlatAvg = 0;
    float TimeChargeWeightedAvg = 0;

    float DeltaT[7] = {-99.};

    
    
     //require photek (for electron selection)
    if( !(photekAmp > 0.05 && photekCharge > 2) ) continue;
    
    //require signal in the central pixel
    if( !(centerAmp > 0.03 && centerCharge > 2) ) continue;
    
    for ( int j = 1; j <= 7; j++)
      {
	if ( amp[j] > 0.02 && integral[j] > 1 )
	  {
	    
	    DeltaT[j-1] = gauspeak[0] - linearTime45[j];
	    
	    NumberOfChannels++;
	    TotalCharge += integral[j];
	    TimeFlatAvg += linearTime45[j];
	    TimeChargeWeightedAvg += linearTime45[j]*integral[j];
	    if (j > 1) {
	      //cout << j << " : " << linearTime45[j] << " " << integral[j] << " : " << RingOneTimeFlatAvg << " " << RingOneTimeChargeWeightedAvg << " " << TimeFlatAvg << " " << TimeChargeWeightedAvg << "\n";
	      NumberOfChannelsInRingOne++;
	      RingOneCharge += integral[j];
	      RingOneTimeFlatAvg += linearTime45[j];
	      RingOneTimeChargeWeightedAvg += linearTime45[j]*integral[j];
	    }
	  }
      }
    
    for(int jj=0; jj<7; jj++) histDeltaT[jj]->Fill(DeltaT[jj]);
    
    // do the time averaging
    float DeltaTCombined = -99.;

    
    //if(amp[0]>0.03 && integral[0]>3. && integral[1]>3. && integral[5]>3. && integral[4]>3. && integral[3]>3)
    //DeltaTCombined = gauspeak[0]-(linearTime45[1] + linearTime45[4] + linearTime45[5]+ linearTime45[3])/4.;
    
    histDeltaTCombined->Fill(DeltaTCombined);    
    
    RingOneTimeFlatAvg = RingOneTimeFlatAvg / NumberOfChannelsInRingOne;
    RingOneTimeChargeWeightedAvg = RingOneTimeChargeWeightedAvg / RingOneCharge;
    TimeFlatAvg = TimeFlatAvg / NumberOfChannels;
    TimeChargeWeightedAvg = TimeChargeWeightedAvg / TotalCharge;

    if (NumberOfChannelsInRingOne==3) {
      histTOFRingOneFlatAvg->Fill( RingOneTimeFlatAvg - photekTimeGauss);
      histTOFRingOneChargeWeightedAvg->Fill( RingOneTimeChargeWeightedAvg - photekTimeGauss);      
    }

    //Fill histograms
    histTotalCharge->Fill(TotalCharge);
    histNChannelRingOne->Fill(NumberOfChannelsInRingOne);
    histChargeRingOne->Fill(RingOneCharge);
    histChargeCenterOverTotalCharge->Fill(centerCharge / (centerCharge + RingOneCharge));
    histChargeRingOneOverTotalCharge->Fill(RingOneCharge / (centerCharge + RingOneCharge));

    histTOFCenter->Fill(centerTime - photekTimeGauss);
    histTOFFlatAvg->Fill( TimeFlatAvg - photekTimeGauss);
    histTOFChargeWeightedAvg->Fill( TimeChargeWeightedAvg - photekTimeGauss);
    
  }
  
  float meanT[7];
  for ( int i = 0; i < 7; i++ )
    {
      meanT[i] = histDeltaT[i]->GetMean();
      std::cout << "MEAN " << i << "-->" << meanT[i] << std::endl;
    }


  TH1F *histTOFCenter_C;
  TH1F *histTOFFlatAvg_C;
  TH1F *histTOFChargeWeightedAvg_C;
  TH1F *histTOFRingOneFlatAvg_C;
  TH1F *histTOFRingOneChargeWeightedAvg_C;
  TH1F *histDeltaTCombined_C;

  TH1F *histDeltaT_C[7];

  for(int j=0; j < 7; j++) histDeltaT_C[j]= new TH1F(Form("histDeltaT_C_%d",j),"; Time [ns];Number of Events", 200, -1,1);
  
  histTOFCenter_C = new TH1F("histTOFCenter_C","; Time [ns];Number of Events", 200, -6,-4);
  TH1F* histTOF2Pixel_C = new TH1F("histTOF__Pixel_C","; Time [ns];Number of Events", 200, -1,1);
  
  histTOFFlatAvg_C = new TH1F("histTOFFlatAvg_C","; Time [ns];Number of Events", 200, -6,-4);
  histTOFChargeWeightedAvg_C = new TH1F("histTOFChargeWeightedAvg_C","; Time [ns];Number of Events", 200, -1,-1);
  histTOFRingOneFlatAvg_C = new TH1F("histTOFRingOneFlatAvg_C","; Time [ns];Number of Events", 200, -10,10);
  histTOFRingOneChargeWeightedAvg_C = new TH1F("histTOFRingOneChargeWeightedAvg_C","; Time [ns];Number of Events", 200, -10,10);
  histDeltaTCombined_C = new TH1F("histDeltaTCombined_C","; Time [ns];Number of Events", 100, 4,5.5);


  for (Long64_t iEntry=0;iEntry<nentries;iEntry++)
    {
      if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
      tree->GetEntry(iEntry);    
      float DeltaT[7] = {-99.};
      //require photek (for electron selection)
      float photekTimeGauss = gauspeak[0];
      float photekAmp = amp[0];
      float photekCharge = integral[0];

      float centerAmp = amp[1];
      float centerCharge = integral[1];
      //float Charge4 = integral[4];
      float centerTime = linearTime45[1];
    
      if( !(photekAmp > 0.03 && photekCharge > 3) ) continue;
      //require signal in the central pixel
      if( centerCharge < 3 ) continue;

      std::vector< Pixel > vect;
      Pixel pixel;
      for ( int j = 1; j <= 7; j++)
	{
	  if ( integral[j] > 3 ) histDeltaT_C[j-1]->Fill(gauspeak[0] - linearTime45[j] - meanT[j-1]);
	  pixel.index = j;
	  if ( j == 1 ) pixel.charge = 3.3*integral[j];
	  else pixel.charge = integral[j];
	  pixel.time = linearTime45[j] - meanT[j-1];
	  vect.push_back( pixel ) ;
	}

      auto sortPixel = []( Pixel a, Pixel b ) {return a.charge > b.charge ?  true : false;};
      //std::sort( vect.begin(), vect.end(), sortPixel );
      float average2;
      if ( vect[3].charge > 3 )
	{
	  //average2 = gauspeak[0] - ( vect[0].charge*vect[0].time + vect[3].charge*vect[3].time + vect[4].charge*vect[4].time + vect[2].charge*vect[2].time)/(vect[0].charge+vect[3].charge+vect[4].charge+vect[2].charge);
	  //average2 = gauspeak[0] - ( vect[0].charge*vect[0].time + vect[3].charge*vect[3].time + vect[4].charge*vect[4].time )/(vect[0].charge+vect[3].charge+vect[4].charge);
	  average2 = gauspeak[0] - (vect[0].charge*vect[0].time + vect[3].charge*vect[3].time )/(vect[0].charge+vect[3].charge);
	  histTOF2Pixel_C->Fill( average2 );
	}
      std::cout << "event: " << iEntry << "-->" << vect[0].charge << " " << vect[1].charge << " " << vect[2].charge << std::endl;
    }
  
  /*for ( int i = 0; i < 7; i++)
    {
      for ( int k = 1; k < histDeltaT[0]->GetNbinsX(); k++ )
	{
	  histDeltaT_C[i]->Fill( histDeltaT_C[i]->GetBinContent(k)-meanT[i] );
	}
    }
  */
  //Normalize hists
  histTotalCharge = NormalizeHist(histTotalCharge);
  histChargeRingOne = NormalizeHist(histChargeRingOne);
  histNChannelRingOne = NormalizeHist(histNChannelRingOne);
  histChargeCenterOverTotalCharge = NormalizeHist(histChargeCenterOverTotalCharge);
  histChargeRingOneOverTotalCharge = NormalizeHist(histChargeRingOneOverTotalCharge);

  // Do Gaussian fit of delta T distributions
  for(int j=0; j<7; j++) {
    double mean = histDeltaT[j]->GetMean();
    double rms = histDeltaT[j]->GetRMS();
    double xmin = mean-2.0*rms;
    double xmax = mean+2.0*rms;
    cout << "\nFitting Channel #" << j << ":\n" << endl;
    histDeltaT[j]->Fit("gaus","MLES","",xmin,xmax);
  }


  TFile *file = TFile::Open(outputFilename.c_str(), "RECREATE");
  file->cd();
  file->WriteTObject(histTotalCharge,"histTotalCharge", "WriteDelete");
  file->WriteTObject(histTOFCenter,"histTOFCenter", "WriteDelete");
  file->WriteTObject(histTOFFlatAvg,"histTOFFlatAvg", "WriteDelete");
  file->WriteTObject(histTOFChargeWeightedAvg,"histTOFChargeWeightedAvg", "WriteDelete");
  file->WriteTObject(histTOFRingOneFlatAvg,"histTOFRingOneFlatAvg", "WriteDelete");
  file->WriteTObject(histTOFRingOneChargeWeightedAvg,"histTOFRingOneChargeWeightedAvg", "WriteDelete");
  file->WriteTObject(histNChannelRingOne,"histNChannelRingOne", "WriteDelete");
  file->WriteTObject(histChargeRingOne,"histChargeRingOne", "WriteDelete");
  file->WriteTObject(histChargeCenterOverTotalCharge,"histChargeCenterOverTotalCharge", "WriteDelete");
  file->WriteTObject(histChargeRingOneOverTotalCharge,"histChargeRingOneOverTotalCharge", "WriteDelete");
  file->WriteTObject(histDeltaTCombined,"histDeltaTCombined", "WriteDelete");

  file->WriteTObject(histTOF2Pixel_C,"Average2Pixels", "WriteDelete");
      
  for( int i = 0; i < 7; i++ )
    {
      histDeltaT_C[i]->Write( Form("deltaT_%d_corr", i) );
      histDeltaT[i]->Write( Form("deltaT_%d", i) );
    }
  
  file->Close();
  delete file;
  
}

void MultiChannelStudy_TimingMethod1() {

  // DoMultiChannelStudy("t1065-jun-2016-90.dat-full.root","output.90.root");
  // DoMultiChannelStudy("t1065-jun-2016-94.dat-full.root","output.94.root");
  // DoMultiChannelStudy("t1065-jun-2016-81.dat-full.root","output.81.root");
  DoMultiChannelStudy("test_timing.root","output.TungstenAll.root");


}
