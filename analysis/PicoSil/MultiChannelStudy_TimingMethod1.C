
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

          if (j > 1) 
          {
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

    if (NumberOfChannelsInRingOne==3) 
    {
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
  // these are histograms to compare channels pairwise
  TH1F* histTOF2Pixel12_C = new TH1F("histTOF__Pixel12_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOF2Pixel13_C = new TH1F("histTOF__Pixel13_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOF2Pixel14_C = new TH1F("histTOF__Pixel14_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOF2Pixel15_C = new TH1F("histTOF__Pixel15_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOF2Pixel16_C = new TH1F("histTOF__Pixel16_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOF2Pixel17_C = new TH1F("histTOF__Pixel17_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOF2Pixel_All_C = new TH1F("histTOF__Pixel_All_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOF2Pixel_167_C = new TH1F("histTOF__Pixel_167_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOF2Pixel_172_C = new TH1F("histTOF__Pixel_172_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOF2Pixel_123_C = new TH1F("histTOF__Pixel_123_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOF2Pixel_134_C = new TH1F("histTOF__Pixel_134_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOF2Pixel_145_C = new TH1F("histTOF__Pixel_145_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOF2Pixel_156_C = new TH1F("histTOF__Pixel_156_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOFPixellargest_C = new TH1F("histTOF__Pixellargest_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOFPixel2largest_C = new TH1F("histTOF__Pixel2largest_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOFPixel3largest_C = new TH1F("histTOF__Pixel3largest_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOFPixel4largest_C = new TH1F("histTOF__Pixel4largest_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOFPixel5largest_C = new TH1F("histTOF__Pixel5largest_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOFPixel6largest_C = new TH1F("histTOF__Pixel6largest_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histTOFPixel7largest_C = new TH1F("histTOF__Pixel7largest_C","; Time [ns];Number of Events", 300, -0.5,0.5);
  TH1F* histMaxIndex_C = new TH1F("histMaxIndex","; Index;Number of Events", 7, 0.5,7.5);
  TH1F* histEnergy1_C = new TH1F("histEnergy1","; Percentage of total energy in event with most energy;Number of Events", 200, 0.5,100);
  TH1F* histEnergy2_C = new TH1F("histEnergy2","; Percentage of total energy in two events with most energy;Number of Events", 200, 0.5,100);
  TH1F* histEnergy3_C = new TH1F("histEnergy3","; Percentage of total energy in three events with most energy;Number of Events", 200, 0.5,100);
  TH1F* histEnergy4_C = new TH1F("histEnergy4","; Percentage of total energy in four events with most energy;Number of Events", 200, 0.5,100);
  TH1F* histEnergy5_C = new TH1F("histEnergy5","; Percentage of total energy in five events with most energy;Number of Events", 200, 0.5,100);
  TH1F* histEnergy6_C = new TH1F("histEnergy6","; Percentage of total energy in six events with most energy;Number of Events", 200, 0.5,100);
  TH1F* histEnergy7_C = new TH1F("histEnergy7","; Percentage of total energy in seven events with most energy;Number of Events", 200, 90,110);


  histTOFFlatAvg_C = new TH1F("histTOFFlatAvg_C","; Time [ns];Number of Events", 200, -6,-4);
  histTOFChargeWeightedAvg_C = new TH1F("histTOFChargeWeightedAvg_C","; Time [ns];Number of Events", 200, -1,-1);
  histTOFRingOneFlatAvg_C = new TH1F("histTOFRingOneFlatAvg_C","; Time [ns];Number of Events", 200, -10,10);
  histTOFRingOneChargeWeightedAvg_C = new TH1F("histTOFRingOneChargeWeightedAvg_C","; Time [ns];Number of Events", 200, -10,10);
  histDeltaTCombined_C = new TH1F("histDeltaTCombined_C","; Time [ns];Number of Events", 100, 4,5.5);


  for (Long64_t iEntry=0;iEntry<nentries;iEntry++)
    {
      if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
      tree->GetEntry(iEntry);    
      //float DeltaT[7] = {-99.};
      //require photek (for electron selection)
      float photekTimeGauss = gauspeak[0];
      float photekAmp = amp[0];
      float photekCharge = integral[0];

      float totalEnergy = 0;
      float Energy = 0;
      float percent = 0;

      float largestIndex = 0;
      float largestIndexIntegral = 0;

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
        pixel.time = gauspeak[0]-(linearTime45[j] + meanT[j-1]);
        vect.push_back( pixel ) ;
      }

      auto sortPixel = []( Pixel a, Pixel b ) {return a.charge > b.charge ?  true : false;};
      //std::sort( vect.begin(), vect.end(), sortPixel );
      float average ;
      if ( vect[3].charge > 3 )
      {
        //average2 = gauspeak[0] - ( vect[0].charge*vect[0].time + vect[3].charge*vect[3].time + vect[4].charge*vect[4].time + vect[2].charge*vect[2].time)/(vect[0].charge+vect[3].charge+vect[4].charge+vect[2].charge);
        //average2 = gauspeak[0] - ( vect[0].charge*vect[0].time + vect[3].charge*vect[3].time + vect[4].charge*vect[4].time )/(vect[0].charge+vect[3].charge+vect[4].charge);
        average  = (vect[0].charge*vect[0].time + vect[1].charge*vect[1].time )/(vect[0].charge+vect[1].charge);
        histTOF2Pixel12_C->Fill( average );
        average = (vect[0].charge*vect[0].time + vect[2].charge*vect[2].time )/(vect[0].charge+vect[2].charge);
        histTOF2Pixel13_C->Fill( average );
        average  = (vect[0].charge*vect[0].time + vect[3].charge*vect[3].time )/(vect[0].charge+vect[3].charge);
        histTOF2Pixel14_C->Fill( average );
        average = (vect[0].charge*vect[0].time + vect[4].charge*vect[4].time )/(vect[0].charge+vect[4].charge);
        histTOF2Pixel15_C->Fill( average );
        average  = (vect[0].charge*vect[0].time + vect[5].charge*vect[5].time )/(vect[0].charge+vect[5].charge);
        histTOF2Pixel16_C->Fill( average );
        average = (vect[0].charge*vect[0].time + vect[6].charge*vect[6].time )/(vect[0].charge+vect[6].charge);
        histTOF2Pixel17_C->Fill( average );
        average = (vect[0].charge*vect[0].time + vect[1].charge*vect[1].time + vect[2].charge*vect[2].time + vect[3].charge*vect[3].time + vect[4].charge*vect[4].time + vect[5].charge*vect[5].time + vect[6].charge*vect[6].time )/(vect[0].charge+vect[1].charge +vect[2].charge+vect[3].charge +vect[4].charge+vect[5].charge + vect[6].charge );
        histTOF2Pixel_All_C->Fill( average );
        average = (vect[0].charge*vect[0].time + vect[5].charge*vect[5].time + vect[6].charge*vect[6].time )/(vect[0].charge+vect[5].charge + vect[6].charge );
        histTOF2Pixel_167_C->Fill( average );
        average = (vect[0].charge*vect[0].time + vect[1].charge*vect[1].time + vect[6].charge*vect[6].time )/(vect[0].charge+vect[1].charge + vect[6].charge );
        histTOF2Pixel_172_C->Fill( average );
        average = (vect[0].charge*vect[0].time + vect[1].charge*vect[1].time + vect[2].charge*vect[2].time )/(vect[0].charge+vect[1].charge + vect[2].charge );
        histTOF2Pixel_123_C->Fill( average );
        average = (vect[0].charge*vect[0].time + vect[3].charge*vect[3].time + vect[2].charge*vect[2].time )/(vect[0].charge+vect[3].charge + vect[2].charge );
        histTOF2Pixel_134_C->Fill( average );
        average = (vect[0].charge*vect[0].time + vect[3].charge*vect[3].time + vect[4].charge*vect[4].time )/(vect[0].charge+vect[3].charge + vect[4].charge );
        histTOF2Pixel_145_C->Fill( average );
        average = (vect[0].charge*vect[0].time + vect[5].charge*vect[5].time + vect[4].charge*vect[4].time )/(vect[0].charge+vect[5].charge + vect[4].charge );
        histTOF2Pixel_156_C->Fill( average );
        // 1 is the center channel of the pico sil, 2-7 are the first ring, but 0 indexed here


      }

      // this will find which channel has the maximum charge and plot the results in the MaxIndex histogram. largestIndex and largestIndexIntegral starts as 0
      // channels are 1 indexed for this
      for ( int j = 0; j <= 7; j++)
        {
          if ( integral[j] > largestIndexIntegral )
            {
              largestIndexIntegral = integral[j];
              largestIndex = j;
            }
        }
        histMaxIndex_C->Fill( largestIndex );


      // Sort the vect and plot the combination of the largest energy pixels

      {
        std::sort( vect.begin(), vect.end(), sortPixel );

        // histogram of the percentage of the total energy contained in the most energetic x events
        // 2* to account for the 6db attenuator on the center pixel

        totalEnergy = 2*integral[0]+integral[1]+integral[2]+integral[3]+integral[4]+integral[5]+integral[6];
        Energy = 2*integral[0];
        percent = Energy/totalEnergy * 100;
        histEnergy1_C->Fill( percent );
        Energy = 2*integral[0]+integral[1];
        percent = Energy/totalEnergy * 100;
        histEnergy2_C->Fill( percent );
        Energy = 2*integral[0]+integral[1]+integral[2];
        percent = Energy/totalEnergy * 100;
        histEnergy3_C->Fill( percent );
        Energy = 2*integral[0]+integral[1]+integral[2]+integral[3];
        percent = Energy/totalEnergy * 100;
        histEnergy4_C->Fill( percent );
        Energy = 2*integral[0]+integral[1]+integral[2]+integral[3]+integral[4];
        percent = Energy/totalEnergy * 100;
        histEnergy5_C->Fill( percent );
        Energy = 2*integral[0]+integral[1]+integral[2]+integral[3]+integral[4]+integral[5];
        percent = Energy/totalEnergy * 100;
        histEnergy6_C->Fill( percent );
        Energy = 2*integral[0]+integral[1]+integral[2]+integral[3]+integral[4]+integral[5]+integral[6];
        percent = Energy/totalEnergy * 100;
        histEnergy7_C->Fill( percent );

        // make histogram with pixels with largest charges

        average = (vect[0].charge*vect[0].time )/(vect[0].charge );
        histTOFPixellargest_C->Fill( average );
        average = (vect[0].charge*vect[0].time + vect[1].charge*vect[1].time )/(vect[0].charge+vect[1].charge );
        histTOFPixel2largest_C->Fill( average );
        average = (vect[0].charge*vect[0].time + vect[1].charge*vect[1].time + vect[2].charge*vect[2].time )/(vect[0].charge+vect[1].charge + vect[2].charge );
        histTOFPixel3largest_C->Fill( average );
        average = (vect[0].charge*vect[0].time + vect[1].charge*vect[1].time + vect[2].charge*vect[2].time + vect[3].charge*vect[3].time )/(vect[0].charge+vect[1].charge + vect[2].charge + vect[3].charge );
        histTOFPixel4largest_C->Fill( average );
        average = (vect[0].charge*vect[0].time + vect[1].charge*vect[1].time + vect[2].charge*vect[2].time + vect[3].charge*vect[3].time + vect[4].charge*vect[4].time )/(vect[0].charge+vect[1].charge + vect[2].charge + vect[3].charge + vect[4].charge );
        histTOFPixel5largest_C->Fill( average );
        average = (vect[0].charge*vect[0].time + vect[1].charge*vect[1].time + vect[2].charge*vect[2].time + vect[3].charge*vect[3].time + vect[4].charge*vect[4].time + vect[5].charge*vect[5].time )/(vect[0].charge+vect[1].charge + vect[2].charge + vect[3].charge + vect[4].charge + vect[5].charge );
        histTOFPixel6largest_C->Fill( average );
        average = (vect[0].charge*vect[0].time + vect[1].charge*vect[1].time + vect[2].charge*vect[2].time + vect[3].charge*vect[3].time + vect[4].charge*vect[4].time + vect[5].charge*vect[5].time + vect[6].charge*vect[6].time )/(vect[0].charge+vect[1].charge +vect[2].charge+vect[3].charge +vect[4].charge+vect[5].charge + vect[6].charge );
        histTOFPixel7largest_C->Fill( average );
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

  // creates root file
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

  file->WriteTObject(histTOF2Pixel12_C,"AveragePixels12", "WriteDelete");
  file->WriteTObject(histTOF2Pixel13_C,"AveragePixels13", "WriteDelete");
  file->WriteTObject(histTOF2Pixel14_C,"AveragePixels14", "WriteDelete");
  file->WriteTObject(histTOF2Pixel15_C,"AveragePixels15", "WriteDelete");
  file->WriteTObject(histTOF2Pixel16_C,"AveragePixels16", "WriteDelete");
  file->WriteTObject(histTOF2Pixel17_C,"AveragePixels17", "WriteDelete");
  file->WriteTObject(histTOF2Pixel_All_C,"AveragePixels_All", "WriteDelete");
  file->WriteTObject(histTOF2Pixel_167_C,"AveragePixels_167", "WriteDelete");
  file->WriteTObject(histTOF2Pixel_172_C,"AveragePixels_172", "WriteDelete");
  file->WriteTObject(histTOF2Pixel_123_C,"AveragePixels_123", "WriteDelete");
  file->WriteTObject(histTOF2Pixel_134_C,"AveragePixels_134", "WriteDelete");
  file->WriteTObject(histTOF2Pixel_145_C,"AveragePixels_145", "WriteDelete");
  file->WriteTObject(histTOF2Pixel_156_C,"AveragePixels_156", "WriteDelete");
  file->WriteTObject(histTOFPixellargest_C,"AveragePixels_largest", "WriteDelete");
  file->WriteTObject(histTOFPixel2largest_C,"AveragePixels_2largest", "WriteDelete");
  file->WriteTObject(histTOFPixel3largest_C,"AveragePixels_3largest", "WriteDelete");
  file->WriteTObject(histTOFPixel4largest_C,"AveragePixels_4largest", "WriteDelete");
  file->WriteTObject(histTOFPixel5largest_C,"AveragePixels_5largest", "WriteDelete");
  file->WriteTObject(histTOFPixel6largest_C,"AveragePixels_6largest", "WriteDelete");
  file->WriteTObject(histTOFPixel7largest_C,"AveragePixels_7largest", "WriteDelete");
  file->WriteTObject(histMaxIndex_C,"MaxIndex", "WriteDelete");
  file->WriteTObject(histEnergy1_C,"PercentTotalEnergy1", "WriteDelete");
  file->WriteTObject(histEnergy2_C,"PercentTotalEnergy2", "WriteDelete");
  file->WriteTObject(histEnergy3_C,"PercentTotalEnergy3", "WriteDelete");
  file->WriteTObject(histEnergy4_C,"PercentTotalEnergy4", "WriteDelete");
  file->WriteTObject(histEnergy5_C,"PercentTotalEnergy5", "WriteDelete");
  file->WriteTObject(histEnergy6_C,"PercentTotalEnergy6", "WriteDelete");
  file->WriteTObject(histEnergy7_C,"PercentTotalEnergy7", "WriteDelete");

  for( int i = 0; i < 7; i++ )
    {
      histDeltaT_C[i]->Write( Form("deltaT_%d_corr", i+1) );
      histDeltaT[i]->Write( Form("deltaT_%d", i+1) );
    }

  file->Close();
  delete file;
  
}


void MultiChannelStudy_TimingMethod1() {

  // DoMultiChannelStudy("t1065-jun-2016-90.dat-full.root","output.90.root");
  // DoMultiChannelStudy("t1065-jun-2016-94.dat-full.root","output.94.root");
  // DoMultiChannelStudy("t1065-jun-2016-81.dat-full.root","output.81.root");
  DoMultiChannelStudy("../../raw/combine_32gev_1mm.root","output_32gev_1mm.root");
  DoMultiChannelStudy("../../raw/combine_32gev_1cm.root","output_32gev_1cm.root");
  DoMultiChannelStudy("../../raw/combine_16gev_1mm.root","output_16gev_1mm.root");


}
