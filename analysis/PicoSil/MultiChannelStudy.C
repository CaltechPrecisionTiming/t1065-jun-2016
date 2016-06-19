
#include <iostream>
#include <fstream> 
#include <sstream>
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






void MultiChannelStudy() {


  TFile *inputfile = TFile::Open("t1065-jun-2016-65To80.root","READ");
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

  TH1F *histTOFCenter;
  TH1F *histTOFFlatAvg;
  TH1F *histTOFChargeWeightedAvg;
  TH1F *histTOFRing1FlatAvg;
  TH1F *histTOFRing1ChargeWeightedAvg;
  TH1F *histNChannelRing1;
  TH1F *histChargeRing1;
  TH1F *histChargeCenterOverTotalCharge;
  TH1F *histChargeRing1OverTotalCharge;

  histTOFCenter = new TH1F("histTOFCenter","; Time [ns];Number of Events", 200, -6,-4);
  histTOFFlatAvg = new TH1F("histTOFFlatAvg","; Time [ns];Number of Events", 200, -6,-4);
  histTOFChargeWeightedAvg = new TH1F("histTOFChargeWeightedAvg","; Time [ns];Number of Events", 200, -6,-4);
  histTOFRing1FlatAvg = new TH1F("histTOFRing1FlatAvg","; Time [ns];Number of Events", 200, -6,-4);
  histTOFRing1ChargeWeightedAvg = new TH1F("histTOFRing1ChargeWeightedAvg","; Time [ns];Number of Events", 200, -6,-4);
  histNChannelRing1 = new TH1F("histNChannelRing1","; Time [ns];Number of Events", 7, -0.5,6.5);
  histChargeRing1 = new TH1F("histChargeRing1","; Time [ns];Number of Events", 200, 0,100);
  histChargeCenterOverTotalCharge = new TH1F("histChargeCenterOverTotalCharge","; Time [ns];Number of Events", 50, 0,1);
  histChargeRing1OverTotalCharge = new TH1F("histChargeRing1OverTotalCharge","; Time [ns];Number of Events", 50, 0,1);

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

    float Ring1Charge = 0;
    int NumberOfChannelsInRing1 = 0;
    float Ring1TimeFlatAvg = 0;
    float Ring1TimeChargeWeightedAvg = 0;
    float TotalCharge = 0;
    int NumberOfChannels = 0;
    float TimeFlatAvg = 0;
    float TimeChargeWeightedAvg = 0;
    for (int j=1; j <= 7; j++) {
      if (amp[j]>0.02 && integral[j]>1) {
	NumberOfChannels++;
	TotalCharge += integral[j];
	TimeFlatAvg += linearTime45[j];
	TimeChargeWeightedAvg += linearTime45[j]*integral[j];
	if (j > 1) {
	  NumberOfChannelsInRing1++;
	  Ring1Charge += integral[j];
	  Ring1TimeFlatAvg += linearTime45[j];
	  Ring1TimeChargeWeightedAvg += linearTime45[j]*integral[j];
	}      	
    }
    Ring1TimeFlatAvg = Ring1TimeFlatAvg / NumberOfChannelsInRing1;
    Ring1TimeChargeWeightedAvg = Ring1TimeChargeWeightedAvg / Ring1Charge;
    TimeFlatAvg = TimeFlatAvg / NumberOfChannels;
    TimeChargeWeightedAvg = TimeChargeWeightedAvg / TotalCharge;


    //require photek (for electron selection)
    if( !(photekAmp > 0.05 && photekCharge > 2)) continue;
    
    //require signal in the central pixel
    if( !(centerAmp > 0.03 && centerCharge > 2)) continue;

    //Fill histograms
    histNChannelRing1->Fill(NumberOfChannelsInRing1);
    histChargeRing1->Fill(Ring1Charge);
    histChargeCenterOverTotalCharge->Fill(centerCharge / (centerCharge + Ring1Charge));
    histChargeRing1OverTotalCharge->Fill(Ring1Charge / (centerCharge + Ring1Charge));

    histTOFCenter->Fill(centerTime - photekTimeGauss);
    histTOFFlatAvg->Fill( TimeFlatAvg - photekTimeGauss);
    histTOFChargeWeightedAvg->Fill( TimeChargeWeightedAvg - photekTimeGauss);
    histTOFRing1FlatAvg->Fill( Ring1TimeFlatAvg - photekTimeGauss);
    histTOFRing1ChargeWeightedAvg->Fill( Ring1TimeChargeWeightedAvg - photekTimeGauss);
    
    }
  }

  TFile *file = TFile::Open("output.root", "UPDATE");
  file->cd();
  file->WriteTObject(histTOFCenter,"histTOFCenter", "WriteDelete");
  file->WriteTObject(histTOFFlatAvg,"histTOFFlatAvg", "WriteDelete");
  file->WriteTObject(histTOFChargeWeightedAvg,"histTOFChargeWeightedAvg", "WriteDelete");
  file->WriteTObject(histTOFRing1FlatAvg,"histTOFRing1FlatAvg", "WriteDelete");
  file->WriteTObject(histTOFRing1ChargeWeightedAvg,"histTOFRing1ChargeWeightedAvg", "WriteDelete");
  file->WriteTObject(histNChannelRing1,"histNChannelRing1", "WriteDelete");
  file->WriteTObject(histChargeRing1,"histChargeRing1", "WriteDelete");
  file->WriteTObject(histChargeCenterOverTotalCharge,"histChargeCenterOverTotalCharge", "WriteDelete");
  file->WriteTObject(histChargeRing1OverTotalCharge,"histChargeRing1OverTotalCharge", "WriteDelete");
  file->Close();
  delete file;


}
