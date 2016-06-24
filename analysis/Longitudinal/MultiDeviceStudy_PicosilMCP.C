
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

// This code generates the combined Delta T histogram for the picosil center pixel and the MCP for the re-cabled configuration. The SiPad is not used.

void DoMultiDeviceStudy( string filename, string outputFilename ) {

  TFile *inputfile = TFile::Open(filename.c_str(),"READ");
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

  //Create histogram
  TH1F *histDeltaT = new TH1F("histDeltaT","; Time [ns];Number of Events", 50, 3.2, 4.2); // Weights 33-33-33
  TH1F *histDeltaTWeighted = new TH1F("histDeltaTWeighted","; Time [ns];Number of Events", 50, 2, 6); //Weights each event based on charge
  TH1F *histDeltaTCenter = new TH1F("histDeltaTCenter","; Time [ns];Number of Events", 50, 4, 5); //DeltaT of center picosil pixel
  TH1F *histDeltaTMCP = new TH1F("histDeltaTMCP","; Time [ns];Number of Events", 50, 2, 3); //DeltaT of MCP
  

  float totalCenterCharge = 0;
  float totalMCPCharge = 0;


  //read all entries and fill the histogram
  Long64_t nentries = tree->GetEntries();

  //Loop through every event in .root file
  std::cout<<"Number of events in Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss0 = gauspeak[0];
    float photekAmp0 = amp[0];
    float photekCharge0 = integral[0];

    float photekTimeGauss1 = gauspeak[9];
    float photekAmp1 = amp[9];
    float photekCharge1 = integral[9];

    float centerAmp = amp[1];
    float centerCharge = integral[1];
    float centerTime = linearTime45[1];

    float MCPAmp = amp[11];
    float MCPCharge = integral[11];
    float MCPTime = linearTime45[11];


    // APPLY EVENT CUTS:
    //require photek (for electron selection)
    if( !(photekAmp0 > 0.05 && photekCharge0 > 2) ) continue;
    
    //require signal in the central pixel
    if( !(centerAmp > 0.03 && centerCharge > 2) ) continue;

    //require Si Pad and MCP minimum amplitude. Possibly will include cut for charge at some point...
    if( !( MCPAmp > 0.03) ) continue;

    //Calculates the Delta T's if the event passes the cuts:
    float DeltaTCenter = photekTimeGauss0 - centerTime; // Will be used to fill the histogram at every event
    float DeltaTMCP = photekTimeGauss1 - MCPTime;

    float DeltaT = (DeltaTCenter + DeltaTMCP)/2; //Evenly weighted
    float DeltaTWeighted = (DeltaTCenter*centerCharge + DeltaTMCP*MCPCharge) / (centerCharge+MCPCharge);
    //The above weights each individual event being added by the charge of each event. 

    totalCenterCharge += centerCharge;
    totalMCPCharge += MCPCharge;

    histDeltaT->Fill(DeltaT);
    histDeltaTWeighted->Fill(DeltaTWeighted);
    histDeltaTCenter->Fill(DeltaTCenter);
    histDeltaTMCP->Fill(DeltaTMCP);
    
  }

  // Do Gaussian fit of delta T distributions from (mean-2RMS) to (mean+2RMS)
  double mean = histDeltaT->GetMean();
  double rms = histDeltaT->GetRMS();
  double xmin = mean-2.0*rms;
  double xmax = mean+2.0*rms;
  histDeltaT->Fit("gaus","MLES","",xmin,xmax);

  mean = histDeltaTWeighted->GetMean();
  rms = histDeltaTWeighted->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaTWeighted->Fit("gaus","MLES","",xmin,xmax);

  double meanCenter = histDeltaTCenter->GetMean();
  rms = histDeltaTCenter->GetRMS();
  xmin = meanCenter-2.0*rms;
  xmax = meanCenter+2.0*rms;
  histDeltaTCenter->Fit("gaus","MLES","",xmin,xmax);

  double meanMCP = histDeltaTMCP->GetMean();
  rms = histDeltaTMCP->GetRMS();
  xmin = meanMCP-2.0*rms;
  xmax = meanMCP+2.0*rms;
  histDeltaTMCP->Fit("gaus","MLES","",xmin,xmax);


  TH1F *histDeltaTWeightsCorrected = new TH1F("histDeltaTWeightsCorrected","; Time [ns];Number of Events", 50, -0.5, 0.5);//Corrects for histograms not initially centered at 0.
  TH1F *histDeltaTWeightsCorrected_totalcharge = new TH1F("histDeltaTWeightsCorrected_totalcharge","; Time [ns];Number of Events", 50, -0.5, 0.5);


  //Loop through again. This time, subtracting the respective means. Process is almost the same.
  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss0 = gauspeak[0];
    float photekAmp0 = amp[0];
    float photekCharge0 = integral[0];

    float photekTimeGauss1 = gauspeak[9];
    float photekAmp1 = amp[9];
    float photekCharge1 = integral[9];

    float centerAmp = amp[1];
    float centerCharge = integral[1];
    float centerTime = linearTime45[1];

    float MCPAmp = amp[11];
    float MCPCharge = integral[11];
    float MCPTime = linearTime45[11];

    if( !(photekAmp0 > 0.05 && photekCharge0 > 2) ) continue;
    if( !(centerAmp > 0.03 && centerCharge > 2) ) continue;
    if( !( MCPAmp > 0.03) ) continue;

    //Here are the corrections that center the histograms at 0.
    float DeltaTCenter = photekTimeGauss0 - centerTime - meanCenter; 
    float DeltaTMCP = photekTimeGauss1 - MCPTime - meanMCP;
    float DeltaTWeighted = (DeltaTCenter*centerCharge + DeltaTMCP*MCPCharge) / (centerCharge+MCPCharge);
    float DeltaTWeightsCorrected_totalcharge = (DeltaTCenter*totalCenterCharge + DeltaTMCP*totalMCPCharge) / (totalCenterCharge+totalMCPCharge);

    histDeltaTWeightsCorrected->Fill(DeltaTWeighted); // Uses charge in each device over each individual event as weighting.
    histDeltaTWeightsCorrected_totalcharge->Fill(DeltaTWeightsCorrected_totalcharge); // Uses total charge in each device over whole run as weighting.
  }

  mean = histDeltaTWeightsCorrected->GetMean();
  rms = histDeltaTWeightsCorrected->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaTWeightsCorrected->Fit("gaus","MLES","",xmin,xmax);

  mean = histDeltaTWeightsCorrected_totalcharge->GetMean();
  rms = histDeltaTWeightsCorrected_totalcharge->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaTWeightsCorrected_totalcharge->Fit("gaus","MLES","",xmin,xmax);


  // Creates output root file
  TFile *file = TFile::Open(outputFilename.c_str(), "RECREATE");
  file->cd();
  file->WriteTObject(histDeltaT,"histDeltaT", "WriteDelete");
  file->WriteTObject(histDeltaTWeighted,"histDeltaTWeighted", "WriteDelete");
  file->WriteTObject(histDeltaTCenter,"histDeltaTCenter", "WriteDelete");
  file->WriteTObject(histDeltaTMCP,"histDeltaTMCP", "WriteDelete");
  file->WriteTObject(histDeltaTWeightsCorrected,"histDeltaTWeightsCorrected", "WriteDelete");
  file->WriteTObject(histDeltaTWeightsCorrected_totalcharge,"histDeltaTWeightsCorrected_totalcharge", "WriteDelete");


  file->Close();
  delete file;

}

void MultiDeviceStudy_PicosilMCP() {
  DoMultiDeviceStudy("t1065-jun-2016-66.dat-full.root","output66_PicosilMCP.root");
}
