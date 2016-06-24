
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

// This code generates the combined Delta T histogram for the picosil center pixel, the SiPad, and the MCP for the initial configuration. Runs 1-52 have this configuration.

void DoMultiDeviceStudy_InitialWiring( string filename, string outputFilename ) {

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
  TH1F *histDeltaT = new TH1F("histDeltaT","; Time [ns];Number of Events", 50, -5,5); // Weights 33-33-33
  TH1F *histDeltaTWeighted = new TH1F("histDeltaTWeighted","; Time [ns];Number of Events", 50, -5,5); //Weights each event based on charge
  TH1F *histDeltaTCenter = new TH1F("histDeltaTCenter","; Time [ns];Number of Events", 50, 4, 5); //DeltaT of center picosil pixel
  TH1F *histDeltaTSiPad = new TH1F("histDeltaTSiPad","; Time [ns];Number of Events", 50, -5, 5); //DeltaT of SiPad
  TH1F *histDeltaTMCP = new TH1F("histDeltaTMCP","; Time [ns];Number of Events", 50, 2, 3); //DeltaT of MCP
  

  float totalCenterCharge = 0;
  float totalSiPadCharge = 0;
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

    float photekTimeGauss3 = gauspeak[27];
    float photekAmp3 = amp[27];
    float photekCharge3 = integral[27];

    float centerAmp = amp[34];
    float centerCharge = integral[34];
    float centerTime = linearTime45[34];

    float SiPadAmp = amp[1];
    float SiPadCharge = integral[1];
    float SiPadTime = linearTime45[1];

    float MCPAmp = amp[2];
    float MCPCharge = integral[2];
    float MCPTime = linearTime45[2];


    // APPLY EVENT CUTS:
    //require photek (for electron selection)
    if( !(photekAmp0 > 0.05 && photekCharge0 > 2) ) continue;
    
    //require signal in the central pixel
    if( !(centerAmp > 0.03 && centerCharge > 2) ) continue;

    //require Si Pad and MCP minimum amplitude. Possibly will include cut for charge at some point...
    if( !(SiPadAmp > 0.03 && MCPAmp > 0.03) ) continue;

    //Calculates the Delta T's if the event passes the cuts:
    float DeltaTCenter = photekTimeGauss3 - centerTime; // Will be used to fill the histogram at every event
    float DeltaTSiPad = photekTimeGauss0 - SiPadTime;
    float DeltaTMCP = photekTimeGauss0 - MCPTime;

    float DeltaT = (DeltaTCenter + DeltaTSiPad + DeltaTMCP)/3; //Evenly weighted
    float DeltaTWeighted = (DeltaTCenter*centerCharge + DeltaTSiPad*SiPadCharge + DeltaTMCP*MCPCharge) / (centerCharge+SiPadCharge+MCPCharge);
    //The above weights each individual event being added by the charge of each event. 

    totalCenterCharge += centerCharge;
    totalSiPadCharge += SiPadCharge;
    totalMCPCharge += MCPCharge;

    histDeltaT->Fill(DeltaT);
    histDeltaTWeighted->Fill(DeltaTWeighted);
    histDeltaTCenter->Fill(DeltaTCenter);
    histDeltaTSiPad->Fill(DeltaTSiPad);
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

  double meanSiPad = histDeltaTSiPad->GetMean();
  rms = histDeltaTSiPad->GetRMS();
  xmin = meanSiPad-2.0*rms;
  xmax = meanSiPad+2.0*rms;
  histDeltaTSiPad->Fit("gaus","MLES","",xmin,xmax);

  double meanMCP = histDeltaTMCP->GetMean();
  rms = histDeltaTMCP->GetRMS();
  xmin = meanMCP-2.0*rms;
  xmax = meanMCP+2.0*rms;
  histDeltaTMCP->Fit("gaus","MLES","",xmin,xmax);


  TH1F *histDeltaTWeightsCorrected = new TH1F("histDeltaTWeightsCorrected","; Time [ns];Number of Events", 50, -2,2);//Corrects for histograms not initially centered at 0.
  TH1F *histDeltaTWeightsCorrected_totalcharge = new TH1F("histDeltaTWeightsCorrected_totalcharge","; Time [ns];Number of Events", 50, -2,2);


  //Loop through again. This time, subtracting the respective means. Process is almost the same.
  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss0 = gauspeak[0];
    float photekAmp0 = amp[0];
    float photekCharge0 = integral[0];

    float photekTimeGauss3 = gauspeak[27];
    float photekAmp3 = amp[27];
    float photekCharge3 = integral[27];

    float centerAmp = amp[34];
    float centerCharge = integral[34];
    float centerTime = linearTime45[34];

    float SiPadAmp = amp[1];
    float SiPadCharge = integral[1];
    float SiPadTime = linearTime45[1];

    float MCPAmp = amp[2];
    float MCPCharge = integral[2];
    float MCPTime = linearTime45[2];

    if( !(photekAmp0 > 0.05 && photekCharge0 > 2) ) continue;
    if( !(centerAmp > 0.03 && centerCharge > 2) ) continue;
    if( !(SiPadAmp > 0.03 && MCPAmp > 0.03) ) continue;

    //Here are the corrections that center the histograms at 0.
    float DeltaTCenter = photekTimeGauss3 - centerTime - meanCenter; 
    float DeltaTSiPad = photekTimeGauss0 - SiPadTime - meanSiPad;
    float DeltaTMCP = photekTimeGauss0 - MCPTime - meanMCP;
    float DeltaTWeighted = (DeltaTCenter*centerCharge + DeltaTSiPad*SiPadCharge + DeltaTMCP*MCPCharge) / (centerCharge+SiPadCharge+MCPCharge);
    float DeltaTWeightsCorrected_totalcharge = (DeltaTCenter*totalCenterCharge + DeltaTSiPad*totalSiPadCharge + DeltaTMCP*totalMCPCharge) / (totalCenterCharge+totalSiPadCharge+totalMCPCharge);

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
  file->WriteTObject(histDeltaTSiPad,"histDeltaTSiPad", "WriteDelete");
  file->WriteTObject(histDeltaTMCP,"histDeltaTMCP", "WriteDelete");
  file->WriteTObject(histDeltaTWeightsCorrected,"histDeltaTWeightsCorrected", "WriteDelete");
  file->WriteTObject(histDeltaTWeightsCorrected_totalcharge,"histDeltaTWeightsCorrected_totalcharge", "WriteDelete");


  file->Close();
  delete file;

}

void MultiDeviceStudy_InitialWiring() {
  DoMultiDeviceStudy_InitialWiring("t1065-jun-2016-46.dat-full.root","output46.root");
}
