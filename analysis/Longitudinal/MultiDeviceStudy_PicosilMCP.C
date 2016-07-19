
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
#include <time.h>
#include "TRandom3.h" 

// This code generates the combined Delta T histogram for the picosil (center pixel and first ring of pixels) and the Photonis MCP for the re-cabled configuration **with respect to the Photek**. The SiPad is not used. 
// An update to the code includes use of the first ring pixels in the picosil.
// Another update gives the Delta T histogram between the picosil and the Photonis MCP.
// Author: Daniel Gawerc


// Histogram names start with histDeltaT. Then they say the devices they incorporate: either MCP, PicoSil center pixel, or PicoSil with all pixels. The PicoSil with all the Pixels will either specify equal, total charge, or event charge; this signifies all pixels weighted equally (1/7), or by charge (either weighting differently event-by-event or the same weights using total charge in the each pixel throughout the run). Then the histogram will specify equal, total charge, or event charge, with the same meanings as above, except for weighting the PicoSil delta T against the MCP delta T. Finally, a last option No Shift indicates the component histograms hadn't been shifted to zero by having their means subtracted, which should result in a histogram with a high timing resolution.
// Having At0 indicates the mean has been subtracted from the original hist.

void DoMultiDeviceStudy( string filename ) {

  srand (time(NULL));
  int seed = rand();

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

  //Create histograms
  TH1F *histDeltaT_Center_MCP_Equal = new TH1F("histDeltaT_Center_MCP_Equal","; Time [ns];Number of Events", 150, -.3, .3); // Weights MCP and PicoSil center pixel 50-50
  TH1F *histDeltaT_PicoSilEventCharge_MCP_Equal = new TH1F("histDeltaT_PicoSilEventCharge_MCP_Equal","; Time [ns];Number of Events", 150, -.3, .3);// Picosil delta T found by weighting with event pixel charge, but then overall delta T fund by weighting PicoSil and MCP equally.
  TH1F *histDeltaT_PicoSilTotalCharge_MCP_Equal = new TH1F("histDeltaT_PicoSilTotalCharge_MCP_Equal","; Time [ns];Number of Events", 150, -.3, .3);// Picosil delta T found by weighting with total pixel charge, but then overall delta T fund by weighting PicoSil and MCP equally.
  TH1F *histDeltaT_PicoSilLandauCharge_MCP_Equal = new TH1F("histDeltaT_PicoSilLandauCharge_MCP_Equal","; Time [ns];Number of Events", 150, -.3, .3);
  TH1F *histDeltaT_PicoSilLandauCharge_MCP_Equal_Smear = new TH1F("histDeltaT_PicoSilLandauCharge_MCP_Equal_Smear","; Time [ns];Number of Events", 150, -.3, .3);//Smear
  TH1F *histDeltaT_PicoSilEqual_MCP_Equal = new TH1F("histDeltaT_PicoSilEqual_MCP_Equal","; Time [ns];Number of Events", 150, -.3, .3);//Picosil delta T found by weighting pixels equally, and then overall delta T fund by weighting PicoSil and MCP equally. Thus MCP is weighted by 1/2 and each pixel is weighted 1/14.
  TH1F *histDeltaT_PicoSil_MCP_Equal = new TH1F("histDeltaT_PicoSil_MCP_Equal","; Time [ns];Number of Events", 150, -.3, .3);// delta T found by placing equal emphasis on the pixels as on the MCP. Thus everything is weighted 1/8.
  TH1F *histDeltaT_PicoSil_MCP_EventCharge = new TH1F("histDeltaT_PicoSil_MCP_EventCharge","; Time [ns];Number of Events", 150, -.3, .3);// delta T found by weighting with event charge in each device.
  TH1F *histDeltaT_PicoSil_MCP_TotalCharge = new TH1F("histDeltaT_PicoSil_MCP_TotalCharge","; Time [ns];Number of Events", 150, -.3, .3);// delta T found by constant weighting with charge in device over entire run.
  TH1F *histDeltaT_Center_MCP_EventCharge_NoShift = new TH1F("histDeltaT_Center_MCP_EventCharge_NoShift","; Time [ns];Number of Events", 100, 2, 6); //Weights each event based on charge. Poor resolution.
  TH1F *histDeltaTPicoSil[7];
  TH1F *histDeltaTPicoSilSmear[7];
  for(int i=0; i<7; i++) {
    histDeltaTPicoSil[i] = new TH1F(Form("histDeltaTPicoSil_%d",i),"; Time [ns];Number of Events", 50, 3, 6); //DeltaT of PicoSil pixels
    histDeltaTPicoSilSmear[i] = new TH1F(Form("histDeltaTPicoSilSmear_%d",i),"; Time [ns];Number of Events", 50, 3, 6);
  }
  TH1F *histDeltaTCenter = new TH1F("histDeltaTCenter","; Time [ns];Number of Events", 50, 4, 5); //DeltaT of center picosil pixel
  TH1F *histDeltaTMCP = new TH1F("histDeltaTMCP","; Time [ns];Number of Events", 50, 2, 3); //DeltaT of MCP
  TH1F *histDeltaTPicoSilAt0TotalCharge = new TH1F("histDeltaTPicoSilAt0TotalCharge","; Time [ns];Number of Events", 150, -0.3, 0.3); //All pixels combined
  TH1F *histDeltaTPicoSilAt0EventCharge = new TH1F("histDeltaTPicoSilAt0EventCharge","; Time [ns];Number of Events", 150, -0.3, 0.3);
  TH1F *histDeltaTPicoSilAt0LandauCharge = new TH1F("histDeltaTPicoSilAt0LandauCharge","; Time [ns];Number of Events", 150, -0.3, 0.3);// NEW: uses charge MPV
  TH1F *histDeltaTPicoSilAt0LandauChargeSmear = new TH1F("histDeltaTPicoSilAt0LandauChargeSmear", "; Time [ns];Number of Events", 150, -0.3, 0.3);//SKIROC
  TH1F *histDeltaTCenterAt0 = new TH1F("histDeltaTCenterAt0","; Time [ns];Number of Events", 150, -0.3, 0.3); //shifted to be centered at zero
  TH1F *histDeltaTMCPAt0 = new TH1F("histDeltaTMCPAt0","; Time [ns];Number of Events", 150, -0.3, 0.3); //shifted to be centered at zero

  TH1F *histDeltaT_PicoSil_vs_MCP_EventCharge = new TH1F("histDeltaT_PicoSil_vs_MCP_EventCharge","; Time [ns];Number of Events",150,-0.3,0.3);
  TH1F *histDeltaT_PicoSil_vs_MCP_TotalCharge = new TH1F("histDeltaT_PicoSil_vs_MCP_TotalCharge","; Time [ns];Number of Events",150,-0.3,0.3);
  TH1F *histDeltaT_PicoSil_vs_MCP[7];
  for(int i=0; i<7; i++) histDeltaT_PicoSil_vs_MCP[i] = new TH1F(Form("histDeltaT_PicoSil_vs_MCP_%d",i),"; Time [ns];Number of Events", 300, -3.5, -1.5); // DeltaT between PicoSil and MCP instead of Photek.

  TH1F *histCharges[7]; // collects charge values for picosil pixels in every event in which they pass the cuts.
  for(int i=0; i<7; i++) histCharges[i] = new TH1F( Form("histCharges_%d",i),"; Charge [pC];Number of Events", 160, 0, 80);
 
  //TH2F *histDeltaT_vs_Charge_PicoSil = new TH2F("histDeltaT_vs_Charge_PicoSil","; Event Charge [pC]; Time [ns];Number of Events", 200, 0.1, 80, 300, -0.3, 0.3);
  //TH2F *histDeltaT_vs_Charge_MCP = new TH2F("histDeltaT_vs_Charge_MCP","; Event Charge [pC]; Time [ns];Number of Events", 200, 0.1, 80, 300, -0.3, 0.3);
  // THE ABOVE TWO HISTOGRAMS ARE INITIALIZED BUT UNUSED CURRENTLY. -- Commented out.


  float totalPicoSilCharge[7] = {0.};
  float totalCenterCharge = 0; // Should be same value as totalPicoSilCharge[0]
  float totalMCPCharge = 0;

  float photekAmpCut = sqrt(10)*0.1; //THESE ARE THE CUT VALUES AFTER ADJUSTING FOR ATTENUATORS
  float photekChargeCut = sqrt(10)*2;
  float centerAmpCut = 2*0.15;
  float centerChargeCut = 2*11;
  float MCPAmpCut = 0.08;

  //read all entries and fill the histogram
  Long64_t nentries = tree->GetEntries();

  //Loop through every event in .root file
  TRandom3 *rando = new TRandom3(seed);
  std::cout<<"Number of events in Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss0 = gauspeak[0];
    float photekAmp0 = sqrt(10)*amp[0]; //accounts for 10dB attenuator
    float photekCharge0 = sqrt(10)*integral[0];

    float photekTimeGauss1 = gauspeak[9];
    float photekAmp1 = sqrt(10)*amp[9];
    float photekCharge1 = sqrt(10)*integral[9];

    float centerAmp = 2*amp[1]; //accounts for 6dB attenuator
    float centerCharge = 2*integral[1];
    float centerTime = linearTime45[1];

    float MCPAmp = amp[11];
    float MCPCharge = integral[11];
    float MCPTime = linearTime45[11];


    // APPLY EVENT CUTS:
    //require photek (for electron selection)
    if( !(photekAmp0 > photekAmpCut && photekCharge0 > photekChargeCut) ) continue;
    //require signal in the central pixel
    if( !(centerAmp > centerAmpCut && centerCharge > centerChargeCut) ) continue;
    //require MCP minimum amplitude. Possibly will include cut for charge at some point...
    if( !( MCPAmp > MCPAmpCut) ) continue;

    
    double linearTime45Smear[7];
    for (int j = 0; j < 7; j++)  linearTime45Smear[j] = rando->Gaus(linearTime45[j+1], 0.050); //Samples from smear


    //Calculates the Delta T's if the event passes the cuts:
    // (Will be used to fill the histogram at every event)
    float DeltaTCenter = photekTimeGauss0 - centerTime; 
    float DeltaTMCP = photekTimeGauss1 - MCPTime;
    float DeltaTPicoSil[7]; 
    float DeltaTPicoSilSmear[7];
    float DeltaTPicoSil_vs_MCP[7];
    std::fill(DeltaTPicoSil, DeltaTPicoSil+7, -99);
    std::fill(DeltaTPicoSilSmear, DeltaTPicoSilSmear+7, -99);
    std::fill(DeltaTPicoSil_vs_MCP, DeltaTPicoSil_vs_MCP+7, -99);
    for ( int j = 1; j <= 7; j++){
      if ( (amp[j] > 0.01 && integral[j] > 1) || j == 1 ) {
	DeltaTPicoSil[j-1] = photekTimeGauss0 - linearTime45[j];
        DeltaTPicoSilSmear[j-1] = photekTimeGauss0 - linearTime45Smear[j-1];
	DeltaTPicoSil_vs_MCP[j-1] = (linearTime45[j] - photekTimeGauss0) - (MCPTime - photekTimeGauss1); //subtracting photek accts for MCP and HGC in diff groups
	if ( j == 1 ) {
	  totalPicoSilCharge[j-1] += centerCharge;
	  histCharges[j-1]->Fill(centerCharge);
	}
	else {
	  totalPicoSilCharge[j-1] += integral[j];
	  histCharges[j-1]->Fill(integral[j]);
	}
      }
    }

    float DeltaT_Center_MCP_EventCharge_NoShift = (DeltaTCenter*centerCharge + DeltaTMCP*MCPCharge) / (centerCharge+MCPCharge);
    //The above weights each individual event being added by the charge of each event. The resolution will be large, as the mean should be subtracted from each of the events prior to weighting them by charge.

    totalCenterCharge += centerCharge;
    totalMCPCharge += MCPCharge;

    histDeltaT_Center_MCP_EventCharge_NoShift->Fill(DeltaT_Center_MCP_EventCharge_NoShift);
    histDeltaTCenter->Fill(DeltaTCenter);
    histDeltaTMCP->Fill(DeltaTMCP);
    for(int k=0; k<7; k++) {
      if (DeltaTPicoSil[k] != -99.) {
	histDeltaTPicoSil[k]->Fill(DeltaTPicoSil[k]);
        histDeltaTPicoSilSmear[k]->Fill(DeltaTPicoSilSmear[k]);
	histDeltaT_PicoSil_vs_MCP[k]->Fill(DeltaTPicoSil_vs_MCP[k]);
      }
    }
  }


  float ringChargeTotal = 0;
  for(int i=1; i<=6; i++) ringChargeTotal += totalPicoSilCharge[i];

  // Do Gaussian fit of delta T distributions from (mean-2RMS) to (mean+2RMS)
  double mean = histDeltaT_Center_MCP_EventCharge_NoShift->GetMean();
  double rms = histDeltaT_Center_MCP_EventCharge_NoShift->GetRMS();
  double xmin = mean-2.0*rms;
  double xmax = mean+2.0*rms;
  histDeltaT_Center_MCP_EventCharge_NoShift->Fit("gaus","QMLES","",xmin,xmax); // Q suppresses fit results

  double meanCenter = histDeltaTCenter->GetMean();
  rms = histDeltaTCenter->GetRMS();
  xmin = meanCenter-2.0*rms;
  xmax = meanCenter+2.0*rms;
  histDeltaTCenter->Fit("gaus","QMLES","",xmin,xmax);

  double meanMCP = histDeltaTMCP->GetMean();
  rms = histDeltaTMCP->GetRMS();
  xmin = meanMCP-2.0*rms;
  xmax = meanMCP+2.0*rms;
  histDeltaTMCP->Fit("gaus","QMLES","",xmin,xmax);


  TF1 *flandau[7];
  double MPVlandau[7]; // Seeing if we can get better results by using Landau mean as weighting charge value.
  for(int i=1; i<=6; i++) {
  mean = histCharges[i]->GetMean();
  rms = histCharges[i]->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  flandau[i] = new TF1( Form("flandau_%d",i), "landau", xmin, xmax); // 1-D landau func
  histCharges[i]->Fit( Form("flandau_%d",i), "QMLES","", xmin, xmax);
  gStyle->SetOptFit(1);
  MPVlandau[i] = flandau[i]->GetMaximumX(); // In order to find MPV.
  }
  mean = histCharges[0]->GetMean();
  rms = histCharges[0]->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  flandau[0] = new TF1("flandau_0","gaus", xmin, xmax); // storing a Gaus fit for the central pixel with the other landau fits
  histCharges[0]->Fit("flandau_0","QMLES","", xmin, xmax);
  gStyle->SetOptFit(1);
  MPVlandau[0] = flandau[0]->GetParameter(1);


  double meanPicoSil[7];
  double meanPicoSilSmear[7];
  double meanPicoSil_vs_MCP[7];
  for (int i = 0; i < 7; i++){
    meanPicoSil[i] = histDeltaTPicoSil[i]->GetMean();
    meanPicoSilSmear[i] = histDeltaTPicoSilSmear[i]->GetMean();
    meanPicoSil_vs_MCP[i] = histDeltaT_PicoSil_vs_MCP[i]->GetMean();
  }

  TRandom3 *rando2 = new TRandom3(seed);
  //Loop through again. This time, subtracting the respective means. Process is almost the same.
  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss0 = gauspeak[0];
    float photekAmp0 = sqrt(10)*amp[0]; //accounts for 10dB attenuator
    float photekCharge0 = sqrt(10)*integral[0];

    float photekTimeGauss1 = gauspeak[9];
    float photekAmp1 = sqrt(10)*amp[9];
    float photekCharge1 = sqrt(10)*integral[9];

    float centerAmp = 2*amp[1]; //accounts for 6dB attenuator
    float centerCharge = 2*integral[1];
    float centerTime = linearTime45[1];

    float MCPAmp = amp[11];
    float MCPCharge = integral[11];
    float MCPTime = linearTime45[11];

    if( !(photekAmp0 > photekAmpCut && photekCharge0 > photekChargeCut) ) continue;
    if( !(centerAmp > centerAmpCut && centerCharge > centerChargeCut) ) continue;
    if( !( MCPAmp > MCPAmpCut) ) continue;

    double linearTime45Smear[7];
    for (int j = 0; j < 7; j++)  linearTime45Smear[j] = rando2->Gaus(linearTime45[j+1],0.050);

    float DeltaTPicoSil[7] = {0.}; 
    float DeltaTPicoSilSmear[7] = {0.};
    float DeltaTPicoSil_vs_MCP[7] = {0.};
    for ( int j = 1; j <= 7; j++){
      if ( (amp[j] > 0.01 && integral[j] > 1) || j == 1 ) {
	DeltaTPicoSil[j-1] = photekTimeGauss0 - linearTime45[j] - meanPicoSil[j-1];
        DeltaTPicoSilSmear[j-1] = photekTimeGauss0 - linearTime45Smear[j-1] - meanPicoSilSmear[j-1];
	DeltaTPicoSil_vs_MCP[j-1] = (linearTime45[j] - photekTimeGauss0) - (MCPTime - photekTimeGauss1) - meanPicoSil_vs_MCP[j-1];
      }
    }


    float ringWeightTotal = 0;
    float ringWeightEvent = 0;
    float ringChargeEvent = 0;
    for(int ii = 2; ii <= 7; ii++){
      ringWeightTotal += totalPicoSilCharge[ii-1] * DeltaTPicoSil[ii-1];
      ringWeightEvent += integral[ii] * DeltaTPicoSil[ii-1];
      if(DeltaTPicoSil[ii-1] != 0.) ringChargeEvent += integral[ii];
    }


    float PicoSil_vs_MCP_WeightTotal = 0;
    float PicoSil_vs_MCP_WeightEvent = 0;
    float PicoSil_vs_MCP_ChargeTotal = 0;
    float PicoSil_vs_MCP_ChargeEvent = 0;
    for (int jj = 1; jj <= 7; jj++) {
      PicoSil_vs_MCP_WeightTotal += totalPicoSilCharge[jj-1] * DeltaTPicoSil_vs_MCP[jj-1];
      if (jj != 1) PicoSil_vs_MCP_WeightEvent += integral[jj] * DeltaTPicoSil_vs_MCP[jj-1];
      else PicoSil_vs_MCP_WeightEvent += centerCharge * DeltaTPicoSil_vs_MCP[jj-1];
      if (DeltaTPicoSil_vs_MCP[jj-1] != 0.) {
      	if ( jj != 1 ) PicoSil_vs_MCP_ChargeEvent += integral[jj];
      	else PicoSil_vs_MCP_ChargeEvent += centerCharge;
      	PicoSil_vs_MCP_ChargeTotal += totalPicoSilCharge[jj-1];
      }
    }

    
    //Here are the corrections that center the histograms at 0.
    float DeltaTCenter = photekTimeGauss0 - centerTime - meanCenter; 
    float DeltaTMCP = photekTimeGauss1 - MCPTime - meanMCP;

    float DeltaT_Center_MCP_Equal = (DeltaTCenter + DeltaTMCP)/2; //MCP and center pixel evenly weighted
    float DeltaT_PicoSil_MCP_EventCharge = (DeltaTCenter*centerCharge + DeltaTMCP*MCPCharge + ringWeightEvent) / (centerCharge+MCPCharge+ringChargeEvent);
    float DeltaT_PicoSil_MCP_TotalCharge = (DeltaTCenter*totalCenterCharge + DeltaTMCP*totalMCPCharge + ringWeightTotal) / (totalCenterCharge+totalMCPCharge+ringChargeTotal);

    float DeltaT_PicoSilTotalCharge_MCP_Equal = 0.5*DeltaTMCP; //Add PicoSil on next line.
    float temp_pstotalweight = 0;
    for (int j = 0; j <= 6; j++) temp_pstotalweight += DeltaTPicoSil[j]*totalPicoSilCharge[j];
    float temp_pstotalcharge = 0;
    for (int j = 0; j <= 6; j++) {if (DeltaTPicoSil[j] != 0) temp_pstotalcharge += totalPicoSilCharge[j];}
    DeltaT_PicoSilTotalCharge_MCP_Equal += 0.5*temp_pstotalweight/temp_pstotalcharge;

    float DeltaT_PicoSilEventCharge_MCP_Equal = 0.5*DeltaTMCP;
    float temp_pseventweight = 0;
    for (int j = 1; j <= 6; j++) temp_pseventweight += DeltaTPicoSil[j]*integral[j+1];
    temp_pseventweight += DeltaTPicoSil[0]*centerCharge;
    float temp_pseventcharge = 0;
    for (int j = 1; j <= 6; j++) {if (DeltaTPicoSil[j] != 0) temp_pseventcharge += integral[j+1];}
    if (DeltaTPicoSil[0] != 0) temp_pseventcharge += centerCharge;
    DeltaT_PicoSilEventCharge_MCP_Equal += 0.5*temp_pseventweight/temp_pseventcharge;

    float DeltaT_PicoSilLandauCharge_MCP_Equal = 0.5*DeltaTMCP;
    float temp_pslandauweight = 0;
    float temp_pslandaucharge = 0;
    for (int j = 0; j <= 6; j++) {
      temp_pslandauweight += DeltaTPicoSil[j]*MPVlandau[j];
      if (DeltaTPicoSil[j] != 0) temp_pslandaucharge += MPVlandau[j];
    }
    DeltaT_PicoSilLandauCharge_MCP_Equal += 0.5*temp_pslandauweight/temp_pslandaucharge;

    float DeltaT_PicoSilLandauCharge_MCP_Equal_Smear = 0.5*DeltaTMCP; // MCP, photek not smeared
    float temp_pslandauweight_smear = 0;
    float temp_pslandaucharge_smear = 0;
    for (int j = 0; j <= 6; j++) {
      temp_pslandauweight_smear += DeltaTPicoSilSmear[j]*MPVlandau[j];
      if (DeltaTPicoSilSmear[j] != 0) temp_pslandaucharge_smear += MPVlandau[j];
    }
    DeltaT_PicoSilLandauCharge_MCP_Equal_Smear += 0.5*temp_pslandauweight_smear/temp_pslandaucharge_smear;

    float DeltaT_PicoSil_MCP_Equal = DeltaTMCP;
    int inc = 1; //Divide by inc at the end, because some PicoSil pixels may not have passed cuts.
    for (int j = 0; j <= 6; j++) {
      DeltaT_PicoSil_MCP_Equal += DeltaTPicoSil[j];
      if (DeltaTPicoSil[j] != 0.) inc += 1;
    }
    DeltaT_PicoSil_MCP_Equal /= inc;

    float DeltaT_PicoSilEqual_MCP_Equal = 0.5*DeltaTMCP;
    inc = 1;
    float temp_psdeltat = 0;
    for (int j = 0; j <= 6; j++) {
      temp_psdeltat += DeltaTPicoSil[j];
      if (DeltaTPicoSil[j] != 0.) inc += 1;
    }
    temp_psdeltat /= inc;
    DeltaT_PicoSilEqual_MCP_Equal += 0.5*temp_psdeltat;

    float DeltaT_PicoSil_vs_MCP_TotalCharge = PicoSil_vs_MCP_WeightTotal / PicoSil_vs_MCP_ChargeTotal;
    float DeltaT_PicoSil_vs_MCP_EventCharge = PicoSil_vs_MCP_WeightEvent / PicoSil_vs_MCP_ChargeEvent;


    histDeltaT_Center_MCP_Equal->Fill(DeltaT_Center_MCP_Equal);
    histDeltaT_PicoSil_MCP_EventCharge->Fill(DeltaT_PicoSil_MCP_EventCharge); 
    histDeltaT_PicoSil_MCP_TotalCharge->Fill(DeltaT_PicoSil_MCP_TotalCharge); 
    histDeltaT_PicoSilTotalCharge_MCP_Equal->Fill(DeltaT_PicoSilTotalCharge_MCP_Equal);
    histDeltaT_PicoSilEventCharge_MCP_Equal->Fill(DeltaT_PicoSilEventCharge_MCP_Equal);
    histDeltaT_PicoSilLandauCharge_MCP_Equal->Fill(DeltaT_PicoSilLandauCharge_MCP_Equal);
    histDeltaT_PicoSilLandauCharge_MCP_Equal_Smear->Fill(DeltaT_PicoSilLandauCharge_MCP_Equal_Smear);
    histDeltaT_PicoSil_MCP_Equal->Fill(DeltaT_PicoSil_MCP_Equal);
    histDeltaT_PicoSilEqual_MCP_Equal->Fill(DeltaT_PicoSilEqual_MCP_Equal);

    histDeltaTPicoSilAt0TotalCharge->Fill(temp_pstotalweight/temp_pstotalcharge);
    histDeltaTPicoSilAt0EventCharge->Fill(temp_pseventweight/temp_pseventcharge);
    histDeltaTPicoSilAt0LandauCharge->Fill(temp_pslandauweight/temp_pslandaucharge);
    histDeltaTPicoSilAt0LandauChargeSmear->Fill(temp_pslandauweight_smear/temp_pslandaucharge_smear);
    histDeltaTCenterAt0->Fill(DeltaTCenter);
    histDeltaTMCPAt0->Fill(DeltaTMCP);

    histDeltaT_PicoSil_vs_MCP_TotalCharge->Fill(DeltaT_PicoSil_vs_MCP_TotalCharge);
    histDeltaT_PicoSil_vs_MCP_EventCharge->Fill(DeltaT_PicoSil_vs_MCP_EventCharge);
  }

  mean = histDeltaT_Center_MCP_Equal->GetMean();
  rms = histDeltaT_Center_MCP_Equal->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_Center_MCP_Equal->Fit("gaus","QMLES","",xmin,xmax);
  gStyle->SetOptFit(1);

  mean = histDeltaT_PicoSil_MCP_Equal->GetMean();
  rms = histDeltaT_PicoSil_MCP_Equal->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_PicoSil_MCP_Equal->Fit("gaus","QMLES","",xmin,xmax);
  gStyle->SetOptFit(1);

  mean = histDeltaT_PicoSilEqual_MCP_Equal->GetMean();
  rms = histDeltaT_PicoSilEqual_MCP_Equal->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_PicoSilEqual_MCP_Equal->Fit("gaus","QMLES","",xmin,xmax);
  gStyle->SetOptFit(1);

  mean = histDeltaT_PicoSil_MCP_EventCharge->GetMean();
  rms = histDeltaT_PicoSil_MCP_EventCharge->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_PicoSil_MCP_EventCharge->Fit("gaus","QMLES","",xmin,xmax);
  gStyle->SetOptFit(1);

  mean = histDeltaT_PicoSil_MCP_TotalCharge->GetMean();
  rms = histDeltaT_PicoSil_MCP_TotalCharge->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_PicoSil_MCP_TotalCharge->Fit("gaus","QMLES","",xmin,xmax);
  gStyle->SetOptFit(1);

  mean = histDeltaT_PicoSilTotalCharge_MCP_Equal->GetMean();
  rms = histDeltaT_PicoSilTotalCharge_MCP_Equal->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_PicoSilTotalCharge_MCP_Equal->Fit("gaus","QMLES","",xmin,xmax);
  gStyle->SetOptFit(1);

  mean = histDeltaT_PicoSilEventCharge_MCP_Equal->GetMean();
  rms = histDeltaT_PicoSilEventCharge_MCP_Equal->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_PicoSilEventCharge_MCP_Equal->Fit("gaus","QMLES","",xmin,xmax);
  gStyle->SetOptFit(1);

  mean = histDeltaT_PicoSilLandauCharge_MCP_Equal->GetMean();
  rms = histDeltaT_PicoSilLandauCharge_MCP_Equal->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_PicoSilLandauCharge_MCP_Equal->Fit("gaus","QMLES","",xmin,xmax);
  gStyle->SetOptFit(1);

  mean = histDeltaT_PicoSilLandauCharge_MCP_Equal_Smear->GetMean();
  rms = histDeltaT_PicoSilLandauCharge_MCP_Equal_Smear->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_PicoSilLandauCharge_MCP_Equal_Smear->Fit("gaus","QMLES","",xmin,xmax);
  gStyle->SetOptFit(1);

  mean = histDeltaTCenterAt0->GetMean();
  rms = histDeltaTCenterAt0->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaTCenterAt0->Fit("gaus","QMLES","",xmin,xmax);
  gStyle->SetOptFit(1);

  mean = histDeltaTMCPAt0->GetMean();
  rms = histDeltaTMCPAt0->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaTMCPAt0->Fit("gaus","QMLES","",xmin,xmax);
  gStyle->SetOptFit(1);

  mean = histDeltaTPicoSilAt0TotalCharge->GetMean();
  rms = histDeltaTPicoSilAt0TotalCharge->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaTPicoSilAt0TotalCharge->Fit("gaus","QMLES","",xmin,xmax);
  gStyle->SetOptFit(1);

  mean = histDeltaTPicoSilAt0EventCharge->GetMean();
  rms = histDeltaTPicoSilAt0EventCharge->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaTPicoSilAt0EventCharge->Fit("gaus","QMLES","",xmin,xmax);
  gStyle->SetOptFit(1);

  mean = histDeltaTPicoSilAt0LandauCharge->GetMean();
  rms = histDeltaTPicoSilAt0LandauCharge->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaTPicoSilAt0LandauCharge->Fit("gaus","QMLES","",xmin,xmax);
  gStyle->SetOptFit(1);

  mean = histDeltaTPicoSilAt0LandauChargeSmear->GetMean();
  rms = histDeltaTPicoSilAt0LandauChargeSmear->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaTPicoSilAt0LandauChargeSmear->Fit("gaus","QMLES","",xmin,xmax);
  gStyle->SetOptFit(1);

  mean = histDeltaT_PicoSil_vs_MCP_TotalCharge->GetMean();
  rms = histDeltaT_PicoSil_vs_MCP_TotalCharge->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_PicoSil_vs_MCP_TotalCharge->Fit("gaus","QMLES","",xmin,xmax);
  gStyle->SetOptFit(1);

  mean = histDeltaT_PicoSil_vs_MCP_EventCharge->GetMean();
  rms = histDeltaT_PicoSil_vs_MCP_EventCharge->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_PicoSil_vs_MCP_EventCharge->Fit("gaus","QMLES","",xmin,xmax);
  gStyle->SetOptFit(1);


  // Creates output root file
  TFile *file = TFile::Open(("output"+filename).c_str(), "RECREATE");
  file->cd();
  file->WriteTObject(histDeltaTCenterAt0,"histDeltaTCenter", "WriteDelete");
  file->WriteTObject(histDeltaTPicoSilAt0EventCharge,"histDeltaTPicoSilEventCharge", "WriteDelete");
  file->WriteTObject(histDeltaTPicoSilAt0TotalCharge,"histDeltaTPicoSilTotalCharge", "WriteDelete");
  file->WriteTObject(histDeltaTPicoSilAt0LandauCharge,"histDeltaTPicoSilLandauCharge", "WriteDelete");
  file->WriteTObject(histDeltaTPicoSilAt0LandauChargeSmear,"histDeltaTPicoSilLandauChargeSmear", "WriteDelete");
  file->WriteTObject(histDeltaTMCPAt0,"histDeltaTMCP", "WriteDelete");
  file->WriteTObject(histDeltaT_Center_MCP_EventCharge_NoShift,"histDeltaT_Center_MCP_EventCharge_NoShift", "WriteDelete");
  file->WriteTObject(histDeltaT_Center_MCP_Equal,"histDeltaT_Center_MCP_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSil_MCP_Equal,"histDeltaT_PicoSil_MCP_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSilEqual_MCP_Equal,"histDeltaT_PicoSilEqual_MCP_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSilEventCharge_MCP_Equal,"histDeltaT_PicoSilEventCharge_MCP_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSilTotalCharge_MCP_Equal,"histDeltaT_PicoSilTotalCharge_MCP_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSilLandauCharge_MCP_Equal,"histDeltaT_PicoSilLandauCharge_MCP_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSilLandauCharge_MCP_Equal_Smear,"histDeltaT_PicoSilLandauCharge_MCP_Equal_Smear", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSil_MCP_EventCharge,"histDeltaT_PicoSil_MCP_EventCharge", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSil_MCP_TotalCharge,"histDeltaT_PicoSil_MCP_TotalCharge", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSil_vs_MCP_TotalCharge,"histDeltaT_PicoSil_vs_MCP_TotalCharge", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSil_vs_MCP_EventCharge,"histDeltaT_PicoSil_vs_MCP_EventCharge", "WriteDelete");
  for(int i=0; i<=6; i++) {
    file->WriteTObject(histDeltaT_PicoSil_vs_MCP[i],Form("histDeltaT_PicoSil_vs_MCP[%d]",i),"WriteDelete");
  }  
  for(int i=0; i<=6; i++) {
    file->WriteTObject(histCharges[i],Form("histCharges[%d]",i),"WriteDelete");
  }


  TH1F *histPhotekAmpCut = new TH1F("histPhotekAmpCut","; Amp;Number of Events", 400, 0, 2.5);
  TH1F *histPhotekChargeCut = new TH1F("histPhotekChargeCut","; Charge;Number of Events", 400, 0, 30);
  TH1F *histCenterAmpCut = new TH1F("histCenterAmpCut","; Amp;Number of Events", 200, 0, 1.5);
  TH1F *histCenterChargeCut = new TH1F("histCenterChargeCut","; Charge;Number of Events", 400, 0, 60);
  TH1F *histMCPAmpCut = new TH1F("histMCPAmpCut","; Amp;Number of Events", 100, 0, 0.75);

  tree->Draw("sqrt(10)*amp[0]>>histPhotekAmpCut", Form("sqrt(10)*amp[0]>%f",photekAmpCut) );
  tree->Draw("sqrt(10)*int[0]>>histPhotekChargeCut", Form("sqrt(10)*int[0]>%f",photekChargeCut));
  tree->Draw("2*amp[1]>>histCenterAmpCut", Form("2*amp[1]>%f",centerAmpCut));
  tree->Draw("2*int[1]>>histCenterChargeCut", Form("2*int[1]>%f",centerChargeCut));
  tree->Draw("amp[11]>>histMCPAmpCut", Form("amp[11]>%f",MCPAmpCut));

  TH1F *histPhotekAmp = new TH1F("histPhotekAmp","; Amp;Number of Events", 400, 0, 2.5);
  TH1F *histPhotekCharge = new TH1F("histPhotekCharge","; Charge;Number of Events", 400, 0, 30);
  TH1F *histCenterAmp = new TH1F("histCenterAmp","; Amp;Number of Events", 200, 0, 1.5);
  TH1F *histCenterCharge = new TH1F("histCenterCharge","; Charge;Number of Events", 400, 0, 60);
  TH1F *histMCPAmp = new TH1F("histMCPAmp","; Amp;Number of Events", 100, 0, 0.75);

  tree->Draw("sqrt(10)*amp[0]>>histPhotekAmp" );
  tree->Draw("sqrt(10)*int[0]>>histPhotekCharge");
  tree->Draw("2*amp[1]>>histCenterAmp");
  tree->Draw("2*int[1]>>histCenterCharge");
  tree->Draw("amp[11]>>histMCPAmp");

  file->WriteTObject(histPhotekAmp, "Photek Amp", "WriteDelete");
  file->WriteTObject(histPhotekAmpCut, "Cut on Photek Amp", "WriteDelete");
  file->WriteTObject(histPhotekCharge, "Photek Charge", "WriteDelete");
  file->WriteTObject(histPhotekChargeCut, "Cut on Photek Charge", "WriteDelete");
  file->WriteTObject(histCenterAmp, "Center Amp", "WriteDelete");
  file->WriteTObject(histCenterAmpCut, "Cut on Center Amp", "WriteDelete");
  file->WriteTObject(histCenterCharge, "Center Charge", "WriteDelete");
  file->WriteTObject(histCenterChargeCut, "Cut on Center Charge", "WriteDelete");
  file->WriteTObject(histMCPAmp, "MCP Amp", "WriteDelete");
  file->WriteTObject(histMCPAmpCut, "Cut on MCP Amp", "WriteDelete");

  file->Close();
  delete file;

}




void PlotDeltaTPDF(TCanvas *c, TLatex *tex, TH1F *hist, string outfile) {
  hist->Draw();
  gStyle->SetOptFit(0); //Hides the parameter box
  gStyle->SetOptStat(0);
  double mean = hist->GetMean();
  double rms = hist->GetRMS();
  TF1 *gausfit = new TF1("gausfit","gaus", mean - 2.0*rms, mean + 2.0*rms);//1-D gaus function defined around hist peak
  hist->Fit("gausfit","QMLES","", mean - 2.0*rms, mean + 2.0*rms);// Fit the hist; Q-quiet, L-log likelihood method, E-Minos errors technique, M-improve fit results
  hist->GetXaxis()->SetTitle("Time Resolution [ns]");
  tex->DrawLatex(0.6, 0.8, Form("#sigma = %.1f #pm %.1f ps", 1000*gausfit->GetParameter(2), 1000*gausfit->GetParError(2)));
  c->SaveAs(outfile.c_str()); //outfile should end in .pdf
}


void makeTimeResolution( string filename ) {

  DoMultiDeviceStudy( filename.c_str() );

  TFile *_file = TFile::Open( ("output"+filename).c_str() ); //Should be .root

  //Create variables containing hists:
  TH1F *histDeltaTCenter = (TH1F*)_file->Get("histDeltaTCenter"); //Picosil center pixel
  TH1F *histDeltaTPicoSilTotalCharge = (TH1F*)_file->Get("histDeltaTPicoSilTotalCharge");//Picosil, total charge
  TH1F *histDeltaTPicoSilEventCharge = (TH1F*)_file->Get("histDeltaTPicoSilEventCharge");//Picosil, event charge
  TH1F *histDeltaTPicoSilLandauCharge= (TH1F*)_file->Get("histDeltaTPicoSilLandauCharge");//Picosil, landau charge
  TH1F *histDeltaTPicoSilLandauChargeSmear = (TH1F*)_file->Get("histDeltaTPicoSilLandauChargeSmear");//Picosil, landau charge
  TH1F *histDeltaTMCP = (TH1F*)_file->Get("histDeltaTMCP"); // MCP
  TH1F *histDeltaT_Center_MCP_Equal = (TH1F*)_file->Get("histDeltaT_Center_MCP_Equal"); // each device weighted equally
  TH1F *histDeltaT_PicoSil_MCP_Equal = (TH1F*)_file->Get("histDeltaT_PicoSil_MCP_Equal");
  TH1F *histDeltaT_PicoSilEqual_MCP_Equal = (TH1F*)_file->Get("histDeltaT_PicoSilEqual_MCP_Equal");
  TH1F *histDeltaT_PicoSilEventCharge_MCP_Equal = (TH1F*)_file->Get("histDeltaT_PicoSilEventCharge_MCP_Equal");
  TH1F *histDeltaT_PicoSilTotalCharge_MCP_Equal = (TH1F*)_file->Get("histDeltaT_PicoSilTotalCharge_MCP_Equal");
  TH1F *histDeltaT_PicoSilLandauCharge_MCP_Equal= (TH1F*)_file->Get("histDeltaT_PicoSilLandauCharge_MCP_Equal");
  TH1F *histDeltaT_PicoSilLandauCharge_MCP_Equal_Smear = (TH1F*)_file->Get("histDeltaT_PicoSilLandauCharge_MCP_Equal_Smear");
  TH1F *histDeltaT_PicoSil_MCP_EventCharge = (TH1F*)_file->Get("histDeltaT_PicoSil_MCP_EventCharge"); //Combination of device Delta T's after shifting distributions around 0 and then weighting event by event.
  TH1F *histDeltaT_PicoSil_MCP_TotalCharge = (TH1F*)_file->Get("histDeltaT_PicoSil_MCP_TotalCharge"); //Combination after shifting around 0 and weighting with total charge.

  TCanvas *c = new TCanvas ("c","c",800, 600); 
  TLatex *tex = new TLatex();
  tex->SetNDC(); // Sets coords such that (0,0) is bottom left & (1,1) is top right.
  tex->SetTextSize(0.060);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);

  c->cd();

  histDeltaTCenter->SetTitle("HGC Center Pixel: TOF");
  histDeltaTPicoSilEventCharge->SetTitle("HGC: TOF w/ Event Charge Weighting");
  histDeltaTPicoSilTotalCharge->SetTitle("HGC: TOF w/ Total Charge Weighting");
  histDeltaTPicoSilLandauCharge->SetTitle("HGC: TOF w/ Landau MPV Charge Weighting");
  histDeltaTPicoSilLandauChargeSmear->SetTitle("SKIROC Emulation: HGC TOF w/ Landau MPV Charge Weighting");
  histDeltaTMCP->SetTitle("MCP: TOF");
  histDeltaT_Center_MCP_Equal->SetTitle("HGC Center Pixel, MCP: TOF w/ Equal Weighting");
  histDeltaT_PicoSil_MCP_Equal->SetTitle("HGC Inner Ring and Center Pixel, MCP: TOF w/ Weighting 1/8");
  histDeltaT_PicoSilEqual_MCP_Equal->SetTitle("HGC Inner Ring and Center Pixel Weighted 1/14, MCP Weighted 1/2: TOF");
  histDeltaT_PicoSil_MCP_EventCharge->SetTitle("HGC, MCP: TOF w/ Event Charge Weighting");
  histDeltaT_PicoSil_MCP_TotalCharge->SetTitle("HGC, MCP: TOF w/ Total Charge Weighting");
  histDeltaT_PicoSilEventCharge_MCP_Equal->SetTitle("HGC w/ Event Charge Weighting then *1/2, MCP Weighted 1/2: TOF");
  histDeltaT_PicoSilTotalCharge_MCP_Equal->SetTitle("HGC w/ Total Charge Weighting then *1/2, MCP Weighted 1/2: TOF");
  histDeltaT_PicoSilLandauCharge_MCP_Equal->SetTitle("HGC w/ Landau MPV Charge Weighting then *1/2, MCP Weighted 1/2: TOF");
  histDeltaT_PicoSilLandauCharge_MCP_Equal_Smear->SetTitle("SKIROC Emulation: HGC w/ Landau MPV Charge Weighting then *1/2, MCP Weighted 1/2: TOF");


  PlotDeltaTPDF(c, tex, histDeltaTCenter, "deltaTCenter.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTPicoSilEventCharge, "deltaTPicoSilEventCharge.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTPicoSilTotalCharge, "deltaTPicoSilTotalCharge.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTPicoSilLandauCharge, "deltaTPicoSilLandauCharge.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTPicoSilLandauChargeSmear, "deltaTPicoSilLandauChargeSKIROC.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTMCP, "deltaTMCP.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_Center_MCP_Equal, "deltaT_Center_MCP_Equal.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSil_MCP_Equal, "deltaT_PicoSil_MCP_Equal.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSilEqual_MCP_Equal, "deltaT_PicoSilEqual_MCP_Equal.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSil_MCP_EventCharge, "deltaT_PicoSil_MCP_EventCharge.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSil_MCP_TotalCharge, "deltaT_PicoSil_MCP_TotalCharge.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSilEventCharge_MCP_Equal, "deltaT_PicoSilEventCharge_MCP_Equal.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSilTotalCharge_MCP_Equal, "deltaT_PicoSilTotalCharge_MCP_Equal.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSilLandauCharge_MCP_Equal, "deltaT_PicoSilLandauCharge_MCP_Equal.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSilLandauCharge_MCP_Equal_Smear, "deltaT_PicoSilLandauCharge_MCP_Equal_SKIROC.pdf");
}


void MultiDeviceStudy_PicosilMCP() {
  makeTimeResolution("104-116except111-114.root"); // Outputs PDFs with histograms
}
