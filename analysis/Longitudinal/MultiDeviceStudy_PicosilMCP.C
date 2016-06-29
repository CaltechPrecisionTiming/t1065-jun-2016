
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
// An update to the code includes use of the first ring pixels in the picosil. -Daniel

//Histogram names start with histDeltaT. Then they say the devices they incorporate: either MCP, PicoSil center pixel, or PicoSil with all pixels. The PicoSil with all the Pixels will either specify equal, total charge, or event charge; this signifies all pixels weighted equally (1/7), or by charge (either weighting differently event-by-event or the same weights using total charge in the each pixel throughout the run). Then the histogram will specify equal, total charge, or event charge, with the same meanings as above, except for weighting the PicoSil delta T against the MCP delta T. Finally, a last option No Shift indicates the component histograms hadn't been shifted to zero by having their means subtracted, which should result in a histogram with a high timing resolution.
// Having At0 indicates the mean has been subtracted from the original hist.

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

  //Create histograms
  TH1F *histDeltaT_Center_MCP_Equal = new TH1F("histDeltaT_Center_MCP_Equal","; Time [ns];Number of Events", 50, -.3, .3); // Weights MCP and PicoSil center pixel 50-50
  TH1F *histDeltaT_PicoSilEventCharge_MCP_Equal = new TH1F("histDeltaT_PicoSilEventCharge_MCP_Equal","; Time [ns];Number of Events", 50, -.3, .3);// Picosil delta T found by weighting with event pixel charge, but then overall delta T fund by weighting PicoSil and MCP equally.
  TH1F *histDeltaT_PicoSilTotalCharge_MCP_Equal = new TH1F("histDeltaT_PicoSilTotalCharge_MCP_Equal","; Time [ns];Number of Events", 50, -.3, .3);// Picosil delta T found by weighting with total pixel charge, but then overall delta T fund by weighting PicoSil and MCP equally.
  TH1F *histDeltaT_PicoSilEqual_MCP_Equal = new TH1F("histDeltaT_PicoSilEqual_MCP_Equal","; Time [ns];Number of Events", 50, -.3, .3);//Picosil delta T found by weighting pixels equally, and then overall delta T fund by weighting PicoSil and MCP equally. Thus MCP is weighted by 1/2 and each pixel is weighted 1/14.
  TH1F *histDeltaT_PicoSil_MCP_Equal = new TH1F("histDeltaT_PicoSil_MCP_Equal","; Time [ns];Number of Events", 50, -.3, .3);// delta T found by placing equal emphasis on the pixels as on the MCP. Thus everything is weighted 1/8.
  TH1F *histDeltaT_PicoSil_MCP_EventCharge = new TH1F("histDeltaT_PicoSil_MCP_EventCharge","; Time [ns];Number of Events", 50, -.3, .3);// delta T found by weighting with event charge in each device.
  TH1F *histDeltaT_PicoSil_MCP_TotalCharge = new TH1F("histDeltaT_PicoSil_MCP_TotalCharge","; Time [ns];Number of Events", 50, -.3, .3);// delta T found by constant weighting with charge in device over entire run.
  TH1F *histDeltaT_Center_MCP_EventCharge_NoShift = new TH1F("histDeltaT_Center_MCP_EventCharge_NoShift","; Time [ns];Number of Events", 100, 2, 6); //Weights each event based on charge. Poor resolution.
  TH1F *histDeltaTPicoSil[7];
  for(int i=0; i<7; i++) histDeltaTPicoSil[i] = new TH1F(Form("histDeltaTPicoSil_%d",i),"; Time [ns];Number of Events", 50, 3, 6); //DeltaT of PicoSil pixels
  TH1F *histDeltaTCenter = new TH1F("histDeltaTCenter","; Time [ns];Number of Events", 50, 4, 5); //DeltaT of center picosil pixel
  TH1F *histDeltaTMCP = new TH1F("histDeltaTMCP","; Time [ns];Number of Events", 50, 2, 3); //DeltaT of MCP
  TH1F *histDeltaTPicoSilAt0 = new TH1F("histDeltaTPicoSilAt0","; Time [ns];Number of Events", 50, -0.3, 0.3); //All pixels combined
  TH1F *histDeltaTCenterAt0 = new TH1F("histDeltaTCenterAt0","; Time [ns];Number of Events", 50, -0.3, 0.3); //shifted to be centered at zero
  TH1F *histDeltaTMCPAt0 = new TH1F("histDeltaTMCPAt0","; Time [ns];Number of Events", 50, -0.3, 0.3); //shifted to be centered at zero
  

  float totalPicoSilCharge[7] = {0.};
  float totalCenterCharge = 0; // Should be same value as totalPicoSilCharge[0]
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
    // (Will be used to fill the histogram at every event)
    float DeltaTCenter = photekTimeGauss0 - centerTime; 
    float DeltaTMCP = photekTimeGauss1 - MCPTime;
    float DeltaTPicoSil[7] = {-99.}; 
    for ( int j = 1; j <= 7; j++){
      if ( amp[j] > 0.02 && integral[j] > 1 ) {
	DeltaTPicoSil[j-1] = photekTimeGauss0 - linearTime45[j];
	totalPicoSilCharge[j-1] += integral[j];
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
      if (DeltaTPicoSil[k] != -99.) histDeltaTPicoSil[k]->Fill(DeltaTPicoSil[k]);
    }
  }

  totalPicoSilCharge[0] *= 2; //Account for 6dB attenuator
  totalCenterCharge *= 2;

  float ringChargeTotal = 0;
  for(int i=1; i<=6; i++) ringChargeTotal += totalPicoSilCharge[i];

  // Do Gaussian fit of delta T distributions from (mean-2RMS) to (mean+2RMS)
  double mean = histDeltaT_Center_MCP_EventCharge_NoShift->GetMean();
  double rms = histDeltaT_Center_MCP_EventCharge_NoShift->GetRMS();
  double xmin = mean-2.0*rms;
  double xmax = mean+2.0*rms;
  histDeltaT_Center_MCP_EventCharge_NoShift->Fit("gaus","MLES","",xmin,xmax);

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

  double meanPicoSil[7];
  for (int i = 0; i < 7; i++){
    meanPicoSil[i] = histDeltaTPicoSil[i]->GetMean();
  }


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

    float DeltaTPicoSil[7] = {0.}; 
    for ( int j = 1; j <= 7; j++){
      if ( amp[j] > 0.02 && integral[j] > 1 ) {
	DeltaTPicoSil[j-1] = photekTimeGauss0 - linearTime45[j] - meanPicoSil[j-1];
      }
    }

    float ringWeightTotal = 0;
    for(int ii = 1; ii <= 6; ii++)ringWeightTotal += totalPicoSilCharge[ii] * DeltaTPicoSil[ii];

    float ringWeightEvent = 0;
    for (int jj = 2; jj <= 7; jj++) ringWeightEvent += integral[jj] * DeltaTPicoSil[jj-1];

    float ringChargeEvent = 0;
    for (int kk = 2; kk<=7; kk++) {if(DeltaTPicoSil[kk-1] != 0.) ringChargeEvent += integral[kk];}

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
    for (int j = 0; j <= 6; j++) temp_pseventweight += DeltaTPicoSil[j]*integral[j+1];
    float temp_pseventcharge = 0;
    for (int j = 0; j <= 6; j++) {if (DeltaTPicoSil[j] != 0) temp_pseventcharge += integral[j+1];}
    DeltaT_PicoSilEventCharge_MCP_Equal += 0.5*temp_pseventweight/temp_pseventcharge;

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

    histDeltaT_Center_MCP_Equal->Fill(DeltaT_Center_MCP_Equal);
    histDeltaT_PicoSil_MCP_EventCharge->Fill(DeltaT_PicoSil_MCP_EventCharge); 
    histDeltaT_PicoSil_MCP_TotalCharge->Fill(DeltaT_PicoSil_MCP_TotalCharge); 
    histDeltaT_PicoSilTotalCharge_MCP_Equal->Fill(DeltaT_PicoSilTotalCharge_MCP_Equal);
    histDeltaT_PicoSilEventCharge_MCP_Equal->Fill(DeltaT_PicoSilEventCharge_MCP_Equal);
    histDeltaT_PicoSil_MCP_Equal->Fill(DeltaT_PicoSil_MCP_Equal);
    histDeltaT_PicoSilEqual_MCP_Equal->Fill(DeltaT_PicoSilEqual_MCP_Equal);

    histDeltaTPicoSilAt0->Fill(ringWeightTotal/ringChargeTotal);
    histDeltaTCenterAt0->Fill(DeltaTCenter);
    histDeltaTMCPAt0->Fill(DeltaTMCP);
  }

  mean = histDeltaT_Center_MCP_Equal->GetMean();
  rms = histDeltaT_Center_MCP_Equal->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_Center_MCP_Equal->Fit("gaus","MLES","",xmin,xmax);

  mean = histDeltaT_PicoSil_MCP_Equal->GetMean();
  rms = histDeltaT_PicoSil_MCP_Equal->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_PicoSil_MCP_Equal->Fit("gaus","MLES","",xmin,xmax);

  mean = histDeltaT_PicoSilEqual_MCP_Equal->GetMean();
  rms = histDeltaT_PicoSilEqual_MCP_Equal->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_PicoSilEqual_MCP_Equal->Fit("gaus","MLES","",xmin,xmax);

  mean = histDeltaT_PicoSil_MCP_EventCharge->GetMean();
  rms = histDeltaT_PicoSil_MCP_EventCharge->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_PicoSil_MCP_EventCharge->Fit("gaus","MLES","",xmin,xmax);

  mean = histDeltaT_PicoSil_MCP_TotalCharge->GetMean();
  rms = histDeltaT_PicoSil_MCP_TotalCharge->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_PicoSil_MCP_TotalCharge->Fit("gaus","MLES","",xmin,xmax);

  mean = histDeltaT_PicoSilTotalCharge_MCP_Equal->GetMean();
  rms = histDeltaT_PicoSilTotalCharge_MCP_Equal->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_PicoSilTotalCharge_MCP_Equal->Fit("gaus","MLES","",xmin,xmax);

  mean = histDeltaT_PicoSilEventCharge_MCP_Equal->GetMean();
  rms = histDeltaT_PicoSilEventCharge_MCP_Equal->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaT_PicoSilEventCharge_MCP_Equal->Fit("gaus","MLES","",xmin,xmax);

  mean = histDeltaTCenterAt0->GetMean();
  rms = histDeltaTCenterAt0->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaTCenterAt0->Fit("gaus","MLES","",xmin,xmax);

  mean = histDeltaTMCPAt0->GetMean();
  rms = histDeltaTMCPAt0->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaTMCPAt0->Fit("gaus","MLES","",xmin,xmax);

  mean = histDeltaTPicoSilAt0->GetMean();
  rms = histDeltaTPicoSilAt0->GetRMS();
  xmin = mean-2.0*rms;
  xmax = mean+2.0*rms;
  histDeltaTPicoSilAt0->Fit("gaus","MLES","",xmin,xmax);


  // Creates output root file
  TFile *file = TFile::Open(outputFilename.c_str(), "RECREATE");
  file->cd();
  file->WriteTObject(histDeltaTCenterAt0,"histDeltaTCenter", "WriteDelete");
  file->WriteTObject(histDeltaTPicoSilAt0,"histDeltaTPicoSil", "WriteDelete");
  file->WriteTObject(histDeltaTMCPAt0,"histDeltaTMCP", "WriteDelete");
  file->WriteTObject(histDeltaT_Center_MCP_EventCharge_NoShift,"histDeltaT_Center_MCP_EventCharge_NoShift", "WriteDelete");
  file->WriteTObject(histDeltaT_Center_MCP_Equal,"histDeltaT_Center_MCP_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSil_MCP_Equal,"histDeltaT_PicoSil_MCP_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSilEqual_MCP_Equal,"histDeltaT_PicoSilEqual_MCP_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSilEventCharge_MCP_Equal,"histDeltaT_PicoSilEventCharge_MCP_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSilTotalCharge_MCP_Equal,"histDeltaT_PicoSilTotalCharge_MCP_Equal", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSil_MCP_EventCharge,"histDeltaT_PicoSil_MCP_EventCharge", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSil_MCP_TotalCharge,"histDeltaT_PicoSil_MCP_TotalCharge", "WriteDelete");

  file->Close();
  delete file;

}




void makeTimeResolution( string filename ) {

  TFile *_file = TFile::Open( filename.c_str() ); //Should be .root

  //Create variables containing hists:
  TH1F* histDeltaT = (TH1F*)_file->Get("histDeltaT_Center_MCP_Equal"); // each device weighted equally
  // Not doing histDeltaTWeighted, since it gives bad results.
  TH1F* histDeltaTCenter = (TH1F*)_file->Get("histDeltaTCenter"); //Picosil center pixel
  TH1F* histDeltaTMCP = (TH1F*)_file->Get("histDeltaTMCP"); // MCP
  TH1F* histDeltaTWeightsCorrected = (TH1F*)_file->Get("histDeltaT_PicoSil_MCP_EventCharge"); //Combination of device Delta T's after shifting distributions around 0 and then weighting event by event.
  TH1F* histDeltaTWeightsCorrected_totalcharge = (TH1F*)_file->Get("histDeltaT_PicoSil_MCP_TotalCharge"); //Combination after shifting around 0 and weighting with total charge.

  TCanvas *c = new TCanvas ("c","c",800, 600); 
  TLatex *tex = new TLatex();
  tex->SetNDC(); // Sets coords such that (0,0) is bottom left & (1,1) is top right.
  tex->SetTextSize(0.060);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);

  c->cd();
  gStyle->SetOptStat(0); // Hides the parameter box


  histDeltaT->Draw();
  TF1* fgausEqual = new TF1("fgausEqual","gaus", histDeltaT->GetMean() - 2.0*histDeltaT->GetRMS(), histDeltaT->GetMean() + 2.0*histDeltaT->GetRMS()); // 1-D gaus func defined around hist peak
  histDeltaT->Fit("fgausEqual","QMLE","", histDeltaT->GetMean() - 2.0*histDeltaT->GetRMS(), histDeltaT->GetMean() + 2.0*histDeltaT->GetRMS()); // Fit the hist; Q-quiet, L-log likelihood method, E-Minos errors technique, M-improve fit results
  histDeltaT->GetXaxis()->SetTitle("Time Resolution [ps]");
  histDeltaT->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgausEqual->GetParameter(2), 1000*fgausEqual->GetParError(2)));
  c->SaveAs("deltaTEqual.pdf");

  histDeltaTCenter->Draw();
  TF1* fgausCenter = new TF1("fgausCenter","gaus", histDeltaTCenter->GetMean() - 2.0*histDeltaTCenter->GetRMS(), histDeltaTCenter->GetMean() + 2.0*histDeltaTCenter->GetRMS()); 
  histDeltaTCenter->Fit("fgausCenter","QMLE","", histDeltaTCenter->GetMean() - 2.0*histDeltaTCenter->GetRMS(), histDeltaTCenter->GetMean() + 2.0*histDeltaTCenter->GetRMS());
  histDeltaTCenter->GetXaxis()->SetTitle("Time Resolution [ps]");
  histDeltaTCenter->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgausCenter->GetParameter(2), 1000*fgausCenter->GetParError(2)));
  c->SaveAs("deltaTCenter.pdf");

  histDeltaTMCP->Draw();
  TF1* fgausMCP = new TF1("fgausMCP","gaus", histDeltaTMCP->GetMean() - 2.0*histDeltaTMCP->GetRMS(), histDeltaTMCP->GetMean() + 2.0*histDeltaTMCP->GetRMS()); 
  histDeltaTMCP->Fit("fgausMCP","QMLE","", histDeltaTMCP->GetMean() - 2.0*histDeltaTMCP->GetRMS(), histDeltaTMCP->GetMean() + 2.0*histDeltaTMCP->GetRMS());
  histDeltaTMCP->GetXaxis()->SetTitle("Time Resolution [ps]");
  histDeltaTMCP->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgausMCP->GetParameter(2), 1000*fgausMCP->GetParError(2)));
  c->SaveAs("deltaTMCP.pdf");
  
  histDeltaTWeightsCorrected->Draw();
  TF1* fgausWC = new TF1("fgausWC","gaus", histDeltaTWeightsCorrected->GetMean() - 2.0*histDeltaTWeightsCorrected->GetRMS(), histDeltaTWeightsCorrected->GetMean() + 2.0*histDeltaTWeightsCorrected->GetRMS()); 
  histDeltaTWeightsCorrected->Fit("fgausWC","QMLE","", histDeltaTWeightsCorrected->GetMean() - 2.0*histDeltaTWeightsCorrected->GetRMS(), histDeltaTWeightsCorrected->GetMean() + 2.0*histDeltaTWeightsCorrected->GetRMS());
  histDeltaTWeightsCorrected->GetXaxis()->SetTitle("Time Resolution [ps]");
  histDeltaTWeightsCorrected->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgausWC->GetParameter(2), 1000*fgausWC->GetParError(2)));
  c->SaveAs("deltaTWeightsCorrected.pdf");

  histDeltaTWeightsCorrected_totalcharge->Draw();
  TF1* fgausWCtc = new TF1("fgausWCtc","gaus", histDeltaTWeightsCorrected_totalcharge->GetMean() - 2.0*histDeltaTWeightsCorrected_totalcharge->GetRMS(), histDeltaTWeightsCorrected_totalcharge->GetMean() + 2.0*histDeltaTWeightsCorrected_totalcharge->GetRMS()); 
  histDeltaTWeightsCorrected_totalcharge->Fit("fgausWCtc","QMLE","", histDeltaTWeightsCorrected_totalcharge->GetMean() - 2.0*histDeltaTWeightsCorrected_totalcharge->GetRMS(), histDeltaTWeightsCorrected_totalcharge->GetMean() + 2.0*histDeltaTWeightsCorrected_totalcharge->GetRMS());
  histDeltaTWeightsCorrected_totalcharge->GetXaxis()->SetTitle("Time Resolution [ps]");
  histDeltaTWeightsCorrected_totalcharge->SetTitle("");
  tex->DrawLatex(0.6, 0.80, Form("#sigma = %.0f #pm %.2f ps", 1000*fgausWCtc->GetParameter(2), 1000*fgausWCtc->GetParError(2)));
  c->SaveAs("deltaTWeightsCorrected_totalcharge.pdf");
  
}


void MultiDeviceStudy_PicosilMCP() {
  DoMultiDeviceStudy("65-83.root","output65-83.root");
  makeTimeResolution("output65-83.root"); // Outputs PDFs with histograms
}
