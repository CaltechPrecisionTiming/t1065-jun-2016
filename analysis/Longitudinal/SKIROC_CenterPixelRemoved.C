
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
#include <algorithm>
#include <functional>
#include "TRandom3.h"

// SKIROC emulation with center HGC pixel removed/not taken into account.
// Author: Daniel Gawerc

string pixels_added[6];


void Fitter(TH1F *hist) {
  //Helper function for fitting Gaussian
  double xmin = hist->GetMean() - 2.0*hist->GetRMS();
  double xmax = hist->GetMean() + 2.0*hist->GetRMS();
  hist->Fit("gaus","QMLES","",xmin,xmax); // Q suppresses fit results
  gStyle->SetOptFit(1);
}

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
  float smearWidth = 1.5;
  float smearBins = 45;

  TH1F *histDeltaT_PicoSilEqual_MCP_Equal_BothSmear = new TH1F("histDeltaT_PicoSilEqual_MCP_Equal_BothSmear","; Time [ns];Number of Events", smearBins, -smearWidth, smearWidth);
  TH1F *histDeltaTPicoSilSmear[6];
  TH1F *histDeltaTPicoSilSmearAt0[6];
  for(int i=0; i<6; i++) {
    histDeltaTPicoSilSmear[i] = new TH1F(Form("histDeltaTPicoSilSmear_%d",i),"; Time [ns];Number of Events", 50, 3, 6);
    histDeltaTPicoSilSmearAt0[i] = new TH1F(Form("histDeltaTPicoSilSmearAt0_%d",i),"; Time [ns];Number of Events", smearBins, -smearWidth, smearWidth); 
  }
  TH1F *histDeltaTMCPSmear = new TH1F("histDeltaTMCPSmear","; Time [ns];Number of Events", 50, 2, 3); 
  TH1F *histDeltaTMCPAt0Smear = new TH1F("histDeltaTMCPAt0Smear","; Time [ns];Number of Events", smearBins, -smearWidth, smearWidth); //shifted to be centered at zero
  TH1F *histDeltaTPicoSilAt0EqualSmear = new TH1F("histDeltaTPicoSilAt0EqualSmear", "; Time [ns];Number of Events", smearBins, -smearWidth, smearWidth);
  TH1F *histDeltaTPicoSilAt0EqualSmear_nEventsCombine[6];
  for(int i=0; i<6; i++) histDeltaTPicoSilAt0EqualSmear_nEventsCombine[i] = 
      new TH1F(Form("histDeltaTPicoSilAt0EqualSmear_nEventsCombine_%dPixels",i+1), "; Time [ns];Number of Events", smearBins, -smearWidth, smearWidth);

  float photekAmpCut = sqrt(10)*0.1; //THESE ARE THE CUT VALUES AFTER ADJUSTING FOR ATTENUATORS
  float photekChargeCut = sqrt(10)*2;
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

    float MCPAmp = amp[11];
    float MCPCharge = integral[11];
    float MCPTime = linearTime45[11];


    // APPLY EVENT CUTS:
    //require photek (for electron selection)
    if( !(photekAmp0 > photekAmpCut && photekCharge0 > photekChargeCut) ) continue;
    //require MCP minimum amplitude. 
    if( !( MCPAmp > MCPAmpCut) ) continue;


    double linearTime45Smear[6];
    for (int j = 0; j < 6; j++)  linearTime45Smear[j] = rando->Gaus(linearTime45[j+2], 0.500); //Samples from smear
    double MCPTimeSmear = rando->Gaus(MCPTime, 0.350);


    //Calculates the Delta T's if the event passes the cuts:
    for ( int j = 0; j < 6; j++){
      if ( amp[j+2] > 0.01 && integral[j+2] > 1 ) {
        histDeltaTPicoSilSmear[j]->Fill( photekTimeGauss0 - linearTime45Smear[j] );
      }
    }

    histDeltaTMCPSmear->Fill(photekTimeGauss1 - MCPTimeSmear);
  }


  double meanMCPSmear = histDeltaTMCPSmear->GetMean();
  double meanPicoSilSmear[6];
  int DeltaTPicoSilSmear_Events[6];
  for (int i = 0; i < 6; i++){
    meanPicoSilSmear[i] = histDeltaTPicoSilSmear[i]->GetMean();
    DeltaTPicoSilSmear_Events[i] = histDeltaTPicoSilSmear[i]->GetEntries();
  }
  int DeltaTPicoSilSmear_Events_Sorted[6];
  int DeltaTPicoSilSmear_Events_SortedIndices[6];
  std::copy(DeltaTPicoSilSmear_Events, DeltaTPicoSilSmear_Events+6, DeltaTPicoSilSmear_Events_Sorted);
  std::sort(DeltaTPicoSilSmear_Events_Sorted, DeltaTPicoSilSmear_Events_Sorted+6, std::greater<int>() );
  for (int i = 0; i < 6; i++){
    for (int j = 0; j < 6; j++) {
      if( DeltaTPicoSilSmear_Events[j] == DeltaTPicoSilSmear_Events_Sorted[i] ) {
        DeltaTPicoSilSmear_Events_SortedIndices[i] = j;
        break;
      }
    }
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

    float MCPAmp = amp[11];
    float MCPCharge = integral[11];
    float MCPTime = linearTime45[11];

    if( !(photekAmp0 > photekAmpCut && photekCharge0 > photekChargeCut) ) continue;
    if( !( MCPAmp > MCPAmpCut) ) continue;

    double linearTime45Smear[6];
    for (int j = 0; j < 6; j++)  linearTime45Smear[j] = rando2->Gaus(linearTime45[j+2],0.500);
    double MCPTimeSmear = rando2->Gaus(MCPTime, 0.350);

    float DeltaTPicoSilSmear[6];
    std::fill(DeltaTPicoSilSmear, DeltaTPicoSilSmear+6, -99);
    for ( int j = 0; j < 6; j++){
      if ( amp[j+2] > 0.01 && integral[j+2] > 1 ) {
        DeltaTPicoSilSmear[j] = photekTimeGauss0 - linearTime45Smear[j] - meanPicoSilSmear[j];
        histDeltaTPicoSilSmearAt0[j]->Fill(DeltaTPicoSilSmear[j]);
      }
    }

    
    //Here are the corrections that center the histograms at 0.
    float DeltaTMCPSmear = photekTimeGauss1 - MCPTimeSmear - meanMCPSmear;


    float DeltaTPicoSilAt0EqualSmear = 0;
    int inc = 0;
    for (int j = 0; j < 6; j++) {
      if (DeltaTPicoSilSmear[j] != -99.) {
        DeltaTPicoSilAt0EqualSmear += DeltaTPicoSilSmear[j];
        inc += 1;
      }
    }
    if (inc != 0) DeltaTPicoSilAt0EqualSmear /= inc;
    else DeltaTPicoSilAt0EqualSmear = -99;


    // For adding in each pixel in descending order by number of events.
    for(int j = 0; j < 6; j++) {
      double fill = 0;
      inc = 0;
      for(int k = 0; k <= j; k++) {
        if( DeltaTPicoSilSmear[ DeltaTPicoSilSmear_Events_SortedIndices[k] ] != -99. ) {
          fill += DeltaTPicoSilSmear[ DeltaTPicoSilSmear_Events_SortedIndices[k] ];
          inc += 1;
        }
      }
      if (inc != 0) {
        fill /= inc;
        histDeltaTPicoSilAt0EqualSmear_nEventsCombine[j]->Fill( fill );
      }
    }


    if(DeltaTPicoSilAt0EqualSmear != -99) {
      histDeltaT_PicoSilEqual_MCP_Equal_BothSmear->Fill(0.5*DeltaTPicoSilAt0EqualSmear + 0.5*DeltaTMCPSmear);
      histDeltaTPicoSilAt0EqualSmear->Fill(DeltaTPicoSilAt0EqualSmear);
    }
    histDeltaTMCPAt0Smear->Fill(DeltaTMCPSmear);
  }

  // Add Gaussian fit
  Fitter(histDeltaT_PicoSilEqual_MCP_Equal_BothSmear);
  Fitter(histDeltaTMCPAt0Smear);
  Fitter(histDeltaTPicoSilAt0EqualSmear);

  for (int i = 0; i < 6; i++) {
    Fitter(histDeltaTPicoSilSmearAt0[i]);
    Fitter(histDeltaTPicoSilAt0EqualSmear_nEventsCombine[i]);
  }


  // Creates output root file
  TFile *file = TFile::Open(("output_nocenter"+filename).c_str(), "RECREATE");
  file->cd();
  file->WriteTObject(histDeltaTPicoSilAt0EqualSmear,"histDeltaTPicoSilEqualSmear", "WriteDelete");
  file->WriteTObject(histDeltaTMCPAt0Smear,"histDeltaTMCPSmear", "WriteDelete");
  file->WriteTObject(histDeltaT_PicoSilEqual_MCP_Equal_BothSmear,"histDeltaT_PicoSilEqual_MCP_Equal_BothSmear", "WriteDelete");
  for(int i=0; i<6; i++) {
    file->WriteTObject(histDeltaTPicoSilSmearAt0[i],Form("histDeltaTPicoSilSmear[%d]",i+1),"WriteDelete");
  }

  pixels_added[0] = Form("#Deltat of Equal-Weight HGC Ring 1 Pixels in Order of # Events. Smeared Pixels: %d", DeltaTPicoSilSmear_Events_SortedIndices[0]+1);
  for(int i=1; i<6; i++) pixels_added[i] = pixels_added[i-1] + Form(",%d",DeltaTPicoSilSmear_Events_SortedIndices[i]+1);
  for(int i=0; i<6; i++) {
    file->WriteTObject(histDeltaTPicoSilAt0EqualSmear_nEventsCombine[i], pixels_added[i].c_str(), "WriteDelete");
  }


  TH1F *histPhotekAmpCut = new TH1F("histPhotekAmpCut","; Amp;Number of Events", 400, 0, 2.5);
  TH1F *histPhotekChargeCut = new TH1F("histPhotekChargeCut","; Charge;Number of Events", 400, 0, 30);
  TH1F *histMCPAmpCut = new TH1F("histMCPAmpCut","; Amp;Number of Events", 100, 0, 0.75);

  tree->Draw("sqrt(10)*amp[0]>>histPhotekAmpCut", Form("sqrt(10)*amp[0]>%f",photekAmpCut) );
  tree->Draw("sqrt(10)*int[0]>>histPhotekChargeCut", Form("sqrt(10)*int[0]>%f",photekChargeCut));
  tree->Draw("amp[11]>>histMCPAmpCut", Form("amp[11]>%f",MCPAmpCut));

  TH1F *histPhotekAmp = new TH1F("histPhotekAmp","; Amp;Number of Events", 400, 0, 2.5);
  TH1F *histPhotekCharge = new TH1F("histPhotekCharge","; Charge;Number of Events", 400, 0, 30);
  TH1F *histMCPAmp = new TH1F("histMCPAmp","; Amp;Number of Events", 100, 0, 0.75);

  tree->Draw("sqrt(10)*amp[0]>>histPhotekAmp" );
  tree->Draw("sqrt(10)*int[0]>>histPhotekCharge");
  tree->Draw("amp[11]>>histMCPAmp");

  file->WriteTObject(histPhotekAmp, "Photek Amp", "WriteDelete");
  file->WriteTObject(histPhotekAmpCut, "Cut on Photek Amp", "WriteDelete");
  file->WriteTObject(histPhotekCharge, "Photek Charge", "WriteDelete");
  file->WriteTObject(histPhotekChargeCut, "Cut on Photek Charge", "WriteDelete");
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
  tex->DrawLatex(0.6, 0.8, Form("#sigma = %.0f #pm %.0f ps", 1000*gausfit->GetParameter(2), 1000*gausfit->GetParError(2)));
  c->SaveAs(outfile.c_str()); //outfile should end in .pdf
}



void makeTimeResolution( string filename ) {

  DoMultiDeviceStudy( filename.c_str() );

  TFile *_file = TFile::Open( ("output_nocenter"+filename).c_str() ); //Should be .root

  //Create variables containing hists:
  TH1F *histDeltaTPicoSilEqualSmear = (TH1F*)_file->Get("histDeltaTPicoSilEqualSmear");
  TH1F *histDeltaTMCPSmear = (TH1F*)_file->Get("histDeltaTMCPSmear"); // MCP
  TH1F *histDeltaT_PicoSilEqual_MCP_Equal_BothSmear = (TH1F*)_file->Get("histDeltaT_PicoSilEqual_MCP_Equal_BothSmear");
  TH1F *combo[6];
  for (int i=0; i<6; i++) combo[i] = (TH1F*)_file->Get( pixels_added[i].c_str() );

  TCanvas *c = new TCanvas ("c","c",800, 600); 
  TLatex *tex = new TLatex();
  tex->SetNDC(); // Sets coords such that (0,0) is bottom left & (1,1) is top right.
  tex->SetTextSize(0.060);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);

  c->cd();

  histDeltaTPicoSilEqualSmear->SetTitle("SKIROC Emulation: HGC Ring 1 TOF w/ Equal Weighting");
  histDeltaTMCPSmear->SetTitle("SKIROC Emulation: MCP TOF");
  histDeltaT_PicoSilEqual_MCP_Equal_BothSmear->SetTitle("1/12 HGC Smeared Ring 1 Pixels, 1/2 MCP Smeared: TOF");
  for (int i=0; i<6; i++) combo[i]->SetTitle( pixels_added[i].c_str() );

  PlotDeltaTPDF(c, tex, histDeltaTPicoSilEqualSmear, "deltaTPicoSilEqualSmear_NoCenter.pdf");
  PlotDeltaTPDF(c, tex, histDeltaTMCPSmear, "deltaTMCPSmear_NoCenter.pdf");
  PlotDeltaTPDF(c, tex, histDeltaT_PicoSilEqual_MCP_Equal_BothSmear, "deltaT_PicoSilEqual_MCP_Equal_BothSmear_NoCenter.pdf");
  for (int i=0; i<6; i++) PlotDeltaTPDF(c, tex, combo[i], Form("SKIROC_%d_Pixels_NoCenter.pdf",i+1) );
}


void SKIROC_CenterPixelRemoved() {
  makeTimeResolution("104-116except111-114.root"); // Outputs PDFs with histograms
}
