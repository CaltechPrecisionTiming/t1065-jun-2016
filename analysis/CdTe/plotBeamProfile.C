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
//#include "SiliconPadUtils.h"




void plotBeamProfile() {

  TFile *file_100GeV = TFile::Open( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v5_CdTe/analysis_5568.root" , "READ");
  TTree *tree_100GeV = (TTree*)(file_100GeV->Get("t1065"));

  TH2F *beam_100GeV = new TH2F("beam_100GeV", " ; Beam X Position [mm] ; Beam Y Position [mm]; Number of Events", 35,-20,15,35,-10,25);
  TH2F *beamSensor_100GeV = new TH2F("beamSensor_100GeV", " ; Beam X Position [mm] ; Beam Y Position [mm]; Amplitude Weighted Events", 35,-20,15,35,-10,25);
  TH2F *beamXVsAmp_100GeV = new TH2F("beamXVsAmp_100GeV", " ; Beam X Position [mm] ; Signal Amplitude [V]; Number of Events", 30,-5,25,40,0,1.0);
  TH2F *beamYVsAmp_100GeV = new TH2F("beamYVsAmp_100GeV", " ; Beam Y Position [mm] ; Signal Amplitude [V]; Number of Events", 30,-15,15,40,0,1.0);
  tree_100GeV->Draw("TDCx:TDCy>>beam_100GeV","","colz");
  tree_100GeV->Draw("TDCx:TDCy>>beamSensor_100GeV","amp[1]","colz");
  tree_100GeV->Draw("amp[1]*3.16228*(1.0/63.0957):TDCx>>beamXVsAmp_100GeV","","colz");
  tree_100GeV->Draw("amp[1]*3.16228*(1.0/63.0957):TDCy>>beamYVsAmp_100GeV","","colz");


  TCanvas *cv = new TCanvas("cv","cv", 800, 800);
  beam_100GeV->Draw("colz");
  beam_100GeV->SetStats(false);
  beam_100GeV->GetYaxis()->SetTitleOffset(1.5);
  beam_100GeV->GetXaxis()->SetTitleOffset(1.2);

  cv->SetRightMargin(0.15);
  cv->SetLeftMargin(0.12);
  cv->SetBottomMargin(0.12);

  cv->SaveAs("BeamProfile.png");
  cv->SaveAs("BeamProfile.pdf");

  cv = new TCanvas("cv","cv", 800, 800);
  beamSensor_100GeV->Draw("colz");
  beamSensor_100GeV->SetStats(false);
  beamSensor_100GeV->GetYaxis()->SetTitleOffset(1.5);
  beamSensor_100GeV->GetXaxis()->SetTitleOffset(1.2);

  cv->SetRightMargin(0.15);
  cv->SetLeftMargin(0.12);
  cv->SetBottomMargin(0.12);

  cv->SaveAs("BeamSensorProfile.png");
  cv->SaveAs("BeamSensorProfile.pdf");

  cv = new TCanvas("cv","cv", 800, 800);
  TProfile *beamXProfile_100GeV = beamXVsAmp_100GeV->ProfileX();
  beamXProfile_100GeV->Draw();
  beamXProfile_100GeV->SetLineWidth(2);
  beamXProfile_100GeV->SetStats(false);
  beamXProfile_100GeV->GetYaxis()->SetTitleOffset(2.5);
  beamXProfile_100GeV->GetYaxis()->SetTitle("Mean Signal Amplitude [V]");
  beamXProfile_100GeV->GetXaxis()->SetTitleOffset(1.2);

  cv->SetLeftMargin(0.2);
  cv->SetBottomMargin(0.12);

  cv->SaveAs("SensorXProfile.png");
  cv->SaveAs("SensorXProfile.pdf");

  cv = new TCanvas("cv","cv", 800, 800);
  TProfile *beamYProfile_100GeV = beamYVsAmp_100GeV->ProfileX();
  beamYProfile_100GeV->Draw();
  beamYProfile_100GeV->SetLineWidth(2);
  beamYProfile_100GeV->SetStats(false);
  beamYProfile_100GeV->GetYaxis()->SetTitleOffset(2.5);
  beamYProfile_100GeV->GetYaxis()->SetTitle("Mean Signal Amplitude [V]");
  beamYProfile_100GeV->GetXaxis()->SetTitleOffset(1.2);

  cv->SetLeftMargin(0.2);
  cv->SetBottomMargin(0.12);

  cv->SaveAs("SensorYProfile.png");
  cv->SaveAs("SensorYProfile.pdf");



}
