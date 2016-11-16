//#ifdef __MAKECINT__
//#pragma link C++ class vector<vector<float> >+;
//#endif

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
#include "SiliconPadUtils.h"




void plotPulses() {

  TFile *file_100GeV = TFile::Open( "eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v5/analysis_5568.root" , "READ");
  TTree *tree_100GeV = (TTree*)(file_100GeV->Get("t1065"));

  TH1F *pulse_100GeV = new TH1F("pulse_100GeV", ";time;amplitude", 1024, -0.5, 1023.5);
  tree_100GeV->Draw("t0>>pulse_100GeV", "raw[17]*(event==100)");
  
  TCanvas *cv = new TCanvas("cv","cv", 800, 800);
  pulse_100GeV->Draw("hist");



}
