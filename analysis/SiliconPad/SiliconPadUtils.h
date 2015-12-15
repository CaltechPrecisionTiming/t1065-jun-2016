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

void MakeChargePlot(string filename, string plotname, double ampCutOnPhotek, double attenuationFactor, 
		    double xmin, double xmax, 
		    double fitmin, double fitmax) {
  // Get the tree

  TFile *inputfile = TFile::Open(filename.c_str(),"READ");
  TTree *tree = (TTree*)inputfile->Get("pulse");

  // get the variables from the ntuple
  float amp[36];
  float integral[36];
  float gauspeak[36];

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("gauspeak",1);
  tree->SetBranchStatus("amp",1);
  tree->SetBranchStatus("int",1);
  tree->SetBranchAddress("gauspeak",gauspeak);
  tree->SetBranchAddress("amp",amp);
  tree->SetBranchAddress("int",integral);

  //create histograms
  TH1F *histIntCharge;
  histIntCharge = new TH1F("histIntCharge","; Integrated Charge [pC];Number of Events", 50, xmin,xmax);
  
  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();

  std::cout<<"Number of events in Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
    if (iEntry %1000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
  
    float photekTimeGauss = gauspeak[18];
    float siliconTimeGauss = gauspeak[21];
    float photekAmp = amp[18];
    float siliconAmp = amp[21];
    float photekIntegral = integral[18];
    float siliconIntegral = integral[21];
       
    //use photek amplitude cut for electron ID
    //cout << "test: " << photekAmp << " " << siliconIntegral << "\n";
    if( !(photekAmp > ampCutOnPhotek)) continue;

    
    //Fill histogram
    histIntCharge->Fill(siliconIntegral * attenuationFactor);

  }


  TCanvas * c = 0;


  //Energy plot
  c = new TCanvas("c","c",600,600);  
  histIntCharge->SetAxisRange(xmin,xmax,"X");
  histIntCharge->SetTitle("");
  histIntCharge->GetXaxis()->SetTitle("Integrated Charge [pC]");
  histIntCharge->GetYaxis()->SetTitle("Number of Events");
  histIntCharge->GetYaxis()->SetTitleOffset(1.3);
  histIntCharge->SetMaximum(1.2*histIntCharge->GetMaximum());
  histIntCharge->Draw();
  histIntCharge->SetStats(0);
  histIntCharge->Fit("gaus","","",fitmin,fitmax);
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.040);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.45, 0.85, Form("Mean = %.2f #pm %.2f %s",fitter->GetParameter(1),TMath::Max(0.01,fitter->GetParError(1)),"pC"));
  tex->DrawLatex(0.45, 0.80, Form("#sigma = %.2f #pm %.2f %s",fitter->GetParameter(2),TMath::Max(0.01,fitter->GetParError(2)),"pC"));
  //tex->DrawLatex(0.15, 0.92, Form("Attenuation Factor = %.3f",attenuationFactor));
  
  c->SaveAs( Form("%s_charge.gif", plotname.c_str()) );
  c->SaveAs( Form("%s_charge.pdf", plotname.c_str()) );
 
}


