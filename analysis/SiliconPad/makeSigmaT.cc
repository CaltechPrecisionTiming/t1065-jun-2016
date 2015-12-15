#ifdef __MAKECINT__
#pragma link C++ class vector<vector<float> >+;
#endif

#include <iostream>
#include <string>
#include <fstream> 
#include <sstream>
#include <vector>
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
#include "TGaxis.h"
#include "TPad.h"
#include <math.h> 

void MakeAmplitudePlot(std::string filename, std::string plotname, float electronIDcut) {
  // Get the tree


  TFile *inputfile = new TFile(filename.c_str(),"READ");
  TTree *tree = (TTree*)inputfile->Get("pulse");

  // get the variables from the ntuple
  Float_t amp[36];
  Float_t integral[36];
  Float_t gauspeak[36];

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("amp",1);
  tree->SetBranchStatus("int",1);
  tree->SetBranchStatus("gauspeak",1);

  tree->SetBranchAddress("amp",&amp);
  tree->SetBranchAddress("int",&integral);
  tree->SetBranchAddress("gauspeak",&gauspeak);

  //create histograms
  const int nBinsX = 75;
  const float Xrange = 6.5*nBinsX;

  TProfile * DeltaT_vs_Charge = new TProfile("DeltaT_vs_Charge", "DeltaT_vs_Charge", nBinsX, 0, Xrange, 2, 4);
  TH2F *DeltaT_vs_Charge_Corrected = new TH2F("DeltaT_vs_Charge_Corrected","; X; Y", nBinsX, 0, Xrange, 40, 2, 4);

  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();
  std::cout<<"Number of events in Physics Sample: "<<nentries<<std::endl;  

  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
    
    tree->GetEntry(iEntry);    

    if ( iEntry%1000 == 0 ) std::cout << "event: " << iEntry << std::endl;
    
    if(amp[18] > electronIDcut)
      DeltaT_vs_Charge->Fill(10*integral[21], gauspeak[18]-gauspeak[21]);
  }
  
  std::cout<<"Done Filling the raw plots... "<<std::endl;
  
  TF1* fslopecorr = new TF1("fslopecorr","pol1", 0, Xrange);
  DeltaT_vs_Charge->Fit("fslopecorr","Q","", 0, Xrange);
  
  float intercept = fslopecorr->GetParameter(0);
  float slope = fslopecorr->GetParameter(1);
  
  // correct the slope
  std::cout<<"Correcting the slope"<<std::endl;

  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) 
    {      
      tree->GetEntry(iEntry);    
      
      if ( iEntry%1000 == 0 ) std::cout << "event: " << iEntry << std::endl;

      float a = fslopecorr->GetParameter(0);
      float b = fslopecorr->GetParameter(1);

      if(amp[18] > electronIDcut)
	{
	  DeltaT_vs_Charge_Corrected->Fill(10*integral[21], gauspeak[18]-gauspeak[21] + (a+b*0.02) - (a+b*10*integral[21]));

	}
   }

  
  // make the sigmaT vs integral plot

  float res[nBinsX]; 
  float charge[nBinsX];
  float errorX[nBinsX], errorY[nBinsX];

  for(int i = 0; i < DeltaT_vs_Charge_Corrected->GetNbinsX(); i++)
    {
      TH1F* fHprojy = (TH1F*)DeltaT_vs_Charge_Corrected->ProjectionY("fHprojy", i, i+1);
 
      TF1* fgaus = new TF1("fgaus","gaus", fHprojy->GetMean() - 2*fHprojy->GetRMS(), fHprojy->GetMean() + 2*fHprojy->GetRMS());
      fHprojy->Fit("fgaus","Q","", fHprojy->GetMean() - 2*fHprojy->GetRMS(), fHprojy->GetMean() + 2*fHprojy->GetRMS());

      TCanvas *cv = new TCanvas("cv","cv", 600, 600);
      
      fHprojy->Draw();
      cv->SaveAs(Form("plots/fHprojy_%d.pdf",i));
      
      res[i] = fgaus->GetParameter(2);
      charge[i] = DeltaT_vs_Charge_Corrected->GetXaxis()->GetBinLowEdge(i);
      errorX[i] = 0.0;
      errorY[i] = fgaus->GetParError(2);
   }
  
  TGraphErrors* gr = new TGraphErrors( nBinsX, charge,res, errorX, errorY );
  gr->Draw("AC*"); 

  DeltaT_vs_Charge->SaveAs("test.root");
  DeltaT_vs_Charge_Corrected->SaveAs("test3.root");
  gr->SaveAs("test45.root");
}

int main(int argc, char **argv)
{
  
  MakeAmplitudePlot(argv[1], argv[2], atof(argv[3]));
}
