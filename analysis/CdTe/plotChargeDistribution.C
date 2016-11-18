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




void makeChargeDistribution(string filename, string plotname, double ampCutOnPhotek, 
			    double beamXMin, double beamXMax, double beamYMin, double beamYMax,
			    double xmin, double xmax, double fitmin, double fitmax) {

  // double xmin = 0;
  // double xmax = 100*3.16228*(1.0/63.0957);
  // double fitmin = 40*3.16228*(1.0/63.0957);
  // double fitmax = 80*3.16228*(1.0/63.0957);

  TFile *inputfile = TFile::Open(filename.c_str(),"READ");
  
  TTree *tree = (TTree*)inputfile->Get("t1065");

  // get the variables from the ntuple
  float amp[36];
  float integral[36];
  float gauspeak[36];
  float linearTime30[36];
  float beamX;
  float beamY;

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("gauspeak",1);
  tree->SetBranchStatus("amp",1);
  tree->SetBranchStatus("intfull",1);
  tree->SetBranchStatus("linearTime30",1);
  tree->SetBranchStatus("TDCx",1);
  tree->SetBranchStatus("TDCy",1);
  tree->SetBranchAddress("gauspeak",gauspeak);
  tree->SetBranchAddress("amp",amp);
  tree->SetBranchAddress("intfull",integral);
  tree->SetBranchAddress("linearTime30",linearTime30);
  tree->SetBranchAddress("TDCx",&beamX);
  tree->SetBranchAddress("TDCy",&beamY);

  //create histograms
  TH1F *histIntCharge;
  histIntCharge = new TH1F("histIntCharge","; Integrated Charge [fC];Number of Events", 25, xmin, xmax);

  
  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();

  std::cout<<"Number of events in Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
    if (iEntry %1000 == 0) 
      cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
    // cout << "here1\n";
    float photekTimeGauss = gauspeak[0];
    float CdTeTime = linearTime30[1];
    float photekAmp = amp[0]*(10.0);
    float CdTeAmp = amp[1]*3.16228*(1.0/63.0957);
    float photekCharge = integral[0]*(10.0);
    float CdTeCharge = integral[1]*3.16228*(1.0/63.0957);
    // cout << "here2\n";
       
    //use photek amplitude cut for electron ID
    //cout << "test: " << photekAmp << " " << siliconIntegral << "\n";
    if( !(photekAmp > 0.5)) continue;
    if(!(beamX > 4.5 && beamX < 12.5)) continue;
    if(!(beamY > -6.5 && beamY < 1.5)) continue;
    // cout << "here3\n";

    //don't fill overflow bins
    //if (1000* siliconIntegral * attenuationFactor / amplificationFactor > xmax) continue;
    
    //histIntCharge->Fill( CdTeCharge );
    histIntCharge->Fill( CdTeAmp );

    cout << CdTeCharge << " " << beamX << " " << beamY << " " << CdTeAmp << "\n";

    //cout << 1000* amp[21] << " : " << amplificationFactor << " : " << siliconIntegral * attenuationFactor / amplificationFactor << "\n";
 
  }


  TCanvas * c = 0;


  //Energy plot
  c = new TCanvas("c","c",600,600);  
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.17);
  histIntCharge->SetAxisRange(xmin,xmax,"X");
  histIntCharge->SetTitle("");
  histIntCharge->GetXaxis()->SetTitle("Integrated Charge [fC]");
  histIntCharge->GetXaxis()->SetTitleSize(0.045);
  histIntCharge->GetXaxis()->SetLabelSize(0.045);
  histIntCharge->GetYaxis()->SetTitle("Number of Events");
  histIntCharge->GetYaxis()->SetTitleOffset(1.3);
  histIntCharge->GetYaxis()->SetTitleSize(0.05);
  histIntCharge->GetYaxis()->SetLabelSize(0.045);
  histIntCharge->GetYaxis()->SetLabelOffset(0.015);
  histIntCharge->GetYaxis()->SetTitleOffset(1.7);
  histIntCharge->SetMaximum(1.2*histIntCharge->GetMaximum());
  histIntCharge->Draw();
  histIntCharge->SetStats(0);
  histIntCharge->Fit("gaus","","",fitmin,fitmax);
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.050);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  /* tex->DrawLatex(0.45, 0.85, Form("Mean = %.1f #pm %.1f %s",fitter->GetParameter(1),TMath::Max(0.01,fitter->GetParError(1)),"fC")); */
  /* tex->DrawLatex(0.45, 0.80, Form("#sigma = %.1f #pm %.1f %s",fitter->GetParameter(2),TMath::Max(0.01,fitter->GetParError(2)),"fC")); */
  tex->DrawLatex(0.45, 0.85, Form("Mean = %.3f %s",fitter->GetParameter(1),"fC"));
  tex->DrawLatex(0.45, 0.80, Form("#sigma = %.3f %s",fitter->GetParameter(2),"fC"));
  //tex->DrawLatex(0.15, 0.92, Form("Attenuation Factor = %.3f",attenuationFactor));
  
  c->SaveAs( Form("%s_charge.gif", plotname.c_str()) );
  c->SaveAs( Form("%s_charge.pdf", plotname.c_str()) );
 

}

void plotChargeDistribution() {

  // makeChargeDistribution( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v5/analysis_5568.root", 
  // 			  "100GeV", 0.5,
  //                      4.5,12.5,-6.5,1.5,
  // 			  //0, 10, 2, 4
  // 			  0, 0.05, 0.015, 0.028
  // 			  );

   // makeChargeDistribution( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v5/analysis_5570.root", 
   // 			   "200GeV", 0.5,
   // 			   2.5,13.5,-8,2,
   // 			   //0, 10, 3,6.5
   // 			   0, 0.08, 0.025,0.045
   // 			   );

  makeChargeDistribution( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Nov2016CERN/ntuples_v5/analysis_5571-5576.root", 
  			  "50GeV", 0.5,
			  2.5,13.5,-8.0,2.0,
  			  //0, 5,1,2.5
			  0,0.03, 0.007, 0.017
  			  );
 

}
