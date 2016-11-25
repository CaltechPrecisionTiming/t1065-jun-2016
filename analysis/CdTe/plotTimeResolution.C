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



void MakeTimeResolutionVsBeamEnergyGraph() {

  //use beam energy for xaxis
  const int nPoints_noCorr = 7;
  float x_noCorr[nPoints_noCorr] = { 2, 3.5, 5, 7, 50.0 , 100.0, 200.0 };
  float xerr_noCorr[nPoints_noCorr] = { 0, 0, 0, 0, 0, 0, 0 };
  float y_charge_noCorr[nPoints_noCorr] = { 113 , 99 , 92 , 86 , 73 , 43, 43 }; 
  float yerr_charge_noCorr[nPoints_noCorr] = { 13, 7, 5 ,  7 , 8, 4, 3 };


  //We decide not to use these because the cuts made on the beam position
  //are not consistent between the different energy points
  // const int nPoints_Corr = 3;
  // float x_Corr[nPoints_Corr] = { 50.0 , 100.0, 200.0 };
  // float xerr_Corr[nPoints_Corr] = { 0, 0, 0 };
  // float y_charge_Corr[nPoints_Corr] = { 36 , 23, 21 }; 
  // float yerr_charge_Corr[nPoints_Corr] = { 4, 4, 3 };



  TGraphErrors *graphChargeVsEnergyAt6X0_Resolution = new TGraphErrors(nPoints_noCorr,x_noCorr,y_charge_noCorr,xerr_noCorr,yerr_charge_noCorr);
  graphChargeVsEnergyAt6X0_Resolution->SetMarkerStyle(0);
  graphChargeVsEnergyAt6X0_Resolution->SetLineColor(kBlue+1);
  graphChargeVsEnergyAt6X0_Resolution->SetLineWidth(2);


  TCanvas *c = 0;
  TVirtualFitter *fitter = 0;
  TLatex *tex = 0;

  c = new TCanvas("c","c",800,800);
  c->SetLeftMargin(0.18);
  c->SetBottomMargin(0.15);
  c->SetLogx();

  graphChargeVsEnergyAt6X0_Resolution->Draw("AP");
  graphChargeVsEnergyAt6X0_Resolution->SetTitle("");
  graphChargeVsEnergyAt6X0_Resolution->GetXaxis()->SetTitle("Electron Beam Energy [GeV]");
  graphChargeVsEnergyAt6X0_Resolution->GetXaxis()->SetTitleSize(0.045);
  graphChargeVsEnergyAt6X0_Resolution->GetXaxis()->SetLabelSize(0.045);
  graphChargeVsEnergyAt6X0_Resolution->GetXaxis()->SetTitleOffset(1.2);
  graphChargeVsEnergyAt6X0_Resolution->GetYaxis()->SetTitle("Time Resolution [ps]");
  graphChargeVsEnergyAt6X0_Resolution->GetYaxis()->SetTitleOffset(1.4);
  graphChargeVsEnergyAt6X0_Resolution->GetYaxis()->SetTitleSize(0.05);
  graphChargeVsEnergyAt6X0_Resolution->GetYaxis()->SetLabelSize(0.045);
  //graphChargeVsEnergyAt6X0_Resolution->GetXaxis()->SetRangeUser(0,40);
  graphChargeVsEnergyAt6X0_Resolution->SetMarkerStyle(20);
  graphChargeVsEnergyAt6X0_Resolution->SetMarkerSize(1);
  graphChargeVsEnergyAt6X0_Resolution->GetYaxis()->SetRangeUser(0,150);

  // graphChargeVsEnergyAt6X0_Resolution->Fit("pol1","","");
  // fitter = TVirtualFitter::GetFitter();
  // c->SetLeftMargin(0.15);
  // c->SetBottomMargin(0.12);

  // tex = new TLatex();
  // tex->SetNDC();
  // tex->SetTextSize(0.050);
  // tex->SetTextFont(42);
  // tex->SetTextColor(kBlack);
  // tex->DrawLatex(0.35, 0.93, "Absorber : 6 X_{0} Tungsten");


  c->SaveAs( "TimeResolutionVsEnergy.gif" );
  c->SaveAs( "TimeResolutionVsEnergy.pdf" );

  return;


  //use beam energy for xaxis
  const int nPoints_T9 = 4;
  float x_T9[nPoints_T9] = { 2, 3.5, 5, 7 };
  float xerr_T9[nPoints_T9] = { 0, 0, 0, 0 };
  float y_charge_T9[nPoints_T9] = {  0.16, 0.26, 0.26, 0.40 }; 
  float yerr_charge_T9[nPoints_T9] = { 0.008, 0.015, 0.03, 0.026 };
  float yResolution_charge_T9[nPoints_T9] = { 0.08, 0.15, 0.17, 0.20 };


  TGraphErrors *graphChargeVsEnergyAt2X0_Resolution = new TGraphErrors(nPoints_T9,x_T9,y_charge_T9,xerr_T9,yResolution_charge_T9);
  TGraphErrors *graphChargeVsEnergyAt2X0 = new TGraphErrors(nPoints_T9,x_T9,y_charge_T9,xerr_T9,yerr_charge_T9);
  graphChargeVsEnergyAt2X0_Resolution->SetMarkerStyle(20);
  graphChargeVsEnergyAt2X0_Resolution->SetLineWidth(3);
  graphChargeVsEnergyAt2X0_Resolution->SetLineColor(kGreen+1);
  graphChargeVsEnergyAt2X0->SetLineColor(kBlue);
  graphChargeVsEnergyAt2X0->SetLineWidth(2);

  c = new TCanvas("c","c",800,800);
  graphChargeVsEnergyAt2X0_Resolution->Draw("AP");
  graphChargeVsEnergyAt2X0_Resolution->SetTitle("");
  graphChargeVsEnergyAt2X0_Resolution->GetXaxis()->SetTitle("Electron Beam Energy [GeV/c^{2}]");
  graphChargeVsEnergyAt2X0_Resolution->GetXaxis()->SetTitleSize(0.045);
  graphChargeVsEnergyAt2X0_Resolution->GetXaxis()->SetLabelSize(0.045);
  graphChargeVsEnergyAt2X0_Resolution->GetXaxis()->SetTitleOffset(1.05);
  graphChargeVsEnergyAt2X0_Resolution->GetYaxis()->SetTitle("Integrated Charge [pC]");
  graphChargeVsEnergyAt2X0_Resolution->GetYaxis()->SetTitleOffset(1.2);
  graphChargeVsEnergyAt2X0_Resolution->GetYaxis()->SetTitleSize(0.05);
  graphChargeVsEnergyAt2X0_Resolution->GetYaxis()->SetLabelSize(0.045);
  //graphChargeVsEnergyAt2X0_Resolution->GetXaxis()->SetRangeUser(0,40);
  graphChargeVsEnergyAt2X0_Resolution->GetYaxis()->SetRangeUser(0,0.75);

  graphChargeVsEnergyAt2X0->Draw("PE1same");
  graphChargeVsEnergyAt2X0->Fit("pol1","","");
  fitter = TVirtualFitter::GetFitter();

  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.12);

  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.050);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.35, 0.93, "Absorber : 2 X_{0} Lead");

  c->SaveAs( "ChargeVsEnergyAt2X0.gif" );
  c->SaveAs( "ChargeVsEnergyAt2X0.pdf" );


}



void makeTimeResolutionDistributionT9(string filename, string plotname, string plotTitle,
				      double ampCutOnMCP, double ampCutOnLYSO, double ampCutOnTrigger,
				      double beamXMin, double beamXMax, double beamYMin, double beamYMax,
				      int nbins, double xmin, double xmax, double fitmin, double fitmax) {
  

  TFile *inputfile = TFile::Open(filename.c_str(),"READ");
  
  TTree *tree = (TTree*)inputfile->Get("TBtime");

  // get the variables from the ntuple
  float amp_ref;
  float amp_CdTe;
  float amp_LYSO;
  float amp_Trigger;
  float time_ref;
  float time_CdTe;
  float beamX;
  float beamY;

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("time_second",1);
  tree->SetBranchStatus("time_ref2",1);
  tree->SetBranchStatus("max_ref2",1);
  tree->SetBranchStatus("max_second",1);
  tree->SetBranchStatus("max_third",1);
  tree->SetBranchStatus("max_fourth",1);
  tree->SetBranchStatus("TDCx",1);
  tree->SetBranchStatus("TDCy",1);
  tree->SetBranchAddress("time_second",&time_CdTe);
  tree->SetBranchAddress("time_ref2",&time_ref);
  tree->SetBranchAddress("max_ref2",&amp_ref);
  tree->SetBranchAddress("max_second",&amp_CdTe);
  tree->SetBranchAddress("max_third",&amp_LYSO);
  tree->SetBranchAddress("max_fourth",&amp_Trigger);
  tree->SetBranchAddress("TDCx",&beamX);
  tree->SetBranchAddress("TDCy",&beamY);

  //create histograms
  TH1F *histDeltaT;
  histDeltaT = new TH1F("histDeltaT","; #Delta t [ns];Number of Events", nbins, xmin, xmax);

  
  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();

  std::cout<<"Number of events in Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
    if (iEntry %1000 == 0) 
      cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);    
    // cout << "here1\n";
    float MCPTimeGauss = time_ref;
    float CdTeTime = time_CdTe;
    float MCPAmp = amp_ref;
    float LYSOAmp = amp_LYSO;
    float TriggerAmp = amp_Trigger;
    float CdTeAmp = amp_CdTe*(1.0/63.0957);
    // cout << "here2\n";
       
    //use MCP amplitude cut for electron ID
    //cout << "test: " << MCPAmp << " " << TriggerAmp << " " << LYSOAmp  << "\n";
    if( !(MCPAmp > ampCutOnMCP)) continue;
    if( !(TriggerAmp > ampCutOnTrigger)) continue;
    if( !(LYSOAmp > ampCutOnLYSO)) continue;
    if(!(beamX > beamXMin && beamX < beamXMax)) continue;
    if(!(beamY > beamYMin && beamY < beamYMax)) continue;
     // cout << "here3\n";

    //don't fill overflow bins
    //if (1000* siliconIntegral * attenuationFactor / amplificationFactor > xmax) continue;
    
    // cout << "time:" << CdTeTime << " " << MCPTimeGauss << "\n";
    histDeltaT->Fill( 0.2*(CdTeTime - MCPTimeGauss) );

    //cout << CdTeCharge << " " << beamX << " " << beamY << " " << CdTeAmp << "\n";

    //cout << 1000* amp[21] << " : " << amplificationFactor << " : " << siliconIntegral * attenuationFactor / amplificationFactor << "\n";
 
  }


  TCanvas * c = 0;


  //Energy plot
  c = new TCanvas("c","c",600,600);  
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.17);
  histDeltaT->SetAxisRange(xmin,xmax,"X");
  histDeltaT->SetTitle("");
  histDeltaT->GetXaxis()->SetTitle("#Delta t [ns]");
  histDeltaT->GetXaxis()->SetTitleSize(0.045);
  histDeltaT->GetXaxis()->SetLabelSize(0.045);
  histDeltaT->GetYaxis()->SetTitle("Number of Events");
  histDeltaT->GetYaxis()->SetTitleOffset(1.3);
  histDeltaT->GetYaxis()->SetTitleSize(0.05);
  histDeltaT->GetYaxis()->SetLabelSize(0.045);
  histDeltaT->GetYaxis()->SetLabelOffset(0.015);
  histDeltaT->GetYaxis()->SetTitleOffset(1.7);
  histDeltaT->SetMaximum(1.2*histDeltaT->GetMaximum());
  histDeltaT->Draw();
  histDeltaT->SetStats(0);
  histDeltaT->Fit("gaus","","",fitmin,fitmax);
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.050);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  //tex->DrawLatex(0.45, 0.85, Form("Mean = %.2f %s",fitter->GetParameter(1),"pC"));
  tex->DrawLatex(0.45, 0.80, Form("#sigma = %.3f #pm %.3f %s",fitter->GetParameter(2),fitter->GetParError(2),"ns"));

  tex->DrawLatex(0.18, 0.93, Form("%s", plotTitle.c_str()));

  c->SaveAs( Form("%s_charge.gif", plotname.c_str()) );
  c->SaveAs( Form("%s_charge.pdf", plotname.c_str()) );
 

}





void plotTimeResolution(double energy = -1) {

  if (energy == 2) {
    makeTimeResolutionDistributionT9( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Jun2016CERN/treeBaseFall/treeBaseFall_4627.root", 
				      "2GeV", "2 GeV Electrons, 2 X_{0} Lead Absorber", 
				      5, 700, 150,
				      -20,20, -10,10,
				      30, 3,4.5,3.5,4.0
				      );
  }
  
 if (energy == 3.5) {
    makeTimeResolutionDistributionT9( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Jun2016CERN/treeBaseFall/treeBaseFall_4626.root", 
				      "3p5GeV", "3.5 GeV Electrons, 2 X_{0} Lead Absorber", 
				      5, 700, 150,
				      -20,20, -10,10,
				      30, 3, 4.5, 3.5, 4.0
				      );
  }

  if (energy == 5) {
    makeTimeResolutionDistributionT9( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Jun2016CERN/treeBaseFall/treeBaseFall_4629.root", 
				      "5GeV", "5 GeV Electrons, 2 X_{0} Lead Absorber", 
				      5, 700, 150,
				      -20,20, -10,10,
				      30, 3, 4.5, 3.5, 3.9
				      );
  }

  if (energy == 7) {
    makeTimeResolutionDistributionT9( "/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/Jun2016CERN/treeBaseFall/treeBaseFall_4632.root", 
				      "7GeV", "7 GeV Electrons, 2 X_{0} Lead Absorber", 
				      5, 700, 150,
				      -20,20, -10,10,
				      60, 3, 4.5, 3.6, 3.9
				      );
  }

 
  if (energy == -1) {

    MakeTimeResolutionVsBeamEnergyGraph();
  }

}
