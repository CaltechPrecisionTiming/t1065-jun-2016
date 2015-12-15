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


// void MakePlotAlternativeFormat(string filename, string plotname, double scalefactor, double fitmin, double fitmax) {
//   // Get the tree


//   TFile *inputfile = new TFile(filename.c_str(),"READ");
//   TTree *tree = (TTree*)inputfile->Get("tree");

//   // get the variables from the ntuple
//   float t1 = 0;
//   float t2 = 0;
//   float t3 = 0;
//   float t4 = 0;
//   float ch1Amp = 0;
//   float ch2Amp = 0;
//   float ch3Amp = 0;
//   float ch4Amp = 0;

//   tree->SetBranchAddress("ch1Time",&t1);
//   tree->SetBranchAddress("ch2Time",&t2);
//   tree->SetBranchAddress("ch3Time",&t3);
//   tree->SetBranchAddress("ch4Time",&t4);
//   tree->SetBranchAddress("ch1Amp",&ch1Amp);
//   tree->SetBranchAddress("ch2Amp",&ch2Amp);
//   tree->SetBranchAddress("ch3Amp",&ch3Amp);
//   tree->SetBranchAddress("ch4Amp",&ch4Amp);

//   //create histograms
//   TH1F *dt;
//   TH1F *histAmplitude;
//   histAmplitude = new TH1F("histAmplitude","; Amplitude [V];Number of Events",50,0,0.5);
//   dt = new TH1F("dt","; #Delta t [ns]; Number of Events", 400, -1,1);
  
//   //read all entries and fill the histograms
//   Long64_t nentries = tree->GetEntries();
//   std::cout<<"Number of events in Physics Sample: "<<nentries<<std::endl;  
//   for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
//       tree->GetEntry(iEntry);    
      
//       //require cherenkov and signal in the front MCP
//       if (ch1Amp > 20 && ch4Amp > 100) {
// 	histAmplitude->Fill(ch2Amp/1000);
// 	if (ch2Amp > 20 
// 	    && ch1Amp < 490 && ch2Amp < 490
// 	    ) {
// 	  dt->Fill(t2 - t1);
// 	}
//       }
//   }

//   TCanvas * c = 0;


//   //Energy plot
//   c = new TCanvas("c","c",600,600);
  
//   histAmplitude->SetAxisRange(0.0,0.5,"X");
//   histAmplitude->SetTitle("");
//   histAmplitude->GetXaxis()->SetTitle("Pulse Height [V]");
//   histAmplitude->GetYaxis()->SetTitle("Number of Events");
//   histAmplitude->GetYaxis()->SetTitleOffset(1.3);
//   histAmplitude->SetMaximum(1.2*histAmplitude->GetMaximum());
//   histAmplitude->Draw();
//   histAmplitude->SetStats(0);
//   histAmplitude->Fit("gaus","","",fitmin,fitmax);
//   TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  
//   TLatex *tex = new TLatex();
//   tex->SetNDC();
//   tex->SetTextSize(0.040);
//   tex->SetTextFont(42);
//   tex->SetTextColor(kBlack);
//   //tex->DrawLatex(0.55, 0.80, Form("#sigma/#mu = %.0f %s",100*fitter->GetParameter(2)/fitter->GetParameter(1),"%"));
//   tex->DrawLatex(0.55, 0.85, Form("Mean = %.2f #pm %.2f %s",fitter->GetParameter(1),TMath::Max(0.01,fitter->GetParError(1)),"V"));
//   tex->DrawLatex(0.15, 0.92, Form("Attenuation Factor = %.3f",scalefactor));
  
//   c->SaveAs( Form("%s_amplitude.gif", plotname.c_str()) );
//   c->SaveAs( Form("%s_amplitude.pdf", plotname.c_str()) );


//   //time resolution plot
//   c = new TCanvas("c","c",600,600);
  
//   dt->SetAxisRange(0.1,0.4,"X");
//   dt->SetTitle("");
//   dt->GetXaxis()->SetTitle("#Delta t [ns]");
//   dt->GetXaxis()->SetTitleSize(0.045);
//   dt->GetXaxis()->SetLabelSize(0.040);
//   dt->GetXaxis()->SetLabelOffset(0.012);
//   dt->GetXaxis()->SetTitleOffset(1.00);
//   dt->GetYaxis()->SetTitle("Number of Events");
//   dt->GetYaxis()->SetTitleOffset(1.02);
//   dt->GetYaxis()->SetTitleSize(0.045);
//   dt->GetYaxis()->SetLabelSize(0.040);
//   dt->SetMaximum(1.2*dt->GetMaximum());
//   dt->Draw();
//   dt->SetStats(0);
//   dt->Fit("gaus","","",0.2,0.3);
//   TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  
//   TLatex *tex = new TLatex();
//   tex->SetNDC();
//   tex->SetTextSize(0.050);
//   tex->SetTextFont(42);
//   tex->SetTextColor(kBlack);
//   tex->DrawLatex(0.50, 0.80, Form("#sigma = %.1f #pm %.1f ps",1000*fitter->GetParameter(2),1000*fitter->GetParError(2)));
//   //tex->DrawLatex(0.55, 0.85, Form("Mean = %.2f #pm %.2f %s",fitter->GetParameter(1),fitter->GetParError(1),"V"));
//   tex->DrawLatex(0.50, 0.85, Form("Mean = %.2f ns",fitter->GetParameter(1)));
  
//   c->SaveAs( Form("%s_dt.gif", plotname.c_str()) );
//   c->SaveAs( Form("%s_dt.pdf", plotname.c_str()) );


  
// }


// Double_t TOFResolutionFunction(Double_t *x,Double_t *par)
// {
//   Double_t arg = 0;
// //if (par[2] != 0) arg = (x[0] - par[1])/par[2];
//   arg = x[0];
//   Double_t fitval = par[0]*(1.0 / sqrt(arg)) + par[1];
//   return fitval;
// }


void MakeAmplitudeVsBeamEnergyGraph() {

  //use beam energy for xaxis
  const int nPoints = 4;
  float x[nPoints] = { 4.0, 8.0, 16.0, 32.0 };
  float xerr[nPoints] = { 0.027*4.0, 0.023*8.0, 0.045*16.0, 0.05*32.0 };
  float y_charge[nPoints] = {  13.7, 33.6, 77.1, 176 }; 
  float yerr_charge[nPoints] = { 0.3, 0.3, 0.5, 3 };
  float y_MIP[nPoints] = {  0.0, 0.0, 0.0, 0.0 }; 
  float yerr_MIP[nPoints] = { 0.0, 0.0, 0.0, 0.0 };

  double chargePerMIP = 1.3;
  for (int i=0; i<nPoints; ++i) {
    y_MIP[i] = y_charge[i] / chargePerMIP;
    yerr_MIP[i] =  yerr_charge[i] / chargePerMIP;
  }

  TGraphErrors *graphChargeVsEnergy = new TGraphErrors(nPoints,x,y_charge,xerr,yerr_charge);
  graphChargeVsEnergy->SetLineWidth(3);
  TGraphErrors *graphMIPVsEnergy = new TGraphErrors(nPoints,x,y_MIP,xerr,yerr_MIP);
  graphMIPVsEnergy->SetLineWidth(3);

  TCanvas *c = 0;
  TVirtualFitter *fitter = 0;

  c = new TCanvas("c","c",800,600);
  graphChargeVsEnergy->Draw("AP");
  graphChargeVsEnergy->SetTitle("");
  graphChargeVsEnergy->GetXaxis()->SetTitle("Electron Beam Energy [GeV/c^{2}]");
  graphChargeVsEnergy->GetXaxis()->SetTitleSize(0.045);
  graphChargeVsEnergy->GetXaxis()->SetLabelSize(0.045);
  graphChargeVsEnergy->GetXaxis()->SetTitleOffset(1.0);
  graphChargeVsEnergy->GetYaxis()->SetTitle("Integrated Charge [pC]");
  graphChargeVsEnergy->GetYaxis()->SetTitleOffset(1.02);
  graphChargeVsEnergy->GetYaxis()->SetTitleSize(0.045);
  graphChargeVsEnergy->GetYaxis()->SetLabelSize(0.045);
  graphChargeVsEnergy->GetXaxis()->SetRangeUser(0,40);
  graphChargeVsEnergy->GetYaxis()->SetRangeUser(0,200);

  graphChargeVsEnergy->Fit("pol1","","");
  fitter = TVirtualFitter::GetFitter();
  
  c->SaveAs( "ChargeVsEnergyAt6X0.gif" );
  c->SaveAs( "ChargeVsEnergyAt6X0.pdf" );


  c = new TCanvas("c","c",800,600);
  graphMIPVsEnergy->Draw("AP");
  graphMIPVsEnergy->SetTitle("");
  graphMIPVsEnergy->GetXaxis()->SetTitle("Electron Beam Energy [GeV/c^{2}]");
  graphMIPVsEnergy->GetXaxis()->SetTitleSize(0.045);
  graphMIPVsEnergy->GetXaxis()->SetLabelSize(0.045);
  graphMIPVsEnergy->GetXaxis()->SetTitleOffset(1.0);
  graphMIPVsEnergy->GetYaxis()->SetTitle("Charge / Charge per MIP");
  graphMIPVsEnergy->GetYaxis()->SetTitleOffset(1.02);
  graphMIPVsEnergy->GetYaxis()->SetTitleSize(0.045);
  graphMIPVsEnergy->GetYaxis()->SetLabelSize(0.045);
  graphMIPVsEnergy->GetXaxis()->SetRangeUser(0,40);
  graphMIPVsEnergy->GetYaxis()->SetRangeUser(0,200);

  graphMIPVsEnergy->Fit("pol1","","");
  fitter = TVirtualFitter::GetFitter();
  
  c->SaveAs( "MIPVsEnergyAt6X0.gif" );
  c->SaveAs( "MIPVsEnergyAt6X0.pdf" );


}

void MakeAmplitudeVsAbsorberThicknessGraph() {

  //use beam energy for xaxis
  const int nPoints = 5;
  float x[nPoints] = { 0.0, 1.0, 2.0, 4.0, 6.0 };
  float xerr[nPoints] = { 0.1, 0.25, 0.25, 0.25, 0.25 };
  float y_charge[nPoints] = {  1.4, 5.5, 14.0, 34.4, 32.5 }; 
  float yerr_charge[nPoints] = { 0.1, 1.0, 0.3, 0.4, 0.4 };
  float y_MIP[nPoints] = {  0.0, 0.0, 0.0, 0.0, 0.0 }; 
  float yerr_MIP[nPoints] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

  double chargePerMIP = 1.3;
  for (int i=0; i<nPoints; ++i) {
    y_MIP[i] = y_charge[i] / chargePerMIP;
    yerr_MIP[i] =  yerr_charge[i] / chargePerMIP;
  }

  TGraphErrors *graphChargeVsAbsorber = new TGraphErrors(nPoints,x,y_charge,xerr,yerr_charge);
  graphChargeVsAbsorber->SetLineWidth(3);
  TGraphErrors *graphMIPVsAbsorber = new TGraphErrors(nPoints,x,y_MIP,xerr,yerr_MIP);
  graphMIPVsAbsorber->SetLineWidth(3);

  TCanvas *c = 0;
  TVirtualFitter *fitter = 0;

  c = new TCanvas("c","c",800,600);
  graphChargeVsAbsorber->Draw("AP");
  graphChargeVsAbsorber->SetTitle("");
  graphChargeVsAbsorber->GetXaxis()->SetTitle("Tungsten Absorber Thickness [X_{0}]");
  graphChargeVsAbsorber->GetXaxis()->SetTitleSize(0.045);
  graphChargeVsAbsorber->GetXaxis()->SetLabelSize(0.045);
  graphChargeVsAbsorber->GetXaxis()->SetTitleOffset(1.0);
  graphChargeVsAbsorber->GetYaxis()->SetTitle("Integrated Charge [pC]");
  graphChargeVsAbsorber->GetYaxis()->SetTitleOffset(1.02);
  graphChargeVsAbsorber->GetYaxis()->SetTitleSize(0.045);
  graphChargeVsAbsorber->GetYaxis()->SetLabelSize(0.045);
  graphChargeVsAbsorber->GetXaxis()->SetRangeUser(0.0,8.0);
  graphChargeVsAbsorber->GetYaxis()->SetRangeUser(0,40);

  // graphChargeVsAbsorber->Fit("pol2","","");
  // fitter = TVirtualFitter::GetFitter();
  
  c->SaveAs( "ChargeVsAbsorberAt8GeV.gif" );
  c->SaveAs( "ChargeVsAbsorberAt8GeV.pdf" );


  c = new TCanvas("c","c",800,600);
  graphMIPVsAbsorber->Draw("AP");
  graphMIPVsAbsorber->SetTitle("");
  graphMIPVsAbsorber->GetXaxis()->SetTitle("Tungsten Absorber Thickness [X_{0}]");
  graphMIPVsAbsorber->GetXaxis()->SetTitleSize(0.045);
  graphMIPVsAbsorber->GetXaxis()->SetLabelSize(0.045);
  graphMIPVsAbsorber->GetXaxis()->SetTitleOffset(1.0);
  graphMIPVsAbsorber->GetYaxis()->SetTitle("Charge / Charge per MIP");
  graphMIPVsAbsorber->GetYaxis()->SetTitleOffset(1.02);
  graphMIPVsAbsorber->GetYaxis()->SetTitleSize(0.045);
  graphMIPVsAbsorber->GetYaxis()->SetLabelSize(0.045);
  graphMIPVsAbsorber->GetXaxis()->SetRangeUser(0.0, 8.0);
  graphMIPVsAbsorber->GetYaxis()->SetRangeUser(0, 35);

  // graphMIPVsAbsorber->Fit("pol1","","");
  // fitter = TVirtualFitter::GetFitter();
  
  c->SaveAs( "MIPVsAbsorberAt8GeV.gif" );
  c->SaveAs( "MIPVsAbsorberAt8GeV.pdf" );


}



void MakeAmplitudeVsBiasVoltageGraph() {

  //use beam energy for xaxis
  const int nPoints = 5;
  float x[nPoints] = { 200, 300, 400, 500, 600 };
  float xerr[nPoints] = { 0.01*200, 0.01*300, 0.01*400, 0.01*500, 0.01*600 };
  float y_charge[nPoints] = {  58.4, 68.0, 73.8, 77.1, 78.6 }; 
  float yerr_charge[nPoints] = { 1.0, 1.4, 1.3, 0.5, 1.4 };
  float y_MIP[nPoints] = {  0.0, 0.0, 0.0, 0.0, 0.0 }; 
  float yerr_MIP[nPoints] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

  double chargePerMIP = 1.3;
  for (int i=0; i<nPoints; ++i) {
    y_MIP[i] = y_charge[i] / chargePerMIP;
    yerr_MIP[i] =  yerr_charge[i] / chargePerMIP;
  }

  TGraphErrors *graphChargeVsBiasVoltage = new TGraphErrors(nPoints,x,y_charge,xerr,yerr_charge);
  graphChargeVsBiasVoltage->SetLineWidth(3);
  TGraphErrors *graphMIPVsBiasVoltage = new TGraphErrors(nPoints,x,y_MIP,xerr,yerr_MIP);
  graphMIPVsBiasVoltage->SetLineWidth(3);

  TCanvas *c = 0;
  TVirtualFitter *fitter = 0;

  c = new TCanvas("c","c",800,600);
  graphChargeVsBiasVoltage->Draw("AP");
  graphChargeVsBiasVoltage->SetTitle("");
  graphChargeVsBiasVoltage->GetXaxis()->SetTitle("Bias Voltage [V]");
  graphChargeVsBiasVoltage->GetXaxis()->SetTitleSize(0.045);
  graphChargeVsBiasVoltage->GetXaxis()->SetLabelSize(0.045);
  graphChargeVsBiasVoltage->GetXaxis()->SetTitleOffset(1.0);
  graphChargeVsBiasVoltage->GetYaxis()->SetTitle("Integrated Charge [pC]");
  graphChargeVsBiasVoltage->GetYaxis()->SetTitleOffset(1.02);
  graphChargeVsBiasVoltage->GetYaxis()->SetTitleSize(0.045);
  graphChargeVsBiasVoltage->GetYaxis()->SetLabelSize(0.045);
  graphChargeVsBiasVoltage->GetXaxis()->SetRangeUser(0.0,700);
  graphChargeVsBiasVoltage->GetYaxis()->SetRangeUser(0,90);

  // graphChargeVsBiasVoltage->Fit("pol2","","");
  // fitter = TVirtualFitter::GetFitter();
  
  c->SaveAs( "ChargeVsBiasVoltageAt6X0_16GeV.gif" );
  c->SaveAs( "ChargeVsBiasVoltageAt6X0_16GeV.pdf" );


  c = new TCanvas("c","c",800,600);
  graphMIPVsBiasVoltage->Draw("AP");
  graphMIPVsBiasVoltage->SetTitle("");
  graphMIPVsBiasVoltage->GetXaxis()->SetTitle("Bias Voltage [V]");
  graphMIPVsBiasVoltage->GetXaxis()->SetTitleSize(0.045);
  graphMIPVsBiasVoltage->GetXaxis()->SetLabelSize(0.045);
  graphMIPVsBiasVoltage->GetXaxis()->SetTitleOffset(1.0);
  graphMIPVsBiasVoltage->GetYaxis()->SetTitle("Charge / Charge per MIP");
  graphMIPVsBiasVoltage->GetYaxis()->SetTitleOffset(1.02);
  graphMIPVsBiasVoltage->GetYaxis()->SetTitleSize(0.045);
  graphMIPVsBiasVoltage->GetYaxis()->SetLabelSize(0.045);
  graphMIPVsBiasVoltage->GetXaxis()->SetRangeUser(0.0, 700);
  graphMIPVsBiasVoltage->GetYaxis()->SetRangeUser(0, 70);

  // graphMIPVsBiasVoltage->Fit("pol1","","");
  // fitter = TVirtualFitter::GetFitter();
  
  c->SaveAs( "MIPVsBiasVoltageAt6X0_16GeV.gif" );
  c->SaveAs( "MIPVsBiasVoltageAt6X0_16GeV.pdf" );


}



void makeEMShowerEnergyPlots() {

  //*************************************
  // Charge Vs Energy
  //*************************************

  //Noise Runs
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run120-0.dat-full.root", "NoiseNoBeam", -99, sqrt(10.0), -10, 10, -2,2 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run121-0.dat-full.root", "NoiseNoBeamNoAttenuator", -99, 1, -10, 10, -2,2 );


  //Protons
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run125.root", "Proton", -99, 1, -2, 8, 0,3 );


  //Electrons
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run123-124.root", "Electron_0X0", -99, 1, -2, 8, 0,3 );
  
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run112-500V_123.root", "Electron_6X0_4GeV", 0.02, 10.0, 0, 80, 5, 25 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run109.root","Electron_6X0_8GeV", 0.05, 10.0, 0, 300, 10, 60);
  //MakeChargePlot("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run110.root","Electron_6X0_16GeV", 0.1, 10.0, 0, 300, 30, 140);
  // MakeChargePlot("/afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run113-114.root","Electron_6X0_32GeV", 0.1, 10.0, 0, 500, 100, 300);
 
  //MakeAmplitudeVsBeamEnergyGraph();
  

  //*************************************
  // Charge Vs Absorber Thickness
  //*************************************
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run122.root", "Electron_1X0_8GeV", 0.05, 1.0, 0, 40, 2, 12 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run119.root", "Electron_1X0_8GeV_10dbAttenuator", 0.02, sqrt(10.0), 0, 40, 2, 12 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run118.root", "Electron_2X0_8GeV", 0.03, sqrt(10.0), 0, 80, 4, 30 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run117.root", "Electron_4X0_8GeV", 0.05, sqrt(10.0), 0, 150, 10, 60 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run116.root", "Electron_6X0_8GeV", 0.05, sqrt(10.0), 0, 150, 10, 60 );

  //MakeAmplitudeVsAbsorberThicknessGraph();


  //*************************************
  // Bias Voltage Scan
  //*************************************
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run111_600V.root", "Electron_6X0_8GeV_600V", 0.1, 10.0, 0, 300, 30, 150 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run110.root", "Electron_6X0_8GeV_500V", 0.1, 10.0, 0, 300, 30, 150 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run111_400V.root", "Electron_6X0_8GeV_400V", 0.1, 10.0, 0, 300, 30, 140 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run111_300V.root", "Electron_6X0_8GeV_300V", 0.1, 10.0, 0, 300, 30, 130 );
  //MakeChargePlot("root://eoscms:///eos/cms/store/group/phys_susy/razor/Timing/t1065-dec-2015/reco/v1/t1065_dec2015_run111_200V.root", "Electron_6X0_8GeV_200V", 0.1, 10.0, 0, 300, 20, 100 );

  MakeAmplitudeVsBiasVoltageGraph();

}
