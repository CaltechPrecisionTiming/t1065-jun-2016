#ifndef Aux_HH
#define Aux_HH

#include <TGraphErrors.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>


TGraphErrors* GetTGraph( short* channel, float* time );
double GetGaussTime( TGraphErrors* pulse );
int FindMin( int n, short *a);
float GausFit_MeanTime(TGraphErrors * pulse, const float index_first, const float index_last);
float RisingEdgeFitTime(TGraphErrors * pulse, const float index_min, TString fname);
float GausFit_MeanTime(TGraphErrors* pulse, const float index_first, const float index_last, TString fname);
float GetBaseline( int peak, short *a );

TGraphErrors* GetTGraphFilter( short* channel, float* time, TString pulseName );

#endif
