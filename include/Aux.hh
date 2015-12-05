#ifndef Aux_HH
#define Aux_HH

#include <TGraphErrors.h>
#include <TH1F.h>


TGraphErrors* GetTGraph( short* channel, float* time );
double GetGaussTime( TGraphErrors* pulse );
#endif
