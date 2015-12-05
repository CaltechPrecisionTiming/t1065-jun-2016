#include "Aux.hh"

TGraphErrors* GetTGraph(  short* channel, float* time )
{		
  //Setting Errors
  float errorX[1024], errorY[1024], channelFloat[1024];
  float _errorY = 0.0; //5%error on Y
  for ( int i = 0; i < 1024; i++ )
    {
      errorX[i]       = .0;
      errorY[i]       = _errorY*channel[i];
      channelFloat[i] = -channel[i];
    }
  TGraphErrors* tg = new TGraphErrors( 1024, time, channelFloat, errorX, errorY );
  return tg;
};
double GetGaussTime( TGraphErrors* pulse )
{
  return 0;
};
