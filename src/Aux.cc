#include "Aux.hh"

TGraphErrors* GetTGraph(  short* channel, float* time )
{		
  //Setting Errors
  float errorX[1024], errorY[1024], channelFloat[1024];
  float _errorY = 0.00; //5%error on Y
  for ( int i = 0; i < 1024; i++ )
    {
      errorX[i]       = .0;
      errorY[i]       = _errorY*channel[i];
      channelFloat[i] = -channel[i];
    }
  TGraphErrors* tg = new TGraphErrors( 1024, time, channelFloat, errorX, errorY );
  return tg;
};

////////////////////////////////////////////
// find minimum of the pulse
// aa added protection against pulses with single high bin
////////////////////////////////////////////
int FindMin( int n, short *a) {
  
  if (n <= 0 || !a) return -1;
  float xmin = a[5];
  int loc = 0;
  for  (int i = 5; i < n-5; i++) {
    if (xmin > a[i] && a[i+1] < 0.5*a[i])  {
      xmin = a[i];
      loc = i;
    }
  }
  
  return loc;
}

// find the mean time from gaus fit
float GausFit_MeanTime(TGraphErrors* pulse, const float index_first, const float index_last)
{
  TF1* fpeak = new TF1("fpeak","gaus", index_first, index_last);
  pulse->Fit("fpeak","Q","", index_first, index_last);
  
  float timepeak = fpeak->GetParameter(1);
  delete fpeak;
  
  return timepeak;
}

float GausFit_MeanTime(TGraphErrors* pulse, const float index_first, const float index_last, TString fname)
{
  TF1* fpeak = new TF1("fpeak","gaus", index_first, index_last);
  pulse->Fit("fpeak","Q","", index_first, index_last);
  
  TCanvas* c = new TCanvas("canvas","canvas",800,400) ;
  float timepeak = fpeak->GetParameter(1);
  pulse->GetXaxis()->SetLimits( timepeak-10, timepeak+10);
  pulse->SetMarkerSize(1);
  pulse->SetMarkerStyle(20);
  pulse->Draw("AP");
  c->SaveAs(fname+".pdf");
  delete fpeak;
  
  return timepeak;
}

double GetGaussTime( TGraphErrors* pulse )
{
  return 0;
};


float GetBaseline( int peak, short *a ) {

  float tmpsum = 0;
  float tmpcount = 0;
  //std::cout << "GGG\n";
  if (peak < 300) {
    for  (int i = peak + 200; i < 1000; i++) {
      // std::cout << i << " : " << a[i] << "\n";
      tmpsum += a[i];
      tmpcount += 1.0;
    }
  } else {
    for  (int i = 5; i < peak-200; i++) {
      // std::cout << i << " : " << a[i] << "\n";
      tmpsum += a[i];
      tmpcount += 1.0;
    }
  }
  // std::cout << tmpsum / tmpcount << "\n";

  return tmpsum / tmpcount;
}


float GetPulseIntegral(int peak, short *a) 
{
  float integral = 0.;

  //for sharp gaussian peak type pulses
  for (int i=245; i < 295; i++) {
    integral += a[i] * 0.2 * 1e-9 * (1.0/4096.0) * (1.0/50.0) * 1e12; //in units of fC, for 50Ohm termination
  }

  return -1.0 * integral;
}
