#ifndef ImageFactory_hh
#define ImageFactory_hh

#include "TH2.h"
#include "TGraph.h"
#include "TDirectory.h"
#include "TApplication.h"
#include "TCanvas.h"



class ImageFactory
{

public:
  ImageFactory(){}
  ~ImageFactory(){}


private:

  const static int NBINX = 512;
  const static int NBINY = 512;
  const static int OPT = 1;
  const static int X_EXP = 9;
  const static int Y_EXP = 9;

  TH2F *orgimage;
  TH2F *fft_orgimage;
  TH2F *fft_realpart;
  TH2F *fft_impart;
  TH2F *filtered_image;

  TH2F *lowpassfilter;
  TH2F *highpassfilter;
  TH2F *bandpassfilter;

  TH2F *psf;
  TH2F *psf_fft;
  TH2F *psf_realpart;
  TH2F *psf_impart;

  TH2F *power;
  //  TH1D *powerspectrum;
  TGraph *powerspectrum;
  std::vector<float> x_reg1,y_reg1;
  //  float *x_reg1_pointer;
  //  float *y_reg1_pointer;
  bool flag;

public:

  void Init(TApplication *app);  
  void SetInitImage(TH2F *his);
  void SetPsfImage(TH2F *his);
  void Run();
  void Save(TDirectory *dir);

private:

  bool FFT2();
  bool FFT2Inv();
  bool FFT2LowPassFilter();
  bool FFT2HighPassFilter();
  bool FFT2BandPassFilter();
  bool FFT2Deconvolved();

  void FFT2core(float a_rl[NBINY][NBINX], float a_im[NBINY][NBINX],int inv);
  void FFT1core(float *a_rl, float *a_im,
		int length, int ex,
		float *sin_tbl, float *cos_tbl ,float *buf);
  void cstb(int length, int inv, float *sin_tbl, float *cos_tbl);
  void birv(float *a,int length,int ex,float *b);
  void rvmtx1(float a[NBINY][NBINX],float b[NBINX][NBINY],int xsize,int ysize);
  void rvmtx2(float a[NBINX][NBINY],float b[NBINY][NBINX],int xsize,int ysize);

  void palette_b(bool yesorno);
  void palette_rainbow(bool yesorno);
  void Print();
  TApplication *theApp;
  
  TCanvas *c1 ;
  TCanvas *c2 ;
  TCanvas *c3 ;
  TCanvas *c4 ;
  
};

#endif

