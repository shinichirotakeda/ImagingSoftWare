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

  const static int NBINX = 128;
  const static int NBINY = 128;
  const static int OPT = 1;
  const static int X_EXP = 7;
  const static int Y_EXP = 7;

  TH2D *orgimage;
  TH2D *fft_orgimage;
  TH2D *fft_realpart;
  TH2D *fft_impart;
  TH2D *filtered_image;

  TH2D *lowpassfilter;
  TH2D *highpassfilter;
  TH2D *bandpassfilter;

  TH2D *psf;
  TH2D *psf_fft;
  TH2D *psf_realpart;
  TH2D *psf_impart;

  TH2D *power;
  //  TH1D *powerspectrum;
  TGraph *powerspectrum;
  std::vector<double> x_reg1,y_reg1;
  //  double *x_reg1_pointer;
  //  double *y_reg1_pointer;
  bool flag;

public:

  void Init(TApplication *app);  
  void SetInitImage(TH2D *his);
  void SetPsfImage(TH2D *his);
  void Run();
  void Save(TDirectory *dir);

private:

  bool FFT2();
  bool FFT2Inv();
  bool FFT2LowPassFilter();
  bool FFT2HighPassFilter();
  bool FFT2BandPassFilter();
  bool FFT2Deconvolved();

  void FFT2core(double a_rl[NBINY][NBINX], double a_im[NBINY][NBINX],int inv);
  void FFT1core(double *a_rl, double *a_im,
		int length, int ex,
		double *sin_tbl, double *cos_tbl ,double *buf);
  void cstb(int length, int inv, double *sin_tbl, double *cos_tbl);
  void birv(double *a,int length,int ex,double *b);
  void rvmtx1(double a[NBINY][NBINX],double b[NBINX][NBINY],int xsize,int ysize);
  void rvmtx2(double a[NBINX][NBINY],double b[NBINY][NBINX],int xsize,int ysize);

  void palette_b(bool yesorno);
  void Print();
  TApplication *theApp;
  
  TCanvas *c1 ;
  TCanvas *c2 ;
  TCanvas *c3 ;
  TCanvas *c4 ;
  
};

#endif

