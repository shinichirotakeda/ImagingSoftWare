#include "ImageFactory.hh"

#include <iostream>

//#include "com.h"
//#include "cli.h"

#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"
#include "TF1.h"
#include "TPad.h"
#include "NextCLI.hh"


using namespace anl;


void ImageFactory::Init(TApplication *app)
{
  flag = false;
  std::cout << "Init ......." << std::endl;
  c1 = new TCanvas("c1","c1",0,0,500,500);
  c2 = new TCanvas("c2","c2",0,500,500,500);
  c3 = new TCanvas("c3","c3",500,0,500,500);
  c4 = new TCanvas("c4","c4",500,500,500,500);
  c1->Divide(2,2);
  c2->Divide(2,2);
  c3->Divide(2,2);
  c4->Divide(2,2);
  theApp = app;
}

void ImageFactory::SetInitImage(TH2D *his)
{
  orgimage = his;
  fft_orgimage = new TH2D("FFTPower","FFTPower",NBINX,-NBINX/2,NBINX/2,NBINY,-NBINY/2,NBINY/2);
  fft_realpart = new TH2D("FFTREAL","FFTREAL",NBINX,-NBINX/2,NBINX/2,NBINY,-NBINY/2,NBINY/2);
  fft_impart = new TH2D("FFTIM","FFTIM",NBINX,-NBINX/2,NBINX/2,NBINY,-NBINY/2,NBINY/2);
  filtered_image = new TH2D("FilteredImage","Filtered Image",NBINX,orgimage->GetXaxis()->GetXmin(),orgimage->GetXaxis()->GetXmax(),
			    NBINY,orgimage->GetYaxis()->GetXmin(),orgimage->GetYaxis()->GetXmax());

  lowpassfilter = new TH2D("lowpassfilter","lowpassfilter",NBINX,-NBINX/2,NBINX/2,NBINY,-NBINY/2,NBINY/2);
  highpassfilter = new TH2D("highpassfilter","highpassfilter",NBINX,-NBINX/2,NBINX/2,NBINY,-NBINY/2,NBINY/2);
  bandpassfilter = new TH2D("bandpassfilter","bandpassfilter",NBINX,-NBINX/2,NBINX/2,NBINY,-NBINY/2,NBINY/2);

  power = new TH2D("power","power",NBINX,-NBINX/2,NBINX/2,NBINY,-NBINY/2,NBINY/2);

  //  powerspectrum = new TH1D("powerspectrum","powerspectrum",100,0,100);

}

void ImageFactory::SetPsfImage(TH2D *his)
{
  psf  = his;
  psf_fft = new TH2D("PSFPower","PSFPower",NBINX,-NBINX/2,NBINX/2,NBINY,-NBINY/2,NBINY/2);
  psf_realpart = new TH2D("PSFREAL","PSFREAL",NBINX,-NBINX/2,NBINX/2,NBINY,-NBINY/2,NBINY/2);
  psf_impart = new TH2D("PSFIM","PSFTIM",NBINX,-NBINX/2,NBINX/2,NBINY,-NBINY/2,NBINY/2);

}

void ImageFactory::Run()
{
  int num;

  std::cout << "Run ...... select number" << std::endl;

 runstart:


  std::cout << " FFT and Reverse (1)" << std::endl;
  std::cout << " FFT and LowPass Filtering (2)" << std::endl;
  std::cout << " FFT and HighPass Filtering (3)" << std::endl;
  std::cout << " FFT and BandPass Filtering (4)" << std::endl;
  std::cout << " FFT and Deconvolve d (5)" << std::endl;
  std::cout << " Print (6)" << std::endl;
  std::cout << " EXIT (0)" << std::endl;
  std::cout << " > " ;
  std::cin >> num ;

  if(num == 6){
    Print();
    goto runstart;
  }else if(num==1){
    FFT2();
    FFT2Inv();
    Print();
    goto runstart;
  }else if(num == 2){
    FFT2();
    FFT2LowPassFilter();
    Print();
    goto runstart;
  }else if(num == 3){
    FFT2();
    FFT2HighPassFilter();
    Print();
    goto runstart;

  }else if(num == 4){
    FFT2();
    FFT2BandPassFilter();
    Print();
    goto runstart;
  }else if(num == 5){
    FFT2();
    FFT2Deconvolved();
    Print();
    goto runstart;

  }else if(num == 0){

  }else{
    goto runstart;
  }


}

void ImageFactory::Print()
{
  std::cout << "Print ......." << std::endl;
  palette_b(false);

  c1->cd(1);
  orgimage->Draw("colz");
  c1->cd(2);
  fft_orgimage->Draw("colz");
  c1->cd(3);
  fft_realpart->Draw("colz");
  c1->cd(4);
  fft_impart->Draw("colz");

  c2->cd(1);
  psf->Draw("colz");
  c2->cd(2);
  psf_fft->Draw("colz");
  c2->cd(3);
  psf_realpart->Draw("colz");
  c2->cd(4);
  psf_impart->Draw("colz");

  c3->cd(1);
  lowpassfilter->Draw("colz");
  c3->cd(2);
  highpassfilter->Draw("colz");
  c3->cd(3);
  bandpassfilter->Draw("colz");
  c3->cd(4);
  filtered_image->Draw("colz");

  c4->cd(1);
  power->Draw("colz");
  c4->cd(2);
  if(flag)  powerspectrum->Draw("AP");

  c1->Update();  
  c2->Update();  
  c3->Update();  
  c4->Update();
  theApp->Run(kTRUE);

}


void ImageFactory::Save(TDirectory *dir)
{
  std::cout << "data save" << std::endl;
  dir->cd();
  std::cout << "dir ->cd ()" << std::endl;
  orgimage->Write();
  fft_orgimage->Write();
  fft_realpart->Write();
  fft_impart->Write();
  filtered_image->Write();
  lowpassfilter->Write();
  highpassfilter->Write();
  bandpassfilter->Write();

  psf->Write();
  psf_fft->Write();
  psf_realpart->Write();
  psf_impart->Write();

  if(flag)  powerspectrum->Write();
  power->Write();
}


bool ImageFactory::FFT2()
{

  int n_xbin = orgimage->GetXaxis()->GetNbins();
  int n_ybin = orgimage->GetYaxis()->GetNbins();

  if(n_xbin != NBINX || n_ybin != NBINY){
    std::cout << "Number of Bins shold be " << NBINX << " : " << NBINY << std::endl; 
    return false;
  }

  double norm,max;
  // double data;

  double a_rl[NBINY][NBINX];
  double a_im[NBINY][NBINX];

  for(int i=1;i<=n_ybin;i++){
    for(int j=1;j<=n_xbin;j++){
      a_rl[i-1][j-1] = orgimage->GetBinContent(j,i);
      a_im[i-1][j-1] = 0.;
    }
  }

  FFT2core(a_rl,a_im,1);

  max = 0.;
  for(int i=0;i<NBINY;i++){
    for(int j=0;j<NBINX;j++){
      fft_realpart ->SetBinContent(j+1,i+1,a_rl[i][j]);
      fft_impart ->SetBinContent(j+1,i+1,a_im[i][j]);
      norm = a_rl[i][j]*a_rl[i][j] + a_im[i][j]*a_im[i][j];
      //      if(norm != 0.0){
      //	norm = TMath::Log(norm) / 2.0;
      //      }else{
      //	norm = 0.0;
      //      }

      fft_orgimage->SetBinContent(j+1,i+1,norm);
      if(norm>max) max = norm;

    }
  }

  return true;

}

bool ImageFactory::FFT2Inv()
{

  double a_rl[NBINY][NBINX];
  double a_im[NBINY][NBINX];

  for(int i=1;i<=NBINY;i++){
    for(int j=1;j<=NBINX;j++){
      a_rl[i-1][j-1] = fft_realpart->GetBinContent(j,i);
      a_im[i-1][j-1] = fft_impart->GetBinContent(j,i);
    }
  }
  FFT2core(a_rl,a_im,-1);
  for(int i=0;i<NBINY;i++){
    for(int j=0;j<NBINX;j++){
      filtered_image->SetBinContent(j+1,i+1,a_rl[i][j]);
    }
  }

  return true;
}




bool ImageFactory::FFT2LowPassFilter()
{

  double range = (double ) NBINX/2/2;
  CLread("Filter Range (Default is NBINX/4) ", &range);

  for(int i=1;i<=NBINX;i++){
    for(int j=1;j<=NBINX;j++){
      double xcenter = lowpassfilter->GetXaxis()->GetBinCenter(i);
      double ycenter = lowpassfilter->GetYaxis()->GetBinCenter(j);
      double distance = TMath::Sqrt(xcenter*xcenter + ycenter*ycenter);
      if(distance <= range){
	lowpassfilter-> SetBinContent(i,j,1.0);
      }else{
	lowpassfilter-> SetBinContent(i,j,0.0);
      }
    }
  }


  double a_rl[NBINY][NBINX];
  double a_im[NBINY][NBINX];

  for(int i=1;i<=NBINY;i++){
    for(int j=1;j<=NBINX;j++){
      a_rl[i-1][j-1] = fft_realpart->GetBinContent(j,i)*lowpassfilter->GetBinContent(j,i);
      a_im[i-1][j-1] = fft_impart->GetBinContent(j,i)*lowpassfilter->GetBinContent(j,i) ;
    }
  }

  FFT2core(a_rl,a_im,-1);
  for(int i=0;i<NBINY;i++){
    for(int j=0;j<NBINX;j++){
      filtered_image->SetBinContent(j+1,i+1,a_rl[i][j]);
    }
  }

  return true;

}

bool ImageFactory::FFT2HighPassFilter()
{

 double range = (double ) NBINX/2/2;
  CLread("Filter Range (Default is NBINX/4) ", &range);

  for(int i=1;i<=NBINX;i++){
    for(int j=1;j<=NBINX;j++){
      double xcenter = lowpassfilter->GetXaxis()->GetBinCenter(i);
      double ycenter = lowpassfilter->GetYaxis()->GetBinCenter(j);
      double distance = TMath::Sqrt(xcenter*xcenter + ycenter*ycenter);
      if(distance >= range){
	highpassfilter-> SetBinContent(i,j,1.0);
      }else{
	highpassfilter-> SetBinContent(i,j,0.0);
      }
    }
  }


  double a_rl[NBINY][NBINX];
  double a_im[NBINY][NBINX];

  for(int i=1;i<=NBINY;i++){
    for(int j=1;j<=NBINX;j++){
      a_rl[i-1][j-1] = fft_realpart->GetBinContent(j,i)*highpassfilter->GetBinContent(j,i);
      a_im[i-1][j-1] = fft_impart->GetBinContent(j,i)*highpassfilter->GetBinContent(j,i) ;
    }
  }
  FFT2core(a_rl,a_im,-1);
  for(int i=0;i<NBINY;i++){
    for(int j=0;j<NBINX;j++){
      filtered_image->SetBinContent(j+1,i+1,a_rl[i][j]);
    }
  }

  return true;

}



bool ImageFactory::FFT2BandPassFilter()
{

  double ld = (double ) NBINX/32;
  CLread("Filter Range LD (Default is NBINX/32) ", &ld);

  double ud = (double ) NBINX/8;
  CLread("Filter Range UD (Default is NBINX/8) ", &ud);

  for(int i=1;i<=NBINX;i++){
    for(int j=1;j<=NBINX;j++){
      double xcenter = bandpassfilter->GetXaxis()->GetBinCenter(i);
      double ycenter = bandpassfilter->GetYaxis()->GetBinCenter(j);
      double distance = TMath::Sqrt(xcenter*xcenter + ycenter*ycenter);
      if(distance >= ld && ud >= distance){
	bandpassfilter-> SetBinContent(i,j,1.0);
      }else{
	bandpassfilter-> SetBinContent(i,j,0.0);
      }
    }
  }


  double a_rl[NBINY][NBINX];
  double a_im[NBINY][NBINX];

  for(int i=1;i<=NBINY;i++){
    for(int j=1;j<=NBINX;j++){
      a_rl[i-1][j-1] = fft_realpart->GetBinContent(j,i)*bandpassfilter->GetBinContent(j,i);
      a_im[i-1][j-1] = fft_impart->GetBinContent(j,i)*bandpassfilter->GetBinContent(j,i) ;
    }
  }
  FFT2core(a_rl,a_im,-1);
  for(int i=0;i<NBINY;i++){
    for(int j=0;j<NBINX;j++){
      filtered_image->SetBinContent(j+1,i+1,a_rl[i][j]);
    }
  }

  return true;

}

bool ImageFactory::FFT2Deconvolved()
{
  flag=true;
  // FFT OF PSF

  int n_xbin = psf->GetXaxis()->GetNbins();
  int n_ybin = psf->GetYaxis()->GetNbins();

  if(n_xbin != NBINX || n_ybin != NBINY){
    std::cout << "Number of Bins should be " << NBINX << " : " << NBINY << std::endl; 
    return false;
  }

  double norm,max;

  double a_rl[NBINY][NBINX];
  double a_im[NBINY][NBINX];

  for(int i=1;i<=n_ybin;i++){
    for(int j=1;j<=n_xbin;j++){
      a_rl[i-1][j-1] = psf->GetBinContent(j,i);
      a_im[i-1][j-1] = 0.;
    }
  }

  FFT2core(a_rl,a_im,1);

  max = 0.;
  for(int i=0;i<NBINY;i++){
    for(int j=0;j<NBINX;j++){
      psf_realpart ->SetBinContent(j+1,i+1,a_rl[i][j]);
      psf_impart ->SetBinContent(j+1,i+1,a_im[i][j]);
      norm = a_rl[i][j]*a_rl[i][j] + a_im[i][j]*a_im[i][j];
      //      if(norm != 0.0){
      //	norm = TMath::Log(norm) / 2.0;
      //      }else{
      //	norm = 0.0;
      //      }

      psf_fft->SetBinContent(j+1,i+1,norm);
      if(norm>max) max = norm;

    }
  }


  // Deconvolved
  


  double aa_rl[NBINY][NBINX];
  double aa_im[NBINY][NBINX];

  double ld = 0.;
  double ud  =10.;
  CLread("LD is ",&ld);
  CLread("UD is ",&ud);

  double cutoff = 6.0;
  CLread("Cut off is  ",&cutoff);

  int n_zi = 3;

  TF1 *filt = new TF1("filt","1./(1.+TMath::Power(x/[0],2*[1]))",0,100);
  filt->SetParameter(0,cutoff);
  filt->SetParameter(1,n_zi);

  for(int i=1;i<=NBINY;i++){
    for(int j=1;j<=NBINX;j++){

      double a = fft_realpart->GetBinContent(j,i);
      double b = fft_impart->GetBinContent(j,i);

      double ppp = TMath::Power(a,2)+ TMath::Power(b,2);
      double k =  TMath::Sqrt( TMath::Power (NBINX/2 -j,2) + TMath::Power(NBINY/2 -i ,2));

      x_reg1.push_back(k);
      y_reg1.push_back(ppp);

      double c = psf_realpart->GetBinContent(j,i);
      double d = psf_impart->GetBinContent(j,i);

      if(k<ud){
	if(c*c+d*d!=0.){
	  aa_rl[i-1][j-1] = (a*c+b*d)/(c*c+d*d);
	  aa_im[i-1][j-1] = (b*c-a*d)/(c*c+d*d);	
	}else{
	  std::cout << "error, devided by zero " << std::endl;
	  aa_rl[i-1][j-1] =0;
	  aa_im[i-1][j-1] =0;
	}
      }else{
	aa_rl[i-1][j-1] =0;
	aa_im[i-1][j-1] =0;
      }

      double factor = filt->Eval(k);

      aa_rl[i-1][j-1] = aa_rl[i-1][j-1] * factor;
      aa_im[i-1][j-1] = aa_im[i-1][j-1] * factor;


      power->SetBinContent(j,i,TMath::Power(aa_rl[i-1][j-1],2)+ TMath::Power(aa_im[i-1][j-1],2));

      if( !(ld<k  &&  k< ud)  ){
      	aa_rl[i-1][j-1]=0;
      	aa_im[i-1][j-1]=0;
      }

    }
  }
  
  delete filt;

  double *x_reg1_pointer =&(x_reg1.at(0));
  double *y_reg1_pointer =&(y_reg1.at(0));
  powerspectrum = new TGraph(x_reg1.size(),x_reg1_pointer,y_reg1_pointer);

  FFT2core(aa_rl,aa_im,-1);

  for(int i=0;i<NBINY;i++){
    for(int j=0;j<NBINX;j++){
      filtered_image->SetBinContent(j+1,i+1,aa_rl[i][j]);
    }
  }

  return true;

}


void ImageFactory::FFT2core(double a_rl[NBINY][NBINX], double a_im[NBINY][NBINX],int inv)
{

  double Xsin_tbl[NBINX];
  double Xcos_tbl[NBINX];
  double Ysin_tbl[NBINY];
  double Ycos_tbl[NBINY];

  double b_rl[NBINX][NBINY];
  double b_im[NBINX][NBINY];

  double buf_x[NBINX];
  double buf_y[NBINY];

  cstb(NBINX,inv,Xsin_tbl,Xcos_tbl);
  cstb(NBINY,inv,Ysin_tbl,Ycos_tbl);



  for(int i=0;i<NBINY;i++){
    FFT1core(&a_rl[i][0],&a_im[i][0],NBINX,X_EXP,Xsin_tbl,Xcos_tbl,buf_x);
  }

  rvmtx1(a_rl,b_rl,NBINX,NBINY);
  rvmtx1(a_im,b_im,NBINX,NBINY);



  for(int i=0;i<NBINX;i++){
    FFT1core(&b_rl[i][0], &b_im[i][0],NBINY,Y_EXP,Ysin_tbl,Ycos_tbl,buf_y);
  }


  rvmtx2(b_rl,a_rl,NBINX,NBINY);
  rvmtx2(b_im,a_im,NBINX,NBINY);
}


void ImageFactory::cstb(int length, int inv, double *sin_tbl, double *cos_tbl)
{

  double xx, arg;
  xx = -TMath::Pi() * 2.0 /(double)length;

  if(inv < 0) xx = -xx;
  for(int i=0;i<length;i++){
    arg = i*xx;
    sin_tbl[i] = TMath::Sin(arg);
    cos_tbl[i] = TMath::Cos(arg);
  }


}

void ImageFactory::FFT1core(double *a_rl, double *a_im,
			    int length, int ex,
			    double *sin_tbl, double *cos_tbl,double *buf)
{
  int i,j,k,w,j1,j2;
  int numb, lenb, timb;
  double xr,xi,yr,yi,nrml;

  if(OPT == 1){
    for(i=1;i<length;i+=2){
      a_rl[i] = -a_rl[i];
      a_im[i] = -a_im[i];
    }
  }

  numb = 1;
  lenb = length;
  for(i = 0;i<ex;i++){
    lenb /= 2;
    timb = 0;
    for(j=0;j<numb;j++){
      w=0;
      for(k=0;k<lenb;k++){
	j1 = timb + k;
	j2 = j1 + lenb;
	xr = a_rl[j1];
	xi = a_im[j1];
	yr = a_rl[j2];
	yi = a_im[j2];
	a_rl[j1] = xr + yr;
	a_im[j1] = xi + yi;
	xr = xr - yr;
	xi = xi - yi;
	a_rl[j2] = xr*cos_tbl[w] - xi*sin_tbl[w];
	a_im[j2] = xr*sin_tbl[w] + xi*cos_tbl[w];
	w += numb;
      }
      timb += (2*lenb);
    }
    numb *=2;
  }

  birv(a_rl,length,ex,buf);
  birv(a_im,length,ex,buf);
  if(OPT == 1){
    for(i=1;i<length;i+=2){
      a_rl[i] = -a_rl[i];
      a_im[i] = -a_im[i];
    }
  }
  nrml = 1.0 / TMath::Sqrt((double)length);
  for(i=0;i<length;i++){
    a_rl[i] *= nrml;
    a_im[i] *= nrml;
  }

}

void ImageFactory::birv(double *a,int length,int ex,double *b)
{

  int i,ii,k,bit;
  for(i=0;i<length;i++){
    for(k=0,ii=i,bit=0; ;bit<<=1,ii>>=1){
      bit = (ii & 1)| bit;
      if(++k == ex) break;
    }
    b[i] = a[bit];
  }

  for(i=0;i<length;i++) a[i]=b[i];


}

void ImageFactory::rvmtx1(double a[NBINY][NBINX],double b[NBINX][NBINY],int xsize,int ysize)
{

  int i,j;

  for(i=0;i<ysize;i++){
    for(j=0;j<xsize;j++){
      //std::cout<< i << " " << j << " " << a[i][j] << std::endl;
      b[j][i] = a[i][j];
    }
  }

}


void ImageFactory::rvmtx2(double a[NBINX][NBINY],double b[NBINY][NBINX],int xsize,int ysize)
{

  int i,j;

  for(i=0;i<ysize;i++){
    for(j=0;j<xsize;j++){
      //std::cout<< i << " " << j << " " << a[i][j] << std::endl;
      b[i][j]  = a[j][i];
    }
  }

}

void ImageFactory::palette_b(bool yesorno){

  const int colNum = 2000;
  int palette[colNum];
  for(int i=0; i<colNum; i++)
    {
      
      float red;
      float green;
      float blue;
      if(i<500)
	{
	  blue = (float)i/500.;
	}
      else if(i<1000)
	{
	  blue = 1.0 - (float)(i-500)/500.;
	}
      else if(i<1500)
	{
	  blue = 0.0;
	}
      else
	{
	  blue = (float)(i-1500)/500.;
	}
      
      if(i<500)
	{
	  red = 0.0;
	}
      else if(i<1000)
	{
	  red = (float)(i-500)/500.;
	}
      else  
	{
	  red = 1.0;
	}

      if(i<1000)
	{
	  green = 0.0;
	}
      else if( i < 1500)
	{
	  green = (float)(i-1000)/500.;
	}
      else
	{
	  green = 1.0;
	}

      if(! gROOT->GetColor(230+i))
	{
	  TColor* color = new TColor(230+i, red, green, blue, "");
	}
      else
	{
	  TColor* color = gROOT->GetColor(230+i);
	  color->SetRGB(red, green, blue);
	}

      if(yesorno)
	{
	  palette[i] = 230 + colNum - i;
	}
      else
	{
	  palette[i] = 230 + i;
	}    

    }

  gStyle->SetPalette(colNum, palette);

}

