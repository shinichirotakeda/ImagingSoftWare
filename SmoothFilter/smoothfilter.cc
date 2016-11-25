#include <iostream>
#include "NextCLI.hh"
#include "TFile.h"
#include "TH2.h"

int main(){


  std::cout << "Hellow " << std::endl;

  std::string inFileName ="input.root";
  anl::CLread("Input file : ", &inFileName);
  std::string imageName ="image_000";
  anl::CLread("Image Name : ", &imageName);
  std::string inFileKernel ="kernel.root";
  anl::CLread("Input file (Kernel): ", &inFileKernel);
  std::string outFileName ="output.root";
  anl::CLread("Output file : ", &outFileName);

  
  TFile* inFile = new TFile(inFileName.c_str());
  TH2D *org = (TH2D *)inFile->Get(imageName.c_str());
  TH2D *smoothed = (TH2D *)org->Clone("Smoothed");
  

  TFile* inKernel = new TFile(inFileKernel.c_str());
  TH2D *kernel = (TH2D *)inKernel->Get("kernel");

  TFile *fout = new TFile(outFileName.c_str(),"recreate");

  int nbinx = org->GetXaxis()->GetNbins();
  int nbiny = org->GetYaxis()->GetNbins();

  int kernelx = kernel->GetXaxis()->GetNbins();
  int kernely = kernel->GetYaxis()->GetNbins();
  if(kernelx%2==0 || kernely%2==0){
    std::cout << "Kernel size should be a odd number" << std::endl;
    return 0;
  }

  for(int i=1;i<=nbinx;i++){
    std::cout << i << std::endl;
    for(int j=1;j<=nbiny;j++){

      double value = 0.;
      int xld = i - (kernelx-1)/2;
      int xud = i + (kernelx-1)/2;
      int yld = j - (kernely-1)/2;
      int yud = j + (kernely-1)/2;

      if(xld<=0) xld = 1;
      if(nbinx<xud) xud = nbinx;
      if(yld<=0) yld = 1;
      if(nbiny<yud) yud = nbiny;

      for(int ii=xld;ii<=xud;ii++){
	for(int jj=yld;jj<=yud;jj++){

	  int dx = ii-i;
	  int dy = jj-j;
	  
	  double data = org->GetBinContent(ii,jj);
	  double factor = kernel->GetBinContent((kernelx+1)/2+dx,(kernely+1)/2+dy);
	  value = value + factor*data;
	}
      }

      smoothed->SetBinContent(i,j,value);
  
    }
  }

  fout->cd();
  smoothed->Write();
  fout->Close();
  return 0;
  
}
