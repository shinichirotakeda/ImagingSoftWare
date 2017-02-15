#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TH2.h"
#include "TDirectory.h"

#include "ImageFactory.hh"

#include "NextCLI.hh"

//#include "com.h"
//#include "cli.h"

using namespace anl;

int main(int argc, char **argv){

  TApplication theApp("App", &argc, argv);
  ImageFactory *image_factory = new ImageFactory();
  
  image_factory->Init(&theApp);
  
  
  std::string name="input.root";
  CLread("INPUT FILE NAME ? ", &name);
  TFile *fin = new TFile(name.c_str());
  
  std::string h_name="image";
  CLread("IMAGE NAME ? ", &h_name);
  TH2F *orgimage = (TH2F *)fin->Get(h_name.c_str());
  
  std::string name_psf="psf.root";
  CLread("PSF FILE NAME ? ", &name_psf);
  TFile *fpsf = new TFile(name_psf.c_str());
  
  std::string psf_name="psf";
  CLread("PSF NAME ? ", &psf_name);
  TH2F *psf = (TH2F *)fpsf->Get(psf_name.c_str());

  std::string n_fout = "output.root";
  CLread("OUTPUT ROOT FILE NAME ?",&n_fout);
  TFile *fout = new TFile(n_fout.c_str(),"recreate");
  TDirectory *dir = fout->GetDirectory(0);
  
  image_factory->SetInitImage(orgimage);
  image_factory->SetPsfImage(psf);
  
  image_factory->Run();
  image_factory->Save(dir);



  return 0;

}


