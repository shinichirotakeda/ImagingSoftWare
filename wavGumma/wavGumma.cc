#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TH2.h"
#include "TDirectory.h"

//#include "com.h"
//#include "cli.h"

#include "PSdetection.hh"
//#include "iostream.h"

#include "NextCLI.hh"

using namespace anl;

int main(int argc, char **argv){
  
  //  #pragma omp parallel
  //  {
  //    std::cout << "Hellow OpneMP" << std::endl;
  //  }
  
  
  TApplication theApp("App", &argc, argv);
  PSdetection *ps_detection = new PSdetection();
  ps_detection->Init();
  
  std::string name="input.root";
  CLread("ROOT FILE NAME ? ", &name);
  std::string image_name="image";
  CLread("Image NAME ? ", &image_name);
  TFile *fin = new TFile(name.c_str());
  TH2D *orgimage = (TH2D *)fin->Get(image_name.c_str());
  //  TH2D *orgimage = (TH2D *)fin->Get("Simple Back Projection");
  ps_detection->SetInitImage(orgimage);
  
  
  double x_sigma = 0.02;
  double y_sigma = 0.02;
  CLread("X sigma ?", &x_sigma);
  CLread("Y sigma ?", &y_sigma);
  ps_detection->SetSigma(x_sigma,y_sigma);
  
  char nn[256];
  sprintf(nn,"output_x%.3f_y%.3f.root",x_sigma,y_sigma);
  std::string n_fout(nn);
  
  CLread("OUTPUT ROOT FILE NAME ?",&n_fout);
  TFile *fout = new TFile(n_fout.c_str(),"recreate");
  TDirectory *dir = fout->GetDirectory(0);
  
  
  
  ps_detection->Run();
  ps_detection->Save(dir);
  return 0;
  
}


