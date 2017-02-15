#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>

#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH2.h"
#include "TVector3.h"

#include "TwoHitComptonEvent.hh"
#include "OneEventResponse.hh"

class SBPComptonImaging
{
public:
  SBPComptonImaging(int nEvent, int nIteration = 10,int initnum = 0);
  virtual ~SBPComptonImaging();

private:
  SBPComptonImaging();

private:
  const int NUM_EVENT;
  const int NUM_ITERATION;
  int NPixelX;
  int NPixelY;

  double m_ImageWidthX;
  double m_ImageWidthY;

  double m_ImageCenterX;
  double m_ImageCenterY;

  double voigt_sigma;
  double voigt_ld;
  double voigt_fwhm;
  double gauss_fwhm;

  int init_itenum;
  
  TH2F** m_ImageArray; // label:l,j,k
  TH2* m_Efficiency; // label:j
  double zcenter;
  double weight;

  OneEventResponse* m_Response;
  double *m_Visibility;
  

  // double ExpectedDataCount(double* image, int i);
  //  double ImageNewFactor(double* image, int j);
  void improveImage(TH2* newimage, const TH2* orgimage, int ite);

public:
  void ana();
  bool setEfficiencyMap(TH2* h_eff);
  void createImageHist(TFile* file);
  bool makeResponce(TTree* eventtree);
  bool backProjection(TTree* eventtree);
  void SetInitImage(TH2 *init);
  void SetResolutionFunc(double fwhm){gauss_fwhm=fwhm;}
  double calcvoigtfwhm(double sigmaG, double gammaL);
  void setZcenter(double zc){zcenter=zc;}
  void setWeight(double w){weight=w;}

};

