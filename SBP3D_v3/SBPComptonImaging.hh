#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>

#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH3.h"
#include "TVector3.h"

#include "TwoHitComptonEvent.hh"
#include "OneEventResponse.hh"

class SBPComptonImaging
{
public:
  SBPComptonImaging(int nEvent, int nIteration = 10,int initnum = 0);
  //  SBPComptonImaging(int nEvent, int nIteration = 10, int initnum);
  virtual ~SBPComptonImaging();

private:
  SBPComptonImaging();

private:
  const int NUM_EVENT;
  const int NUM_ITERATION;
  int NPixelX;
  int NPixelY;
  int NPixelZ;

  //  static const int IMAGE_SIZE = NPixel*NPixel*NPixel;

  double m_ImageWidthX;
  double m_ImageWidthY;
  double m_ImageWidthZ;

  double m_ImageCenterX;
  double m_ImageCenterY;
  double m_ImageCenterZ;


  //  double resolution;
  double voigt_sigma;
  double voigt_ld;
  double voigt_fwhm;
  double gauss_fwhm;

  int init_itenum;
  
  // Image space

  TH3F** m_ImageArray; // label:l,j,k
  TH3F* m_Efficiency; // label:j

  // Response

  OneEventResponse* m_Response;
  double *m_Visibility;


  // double ExpectedDataCount(double* image, int i);
  //  double ImageNewFactor(double* image, int j);
  void improveImage(TH3F* newimage, const TH3F* orgimage, int ite);

public:
  void ana();
  bool setEfficiencyMap(TH3F* h_eff);
  void createImageHist(TFile* file);
  bool makeResponce(TTree* eventtree);
  bool backProjection(TTree* eventtree, double weight);
  void SetInitImage(TH3F *init);
  void SetResolutionFunc(double fwhm){gauss_fwhm=fwhm;}
  double calcvoigtfwhm(double sigmaG, double gammaL);


};

