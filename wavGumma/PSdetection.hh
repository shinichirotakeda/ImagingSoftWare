#ifndef PSdetection_hh
#define PSdetection_hh

#include "TDirectory.h"
#include "TH2.h"

class PSdetection
{

public:
  PSdetection(){}
  ~PSdetection(){}

private:

  double x_sigma;
  double y_sigma;

  double x_range;
  double y_range;

  double xbinw;
  double ybinw;

  //  double confidence;

  TH2D *orgimage;
  TH2D *newimage;
  TH2D *correlationdata;
  TH2D *correlationdata2;
  TH2D *orgbackground;
  TH2D *background;
  TH2D *sourceposition3;
  TH2D *sourceposition4;
  TH2D *sourceposition5;
  TH2D *sourceposition6;

  void CorrelationWithMHFunctionData();
  void BackGroundEstimation();
  void CalcNewImage();
  void ExtimatedSorcePosition();

public:
  void Init();
  void Save(TDirectory *dir);
  void Run();
  void SetInitImage(TH2D *his);
  void SetSigma(double dx,double dy);

};

#endif
