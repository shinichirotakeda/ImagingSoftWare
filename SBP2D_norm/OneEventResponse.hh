#ifndef OneEventResponse_HH
#define OneEventResponse_HH

#include <iostream>
#include <map>
//#include <vector>
#include <algorithm>

#include "TH2.h"
#include "TwoHitComptonEvent.hh"

class OneEventResponse
{

public:
  OneEventResponse();
  ~OneEventResponse(){};

public:
  void MakeOneEventResponse(TwoHitComptonEvent& twohit, TH2 *Tmp, TH2* eff, TH2 *backprojection, double *vis, double zcenter, double gauss_fwhm, double weight);
  float SearchNonZeroElement(int gbin);
  float ExpectionOfDataCount(const TH2* orgimage);

private:
  //  static const int Nspace =6;
  std::map<int, float> Response;
  //  std::vector<int> Response;



};

#endif
