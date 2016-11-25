#ifndef OneEventResponse_HH
#define OneEventResponse_HH

#include <iostream>
//#include <map>
#include <vector>
#include <algorithm>

#include "TH3.h"
#include "TwoHitComptonEvent.hh"

class OneEventResponse
{

public:
  OneEventResponse();
  ~OneEventResponse(){};

public:
  void MakeOneEventResponse(TwoHitComptonEvent& twohit, TH3F *Tmp, TH3F* eff, TH3F *backprojection, double *vis);
  float SearchNonZeroElement(int gbin);
  float ExpectionOfDataCount(const TH3F* orgimage);

private:
  //  static const int Nspace =6;
  //  std::map<int, float> Response;
  std::vector<int> Response;



};

#endif
