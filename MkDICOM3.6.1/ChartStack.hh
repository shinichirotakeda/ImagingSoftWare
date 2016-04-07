#ifndef ChartStack_hh
#define ChartStack_hh

#include "TH2.h"

class ChartStack
{


public:

  void Set2DHist(TH2D *h){hist = h;};
  void SetSliceLocation(double a){slicelocation  = a ;};
  void Print2DHist();

  TH2D *Get2DHist() const {return hist;};
  double GetSlicePos() const {return slicelocation;};
private:

  double slicelocation;
  TH2D *hist;





};


#endif
