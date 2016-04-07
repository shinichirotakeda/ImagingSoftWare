#include "TROOT.h"
#include "TApplication.h"

#include "MkDICOM.hh"
#include "ChartStack.hh"





int main(int argc, char **argv){

  TApplication theApp("App", &argc, argv);

  MkDICOM *DICOMFactory = new MkDICOM();


  DICOMFactory->Init();
  DICOMFactory->MkChart();
  DICOMFactory->DICOMFileGen();


  return 0;

}


