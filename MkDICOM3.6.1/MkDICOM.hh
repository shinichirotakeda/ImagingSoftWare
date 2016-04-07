#ifndef MkDICOM_hh
#define MkDICOM_hh

#include <vector>
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>

#include "TCanvas.h"
#include "TH3.h"
#include "TFile.h"

#include "ChartStack.hh"




class MkDICOM
{

public:
  MkDICOM(){}
  ~MkDICOM(){}
  


private:

  TH3 *sbp;
  double h_maxvalue;
  double h_minvalue;

  std::string p_name;
  std::string base_name;
  std::string n_study;

  std::string slice;
  int sli_dir;
  std::string name_folder;
  //  int number;
  //  int sift;
  //  int bit;

  std::vector<ChartStack *> ChartStackVec;
  void DICOMFileFactory(ChartStack *s);
  void AddZeroHedder(ofstream *fout,const ChartStack *s);
  void Add0002Tag(ofstream *fout,const ChartStack *s);
  void Add0008Tag(ofstream *fout,const ChartStack *s);
  void Add0010Tag(ofstream *fout,const ChartStack *s);
  void Add0018Tag(ofstream *fout,const ChartStack *s);
  void Add0020Tag(ofstream *fout,const ChartStack *s);
  void Add0028Tag(ofstream *fout,const ChartStack *s);
  void Add0032Tag(ofstream *fout,const ChartStack *s);
  void Add0040Tag(ofstream *fout,const ChartStack *s);
  void AddPixelData(ofstream *fout,const ChartStack *s);

  void Set3DHist(TH3 *h3);

public:

  void Init(); 
  void MkChart();
  void DICOMFileGen();
  void PrintChart();


};

#endif


