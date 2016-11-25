#include "SBPComptonImaging.hh"

//#include "com.h"
//#include "cli.h"

#include "NextCLI.hh"

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TMath.h"

#include <iostream>

using namespace anl;

int main()
{
  int nIteration = 0; // Back Projection
  //  CLintrd("Number of Iteration : ", &nIteration);

  double weight = 0.0;

  double fwhm = 2.4;
  std::string inFileName ="twohit.root";
  std::string effFileName="efficiency.root";
  std::string outFileName="image.root";
  //  strcpy(inFileName, "twohit.root");
  //  strcpy(effFileName, "efficiency.root");
  //  strcpy(outFileName, "image.root");
  CLread("Input file : ", &inFileName);
  CLread("efficiency file name : ", &effFileName);
  CLread("Output file : ", &outFileName);
  CLread("Angular resolution FWHM [degree]", &fwhm);
  CLread("BackProjection weight UEP:0.0 SEP:2.0 ", &weight);

  TFile* inFile = new TFile(inFileName.c_str());
  TFile* effFile = new TFile(effFileName.c_str());
  TFile* outFile = new TFile(outFileName.c_str(), "recreate");
  TTree* tree = (TTree*)inFile->Get("comptree");
  TH3F *s_eff = (TH3F *)effFile->Get("efficiency");
  int nEvent = tree->GetEntries();


  SBPComptonImaging img(nEvent, nIteration);

  img.setEfficiencyMap(s_eff);
  img.createImageHist(outFile);
  img.SetResolutionFunc(fwhm);
  img.backProjection(tree,weight);

  outFile->Write();
  outFile->Close();
  inFile->Close();

  return 0;
}
