#include "PSdetection.hh"
#include <iostream>

#include "TMath.h"
//#include "omp.h"

void PSdetection::Init(){

  x_sigma = 0.02;
  y_sigma = 0.02;
  //  confidence = 5.;
}

void PSdetection::SetSigma(double dx,double dy){

  x_sigma = dx;
  y_sigma = dy;

  x_sigma = x_sigma/xbinw; // pixel dimension
  y_sigma = y_sigma/ybinw; // pixel dimension

  x_range = 5.*x_sigma;
  y_range = 5.*y_sigma;


}

void PSdetection::Save(TDirectory *dir)
{
  std::cout << "data save" << std::endl;
  dir->cd();
  orgimage->Write();
  correlationdata->Write();
  background->Write();

}

void PSdetection::Run()
{

  CorrelationWithMHFunctionData(); // pixel dimension
  BackGroundEstimation(); // backgound 
  // ExtimatedSorcePosition();
  // CalcNewImage();
}

void PSdetection::SetInitImage(TH2D *his)
{
  orgimage = his;

  int nbinx = orgimage->GetXaxis()->GetNbins();
  int nbiny = orgimage->GetYaxis()->GetNbins();

  double xmin = orgimage->GetXaxis()->GetXmin();
  double xmax = orgimage->GetXaxis()->GetXmax();

  double ymin = orgimage->GetYaxis()->GetXmin();
  double ymax = orgimage->GetYaxis()->GetXmax();

  xbinw = orgimage->GetXaxis()->GetBinWidth(1);
  ybinw = orgimage->GetYaxis()->GetBinWidth(1);

  correlationdata = new TH2D("correlationdata","correlationdata",nbinx,xmin,xmax,nbiny,ymin,ymax);
  background = new TH2D("background","background",nbinx,xmin,xmax,nbiny,ymin,ymax);

}


void PSdetection::CorrelationWithMHFunctionData()
{

  std::cout << "calculating WT of orgimage ...." << std::endl;

  int nbinx = orgimage->GetXaxis()->GetNbins();
  int nbiny = orgimage->GetYaxis()->GetNbins();




  std::cout << x_sigma << " pixel " << std::endl; 
  std::cout << y_sigma << " pixel " << std::endl; 

#pragma omp parallel num_threads(2)
  {
#pragma omp for
  for(int i=1;i<=nbinx;i++){
    for(int j=1;j<=nbiny;j++){
      int gbin = orgimage->GetBin(i,j);

      // Correlation 
      double corre = 0.;
      for(int ii=1;ii<=nbinx;ii++){

	double x_dist = TMath::Abs((double)i-(double)ii);
	if(x_dist>x_range) continue;

	for(int jj=1;jj<=nbiny;jj++){

	  double y_dist = TMath::Abs((double)j-(double)jj);
	  if(y_dist>y_range) continue;

	  double value = orgimage->GetBinContent(ii,jj);
	  double x1 = i-ii - 0.5;
	  double x2 = i-ii + 0.5;
	  double y1 = j-jj - 0.5;
	  double y2 = j-jj + 0.5;

	  double x_erf = TMath::Erf(x2/TMath::Sqrt(2)/x_sigma) - TMath::Erf(x1/TMath::Sqrt(2)/x_sigma); 
	  double y_erf = TMath::Erf(y2/TMath::Sqrt(2)/y_sigma) - TMath::Erf(y1/TMath::Sqrt(2)/y_sigma); 

	  double x_diff = x1*TMath::Exp(-x1*x1/2./x_sigma/x_sigma) - x2*TMath::Exp(-x2*x2/2./x_sigma/x_sigma);
	  double y_diff = y1*TMath::Exp(-y1*y1/2./y_sigma/y_sigma) - y2*TMath::Exp(-y2*y2/2./y_sigma/y_sigma);


	  double w = TMath::Pi()*x_sigma*y_sigma*x_erf*y_erf - 
	    TMath::Sqrt(TMath::Pi()/2.)*y_sigma*y_erf*(x_diff + TMath::Sqrt(TMath::Pi()/2.)*x_sigma*x_erf) -
	    TMath::Sqrt(TMath::Pi()/2.)*x_sigma*x_erf*(y_diff + TMath::Sqrt(TMath::Pi()/2.)*y_sigma*y_erf);

	  corre = corre + w*value;


	}
      }

      correlationdata->SetBinContent(gbin,corre);

    }
  }

  }

}


void PSdetection::BackGroundEstimation()
{
  std::cout << "Estimaging Background File...." << std::endl;

  int nbinx = orgimage->GetXaxis()->GetNbins();
  int nbiny = orgimage->GetYaxis()->GetNbins();

  std::cout << x_sigma << " pixel " << std::endl; 
  std::cout << y_sigma << " pixel " << std::endl; 

#pragma omp parallel num_threads(2)
  {
#pragma omp for
  for(int i=1;i<=nbinx;i++){
    for(int j=1;j<=nbiny;j++){
      int gbin = orgimage->GetBin(i,j);

      // Correlation with NW
      double corre = 0.;
      for(int ii=1;ii<=nbinx;ii++){

	double x_dist = TMath::Abs((double)i-(double)ii);
	if(x_dist>x_range) continue;

	for(int jj=1;jj<=nbiny;jj++){

	  double y_dist = TMath::Abs((double)j-(double)jj);
	  if(y_dist>y_range) continue;

	  double value = orgimage->GetBinContent(ii,jj);
	  double x1 = i-ii - 0.5;
	  double x2 = i-ii + 0.5;
	  double y1 = j-jj - 0.5;
	  double y2 = j-jj + 0.5;

	  double x_erf = TMath::Erf(x2/TMath::Sqrt(2)/x_sigma) - TMath::Erf(x1/TMath::Sqrt(2)/x_sigma); 
	  double y_erf = TMath::Erf(y2/TMath::Sqrt(2)/y_sigma) - TMath::Erf(y1/TMath::Sqrt(2)/y_sigma); 

	  double x_diff = x1*TMath::Exp(-x1*x1/2./x_sigma/x_sigma) - x2*TMath::Exp(-x2*x2/2./x_sigma/x_sigma);
	  double y_diff = y1*TMath::Exp(-y1*y1/2./y_sigma/y_sigma) - y2*TMath::Exp(-y2*y2/2./y_sigma/y_sigma);


	  double w = TMath::Pi()*x_sigma*y_sigma*x_erf*y_erf - 
	    TMath::Sqrt(TMath::Pi()/2.)*y_sigma*y_erf*(x_diff + TMath::Sqrt(TMath::Pi()/2.)*x_sigma*x_erf) -
	    TMath::Sqrt(TMath::Pi()/2.)*x_sigma*x_erf*(y_diff + TMath::Sqrt(TMath::Pi()/2.)*y_sigma*y_erf);

	  if(w<0.)corre = corre - w*value;


	}
      }

      background->SetBinContent(gbin,TMath::Exp(1)/(4.*TMath::Pi()*x_sigma*y_sigma)*corre);

    }
  }
  }

}

/*
void PSdetection::CalcNewImage(){

  int nbinx = orgimage->GetXaxis()->GetNbins();
  int nbiny = orgimage->GetYaxis()->GetNbins();

  std::cout << x_sigma << " pixel " << std::endl; 
  std::cout << y_sigma << " pixel " << std::endl; 

  for(int i=1;i<=nbinx;i++){
    for(int j=1;j<=nbiny;j++){
      int gbin = orgimage->GetBin(i,j);
      newimage->SetBinContent(gbin,(orgimage->GetBinContent(gbin)) - (background->GetBinContent(gbin)));
    }
  }


}
*/
/*
void PSdetection::ExtimatedSorcePosition()
{
  int nbinx = orgimage->GetXaxis()->GetNbins();
  int nbiny = orgimage->GetYaxis()->GetNbins();
  double minimum =100000.;
  for(int i=1;i<=nbinx;i++){
    for(int j=1;j<=nbiny;j++){

      sourceposition3->SetBinContent(i,j,correlationdata->GetBinContent(i,j));
      sourceposition4->SetBinContent(i,j,correlationdata->GetBinContent(i,j));
      sourceposition5->SetBinContent(i,j,correlationdata->GetBinContent(i,j));
      sourceposition6->SetBinContent(i,j,correlationdata->GetBinContent(i,j));
      double x_center = correlationdata->GetXaxis()->GetBinCenter(i);
      double y_center = correlationdata->GetYaxis()->GetBinCenter(j);

      if(TMath::Sqrt(x_center*x_center + y_center*y_center)<=1.){

	double value = correlationdata->GetBinContent(i,j);
	if(value<minimum) minimum = value;
      }
    }
  }

  std::cout << "Minimum value is " << minimum << std::endl;

  double confidence = - (3./4.1 * minimum);
  sourceposition3->SetMinimum(confidence);
  confidence = - (4./4.1 * minimum);
  sourceposition4->SetMinimum(confidence);
  confidence = - (5./4.1 * minimum);
  sourceposition5->SetMinimum(confidence);
  confidence = - (6./4.1 * minimum);
  sourceposition6->SetMinimum(confidence);
}
*/

/*
void PSdetection::CorrelationWithMHFunctionData()
{

  int nbinx = orgimage->GetXaxis()->GetNbins();
  int nbiny = orgimage->GetYaxis()->GetNbins();

  double xbinw = orgimage->GetXaxis()->GetBinWidth(1);
  double ybinw = orgimage->GetYaxis()->GetBinWidth(1);


  for(int i=1;i<=nbinx;i++){
    for(int j=1;j<=nbiny;j++){
      int gbin = orgimage->GetBin(i,j);
      double xcenter = orgimage->GetXaxis()->GetBinCenter(i);
      double ycenter = orgimage->GetYaxis()->GetBinCenter(j);

      std::cout << i << " " << j << std::endl;

      // Correlation 
      double corre = 0.;
      for(int ii=1;ii<=nbinx;ii++){

	double xposition = orgimage->GetXaxis()->GetBinCenter(ii);
	double x_dist = TMath::Abs(xcenter-xposition);
	if(x_dist>x_range) continue;

	for(int jj=1;jj<=nbiny;jj++){
	  double yposition = orgimage->GetYaxis()->GetBinCenter(jj);
	  double y_dist = TMath::Abs(ycenter-yposition);
	  if(y_dist>y_range) continue;

	  double value = orgimage->GetBinContent(ii,jj);
	  double x1 = xposition-xcenter - 0.5*xbinw;
	  double x2 = xposition-xcenter + 0.5*xbinw;
	  double y1 = yposition-ycenter - 0.5*ybinw;
	  double y2 = yposition-ycenter + 0.5*ybinw;

	  double x_erf = TMath::Erf(x2/TMath::Sqrt(2)/x_sigma) - TMath::Erf(x1/TMath::Sqrt(2)/x_sigma); 
	  double y_erf = TMath::Erf(y2/TMath::Sqrt(2)/y_sigma) - TMath::Erf(y1/TMath::Sqrt(2)/y_sigma); 

	  double x_diff = x1*TMath::Exp(-x1*x1/2./x_sigma/x_sigma) - x2*TMath::Exp(-x2*x2/2./x_sigma/x_sigma);
	  double y_diff = y1*TMath::Exp(-y1*y1/2./y_sigma/y_sigma) - y2*TMath::Exp(-y2*y2/2./y_sigma/y_sigma);


	  double w = TMath::Pi()*x_sigma*y_sigma*x_erf*y_erf - 
	    TMath::Sqrt(TMath::Pi()/2.)*y_sigma*y_erf*(x_diff + TMath::Sqrt(TMath::Pi()/2.)*x_sigma*x_erf) -
	    TMath::Sqrt(TMath::Pi()/2.)*x_sigma*x_erf*(y_diff + TMath::Sqrt(TMath::Pi()/2.)*y_sigma*y_erf);

	  corre = corre + w*value;


	}
      }

      correlationdata->SetBinContent(gbin,corre);

    }
  }

}
*/
