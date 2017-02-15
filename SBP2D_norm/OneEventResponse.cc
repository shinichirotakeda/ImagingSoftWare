#include "OneEventResponse.hh"
#include "TMath.h"

OneEventResponse::OneEventResponse()
{
 
}


void OneEventResponse::MakeOneEventResponse(TwoHitComptonEvent& twohit,TH2 *Tmp, TH2* eff, TH2 *backprojection, double *vis, double zcenter, double gauss_fwhm, double weight)
{
  int TmpBinX = Tmp->GetXaxis()->GetNbins();
  int TmpBinY = Tmp->GetYaxis()->GetNbins();


  double binwidth = eff->GetXaxis()->GetBinWidth(1);
  double sigma = TMath::Pi()/180.* (gauss_fwhm/2.355);

  double  thetaE =  TMath::ACos(twohit.CosThetaE());
  TVector3 axis = twohit.ConeAxis();
  axis = -axis;
  TVector3 ver = twohit.ConeVertex();
  
  TVector3 pi(0.0,0.0,0.0);
  TVector3 diri(0.0,0.0,0.0);
  double rcrit = 0.0;
  double mi = 0.;
  double phii = 0.;
  double chii = 0.;
 
  for(int i=1;i<=TmpBinX;i++){
    for(int j=1;j<=TmpBinY;j++){
 
        pi.SetXYZ(eff->GetXaxis()->GetBinCenter(i),eff->GetYaxis()->GetBinCenter(j), zcenter);
        diri = pi-ver;
        phii = axis.Angle(diri);
        chii = TMath::Abs(TMath::Pi()-thetaE-phii);
        mi = (diri.Mag()) * TMath::Sin(chii) - (diri.Mag()) * TMath::Cos(chii)* TMath::Tan(3.*sigma);
        rcrit = (binwidth) * 0.866;
	
	if(mi<rcrit){
          int gbin = eff->GetBin(i,j);
          double aa = (diri.Mag()) * TMath::Cos(chii)* TMath::Tan(sigma);
	  double value = 1.0/(TMath::Power((diri.Mag()) * TMath::Cos(chii),weight))*TMath::Gaus((diri.Mag()) * TMath::Sin(chii),0,aa);
	  backprojection->AddBinContent(gbin,value);
	  Tmp->AddBinContent(gbin,value);
	  *vis += value;
	  //	  Response.insert( std::map<int, float>::value_type(gbin,(float)value*eff->GetBinContent(gbin)) );
	  Response.insert( std::map<int, float>::value_type(gbin,(float)value) );
	}

            
    }
  }

}

float OneEventResponse::SearchNonZeroElement(int gbin)
{
  if(Response.size()==0){
    return -1.0;
  }

  std::map<int, float>::iterator it = Response.begin();

  it = Response.find(gbin);
  if(it==Response.end()){
    return -1.0;
  }else{
    return (*it).second;
  }

}

float OneEventResponse::ExpectionOfDataCount(const TH2* orgimage)
{


  float sum =0.;
  std::map<int, float>::iterator it = Response.begin();
  
  while(it!=Response.end()){
    int gbin = (*it).first;
    float value = (*it).second;
    sum += (orgimage->GetBinContent(gbin)) * value;
    it++;
  }

  return sum; 
 
}

