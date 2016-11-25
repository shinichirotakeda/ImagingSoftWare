#include "OneEventResponse.hh"
#include "TMath.h"

OneEventResponse::OneEventResponse()
{
 
}


void OneEventResponse::MakeOneEventResponse(TwoHitComptonEvent& twohit,TH3F *Tmp, TH3F* eff, TH3F *backprojection, double *vis)
{
  int TmpBinX = Tmp->GetXaxis()->GetNbins();
  int TmpBinY = Tmp->GetYaxis()->GetNbins();
  int TmpBinZ = Tmp->GetZaxis()->GetNbins();


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
      for(int k=1;k<=TmpBinZ;k++){

 
	pi.SetXYZ(Tmp->GetXaxis()->GetBinCenter(i),Tmp->GetYaxis()->GetBinCenter(j),Tmp->GetZaxis()->GetBinCenter(k));
	diri = pi-ver;
	phii = axis.Angle(diri);
	chii = TMath::Abs(TMath::Pi()-thetaE-phii);
	mi = (diri.Mag()) * TMath::Sin(chii);
	rcrit = (Tmp->GetXaxis()->GetBinWidth(1) ) * 0.866;
	
	if(mi<rcrit){
	  int gbin = Tmp->GetBin(i,j,k);
	  Tmp->AddBinContent(gbin,1.0);
	  backprojection->AddBinContent(gbin,1.0);
	  //	  Response.insert( std::map<int, float>::value_type(gbin,1.0*eff->GetBinContent(gbin)) );
	  Response.push_back(gbin);
	}

      }      
    }
  }

  *vis = Tmp->Integral(1,TmpBinX,1,TmpBinY,1,TmpBinZ);
  std::sort(Response.begin(),Response.end()); //sort
}

float OneEventResponse::SearchNonZeroElement(int gbin)
{
  if(Response.size()==0){
    return -1.0;
  }

  //std::map<unsigned short int, unsigned char>::iterator it = Response.begin();
  std::vector<int>::iterator it = find(Response.begin(),Response.end(), gbin);

  if(it==Response.end()){
    return -1.0;
  }else{
    return 1.0;
  }

}

float OneEventResponse::ExpectionOfDataCount(const TH3F* orgimage)
{
  float sum =0.;
  //  std::map<unsigned short int, unsigned char >::iterator it = Response.begin();
  std::vector<int>::iterator it = Response.begin();

  while(it!=Response.end()){
    int gbin = (*it);
      //    float value = (float) (*it).second;
    sum += (orgimage->GetBinContent(gbin)) * 1.0;
    it++;
  }

  return sum; 
 
}

