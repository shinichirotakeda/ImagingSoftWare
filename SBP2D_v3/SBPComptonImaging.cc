#include "SBPComptonImaging.hh"

#include "TRandom3.h"
#include "TMath.h"


SBPComptonImaging::SBPComptonImaging(int nEvent, int nIteration, int initnum) :
  NUM_EVENT(nEvent), NUM_ITERATION(nIteration),init_itenum(initnum)
{

  init_itenum = initnum;

  m_ImageArray = new TH2F*[NUM_ITERATION+1];
  m_Response = new OneEventResponse[NUM_EVENT];
  m_Visibility = new double[NUM_EVENT];

}

SBPComptonImaging::~SBPComptonImaging()
{

}


void SBPComptonImaging::createImageHist(TFile* file)
{
  file->cd();
  char name[256];
  for (int l = 0; l<= NUM_ITERATION; l++) {
    sprintf(name, "image_%03d", init_itenum + l);
    m_ImageArray[l] = new TH2F(name, name, NPixelX, m_ImageCenterX-0.5*m_ImageWidthX, m_ImageCenterX+0.5*m_ImageWidthX, 
			       NPixelY, m_ImageCenterY-0.5*m_ImageWidthY, m_ImageCenterY+0.5*m_ImageWidthY);
  }


  for(int i=1;i<=NPixelX;i++){
    for(int j=1;j<=NPixelY;j++){
      int gbin = m_ImageArray[0]->GetBin(i,j);
      m_ImageArray[0]->SetBinContent(gbin,0);
    }
  }
  

}

bool SBPComptonImaging::backProjection(TTree* eventtree)
{
  const int nEvent = NUM_EVENT;
  std::cout << "Number of Events : " << nEvent << std::endl;
  
  if (eventtree->GetEntries() != nEvent) {
    std::cout << "Event tree is invalid." << std::endl;
    return false;
  }



  double binwidth = m_Efficiency->GetXaxis()->GetBinWidth(1);  
  int nbinx = m_Efficiency->GetXaxis()->GetNbins();
  int nbiny = m_Efficiency->GetYaxis()->GetNbins();

  double sigma = TMath::Pi()/180.* (gauss_fwhm/2.355);

  double edep1;
  double posx1, posy1, posz1;
  double edep2;
  double posx2, posy2, posz2;
  eventtree->SetBranchAddress("h1edep",      &edep1);
  eventtree->SetBranchAddress("h1posx",      &posx1);
  eventtree->SetBranchAddress("h1posy",      &posy1);
  eventtree->SetBranchAddress("h1posz",      &posz1);
  eventtree->SetBranchAddress("h2edep",      &edep2);
  eventtree->SetBranchAddress("h2posx",      &posx2);
  eventtree->SetBranchAddress("h2posy",      &posy2);
  eventtree->SetBranchAddress("h2posz",      &posz2);

  
  for (int ev=0; ev<nEvent; ev++) {
    if (ev%100==0) { std::cout << "    event = " << ev << std::endl; }

    eventtree->GetEntry(ev);
    TwoHitComptonEvent twohit;
    twohit.setH1Energy(edep1);
    twohit.setH1Pos(posx1,posy1,posz1);
    twohit.setH2Energy(edep2);
    twohit.setH2Pos(posx2,posy2,posz2);
    
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
 


    for(int i=1;i<=nbinx;i++){
      for(int j=1;j<=nbiny;j++){
        
        pi.SetXYZ(m_Efficiency->GetXaxis()->GetBinCenter(i),m_Efficiency->GetYaxis()->GetBinCenter(j), zcenter);
        diri = pi-ver;
        phii = axis.Angle(diri);
        chii = TMath::Abs(TMath::Pi()-thetaE-phii);
        mi = (diri.Mag()) * TMath::Sin(chii) - (diri.Mag()) * TMath::Cos(chii)* TMath::Tan(3.*sigma);
        rcrit = (binwidth) * 0.866;
        
        if(mi<rcrit){
          int gbin = m_Efficiency->GetBin(i,j);
          double aa = (diri.Mag()) * TMath::Cos(chii)* TMath::Tan(sigma);
	  double value = 1.0/(TMath::Power((diri.Mag()) * TMath::Cos(chii),weight))*TMath::Gaus((diri.Mag()) * TMath::Sin(chii),0,aa);
           m_ImageArray[0]->AddBinContent(gbin,value);
        }
        
      }
      
    }


  }

  return true;

}

double SBPComptonImaging::calcvoigtfwhm(double sigmaG, double gammaL)
{

  double fwhmG = sigmaG * 2.0 * TMath::Sqrt(2.0*TMath::Log(2.));
  double fwhmL = 2.0*gammaL;
  double phi=fwhmL/fwhmG;
  if(phi < 1e-2) return fwhmG;
  if(phi > 100.) return fwhmL;


  double c_0 = 2.0056;
  double c_1 = 1.0593;
  
  return fwhmG * (1.0 - c_0 * c_1 + TMath::Sqrt(phi*phi+2.0*c_1*phi+c_0*c_0*c_1*c_1));
  
}


bool SBPComptonImaging::makeResponce(TTree* eventtree)
{
  const int nEvent = NUM_EVENT;
  std::cout << "Number of Events : " << nEvent << std::endl;
  
  if (eventtree->GetEntries() != nEvent) {
    std::cout << "Event tree is invalid." << std::endl;
    return false;
  }

  double edep1;
  double posx1, posy1, posz1;
  double edep2;
  double posx2, posy2, posz2;
  eventtree->SetBranchAddress("h1edep",      &edep1);
  eventtree->SetBranchAddress("h1posx",      &posx1);
  eventtree->SetBranchAddress("h1posy",      &posy1);
  eventtree->SetBranchAddress("h1posz",      &posz1);
  eventtree->SetBranchAddress("h2edep",      &edep2);
  eventtree->SetBranchAddress("h2posx",      &posx2);
  eventtree->SetBranchAddress("h2posy",      &posy2);
  eventtree->SetBranchAddress("h2posz",      &posz2);

  std::cout << "Make responce " << std::endl;
  TH2F* responceTmp = new TH2F("resp", "resp", NPixelX, m_ImageCenterX-0.5*m_ImageWidthX, m_ImageCenterX+0.5*m_ImageWidthX, 
			       NPixelY, m_ImageCenterY-0.5*m_ImageWidthY, m_ImageCenterY+0.5*m_ImageWidthY);
  
  for (int i=0; i<nEvent; i++) {
    if (i%100==0) { std::cout << "    i = " << i << std::endl; }

    eventtree->GetEntry(i);
    TwoHitComptonEvent twohit;
    twohit.setH1Energy(edep1);
    twohit.setH1Pos(posx1,posy1,posz1);
    twohit.setH2Energy(edep2);
    twohit.setH2Pos(posx2,posy2,posz2);
    
    responceTmp->Reset();
    m_Response[i].MakeOneEventResponse(twohit, responceTmp, m_Efficiency, m_ImageArray[0], &m_Visibility[i],zcenter,gauss_fwhm, weight); 

  }
  
  for(int j=0;j<NUM_EVENT;j++){
    std::cout<< m_Visibility[j] << std::endl;
  }

  responceTmp->Delete();
  return true;
}

bool SBPComptonImaging::setEfficiencyMap(TH2* h_eff)
{

  NPixelX = h_eff->GetXaxis()->GetNbins();
  NPixelY = h_eff->GetYaxis()->GetNbins();

  m_Efficiency = h_eff;

  m_ImageWidthX = ( h_eff->GetXaxis()->GetXmax() ) - ( h_eff->GetXaxis()->GetXmin() ) ;
  m_ImageWidthY = ( h_eff->GetYaxis()->GetXmax() ) - ( h_eff->GetYaxis()->GetXmin() ) ;

  m_ImageCenterX = (( h_eff->GetXaxis()->GetXmax() ) + ( h_eff->GetXaxis()->GetXmin() ) ) /2. ;
  m_ImageCenterY = (( h_eff->GetYaxis()->GetXmax() ) + ( h_eff->GetYaxis()->GetXmin() ) ) /2. ;

  return true;
}

void SBPComptonImaging::SetInitImage(TH2 *init){

  for(int i = 1; i<=NPixelX;i++){
    for(int j = 1; j<=NPixelY;j++){
	int gbin = init->GetBin(i,j);
	m_ImageArray[0] -> SetBinContent(gbin,init->GetBinContent(gbin));
    }
  }



}

void SBPComptonImaging::ana()
{
  std::cout << "LM-EM-ML Imaging start" << std::endl;
  for (int l=0; l<NUM_ITERATION; l++) {
    std::cout << "Iteration l = " << l << std::endl;
    improveImage(m_ImageArray[l+1], m_ImageArray[l],l);
  }
}



void SBPComptonImaging::improveImage(TH2* newimage, const TH2* orgimage, int ite)
{

  int j1, j2;
  double factor=0.;

  for (j1=1; j1<=NPixelX; j1++) {
    for (j2=1; j2<=NPixelY; j2++) {

      factor = 0.0;
      int gbin = orgimage->GetBin(j1,j2);

      for(int k = 0;k<NUM_EVENT;k++){
	float tij = m_Response[k].SearchNonZeroElement(gbin);

	  if(tij > 0.0){
	    tij = tij * m_Visibility[k] / m_Response[k].ExpectionOfDataCount(orgimage);
	    factor += tij;
	  }else{
	    continue;
	  }
	  //	  if(k%100) std::cout << k << std::endl;
	}
	//	std::cout << "calculate ..." << j1 << " " << j2 << " " << j3 << std::endl; 
	newimage->SetBinContent(gbin, orgimage->GetBinContent(gbin) / m_Efficiency->GetBinContent(gbin) * factor);

	
    }
  }
  
}






/*
double SBPComptonImaging::ImageNewFactor(double* image, int j)
{
  double** const t = m_Responce;
  std::vector<double>& v = m_Visibility;
  
  double factor = 0.;
  double* const ex_data = m_ExpectedDataTmp;
  int n = NUM_EVENT;
  for (int i=0; i<n; i++) {
    if (v[i] != 0.0 && ex_data[i] != 0.0) {
      factor += t[i][j] * v[i] / ex_data[i];
    }
  }

  factor = factor/m_Efficiency[j];
  return factor;
}

*/

/*
double SBPComptonImaging::ExpectedDataCount(double* image, int i)
{
  double* const t_i = m_Responce[i];
  
  double count = 0.;
  for (int j=0; j<IMAGE_SIZE; j++) {
  count += t_i[j] * image[j];
  }
  //count += m_Background[i];
  
  return count;
  }
*/
