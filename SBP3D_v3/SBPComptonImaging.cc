#include "SBPComptonImaging.hh"

#include "TRandom3.h"
#include "TMath.h"

/*
SBPComptonImaging::SBPComptonImaging(int nEvent, int nIteration) :
  NUM_EVENT(nEvent), NUM_ITERATION(nIteration)
{
  init_itenum = 0;
  m_ImageWidth = 32.;
  m_ImageArray = new TH3F*[NUM_ITERATION+1];
  m_Response = new OneEventResponse[NUM_EVENT];
  m_Visibility = new double[NUM_EVENT];
}
*/

SBPComptonImaging::SBPComptonImaging(int nEvent, int nIteration, int initnum) :
  NUM_EVENT(nEvent), NUM_ITERATION(nIteration),init_itenum(initnum)
{

  init_itenum = initnum;

  m_ImageArray = new TH3F*[NUM_ITERATION+1];
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
    m_ImageArray[l] = new TH3F(name, name, NPixelX, m_ImageCenterX-0.5*m_ImageWidthX, m_ImageCenterX+0.5*m_ImageWidthX, 
			       NPixelY, m_ImageCenterY-0.5*m_ImageWidthY, m_ImageCenterY+0.5*m_ImageWidthY,
			       NPixelZ, m_ImageCenterZ-0.5*m_ImageWidthZ, m_ImageCenterZ+0.5*m_ImageWidthZ);
  }
}

bool SBPComptonImaging::backProjection(TTree* eventtree, double weight)
{
  const int nEvent = NUM_EVENT;
  std::cout << "Number of Events : " << nEvent << std::endl;
  
  if (eventtree->GetEntries() != nEvent) {
    std::cout << "Event tree is invalid." << std::endl;
    return false;
  }

  TH3F *imagespace[4][7];
  TH3F *Tmp[4];
  double binwidth[8];

  binwidth[0] = m_ImageWidthZ;
  binwidth[1] = m_ImageWidthZ/2.;
  binwidth[2] = m_ImageWidthZ/4.;
  binwidth[3] = m_ImageWidthZ/8.;
  binwidth[4] = m_ImageWidthZ/16.;
  binwidth[5] = m_ImageWidthZ/32.;
  binwidth[6] = m_ImageWidthZ/64.;
  binwidth[7] = m_ImageWidthZ/128.;

  char name[256];
  for(int j=0;j<7;j++){
    sprintf(name,"TmpImageSpace%d%d",0,j);
    imagespace[0][j] = new TH3F(name,name,
				(int)TMath::Power(2,j), m_ImageCenterX-0.5*m_ImageWidthX , m_ImageCenterX, 
				(int)TMath::Power(2,j), m_ImageCenterY-0.5*m_ImageWidthY , m_ImageCenterY,
				(int)TMath::Power(2,j), m_ImageCenterZ-0.5*m_ImageWidthZ , m_ImageCenterZ+0.5*m_ImageWidthZ);
  }

  sprintf(name,"Tmp%d",0);
  Tmp[0] = new TH3F(name,name,
		    128, m_ImageCenterX-0.5*m_ImageWidthX , m_ImageCenterX, 
		    128, m_ImageCenterY-0.5*m_ImageWidthY , m_ImageCenterY,
		    128, m_ImageCenterZ-0.5*m_ImageWidthZ , m_ImageCenterZ+0.5*m_ImageWidthZ);
  


  for(int j=0;j<7;j++){
    sprintf(name,"TmpImageSpace%d%d",1,j);
    imagespace[1][j] = new TH3F(name,name,
				(int)TMath::Power(2,j), m_ImageCenterX, m_ImageCenterX+0.5*m_ImageWidthX, 
				(int)TMath::Power(2,j), m_ImageCenterY-0.5*m_ImageWidthY , m_ImageCenterY,
				(int)TMath::Power(2,j), m_ImageCenterZ-0.5*m_ImageWidthZ , m_ImageCenterZ+0.5*m_ImageWidthZ);
  }
  sprintf(name,"Tmp%d",1);
  Tmp[1] = new TH3F(name,name,
		    128, m_ImageCenterX, m_ImageCenterX+0.5*m_ImageWidthX, 
		    128, m_ImageCenterY-0.5*m_ImageWidthY , m_ImageCenterY,
		    128, m_ImageCenterZ-0.5*m_ImageWidthZ , m_ImageCenterZ+0.5*m_ImageWidthZ);
  
  for(int j=0;j<7;j++){
    sprintf(name,"TmpImageSpace%d%d",2,j);
    imagespace[2][j] = new TH3F(name,name,
				(int)TMath::Power(2,j), m_ImageCenterX, m_ImageCenterX+0.5*m_ImageWidthX, 
				(int)TMath::Power(2,j), m_ImageCenterY, m_ImageCenterY+0.5*m_ImageWidthY,
				(int)TMath::Power(2,j), m_ImageCenterZ-0.5*m_ImageWidthZ , m_ImageCenterZ+0.5*m_ImageWidthZ);
  }

  sprintf(name,"Tmp%d",2);
  Tmp[2] = new TH3F(name,name,
		    128, m_ImageCenterX, m_ImageCenterX+0.5*m_ImageWidthX, 
		    128, m_ImageCenterY, m_ImageCenterY+0.5*m_ImageWidthY,
		    128, m_ImageCenterZ-0.5*m_ImageWidthZ , m_ImageCenterZ+0.5*m_ImageWidthZ);


  for(int j=0;j<7;j++){
    sprintf(name,"TmpImageSpace%d%d",3,j);
    imagespace[3][j] = new TH3F(name,name,
				(int)TMath::Power(2,j), m_ImageCenterX-0.5*m_ImageWidthX , m_ImageCenterX, 
				(int)TMath::Power(2,j), m_ImageCenterY,m_ImageCenterY+0.5*m_ImageWidthY,
				(int)TMath::Power(2,j), m_ImageCenterZ-0.5*m_ImageWidthZ , m_ImageCenterZ+0.5*m_ImageWidthZ);
  }


  sprintf(name,"Tmp%d",3);
  Tmp[3] = new TH3F(name,name,
		    128, m_ImageCenterX-0.5*m_ImageWidthX , m_ImageCenterX, 
		    128, m_ImageCenterY,m_ImageCenterY+0.5*m_ImageWidthY,
		    128, m_ImageCenterZ-0.5*m_ImageWidthZ , m_ImageCenterZ+0.5*m_ImageWidthZ);
  


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
 

    for(int id=0;id<4;id++){

      pi.SetXYZ(imagespace[id][0]->GetXaxis()->GetBinCenter(1),imagespace[id][0]->GetYaxis()->GetBinCenter(1),imagespace[id][0]->GetZaxis()->GetBinCenter(1));
      diri = pi-ver;
      phii = axis.Angle(diri);
      chii = TMath::Abs(TMath::Pi()-thetaE-phii);
      mi = (diri.Mag()) * TMath::Sin(chii) - (diri.Mag()) * TMath::Cos(chii)* TMath::Tan(3.*sigma);
      rcrit = (binwidth[0] ) * 0.866;
      
      if(rcrit<mi){
	continue;
      }


      int nbinx= imagespace[id][1]->GetXaxis()->GetNbins();
      int nbiny= imagespace[id][1]->GetYaxis()->GetNbins();
      int nbinz= imagespace[id][1]->GetZaxis()->GetNbins();
      
      for(int i=1;i<=nbinx;i++){
	for(int j=1;j<=nbiny;j++){
	  for(int k=1;k<=nbinz;k++){
	    pi.SetXYZ(imagespace[id][1]->GetXaxis()->GetBinCenter(i),imagespace[id][1]->GetYaxis()->GetBinCenter(j),imagespace[id][1]->GetZaxis()->GetBinCenter(k));
	    diri = pi-ver;
	    phii = axis.Angle(diri);
	    chii = TMath::Abs(TMath::Pi()-thetaE-phii);
	    mi = (diri.Mag()) * TMath::Sin(chii) - (diri.Mag()) * TMath::Cos(chii)* TMath::Tan(3.*sigma);
	    rcrit = (binwidth[1]) * 0.866;

	    if(mi<rcrit){
	      for(int ii=2*(i-1)+1;ii<=2*i;ii++){
		for(int jj=2*(j-1)+1;jj<=2*j;jj++){
		  for(int kk=2*(k-1)+1;kk<=2*k;kk++){
		    pi.SetXYZ(imagespace[id][2]->GetXaxis()->GetBinCenter(ii),imagespace[id][2]->GetYaxis()->GetBinCenter(jj),imagespace[id][2]->GetZaxis()->GetBinCenter(kk));
		    diri = pi-ver;
		    phii = axis.Angle(diri);
		    chii = TMath::Abs(TMath::Pi()-thetaE-phii);
		    mi = (diri.Mag()) * TMath::Sin(chii) - (diri.Mag()) * TMath::Cos(chii)* TMath::Tan(3.*sigma);
		    rcrit = (binwidth[2] ) * 0.866;


		    if(mi<rcrit){
		      for(int iii=2*(ii-1)+1;iii<=2*ii;iii++){
			for(int jjj=2*(jj-1)+1;jjj<=2*jj;jjj++){
			  for(int kkk=2*(kk-1)+1;kkk<=2*kk;kkk++){
			    pi.SetXYZ(imagespace[id][3]->GetXaxis()->GetBinCenter(iii),imagespace[id][3]->GetYaxis()->GetBinCenter(jjj),imagespace[id][3]->GetZaxis()->GetBinCenter(kkk));
			    diri = pi-ver;
			    phii = axis.Angle(diri);
			    chii = TMath::Abs(TMath::Pi()-thetaE-phii);
			    mi = (diri.Mag()) * TMath::Sin(chii) - (diri.Mag()) * TMath::Cos(chii)* TMath::Tan(3.*sigma);
			    rcrit = (binwidth[3]) * 0.866;

			    if(mi<rcrit){
			      for(int iiii=2*(iii-1)+1;iiii<=2*iii;iiii++){
				for(int jjjj=2*(jjj-1)+1;jjjj<=2*jjj;jjjj++){
				  for(int kkkk=2*(kkk-1)+1;kkkk<=2*kkk;kkkk++){
				    pi.SetXYZ(imagespace[id][4]->GetXaxis()->GetBinCenter(iiii),imagespace[id][4]->GetYaxis()->GetBinCenter(jjjj),imagespace[id][4]->GetZaxis()->GetBinCenter(kkkk));
				    diri = pi-ver;
				    phii = axis.Angle(diri);
				    chii = TMath::Abs(TMath::Pi()-thetaE-phii);
				    mi = (diri.Mag()) * TMath::Sin(chii) - (diri.Mag()) * TMath::Cos(chii)* TMath::Tan(3.*sigma);
				    rcrit = (binwidth[4]) * 0.866;

				    if(mi<rcrit){
				      for(int iiiii=2*(iiii-1)+1;iiiii<=2*iiii;iiiii++){
					for(int jjjjj=2*(jjjj-1)+1;jjjjj<=2*jjjj;jjjjj++){
					  for(int kkkkk=2*(kkkk-1)+1;kkkkk<=2*kkkk;kkkkk++){
					    pi.SetXYZ(imagespace[id][5]->GetXaxis()->GetBinCenter(iiiii),imagespace[id][5]->GetYaxis()->GetBinCenter(jjjjj),imagespace[id][5]->GetZaxis()->GetBinCenter(kkkkk));
					    diri = pi-ver;
					    phii = axis.Angle(diri);
					    chii = TMath::Abs(TMath::Pi()-thetaE-phii);
					    mi = (diri.Mag()) * TMath::Sin(chii)- (diri.Mag()) * TMath::Cos(chii)* TMath::Tan(3.*sigma);
					    rcrit = (binwidth[5]) * 0.866;


					    if(mi<rcrit){
					      for(int iiiiii=2*(iiiii-1)+1;iiiiii<=2*iiiii;iiiiii++){
						for(int jjjjjj=2*(jjjjj-1)+1;jjjjjj<=2*jjjjj;jjjjjj++){
						  for(int kkkkkk=2*(kkkkk-1)+1;kkkkkk<=2*kkkkk;kkkkkk++){
						    pi.SetXYZ(imagespace[id][6]->GetXaxis()->GetBinCenter(iiiiii),imagespace[id][6]->GetYaxis()->GetBinCenter(jjjjjj),imagespace[id][6]->GetZaxis()->GetBinCenter(kkkkkk));
						    diri = pi-ver;
						    phii = axis.Angle(diri);
						    chii = TMath::Abs(TMath::Pi()-thetaE-phii);
						    mi = (diri.Mag()) * TMath::Sin(chii) - (diri.Mag()) * TMath::Cos(chii)* TMath::Tan(3.*sigma);
						    rcrit = (binwidth[6]) * 0.866;

						    if(mi<rcrit){
						      for(int iiiiiii=2*(iiiiii-1)+1;iiiiiii<=2*iiiiii;iiiiiii++){
							for(int jjjjjjj=2*(jjjjjj-1)+1;jjjjjjj<=2*jjjjjj;jjjjjjj++){
							  for(int kkkkkkk=2*(kkkkkk-1)+1;kkkkkkk<=2*kkkkkk;kkkkkkk++){
							    pi.SetXYZ(Tmp[id]->GetXaxis()->GetBinCenter(iiiiiii),Tmp[id]->GetYaxis()->GetBinCenter(jjjjjjj),Tmp[id]->GetZaxis()->GetBinCenter(kkkkkkk));
							    diri = pi-ver;
							    phii = axis.Angle(diri);
							    chii = TMath::Abs(TMath::Pi()-thetaE-phii);
							    mi = (diri.Mag()) * TMath::Sin(chii) - (diri.Mag()) * TMath::Cos(chii)* TMath::Tan(3.*sigma);
							    rcrit = (binwidth[7] ) * 0.866;
							    if(mi<rcrit){
							      int gbin = Tmp[id]->GetBin(iiiiiii,jjjjjjj,kkkkkkk);
							      double aa = (diri.Mag()) * TMath::Cos(chii)* TMath::Tan(sigma);

							      Tmp[id]->AddBinContent(gbin,1.0/(TMath::Power((diri.Mag()) * TMath::Cos(chii),weight))*TMath::Gaus((diri.Mag()) * TMath::Sin(chii),0,aa) );

							      
							    }else{
							      continue;
							    }
							  }
							}
						      }
						    }else{
						      continue;
						    }
						  }
						}
					      }
					    }else{
					      continue;
					    }
					  }
					}
				      }
				    }else{
				      continue;
				    }
				  }
				}
			      }
			      
			    }else{
			      continue;
			    }
			  }
			}
		      }
		    }else{
		      continue;
		    }
		  }
		}
	      }
	    }else{
	      continue;
	    }
	  }
	}
      }


    }
  }



  for(int i=1;i<=128;i++){
    for(int j=1;j<=128;j++){
      for(int k=1;k<=128;k++){
	double value = Tmp[0]->GetBinContent(i,j,k);
	int gbin = m_ImageArray[0]->GetBin(i,j,k);
	m_ImageArray[0]->SetBinContent(gbin,value);

	value = Tmp[1]->GetBinContent(i,j,k);
	gbin = m_ImageArray[0]->GetBin(i+128,j,k);
	m_ImageArray[0]->SetBinContent(gbin,value);

	value = Tmp[2]->GetBinContent(i,j,k);
	gbin = m_ImageArray[0]->GetBin(i+128,j+128,k);
	m_ImageArray[0]->SetBinContent(gbin,value);

	value = Tmp[3]->GetBinContent(i,j,k);
	gbin = m_ImageArray[0]->GetBin(i,j+128,k);
	m_ImageArray[0]->SetBinContent(gbin,value);
      }
    }
  }

  for(int i=0;i<4;i++){
    delete Tmp[i];
  }
  for(int i=0;i<7;i++){
    delete imagespace[0][i];
    delete imagespace[1][i];
    delete imagespace[2][i];
    delete imagespace[3][i];
  }

  return true;
}

/*
void SBPComptonImaging::SetResolutionFunc(double sigma, double ld)
{
  voigt_sigma = sigma;
  voigt_ld = ld;
  voigt_fwhm = calcvoigtfwhm(sigma,ld/2.);
  std::cout << "Voigt FWHM [degree] is " << voigt_fwhm << std::endl;
}
*/

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
  TH3F* responceTmp = new TH3F("resp", "resp", NPixelX, m_ImageCenterX-0.5*m_ImageWidthX, m_ImageCenterX+0.5*m_ImageWidthX, 
			       NPixelY, m_ImageCenterY-0.5*m_ImageWidthY, m_ImageCenterY+0.5*m_ImageWidthY,
			       NPixelZ, m_ImageCenterZ-0.5*m_ImageWidthZ, m_ImageCenterZ+0.5*m_ImageWidthZ);
  
  for (int i=0; i<nEvent; i++) {
    if (i%100==0) { std::cout << "    i = " << i << std::endl; }

    eventtree->GetEntry(i);
    TwoHitComptonEvent twohit;
    twohit.setH1Energy(edep1);
    twohit.setH1Pos(posx1,posy1,posz1);
    twohit.setH2Energy(edep2);
    twohit.setH2Pos(posx2,posy2,posz2);
    
    responceTmp->Reset();
    m_Response[i].MakeOneEventResponse(twohit,responceTmp,m_Efficiency,m_ImageArray[0],&m_Visibility[i]); 

  }
  
  for(int j=0;j<NUM_EVENT;j++){
    std::cout<< m_Visibility[j] << std::endl;
  }

  responceTmp->Delete();
  return true;
}

bool SBPComptonImaging::setEfficiencyMap(TH3F* h_eff)
{

  NPixelX = h_eff->GetXaxis()->GetNbins();
  NPixelY = h_eff->GetYaxis()->GetNbins();
  NPixelZ = h_eff->GetZaxis()->GetNbins();

  if(!(NPixelZ==128 && NPixelX==256 && NPixelY==256)){
    std::cout << "Voxel size is not 256x256x128 " << NPixelX << NPixelY << NPixelZ <<std::endl;
    return false;
  }

  m_Efficiency = h_eff;

  m_ImageWidthX = ( h_eff->GetXaxis()->GetXmax() ) - ( h_eff->GetXaxis()->GetXmin() ) ;
  m_ImageWidthY = ( h_eff->GetYaxis()->GetXmax() ) - ( h_eff->GetYaxis()->GetXmin() ) ;
  m_ImageWidthZ = ( h_eff->GetZaxis()->GetXmax() ) - ( h_eff->GetZaxis()->GetXmin() ) ;

  m_ImageCenterX = (( h_eff->GetXaxis()->GetXmax() ) + ( h_eff->GetXaxis()->GetXmin() ) ) /2. ;
  m_ImageCenterY = (( h_eff->GetYaxis()->GetXmax() ) + ( h_eff->GetYaxis()->GetXmin() ) ) /2. ;
  m_ImageCenterZ = (( h_eff->GetZaxis()->GetXmax() ) + ( h_eff->GetZaxis()->GetXmin() ) ) /2.;

  return true;
}

void SBPComptonImaging::SetInitImage(TH3F *init){

  for(int i = 1; i<=NPixelX;i++){
    for(int j = 1; j<=NPixelY;j++){
      for(int k = 1; k<=NPixelZ;k++){
	int gbin = init->GetBin(i,j,k);
	m_ImageArray[0] -> SetBinContent(gbin,init->GetBinContent(gbin));

      }
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



void SBPComptonImaging::improveImage(TH3F* newimage, const TH3F* orgimage, int ite)
{

  int j1, j2, j3;
  double factor=0.;

  for (j1=1; j1<=NPixelX; j1++) {
    for (j2=1; j2<=NPixelY; j2++) {
      for (j3=1; j3<=NPixelZ; j3++) {
	factor = 0.0;
	int gbin = orgimage->GetBin(j1,j2,j3);

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
