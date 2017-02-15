


void input(){


  gROOT->Reset();

  int nbinx = 128;
  int nbiny = 128;

  TFile *fout = new TFile("input.root","recreate");
  TH2D *image = new TH2D("image","image",nbinx,-102.4,102.4,nbiny,-102.4,102.4);
  //  FILE *fo;
  //  if((fo=fopen("psf.dat","wb"))==NULL){
  //    std::cout << "error " << std::endl;
  //  }

  //  double xcenter = 0.05;
  //  double ycenter =0.05;

  double xcenter = 30.0;
  double ycenter = 30.0;

  for(int i=1;i<=nbinx;i++){
    for(int j=1;j<=nbiny;j++){
      double a = image->GetXaxis()->GetBinCenter(i);
      double b = image->GetYaxis()->GetBinCenter(j);
      double r =  TMath::Sqrt(TMath::Power(xcenter-a,2)+TMath::Power(ycenter-b,2));

      double value =  TMath::Gaus(r,0,5.0);
      //value = value * (1.+0.01*gRandom->Rndm());
  
      image->SetBinContent(i,j,value);

    }
  }


  image->Write();
  fout->Close();
  //  fclose(fo);

}
