


void psf(){


  gROOT->Reset();

  int nbinx = 128;
  int nbiny = 128;

  TFile *fout = new TFile("psf.root","recreate");
  TH2D *image = new TH2D("psf","psf",nbinx,-102.4,102.4,nbiny,-102.4,102.4);
  //  FILE *fo;
  //  if((fo=fopen("psf.dat","wb"))==NULL){
  //    std::cout << "error " << std::endl;
  //  }

  //  double xcenter = 0.05;
  //  double ycenter =0.05;

  double xcenter = 0.0;
  double ycenter =0.0;
  double sum = 0.0;
  for(int i=1;i<=nbinx;i++){
    for(int j=1;j<=nbiny;j++){
      double a = image->GetXaxis()->GetBinCenter(i);
      double b = image->GetYaxis()->GetBinCenter(j);
      double r =  TMath::Sqrt(TMath::Power(xcenter-a,2)+TMath::Power(ycenter-b,2));

      double value =  TMath::Gaus(r,0,8.0,kTRUE);
      //      value = value * (1.+0.01*gRandom->Rndm());
      sum = sum + value;
      image->SetBinContent(i,j,value);

    }
  }

  for(int i=1;i<=nbinx;i++){
    for(int j=1;j<=nbiny;j++){
      image->SetBinContent(i,j,image->GetBinContent(i,j)/sum);
    }
  }

  image->Write();
  fout->Close();
  //  fclose(fo);

}
