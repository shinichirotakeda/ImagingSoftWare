


void psf(){


  gROOT->Reset();

  int nbinx = 512;
  int nbiny = 512;

  TFile *fout = new TFile("psf.root","recreate");
  TH2F *image = new TH2F("psf","psf",nbinx,-250.0,250.0,nbiny,-250.0,250.0);
  //  FILE *fo;
  //  if((fo=fopen("psf.dat","wb"))==NULL){
  //    std::cout << "error " << std::endl;
  //  }

  //  float xcenter = 0.05;
  //  float ycenter =0.05;

  float xcenter = 0.0;
  float ycenter =0.0;
  float sum = 0.0;
  for(int i=1;i<=nbinx;i++){
    for(int j=1;j<=nbiny;j++){
      float a = image->GetXaxis()->GetBinCenter(i);
      float b = image->GetYaxis()->GetBinCenter(j);
      float r =  TMath::Sqrt(TMath::Power(xcenter-a,2)+TMath::Power(ycenter-b,2));

      float value =  TMath::Gaus(r,0,8.0,kTRUE);
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
