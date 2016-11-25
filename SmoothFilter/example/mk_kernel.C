
void mk_kernel(){

  gROOT->Reset();

  int nbin = 11; // odd number

  double fwhm = 2.;
  std::cout << " FWHM : " ;
  std::cin >> fwhm;

  double sigma = fwhm/2.355;
  std::cout << sigma << std::endl;

  TFile *fout = new TFile("kernel.root","recreate");
  
  TH2D *kernel = new TH2D("kernel","kernel",nbin,0.5,nbin+0.5,nbin,0.5,nbin+0.5);

  TF1 *fg = new TF1("fg","(1./TMath::Sqrt(2.*TMath::Pi()*[1]*[1]))*TMath::Exp(-0.5 * ((x-[0])/[1])**2)",-nbin*10,nbin*10);
  fg->SetParameter(0,0);
  fg->SetParameter(1,sigma);


  
  std::cout << fg->Eval(0) << " " <<fg->Eval(sigma) << std::endl;

  double sum = 0.;
  for(int i=1;i<=nbin;i++){
    for(int j=1;j<=nbin;j++){
      double x = (nbin+1)/2. -i;
      double y = (nbin+1)/2. -j;
      double r = TMath::Sqrt(x*x + y*y);
      double value = fg->Eval(r);
      sum = sum + value;
      kernel->SetBinContent(i,j,value);
    }
  }


  for(int i=1;i<=nbin;i++){
    for(int j=1;j<=nbin;j++){
      kernel->SetBinContent(i,j,kernel->GetBinContent(i,j)/sum);
    }
  }

  
  std::cout << sum << std::endl;
  TCanvas *c1 = new TCanvas("c1","c1",0,0,500,500);
  c1->Divide(2,1);
  c1->cd(1);  
  kernel->Draw("colz");
  c1->cd(2);
  fg->Draw();

  fout->cd();
  kernel->Write();
  fout->Close();
}
