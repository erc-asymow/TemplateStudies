{

  int ibin = 139;
  
  using namespace RooFit;

  TFile* f = TFile::Open("massscales_PostVFPMoreBins_Iter0.root");
  TH2D* h2 = (TH2D*)f->Get("h_reco_bin_dm");
  TH2D* h2m = (TH2D*)f->Get("h_reco_bin_m");
  TH1D* h = (TH1D*)h2->ProjectionY("h",ibin,ibin);
  TH1D* hm = (TH1D*)h2m->ProjectionY("hm",ibin,ibin);

  RooRealVar mass("mass", "", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  mass.setRange("r1", -8.0, 8.0);
  
  RooDataHist data("data", "", RooArgList(mass), h );

  RooRealVar x0("x0", "", h->GetMean(), mass.getMin(), mass.getMax() );
  RooRealVar sigmaL("sigmaL", "", h->GetRMS(), h->GetRMS()*0.5, h->GetRMS()*2 );
  RooRealVar sigmaR("sigmaR", "", h->GetRMS(), h->GetRMS()*0.5, h->GetRMS()*2 );
  RooRealVar alphaL("alphaL", "", 1.0, 0.2, +10 );
  RooRealVar alphaR("alphaR", "", 1.0, 0.2, +10 );
  RooRealVar nL("nL", "", 2, 1, 100 );
  RooRealVar nR("nR", "", 2, 1, 100 );
  
  RooCrystalBall pdf("pdf", "", mass, x0, sigmaL, sigmaR, alphaL, nL, alphaR, nR);
  //RooGaussian pdf("pdf", "", mass, x0, sigmaL);
  //RooCBShape pdf("pdf", "", mass, x0, sigmaL, alphaL, nL);
  
  std::unique_ptr<RooFitResult> r{pdf.fitTo(data,
					    InitialHesse(true),
					    Minimizer("Minuit2"),
					    Range("r1"),
					    Save(), SumW2Error(true) )};
  r->Print();

  RooPlot* frame = mass.frame();
  data.plotOn(frame);
  pdf.plotOn(frame, VisualizeError(*r), Range("r1"));
  data.plotOn(frame);
  frame->Draw();  

  
  TCanvas* c = new TCanvas("c", "canvas", 1600, 400);
  c->Divide(3,1);

  TH1D* h_der = new TH1D("h_der", "", h->GetXaxis()->GetNbins()*2, h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax() );
  h_der->Reset();
  TH1D* h_jscale = (TH1D*)h_der->Clone("h_jscale");
  TH1D* h_jwidth = (TH1D*)h_der->Clone("h_jwidth");  
  RooDerivative* der = pdf.derivative(mass, 1, 0.001 ); 
  for(int ib=1; ib<=h_der->GetXaxis()->GetNbins();ib++){
    double x = h_der->GetXaxis()->GetBinCenter(ib);
    mass.setVal( x );
    //RooDerivative* der = pdf.derivative(mass, 1, 0.001 ); 
    double fprime = der->getVal();
    double f = pdf.getVal();
    h_der->SetBinContent(ib, fprime);
    h_jscale->SetBinContent(ib, -fprime/f * hm->GetMean());   
    h_jwidth->SetBinContent(ib, -(1+x*fprime/f) );
  }
  h_jscale->Print("all");
  c->cd(1);
  h_der->Draw();
  c->cd(2);
  h_jscale->Draw();
  c->cd(3);
  h_jwidth->Draw();
  c->Update();
  c->Draw();

}
