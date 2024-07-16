#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "RooRealVar.h"
#include "RooDerivative.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooCrystalBall.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"

using namespace RooFit;

int fit(int ib=-1, float minFrac=0.995, bool savePlots=false){


  if(ib<-1){
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  }
    
  TFile* fout = TFile::Open("cb.root", "RECREATE");
  TTree* tree = new TTree("tree", "");
  float edm, fr, norm;
  int status, bin, flag;
  tree->Branch("edm", &edm, "edm/F");
  tree->Branch("norm", &norm, "norm/F");
  tree->Branch("fr", &fr, "fr/F");
  tree->Branch("status", &status, "status/I");
  tree->Branch("bin", &bin, "bin/I");
  tree->Branch("bin", &flag, "flag/I");
  
  //int ibin = ib;
  

  TFile* f = TFile::Open("massscales_PostVFPMoreBins_Iter0.root");
  TH2D* h2 = (TH2D*)f->Get("h_reco_bin_dm");
  TH2D* h2m = (TH2D*)f->Get("h_reco_bin_m");

  for(int ibin=0; ibin<h2->GetXaxis()->GetNbins(); ibin++){

    if(ib>=0 && ibin!=ib) continue;

    if(ibin%500==0) cout << "Doing bin " << ibin << endl;
    
    TH1D* h = (TH1D*)h2->ProjectionY("h",ibin,ibin);
    TH1D* hm = (TH1D*)h2m->ProjectionY("hm",ibin,ibin);

    edm = -99.;
    status = -99;
    fr = -99.;
    flag = -99;
    bin = ibin;
    norm = h->Integral();

    if(norm<1.){
      tree->Fill();
      continue;
    }          
    
    RooRealVar mass("mass", Form("mass for bin %d", ibin), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    mass.setRange("r1", -10.0, 10.0);
    
    RooDataHist data("data", "", RooArgList(mass), h );
    
    RooRealVar x0("x0", "", h->GetMean(), mass.getMin(), mass.getMax() );
    RooRealVar sigmaL("sigmaL", "", h->GetRMS(), h->GetRMS()*0.5, h->GetRMS()*2 );
    RooRealVar sigmaR("sigmaR", "", h->GetRMS(), h->GetRMS()*0.5, h->GetRMS()*2 );
    RooRealVar alphaL("alphaL", "", 1.0, 0.2, +10 );
    RooRealVar alphaR("alphaR", "", 1.0, 0.2, +10 );
    RooRealVar nL("nL", "", 2, 1, 100 );
    RooRealVar nR("nR", "", 2, 1, 100 );
    
    RooCrystalBall pdf("pdf", "", mass, x0, sigmaL, sigmaR, alphaL, nL, alphaR, nR);

    RooGaussian gaus("pdf", "", mass, x0, sigmaL);
    //RooCBShape pdf("pdf", "", mass, x0, sigmaL, alphaL, nL);
    
    RooRealVar tau("tau", "", 0., -10, 10);
    RooExponential bkg("bkg", "", mass, tau);
    
    RooRealVar frac("frac", "", 0.9, 0., 1.);
    
    RooAddPdf pdfTot("pdfTot", "", {pdf, bkg}, frac);

    TString rname = "r1";
    pdfTot.fixCoefRange( rname.Data() );
    pdfTot.fixCoefNormalization(mass);
    
    std::shared_ptr<RooFitResult> r{pdfTot.fitTo(data,
						 InitialHesse(true),
						 Minimizer("Minuit2"),
						 Range( rname.Data() ),
						 Save(), SumW2Error(true),
						 PrintLevel(ib>=0 ? 1 : -1),
						 Verbose(ib>=0) )};
    edm = r->edm();
    status = r->status();
    fr = frac.getVal();
    flag = 0;

    bool fallBack_CB = false;
    bool fallBack_Gaus = false;
    if(fr>minFrac){
      fallBack_CB = true;
      std::shared_ptr<RooFitResult> rn{pdf.fitTo(data,
						 InitialHesse(true),
						 Minimizer("Minuit2"),
						 Range( rname.Data() ),
						 Save(), SumW2Error(true),
						 PrintLevel(ib>=0 ? 1 : -1),
						 Verbose(ib>=0) )};      
      r = rn;
      flag = 1;
    }
    edm = r->edm();
    status = r->status();

    if(status!=0){
      rname = "r2";
      mass.setRange( rname.Data(), -8.0, 8.0);
      std::shared_ptr<RooFitResult> rn{pdf.fitTo(data,
						InitialHesse(true),
						Minimizer("Minuit2"),
						Range( rname.Data() ),
						Save(), SumW2Error(true),
						PrintLevel(ib>=0 ? 1 : -1),
						Verbose(ib>=0) )};      
      r = rn;
      flag = 2;
    }
    edm = r->edm();
    status = r->status();

    if(status!=0){
      rname = "r3";
      mass.setRange( rname.Data() , -6.0, 6.0);
      std::shared_ptr<RooFitResult> rn{pdf.fitTo(data,
						InitialHesse(true),
						Minimizer("Minuit2"),
						Range( rname.Data() ),
						Save(), SumW2Error(true),
						PrintLevel(ib>=0 ? 1 : -1),
						Verbose(ib>=0) )};      
      r = rn;
      flag = 3;
    }
    edm = r->edm();
    status = r->status();
    
    if(status!=0){
      fallBack_Gaus = true;
      rname = "r3";
      std::shared_ptr<RooFitResult> rn{gaus.fitTo(data,
						 InitialHesse(true),
						 Minimizer("Minuit2"),
						 Range( rname.Data() ),
						 Save(), SumW2Error(true),
						 PrintLevel(ib>=0 ? 1 : -1),
						 Verbose(ib>=0) )};      
      r = rn;
      flag = 4;
    }
    edm = r->edm();
    status = r->status();
    
    tree->Fill();
    
    if(ib>=0 || savePlots){
      if(ib>=0) r->Print();      
      RooPlot* frame = mass.frame();
      data.plotOn(frame);
      if(!fallBack_CB && !fallBack_Gaus){
	pdfTot.plotOn(frame, VisualizeError(*r), Range( rname.Data() ));
	pdfTot.plotOn(frame, Components(pdf), LineColor(kGreen), Range( rname.Data() ));
	pdfTot.plotOn(frame, Components(bkg), LineColor(kRed), Range( rname.Data() ));
	pdfTot.plotOn(frame, LineColor(kBlue), Range( rname.Data() ));
      }
      else if(fallBack_CB){
	pdf.plotOn(frame, VisualizeError(*r), Range( rname.Data() ));
      }
      else if(fallBack_Gaus){
	gaus.plotOn(frame, VisualizeError(*r), Range( rname.Data() ));
      }
      data.plotOn(frame);
      
      TCanvas* c = new TCanvas("c", "canvas", 1200, 400);
      c->Divide(3,1);
      
      TH1D* h_der = new TH1D("h_der", "", h->GetXaxis()->GetNbins()*2, h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax() );
      h_der->Reset();
      h_der->SetStats(0);
      h_der->SetTitle("derivative");
      TH1D* h_jscale = (TH1D*)h_der->Clone("h_jscale");
      h_jscale->SetStats(0);
      h_jscale->SetTitle("scale Jacobian");
      TH1D* h_jwidth = (TH1D*)h_der->Clone("h_jwidth");  
      h_jwidth->SetStats(0);
      h_jwidth->SetTitle("resolution Jacobian");
      RooDerivative* der = 0;
      if(!fallBack_CB && !fallBack_Gaus)
	der = pdfTot.derivative(mass, 1, 0.001 );
      else if(fallBack_CB)
	der = pdf.derivative(mass, 1, 0.001 );
      else if(fallBack_Gaus)
	der = gaus.derivative(mass, 1, 0.001 );
	
      for(int ib=1; ib<=h_der->GetXaxis()->GetNbins();ib++){
	double x = h_der->GetXaxis()->GetBinCenter(ib);
	mass.setVal( x );
	//RooDerivative* der = pdf.derivative(mass, 1, 0.001 ); 
	double fprime = der->getVal();
	double f = 0.;
	if(!fallBack_CB && !fallBack_Gaus)
	  f = pdfTot.getVal();
	else if(fallBack_CB)
	  f = pdf.getVal();
	else if(fallBack_Gaus)
	  f = gaus.getVal();	
	h_der->SetBinContent(ib, fprime);
	h_jscale->SetBinContent(ib, -fprime/f * hm->GetMean());   
	h_jwidth->SetBinContent(ib, -(1+x*fprime/f) );
      }
      //h_jscale->Print("all");
      c->cd(1);
      frame->Draw();
      //c->cd(2);
      //h_der->Draw();
      c->cd(2);
      h_jscale->Draw();
      c->cd(3);
      h_jwidth->Draw();
      c->Update();
      c->Draw();
      if(savePlots) c->SaveAs(Form("plots/cb/cbfit_bin%d.png", ibin));
    }
  }

  fout->cd();
  tree->Write();
  fout->Close();
  
  return 0;

}
