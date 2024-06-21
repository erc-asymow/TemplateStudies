void merge_massloop(TString tag = "SmearRealistic3Loops", TString name = "massloop_in_iter2", bool batch=false, bool savePng=false){

  //TString plotname = "massloop_in_iter2";
  //TString tag = "SmearRealistic3Loops";
  TString plotname = name + TString("_") + tag;
  
  TCanvas* c = new TCanvas("c", "canvas", 1600, 400);
  c->Divide(4,1);

  TFile* fout = TFile::Open("merge_massloop_"+tag+"_"+name+".root", "RECREATE");
  TTree* treeout = new TTree("tree","");
  
  TFile* fsIter0 = TFile::Open("massscales_"+tag+"_Iter0.root", "READ");
  TFile* fsIter1 = TFile::Open("massscales_"+tag+"_Iter1.root", "READ");
  TFile* fsIter2 = TFile::Open("massscales_"+tag+"_Iter2.root", "READ");
  TFile* fsIter3 = TFile::Open("massscales_"+tag+"_Iter3.root", "READ");
  
  TH1D* hmass_smear0 = 0;
  TH1D* hmass_smear1 = 0;
  if(string(plotname.Data()).find("iter0")!=string::npos){
    hmass_smear0 = ((TH2D*)fsIter0->Get("h_smear0_bin_m"))->ProjectionY("smear0");
    hmass_smear1 = ((TH2D*)fsIter0->Get("h_smear1_bin_m"))->ProjectionY("smear1");
  }
  else if(string(plotname.Data()).find("iter1")!=string::npos){
    hmass_smear0 = ((TH2D*)fsIter1->Get("h_smear0_bin_m"))->ProjectionY("smear0");
    hmass_smear1 = ((TH2D*)fsIter1->Get("h_smear1_bin_m"))->ProjectionY("smear1");
  }
  else if(string(plotname.Data()).find("iter2")!=string::npos){
    hmass_smear0 = ((TH2D*)fsIter2->Get("h_smear0_bin_m"))->ProjectionY("smear0");
    hmass_smear1 = ((TH2D*)fsIter2->Get("h_smear1_bin_m"))->ProjectionY("smear1");
  }
  else if(string(plotname.Data()).find("iter3")!=string::npos && fsIter3!=0){
    hmass_smear0 = ((TH2D*)fsIter3->Get("h_smear0_bin_m"))->ProjectionY("smear0");
    hmass_smear1 = ((TH2D*)fsIter3->Get("h_smear1_bin_m"))->ProjectionY("smear1");
  }
  hmass_smear1->Divide(hmass_smear0);
  hmass_smear1->SetLineColor(kBlack);
  hmass_smear1->SetMaximum(1.1);
  hmass_smear1->SetMinimum(0.9);
  c->cd(1);
  hmass_smear1->SetStats(0);
  hmass_smear1->SetTitle("data / nominal");
  hmass_smear1->Draw("HISTE");
  
  
  TFile* fIter0 = TFile::Open("massfit_"+tag+"_Iter0.root", "READ");
  TFile* fIter1 = TFile::Open("massfit_"+tag+"_Iter1.root", "READ");
  TFile* fIter2 = TFile::Open("massfit_"+tag+"_Iter2.root", "READ");
  TFile* fIter3 = TFile::Open("massfit_"+tag+"_Iter3.root", "READ");

  //vector<TString> params = {"Ain", "ein", "Min"};
  vector<TString> params = {"A", "e", "M"};
  if(string(name.Data()).find("in")!=string::npos)
    params = {"Ain", "ein", "Min"};
    
  TTree* t0 = (TTree*) fIter0->Get("tree");
  TTree* t1 = (TTree*) fIter1->Get("tree");
  TTree* t2 = (TTree*) fIter2->Get("tree");
  TTree* t3 = 0;
  if(fIter3!=0)
    t3 = (TTree*) fIter3->Get("tree");
  double fmin0, prob0;
  double fmin1, prob1;
  double fmin2, prob2;
  double fmin3, prob3;
  t0->SetBranchAddress("fmin", &fmin0);
  t0->SetBranchAddress("prob", &prob0);
  t1->SetBranchAddress("fmin", &fmin1);
  t1->SetBranchAddress("prob", &prob1);
  t2->SetBranchAddress("fmin", &fmin2);
  t2->SetBranchAddress("prob", &prob2);
  t0->GetEntry(0);
  t1->GetEntry(0);
  t2->GetEntry(0);
  if(fIter3!=0){
    t3->SetBranchAddress("fmin", &fmin3);
    t3->SetBranchAddress("prob", &prob3);
    t3->GetEntry(0);
  }

  TH1D* hp_nom_template  = (TH1D*) fIter0->Get("h_A_vals_nom");
  double A_val_fits[hp_nom_template->GetNbinsX()];
  double A_val_noms[hp_nom_template->GetNbinsX()];
  double A_val_errs[hp_nom_template->GetNbinsX()];
  double e_val_fits[hp_nom_template->GetNbinsX()];
  double e_val_noms[hp_nom_template->GetNbinsX()];
  double e_val_errs[hp_nom_template->GetNbinsX()];
  double M_val_fits[hp_nom_template->GetNbinsX()];
  double M_val_noms[hp_nom_template->GetNbinsX()];
  double M_val_errs[hp_nom_template->GetNbinsX()];

  for(unsigned int p = 0; p < params.size(); p++){
    TH1D* hp_nom  = (TH1D*) fIter0->Get("h_"+params[p]+"_vals_nom");
    TH1D* hp_fit0 = (TH1D*) fIter0->Get("h_"+params[p]+"_vals_fit");
    TH1D* hp_fit1 = (TH1D*) fIter1->Get("h_"+params[p]+"_vals_fit");
    TH1D* hp_fit2 = (TH1D*) fIter2->Get("h_"+params[p]+"_vals_fit");
    TH1D* hp_fit3 = 0;
    if(fIter3!=0)
      hp_fit3 = (TH1D*) fIter3->Get("h_"+params[p]+"_vals_fit");
    
    hp_fit0->SetTitle(Form("#chi^{2}/ndof = %.2f (prob=%.2f)", 1+fmin0, prob0));
    if( string(plotname.Data()).find("iter1")!=string::npos ){
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit0->SetBinError(ib, hp_fit1->GetBinError(ib) );
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit1->SetBinError(ib, 0.);
      hp_fit0->Add(hp_fit1);
      hp_fit0->SetTitle(Form("#chi^{2}/ndof = %.2f (prob=%.2f)", 1+fmin1, prob1));
    }
    else if( string(plotname.Data()).find("iter2")!=string::npos ){
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit0->SetBinError(ib, hp_fit2->GetBinError(ib) );
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit1->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit2->GetXaxis()->GetNbins();ib++) hp_fit2->SetBinError(ib, 0.);
      hp_fit0->Add(hp_fit1);
      hp_fit0->Add(hp_fit2);
      hp_fit0->SetTitle(Form("#chi^{2}/ndof = %.2f (prob=%.2f)", 1+fmin2, prob2));
    }
    else if( string(plotname.Data()).find("iter3")!=string::npos && fIter3!=0 ){
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit0->SetBinError(ib, hp_fit3->GetBinError(ib) );
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit1->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit2->GetXaxis()->GetNbins();ib++) hp_fit2->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit2->GetXaxis()->GetNbins();ib++) hp_fit3->SetBinError(ib, 0.);
      hp_fit0->Add(hp_fit1);
      hp_fit0->Add(hp_fit2);
      hp_fit0->Add(hp_fit3);
      hp_fit0->SetTitle(Form("#chi^{2}/ndof = %.2f (prob=%.2f)", 1+fmin3, prob3));
    }
    hp_nom->SetLineColor(kBlue);
    hp_nom->SetLineWidth(3);
    hp_fit0->SetLineColor(kBlack);
    hp_fit0->SetMarkerColor(kBlack);
    hp_fit0->SetMarkerStyle(kFullCircle);
    c->cd(p+2);
    hp_nom->SetStats(0);
    //hp_nom->SetMaximum( hp_nom->GetMaximum()*1.5);
    //hp_nom->SetMinimum( hp_nom->GetMinimum()*0.5);
    hp_fit0->SetStats(0);
    //hp_fit0->SetMaximum( hp_fit0->GetMaximum()*1.5);
    //hp_fit0->SetMinimum( hp_fit0->GetMinimum()*0.5);

    //if(p>=1){
    //  hp_fit0->SetMaximum( + (hp_fit0->GetMaximum()+hp_fit0->GetBinError(1))*1.1 );
    //  hp_fit0->SetMinimum( - (hp_fit0->GetMaximum()+hp_fit0->GetBinError(1))*1.1 );
    //}

    if(p==0 && savePng){
      hp_fit0->SetMaximum( +0.002 );
      hp_fit0->SetMinimum( -0.002 );
    }
    if(p==1 && savePng){
      hp_fit0->SetMaximum( +0.0025 );
      hp_fit0->SetMinimum( -0.0025 );
    }
    if(p==2 && savePng){
      hp_fit0->SetMaximum( +0.0015 );
      hp_fit0->SetMinimum( -0.0015 );
    }          
    hp_fit0->Draw("HISTPE");
    hp_nom->Draw("HISTSAME");

    for(int ib = 0; ib<hp_nom->GetNbinsX(); ib++){
      TString b_nom_name(Form("%c%d%s", params[p][0], ib, string(params[p].Data()).find("in")!=string::npos ? Form("_intrue") : Form("_true")   ));
      TString b_fit_name(Form("%c%d%s", params[p][0], ib, string(params[p].Data()).find("in")!=string::npos ? Form("_in") : Form("")   ));
      TString b_err_name(Form("%c%d%s", params[p][0], ib, string(params[p].Data()).find("in")!=string::npos ? Form("_inerr") : Form("_err")   )); 
      if(p==0){
	treeout->Branch( b_nom_name.Data(), &(A_val_noms[ib]), (b_nom_name+TString("/D")).Data()   );
	treeout->Branch( b_fit_name.Data(), &(A_val_fits[ib]), (b_fit_name+TString("/D")).Data()   );
	treeout->Branch( b_err_name.Data(), &(A_val_errs[ib]), (b_err_name+TString("/D")).Data()   );
      }
      else if(p==1){
	treeout->Branch( b_nom_name.Data(), &(e_val_noms[ib]), (b_nom_name+TString("/D")).Data()   );
	treeout->Branch( b_fit_name.Data(), &(e_val_fits[ib]), (b_fit_name+TString("/D")).Data()   );
	treeout->Branch( b_err_name.Data(), &(e_val_errs[ib]), (b_err_name+TString("/D")).Data()   );
      }
      else if(p==2){
	treeout->Branch( b_nom_name.Data(), &(M_val_noms[ib]), (b_nom_name+TString("/D")).Data()   );
	treeout->Branch( b_fit_name.Data(), &(M_val_fits[ib]), (b_fit_name+TString("/D")).Data()   );
	treeout->Branch( b_err_name.Data(), &(M_val_errs[ib]), (b_err_name+TString("/D")).Data()   );
      }
    }
    for(int ib = 0; ib<hp_nom->GetNbinsX(); ib++){
      if(p==0){
	A_val_noms[ib] = hp_nom->GetBinContent(ib+1);
	A_val_fits[ib] = hp_fit0->GetBinContent(ib+1);
	A_val_errs[ib] = hp_fit0->GetBinError(ib+1);
      }
      else if(p==1){
	e_val_noms[ib] = hp_nom->GetBinContent(ib+1);
	e_val_fits[ib] = hp_fit0->GetBinContent(ib+1);
	e_val_errs[ib] = hp_fit0->GetBinError(ib+1);
      }
      else if(p==2){
	M_val_noms[ib] = hp_nom->GetBinContent(ib+1);
	M_val_fits[ib] = hp_fit0->GetBinContent(ib+1);
	M_val_errs[ib] = hp_fit0->GetBinError(ib+1);
      }
    }
    fout->cd();
    hp_nom->Write(Form("h_%c_vals_nom", params[p][0]));
  }

  c->Update();
  c->Draw();
  if(savePng) c->SaveAs(plotname+".png");
  if(batch){
    fIter0->Close();
    fIter1->Close();
    fIter2->Close();
    delete c;
  }
  fout->cd();
  treeout->Fill();
  treeout->Write();
  fout->Close();
}


void hadd_massloop_toys(TString name = "massloop_in"){
  vector<TString> iters = {"iter0", "iter1", "iter2", "iter3"};
  for(unsigned int it=0; it<iters.size(); it++){
    cout << "Doing iter " << iters[it] << endl;
    for(int itoy=0; itoy<100; itoy++){
      cout << "Doing toy " << itoy << endl;
      merge_massloop( TString(Form("SmearRealisticRnd_toy%d", itoy)),
		      name+"_"+iters[it], true, false);
    }
    gSystem->Exec( Form("hadd -f merge_massloop_SmearRealisticRnd_merged_%s_%s.root merge_massloop_SmearRealisticRnd_toy*_%s_%s.root ; rm merge_massloop_SmearRealisticRnd_toy*_%s_%s.root", name.Data(), iters[it].Data(), name.Data(), iters[it].Data(), name.Data(), iters[it].Data() ));
  }
}

void merge_massloop_toys(TString tag = "SmearRealisticRnd_merged", TString name = "massloop_in",
			 bool batch=false,
			 bool savePng=false,
			 bool plotMean=false){
  
  TString plotname = name + TString("_") + tag + TString(plotMean ? "_mean" : "_sigma");

  TH1D* h_pIter0 = 0;
  TH1D* h_pIter1 = 0;
  TH1D* h_pIter2 = 0;
  TH1D* h_pIter3 = 0;

  TCanvas* c = new TCanvas("c", "canvas", 1600, 400);
  c->Divide(3,1);

  TLegend* leg1 = new TLegend(0.45, 0.65, 0.85, 0.90, "","brNDC");
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.045);
  leg1->SetFillColor(10);
  leg1->SetHeader("");

  //TFile* fIter0 = TFile::Open("massfit_"+tag+"_Iter0.root", "READ");
  //TFile* fIter1 = TFile::Open("massfit_"+tag+"_Iter1.root", "READ");
  //TFile* fIter2 = TFile::Open("massfit_"+tag+"_Iter2.root", "READ");
  //TFile* fIter3 = TFile::Open("massfit_"+tag+"_Iter3.root", "READ");

  TFile* fIter0 = TFile::Open("merge_massloop_"+tag+"_"+name+"_iter0.root", "READ");
  TFile* fIter1 = TFile::Open("merge_massloop_"+tag+"_"+name+"_iter1.root", "READ");
  TFile* fIter2 = TFile::Open("merge_massloop_"+tag+"_"+name+"_iter2.root", "READ");
  TFile* fIter3 = TFile::Open("merge_massloop_"+tag+"_"+name+"_iter3.root", "READ");
  
    
  vector<TString> params = {"A", "e", "M"};
  
  for(unsigned int p = 0; p < params.size(); p++){


    TH1D* h_template = (TH1D*)fIter0->Get(Form("h_%s_vals_nom", params[p].Data()));
    int nparams = h_template->GetNbinsX();
    c->cd(p+1);
    h_pIter0 = new TH1D(Form("h_%s_Iter0", params[p].Data()), "", nparams, 0, nparams);   
    h_pIter1 = new TH1D(Form("h_%s_Iter1", params[p].Data()), "", nparams, 0, nparams);
    h_pIter2 = new TH1D(Form("h_%s_Iter2", params[p].Data()), "", nparams, 0, nparams);   
    if(fIter3!=0)
      h_pIter3 = new TH1D(Form("h_%s_Iter3", params[p].Data()), "", nparams, 0, nparams);         
      
    for(int ip = 0; ip<nparams; ip++){
      h_pIter0->GetXaxis()->SetBinLabel(ip+1, h_template->GetXaxis()->GetBinLabel(ip+1));
      h_pIter1->GetXaxis()->SetBinLabel(ip+1, h_template->GetXaxis()->GetBinLabel(ip+1));
      h_pIter2->GetXaxis()->SetBinLabel(ip+1, h_template->GetXaxis()->GetBinLabel(ip+1));
      if(fIter3!=0)
	h_pIter3->GetXaxis()->SetBinLabel(ip+1, h_template->GetXaxis()->GetBinLabel(ip+1));
    }
    
    TTree* t0 = (TTree*) fIter0->Get("tree");
    TTree* t1 = (TTree*) fIter1->Get("tree");
    TTree* t2 = (TTree*) fIter2->Get("tree");
    TTree* t3 = 0;
    if(fIter3!=0)
      t3 = (TTree*) fIter3->Get("tree");

    TH1D* hpulls = new TH1D(Form("hpulls%d",p), "", 40, -8, 8);
    for(int ip = 0; ip<nparams; ip++){
      
      TF1* gf = new TF1("gf","[0]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp( -0.5*(x-[1])*(x-[1])/[2]/[2] )", -8, 8);
      gf->SetNpx(10000);
      gf->SetParameter(0, 100.);
      gf->SetParameter(1, 0.);
      gf->SetParameter(2, 1.0 );      

      TString formula(Form("(%s%d_in - %s%d_intrue)/%s%d_inerr>>hpulls%d", params[p].Data(), ip, params[p].Data(), ip, params[p].Data(), ip, p));
      if( string(name.Data()).find("in")==string::npos )
	formula = TString(Form("(%s%d - %s%d_true)/%s%d_err>>hpulls%d", params[p].Data(), ip, params[p].Data(), ip, params[p].Data(), ip, p));
      float mean = 0.;
      float mean_err = 0.;
      float rms = 0.;
      float rms_err = 0.;
      t0->Draw(formula.Data(), "", "");
      hpulls->Fit("gf", "QR", "", -6, 6);
      mean = gf->GetParameter(1);
      mean_err = gf->GetParError(1);
      rms = TMath::Abs(gf->GetParameter(2));
      rms_err = gf->GetParError(2);
      //rms = hpulls->GetRMS();
      //rms_err = hpulls->GetRMSError();
      h_pIter0->SetBinContent(ip+1, plotMean ? mean : rms);
      h_pIter0->SetBinError(ip+1, plotMean ? mean_err : rms_err);
      hpulls->Reset();
	
      t1->Draw(formula.Data(), "", "");
      gf->SetParameter(0, 100.);
      gf->SetParameter(1, 0.);
      gf->SetParameter(2, 1.0 );      
      hpulls->Fit("gf", "QR", "", -3, 3);
      mean = gf->GetParameter(1);
      mean_err = gf->GetParError(1);
      rms = TMath::Abs(gf->GetParameter(2));
      rms_err = gf->GetParError(2);
      //rms = hpulls->GetRMS();
      //rms_err = hpulls->GetRMSError();
      h_pIter1->SetBinContent(ip+1, plotMean ? mean : rms);
      h_pIter1->SetBinError(ip+1, plotMean ? mean_err : rms_err);
      hpulls->Reset();

      t2->Draw(formula.Data(), "", "");
      gf->SetParameter(0, 100.);
      gf->SetParameter(1, 0.);
      gf->SetParameter(2, 1.0 );      
      hpulls->Fit("gf", "QR", "", -3, 3);
      mean = gf->GetParameter(1);
      mean_err = gf->GetParError(1);
      rms = TMath::Abs(gf->GetParameter(2));
      rms_err = gf->GetParError(2);
      //rms = hpulls->GetRMS();
      //rms_err = hpulls->GetRMSError();
      h_pIter2->SetBinContent(ip+1, plotMean ? mean : rms);
      h_pIter2->SetBinError(ip+1, plotMean ? mean_err : rms_err);
      hpulls->Reset();

      if(fIter3!=0){
	t3->Draw(formula.Data(), "", "");
	gf->SetParameter(0, 100.);
	gf->SetParameter(1, 0.);
	gf->SetParameter(2, 1.0 );      
	hpulls->Fit("gf", "QR", "", -3, 3);
	mean = gf->GetParameter(1);
	mean_err = gf->GetParError(1);
	rms = TMath::Abs(gf->GetParameter(2));
	rms_err = gf->GetParError(2);
	//rms = hpulls->GetRMS();
	//rms_err = hpulls->GetRMSError();
	h_pIter3->SetBinContent(ip+1, plotMean ? mean : rms);
	h_pIter3->SetBinError(ip+1, plotMean ? mean_err : rms_err);
	hpulls->Reset();
      }
      delete gf;
    }

    c->cd(p+1);
    h_pIter0->SetLineColor(kBlue);
    h_pIter1->SetLineColor(kGreen);
    h_pIter2->SetLineColor(kRed);
    h_pIter0->SetMarkerStyle(kFullSquare);
    h_pIter1->SetMarkerStyle(kFullTriangleUp);
    h_pIter2->SetMarkerStyle(kFullCircle);
    h_pIter0->SetMarkerColor(kBlue);
    h_pIter1->SetMarkerColor(kGreen);
    h_pIter2->SetMarkerColor(kRed);
    h_pIter0->SetMarkerSize(1.3);
    h_pIter1->SetMarkerSize(1.3);
    h_pIter2->SetMarkerSize(1.3);
    if(fIter3!=0){
      h_pIter3->SetLineColor(kMagenta);
      h_pIter3->SetMarkerStyle(kOpenCircle);
      h_pIter3->SetMarkerColor(kMagenta);
      h_pIter3->SetMarkerSize(1.3);
    }
    
    h_pIter0->SetTitle(Form("%s of pulls for param %s (%s)",
			    plotMean ? Form("Mean") : Form("#sigma"),
			    params[p].Data(),
			    string(name.Data()).find("in")==string::npos ? Form("out") : Form("in") ));    
    h_pIter0->SetStats(0);
    h_pIter0->SetMinimum(plotMean ? -2.0 : 0.25);
    h_pIter0->SetMaximum(plotMean ? +2.0 : 3.0);
    h_pIter0->Draw("E");
    h_pIter1->Draw("ESAME");
    h_pIter2->Draw("ESAME");
    if(fIter3!=0)
      h_pIter3->Draw("ESAME");
    
    if(p==0){
      leg1->AddEntry(h_pIter0, "Iter 0", "PL");
      leg1->AddEntry(h_pIter1, "Iter 1", "PL");
      leg1->AddEntry(h_pIter2, "Iter 2", "PL");
      if(fIter3!=0)
	leg1->AddEntry(h_pIter3, "Iter 3", "PL");
    }
    leg1->Draw();
    TF1* line = new TF1("line", plotMean ? "0.0" : "1.0", 0, nparams);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->Draw("same");
  }
    
  c->Update();
  c->Draw();
  if(savePng) c->SaveAs(plotname+".png");
  if(batch){
    fIter0->Close();
    fIter1->Close();
    fIter2->Close();
    delete c;
  }
}

