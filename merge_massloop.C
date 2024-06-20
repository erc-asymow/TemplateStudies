void merge_massloop(TString tag = "SmearRealistic3Loops", TString name = "massloop_in_iter2", bool batch=false, bool savePng=false){

  //TString plotname = "massloop_in_iter2";
  //TString tag = "SmearRealistic3Loops";
  TString plotname = name + TString("_") + tag;
  
  TCanvas* c = new TCanvas("c", "canvas", 1600, 400);
  c->Divide(4,1);

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

  vector<TString> params = {"Ain", "ein", "Min"};
  //vector<TString> params = {"A", "e", "M"};

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
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit1->SetBinError(ib, 0.);
      hp_fit0->Add(hp_fit1);
      hp_fit0->SetTitle(Form("#chi^{2}/ndof = %.2f (prob=%.2f)", 1+fmin1, prob1));
    }
    else if( string(plotname.Data()).find("iter2")!=string::npos ){
      for(unsigned int ib=1; ib<=hp_fit1->GetXaxis()->GetNbins();ib++) hp_fit1->SetBinError(ib, 0.);
      for(unsigned int ib=1; ib<=hp_fit2->GetXaxis()->GetNbins();ib++) hp_fit2->SetBinError(ib, 0.);
      hp_fit0->Add(hp_fit1);
      hp_fit0->Add(hp_fit2);
      hp_fit0->SetTitle(Form("#chi^{2}/ndof = %.2f (prob=%.2f)", 1+fmin2, prob2));
    }
    else if( string(plotname.Data()).find("iter3")!=string::npos && fIter3!=0 ){
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

    if(p==0){
      hp_fit0->SetMaximum( +0.002 );
      hp_fit0->SetMinimum( -0.002 );
    }
    if(p==1){
      hp_fit0->SetMaximum( +0.0025 );
      hp_fit0->SetMinimum( -0.0025 );
    }
    if(p==2){
      hp_fit0->SetMaximum( +0.0015 );
      hp_fit0->SetMinimum( -0.0015 );
    }          
    hp_fit0->Draw("HISTPE");
    hp_nom->Draw("HISTSAME");
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
