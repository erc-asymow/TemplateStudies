{

  bool logy = true;
  float norm_lumi = (0.100/0.785);
    
  TCanvas* c = new TCanvas("c", "canvas", 1000, 600);
  c->SetGridy();

  TMultiGraph* mg = new TMultiGraph();
  TLegend* leg1 = new TLegend(0.10, 0.62, 0.70, 0.90, "","brNDC");
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.025);
  leg1->SetFillColor(10);
  leg1->SetNColumns(2);
  leg1->SetHeader("Norm. = 1.0#times10^{8} events, p_{T}/GeV#in[25,60], |#eta|<2.5, no A_{4}");

  TString image_name = "deltaM_vs_N_paper_noA4.eps";
  
  vector<TString> posttags = {
    "jUL_j0_j1_j2_j3_DEBUG_ADDMC_ULA0A1A2A3",
    "jUL_j0_j1_j2_j3_DEBUG_ULA0A1A2A3",
  };

    
  ////////////////////////////////////////////////
  for(unsigned int i = 0; i<posttags.size(); i++){
    TString posttag = posttags[i];
    
    ////// nevents
    vector<TString> files_nevents = {
      "../root/fit_SYMFIX_10M_UL_8_6_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_grid_"+posttag+".root",
      "../root/fit_SYMFIX_100M_UL_8_6_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_grid_"+posttag+".root",
      "../root/fit_SYMFIX_1G_UL_8_6_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_grid_"+posttag+".root",
      "../root/fit_SYMFIX_4G_UL_8_6_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_grid_"+posttag+".root",
      "../root/fit_SYMFIX_10G_UL_8_6_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_grid_"+posttag+".root",
    };

    vector<float> edges_nevents = {
      0.01,
      0.1,
      1.,
      4.,
      10.,
    };
    
    for(auto& e : edges_nevents) e /= norm_lumi;
    if(logy){for(auto& e : edges_nevents) e = TMath::Log10(e);}
    
    unsigned int n_nevents = files_nevents.size();
    double xx_nevents[10], yy_nevents[10], yyMC_nevents[10], exx_nevents[10], eyy_nevents[10];
    for(unsigned int i=0; i<n_nevents; i++){
      TString fn = files_nevents[i];
      TFile* f = TFile::Open(fn, "READ");
      TTree *t = f->Get<TTree>("tree");
      double err;
      t->SetBranchAddress("best_mass_err", &err);
      t->GetEntry(0);
      xx_nevents[i] = edges_nevents[i];
      exx_nevents[i] = 0.0;
      yy_nevents[i] = err;
      eyy_nevents[i] = 0.0;
      f->Close();
    }
    TGraphErrors* fit_nevents   = new TGraphErrors(n_nevents,xx_nevents,yy_nevents,exx_nevents,eyy_nevents);
    if(i==0){
      fit_nevents->SetLineStyle(kDashed);
      fit_nevents->SetMarkerStyle(kOpenCircle);
    }
    else{
      fit_nevents->SetLineStyle(kSolid);
      fit_nevents->SetMarkerStyle(kFullCircle);
    }
    fit_nevents->SetMarkerSize(1.5);
    fit_nevents->SetMarkerColor(1);
    fit_nevents->SetLineColor(1);
    fit_nevents->SetMinimum(0.0);
    fit_nevents->SetMaximum(55.0);
    TString leg_text = "GRID (8#times6#times6), |y|<2.5";
    if(i==0) leg_text += " (w/ MC stat)";
    leg1->AddEntry(fit_nevents, leg_text, "PL");
    mg->Add(fit_nevents);
  }

  ////////////////////////////////////////////////
  for(unsigned int i = 0; i<posttags.size(); i++){
    TString posttag = posttags[i];
    
    ////// nevents
    vector<TString> files_nevents = {
      "../root/fit_SYMFIX_10M_UL_10_4_A0_3_2_A1_3_3_A2_3_2_A3_3_2_A4_3_3_full_"+posttag+".root",
      "../root/fit_SYMFIX_100M_UL_10_4_A0_3_2_A1_3_3_A2_3_2_A3_3_2_A4_3_3_full_"+posttag+".root",
      "../root/fit_SYMFIX_1G_UL_10_4_A0_3_2_A1_3_3_A2_3_2_A3_3_2_A4_3_3_full_"+posttag+".root",
      "../root/fit_SYMFIX_4G_UL_10_4_A0_3_2_A1_3_3_A2_3_2_A3_3_2_A4_3_3_full_"+posttag+".root",
      "../root/fit_SYMFIX_10G_UL_10_4_A0_3_2_A1_3_3_A2_3_2_A3_3_2_A4_3_3_full_"+posttag+".root",
      "../root/fit_SYMFIX_40G_UL_10_4_A0_3_2_A1_3_3_A2_3_2_A3_3_2_A4_3_3_full_"+posttag+".root",
    };

    vector<float> edges_nevents = {
      0.01,
      0.1,
      1.,
      4.,
      10.,
      40.,
    };
    
    for(auto& e : edges_nevents) e /= norm_lumi;
    if(logy){for(auto& e : edges_nevents) e = TMath::Log10(e);}
    
    unsigned int n_nevents = files_nevents.size();
    double xx_nevents[10], yy_nevents[10], yyMC_nevents[10], exx_nevents[10], eyy_nevents[10];
    for(unsigned int i=0; i<n_nevents; i++){
      TString fn = files_nevents[i];
      TFile* f = TFile::Open(fn, "READ");
      TTree *t = f->Get<TTree>("tree");
      double err;
      t->SetBranchAddress("best_mass_err", &err);
      t->GetEntry(0);
      xx_nevents[i] = edges_nevents[i];
      exx_nevents[i] = 0.0;
      yy_nevents[i] = err;
      eyy_nevents[i] = 0.0;
      f->Close();
    }
    TGraphErrors* fit_nevents   = new TGraphErrors(n_nevents,xx_nevents,yy_nevents,exx_nevents,eyy_nevents);
    if(i==0){
      fit_nevents->SetLineStyle(kDashed);
      fit_nevents->SetMarkerStyle(kOpenCircle);
    }
    else{
      fit_nevents->SetLineStyle(kSolid);
      fit_nevents->SetMarkerStyle(kFullCircle);
    }
    fit_nevents->SetMarkerSize(1.5);
    fit_nevents->SetMarkerColor(2);
    fit_nevents->SetLineColor(2);
    fit_nevents->SetMinimum(0.0);
    fit_nevents->SetMaximum(55.0);
    TString leg_text = "FULL (pol_{10,4} #times pol_{3,2})";
    if(i==0) leg_text += " (w/ MC stat)";
    leg1->AddEntry(fit_nevents, leg_text, "PL");
    mg->Add(fit_nevents);
  }

  ////////////////////////////////////////////////
  for(unsigned int i = 0; i<posttags.size(); i++){
    TString posttag = posttags[i];
    
    ////// nevents
    vector<TString> files_nevents = {
      "../root/fit_SYMFIX_10M_UL_10_6_A0_3_2_A1_3_3_A2_3_2_A3_3_2_A4_3_3_full_"+posttag+".root",
      "../root/fit_SYMFIX_100M_UL_10_6_A0_3_2_A1_3_3_A2_3_2_A3_3_2_A4_3_3_full_"+posttag+".root",
      "../root/fit_SYMFIX_1G_UL_10_6_A0_3_2_A1_3_3_A2_3_2_A3_3_2_A4_3_3_full_"+posttag+".root",
      "../root/fit_SYMFIX_4G_UL_10_6_A0_3_2_A1_3_3_A2_3_2_A3_3_2_A4_3_3_full_"+posttag+".root",
      "../root/fit_SYMFIX_10G_UL_10_6_A0_3_2_A1_3_3_A2_3_2_A3_3_2_A4_3_3_full_"+posttag+".root",
    };

    vector<float> edges_nevents = {
      0.01,
      0.1,
      1.,
      4.,
      10.,
    };
    
    for(auto& e : edges_nevents) e /= norm_lumi;
    if(logy){for(auto& e : edges_nevents) e = TMath::Log10(e);}
    
    unsigned int n_nevents = files_nevents.size();
    double xx_nevents[10], yy_nevents[10], yyMC_nevents[10], exx_nevents[10], eyy_nevents[10];
    for(unsigned int i=0; i<n_nevents; i++){
      TString fn = files_nevents[i];
      TFile* f = TFile::Open(fn, "READ");
      TTree *t = f->Get<TTree>("tree");
      double err;
      t->SetBranchAddress("best_mass_err", &err);
      t->GetEntry(0);
      xx_nevents[i] = edges_nevents[i];
      exx_nevents[i] = 0.0;
      yy_nevents[i] = err;
      eyy_nevents[i] = 0.0;
      f->Close();
    }
    TGraphErrors* fit_nevents   = new TGraphErrors(n_nevents,xx_nevents,yy_nevents,exx_nevents,eyy_nevents);
    if(i==0){
      fit_nevents->SetLineStyle(kDashed);
      fit_nevents->SetMarkerStyle(kOpenCircle);
    }
    else{
      fit_nevents->SetLineStyle(kSolid);
      fit_nevents->SetMarkerStyle(kFullCircle);
    }
    fit_nevents->SetMarkerSize(1.5);
    fit_nevents->SetMarkerColor(2);
    fit_nevents->SetLineColor(2);
    fit_nevents->SetLineStyle(kDashed);
    fit_nevents->SetMinimum(0.0);
    fit_nevents->SetMaximum(55.0);
    TString leg_text = "FULL (pol_{10,6} #times pol_{3,2})";
    if(i==0) leg_text += " (w/ MC stat)";
    //leg1->AddEntry(fit_nevents, leg_text, "PL");
    //mg->Add(fit_nevents);
  }

  ////////////////////////////////////////////////
  for(unsigned int i = 0; i<posttags.size(); i++){
    TString posttag = posttags[i];
    
    ////// nevents
    vector<TString> files_nevents = {
      "../root/fit_SYMFIX_10M_UL_12_4_A0_3_2_A1_3_3_A2_3_2_A3_3_2_A4_3_3_full_"+posttag+".root",
      "../root/fit_SYMFIX_100M_UL_12_4_A0_3_2_A1_3_3_A2_3_2_A3_3_2_A4_3_3_full_"+posttag+".root",
      "../root/fit_SYMFIX_1G_UL_12_4_A0_3_2_A1_3_3_A2_3_2_A3_3_2_A4_3_3_full_"+posttag+".root",
      "../root/fit_SYMFIX_4G_UL_12_4_A0_3_2_A1_3_3_A2_3_2_A3_3_2_A4_3_3_full_"+posttag+".root",
      "../root/fit_SYMFIX_10G_UL_12_4_A0_3_2_A1_3_3_A2_3_2_A3_3_2_A4_3_3_full_"+posttag+".root",
    };

    vector<float> edges_nevents = {
      0.01,
      0.1,
      1.,
      4.,
      10.,
    };
    
    for(auto& e : edges_nevents) e /= norm_lumi;
    if(logy){for(auto& e : edges_nevents) e = TMath::Log10(e);}
    
    unsigned int n_nevents = files_nevents.size();
    double xx_nevents[10], yy_nevents[10], yyMC_nevents[10], exx_nevents[10], eyy_nevents[10];
    for(unsigned int i=0; i<n_nevents; i++){
      TString fn = files_nevents[i];
      TFile* f = TFile::Open(fn, "READ");
      TTree *t = f->Get<TTree>("tree");
      double err;
      t->SetBranchAddress("best_mass_err", &err);
      t->GetEntry(0);
      xx_nevents[i] = edges_nevents[i];
      exx_nevents[i] = 0.0;
      yy_nevents[i] = err;
      eyy_nevents[i] = 0.0;
      f->Close();
    }
    TGraphErrors* fit_nevents   = new TGraphErrors(n_nevents,xx_nevents,yy_nevents,exx_nevents,eyy_nevents);
    if(i==0){
      fit_nevents->SetLineStyle(kDashed);
      fit_nevents->SetMarkerStyle(kOpenCircle);
    }
    else{
      fit_nevents->SetLineStyle(kSolid);
      fit_nevents->SetMarkerStyle(kFullCircle);
    }
    fit_nevents->SetMarkerSize(1.5);
    fit_nevents->SetMarkerColor(2);
    fit_nevents->SetLineColor(2);
    fit_nevents->SetLineStyle(kDotted);
    fit_nevents->SetMinimum(0.0);
    fit_nevents->SetMaximum(55.0);
    TString leg_text = "FULL (pol_{12,4} #times pol_{3,2})";
    if(i==0) leg_text += " (w/ MC stat)";
    //leg1->AddEntry(fit_nevents, leg_text, "PL");
    //mg->Add(fit_nevents);
  }

    ////////////////////////////////////////////////
  for(unsigned int i = 0; i<posttags.size(); i++){
    TString posttag = posttags[i];
    
    ////// nevents
    vector<TString> files_nevents = {
      "../root/fit_SYMFIXCORR_10M_UL_5_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_100M_UL_5_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_1G_UL_5_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_4G_UL_5_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_10G_UL_5_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_40G_UL_5_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
    };

    vector<float> edges_nevents = {
      0.01,
      0.1,
      1.,
      4.,
      10.,
      40.,
    };
    
    for(auto& e : edges_nevents) e /= norm_lumi;
    if(logy){for(auto& e : edges_nevents) e = TMath::Log10(e);}
    
    unsigned int n_nevents = files_nevents.size();
    double xx_nevents[10], yy_nevents[10], yyMC_nevents[10], exx_nevents[10], eyy_nevents[10];
    for(unsigned int i=0; i<n_nevents; i++){
      TString fn = files_nevents[i];
      TFile* f = TFile::Open(fn, "READ");
      TTree *t = f->Get<TTree>("tree");
      double err;
      t->SetBranchAddress("best_mass_err", &err);
      t->GetEntry(0);
      xx_nevents[i] = edges_nevents[i];
      exx_nevents[i] = 0.0;
      yy_nevents[i] = err;
      eyy_nevents[i] = 0.0;
      f->Close();
    }
    TGraphErrors* fit_nevents   = new TGraphErrors(n_nevents,xx_nevents,yy_nevents,exx_nevents,eyy_nevents);
    if(i==0){
      fit_nevents->SetLineStyle(kDashed);
      fit_nevents->SetMarkerStyle(kOpenCircle);
    }
    else{
      fit_nevents->SetLineStyle(kSolid);
      fit_nevents->SetMarkerStyle(kFullCircle);
    }
    fit_nevents->SetMarkerSize(1.5);
    fit_nevents->SetMarkerColor(4);
    fit_nevents->SetLineColor(4);
    fit_nevents->SetMinimum(0.0);
    fit_nevents->SetMaximum(55.0);
    TString leg_text = "CORR (pol_{5,2} #times pol_{2,2})";
    if(i==0) leg_text += " (w/ MC stat)";
    leg1->AddEntry(fit_nevents, leg_text, "PL");
    mg->Add(fit_nevents);
  }


  ////////////////////////////////////////////////
  for(unsigned int i = 0; i<posttags.size(); i++){
    TString posttag = posttags[i];
    
    ////// nevents
    vector<TString> files_nevents = {
      "../root/fit_SYMFIXCORR_10M_UL_4_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_100M_UL_4_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_1G_UL_4_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_4G_UL_4_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_10G_UL_4_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_40G_UL_4_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
    };

    vector<float> edges_nevents = {
      0.01,
      0.1,
      1.,
      4.,
      10.,
      40.,
    };
    
    for(auto& e : edges_nevents) e /= norm_lumi;
    if(logy){for(auto& e : edges_nevents) e = TMath::Log10(e);}
    
    unsigned int n_nevents = files_nevents.size();
    double xx_nevents[10], yy_nevents[10], yyMC_nevents[10], exx_nevents[10], eyy_nevents[10];
    for(unsigned int i=0; i<n_nevents; i++){
      TString fn = files_nevents[i];
      TFile* f = TFile::Open(fn, "READ");
      TTree *t = f->Get<TTree>("tree");
      double err;
      t->SetBranchAddress("best_mass_err", &err);
      t->GetEntry(0);
      xx_nevents[i] = edges_nevents[i];
      exx_nevents[i] = 0.0;
      yy_nevents[i] = err;
      eyy_nevents[i] = 0.0;
      f->Close();
    }
    TGraphErrors* fit_nevents   = new TGraphErrors(n_nevents,xx_nevents,yy_nevents,exx_nevents,eyy_nevents);
    if(i==0){
      fit_nevents->SetLineStyle(kDashed);
      fit_nevents->SetMarkerStyle(kOpenCircle);
    }
    else{
      fit_nevents->SetLineStyle(kSolid);
      fit_nevents->SetMarkerStyle(kFullCircle);
    }
    fit_nevents->SetMarkerSize(1.5);
    fit_nevents->SetMarkerColor(5);
    fit_nevents->SetLineColor(5);
    fit_nevents->SetMinimum(0.0);
    fit_nevents->SetMaximum(55.0);
    TString leg_text = "CORR (pol_{4,2} #times pol_{2,2})";
    if(i==0) leg_text += " (w/ MC stat)";
    leg1->AddEntry(fit_nevents, leg_text, "PL");
    mg->Add(fit_nevents);
  }

    ////////////////////////////////////////////////
  for(unsigned int i = 0; i<posttags.size(); i++){
    TString posttag = posttags[i];
    
    ////// nevents
    vector<TString> files_nevents = {
      "../root/fit_SYMFIXCORR_10M_UL_3_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_100M_UL_3_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_1G_UL_3_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_4G_UL_3_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_10G_UL_3_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_40G_UL_3_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr_"+posttag+".root",
    };

    vector<float> edges_nevents = {
      0.01,
      0.1,
      1.,
      4.,
      10.,
      40.,
    };
    
    for(auto& e : edges_nevents) e /= norm_lumi;
    if(logy){for(auto& e : edges_nevents) e = TMath::Log10(e);}
    
    unsigned int n_nevents = files_nevents.size();
    double xx_nevents[10], yy_nevents[10], yyMC_nevents[10], exx_nevents[10], eyy_nevents[10];
    for(unsigned int i=0; i<n_nevents; i++){
      TString fn = files_nevents[i];
      TFile* f = TFile::Open(fn, "READ");
      TTree *t = f->Get<TTree>("tree");
      double err;
      t->SetBranchAddress("best_mass_err", &err);
      t->GetEntry(0);
      xx_nevents[i] = edges_nevents[i];
      exx_nevents[i] = 0.0;
      yy_nevents[i] = err;
      eyy_nevents[i] = 0.0;
      f->Close();
    }
    TGraphErrors* fit_nevents   = new TGraphErrors(n_nevents,xx_nevents,yy_nevents,exx_nevents,eyy_nevents);
    if(i==0){
      fit_nevents->SetLineStyle(kDashed);
      fit_nevents->SetMarkerStyle(kOpenCircle);
    }
    else{
      fit_nevents->SetLineStyle(kSolid);
      fit_nevents->SetMarkerStyle(kFullCircle);
    }
    fit_nevents->SetMarkerSize(1.5);
    fit_nevents->SetMarkerColor(3);
    fit_nevents->SetLineColor(3);
    fit_nevents->SetMinimum(0.0);
    fit_nevents->SetMaximum(55.0);
    TString leg_text = "CORR (pol_{3,2} #times pol_{2,2})";
    if(i==0) leg_text += " (w/ MC stat)";
    leg1->AddEntry(fit_nevents, leg_text, "PL");
    mg->Add(fit_nevents);
  }

  ////////////////////////////////////////////////
  for(unsigned int i = 0; i<posttags.size(); i++){
    TString posttag = posttags[i];
    
    ////// nevents
    vector<TString> files_nevents = {
      "../root/fit_SYMFIXCORR_10M_UL_3_2_A0_1_2_A1_1_2_A2_1_2_A3_1_2_A4_1_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_100M_UL_3_2_A0_1_2_A1_1_2_A2_1_2_A3_1_2_A4_1_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_1G_UL_3_2_A0_1_2_A1_1_2_A2_1_2_A3_1_2_A4_1_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_4G_UL_3_2_A0_1_2_A1_1_2_A2_1_2_A3_1_2_A4_1_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_10G_UL_3_2_A0_1_2_A1_1_2_A2_1_2_A3_1_2_A4_1_2_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_40G_UL_3_2_A0_1_2_A1_1_2_A2_1_2_A3_1_2_A4_1_2_corr_"+posttag+".root",
    };

    vector<float> edges_nevents = {
      0.01,
      0.1,
      1.,
      4.,
      10.,
      40.,
    };
    
    for(auto& e : edges_nevents) e /= norm_lumi;
    if(logy){for(auto& e : edges_nevents) e = TMath::Log10(e);}
    
    unsigned int n_nevents = files_nevents.size();
    double xx_nevents[10], yy_nevents[10], yyMC_nevents[10], exx_nevents[10], eyy_nevents[10];
    for(unsigned int i=0; i<n_nevents; i++){
      TString fn = files_nevents[i];
      TFile* f = TFile::Open(fn, "READ");
      TTree *t = f->Get<TTree>("tree");
      double err;
      t->SetBranchAddress("best_mass_err", &err);
      t->GetEntry(0);
      xx_nevents[i] = edges_nevents[i];
      exx_nevents[i] = 0.0;
      yy_nevents[i] = err;
      eyy_nevents[i] = 0.0;
      f->Close();
    }
    TGraphErrors* fit_nevents   = new TGraphErrors(n_nevents,xx_nevents,yy_nevents,exx_nevents,eyy_nevents);
    if(i==0){
      fit_nevents->SetLineStyle(kDashed);
      fit_nevents->SetMarkerStyle(kOpenCircle);
    }
    else{
      fit_nevents->SetLineStyle(kSolid);
      fit_nevents->SetMarkerStyle(kFullCircle);
    }
    fit_nevents->SetMarkerSize(1.5);
    fit_nevents->SetMarkerColor(6);
    fit_nevents->SetLineColor(6);
    fit_nevents->SetMinimum(0.0);
    fit_nevents->SetMaximum(55.0);
    TString leg_text = "CORR (pol_{3,2} #times pol_{1,2})";
    if(i==0) leg_text += " (w/ MC stat)";
    leg1->AddEntry(fit_nevents, leg_text, "PL");
    mg->Add(fit_nevents);
  }

  ////////////////////////////////////////////////
  for(unsigned int i = 0; i<posttags.size(); i++){
    TString posttag = posttags[i];
    
    ////// nevents
    vector<TString> files_nevents = {
      "../root/fit_SYMFIXCORR_10M_UL_3_2_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_100M_UL_3_2_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_1G_UL_3_2_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_4G_UL_3_2_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_10G_UL_3_2_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_corr_"+posttag+".root",
      "../root/fit_SYMFIXCORR_40G_UL_3_2_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_corr_"+posttag+".root",
    };

    vector<float> edges_nevents = {
      0.01,
      0.1,
      1.,
      4.,
      10.,
      40.,
    };
    
    for(auto& e : edges_nevents) e /= norm_lumi;
    if(logy){for(auto& e : edges_nevents) e = TMath::Log10(e);}
    
    unsigned int n_nevents = files_nevents.size();
    double xx_nevents[10], yy_nevents[10], yyMC_nevents[10], exx_nevents[10], eyy_nevents[10];
    for(unsigned int i=0; i<n_nevents; i++){
      TString fn = files_nevents[i];
      TFile* f = TFile::Open(fn, "READ");
      TTree *t = f->Get<TTree>("tree");
      double err;
      t->SetBranchAddress("best_mass_err", &err);
      t->GetEntry(0);
      xx_nevents[i] = edges_nevents[i];
      exx_nevents[i] = 0.0;
      yy_nevents[i] = err;
      eyy_nevents[i] = 0.0;
      f->Close();
    }
    TGraphErrors* fit_nevents   = new TGraphErrors(n_nevents,xx_nevents,yy_nevents,exx_nevents,eyy_nevents);
    if(i==0){
      fit_nevents->SetLineStyle(kDashed);
      fit_nevents->SetMarkerStyle(kOpenCircle);
    }
    else{
      fit_nevents->SetLineStyle(kSolid);
      fit_nevents->SetMarkerStyle(kFullCircle);
    }
    fit_nevents->SetMarkerSize(1.5);
    fit_nevents->SetMarkerColor(7);
    fit_nevents->SetLineColor(7);
    fit_nevents->SetMinimum(0.0);
    fit_nevents->SetMaximum(55.0);
    TString leg_text = "CORR (pol_{3,2} #times pol_{1,1})";
    if(i==0) leg_text += " (w/ MC stat)";
    leg1->AddEntry(fit_nevents, leg_text, "PL");
    mg->Add(fit_nevents);
  }

  ////////////////////////////////////////////////////////
  mg->GetXaxis()->SetTitle("N_{tot}");
  mg->GetXaxis()->SetTitleSize(0.05);
  mg->GetXaxis()->SetTitleOffset(0.8);
  mg->GetYaxis()->SetTitleSize(0.06);
  mg->GetYaxis()->SetTitleOffset(0.7);
  mg->GetYaxis()->SetTitle("#Delta M_{W} (MeV)");
  if(logy) mg->GetXaxis()->SetTitle("log_{10}(N_{unweigh.}/N_{tot})");
  mg->GetYaxis()->SetRangeUser(4, 1000.);
  mg->GetXaxis()->SetRangeUser(0.0001, 100e+09);
  if(logy) mg->GetXaxis()->SetRangeUser(-4, 5);

  mg->GetYaxis()->SetMoreLogLabels();
  mg->GetYaxis()->CenterTitle(true);
  
  mg->Draw("apl");
  c->SetLogy();
  leg1->Draw();
  c->SaveAs(image_name);
  //c->SaveAs("deltaM_vs_y_finalfix.png");

  //TF1* func = new TF1("func", "[0]*TMath::Log(x)+[1]", 0.0001, 50e+09);
  //fit_nevents->Fit(func);

  return;
}
