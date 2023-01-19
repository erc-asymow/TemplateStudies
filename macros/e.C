{

  TCanvas* c = new TCanvas("c", "canvas", 1000, 600);
  TMultiGraph* mg = new TMultiGraph();
  TLegend* leg1 = new TLegend(0.25, 0.75, 0.85, 0.90, "","brNDC");
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetFillColor(10);
  //leg1->SetHeader("");
  leg1->SetNColumns(3);

  ////// SCALEA0A1A2A3A4
  vector<TString> files_SCALEA0A1A2A3A4 = {
    //"../root/fit_SYSTS_1G_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST.root",
    "../root/fit_jacVsM_4G_y2p75_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3A4.root",
    "../root/fit_jacVsM_4G_y3p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3A4.root",
    "../root/fit_jacVsM_4G_y3p25_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3A4.root",
    "../root/fit_jacVsM_4G_y3p50_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3A4.root",
    "../root/fit_jacVsM_4G_y3p75_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3A4.root",
    "../root/fit_jacVsM_4G_y4p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3A4.root",
    "../root/fit_jacVsM_4G_y4p25_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3A4.root",
    "../root/fit_jacVsM_4G_y4p50_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3A4.root",
    "../root/fit_jacVsM_4G_y5p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3A4.root"
  };
  vector<float> edges_SCALEA0A1A2A3A4 = {
    //2.50,
    2.75,
    3.0,
    3.25,
    3.5,
    3.75,
    4.0,
    4.25,
    4.5,
    5.0
  };
  unsigned int n_SCALEA0A1A2A3A4 = files_SCALEA0A1A2A3A4.size();
  double xx_SCALEA0A1A2A3A4[10], yy_SCALEA0A1A2A3A4[10], yyMC_SCALEA0A1A2A3A4[10], exx_SCALEA0A1A2A3A4[10], eyy_SCALEA0A1A2A3A4[10];
  for(unsigned int i=0; i<n_SCALEA0A1A2A3A4; i++){
    TString fn = files_SCALEA0A1A2A3A4[i];
    TFile* f = TFile::Open(fn, "READ");
    TTree *t = f->Get<TTree>("tree");
    double err;
    t->SetBranchAddress("best_mass_err", &err);
    t->GetEntry(0);
    xx_SCALEA0A1A2A3A4[i] = edges_SCALEA0A1A2A3A4[i];
    exx_SCALEA0A1A2A3A4[i] = 0.0;
    yy_SCALEA0A1A2A3A4[i] = err;
    eyy_SCALEA0A1A2A3A4[i] = 0.0;
    //cout << yy[i] << endl;
  }
  TGraphErrors* fit_SCALEA0A1A2A3A4   = new TGraphErrors(n_SCALEA0A1A2A3A4,xx_SCALEA0A1A2A3A4,yy_SCALEA0A1A2A3A4,exx_SCALEA0A1A2A3A4,eyy_SCALEA0A1A2A3A4);
  fit_SCALEA0A1A2A3A4->SetMarkerStyle(kFullCircle);
  fit_SCALEA0A1A2A3A4->SetMarkerSize(1.5);
  fit_SCALEA0A1A2A3A4->SetMarkerColor(kGreen);
  fit_SCALEA0A1A2A3A4->SetMinimum(0.0);
  fit_SCALEA0A1A2A3A4->SetMaximum(110.0);
  //fit_SCALEA0A1A2A3A4->Draw("apl");
  leg1->AddEntry(fit_SCALEA0A1A2A3A4, "Float A_{0,1,2,3,4}", "P");
  mg->Add(fit_SCALEA0A1A2A3A4);

  ////// SCALEA0A1A2A4
  vector<TString> files_SCALEA0A1A2A4 = {
    //"../root/fit_SYSTS_1G_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST.root",
    "../root/fit_jacVsM_4G_y2p75_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A4.root",
    "../root/fit_jacVsM_4G_y3p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A4.root",
    "../root/fit_jacVsM_4G_y3p25_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A4.root",
    "../root/fit_jacVsM_4G_y3p50_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A4.root",
    "../root/fit_jacVsM_4G_y3p75_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A4.root",
    "../root/fit_jacVsM_4G_y4p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A4.root",
    "../root/fit_jacVsM_4G_y4p25_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A4.root",
    "../root/fit_jacVsM_4G_y4p50_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A4.root",
    "../root/fit_jacVsM_4G_y5p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A4.root"
  };
  vector<float> edges_SCALEA0A1A2A4 = {
    //2.50,
    2.75,
    3.0,
    3.25,
    3.5,
    3.75,
    4.0,
    4.25,
    4.5,
    5.0
  };
  unsigned int n_SCALEA0A1A2A4 = files_SCALEA0A1A2A4.size();
  double xx_SCALEA0A1A2A4[10], yy_SCALEA0A1A2A4[10], yyMC_SCALEA0A1A2A4[10], exx_SCALEA0A1A2A4[10], eyy_SCALEA0A1A2A4[10];
  for(unsigned int i=0; i<n_SCALEA0A1A2A4; i++){
    TString fn = files_SCALEA0A1A2A4[i];
    TFile* f = TFile::Open(fn, "READ");
    TTree *t = f->Get<TTree>("tree");
    double err;
    t->SetBranchAddress("best_mass_err", &err);
    t->GetEntry(0);
    xx_SCALEA0A1A2A4[i] = edges_SCALEA0A1A2A4[i];
    exx_SCALEA0A1A2A4[i] = 0.0;
    yy_SCALEA0A1A2A4[i] = err;
    eyy_SCALEA0A1A2A4[i] = 0.0;
    //cout << yy[i] << endl;
  }
  TGraphErrors* fit_SCALEA0A1A2A4   = new TGraphErrors(n_SCALEA0A1A2A4,xx_SCALEA0A1A2A4,yy_SCALEA0A1A2A4,exx_SCALEA0A1A2A4,eyy_SCALEA0A1A2A4);
  fit_SCALEA0A1A2A4->SetMarkerStyle(kFullCircle);
  fit_SCALEA0A1A2A4->SetMarkerSize(1.5);
  fit_SCALEA0A1A2A4->SetMarkerColor(kMagenta);
  fit_SCALEA0A1A2A4->SetMinimum(0.0);
  fit_SCALEA0A1A2A4->SetMaximum(110.0);
  //fit_SCALEA0A1A2A4->Draw("apl");
  leg1->AddEntry(fit_SCALEA0A1A2A4, "Float A_{0,1,2,4}", "P");
  mg->Add(fit_SCALEA0A1A2A4);

  ////// SCALEA0A1A3A4
  vector<TString> files_SCALEA0A1A3A4 = {
    //"../root/fit_SYSTS_1G_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST.root",
    "../root/fit_jacVsM_4G_y2p75_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A3A4.root",
    "../root/fit_jacVsM_4G_y3p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A3A4.root",
    "../root/fit_jacVsM_4G_y3p25_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A3A4.root",
    "../root/fit_jacVsM_4G_y3p50_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A3A4.root",
    "../root/fit_jacVsM_4G_y3p75_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A3A4.root",
    "../root/fit_jacVsM_4G_y4p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A3A4.root",
    "../root/fit_jacVsM_4G_y4p25_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A3A4.root",
    "../root/fit_jacVsM_4G_y4p50_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A3A4.root",
    "../root/fit_jacVsM_4G_y5p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A3A4.root"
  };
  vector<float> edges_SCALEA0A1A3A4 = {
    //2.50,
    2.75,
    3.0,
    3.25,
    3.5,
    3.75,
    4.0,
    4.25,
    4.5,
    5.0
  };
  unsigned int n_SCALEA0A1A3A4 = files_SCALEA0A1A3A4.size();
  double xx_SCALEA0A1A3A4[10], yy_SCALEA0A1A3A4[10], yyMC_SCALEA0A1A3A4[10], exx_SCALEA0A1A3A4[10], eyy_SCALEA0A1A3A4[10];
  for(unsigned int i=0; i<n_SCALEA0A1A3A4; i++){
    TString fn = files_SCALEA0A1A3A4[i];
    TFile* f = TFile::Open(fn, "READ");
    TTree *t = f->Get<TTree>("tree");
    double err;
    t->SetBranchAddress("best_mass_err", &err);
    t->GetEntry(0);
    xx_SCALEA0A1A3A4[i] = edges_SCALEA0A1A3A4[i];
    exx_SCALEA0A1A3A4[i] = 0.0;
    yy_SCALEA0A1A3A4[i] = err;
    eyy_SCALEA0A1A3A4[i] = 0.0;
    //cout << yy[i] << endl;
  }
  TGraphErrors* fit_SCALEA0A1A3A4   = new TGraphErrors(n_SCALEA0A1A3A4,xx_SCALEA0A1A3A4,yy_SCALEA0A1A3A4,exx_SCALEA0A1A3A4,eyy_SCALEA0A1A3A4);
  fit_SCALEA0A1A3A4->SetMarkerStyle(kFullCircle);
  fit_SCALEA0A1A3A4->SetMarkerSize(1.5);
  fit_SCALEA0A1A3A4->SetMarkerColor(kRed-2);
  fit_SCALEA0A1A3A4->SetMinimum(0.0);
  fit_SCALEA0A1A3A4->SetMaximum(110.0);
  //fit_SCALEA0A1A3A4->Draw("apl");
  leg1->AddEntry(fit_SCALEA0A1A3A4, "Float A_{0,1,3,4}", "P");
  mg->Add(fit_SCALEA0A1A3A4);

  ////// SCALEA0A1A2A3
  vector<TString> files_SCALEA0A1A2A3 = {
    //"../root/fit_SYSTS_1G_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST.root",
    "../root/fit_jacVsM_4G_y2p75_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3.root",
    "../root/fit_jacVsM_4G_y3p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3.root",
    "../root/fit_jacVsM_4G_y3p25_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3.root",
    "../root/fit_jacVsM_4G_y3p50_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3.root",
    "../root/fit_jacVsM_4G_y3p75_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3.root",
    "../root/fit_jacVsM_4G_y4p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3.root",
    "../root/fit_jacVsM_4G_y4p25_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3.root",
    "../root/fit_jacVsM_4G_y4p50_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3.root",
    "../root/fit_jacVsM_4G_y5p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3.root"
  };
  vector<float> edges_SCALEA0A1A2A3 = {
    //2.50,
    2.75,
    3.0,
    3.25,
    3.5,
    3.75,
    4.0,
    4.25,
    4.5,
    5.0
  };
  unsigned int n_SCALEA0A1A2A3 = files_SCALEA0A1A2A3.size();
  double xx_SCALEA0A1A2A3[10], yy_SCALEA0A1A2A3[10], yyMC_SCALEA0A1A2A3[10], exx_SCALEA0A1A2A3[10], eyy_SCALEA0A1A2A3[10];
  for(unsigned int i=0; i<n_SCALEA0A1A2A3; i++){
    TString fn = files_SCALEA0A1A2A3[i];
    TFile* f = TFile::Open(fn, "READ");
    TTree *t = f->Get<TTree>("tree");
    double err;
    t->SetBranchAddress("best_mass_err", &err);
    t->GetEntry(0);
    xx_SCALEA0A1A2A3[i] = edges_SCALEA0A1A2A3[i];
    exx_SCALEA0A1A2A3[i] = 0.0;
    yy_SCALEA0A1A2A3[i] = err;
    eyy_SCALEA0A1A2A3[i] = 0.0;
    //cout << yy[i] << endl;
  }
  TGraphErrors* fit_SCALEA0A1A2A3   = new TGraphErrors(n_SCALEA0A1A2A3,xx_SCALEA0A1A2A3,yy_SCALEA0A1A2A3,exx_SCALEA0A1A2A3,eyy_SCALEA0A1A2A3);
  fit_SCALEA0A1A2A3->SetMarkerStyle(kFullCircle);
  fit_SCALEA0A1A2A3->SetMarkerSize(1.5);
  fit_SCALEA0A1A2A3->SetMarkerColor(kBlue);
  fit_SCALEA0A1A2A3->SetMinimum(0.0);
  fit_SCALEA0A1A2A3->SetMaximum(110.0);
  //fit_SCALEA0A1A2A3->Draw("apl");
  leg1->AddEntry(fit_SCALEA0A1A2A3, "Float A_{0,1,2,3}", "P");
  mg->Add(fit_SCALEA0A1A2A3);

  ////// SCALEA0A1A2
  vector<TString> files_SCALEA0A1A2 = {
    //"../root/fit_SYSTS_1G_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST.root",
    "../root/fit_jacVsM_4G_y2p75_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2.root",
    "../root/fit_jacVsM_4G_y3p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2.root",
    "../root/fit_jacVsM_4G_y3p25_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2.root",
    "../root/fit_jacVsM_4G_y3p50_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2.root",
    "../root/fit_jacVsM_4G_y3p75_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2.root",
    "../root/fit_jacVsM_4G_y4p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2.root",
    "../root/fit_jacVsM_4G_y4p25_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2.root",
    "../root/fit_jacVsM_4G_y4p50_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2.root",
    "../root/fit_jacVsM_4G_y5p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2.root"
  };
  vector<float> edges_SCALEA0A1A2 = {
    //2.50,
    2.75,
    3.0,
    3.25,
    3.5,
    3.75,
    4.0,
    4.25,
    4.5,
    5.0
  };
  unsigned int n_SCALEA0A1A2 = files_SCALEA0A1A2.size();
  double xx_SCALEA0A1A2[10], yy_SCALEA0A1A2[10], yyMC_SCALEA0A1A2[10], exx_SCALEA0A1A2[10], eyy_SCALEA0A1A2[10];
  for(unsigned int i=0; i<n_SCALEA0A1A2; i++){
    TString fn = files_SCALEA0A1A2[i];
    TFile* f = TFile::Open(fn, "READ");
    TTree *t = f->Get<TTree>("tree");
    double err;
    t->SetBranchAddress("best_mass_err", &err);
    t->GetEntry(0);
    xx_SCALEA0A1A2[i] = edges_SCALEA0A1A2[i];
    exx_SCALEA0A1A2[i] = 0.0;
    yy_SCALEA0A1A2[i] = err;
    eyy_SCALEA0A1A2[i] = 0.0;
    //cout << yy[i] << endl;
  }
  TGraphErrors* fit_SCALEA0A1A2   = new TGraphErrors(n_SCALEA0A1A2,xx_SCALEA0A1A2,yy_SCALEA0A1A2,exx_SCALEA0A1A2,eyy_SCALEA0A1A2);
  fit_SCALEA0A1A2->SetMarkerStyle(kFullCircle);
  fit_SCALEA0A1A2->SetMarkerSize(1.5);
  fit_SCALEA0A1A2->SetMarkerColor(kRed);
  fit_SCALEA0A1A2->SetMinimum(0.0);
  fit_SCALEA0A1A2->SetMaximum(110.0);
  //fit_SCALEA0A1A2->Draw("apl");
  leg1->AddEntry(fit_SCALEA0A1A2, "Float A_{0,1,2}", "P");
  mg->Add(fit_SCALEA0A1A2);

  ////// SCALEA0A1
  vector<TString> files_SCALEA0A1 = {
    //"../root/fit_SYSTS_1G_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST.root",
    "../root/fit_jacVsM_4G_y2p75_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1.root",
    "../root/fit_jacVsM_4G_y3p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1.root",
    "../root/fit_jacVsM_4G_y3p25_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1.root",
    "../root/fit_jacVsM_4G_y3p50_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1.root",
    "../root/fit_jacVsM_4G_y3p75_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1.root",
    "../root/fit_jacVsM_4G_y4p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1.root",
    "../root/fit_jacVsM_4G_y4p25_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1.root",
    "../root/fit_jacVsM_4G_y4p50_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1.root",
    "../root/fit_jacVsM_4G_y5p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1.root"
  };
  vector<float> edges_SCALEA0A1 = {
    //2.50,
    2.75,
    3.0,
    3.25,
    3.5,
    3.75,
    4.0,
    4.25,
    4.5,
    5.0
  };
  unsigned int n_SCALEA0A1 = files_SCALEA0A1.size();
  double xx_SCALEA0A1[10], yy_SCALEA0A1[10], yyMC_SCALEA0A1[10], exx_SCALEA0A1[10], eyy_SCALEA0A1[10];
  for(unsigned int i=0; i<n_SCALEA0A1; i++){
    TString fn = files_SCALEA0A1[i];
    TFile* f = TFile::Open(fn, "READ");
    TTree *t = f->Get<TTree>("tree");
    double err;
    t->SetBranchAddress("best_mass_err", &err);
    t->GetEntry(0);
    xx_SCALEA0A1[i] = edges_SCALEA0A1[i];
    exx_SCALEA0A1[i] = 0.0;
    yy_SCALEA0A1[i] = err;
    eyy_SCALEA0A1[i] = 0.0;
    //cout << yy[i] << endl;
  }
  TGraphErrors* fit_SCALEA0A1   = new TGraphErrors(n_SCALEA0A1,xx_SCALEA0A1,yy_SCALEA0A1,exx_SCALEA0A1,eyy_SCALEA0A1);
  fit_SCALEA0A1->SetMarkerStyle(kFullCircle);
  fit_SCALEA0A1->SetMarkerSize(1.5);
  fit_SCALEA0A1->SetMarkerColor(kOrange);
  fit_SCALEA0A1->SetMinimum(0.0);
  fit_SCALEA0A1->SetMaximum(110.0);
  //fit_SCALEA0A1->Draw("apl");
  leg1->AddEntry(fit_SCALEA0A1, "Float A_{0,1}", "P");
  mg->Add(fit_SCALEA0A1);

  ////// FULL
  vector<TString> files_FULL = {
    //"../root/fit_SYSTS_1G_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST.root",
    "../root/fit_jacVsM_4G_y2p75_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_FULL.root",
    "../root/fit_jacVsM_4G_y3p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_FULL.root",
    "../root/fit_jacVsM_4G_y3p25_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_FULL.root",
    "../root/fit_jacVsM_4G_y3p50_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_FULL.root",
    "../root/fit_jacVsM_4G_y3p75_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_FULL.root",
    "../root/fit_jacVsM_4G_y4p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_FULL.root",
    "../root/fit_jacVsM_4G_y4p25_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_FULL.root",
    "../root/fit_jacVsM_4G_y4p50_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_FULL.root",
    "../root/fit_jacVsM_4G_y5p00_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_FULL.root"
  };
  vector<float> edges_FULL = {
    //2.50,
    2.75,
    3.0,
    3.25,
    3.5,
    3.75,
    4.0,
    4.25,
    4.5,
    5.0
  };
  unsigned int n_FULL = files_FULL.size();
  double xx_FULL[10], yy_FULL[10], yyMC_FULL[10], exx_FULL[10], eyy_FULL[10];
  for(unsigned int i=0; i<n_FULL; i++){
    TString fn = files_FULL[i];
    TFile* f = TFile::Open(fn, "READ");
    TTree *t = f->Get<TTree>("tree");
    double err;
    t->SetBranchAddress("best_mass_err", &err);
    t->GetEntry(0);
    xx_FULL[i] = edges_FULL[i];
    exx_FULL[i] = 0.0;
    yy_FULL[i] = err;
    eyy_FULL[i] = 0.0;
    //cout << yy[i] << endl;
  }
  TGraphErrors* fit_FULL   = new TGraphErrors(n_FULL,xx_FULL,yy_FULL,exx_FULL,eyy_FULL);
  fit_FULL->SetMarkerStyle(kFullCircle);
  fit_FULL->SetMarkerSize(1.5);
  fit_FULL->SetMarkerColor(kGray);
  fit_FULL->SetMinimum(0.0);
  fit_FULL->SetMaximum(110.0);
  //fit_FULL->Draw("apl");
  leg1->AddEntry(fit_FULL, "Fit all", "P");
  mg->Add(fit_FULL);

  //////// 
  vector<TString> files_A4new = {
    "../root/fit_jacVsM_4G_y4p00_A4new_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_FULL.root",
    "../root/fit_jacVsM_4G_y4p00_A4new_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1.root",
    "../root/fit_jacVsM_4G_y4p00_A4new_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2.root",
    "../root/fit_jacVsM_4G_y4p00_A4new_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3.root",
  };
  vector<float> edges_A4new = {
    4.0-0.04,
    4.0-0.04,
    4.0-0.04,
    4.0-0.04
  };
  unsigned int n_A4new = files_A4new.size();
  double xx_A4new[10], yy_A4new[10], yyMC_A4new[10], exx_A4new[10], eyy_A4new[10];
  for(unsigned int i=0; i<n_A4new; i++){
    TString fn = files_A4new[i];
    TFile* f = TFile::Open(fn, "READ");
    TTree *t = f->Get<TTree>("tree");
    double err;
    t->SetBranchAddress("best_mass_err", &err);
    t->GetEntry(0);
    xx_A4new[i] = edges_A4new[i];
    exx_A4new[i] = 0.0;
    yy_A4new[i] = err;
    eyy_A4new[i] = 0.0;
    //cout << yy[i] << endl;
  }
  TGraphErrors* fit_A4new   = new TGraphErrors(n_A4new,xx_A4new,yy_A4new,exx_A4new,eyy_A4new);
  fit_A4new->SetMarkerStyle(29);
  fit_A4new->SetMarkerSize(1.5);
  fit_A4new->SetMarkerColor(kBlack);
  fit_A4new->SetMinimum(0.0);
  fit_A4new->SetMaximum(110.0);
  //fit_A4new->Draw("apl");
  leg1->AddEntry(fit_A4new, "A_{4} scaled by 6/20", "P");
  mg->Add(fit_A4new);

  
  //////// 
  vector<TString> files_A4new36 = {
    "../root/fit_jacVsM_4G_y4p00_A4new_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_6_closure_jUL_j0_j1_j2_j3_j4_TEST_FULL.root",
    "../root/fit_jacVsM_4G_y4p00_A4new_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_6_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1.root",
    "../root/fit_jacVsM_4G_y4p00_A4new_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_6_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2.root",
    "../root/fit_jacVsM_4G_y4p00_A4new_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_6_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3.root",
  };
  vector<float> edges_A4new36 = {
    4.0-0.02,
    4.0-0.02,
    4.0-0.02,
    4.0-0.02
  };
  unsigned int n_A4new36 = files_A4new36.size();
  double xx_A4new36[10], yy_A4new36[10], yyMC_A4new36[10], exx_A4new36[10], eyy_A4new36[10];
  for(unsigned int i=0; i<n_A4new36; i++){
    TString fn = files_A4new36[i];
    TFile* f = TFile::Open(fn, "READ");
    TTree *t = f->Get<TTree>("tree");
    double err;
    t->SetBranchAddress("best_mass_err", &err);
    t->GetEntry(0);
    xx_A4new36[i] = edges_A4new36[i];
    exx_A4new36[i] = 0.0;
    yy_A4new36[i] = err;
    eyy_A4new36[i] = 0.0;
    //cout << yy[i] << endl;
  }
  TGraphErrors* fit_A4new36   = new TGraphErrors(n_A4new36,xx_A4new36,yy_A4new36,exx_A4new36,eyy_A4new36);
  fit_A4new36->SetMarkerStyle(kFullCross);
  fit_A4new36->SetMarkerSize(1.5);
  fit_A4new36->SetMarkerColor(kBlack);
  fit_A4new36->SetMinimum(0.0);
  fit_A4new36->SetMaximum(110.0);
  //fit_A4new36->Draw("apl");
  leg1->AddEntry(fit_A4new36, "A_{4} scaled by 6/20, A_{4}(y): 3->6", "P");
  mg->Add(fit_A4new36);

  //////// 
  vector<TString> files_8 = {
    "../root/fit_jacVsM_4G_y4p00_UL_10_8_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_FULL.root",
    "../root/fit_jacVsM_4G_y4p00_UL_10_8_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1.root",
    "../root/fit_jacVsM_4G_y4p00_UL_10_8_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2.root",
    "../root/fit_jacVsM_4G_y4p00_UL_10_8_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3.root",
  };
  vector<float> edges_8 = {
    4.0+0.02,
    4.0+0.02,
    4.0+0.02,
    4.0+0.02
  };
  unsigned int n_8 = files_8.size();
  double xx_8[10], yy_8[10], yyMC_8[10], exx_8[10], eyy_8[10];
  for(unsigned int i=0; i<n_8; i++){
    TString fn = files_8[i];
    TFile* f = TFile::Open(fn, "READ");
    TTree *t = f->Get<TTree>("tree");
    double err;
    t->SetBranchAddress("best_mass_err", &err);
    t->GetEntry(0);
    xx_8[i] = edges_8[i];
    exx_8[i] = 0.0;
    yy_8[i] = err;
    eyy_8[i] = 0.0;
    //cout << yy[i] << endl;
  }
  TGraphErrors* fit_8   = new TGraphErrors(n_8,xx_8,yy_8,exx_8,eyy_8);
  fit_8->SetMarkerStyle(22);
  fit_8->SetMarkerSize(1.5);
  fit_8->SetMarkerColor(kBlack);
  fit_8->SetMinimum(0.0);
  fit_8->SetMaximum(110.0);
  //fit_8->Draw("apl");
  leg1->AddEntry(fit_8, "corr(y): 6->8", "P");
  mg->Add(fit_8);


  //////// 
  vector<TString> files_AtanH = {
    "../root/fit_jacVsM_4G_y4p00_A4TanH_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_4_closure_jUL_j0_j1_j2_j3_j4_TEST_FULL.root",
    "../root/fit_jacVsM_4G_y4p00_A4TanH_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_4_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1.root",
    "../root/fit_jacVsM_4G_y4p00_A4TanH_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_4_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2.root",
    "../root/fit_jacVsM_4G_y4p00_A4TanH_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_4_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3.root",
  };
  vector<float> edges_AtanH = {
    4.0+0.04,
    4.0+0.04,
    4.0+0.04,
    4.0+0.04
  };
  unsigned int n_AtanH = files_AtanH.size();
  double xx_AtanH[10], yy_AtanH[10], yyMC_AtanH[10], exx_AtanH[10], eyy_AtanH[10];
  for(unsigned int i=0; i<n_AtanH; i++){
    TString fn = files_AtanH[i];
    TFile* f = TFile::Open(fn, "READ");
    TTree *t = f->Get<TTree>("tree");
    double err;
    t->SetBranchAddress("best_mass_err", &err);
    t->GetEntry(0);
    xx_AtanH[i] = edges_AtanH[i];
    exx_AtanH[i] = 0.0;
    yy_AtanH[i] = err;
    eyy_AtanH[i] = 0.0;
    //cout << yy[i] << endl;
  }
  TGraphErrors* fit_AtanH   = new TGraphErrors(n_AtanH,xx_AtanH,yy_AtanH,exx_AtanH,eyy_AtanH);
  fit_AtanH->SetMarkerStyle(23);
  fit_AtanH->SetMarkerSize(1.5);
  fit_AtanH->SetMarkerColor(kBlack);
  fit_AtanH->SetMinimum(0.0);
  fit_AtanH->SetMaximum(110.0);
  //fit_AtanH->Draw("apl");
  leg1->AddEntry(fit_AtanH, "A_{4}: pol -> ATanH", "P");
  mg->Add(fit_AtanH);

  //////// 
  vector<TString> files_AtanH_10G = {
    "../root/fit_jacVsM_10G_y4p00_A4TanH_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_4_closure_jUL_j0_j1_j2_j3_j4_TEST_FULL.root",
    "../root/fit_jacVsM_10G_y4p00_A4TanH_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_4_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1.root",
    "../root/fit_jacVsM_10G_y4p00_A4TanH_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_4_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2.root",
    "../root/fit_jacVsM_10G_y4p00_A4TanH_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_4_closure_jUL_j0_j1_j2_j3_j4_TEST_SCALEA0A1A2A3.root",
  };
  vector<float> edges_AtanH_10G = {
    4.0+0.06,
    4.0+0.06,
    4.0+0.06,
    4.0+0.06
  };
  unsigned int n_AtanH_10G = files_AtanH_10G.size();
  double xx_AtanH_10G[10], yy_AtanH_10G[10], yyMC_AtanH_10G[10], exx_AtanH_10G[10], eyy_AtanH_10G[10];
  for(unsigned int i=0; i<n_AtanH_10G; i++){
    TString fn = files_AtanH_10G[i];
    TFile* f = TFile::Open(fn, "READ");
    TTree *t = f->Get<TTree>("tree");
    double err;
    t->SetBranchAddress("best_mass_err", &err);
    t->GetEntry(0);
    xx_AtanH_10G[i] = edges_AtanH_10G[i];
    exx_AtanH_10G[i] = 0.0;
    yy_AtanH_10G[i] = err;
    eyy_AtanH_10G[i] = 0.0;
    //cout << yy[i] << endl;
  }
  TGraphErrors* fit_AtanH_10G   = new TGraphErrors(n_AtanH_10G,xx_AtanH_10G,yy_AtanH_10G,exx_AtanH_10G,eyy_AtanH_10G);
  fit_AtanH_10G->SetMarkerStyle(28);
  fit_AtanH_10G->SetMarkerSize(1.5);
  fit_AtanH_10G->SetMarkerColor(kBlack);
  fit_AtanH_10G->SetMinimum(0.0);
  fit_AtanH_10G->SetMaximum(110.0);
  //fit_AtanH_10G->Draw("apl");
  leg1->AddEntry(fit_AtanH_10G, "A_{4}: pol -> ATanH, 4G -> 10G", "P");
  mg->Add(fit_AtanH_10G);

  
  ////////////////////////////////////////////////////////
  mg->GetXaxis()->SetTitle("y_{max}");
  mg->GetYaxis()->SetTitle("#Delta M_{W} (MeV)");
  mg->GetYaxis()->SetRangeUser(0., 150.);
			     
  mg->Draw("apl");
  leg1->Draw();
  c->SaveAs("deltaM_vs_y.png");
  return;
}
