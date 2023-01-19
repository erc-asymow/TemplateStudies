{

  bool logy = true;
  float norm_lumi = 0.1;
    
  TCanvas* c = new TCanvas("c", "canvas", 1000, 600);
  c->SetGridy();

  TMultiGraph* mg = new TMultiGraph();
  TLegend* leg1 = new TLegend(0.25, 0.75, 0.55, 0.90, "","brNDC");
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetFillColor(10);
  //leg1->SetHeader("");
  leg1->SetNColumns(3);

  TString common_tag = "DEBUG_ADDMC";
  
  vector<TString> posttags = {
    "jUL_j0_j1_j2_j3_j4_"+common_tag+"_ULA0A1A2A3A4",
    //"jUL_j0_j1_j2_j3_j4_"+common_tag+"_ULA0A1A2A3A4",
    //"jUL_j0_j1_j2_j3_"+common_tag+"_ULA0A1A2A3",
    //"jUL_j0_j1_j2_"+common_tag+"_ULA0A1A2",
    //"jUL_j0_j1_"+common_tag+"_ULA0A1",
    //"jUL_j0_"+common_tag+"_ULA0",
    //"jUL_"+common_tag+"_UL"
  };

  vector<TString> legends = {
    "UL,A_{0,1,2,3,4}",
    //"UL,A_{0,1,2,3}",
    //"UL,A_{0,1,2}",
    //"UL,A_{0,1}",
    //"UL,A_{0}",
    //"UL",
  };

  for(unsigned int i = 0; i<posttags.size(); i++){
    TString posttag = posttags[i];
    
    ////// nevents
    vector<TString> files_nevents = {
      "../root/fit_DEV_1M_UL_3_2_A0_1_2_A1_1_1_A2_1_2_A3_1_2_A4_1_1_closure_"+posttag+".root",
      "../root/fit_DEV_10M_UL_3_2_A0_1_2_A1_1_1_A2_1_2_A3_1_2_A4_1_1_closure_"+posttag+".root",
      "../root/fit_DEV_100M_UL_3_2_A0_1_2_A1_1_1_A2_1_2_A3_1_2_A4_1_1_closure_"+posttag+".root",
      "../root/fit_DEV_1G_UL_3_2_A0_1_2_A1_1_1_A2_1_2_A3_1_2_A4_1_1_closure_"+posttag+".root",
      "../root/fit_DEV_4G_UL_3_2_A0_1_2_A1_1_1_A2_1_2_A3_1_2_A4_1_1_closure_"+posttag+".root",
      "../root/fit_DEV_10G_UL_3_2_A0_1_2_A1_1_1_A2_1_2_A3_1_2_A4_1_1_closure_"+posttag+".root",
      "../root/fit_DEV_40G_UL_3_2_A0_1_2_A1_1_1_A2_1_2_A3_1_2_A4_1_1_closure_"+posttag+".root",
      "../root/fit_DEV_80G_UL_3_2_A0_1_2_A1_1_1_A2_1_2_A3_1_2_A4_1_1_closure_"+posttag+".root",
      //"../root/fit_DEV_1M_UL_3_2_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_closure_"+posttag+".root",
      //"../root/fit_DEV_10M_UL_3_2_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_closure_"+posttag+".root",
      //"../root/fit_DEV_100M_UL_3_2_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_closure_"+posttag+".root",
      //"../root/fit_DEV_1G_UL_3_2_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_closure_"+posttag+".root",
      //"../root/fit_DEV_4G_UL_3_2_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_closure_"+posttag+".root",
      //"../root/fit_DEV_10G_UL_3_2_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_closure_"+posttag+".root",
      //"../root/fit_DEV_40G_UL_3_2_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_closure_"+posttag+".root",
      //"../root/fit_DEV_80G_UL_3_2_A0_1_1_A1_1_1_A2_1_1_A3_1_1_A4_1_1_closure_"+posttag+".root",
      //"../root/fit_DEV_1M_UL_3_4_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_3_closure_"+posttag+".root",
      //"../root/fit_DEV_10M_UL_3_4_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_3_closure_"+posttag+".root",
      //"../root/fit_DEV_100M_UL_3_4_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_3_closure_"+posttag+".root",
      //"../root/fit_DEV_1G_UL_3_4_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_3_closure_"+posttag+".root",
      //"../root/fit_DEV_10G_UL_3_4_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_3_closure_"+posttag+".root",
      //"../root/fit_DEV_40G_UL_3_4_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_3_closure_"+posttag+".root",
      //"../root/fit_DEV_80G_UL_3_4_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_3_closure_"+posttag+".root",
      //"../root/fit_DEV_200G_UL_3_4_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_3_closure_"+posttag+".root",
    };
    vector<float> edges_nevents = {
      0.001,
      0.01,
      0.1,
      1.,
      4.,
      10.,
      40.,
      80.,
      //200,
    };
    
    for(auto& e : edges_nevents) e *= 0.710/norm_lumi;
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
      //cout << yy[i] << endl;
    }
    TGraphErrors* fit_nevents   = new TGraphErrors(n_nevents,xx_nevents,yy_nevents,exx_nevents,eyy_nevents);
    fit_nevents->SetMarkerStyle(kFullCircle);
    fit_nevents->SetMarkerSize(1.5);
    fit_nevents->SetMarkerColor(i+1);
    fit_nevents->SetMinimum(0.0);
    fit_nevents->SetMaximum(55.0);
    leg1->AddEntry(fit_nevents, legends[i], "P");
    mg->Add(fit_nevents);
    
  }
  ////////////////////////////////////////////////////////
  mg->GetXaxis()->SetTitle("(N effective)/10^{9} events");
  mg->GetYaxis()->SetTitle("#Delta M_{W} (MeV)");
  if(logy) mg->GetXaxis()->SetTitle("log_{10}(N/norm)");
  mg->GetYaxis()->SetRangeUser(0.0, 75.);
  mg->GetXaxis()->SetRangeUser(0.0001, 100e+09);
  if(logy) mg->GetXaxis()->SetRangeUser(-4, 3);
			     
  mg->Draw("apl");
  leg1->Draw();
  c->SaveAs("deltaM_vs_y_v2.png");

  //TF1* func = new TF1("func", "[0]*TMath::Log(x)+[1]", 0.0001, 50e+09);
  //fit_nevents->Fit(func);

  return;
}
