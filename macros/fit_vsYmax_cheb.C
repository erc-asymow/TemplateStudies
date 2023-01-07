{

  bool logy = false;

  float norm_lumi    = (0.100/0.800);
  TString image_name = "deltaM_vs_Ymax_paper2_cheb_GW1p0.eps";
  TString header     = "L = 1.0#times10^{8} events, p_{T}/GeV#in[25,p_{T}^{(H)}], |#eta|<2.5";

  //float norm_lumi    = (0.100/0.750);
  //TString image_name = "deltaM_vs_Ymax_paper2_cheb_GW1p0Y3p5X0p5.eps";
  //TString header     = "L = 1.0#times10^{8} events, p_{T}/GeV#in[25,p_{T}^{(H)}], |#eta|<2.5, |y|<3.5, x<0.5";

  TCanvas* c = new TCanvas("c", "canvas", 1200, 800);
  c->SetGridy();

  TMultiGraph* mg = new TMultiGraph();
  TLegend* leg1 = new TLegend(0.10, 0.68, 0.70, 0.90, "","brNDC");
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.035);
  leg1->SetFillColor(10);
  leg1->SetNColumns(4);
  leg1->SetHeader(header);

  vector<TString> files_n;
  vector<TString> files_l;
  vector<int>     files_c;
  vector<float>   files_w;
  files_n.emplace_back("UL_10_4_A0_3_2_A1_2_3_A2_3_2_A3_3_4_A4_3_3_full"); files_l.emplace_back("w_{cheb}");files_c.emplace_back(kRed); files_w.emplace_back(3); 
  
  vector<TString> posttags = {
    "jUL_j0_j1_j2_j3_j4_DEBUG_ADDMC_ULA0A1A2A3A4",
    //"jUL_j0_j1_j2_j3_j4_DEBUG_ULA0A1A2A3A4",
  };

  vector<TString> nevents_s;
  vector<float>   nevents_f;

  //nevents_s.emplace_back("Ymax70");  nevents_f.emplace_back(70);
  //nevents_s.emplace_back("Ymax68");  nevents_f.emplace_back(68);
  //nevents_s.emplace_back("Ymax66");  nevents_f.emplace_back(66);
  //nevents_s.emplace_back("Ymax64");  nevents_f.emplace_back(64);
  nevents_s.emplace_back("Ymax62");  nevents_f.emplace_back(62);
  nevents_s.emplace_back("Ymax60");  nevents_f.emplace_back(60);
  nevents_s.emplace_back("Ymax58");  nevents_f.emplace_back(58);
  nevents_s.emplace_back("Ymax56");  nevents_f.emplace_back(56);
  nevents_s.emplace_back("Ymax54");  nevents_f.emplace_back(54);
  nevents_s.emplace_back("Ymax52");  nevents_f.emplace_back(52);

  ////////////////////////////////////////////////
  for(unsigned int i1 = 0; i1<posttags.size(); i1++){
    TString posttag = posttags[i1];

    for(unsigned int i2 = 0; i2<files_n.size(); i2++){

      vector<float> edges_nevents = nevents_f;
      //for(auto& e : edges_nevents) e /= norm_lumi;
      if(logy){for(auto& e : edges_nevents) e = TMath::Log10(e);}

      unsigned int n_nevents = 0;
      double xx_nevents[20], yy_nevents[20], yyMC_nevents[20], exx_nevents[20], eyy_nevents[20];
      for(unsigned int i3=0; i3<nevents_s.size(); i3++){
	TString fn = "../root/fit_NEWA3ZEROSMEARGW1p0_10G_"+files_n[i2]+"_"+posttag+"_"+nevents_s[i3]+".root";
	TFile* f = TFile::Open(fn, "READ");
	if(f==0 || f==nullptr || f->IsZombie()) continue;
	TTree *t = f->Get<TTree>("tree");
	double err;
	t->SetBranchAddress("best_mass_err", &err);
	t->GetEntry(0);
	xx_nevents[n_nevents]  = edges_nevents[i3];
	exx_nevents[n_nevents] = 0.0;
	yy_nevents[n_nevents]  = err;
	eyy_nevents[n_nevents] = 0.0;
	n_nevents++;
	f->Close();
      }
      
    TGraphErrors* fit_nevents = new TGraphErrors(n_nevents,xx_nevents,yy_nevents,exx_nevents,eyy_nevents);
    if(i1==0){
      fit_nevents->SetLineStyle(kDashed);
      fit_nevents->SetMarkerStyle(kOpenCircle);
    }
    else{
      fit_nevents->SetLineStyle(kSolid);
      fit_nevents->SetMarkerStyle(kFullCircle);
    }
    if(files_w[i2]>2)
      fit_nevents->SetMarkerSize(1.5);
    else
      fit_nevents->SetMarkerSize(1.0);
      
    fit_nevents->SetMarkerColor(files_c[i2]);
    fit_nevents->SetLineColor(files_c[i2]);
    fit_nevents->SetLineWidth(files_w[i2]);
    fit_nevents->SetMinimum(0.0);
    fit_nevents->SetMaximum(55.0);
    TString leg_text = files_l[i2];
    if(i1==0){
      leg_text += " (w/ MC stat)";
    }
    if(i1==0) leg1->AddEntry(fit_nevents, leg_text, "PL");
    mg->Add(fit_nevents);
    }
  }
  ////////////////////////////////////////////////////////
  mg->GetXaxis()->SetTitle("p^{(H)}_{T} (GeV)");
  mg->GetXaxis()->SetTitleSize(0.045);
  mg->GetXaxis()->SetTitleOffset(0.9);
  mg->GetYaxis()->SetTitleSize(0.06);
  mg->GetYaxis()->SetTitleOffset(0.7);
  mg->GetYaxis()->SetTitle("#Delta M_{W} (MeV)");
  if(logy) mg->GetXaxis()->SetTitle("log_{10}(N/L)");
  mg->GetYaxis()->SetRangeUser(60, 120.);
  mg->GetXaxis()->SetRangeUser(0.0001, 100e+09);
  if(logy) mg->GetXaxis()->SetRangeUser(-4, 5);

  mg->GetYaxis()->SetMoreLogLabels();
  mg->GetYaxis()->CenterTitle(true);
  
  mg->Draw("apl");
  //c->SetLogy();
  leg1->Draw();
  c->SaveAs(image_name);
  //c->SaveAs("deltaM_vs_y_finalfix.png");

  //TF1* func = new TF1("func", "[0]*TMath::Log(x)+[1]", 0.0001, 50e+09);
  //fit_nevents->Fit(func);

  return;
}

