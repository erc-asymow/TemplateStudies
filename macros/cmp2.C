
{

  int fit_qt_y = 0;
  
  int axis = 1;
  
  bool normalize = true;

  std::vector<TFile*> files{};
  std::vector<TString> tags{};

  //files.emplace_back( TFile::Open("../root/histos_NEWA3ZEROSMEAR_cheb_4G_UL_8_4_A0_3_2_A1_2_3_A2_3_2_A3_3_4_A4_3_3_full.root") );
  //tags.emplace_back("n_{x}=8");
  //files.emplace_back( TFile::Open("../root/histos_NEWA3ZEROSMEAR_cheb_4G_UL_10_4_A0_3_2_A1_2_3_A2_3_2_A3_3_4_A4_3_3_full.root") );
  //tags.emplace_back("n_{x}=10");
  //files.emplace_back( TFile::Open("../root/histos_NEWA3ZEROSMEAR_cheb_4G_UL_12_4_A0_3_2_A1_2_3_A2_3_2_A3_3_4_A4_3_3_full.root") );
  //tags.emplace_back("n_{x}=12");
  //files.emplace_back(TFile::Open( "../root/histos_NEWA3ZEROSMEAR_cheb_4G_UL_10_6_A0_3_2_A1_2_3_A2_3_2_A3_3_4_A4_3_3_full.root") );
  //tags.emplace_back("n_{y}=6");

  //files.emplace_back( TFile::Open("../root/histos_NEWA3ZEROSMEARY3p5X0p5Ymax70GW1p0_40G_UL_3_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr.root"));
  //tags.emplace_back("1p0");
  //files.emplace_back( TFile::Open("../root/histos_NEWA3ZEROSMEARY3p5X0p5Ymax70_40G_UL_3_2_A0_2_2_A1_2_2_A2_2_2_A3_2_2_A4_2_2_corr.root"));
  //tags.emplace_back("2p0");

  //files.emplace_back( TFile::Open("../root/histos_NEWA3ZEROSMEARGW1p0_10G_UL_8_6_A0_8_6_A1_8_6_A2_8_6_A3_8_6_A4_8_6_grid.root"));
  //tags.emplace_back("1p0");
  //files.emplace_back( TFile::Open("../root/histos_NEWA3ZEROSMEAR_10G_UL_8_6_A0_8_6_A1_8_6_A2_8_6_A3_8_6_A4_8_6_grid.root"));
  //tags.emplace_back("2p0");

  //files.emplace_back( TFile::Open("../root/histos_NEWA3ZEROSMEARGW1p0Y3p5X0p5_10G_UL_12_6_A0_3_2_A1_2_3_A2_3_2_A3_3_4_A4_3_3_full.root"));
  //tags.emplace_back("0p5");
  files.emplace_back( TFile::Open("../root/histos_NEWA3ZEROSMEARGW1p0Y3p5_4G_UL_10_6_A0_3_2_A1_2_3_A2_3_2_A3_3_4_A4_3_3_full.root"));
  tags.emplace_back("0p5");
  //files.emplace_back( TFile::Open("../root/histos_NEWA3ZEROSMEARGW1p0_10G_UL_10_4_A0_3_2_A1_2_3_A2_3_2_A3_3_4_A4_3_3_full.root"));
  //tags.emplace_back("0p5");

  std::vector<TH1D*> histos{};
  std::vector<TH1D*> histos_mass{};

  TCanvas* c = new TCanvas("c", "canvas", 1200, 800);
  //c->SetGridy();

  TLegend* leg = new TLegend(0.10, 0.68, 0.70, 0.90, "","brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetFillColor(10);
  leg->SetNColumns(1);
  leg->SetHeader("");
  c->cd();

  int i = 0;
  for(auto tag : tags){
    TH2D* h_truth = files[i]->Get<TH2D>("hMC");
    TH1D* hT = axis==0?h_truth->ProjectionX("hT_py"):h_truth->ProjectionY("hT_py");
    TH2D* h_cheb = files[i]->Get<TH2D>("h");
    TH1D* hC = axis==0?h_cheb->ProjectionX("hC_py"):h_cheb->ProjectionY("hC_py");
    if(i==0){
      files[i]->cd();
      if(files[i]->FindObjectAny("hMC_up")){
	TH2D* h_up = files[i]->Get<TH2D>("hMC_up");
	TH1D* hMup =  axis==0?h_up->ProjectionX("hMup_py"):h_up->ProjectionY("hMup_py");
	hMup->Divide(hT);
	histos_mass.emplace_back((TH1D*)hMup->Clone("hMUp"));
	cout << tag << " pushed back" << endl;
      }
      if(files[i]->FindObjectAny("hMC_down")){
	TH2D* h_down = files[i]->Get<TH2D>("hMC_down");
	TH1D* hMdown =  axis==0?h_down->ProjectionX("hMdown_py"):h_down->ProjectionY("hMdown_py");
	hMdown->Divide(hT);
	histos_mass.emplace_back((TH1D*)hMdown->Clone("hMDown"));
	cout << tag << " pushed back" << endl;
      }
    }

    if(normalize){
      hC->Scale(1./hC->Integral());
      hT->Scale(1./hT->Integral());
    }
    hT->Sumw2();
    hC->Sumw2();
    hC->Divide(hT);
    hC->SetLineColor(i+1);
    histos.emplace_back((TH1D*)hC->Clone(tag));    
    cout << tag << " pushed back" << endl;
    i++;
  }

  c->cd();
  i = 0;
  for(auto h : histos){
    if(i==0){
      if(fit_qt_y){
	h->SetMinimum(0.99);
	h->SetMaximum(1.01);
      }
      else{
	h->SetMinimum(0.99);
	h->SetMaximum(1.01);
      }
      h->SetTitle("");
      h->GetXaxis()->SetTitleSize(0.045);
      h->GetXaxis()->SetTitleOffset(0.9);
      h->GetYaxis()->SetTitleSize(0.06);
      h->GetYaxis()->SetTitleOffset(0.8);
      h->GetXaxis()->SetTitle("p_{T} [GeV]");
      h->GetYaxis()->SetTitle("w_{cheb}/w_{toy}");
      h->SetStats(0);
      h->Draw("HIST");
    }
    else h->Draw("HISTSAME");
    h->SetLineWidth(2);
    leg->AddEntry(h, tags[i], "L");
    i++;
  }

  i=0;
  for(auto h : histos_mass){
    h->SetLineWidth(3);
    h->SetLineStyle(kDashed);
    h->Draw("HISTSAME");
    if(i==0) leg->AddEntry(h, "#pm 10 MeV", "L");
    i++;
  }
  leg->Draw();
  //gPad->SaveAs("plot_cmp_1p0.png");
}
