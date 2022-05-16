
{

  bool normalize = true;

  TString name = "test_corrxy";

  std::vector<TString> tags{};

  std::vector<TFile*> files{};

  std::vector<TH1D*> histos{};
  std::vector<TH1D*> histos_mass{};

  for(auto X : ROOT::TSeqI(11) ){
    for(auto Y : ROOT::TSeqI(11) ){
      if(X<6 || X>8) continue;
      if(Y<6 || Y>8) continue;
      TString tag(Form("8_6_%d_%d",X,Y));
      cout << "file " << tag << endl;
      TFile* f = TFile::Open("root/histos_"+name+"_"+tag+"_closure.root");
      if(f==nullptr || f->IsZombie()) continue;
      files.emplace_back(f);    
      tags.emplace_back(tag);
    }
  }

  TCanvas* c = new TCanvas();
  TLegend* leg = new TLegend();
  c->cd();

  int i = 0;
  for(auto tag : tags){
    TH2D* h_truth = files[i]->Get<TH2D>("wMC");
    TH1D* hT = h_truth->ProjectionY("hT_py");
    TH2D* h_cheb = files[i]->Get<TH2D>("w");
    TH1D* hC = h_cheb->ProjectionY("hC_py");

    if(i==0){
      files[i]->cd();
      if(files[i]->FindObjectAny("wMC_up")){
	TH2D* h_up = files[i]->Get<TH2D>("wMC_up");
	TH1D* hMup =  h_up->ProjectionY("hMup_py");
	hMup->Divide(hT);
	histos_mass.emplace_back((TH1D*)hMup->Clone("hMUp"));
	cout << tag << " pushed back" << endl;
      }
      if(files[i]->FindObjectAny("wMC_down")){
	TH2D* h_down = files[i]->Get<TH2D>("wMC_down");
	TH1D* hMdown =  h_down->ProjectionY("hMdown_py");
	hMdown->Divide(hT);
	histos_mass.emplace_back((TH1D*)hMdown->Clone("hMDown"));
	cout << tag << " pushed back" << endl;
      }
    }

    if(normalize){
      hC->Scale(1./hC->Integral());
      hT->Scale(1./hT->Integral());
    }
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
      h->SetMinimum(0.995);
      h->SetMaximum(1.005);
      h->Draw("HIST");
    }
    else h->Draw("HISTSAME");
    leg->AddEntry(h, tags[i], "L");
    i++;
  }

  for(auto h : histos_mass){
    h->SetLineWidth(3);
    h->Draw("HISTSAME");
    leg->AddEntry(h, "#pm 10 MeV", "L");
  }

  leg->Draw()
}
