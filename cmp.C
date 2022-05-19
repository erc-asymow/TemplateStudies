
{

  int axis = 1;
  
  bool normalize = true;

  TString name = "testy";

  std::vector<TString> tags{};

  std::vector<TFile*> files{};

  std::vector<TH1D*> histos{};
  std::vector<TH1D*> histos_mass{};

  for(auto X : ROOT::TSeqI(16) ){
    if(X<8) continue;
    for(auto Y : ROOT::TSeqI(13) ){
      if(Y>8) continue;
      for(auto CORRX : ROOT::TSeqI(13) ){
	if(CORRX<2) continue;
	for(auto CORRY : ROOT::TSeqI(13) ){
	  if(CORRY>8) continue;
	  TString tag(Form("%d_%d_%d_%d",X,Y,CORRX,CORRY));
	  cout << "file " << tag << endl;
	  TFile* f = TFile::Open("root/histos_"+name+"_"+tag+"_closure.root");
	  if(f==nullptr || f->IsZombie()) continue;
	  files.emplace_back(f);    
	  tags.emplace_back(tag);
	}
      }
    }
  }

  TCanvas* c = new TCanvas();
  TLegend* leg = new TLegend(0.1,0.7,0.4,0.9);
  c->cd();

  int i = 0;
  for(auto tag : tags){
    TH2D* h_truth = files[i]->Get<TH2D>("wMC");
    TH1D* hT = axis==0?h_truth->ProjectionX("hT_py"):h_truth->ProjectionY("hT_py");
    TH2D* h_cheb = files[i]->Get<TH2D>("w");
    TH1D* hC = axis==0?h_cheb->ProjectionX("hC_py"):h_cheb->ProjectionY("hC_py");

    if(i==0){
      files[i]->cd();
      if(files[i]->FindObjectAny("wMC_up")){
	TH2D* h_up = files[i]->Get<TH2D>("wMC_up");
	TH1D* hMup =  axis==0?h_up->ProjectionX("hMup_py"):h_up->ProjectionY("hMup_py");
	hMup->Divide(hT);
	histos_mass.emplace_back((TH1D*)hMup->Clone("hMUp"));
	cout << tag << " pushed back" << endl;
      }
      if(files[i]->FindObjectAny("wMC_down")){
	TH2D* h_down = files[i]->Get<TH2D>("wMC_down");
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
      h->SetMinimum(0.998);
      h->SetMaximum(1.002);
      h->Draw("HIST");
    }
    else h->Draw("HISTSAME");
    leg->AddEntry(h, tags[i], "L");
    i++;
  }

  i=0;
  for(auto h : histos_mass){
    h->SetLineWidth(3);
    h->Draw("HISTSAME");
    if(i==0) leg->AddEntry(h, "#pm 10 MeV", "L");
    i++;
  }

  leg->Draw()
}
