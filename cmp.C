
{

  TString name = "toy2";
  std::vector<TString> tags{};

  std::vector<TFile*> files{};

  for(auto X : ROOT::TSeqI(11) ){
    for(auto Y : ROOT::TSeqI(11) ){
      if(X<8) continue;
      if(Y<6) continue;
      TString tag(Form("%d_%d",X,Y));
      TFile* f = TFile::Open("root/histos_"+name+"_"+tag+"_closure.root");
      if(f==nullptr || f->IsZombie()) continue;
      files.emplace_back(f);    
      tags.emplace_back(tag);
    }
  }

  TCanvas* c = new TCanvas();
  TLegend* leg = new TLegend();
  c->cd();
  std::vector<TH1D*> histos{};

  int i = 0;
  for(auto tag : tags){
    TH2D* h_truth = files[i]->Get<TH2D>("wMC");
    TH1D* hT = h_truth->ProjectionY("hT_py");
    TH2D* h_cheb = files[i]->Get<TH2D>("w");
    TH1D* hC = h_cheb->ProjectionY("hC_py");
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

  leg->Draw()
}
