{

  TString plotname = "massloop_out_iter0";
  
  TCanvas* c = new TCanvas("c", "canvas", 1200, 400);
  c->Divide(3,1);

  TFile* fIter0 = TFile::Open("massfit_SmearRealistic_Iter0.root", "READ");
  TFile* fIter1 = TFile::Open("massfit_SmearRealistic_Iter1.root", "READ");

  //vector<TString> params = {"Ain", "ein", "Min"};
  vector<TString> params = {"A", "e", "M"};

  for(unsigned int p = 0; p < params.size(); p++){
    TH1D* hp_nom  = (TH1D*) fIter0->Get("h_"+params[p]+"_vals_nom");
    TH1D* hp_fit0 = (TH1D*) fIter0->Get("h_"+params[p]+"_vals_fit");
    TH1D* hp_fit1 = (TH1D*) fIter1->Get("h_"+params[p]+"_vals_fit");
    //hp_fit0->Add(hp_fit1);
    hp_nom->SetLineColor(kBlue);
    hp_nom->SetLineWidth(3);
    hp_fit0->SetLineColor(kBlack);
    hp_fit0->SetMarkerColor(kBlack);
    hp_fit0->SetMarkerStyle(kFullCircle);
    c->cd(p+1);
    hp_nom->SetStats(0);
    hp_nom->SetMaximum( hp_nom->GetMaximum()*1.2);
    hp_nom->SetMinimum( hp_nom->GetMinimum()*0.8);
    hp_fit0->SetStats(0);
    hp_fit0->SetMaximum( hp_fit0->GetMaximum()*1.2);
    hp_fit0->SetMinimum( hp_fit0->GetMinimum()*0.8);

    if(p>=1){
      hp_fit0->SetMaximum( + (hp_fit0->GetMaximum()+hp_fit0->GetBinError(1))*1.1 );
      hp_fit0->SetMinimum( - (hp_fit0->GetMaximum()+hp_fit0->GetBinError(1))*1.1 );
    }
      
    hp_fit0->Draw("HISTPE");
    hp_nom->Draw("HISTSAME");
  }

  c->Update();
  c->Draw();
  c->SaveAs(plotname+".png");
}
