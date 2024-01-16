

{

  
  TCanvas* c = new TCanvas("c", "canvas", 1200, 800);
  //c->SetGridy();

  TLegend* leg1 = new TLegend(0.50, 0.66, 0.75, 0.86, "","brNDC");
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.045);
  leg1->SetFillColor(10);
  leg1->SetNColumns(2);
  leg1->SetHeader("N = L = 10^{6}, p = 6");

  
  TFile *f = TFile::Open("root/toy_mctstat2_5_closure.root");

  TH1D* h = (TH1D*)f->Get("h");
  h->SetMinimum(7000.);
  h->SetMaximum(13000.);
  h->SetTitle("");
  h->SetStats(0);
  h->GetYaxis()->SetMaxDigits(3);
  h->GetXaxis()->SetTitle("x");
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetNdivisions(202, kFALSE);
  h->GetXaxis()->SetTitleOffset(0.82);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleOffset(0.7);
  h->GetYaxis()->SetTitle("Events/bin");
  //h->GetXaxis()->SetLimits(0.2, 50);
  //h->GetXaxis()->SetMoreLogLabels();
  //h->GetYaxis()->CenterTitle(true);   
  h->SetLineWidth(2);
  h->SetMarkerStyle(kFullCircle);
  h->SetMarkerSize(1.5);
  leg1->AddEntry(h, "#mu = 0.50", "P");
  h->Draw("PE");
    
  TH1D* h_asy = (TH1D*)f->Get("h_asy");
  h_asy->SetLineWidth(3);
  h_asy->SetLineColor(kBlack);
  leg1->AddEntry(h_asy, "True", "L");
  h_asy->Draw("HISTSAME");

  TH1D* h_mass3 = (TH1D*)f->Get("h_mass2");
  h_mass3->SetLineWidth(2);
  h_mass3->SetLineColor(kRed);
  h_mass3->SetMarkerStyle(kFullSquare);
  h_mass3->SetMarkerSize(1.5);
  h_mass3->SetMarkerColor(kRed);  
  leg1->AddEntry(h_mass3, "#mu = 0.45", "P");
  h_mass3->Draw("PESAME");

  TH1D* h_mass3_asy = (TH1D*)f->Get("h_mass2_asy");
  h_mass3_asy->SetLineWidth(3);
  h_mass3_asy->SetLineColor(kRed);
  leg1->AddEntry(h_mass3_asy, "True", "L");
  h_mass3_asy->Draw("HISTSAME");

  TH1D* h_mass0 = (TH1D*)f->Get("h_mass0");
  h_mass0->SetLineWidth(2);
  h_mass0->SetLineColor(kBlue);
  h_mass0->SetMarkerStyle(kFullTriangleUp);
  h_mass0->SetMarkerSize(1.5);
  h_mass0->SetMarkerColor(kBlue);  
  leg1->AddEntry(h_mass0, "#mu = 0.40", "P");
  h_mass0->Draw("PESAME");

  TH1D* h_mass0_asy = (TH1D*)f->Get("h_mass0_asy");
  h_mass0_asy->SetLineWidth(3);
  h_mass0_asy->SetLineColor(kBlue);
  leg1->AddEntry(h_mass0_asy, "True", "L");
  h_mass0_asy->Draw("HISTSAME");

  c->cd();
  
  leg1->Draw();

  c->Modified();
  c->Update();
  //c->SetLogx();
  //c->SetLogy();

  TString plotname = "toy_paper_y";
  c->SaveAs(plotname+".eps");
  c->SaveAs(plotname+".pdf");
  c->SaveAs(plotname+".png");
  
}
