{

  TCanvas* c = new TCanvas("c", "canvas", 1200, 800);
  //c->SetGridy();

  TLegend* leg = new TLegend(0.50, 0.70, 0.80, 0.90, "","brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.06);
  leg->SetFillColor(10);
  leg->SetNColumns(1);
  leg->SetHeader("");
  c->cd();

  double max_x = 0.5;

  TF1* tf1toy_x = new TF1("toy_x", "[0]*x/TMath::Power(x*x+[1], 1.0)", 0.0, max_x);  
  tf1toy_x->SetNpx(10000);
  double p0_x = +2.35e-03;
  tf1toy_x->SetParameter(0, 1.0);
  tf1toy_x->SetParameter(1, p0_x);
  double int_toy_x = tf1toy_x->Integral(0.0, max_x);
  tf1toy_x->SetParameter(0, 1.0/int_toy_x);
  
  tf1toy_x->SetRange(0., max_x);
  tf1toy_x->SetLineColor(kBlack);
  tf1toy_x->SetLineWidth(3);
  tf1toy_x->Draw();
  tf1toy_x->GetYaxis()->SetRangeUser(0., 8.);  
  tf1toy_x->SetTitle("");
  tf1toy_x->GetXaxis()->SetTitleSize(0.045);
  tf1toy_x->GetXaxis()->SetTitleOffset(0.9);
  tf1toy_x->GetYaxis()->SetTitleSize(0.06);
  tf1toy_x->GetYaxis()->SetTitleOffset(0.6);
  tf1toy_x->GetXaxis()->SetTitle("x");
  tf1toy_x->GetYaxis()->SetTitle("w_{x}");
  leg->AddEntry(tf1toy_x, "nominal", "L");

  TF1* tf1toy_x2 = new TF1("toy_x", "[0]*x/TMath::Power(x*x+[1], 1.25)", 0.0, max_x);  
  tf1toy_x2->SetNpx(10000);
  tf1toy_x2->SetParameter(0, 1.0);
  tf1toy_x2->SetParameter(1, p0_x);
  int_toy_x = tf1toy_x2->Integral(0.0, max_x);
  tf1toy_x2->SetParameter(0, 1.0/int_toy_x);
  leg->AddEntry(tf1toy_x2, "alternate", "L");
  tf1toy_x2->SetLineColor(kRed);
  tf1toy_x2->SetLineWidth(3);
  tf1toy_x2->SetLineStyle(kDashed);
  tf1toy_x2->Draw("same");
  
  leg->Draw();  
  gPad->SaveAs("wxy.eps");
}
