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

  double max_y = 3.5;

  TF1* tf1toy_y = new TF1("toy_y", "[2]/TMath::Sqrt(2*TMath::Pi()*[1])*TMath::Exp(-0.5*(x-[0])*(x-[0])/[1])", -max_y, max_y);
  double sigma2_y = 4.0*4.0;
  tf1toy_y->SetParameter(0, 0.0);
  tf1toy_y->SetParameter(1, sigma2_y);
  tf1toy_y->SetParameter(2, 1.0);
  double int_toy_y = tf1toy_y->Integral(-max_y, max_y);
  tf1toy_y->SetParameter(2, 1.0/int_toy_y);

  tf1toy_y->SetRange(0., max_y);
  tf1toy_y->SetLineColor(kBlack);
  tf1toy_y->SetLineWidth(3);
  tf1toy_y->Draw();
  tf1toy_y->GetYaxis()->SetRangeUser(0., 0.2);  
  tf1toy_y->SetTitle("");
  tf1toy_y->GetXaxis()->SetTitleSize(0.045);
  tf1toy_y->GetXaxis()->SetTitleOffset(0.9);
  tf1toy_y->GetYaxis()->SetTitleSize(0.06);
  tf1toy_y->GetYaxis()->SetTitleOffset(0.6);
  tf1toy_y->GetXaxis()->SetTitle("|y|");
  tf1toy_y->GetYaxis()->SetTitle("w_{y}");
  leg->AddEntry(tf1toy_y, "nominal", "L");
  
  TF1* tf1toy_y2 = new TF1("toy_y2", "[2]/TMath::Sqrt(2*TMath::Pi()*[1])*TMath::Exp(-0.5*(x-[0])*(x-[0])/[1])", -max_y, max_y);
  sigma2_y = 3.0*3.0;
  tf1toy_y2->SetParameter(0, 0.0);
  tf1toy_y2->SetParameter(1, sigma2_y);
  tf1toy_y2->SetParameter(2, 1.0);
  int_toy_y = tf1toy_y2->Integral(-max_y, max_y);
  tf1toy_y2->SetParameter(2, 1.0/int_toy_y);
  leg->AddEntry(tf1toy_y2, "alternate", "L");
  tf1toy_y2->SetLineColor(kRed);
  tf1toy_y2->SetLineWidth(3);
  tf1toy_y2->SetLineStyle(kDashed);
  tf1toy_y2->Draw("same");
  
  leg->Draw();  
  gPad->SaveAs("wy.eps");
}
