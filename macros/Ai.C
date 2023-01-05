{

  TCanvas* c = new TCanvas("c", "canvas", 800, 800);
  c->Divide(2,2);

  //gStyle->SetPalette(kOcean);
  
  TLegend* leg = new TLegend(0.10, 0.70, 0.30, 0.90, "","brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.06);
  leg->SetFillColor(10);
  leg->SetNColumns(1);
  leg->SetHeader("");

  double max_x = 0.5;
  double max_y = 3.5;
  
  TF2* tf2toy_A0 = new TF2("toy_A0", "2*x*x*(1 - 0.01*y*y)",                    0., max_x, 0, max_y);  
  TF2* tf2toy_A1 = new TF2("toy_A1", "(0.5*x + 2*x*x)*( 0.05*y + 0.002*y*y*y)", 0., max_x, 0, max_y);
  TF2* tf2toy_A2 = new TF2("toy_A2", "2*x*x*(1 - 0.01*y*y)",                    0., max_x, 0, max_y);
  TF2* tf2toy_A3 = new TF2("toy_A3", "0.3*(x + x*x + x*x*x)*(0.1*y*y)",         0., max_x, 0, max_y);
  TF2* tf2toy_A4 = new TF2("toy_A4", "(1-x)*(y + y*y*y/10.)/5",                 0., max_x, 0, max_y);

  c->cd(1);
  gPad->SetRightMargin(0.2);
  tf2toy_A0->SetNpx(10000);
  tf2toy_A0->Draw("colz");
  tf2toy_A0->SetTitle("A_{0,2}");
  tf2toy_A0->GetXaxis()->SetTitleSize(0.06);
  tf2toy_A0->GetXaxis()->SetTitleOffset(0.8);
  tf2toy_A0->GetYaxis()->SetTitleSize(0.06);
  tf2toy_A0->GetYaxis()->SetTitleOffset(0.6);
  tf2toy_A0->GetXaxis()->SetTitle("x");
  tf2toy_A0->GetYaxis()->SetTitle("y");

  c->cd(2);
  gPad->SetRightMargin(0.2);
  tf2toy_A1->SetNpx(10000);
  tf2toy_A1->Draw("colz");
  tf2toy_A1->SetTitle("A_{1}");
  tf2toy_A1->GetXaxis()->SetTitleSize(0.06);
  tf2toy_A1->GetXaxis()->SetTitleOffset(0.8);
  tf2toy_A1->GetYaxis()->SetTitleSize(0.06);
  tf2toy_A1->GetYaxis()->SetTitleOffset(0.6);
  tf2toy_A1->GetXaxis()->SetTitle("x");
  tf2toy_A1->GetYaxis()->SetTitle("y");

  c->cd(3);
  gPad->SetRightMargin(0.2);
  tf2toy_A3->SetNpx(10000);
  tf2toy_A3->Draw("colz");
  tf2toy_A3->SetTitle("A_{3}");
  tf2toy_A3->GetXaxis()->SetTitleSize(0.06);
  tf2toy_A3->GetXaxis()->SetTitleOffset(0.8);
  tf2toy_A3->GetYaxis()->SetTitleSize(0.06);
  tf2toy_A3->GetYaxis()->SetTitleOffset(0.6);
  tf2toy_A3->GetXaxis()->SetTitle("x");
  tf2toy_A3->GetYaxis()->SetTitle("y");

  c->cd(4);
  gPad->SetRightMargin(0.2);
  tf2toy_A4->SetNpx(10000);
  tf2toy_A4->Draw("colz");
  tf2toy_A4->SetTitle("A_{4}");
  tf2toy_A4->GetXaxis()->SetTitleSize(0.06);
  tf2toy_A4->GetXaxis()->SetTitleOffset(0.8);
  tf2toy_A4->GetYaxis()->SetTitleSize(0.06);
  tf2toy_A4->GetYaxis()->SetTitleOffset(0.6);
  tf2toy_A4->GetXaxis()->SetTitle("x");
  tf2toy_A4->GetYaxis()->SetTitle("y");
  
  leg->Draw();
  c->cd();
  gPad->SaveAs("Ai.eps");
}
