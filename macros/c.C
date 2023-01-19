{

  TFile* f1 = TFile::Open("../root/histos_NEWA3ZEROSMEARGW1p0_10M_UL_8_6_A0_8_6_A1_8_6_A2_8_6_A3_8_6_A4_8_6_grid.root", "READ");
  TH2D* h1 = f1->Get<TH2D>("jac_12");

  TFile* f2 = TFile::Open("../root/histos_NEWA3ZEROSMEARGW1p0_10G_UL_8_6_A0_8_6_A1_8_6_A2_8_6_A3_8_6_A4_8_6_grid.root", "READ");
  TH2D* h2 = f2->Get<TH2D>("jac_12");

  h1->Scale(1.);
  h2->Scale(1.);
  h1->Divide(h2);
  h1->SetMinimum(0.5);
  h1->SetMaximum(1.5);
  h1->Draw("colz");

}
