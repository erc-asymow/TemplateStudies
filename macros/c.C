{

  TFile* f1 = TFile::Open("root/histos_addmass7_seed0_UL_10_4_grid.root", "READ");
  TH2D* h1 = f1->Get<TH2D>("jac_100");

  TFile* f2 = TFile::Open("root/histos_addmass7_merged_UL_10_4_grid.root", "READ");
  TH2D* h2 = f2->Get<TH2D>("jac_100");

  h1->Divide(h2);
  h1->SetMinimum(0.99);
  h1->SetMaximum(1.01);
  h1->Draw("colz");

}
