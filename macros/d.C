{

  TFile *f1 = TFile::Open("root/histos_jacVsM_4G_y4p50_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure.root");
  TH2D* h1 = f1->Get<TH2D>("jac_0");

  TFile *f2 = TFile::Open("root/histos_jacVsM_4G_y2p75_UL_10_6_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure.root");
  TH2D* h2 = f2->Get<TH2D>("jac_0");

  h2->Divide(h1);
  h2->SetMinimum(2.3);
  h2->SetMaximum(2.5);
  h2->Draw("colz");
  h2->Print("all");  
}
