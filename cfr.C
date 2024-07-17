
int cfr(int ibin=0, bool rms=false){
  
  TFile* f0 = TFile::Open("massscales_PostVFPGausPdf_Iter0.root", "READ");
  TFile* f1 = TFile::Open("massscales_PostVFPCBPdf_Iter0.root", "READ");

  TH2D* h0 = (TH2D*)f0->Get( rms ? "h_smear0_bin_jac_width" : "h_smear0_bin_jac_scale");
  TH2D* h1 = (TH2D*)f1->Get( rms ? "h_smear0_bin_jac_width" : "h_smear0_bin_jac_scale");

  TH1D* s0 = h0->ProjectionY( Form("s0%d_%d", ibin, rms ? 1 : 0), ibin+1, ibin+1);
  TH1D* s1 = h1->ProjectionY( Form("s1%d_%d", ibin, rms ? 1 : 0), ibin+1, ibin+1);

  s1->SetLineColor(kRed);

  s0->SetStats(0);

  //s0->SetMaximum( TMath::Max( s0->GetMaximum(), s1->GetMaximum())*1.1 );
  //s0->SetMinimum( TMath::Min( s0->GetMinimum(), s1->GetMinimum())*1.1 );
  s0->Draw("HISTE");
  s1->Draw("HISTESAME");

  return 0;
}
