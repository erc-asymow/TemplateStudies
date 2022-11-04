{

  TLegend* leg1 = new TLegend(0.15 ,0.75, 0.75, 0.88, "","brNDC");
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetFillColor(10);
 
  TCanvas* c = new TCanvas("c", "canvas", 600, 600);

  TString tag = "jacsvsmass100M_50shift";
  
  TFile* f_shift_data = TFile::Open("root/fit_"+tag+"_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_TEST.root", "READ");
  TGraphErrors* chi2_vs_mass = f_shift_data->Get<TGraphErrors>("chi2_vs_mass"); 

  TFile* f_shift_templ_nom = TFile::Open("root/fit_"+tag+"_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_idx0.root", "READ");
  TTree* tree2_idx0 = f_shift_templ_nom->Get<TTree>("tree2");
  double chi2_start_idx0;
  double chi2_min_idx0;
  double mass_test_idx0;
  tree2_idx0->SetBranchAddress("chi2_start", &chi2_start_idx0);
  tree2_idx0->SetBranchAddress("chi2_min", &chi2_min_idx0);
  tree2_idx0->SetBranchAddress("mass_test", &mass_test_idx0);
  tree2_idx0->GetEntry(0);

  TFile* f_shift_templ_up = TFile::Open("root/fit_"+tag+"_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_idx1.root", "READ");
  TTree* tree2_idx1 = f_shift_templ_up->Get<TTree>("tree2");
  double chi2_start_idx1;
  double chi2_min_idx1;
  double mass_test_idx1;
  tree2_idx1->SetBranchAddress("chi2_start", &chi2_start_idx1);
  tree2_idx1->SetBranchAddress("chi2_min", &chi2_min_idx1);
  tree2_idx1->SetBranchAddress("mass_test", &mass_test_idx1);
  tree2_idx1->GetEntry(0);

  TFile* f_shift_templ_down = TFile::Open("root/fit_"+tag+"_UL_10_4_A0_3_3_A1_3_3_A2_3_3_A3_3_3_A4_3_3_closure_jUL_j0_j1_j2_j3_j4_idx2.root", "READ");
  TTree* tree2_idx2 = f_shift_templ_down->Get<TTree>("tree2");
  double chi2_start_idx2;
  double chi2_min_idx2;
  double mass_test_idx2;
  tree2_idx2->SetBranchAddress("chi2_start", &chi2_start_idx2);
  tree2_idx2->SetBranchAddress("chi2_min", &chi2_min_idx2);
  tree2_idx2->SetBranchAddress("mass_test", &mass_test_idx2);
  tree2_idx2->GetEntry(0);

  double xx[3]  = {mass_test_idx0, mass_test_idx1,mass_test_idx2};
  double yy[3]  = {chi2_min_idx0,chi2_min_idx1,chi2_min_idx2};
  double exx[3] = {0.,0.,0.};
  double eyy[3] = {0.,0.,0.};	
  TGraphErrors* templ_shift = new TGraphErrors(3,xx,yy,exx,eyy);

  c->cd();
  chi2_vs_mass->SetMaximum(3);
  chi2_vs_mass->SetLineColor(kRed);
  chi2_vs_mass->SetLineWidth(2);
  chi2_vs_mass->SetMarkerStyle(3);
  chi2_vs_mass->SetTitle("");
  chi2_vs_mass->Draw();
  templ_shift->SetMarkerStyle(kFullCircle);
  templ_shift->SetMarkerColor(kBlue);
  templ_shift->SetMarkerSize(2);
  templ_shift->Draw("PSAME");

  leg1->AddEntry(chi2_vs_mass, "#chi^{2}_{min} by shifting the data", "L");
  leg1->AddEntry(templ_shift, "#chi^{2}_{min} by shifting the templates", "P");
  leg1->Draw();

  c->SaveAs("root/cmp_"+tag+".png")
}
