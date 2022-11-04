{

  TFile* fUL = TFile::Open("./root/histos_addmass7_seed0_UL_10_4_grid.root", "READ");
  vector<TString> histos_from_UL = {"h", "hMC", "hMC_up", "hMC_down"};
  for(int j=0; j<40; j++) histos_from_UL.emplace_back( TString(Form("jac_%d",j)) );
  //for(int m=0; m<20; m++) histos_from_UL.emplace_back( TString(Form("hMC_mass%d",m)) );
  
  TFile* fA0 = TFile::Open("./root/histos_addmass7_seed1_UL_10_4_grid.root", "READ");
  vector<TString> histos_from_A0 = {};
  for(int j=40; j<80; j++) histos_from_A0.emplace_back( TString(Form("jac_%d",j)) );
  //for(int m=0; m<20; m++) histos_from_A0.emplace_back( TString(Form("hMC_mass%d",m)) );
 
  TFile* fA1 = TFile::Open("./root/histos_addmass7_seed2_UL_10_4_grid.root", "READ");
  vector<TString> histos_from_A1 = {};
  for(int j=80; j<120; j++) histos_from_A1.emplace_back( TString(Form("jac_%d",j)) );
  //for(int m=0; m<20; m++) histos_from_A1.emplace_back( TString(Form("hMC_mass%d",m)) );
  
  TFile* fA2 = TFile::Open("./root/histos_addmass7_seed3_UL_10_4_grid.root", "READ");
  vector<TString> histos_from_A2 = {};
  for(int j=120; j<160; j++) histos_from_A2.emplace_back( TString(Form("jac_%d",j)) );
  //for(int m=0; m<20; m++) histos_from_A2.emplace_back( TString(Form("hMC_mass%d",m)) );
  
  TFile* fA3 = TFile::Open("./root/histos_addmass7_seed4_UL_10_4_grid.root", "READ");
  vector<TString> histos_from_A3 = {};
  for(int j=160; j<200; j++) histos_from_A3.emplace_back( TString(Form("jac_%d",j)) );
  //for(int m=0; m<20; m++) histos_from_A3.emplace_back( TString(Form("hMC_mass%d",m)) );

  TFile* fA4 = TFile::Open("./root/histos_addmass7_seed5_UL_10_4_grid.root", "READ");
  vector<TString> histos_from_A4 = {};
  for(int j=200; j<240; j++) histos_from_A4.emplace_back( TString(Form("jac_%d",j)) );
  for(int m=0; m<20; m++) histos_from_A4.emplace_back( TString(Form("hMC_mass%d",m)) );
 
  TFile* fout = TFile::Open("./root/histos_addmass7_merged5_UL_10_4_grid.root", "RECREATE");
  fout->cd();
  for(auto h : histos_from_UL) (fUL->Get<TH2D>(h))->Write("", TObject::kOverwrite);
  for(auto h : histos_from_A0) (fA0->Get<TH2D>(h))->Write("", TObject::kOverwrite);
  for(auto h : histos_from_A1) (fA1->Get<TH2D>(h))->Write("", TObject::kOverwrite);
  for(auto h : histos_from_A2) (fA2->Get<TH2D>(h))->Write("", TObject::kOverwrite);
  for(auto h : histos_from_A3) (fA3->Get<TH2D>(h))->Write("", TObject::kOverwrite);
  for(auto h : histos_from_A4) (fA4->Get<TH2D>(h))->Write("", TObject::kOverwrite);
  TTree* outtree = fUL->Get<TTree>("outtree");
  fout->cd();
  outtree->CloneTree()->Write("", TObject::kOverwrite);
  fout->Close();
}
