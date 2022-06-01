{

  TString name = "norms";
  for(auto X : ROOT::TSeqI(15) ){
    if(X<1) continue;
    TString tag(Form("%d_%d_%d_%d",X,2,2,2));
    //cout << "file " << tag << endl;
    TFile* f = TFile::Open("root/histos_"+name+"_"+tag+"_closure.root");
    if(f==nullptr || f->IsZombie()) continue;
    TTree* outtree = f->Get<TTree>("outtree");    
    int n_pdfx;
    double norms_pdfx[20];
    outtree->SetBranchAddress("n_pdfx", &n_pdfx);
    outtree->SetBranchAddress("norms_pdfx", &norms_pdfx);
    outtree->GetEntry(0);
    double norms_pdfx_update[20];
    double sum = 0.0;
    for(unsigned int i = 0; i<n_pdfx; i++) {
      norms_pdfx_update[i] = 0.5*(norms_pdfx[i]+norms_pdfx[n_pdfx-1-i]) ;
    }
    for(unsigned int i = 0; i<n_pdfx; i++) sum += norms_pdfx_update[i];
    //cout << "Sum=" << sum << endl;
    for(unsigned int i = 0; i<n_pdfx; i++) norms_pdfx_update[i] /= sum;

    cout << "std::vector<double> norms_cheb" << X << " = {";
    for(unsigned int i = 0; i<n_pdfx; i++){
      cout << norms_pdfx_update[i] ;
      if(i<n_pdfx-1) cout << ", ";
      else cout << "};" << endl;
    }
  }

}
