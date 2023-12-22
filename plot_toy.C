

{

  TCanvas* c = new TCanvas("c", "canvas", 1200, 800);
  //c->SetGridy();

  TMultiGraph* mg = new TMultiGraph();
  TLegend* leg1 = new TLegend(0.65, 0.7, 0.90, 0.90, "","brNDC");
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.035);
  leg1->SetFillColor(10);
  leg1->SetNColumns(3);
  leg1->SetHeader("");

  
  vector<int> degs = {2,4,6,8};
  for(unsigned int i = 0; i < degs.size(); i++){
    TChain* ch = new TChain("tree");
    ch->Add(Form("root/toy_mctstat*_%d_closure.root", degs[i]));

    double xx[20], yy[20], yy_jacasy[20], yy_basy[20], exx[20], eyy[20];
    
    float sigma, sigma_asy, sigma_jacasy, sigma_basy, n_effmc, n_data;
    ch->SetBranchAddress("sigma", &sigma);
    ch->SetBranchAddress("sigma_asy", &sigma_asy);
    ch->SetBranchAddress("sigma_jacasy", &sigma_jacasy);
    ch->SetBranchAddress("sigma_basy", &sigma_basy);
    ch->SetBranchAddress("n_data", &n_data);
    ch->SetBranchAddress("n_effmc", &n_effmc);
    cout << ch->GetEntries() << endl;
    for(unsigned int ie = 0 ; ie < ch->GetEntries(); ie++){
      ch->GetEntry(ie);
      xx[ie] = n_effmc/n_data;
      yy[ie] = sigma/sigma_asy;
      yy_jacasy[ie] = sigma_jacasy/sigma_asy;
      yy_basy[ie] = sigma_basy/sigma_asy;
      exx[ie] = 0.;
      eyy[ie] = 0.;
      cout << yy[ie] << endl;
      cout << yy_jacasy[ie] << endl;
      cout << yy_basy[ie] << endl;
    }
    TGraphErrors* h = new TGraphErrors(ch->GetEntries(),xx,yy,exx,eyy);
    TGraphErrors* h_jacasy = new TGraphErrors(ch->GetEntries(),xx,yy_jacasy,exx,eyy);
    TGraphErrors* h_basy = new TGraphErrors(ch->GetEntries(),xx,yy_basy,exx,eyy);
    //h->SetMinimum(0.01);
    //h_jacasy->SetMinimum(0.01);
    //h_basy->SetMinimum(0.01);
    h->SetMarkerStyle(kFullCircle);
    h->SetMarkerSize(2.0);
    h->SetLineWidth(3);
    leg1->AddEntry(h, Form("p = %d", degs[i]+1), "PL");
    h_jacasy->SetMarkerStyle(kOpenCircle);
    h_jacasy->SetLineStyle(kDashed);
    h_jacasy->SetMarkerSize(1.5);
    h_jacasy->SetLineWidth(2);
    leg1->AddEntry(h_jacasy, "U_{#infty}", "PL");
    h_basy->SetMarkerStyle(kOpenSquare);
    h_basy->SetLineStyle(kDotted);
    h_basy->SetMarkerSize(1.5);
    h_basy->SetLineWidth(2);
    leg1->AddEntry(h_basy, "b_{#infty}", "PL");
    if(i==0){
      h->SetMarkerColor(kRed);
      h->SetLineColor(kRed);
      h_jacasy->SetMarkerColor(kRed+1);
      h_jacasy->SetLineColor(kRed+1);
      h_basy->SetMarkerColor(kRed+2);
      h_basy->SetLineColor(kRed+2);
    }
    else if(i==1){
      h->SetMarkerColor(kOrange);
      h->SetLineColor(kOrange);
      h_jacasy->SetMarkerColor(kOrange+1);
      h_jacasy->SetLineColor(kOrange+1);
      h_basy->SetMarkerColor(kOrange+2);
      h_basy->SetLineColor(kOrange+2);
    }
    else if(i==2){
      h->SetMarkerColor(kGreen);
      h->SetLineColor(kGreen);
      h_jacasy->SetMarkerColor(kGreen+1);
      h_jacasy->SetLineColor(kGreen+1);
      h_basy->SetMarkerColor(kGreen+2);
      h_basy->SetLineColor(kGreen+2);
    }
    else if(i==3){
      h->SetMarkerColor(kBlue);
      h->SetLineColor(kBlue);
      h_jacasy->SetMarkerColor(kBlue+1);
      h_jacasy->SetLineColor(kBlue+1);
      h_basy->SetMarkerColor(kBlue+2);
      h_basy->SetLineColor(kBlue+2);
    }
    mg->Add(h);
    mg->Add(h_jacasy);
    mg->Add(h_basy);
  }

  c->cd();
  mg->SetMinimum(0.01);
  mg->SetMaximum(4.0);
  mg->GetXaxis()->SetTitle("N/L");
  mg->GetXaxis()->SetTitleSize(0.045);
  mg->GetXaxis()->SetTitleOffset(0.9);
  mg->GetYaxis()->SetTitleSize(0.06);
  mg->GetYaxis()->SetTitleOffset(0.7);
  mg->GetYaxis()->SetTitle("#hat{#sigma}/#hat{#sigma}_{#infty}");
  mg->GetYaxis()->SetLimits(0.01, 4);
  mg->GetXaxis()->SetLimits(0.1, 50);
  mg->GetXaxis()->SetMoreLogLabels();
  mg->GetYaxis()->CenterTitle(true);
  
  mg->Draw("apl");
  leg1->Draw();

  TF1* line = new TF1("line", "1.0", 0.1, 50 );
  line->SetLineWidth(2);
  line->SetLineStyle(9);
  line->SetLineColor(kBlack);
  line->Draw("same");
  
  mg->GetHistogram()->GetYaxis()->SetRangeUser(0.03, 5);
  //mg->GetHistogram()->GetXaxis()->SetLimits(0.1, 40);
  c->Modified();
  c->Update();
  c->SetLogx();
  c->SetLogy();

  c->SaveAs("toy_paper.eps");
  c->SaveAs("toy_paper.pdf");
  c->SaveAs("toy_paper.png");
  
}
