{


  float max_x = 0.4;
  float max_y = 2.5;
  int nbinX = 8;
  int nbinY = 6;
  int nbinsXY = nbinX*nbinY;
  TF2* tf2toy_A0 = new TF2("toy_A0", "2*x*x*(1 - 0.01*y*y)",                    0., max_x, -max_y, max_y);
  TF2* tf2toy_A1 = new TF2("toy_A1", "(0.5*x + 2*x*x)*( 0.05*y + 0.002*y*y*y)", 0., max_x, -max_y, max_y);
  TF2* tf2toy_A2 = new TF2("toy_A2", "2*x*x*(1 - 0.01*y*y)",                    0., max_x, -max_y, max_y);
  TF2* tf2toy_A3 = new TF2("toy_A3", "0.3*(x + x*x + x*x*x)*(0.1*y*y)",         0., max_x, -max_y, max_y);
  TF2* tf2toy_A4 = new TF2("toy_A4", "(1-x)*(y + y*y*y/10.)/5",                 0., max_x, -max_y, max_y);

  std::map<TString, TF2*> map;
  map.insert(std::make_pair("A0", tf2toy_A0));
  map.insert(std::make_pair("A1", tf2toy_A1));
  map.insert(std::make_pair("A2", tf2toy_A2));
  map.insert(std::make_pair("A3", tf2toy_A3));
  map.insert(std::make_pair("A4", tf2toy_A4));

  TCanvas* c = new TCanvas("c", "canvas", 1400, 600);
  c->SetLogy();
  TH1D* hdelta = new TH1D("hdelta", ";;#Delta A_{i}/A_{i}", nbinsXY*5, 0, nbinsXY*5);
  hdelta->GetXaxis()->LabelsOption("d");
  //hdelta->SetMaximum(1.0);
  hdelta->SetStats(0);
  //hdelta->SetTitle("50% prior");
  hdelta->SetTitle("50% aligned prior");
  
  TFile* f = TFile::Open("root/fit_NEWA3ZEROSMEARGW1p0NEWTOYY_4G_UL_8_6_grid_jUL_j0_j1_j2_j3_j4_DEBUG.root", "READ");
  TH2D* Vp0 = (TH2D*)f->Get("Vp0");
  for(int i=0; i<Vp0->GetXaxis()->GetNbins(); i++){
    if(i<nbinsXY) continue;
    int Ai = i/nbinsXY - 1;
    int xy = i%nbinsXY;
    int ix = xy/nbinY;
    int iy = xy%nbinY;
    float x = max_x/nbinX*(ix + 0.5);
    float y = max_y/nbinY*(iy + 0.5);
    //float delta = map[Form("A%d",Ai)]->Eval(x,y)*TMath::Sqrt(Vp0->GetBinContent(i+1,i+1));
    float delta = TMath::Sqrt(Vp0->GetBinContent(i+1,i+1));
    //float delta = map[Form("A%d",Ai)]->Eval(x,y)*0.5;
    cout << i << ": A" << Ai << " [" << ix << "," << iy << "]: x=" << x << ", y=" << y << ":   (" << map[Form("A%d",Ai)]->Eval(x,y) << " +/- "  << delta << ")" << endl;
    hdelta->SetBinContent(i+1 - nbinsXY, delta);
    hdelta->GetXaxis()->SetBinLabel(i+1 - nbinsXY, Form("A%d, [%d,%d]", Ai, ix,iy));
  }
  c->Draw();
  hdelta->Draw("HIST");
  c->SaveAs("prior_0p5alignedRel.png");
}
