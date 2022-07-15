
{

  TString name = "dev0";

  int fit_qt_y = 0;
  
  int axis = 1;
  
  bool normalize = true;

  std::vector<TString> tags{};

  std::vector<TFile*> files{};

  std::vector<TH1D*> histos{};
  std::vector<TH1D*> histos_mass{};

  for(auto CORRX : ROOT::TSeqI(13) ){
    if( !(CORRX==12 || CORRX==10) ) continue;
    for(auto CORRY : ROOT::TSeqI(8) ){
      if(CORRY!=4) continue;
      for(auto A0X : ROOT::TSeqI(6) ){
	if(A0X!=3) continue;
	for(auto A0Y : ROOT::TSeqI(6) ){
	  if(A0Y!=3) continue;	      
	  for(auto A1X : ROOT::TSeqI(6) ){
	    if(A1X!=3) continue;
	    for(auto A1Y : ROOT::TSeqI(6) ){
	      if(A1Y!=3) continue;	      
	      for(auto A2X : ROOT::TSeqI(6) ){
		if(A2X!=3) continue;
		for(auto A2Y : ROOT::TSeqI(6) ){
		  if(A2Y!=3) continue;	      
		  for(auto A3X : ROOT::TSeqI(6) ){
		    if( !(A3X==3 || A3X==-1) ) continue;
		    for(auto A3Y : ROOT::TSeqI(6) ){
		      if(A3Y!=3) continue;	      
		      for(auto A4X : ROOT::TSeqI(6) ){
			if( !(A4X==3 || A4X==-1) ) continue;
			for(auto A4Y : ROOT::TSeqI(6) ){
			  if(A4Y!=3) continue;	      
			  TString tag(Form("UL_%d_%d_A0_%d_%d_A1_%d_%d_A2_%d_%d_A3_%d_%d_A4_%d_%d",CORRX,CORRY,A0X,A0Y,A1X,A1Y,A2X,A2Y,A3X,A3Y,A4X,A4Y));
			  cout << "file " << tag << endl;
			  TFile* f = TFile::Open("root/histos_"+name+"_"+tag+"_closure.root");
			  if(f==nullptr || f->IsZombie()) continue;
			  files.emplace_back(f);    
			  tags.emplace_back(tag);
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  TCanvas* c = new TCanvas();
  TLegend* leg = new TLegend(0.1,0.7,0.6,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
  c->cd();

  int i = 0;
  for(auto tag : tags){
    TH2D* h_truth = files[i]->Get<TH2D>("wMC");
    TH1D* hT = axis==0?h_truth->ProjectionX("hT_py"):h_truth->ProjectionY("hT_py");
    TH2D* h_cheb = files[i]->Get<TH2D>("w");
    TH1D* hC = axis==0?h_cheb->ProjectionX("hC_py"):h_cheb->ProjectionY("hC_py");

    if(i==0){
      files[i]->cd();
      if(files[i]->FindObjectAny("wMC_up")){
	TH2D* h_up = files[i]->Get<TH2D>("wMC_up");
	TH1D* hMup =  axis==0?h_up->ProjectionX("hMup_py"):h_up->ProjectionY("hMup_py");
	hMup->Divide(hT);
	histos_mass.emplace_back((TH1D*)hMup->Clone("hMUp"));
	cout << tag << " pushed back" << endl;
      }
      if(files[i]->FindObjectAny("wMC_down")){
	TH2D* h_down = files[i]->Get<TH2D>("wMC_down");
	TH1D* hMdown =  axis==0?h_down->ProjectionX("hMdown_py"):h_down->ProjectionY("hMdown_py");
	hMdown->Divide(hT);
	histos_mass.emplace_back((TH1D*)hMdown->Clone("hMDown"));
	cout << tag << " pushed back" << endl;
      }
    }

    if(normalize){
      hC->Scale(1./hC->Integral());
      hT->Scale(1./hT->Integral());
    }
    hC->Divide(hT);
    hC->SetLineColor(i+1);
    histos.emplace_back((TH1D*)hC->Clone(tag));    
    cout << tag << " pushed back" << endl;
    i++;
  }

  c->cd();
  i = 0;
  for(auto h : histos){
    if(i==0){
      if(fit_qt_y){
	h->SetMinimum(0.99);
	h->SetMaximum(1.01);
      }
      else{
	h->SetMinimum(0.998);
	h->SetMaximum(1.002);
      }
      h->SetTitle("");
      h->SetStats(0);
      h->Draw("HIST");
    }
    else h->Draw("HISTSAME");
    leg->AddEntry(h, tags[i], "L");
    i++;
  }

  i=0;
  for(auto h : histos_mass){
    h->SetLineWidth(3);
    h->Draw("HISTSAME");
    if(i==0) leg->AddEntry(h, "#pm 10 MeV", "L");
    i++;
  }
  leg->Draw();
  gPad->SaveAs("plot_cmp_"+name+".png");
}
