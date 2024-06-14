
{

  TString phase_space = "x0p50_y4p00";
  TString input_tag = "DEBUG";
  TString opt = "z";

  bool doPlots = true;

  bool extendX = true;

  bool keepYdependence = true;
  
  TCanvas* c = 0;
  if(doPlots){
    c = new TCanvas("c", "canvas", 1200, 600);
    c->Divide(1,3);
  }
  
  vector<std::pair<TString,int>> scales;
  vector<TString> scales_names;
  /*
  scales.emplace_back( std::make_pair<TString,int>("twotwo", 1) );
  scales_names.emplace_back("muRmuFUp");
  scales.emplace_back( std::make_pair<TString,int>("point5point5", 2) );
  scales_names.emplace_back("muRmuFDown");
  scales.emplace_back( std::make_pair<TString,int>("onepoint5", 3) );
  scales_names.emplace_back("muFDown");
  scales.emplace_back( std::make_pair<TString,int>("onetwo", 6) );
  scales_names.emplace_back("muFUp");
  scales.emplace_back( std::make_pair<TString,int>("point5one", 4) );
  scales_names.emplace_back("muRDown");
  scales.emplace_back( std::make_pair<TString,int>("twoone", 5) );
  scales_names.emplace_back("muRUp");
  */
  scales.emplace_back( std::make_pair<TString,int>("twotwo_oneone", 7) );
  scales_names.emplace_back("muRmuFUp");
  scales.emplace_back( std::make_pair<TString,int>("point5point5_oneone", 11) );
  scales_names.emplace_back("muRmuFDown");
  scales.emplace_back( std::make_pair<TString,int>("onepoint5_oneone", 15) );
  scales_names.emplace_back("muFDown");
  scales.emplace_back( std::make_pair<TString,int>("onetwo_oneone", 27) );
  scales_names.emplace_back("muFUp");
  scales.emplace_back( std::make_pair<TString,int>("point5one_oneone", 19) );
  scales_names.emplace_back("muRDown");
  scales.emplace_back( std::make_pair<TString,int>("twoone_oneone", 23) );
  scales_names.emplace_back("muRUp");
  
  vector<TString> procs = {"A0",
			   "A1",
			   "A2",
			   "A3",
			   "A4",
			   "A5",
			   "A6",
			   "A7"
  };

  TFile* fout = TFile::Open("fout_systTrad_"+opt+"_"+phase_space+"_"+input_tag+".root", "RECREATE");
  cout << fout->GetName() <<  " created" << endl;
  
  TFile* fin = TFile::Open("root/fout_syst_"+opt+"_"+phase_space+"_"+input_tag+".root", "READ");
  cout << fin->GetName() << " opened" << endl;
    
  TH2D* h_syst_i = 0;
  TH2D* h_systenvUp_i = 0;
  TH2D* h_systenvDown_i = 0;

  // loop over procs
  for(unsigned int iproc = 0; iproc<procs.size(); iproc++){
    
    TString proc = procs[iproc];
    cout << proc << endl;
    
    fout->mkdir(proc);
    
    TH1D* h_info = (TH1D*)fin->Get(proc+"/h_info_"+proc);
    int nfpx = int(h_info->GetBinContent(8));
    int nfpy = int(h_info->GetBinContent(9)); 
    double shift = h_info->GetBinContent(10);
    double xf_max = h_info->GetBinContent(11);
    double yf_max = h_info->GetBinContent(12);    
    int nfp = nfpx*nfpy;
    
    TH2D* h_pdf  = (TH2D*)fin->Get(proc+"/h_pdf_"+proc);
    fout->cd(proc);
    h_info->Write("h_info_"+proc);
    if(!extendX) h_pdf->Write();
    else{
      TFile* f_ext = TFile::Open("/scratch/tanmay/Output2016/Final_Uses/V5_extraBinadded/root_files_ang_coeff_qtbyQ_2d.root", "READ");
      TH2D* h_data_ext = (TH2D*)f_ext->Get("ang_coeff_"+opt+"_2d_qtbyQ_vs_absy_A_"+TString(Form("%d", iproc)));
      const Double_t* binsX = h_pdf->GetXaxis()->GetXbins()->GetArray();
      int nbinsX_ext = h_pdf->GetXaxis()->GetNbins()+1;
      vector<double> vbinsX_ext;
      vbinsX_ext.reserve(nbinsX_ext+1);
      for(unsigned int ib = 0; ib<=h_pdf->GetXaxis()->GetNbins(); ib++ ) vbinsX_ext.push_back( binsX[ib] );
      vbinsX_ext.push_back( 5.0 );
      TH2D* h_pdf_ext = new TH2D("h_pdf_ext_"+proc, "", nbinsX_ext, vbinsX_ext.data(), h_pdf->GetYaxis()->GetNbins(), h_pdf->GetYaxis()->GetXbins()->GetArray());
      for(unsigned int ibx=1; ibx<=h_pdf_ext->GetXaxis()->GetNbins();ibx++){
	for(unsigned int iby=1; iby<=h_pdf_ext->GetXaxis()->GetNbins();iby++){
	  if(ibx<h_pdf_ext->GetXaxis()->GetNbins())
	    h_pdf_ext->SetBinContent(ibx,iby, h_pdf->GetBinContent(ibx,iby));
	  else{ 
	    int bin_data = h_data_ext->FindBin( h_pdf_ext->GetXaxis()->GetBinCenter(ibx), h_pdf_ext->GetYaxis()->GetBinCenter(iby)  );
	    h_pdf_ext->SetBinContent(ibx,iby, h_data_ext->GetBinContent(bin_data) );
	  }
	}
      }
      fout->cd(proc);
      h_pdf_ext->Write("h_pdf_"+proc);
      f_ext->Close();
    }

    for(int isyst=0; isyst<nfp; isyst++ ){
      cout << "\tSyst " << isyst << endl;
      
      TH1D* h_node_i = (TH1D*)fin->Get(proc+"/h_nodes_"+proc+"_syst"+TString(Form("%d",isyst)));
      double node_x = h_node_i->GetBinContent(1);
      double node_y = TMath::Abs(h_node_i->GetBinContent(2));
      if(node_y>=yf_max) node_y -= 1e-03;
      if(node_x>=xf_max) node_x -= 1e-03;
      cout << "\tNodes: (" <<  node_x << ", " << node_y << ")" << endl;
      TH2D* h_syst_in_input = (TH2D*)fin->Get(proc+"/h_pdf_"+proc+"_syst"+TString(Form("%d",isyst))+"Up");    
      
      h_syst_i        = (TH2D*)h_syst_in_input->Clone(proc+"_"+TString(Form("%d", isyst)));
      h_systenvUp_i   = (TH2D*)h_syst_in_input->Clone("h_pdf_"+proc+"_syst"+TString(Form("%d",isyst))+"_scaleUp");
      h_systenvDown_i = (TH2D*)h_syst_in_input->Clone("h_pdf_"+proc+"_syst"+TString(Form("%d",isyst))+"_scaleDown");
      double envelope = 0.;
      int env_scale = -1;

      vector<float> envelopes(h_syst_i->GetNbinsY());
      
      // loop over scales	
      for( unsigned int iscale=0; iscale<scales.size(); iscale++){
	cout << "\t\tScale " << iscale << endl;
	TString scale_name = scales[iscale].first;
	int scale_idx  = scales[iscale].second;
	//TFile* fin_iscale = TFile::Open("root/fout_fit_scale_"+proc+"_"+phase_space+"_"+opt+"_"+proc+"_x1_y2.root", "READ");
	TFile* fin_iscale = TFile::Open("root/fout_fit_scale_31point_"+proc+"_"+phase_space+"_"+opt+"_"+proc+"_x1_y2.root", "READ");
	TH2D* h_data_nom    = (TH2D*)fin_iscale->Get("h_data_0");
	TH2D* h_data_iscale = (TH2D*)fin_iscale->Get("h_data_"+TString(Form("%d", scale_idx)));
	int ibin_x = h_data_nom->GetXaxis()->FindBin(node_x);	
	int ibin_y = h_data_nom->GetYaxis()->FindBin(node_y);

	if(isyst==0){
	  cout << "\t\tReading systematics from following x-node to prevent large errors" << endl;
	  ibin_x = h_data_nom->GetXaxis()->FindBin(((TH1D*)fin->Get(proc+"/h_nodes_"+proc+"_syst"+TString(Form("%d",isyst+1))))->GetBinContent(1));
	}
	
	//h_data_nom->Print("all");
	//cout << "Bin:" << ibin << endl;
	//cout << "\t\tInputs:" << h_data_iscale->GetBinContent( ibin_x, ibin_y ) << " / " << h_data_nom->GetBinContent( ibin_x, ibin_y ) << endl;	
	double binned_ratio =
	  h_data_iscale->GetBinContent( ibin_x, ibin_y ) /
	  h_data_nom->GetBinContent( ibin_x, ibin_y );

	      
	double avg_nom = 0.;
	double avg_scale = 0.;
	for(unsigned int k = 1 ; k<=h_data_nom->GetYaxis()->GetNbins(); k++){
	  //avgy_binned_ratio += h_data_iscale->GetBinContent( ibin_x, k ) / h_data_nom->GetBinContent( ibin_x, k );
	  avg_nom   += h_data_nom->GetBinContent( ibin_x, k );
	  avg_scale += h_data_iscale->GetBinContent( ibin_x, k );
	}

	vector<float> rescales(h_systenvUp_i->GetNbinsY());
	for(unsigned int iy=1; iy<=h_systenvUp_i->GetNbinsY(); iy++){
	  double y_i = h_systenvUp_i->GetYaxis()->GetBinCenter(iy);
	  int ibin = h_data_nom->FindBin(node_x, y_i);
	  float ratio_iy = TMath::Abs(h_data_nom->GetBinContent( ibin ))>0. ? h_data_iscale->GetBinContent( ibin ) / h_data_nom->GetBinContent( ibin ) : 1.0;
	  if(TMath::Abs(ratio_iy-1.0)>=TMath::Abs(envelopes.at(iy-1))){
	    rescales[iy-1] = (ratio_iy - 1.0)/shift;
	    envelopes[iy-1] = ratio_iy-1.0;
	    //if(isyst==5 && iproc==0) cout << y_i << ", ibin " << ibin << ": ratio=" << ratio_iy << ", rescale=" <<  rescales[iy-1] << ", envelope=" << envelopes.at(iy-1) << endl;
	  }
	  else{
	    rescales[iy-1] = envelopes[iy-1]/shift;
	  }
	}

	double avgy_binned_ratio = avg_scale/avg_nom;
	//avgy_binned_ratio /= h_data_nom->GetYaxis()->GetNbins();
	cout << "\t\tbinned_ratio at the y node: " << binned_ratio << ", avg = " << avgy_binned_ratio << endl;

	// take the average
	binned_ratio = avgy_binned_ratio;
	double rescale = (binned_ratio-1.0)/shift;
	cout << "\t\tRescale by (" << binned_ratio << " - 1.0)" << "/" << shift << " = " << rescale << endl;

	// consider only those scale variations that don't change the sign
	if( (TMath::Abs(binned_ratio - 1)>envelope && binned_ratio > 0.) || keepYdependence){
 	  for(unsigned int ix=1; ix<=h_systenvUp_i->GetNbinsX(); ix++){
	    if(h_systenvUp_i->GetXaxis()->GetBinUpEdge(ix)>xf_max) continue;
	    for(unsigned int iy=1; iy<=h_systenvUp_i->GetNbinsY(); iy++){
	      if(h_systenvUp_i->GetYaxis()->GetBinUpEdge(iy)>yf_max) continue;
	      double nom_val = h_pdf->GetBinContent(ix,iy);
	      double syst_val = h_syst_i->GetBinContent(ix,iy);
	      double corr = nom_val!=0. ? syst_val/nom_val - 1.0 : 0.;
	      
	      if(keepYdependence){
		corr *=  TMath::Abs(rescales.at(iy-1));		  
		h_systenvUp_i->SetBinContent(ix,iy, nom_val*(1 + corr));
		h_systenvDown_i->SetBinContent(ix,iy, nom_val*(1 - corr));
		continue;	      
	      }
	
	      corr *= rescale;

	      //cout << "new envelope at " << iy-1 << ": rescales=" << rescale << " envelope=" << envelope << " => corr=" << corr  << endl;
	      //double clip_corr =  std::max( -1.0, std::min(corr,1.0) );
	      //if(clip_corr!=corr)
	      //cout << "\t\tClipping corr: " << corr << " --> " << clip_corr << endl;
	      //corr = clip_corr;
	      //cout << "(" << ix << "," << iy << ") => " << 1+corr << endl;
	      h_systenvUp_i->SetBinContent(ix,iy, nom_val*(1 + corr));
	      h_systenvDown_i->SetBinContent(ix,iy, nom_val*(1 - corr));
	    }
	  }	  
	  envelope = TMath::Abs(binned_ratio - 1);
	  cout << "\t\tsyst " << isyst << ": new envelope is |" << binned_ratio  << " - 1.0| = " << envelope << endl;  
	  env_scale = iscale;
	}
	else{
	  cout << "\t\tsyst " << isyst << ": skip because " << TMath::Abs(binned_ratio - 1) << " <= " << envelope << endl;
	}

	//fout->cd(proc);
	TH2D* h_new_syst_i = (TH2D*)h_syst_i->Clone("h_pdf_"+proc+"_syst"+TString(Form("%d",isyst))+"_"+scales_names[iscale]);
	for(unsigned int ix=1; ix<=h_new_syst_i->GetNbinsX(); ix++){
 	  if(h_new_syst_i->GetXaxis()->GetBinUpEdge(ix)>xf_max) continue;
	  for(unsigned int iy=1; iy<=h_new_syst_i->GetNbinsY(); iy++){
 	    if(h_new_syst_i->GetYaxis()->GetBinUpEdge(iy)>yf_max) continue;
	    double nom_val = h_pdf->GetBinContent(ix,iy);
	    double old_val = h_syst_i->GetBinContent(ix,iy);
	    double corr = old_val>0. ? nom_val/old_val - 1.0 : 0.;
	    corr *= rescale;
	    h_new_syst_i->SetBinContent(ix,iy, nom_val*(1+corr));
	  }
	}
	fout->cd(proc);

	if(!extendX)
	  h_new_syst_i->Write();

	/*
	else{
	  const Double_t* binsX = h_new_syst_i->GetXaxis()->GetXbins()->GetArray();
	  int nbinsX_ext = h_new_syst_i->GetXaxis()->GetNbins()+1;
	  vector<double> vbinsX_ext;
	  vbinsX_ext.reserve(nbinsX_ext+1);
	  for(unsigned int ib = 0; ib<=h_new_syst_i->GetXaxis()->GetNbins(); ib++ ){
	    vbinsX_ext.push_back( binsX[ib] );
	  }
	  vbinsX_ext.push_back( 5.0 );
	  TH2D* h_new_syst_i_ext = new TH2D( TString(h_new_syst_i->GetName())+"_ext", "", nbinsX_ext,
					     vbinsX_ext.data(),
					     h_new_syst_i->GetYaxis()->GetNbins(), h_new_syst_i->GetYaxis()->GetXbins()->GetArray());
	  for(unsigned int ibx=1; ibx<=h_new_syst_i_ext->GetXaxis()->GetNbins();ibx++){
	    for(unsigned int iby=1; iby<=h_new_syst_i_ext->GetXaxis()->GetNbins();iby++){
	      if(ibx<h_new_syst_i_ext->GetXaxis()->GetNbins())
		h_new_syst_i_ext->SetBinContent(ibx,iby, h_new_syst_i->GetBinContent(ibx,iby));
	      else 
		h_new_syst_i_ext->SetBinContent(ibx,iby, h_new_syst_i->GetBinContent(ibx-1,iby));
	    }
	  }
	  h_new_syst_i_ext->Write(proc+"/"+TString(h_new_syst_i->GetName()));
	}
	*/
	
	fin_iscale->Close();
      }

      fout->cd(proc);
      int node_bin = h_pdf->FindBin(node_x, node_y);
      bool flip = false;

      TString name_up;
      TString name_down;
      if( h_systenvUp_i->GetBinContent( node_bin)/h_pdf->GetBinContent( node_bin) > h_systenvDown_i->GetBinContent( node_bin)/h_pdf->GetBinContent( node_bin) ){
	name_up   =  "h_pdf_"+proc+"_syst"+TString(Form("%d",isyst))+"Up";
	name_down =  "h_pdf_"+proc+"_syst"+TString(Form("%d",isyst))+"Down";
      }
      else{
	name_up   =  "h_pdf_"+proc+"_syst"+TString(Form("%d",isyst))+"Down";
	name_down =  "h_pdf_"+proc+"_syst"+TString(Form("%d",isyst))+"Up";
	flip = true;
      }

      if(!extendX){
	h_systenvUp_i->Write( name_up );
	h_systenvDown_i->Write( name_down );
      }
      else{

	TH2D* h_pdf_ext = (TH2D*)fout->Get(proc+"/h_pdf_"+proc);
	
	const Double_t* binsXUp = h_systenvUp_i->GetXaxis()->GetXbins()->GetArray();
	int nbinsXUp_ext = h_systenvUp_i->GetXaxis()->GetNbins()+1;
	vector<double> vbinsXUp_ext;
	vbinsXUp_ext.reserve(nbinsXUp_ext+1);
	for(unsigned int ib = 0; ib<=h_systenvUp_i->GetXaxis()->GetNbins(); ib++ ){
	  vbinsXUp_ext.push_back( binsXUp[ib] );
	  //cout << vbinsXUp_ext[ib] << endl; 
	}
	vbinsXUp_ext.push_back( 5.0 );
	//cout << vbinsXUp_ext[nbinsXUp_ext-1] << endl;
	TH2D* h_systenvUp_i_ext = new TH2D( TString(h_systenvUp_i->GetName())+"_ext", "", nbinsXUp_ext,
					    vbinsXUp_ext.data(),
					    h_systenvUp_i->GetYaxis()->GetNbins(), h_systenvUp_i->GetYaxis()->GetXbins()->GetArray());
	for(unsigned int ibx=1; ibx<=h_systenvUp_i_ext->GetXaxis()->GetNbins();ibx++){
	  for(unsigned int iby=1; iby<=h_systenvUp_i_ext->GetXaxis()->GetNbins();iby++){
	    if(ibx<h_systenvUp_i_ext->GetXaxis()->GetNbins())
	      h_systenvUp_i_ext->SetBinContent(ibx,iby, h_systenvUp_i->GetBinContent(ibx,iby));
	    else{ 
	      //TH2D* h_pdf_ext = (TH2D*)fout->Get(proc+"/h_pdf_"+proc);
	      h_systenvUp_i_ext->SetBinContent(ibx,iby, h_systenvUp_i->GetBinContent(ibx-1,iby)/h_pdf_ext->GetBinContent(ibx-1,iby)*h_pdf_ext->GetBinContent(ibx,iby) );
	    }
	  }
	}
	h_systenvUp_i_ext->Write(name_up);

	const Double_t* binsXDown = h_systenvDown_i->GetXaxis()->GetXbins()->GetArray();
	int nbinsXDown_ext = h_systenvDown_i->GetXaxis()->GetNbins()+1;
	vector<double> vbinsXDown_ext;
	vbinsXDown_ext.reserve(nbinsXDown_ext+1);
	for(unsigned int ib = 0; ib<=h_systenvDown_i->GetXaxis()->GetNbins(); ib++ ){
	  vbinsXDown_ext.push_back( binsXDown[ib] );
	  //cout << vbinsXDown_ext[ib] << endl; 
	}
	vbinsXDown_ext.push_back( 5.0 );
	//cout << vbinsXDown_ext[nbinsXDown_ext-1] << endl;
	TH2D* h_systenvDown_i_ext = new TH2D( TString(h_systenvDown_i->GetName())+"_ext", "", nbinsXDown_ext,
					    vbinsXDown_ext.data(),
					    h_systenvDown_i->GetYaxis()->GetNbins(), h_systenvDown_i->GetYaxis()->GetXbins()->GetArray());
	for(unsigned int ibx=1; ibx<=h_systenvDown_i_ext->GetXaxis()->GetNbins();ibx++){
	  for(unsigned int iby=1; iby<=h_systenvDown_i_ext->GetXaxis()->GetNbins();iby++){
	    if(ibx<h_systenvDown_i_ext->GetXaxis()->GetNbins())
	      h_systenvDown_i_ext->SetBinContent(ibx,iby, h_systenvDown_i->GetBinContent(ibx,iby));
	    else {
	      //TH2D* h_pdf_ext = (TH2D*)fout->Get(proc+"/h_pdf_"+proc);
	      h_systenvDown_i_ext->SetBinContent(ibx,iby, h_systenvDown_i->GetBinContent(ibx-1,iby)/h_pdf_ext->GetBinContent(ibx-1,iby)*h_pdf_ext->GetBinContent(ibx,iby) );
	    }
	  }
	}
	h_systenvDown_i_ext->Write(name_down);

      }

      if(doPlots){
	c->cd(1);
 	TH2D* h_pdf_ext = (TH2D*)fout->Get(proc+"/h_pdf_"+proc);
	h_pdf_ext->SetStats(0);
	h_pdf_ext->Draw("COLZ");
	c->cd(2);	
	TH2D* h_up = (TH2D*)fout->Get(proc+"/h_pdf_"+proc+"_syst"+TString(Form("%d",isyst))+"Up");
	h_up->SetTitle(opt+", "+phase_space+", "+proc+", syst"+TString(Form("%d at node (x=%.2f, y=%.2f)",isyst,node_x, node_y)));
	h_up->Divide(h_pdf_ext);
	h_up->SetStats(0);
	h_up->Draw("COLZ");
	c->cd(3);
	TH2D* h_down = (TH2D*)fout->Get(proc+"/h_pdf_"+proc+"_syst"+TString(Form("%d",isyst))+"Down");
	h_down->SetTitle(opt+", "+phase_space+", "+proc+", syst"+TString(Form("%d at node (x=%.2f, y=%.2f)",isyst,node_x, node_y)));
	h_down->Divide(h_pdf_ext);
	h_down->SetStats(0);
	h_down->Draw("COLZ");
	c->SaveAs( "plots/ratio_scale_"+opt+"_"+phase_space+"_"+proc+"_syst"+TString(Form("%d",isyst))+".png" );
      }
      
      cout << "\tMax scale triggered by " << scales_names[env_scale] << endl;
    }
    
  }


  fin->Close();
  fout->Close();

  return;

}
