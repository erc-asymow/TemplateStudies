
{

  TString phase_space = "x0p50_y4p00";
  TString input_tag = "DEBUG";
  TString opt = "z";

  bool doPlots = true;

  TCanvas* c = 0;
  if(doPlots){
    c = new TCanvas("c", "canvas", 1200, 600);
    c->Divide(1,2);
  }
  
  vector<std::pair<TString,int>> scales;
  vector<TString> scales_names;
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

  
  vector<TString> procs = {"A0",
			   "A1",
			   "A2",
			   "A3",
			   "A4"
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
    h_pdf->Write();
    h_info->Write();

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
      
      // loop over scales	
      for( unsigned int iscale=0; iscale<scales.size(); iscale++){
	cout << "\t\tScale " << iscale << endl;
	TString scale_name = scales[iscale].first;
	int scale_idx  = scales[iscale].second;
	TFile* fin_iscale = TFile::Open("root/fout_fit_scale_"+proc+"_"+phase_space+"_"+opt+"_"+proc+"_x1_y2.root", "READ");
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

	double avgy_binned_ratio = avg_scale/avg_nom;
	//avgy_binned_ratio /= h_data_nom->GetYaxis()->GetNbins();
	cout << "\t\tbinned_ratio at the y node: " << binned_ratio << ", avg = " << avgy_binned_ratio << endl;

	// take the average
	binned_ratio = avgy_binned_ratio;
	double rescale = (binned_ratio-1.0)/shift;
	cout << "\t\tRescale by (" << binned_ratio << " - 1.0)" << "/" << shift << " = " << rescale << endl;

	// consider only those scale variations that don't change the sign
	if(TMath::Abs(binned_ratio - 1)>envelope && binned_ratio > 0.){
 	  for(unsigned int ix=1; ix<=h_systenvUp_i->GetNbinsX(); ix++){
	    if(h_systenvUp_i->GetXaxis()->GetBinUpEdge(ix)>xf_max) continue;
	    for(unsigned int iy=1; iy<=h_systenvUp_i->GetNbinsY(); iy++){
	      if(h_systenvUp_i->GetYaxis()->GetBinUpEdge(iy)>yf_max) continue;
	      double nom_val = h_pdf->GetBinContent(ix,iy);
	      double syst_val = h_syst_i->GetBinContent(ix,iy);
	      double corr = nom_val!=0. ? syst_val/nom_val - 1.0 : 0.;
	      corr *= rescale;
	      //double clip_corr =  std::max( -1.0, std::min(corr,1.0) );
	      //if(clip_corr!=corr)
	      //cout << "\t\tClipping corr: " << corr << " --> " << clip_corr << endl;
	      //corr = clip_corr;
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
	h_new_syst_i->Write();
	
	fin_iscale->Close();
      }

      fout->cd(proc);
      int node_bin = h_pdf->FindBin(node_x, node_y);
      bool flip = false;
      if( h_systenvUp_i->GetBinContent( node_bin)/h_pdf->GetBinContent( node_bin) > h_systenvDown_i->GetBinContent( node_bin)/h_pdf->GetBinContent( node_bin) ){
	h_systenvUp_i->Write("h_pdf_"+proc+"_syst"+TString(Form("%d",isyst))+"Up");
	h_systenvDown_i->Write("h_pdf_"+proc+"_syst"+TString(Form("%d",isyst))+"Down");
      }
      else{
	h_systenvUp_i->Write("h_pdf_"+proc+"_syst"+TString(Form("%d",isyst))+"Down");
	h_systenvDown_i->Write("h_pdf_"+proc+"_syst"+TString(Form("%d",isyst))+"Up");
	flip = true;
      }

      if(doPlots){
	c->cd( flip ? 2 : 1);
	h_systenvUp_i->SetTitle(opt+", "+phase_space+", "+proc+", syst"+TString(Form("%d at node (x=%.2f, y=%.2f)",isyst,node_x, node_y)));
	h_systenvUp_i->Divide(h_pdf);
	h_systenvUp_i->SetStats(0);
	h_systenvUp_i->Draw("COLZ");
	c->cd(flip ? 1 : 2);
	h_systenvDown_i->SetTitle(opt+", "+phase_space+", "+proc+", syst"+TString(Form("%d at node (x=%.2f, y=%.2f)",isyst,node_x, node_y)));
	h_systenvDown_i->Divide(h_pdf);
	h_systenvDown_i->SetStats(0);
	h_systenvDown_i->Draw("COLZ");
	c->SaveAs( "plots/ratio_scale_"+opt+"_"+phase_space+"_"+proc+"_syst"+TString(Form("%d",isyst))+".png" );
      }
      
      cout << "\tMax scale triggered by " << scales_names[env_scale] << endl;
    }
    
  }


  fin->Close();
  fout->Close();

  return;

}
