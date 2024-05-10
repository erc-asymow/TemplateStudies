
{

  TString phase_space = "x0p50_y4p00";
  TString input_tag = "DEBUG";
  TString opt = "wp";

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
			   "A1", "A2", "A3", "A4"
  };
  
  TFile* fout = TFile::Open("fout_systTrad_"+opt+"_"+input_tag+".root", "RECREATE");

  TFile* fin = TFile::Open("root/fout_syst_"+opt+"_"+phase_space+"_"+input_tag+".root", "READ");

  TH2D* h_syst_i = 0;

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

    for(int isyst=0; isyst<nfp; isyst++ ){
      cout << "Syst " << isyst << endl;
      
      TH1D* h_node_i = (TH1D*)fin->Get(proc+"/h_nodes_"+proc+"_syst"+TString(Form("%d",isyst)));
      double node_x = h_node_i->GetBinContent(1);
      double node_y = TMath::Abs(h_node_i->GetBinContent(2));
      if(node_y>=yf_max) node_y -= 1e-03;
      if(node_x>=xf_max) node_x -= 1e-03;
      cout << node_x << ", " << node_y << endl;
      TH2D* h_syst_in_input = (TH2D*)fin->Get(proc+"/h_pdf_"+proc+"_syst"+TString(Form("%d",isyst))+"Up");
      
      h_syst_i = (TH2D*)h_syst_in_input->Clone(proc+"_"+TString(Form("%d", isyst)));
      //TH2D* h_systDown_i = (TH2D*)fin->Get(proc+"/h_pdf_"+proc+"_systDown"+TString(Form("%d",isyst)));

      // loop over scales
      for( unsigned int iscale=0; iscale<scales.size(); iscale++){
	cout << "Scale " << iscale << endl;
	TString scale_name = scales[iscale].first;
	int scale_idx  = scales[iscale].second;
	TFile* fin_iscale = TFile::Open("root/fout_fit_scale_"+proc+"_"+phase_space+"_"+opt+"_"+proc+"_x1_y2.root", "READ");
	TH2D* h_data_nom    = (TH2D*)fin_iscale->Get("h_data_0");
	TH2D* h_data_iscale = (TH2D*)fin_iscale->Get("h_data_"+TString(Form("%d", scale_idx)));
	int ibin_x = h_data_nom->GetXaxis()->FindBin(node_x);
	int ibin_y = h_data_nom->GetYaxis()->FindBin(node_y);
	//h_data_nom->Print("all");
	//cout << "Bin:" << ibin << endl;
	cout << h_data_iscale->GetBinContent( ibin_x, ibin_y ) << " / " << h_data_nom->GetBinContent( ibin_x, ibin_y ) << endl;	
	double binned_ratio =
	  h_data_iscale->GetBinContent( ibin_x, ibin_y ) /
	  h_data_nom->GetBinContent( ibin_x, ibin_y );
	double avgy_binned_ratio = 0.;
	for(unsigned int k = 1 ; k<=h_data_nom->GetYaxis()->GetNbins(); k++){
	  avgy_binned_ratio += h_data_iscale->GetBinContent( ibin_x, k ) / h_data_nom->GetBinContent( ibin_x, k );	  
	}
	avgy_binned_ratio /= h_data_nom->GetYaxis()->GetNbins();
	cout << "binned_ratio at the y node: " << binned_ratio << ", avg = " << avgy_binned_ratio << endl;
	binned_ratio = avgy_binned_ratio;
	
	double rescale = (binned_ratio-1.0)/shift;
	cout << "rescale by " << (binned_ratio-1.0) << "/" << shift << " = " << rescale << endl;
	
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
    }
    
  }


  fin->Close();
  fout->Close();
  return;

}
