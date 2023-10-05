{

  TString systname = "scale_A0";

  TString proc = "A0";
  
  vector<TString> phasespace = {"x0p30_y3p00",
				"x0p40_y3p50"
  };
  
  vector<TString> runs = {"wp",
			  "wm", "z"
  };
  vector<TString> dx = {"1", "2", "3", "4"
			//"6",
			//"7","8","9","10","11","12"
  };
  vector<TString> dy = {"2", "4", "6"
			//"2",
			//"4"
  };

  for(int ip = 0; ip<phasespace.size(); ip++ ){
    for(int ip2 = 0; ip2<runs.size(); ip2++ ){
          for(int idx = 0; idx<dx.size(); idx++ ){
	    for(int idy = 0; idy<dy.size(); idy++ ){
	      TFile* fin = TFile::Open("fout_fit_"+systname+"_"+phasespace[ip]+"_"+runs[ip2]+"_"+proc+"_x"+dx[idx]+"_y"+dy[idy]+".root", "READ");
	      TH1D* hinfo = (TH1D*)fin->Get("h_info");
	      TH1D* hstart = (TH1D*)fin->Get("h_start");
	      int ns = hinfo->GetBinContent(3);
	      TCanvas* c = new TCanvas("c", "canvas", 1200, 600);
	      c->Divide(2,1);
	      for(int i = 0; i < ns; i++){
		TH2D* e_i = (TH2D*)fin->Get(Form("h_exp_%d",i));
		TH2D* d_i = (TH2D*)fin->Get(Form("h_data_%d",i));
		TH2D* d_i_clone = (TH2D*) d_i->Clone(phasespace[ip]+"_"+runs[ip2]+"_"+proc+"_x"+dx[idx]+"_y"+dy[idy]+"_"+TString(Form("%d",i)));
		d_i_clone->Divide(hstart);
		d_i_clone->SetMinimum(0.8);
		d_i_clone->SetMaximum(1.2);
		d_i_clone->SetStats(0);
		c->cd(1);		
		d_i_clone->Draw("colz");		

		e_i->Divide(d_i);
		e_i->SetStats(0);
		e_i->SetMinimum(0.97);
		e_i->SetMaximum(1.03);

		c->cd(2);
		gPad->SetRightMargin(0.15);
		e_i->Draw("COLZ");
		c->SaveAs(systname+"_"+TString(Form("_var%d", i))+"_"+phasespace[ip]+"_"+runs[ip2]+"_"+proc+"_x"+dx[idx]+"_y"+dy[idy]+".png");
		return;
	      }
	      fin->Close();
	    }
	  }
    }
  }
  
  return; 
}
