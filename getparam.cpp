#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TH2D.h"
#include "TF1.h"
#include "TF2.h"
#include "TFitResultPtr.h"
#include <TStopwatch.h>
#include <TGraphErrors.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>
#include <boost/program_options.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

using namespace boost::program_options;


double cheb_fct(double *var, double *par){
  double x=var[0];
  double den = 0.;
  double num = 0.;
  for(int i = 0; i <= par[0] ; i++){ // par[0]=n
    int sign = i%2==0 ? +1 :-1;
    double xj = TMath::Cos((par[0]-i)*TMath::Pi()/par[0]);
    xj *= par[2];
    xj += par[1];
    if(x==xj) return 1.0;// protect from nan      
    double val = sign/(x-xj);
    if(i==0 || i==par[0]) val *= 0.5;
    den += val;
    if(i==par[3]) num = val; // par[3]=m
  }                                             
  return num/den;
}

int main(int argc, char* argv[])
{
  
  TStopwatch sw;
  sw.Start();

  variables_map vm;
  try
    {
      options_description desc{"Options"};
      desc.add_options()
	("help,h",      "Help screen")
	("nevents", value<long>()->default_value(1000), "number of events")
	("ntoys",   value<int>()->default_value(-1), "number of toys")
	("extrabinsX",   value<int>()->default_value(10), "pdf x granularity")
	("extrabinsY",   value<int>()->default_value(10), "pdf y granularity")
	("dULx",   value<int>()->default_value(2), "max degree in x of corrxy")
	("dULy",   value<int>()->default_value(2), "max degree in y of corrxy")
	("dA0x",   value<int>()->default_value(2), "max degree in x for A0")
	("dA0y",   value<int>()->default_value(2), "max degree in y for A0")
	("dA1x",   value<int>()->default_value(2), "max degree in x for A1")
	("dA1y",   value<int>()->default_value(2), "max degree in y for A1")
	("dA2x",   value<int>()->default_value(2), "max degree in x for A2")
	("dA2y",   value<int>()->default_value(2), "max degree in y for A2")
	("dA3x",   value<int>()->default_value(2), "max degree in x for A3")
	("dA3y",   value<int>()->default_value(2), "max degree in y for A3")
	("dA4x",   value<int>()->default_value(2), "max degree in x for A4")
	("dA4y",   value<int>()->default_value(2), "max degree in y for A4")
	("cULx",   value<int>()->default_value(1), "constraint in x of corrxy")
	("cULy",   value<int>()->default_value(0), "constraint in y of corrxy")
	("cA0x",   value<int>()->default_value(1), "constraint to 0 in x=0 for A0")
	("cA0y",   value<int>()->default_value(0), "constraint to 0 in y=0 for A0")
	("cA1x",   value<int>()->default_value(1), "constraint to 0 in x=0 for A1")
	("cA1y",   value<int>()->default_value(0), "constraint to 0 in y=0 for A1")
	("cA2x",   value<int>()->default_value(1), "constraint to 0 in x=0 for A2")
	("cA2y",   value<int>()->default_value(0), "constraint to 0 in y=0 for A2")
	("cA3x",   value<int>()->default_value(1), "constraint to 0 in x=0 for A3")
	("cA3y",   value<int>()->default_value(1), "constraint to 0 in y=0 for A3")
	("cA4x",   value<int>()->default_value(0), "constraint to 0 in x=0 for A4")
	("cA4y",   value<int>()->default_value(0), "constraint to 0 in y=0 for A4")
	("fULx",   value<int>()->default_value(-1), "max degree of modifier in x of corrxy")
	("fULy",   value<int>()->default_value(-1), "max degree of modifier in y of corrxy")
	("fA0x",   value<int>()->default_value(-1), "max degree of modifier in x for A0")
	("fA0y",   value<int>()->default_value(-1), "max degree of modifier in y for A0")
	("fA1x",   value<int>()->default_value(-1), "max degree of modifier in x for A1")
	("fA1y",   value<int>()->default_value(-1), "max degree of modifier in y for A1")
	("fA2x",   value<int>()->default_value(-1), "max degree of modifier in x for A2")
	("fA2y",   value<int>()->default_value(-1), "max degree of modifier in y for A2")
	("fA3x",   value<int>()->default_value(-1), "max degree of modifier in x for A3")
	("fA3y",   value<int>()->default_value(-1), "max degree of modifier in y for A3")
	("fA4x",   value<int>()->default_value(-1), "max degree of modifier in x for A4")
	("fA4y",   value<int>()->default_value(-1), "max degree of modifier in y for A4")
	("x_max",   value<double>()->default_value(-1.0), "max x value for fit")
	("y_max",  value<double>()->default_value(-1.0), "max y value for fit")
	("xf_max",  value<double>()->default_value(-1.0), "max x value for syst")
	("yf_max",  value<double>()->default_value(-1.0), "max y value for syst")
	("outtag",     value<std::string>()->default_value("default"), "tag name")
	("intag",     value<std::string>()->default_value("default"), "tag name")
	("run",     value<std::string>()->default_value("wp"), "process name")
	("xvar",    value<std::string>()->default_value("qtbyQ"), "variable x name")
      	("doA0",    bool_switch()->default_value(false), "")
	("doA1",    bool_switch()->default_value(false), "")
	("doA2",    bool_switch()->default_value(false), "")
	("doA3",    bool_switch()->default_value(false), "")
	("doA4",    bool_switch()->default_value(false), "")	
      	("interpolate",   bool_switch()->default_value(false), "")
	("savePdf",   bool_switch()->default_value(false), "")
	("savePdf2data",  bool_switch()->default_value(false), "")
	("saveJac",   bool_switch()->default_value(false), "")
	("saveSyst",  bool_switch()->default_value(false), "")
      	("runfit",   bool_switch()->default_value(false), "")
      	("verbose",   bool_switch()->default_value(false), "")
	("syst_scet",   bool_switch()->default_value(false), "")
	("syst_pdf",   bool_switch()->default_value(false), "")
	("syst_scale",   bool_switch()->default_value(false), "")
	("syst_altpdf",   bool_switch()->default_value(false), "")
	("shift_UL",   value<double>()->default_value(0.1), "max shift")
	("shift_A0",   value<double>()->default_value(0.1), "max shift")
	("shift_A1",   value<double>()->default_value(0.1), "max shift")
	("shift_A2",   value<double>()->default_value(0.1), "max shift")
	("shift_A3",   value<double>()->default_value(0.1), "max shift")
	("shift_A4",   value<double>()->default_value(0.1), "max shift")
	("debug",   bool_switch()->default_value(false), "");
      
      store(parse_command_line(argc, argv, desc), vm);
      notify(vm);
      if (vm.count("help")){
	std::cout << desc << '\n';
	return 0;
      }
    }
  catch (const error &ex)
    {
      std::cerr << ex.what() << '\n';
    }

  long nevents    = vm["nevents"].as<long>();
  int ntoys       = vm["ntoys"].as<int>();
  int extrabinsX  = vm["extrabinsX"].as<int>();
  int extrabinsY  = vm["extrabinsY"].as<int>();
  std::string outtag = vm["outtag"].as<std::string>();
  std::string intag = vm["intag"].as<std::string>();
  std::string run = vm["run"].as<std::string>();
  std::string xvar = vm["xvar"].as<std::string>();
  TString run2 = "";
  if(run=="wp")      run2="wplus";
  else if(run=="wm") run2="wminus";
  else if(run=="z")  run2="z";
  
  int verbose      = vm["verbose"].as<bool>();

  int syst_scet   = vm["syst_scet"].as<bool>();
  int syst_pdf    = vm["syst_pdf"].as<bool>();
  int syst_scale  = vm["syst_scale"].as<bool>();
  int syst_altpdf = vm["syst_altpdf"].as<bool>();
  
  int runfit       = vm["runfit"].as<bool>();
  int interpolate  = vm["interpolate"].as<bool>();
  int savePdf      = vm["savePdf"].as<bool>();
  int savePdf2data = vm["savePdf2data"].as<bool>();
  int saveJac     = vm["saveJac"].as<bool>();
  int saveSyst    = vm["saveSyst"].as<bool>();
  int debug       = vm["debug"].as<bool>();
  int doA0        = vm["doA0"].as<bool>();
  int doA1        = vm["doA1"].as<bool>();
  int doA2        = vm["doA2"].as<bool>();
  int doA3        = vm["doA3"].as<bool>();
  int doA4        = vm["doA4"].as<bool>();

  int dULx = vm["dULx"].as<int>();
  int dULy = vm["dULy"].as<int>();
  int dA0x   = vm["dA0x"].as<int>();
  int dA0y   = vm["dA0y"].as<int>();
  int dA1x   = vm["dA1x"].as<int>();
  int dA1y   = vm["dA1y"].as<int>();
  int dA2x   = vm["dA2x"].as<int>();
  int dA2y   = vm["dA2y"].as<int>();
  int dA3x   = vm["dA3x"].as<int>();
  int dA3y   = vm["dA3y"].as<int>();
  int dA4x   = vm["dA4x"].as<int>();
  int dA4y   = vm["dA4y"].as<int>();

  int cULx   = vm["cULx"].as<int>();
  int cULy   = vm["cULy"].as<int>();
  int cA0x   = vm["cA0x"].as<int>();
  int cA0y   = vm["cA0y"].as<int>();
  int cA1x   = vm["cA1x"].as<int>();
  int cA1y   = vm["cA1y"].as<int>();
  int cA2x   = vm["cA2x"].as<int>();
  int cA2y   = vm["cA2y"].as<int>();
  int cA3x   = vm["cA3x"].as<int>();
  int cA3y   = vm["cA3y"].as<int>();
  int cA4x   = vm["cA4x"].as<int>();
  int cA4y   = vm["cA4y"].as<int>();

  int fULx   = vm["fULx"].as<int>();
  int fULy   = vm["fULy"].as<int>();
  int fA0x   = vm["fA0x"].as<int>();
  int fA0y   = vm["fA0y"].as<int>();
  int fA1x   = vm["fA1x"].as<int>();
  int fA1y   = vm["fA1y"].as<int>();
  int fA2x   = vm["fA2x"].as<int>();
  int fA2y   = vm["fA2y"].as<int>();
  int fA3x   = vm["fA3x"].as<int>();
  int fA3y   = vm["fA3y"].as<int>();
  int fA4x   = vm["fA4x"].as<int>();
  int fA4y   = vm["fA4y"].as<int>();

  double x_max   = vm["x_max"].as<double>();
  double y_max   = vm["y_max"].as<double>();
  double xf_max  = vm["xf_max"].as<double>();
  double yf_max  = vm["yf_max"].as<double>();

  double shift_UL  = vm["shift_UL"].as<double>();
  double shift_A0  = vm["shift_A0"].as<double>();
  double shift_A1  = vm["shift_A1"].as<double>();
  double shift_A2  = vm["shift_A2"].as<double>();
  double shift_A3  = vm["shift_A3"].as<double>();
  double shift_A4  = vm["shift_A4"].as<double>();

  std::map<TString, double> shift_map;
  shift_map.insert( std::make_pair<TString, double >("UL", std::move(shift_UL)) );
  shift_map.insert( std::make_pair<TString, double >("A0", std::move(shift_A0)) );
  shift_map.insert( std::make_pair<TString, double >("A1", std::move(shift_A1)) );
  shift_map.insert( std::make_pair<TString, double >("A2", std::move(shift_A2)) );
  shift_map.insert( std::make_pair<TString, double >("A3", std::move(shift_A3)) );
  shift_map.insert( std::make_pair<TString, double >("A4", std::move(shift_A4)) );
  
  std::vector<TString> proc = {"UL"};
  if(doA0) proc.emplace_back("A0");
  if(doA1) proc.emplace_back("A1");
  if(doA2) proc.emplace_back("A2");
  if(doA3) proc.emplace_back("A3");
  if(doA4) proc.emplace_back("A4");

  if(runfit){
						    
    //TFile* fin_nom  = TFile::Open("root/file_qtbyQ_and_qt_vs_absy_v3.root", "READ");
    TFile* fin_nom  = TFile::Open(("/scratchnvme/tanmay/OutPut_2016/Final_Uses/Plot_root_Files_ang_coeff_"+xvar+"_2d/root_files_ang_coeff_"+xvar+"_2d.root").c_str(), "READ");
    if(fin_nom==0){
      cout << "Cannot find nominal file" << endl;
      return 0;
    }
    
    // scetlib variations on UL
    if(syst_scet){
      
      cout << "Doing syst_scet on UL" << endl;
      
      //TString hname =  run2+(xvar=="qtbyQ" ? "_ptqVgen" : "_ptVgen")+"_2d_absY_vs_"+TString(xvar.c_str())+"_differential_ul";
      TString hname =  "ul_"+TString(run.c_str())+"_2d_"+TString(xvar.c_str())+"_vs_absy";
      TH2D* h_nom = (TH2D*)fin_nom->Get(hname);    
      if(h_nom==0){
	cout << "Nominal histo not found. Continue." << endl;
	return 0;
      }
      
      // X-axis
      int X_nbins  = h_nom->GetXaxis()->GetNbins();
      double X_max = h_nom->GetXaxis()->GetXmax();
      if(xf_max>0.){
	X_nbins = h_nom->GetXaxis()->FindBin(xf_max);
	X_max = h_nom->GetXaxis()->GetBinUpEdge(X_nbins);
      }
      double X_edges[X_nbins+1];
      for(int ib = 0; ib<X_nbins; ib++) X_edges[ib] = h_nom->GetXaxis()->GetBinLowEdge(ib+1);
      X_edges[X_nbins] = h_nom->GetXaxis()->GetBinUpEdge(X_nbins);
      
      // Y-axis
      int Y_nbins  = h_nom->GetYaxis()->GetNbins();
      double Y_max = h_nom->GetYaxis()->GetXmax();
      if(yf_max>0.){
	Y_nbins = h_nom->GetYaxis()->FindBin(yf_max);
	Y_max = h_nom->GetYaxis()->GetBinUpEdge(Y_nbins);
      }
      double Y_edges[Y_nbins+1];
      for(int ib = 0; ib<Y_nbins; ib++) Y_edges[ib] = h_nom->GetYaxis()->GetBinLowEdge(ib+1);
      Y_edges[Y_nbins] = h_nom->GetYaxis()->GetBinUpEdge(Y_nbins) ;
      
      int nd = X_nbins*Y_nbins;
      cout << "Total bins in acceptance: " << X_nbins << "*" << Y_nbins << " = " << nd << endl;
      
      TH2D* hstart = new TH2D("h_start", h_nom->GetTitle(), X_nbins, X_edges, Y_nbins, Y_edges );
      TH1D* hinfo = new TH1D("h_info", "", 3, 0,3);
      
      TFile *fin_jac = TFile::Open(("fout_"+intag+".root").c_str(), "READ");
      if(fin_jac==0){
	cout << "Cannot find jac file" << endl;
	return 0;
      }
      cout << "Jac file " << fin_jac->GetName() << " opened.";
      TH1D* h_info = (TH1D*)fin_jac->Get("UL/h_info_UL");
      int nfpx = h_info->GetBinContent(8);
      int nfpy = h_info->GetBinContent(9);
      int np   = nfpx*nfpy;
      cout << " " << np << " parameters found" << endl;
      
      hinfo->SetBinContent(1, nd);
      hinfo->SetBinContent(2, np);
    
      MatrixXd J = MatrixXd::Zero(nd,np);
      MatrixXd V_inv_sqrt = MatrixXd::Zero(nd, nd);
      VectorXd y(nd);
      VectorXd y0(nd);
      VectorXd y_err(nd);
      
      int counter = 0;
      for(unsigned int idx=1; idx<=X_nbins; idx++){
	for(unsigned int idy=1; idy<=Y_nbins; idy++){
	  y0(counter)    = h_nom->GetBinContent(idx,idy);
	  y_err(counter) = h_nom->GetBinError(idx,idy);
	  for(unsigned int isyst=0; isyst<np; isyst++){
	    TH2D* h_jac = (TH2D*)fin_jac->Get(Form("UL/h_pdf2data_UL_jac%d", isyst));
	    J(counter,isyst) = h_jac->GetBinContent(idx,idy);
	  }
	  V_inv_sqrt(counter,counter) = 1.0/y_err(counter);
	  hstart->SetBinContent(idx,idy,  y0(counter));
	  hstart->SetBinError(idx,idy,  y_err(counter));
	  counter++;
	}
      }
						    
      TFile* fin_syst = TFile::Open(("/scratchnvme/tanmay/OutPut_2016/Final_Uses/Plot_root_Files_ul_"+xvar+"_scet_vars_2d/root_files_ul_"+xvar+"_scet_vars_2d.root").c_str(), "READ");
      if(fin_syst==0){
	cout << "Cannot find syst file" << endl;
	return 0;
      }
      
      vector<TString> scetlib_syst_names = {
	"pdf0",
	"omega_nu0.5",
	"c_nu-0.1-omega_nu0.5",
	"c_nu0.2-omega_nu0.5",
	"c_nu0.01-omega_nu0.5",
	"c_nu-0.2-omega_nu0.5",
	"c_nu-0.01",
	"Lambda20.25",
	"Lambda2-0.25",
	"Lambda4.01",
	"Lambda4.16",
	"Delta_Lambda2-0.02",
	"Delta_Lambda20.02",
	"c_nu-0.1-omega_nu0.707",
	"omega_nu0.707",
	"c_nu-0.1-omega_nu0.25-Omega0.12-Delta_Omega0.0158",
	"c_nu-0.1-omega_nu0.25-Omega0.13-Delta_Omega0.0",
	"gamma_cusp+1",
	"gamma_cusp+5",
	"gamma_mu_q+1",
	"gamma_mu_q+5",
	"gamma_nu+1",
	"gamma_nu+5",
	"h_qqV-0.5",
	"h_qqV-2.0",
	"s+1",
	"s+5",
	"b_qqV+1",
	"b_qqV+5",
	"b_qqbarV+1",
	"b_qqbarV+5",
	"b_qqS+1",
	"b_qqS+5",
	"b_qqDS+1",
	"b_qqDS+5",
	"b_qg+1",
	"b_qg+5",
	"transition_points0.4_0.75_1.1",
	"transition_points0.2_0.45_0.7",
	"transition_points0.4_0.55_0.7",
	"transition_points0.2_0.65_1.1"
      };
      
      hinfo->SetBinContent(3, scetlib_syst_names.size());

      TFile *fout = TFile::Open(("fout_fit_scet_UL_"+outtag+".root").c_str(), "RECREATE");
      
      for(unsigned int isyst=0; isyst<scetlib_syst_names.size(); isyst++){
	
	TString iproc = scetlib_syst_names[isyst];
	TH2D* hdummy = new TH2D(Form("hdummy_%d",isyst), iproc+";q_{T}/Q or q_{T};|y|", X_nbins, X_edges, Y_nbins, Y_edges );

	TH2D* h_syst = (TH2D*)fin_syst->Get("ul_"+run+"_2d_"+TString(xvar.c_str())+"_vs_absy_scetvar_"+iproc );    
	if(h_syst==0){
	  cout << "Histo not found. Continue." << endl;
	  continue;
	}
	cout << "Histo " << h_syst->GetName() << " found" << endl;
	
	counter = 0;
	for(unsigned int idx=1; idx<=X_nbins; idx++){
	  for(unsigned int idy=1; idy<=Y_nbins; idy++){
	    y(counter)  = h_syst->GetBinContent(idx,idy);
	    counter++;
	  }
	}
	
	MatrixXd A = V_inv_sqrt*J;
	MatrixXd b = V_inv_sqrt*(y-y0);
	VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	VectorXd delta = y-y0-J*x;
	VectorXd pull  = b-A*x;
	if(verbose) cout << delta << endl;
	MatrixXd chi2start = b.transpose()*b;
	MatrixXd chi2 = pull.transpose()*pull;
	double chi2startval = chi2start(0,0);
	double chi2val = chi2(0,0);
	int ndof = y.size() - np;
	cout <<  "Chi2 start = " << chi2startval << " --> chi2/ndof = " << chi2val << " / " << ndof << " = " << chi2val/ndof  << endl;

	TH2D* hdelta = (TH2D*)hdummy->Clone(Form("h_delta_%d",isyst));
	TH2D* hpull  = (TH2D*)hdummy->Clone(Form("h_pull_%d",isyst));
	TH2D* hratio = (TH2D*)hdummy->Clone(Form("h_ratio_%d",isyst));
	TH2D* hdata  = (TH2D*)hdummy->Clone(Form("h_data_%d",isyst));
	TH2D* hexp   = (TH2D*)hdummy->Clone(Form("h_exp_%d",isyst));

	int count_d = 0;
	for(unsigned int idx = 0; idx<X_nbins; idx++){
	  for(unsigned int idy = 0; idy<Y_nbins; idy++){	
	    hdelta->SetBinContent(idx+1,idy+1, delta(count_d));
	    hpull->SetBinContent(idx+1,idy+1,  pull(count_d));
	    hratio->SetBinContent(idx+1,idy+1, y(count_d)!=0. ? (y(count_d)-delta(count_d))/y(count_d) : 1.0 );	
	    hratio->SetBinError(idx+1,idy+1, y(count_d)!=0. ? y_err(count_d)/y(count_d) : 0.0 );	
	    hdata->SetBinContent(idx+1,idy+1,  y(count_d));
	    hdata->SetBinError(idx+1,idy+1,  y_err(count_d));
	    hexp->SetBinContent(idx+1,idy+1,  y(count_d)-delta(count_d));
	    //double exp_err = (J.row(count_d)*W*J.row(count_d).transpose())(0,0);
	    //hexp->SetBinError(idx+1,idy+1, exp_err>0. ? TMath::Sqrt(exp_err) : 0.0);
	    count_d++;
	  }
	}
	fout->cd();
	hpull->Write();
	hratio->Write();
	hdelta->Write();
	hdata->Write();
	hexp->Write();    
      }
      fin_syst->Close();
      
      fout->cd();
      hinfo->Write();
      hstart->Write();    
      fout->Close();
    }

    if(syst_pdf){

      TString pr = proc.size()>1 ?  proc[1] : proc[0];
      cout << "Doing syst_pdf on " << pr << endl;      

      TString hname = "ang_coeff_"+TString(run.c_str())+"_2d_"+TString(xvar.c_str())+"_vs_absy_A_"+TString(pr[1]);
      if(pr=="UL")
	hname =  "ul_"+TString(run.c_str())+"_2d_"+TString(xvar.c_str())+"_vs_absy";
      
      TH2D* h_nom = (TH2D*)fin_nom->Get(hname);    
      if(h_nom==0){
	cout << "Nominal histo not found. Continue." << endl;
	return 0;
      }
      
      // X-axis
      int X_nbins  = h_nom->GetXaxis()->GetNbins();
      double X_max = h_nom->GetXaxis()->GetXmax();
      if(xf_max>0.){
	X_nbins = h_nom->GetXaxis()->FindBin(xf_max);
	X_max = h_nom->GetXaxis()->GetBinUpEdge(X_nbins);
      }
      double X_edges[X_nbins+1];
      for(int ib = 0; ib<X_nbins; ib++) X_edges[ib] = h_nom->GetXaxis()->GetBinLowEdge(ib+1);
      X_edges[X_nbins] = h_nom->GetXaxis()->GetBinUpEdge(X_nbins);
      
      // Y-axis
      int Y_nbins  = h_nom->GetYaxis()->GetNbins();
      double Y_max = h_nom->GetYaxis()->GetXmax();
      if(yf_max>0.){
	Y_nbins = h_nom->GetYaxis()->FindBin(yf_max);
	Y_max = h_nom->GetYaxis()->GetBinUpEdge(Y_nbins);
      }
      double Y_edges[Y_nbins+1];
      for(int ib = 0; ib<Y_nbins; ib++) Y_edges[ib] = h_nom->GetYaxis()->GetBinLowEdge(ib+1);
      Y_edges[Y_nbins] = h_nom->GetYaxis()->GetBinUpEdge(Y_nbins) ;
      
      int nd = X_nbins*Y_nbins;
      cout << "Total bins in acceptance: " << X_nbins << "*" << Y_nbins << " = " << nd << endl;
      
      TH2D* hstart = new TH2D("h_start", h_nom->GetTitle(), X_nbins, X_edges, Y_nbins, Y_edges );
      TH1D* hinfo = new TH1D("h_info", "", 3, 0,3);
      
      TFile *fin_jac = TFile::Open(("fout_"+intag+".root").c_str(), "READ");
      if(fin_jac==0){
	cout << "Cannot find jac file" << endl;
	return 0;
      }
      cout << "Jac file " << fin_jac->GetName() << " opened.";
      TH1D* h_info = (TH1D*)fin_jac->Get(pr+"/h_info_"+pr);
      int nfpx = h_info->GetBinContent(8);
      int nfpy = h_info->GetBinContent(9);
      int np   = nfpx*nfpy;
      cout << " " << np << " parameters found" << endl;
      
      hinfo->SetBinContent(1, nd);
      hinfo->SetBinContent(2, np);
    
      MatrixXd J = MatrixXd::Zero(nd,np);
      MatrixXd V_inv_sqrt = MatrixXd::Zero(nd, nd);
      VectorXd y(nd);
      VectorXd y0(nd);
      VectorXd y_err(nd);
      
      int counter = 0;
      for(unsigned int idx=1; idx<=X_nbins; idx++){
	for(unsigned int idy=1; idy<=Y_nbins; idy++){
	  y0(counter)    = h_nom->GetBinContent(idx,idy);
	  y_err(counter) = h_nom->GetBinError(idx,idy);
	  for(unsigned int isyst=0; isyst<np; isyst++){	    
	    TH2D* h_jac = (TH2D*)fin_jac->Get(pr+"/h_pdf2data_"+pr+"_jac"+TString(Form("%d", isyst)));
	    J(counter,isyst) = h_jac->GetBinContent(idx,idy);
	  }
	  V_inv_sqrt(counter,counter) = 1.0/y_err(counter);
	  hstart->SetBinContent(idx,idy,  y0(counter));
	  hstart->SetBinError(idx,idy,  y_err(counter));
	  counter++;
	}
      }
      
      TFile* fin_syst = TFile::Open( ("/scratchnvme/tanmay/OutPut_2016/Final_Uses/Plot_root_Files_ang_coeff_"+xvar+"_pdf_msht20_vars_2d/root_files_ang_coeff_"+xvar+"_pdf_msht20_vars_2d.root").c_str(), "READ");
      if(fin_syst==0){
	cout << "Cannot find syst file" << endl;
	return 0;
      }

      vector<TString> pdf_syst_names = {
	"pdf0MSHT20"
      };
      for(int s = 1; s<=32; s++){
	pdf_syst_names.emplace_back(TString( Form("pdf%dMSHT20Down", s) ));
	pdf_syst_names.emplace_back(TString( Form("pdf%dMSHT20Up", s) ));
      }
      pdf_syst_names.emplace_back("as0116");
      pdf_syst_names.emplace_back("as0120");
      
      hinfo->SetBinContent(3, pdf_syst_names.size());

      TFile *fout = TFile::Open("fout_fit_pdf_"+pr+"_"+TString(outtag.c_str())+".root", "RECREATE");

      for(unsigned int isyst=0; isyst<pdf_syst_names.size(); isyst++){
	
	TString iproc = pdf_syst_names[isyst];
	TH2D* hdummy = new TH2D(Form("hdummy_%d",isyst), iproc+";q_{T}/Q or q_{T};|y|", X_nbins, X_edges, Y_nbins, Y_edges );

	TString hname = "ang_coeff_"+TString(run.c_str())+"_2d_"+TString(xvar.c_str())+"_vs_absy_"+iproc+"_A_"+TString(pr[1]);
	if(pr=="UL")
	  hname = "ul_"+TString(run.c_str())+"_2d_"+TString(xvar.c_str())+"_vs_absy_"+iproc;

	TH2D* h_syst = (TH2D*)fin_syst->Get( hname );    
	if(h_syst==0){
	  cout << "Histo " << hname << " not found. Continue." << endl;
	  continue;
	}
	cout << "Histo " << h_syst->GetName() << " found" << endl;

	counter = 0;
	for(unsigned int idx=1; idx<=X_nbins; idx++){
	  for(unsigned int idy=1; idy<=Y_nbins; idy++){
	    y(counter)  = h_syst->GetBinContent(idx,idy);
	    counter++;
	  }
	}
	
	MatrixXd A = V_inv_sqrt*J;
	MatrixXd b = V_inv_sqrt*(y-y0);
	VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	VectorXd delta = y-y0-J*x;
	VectorXd pull  = b-A*x;
	if(verbose) cout << delta << endl;
	MatrixXd chi2start = b.transpose()*b;
	MatrixXd chi2 = pull.transpose()*pull;
	double chi2startval = chi2start(0,0);
	double chi2val = chi2(0,0);
	int ndof = y.size() - np;
	cout <<  "Chi2 start = " << chi2startval << " --> chi2/ndof = " << chi2val << " / " << ndof << " = " << chi2val/ndof  << endl;

	TH2D* hdelta = (TH2D*)hdummy->Clone(Form("h_delta_%d",isyst));
	TH2D* hpull  = (TH2D*)hdummy->Clone(Form("h_pull_%d",isyst));
	TH2D* hratio = (TH2D*)hdummy->Clone(Form("h_ratio_%d",isyst));
	TH2D* hdata  = (TH2D*)hdummy->Clone(Form("h_data_%d",isyst));
	TH2D* hexp   = (TH2D*)hdummy->Clone(Form("h_exp_%d",isyst));

	int count_d = 0;
	for(unsigned int idx = 0; idx<X_nbins; idx++){
	  for(unsigned int idy = 0; idy<Y_nbins; idy++){	
	    hdelta->SetBinContent(idx+1,idy+1, delta(count_d));
	    hpull->SetBinContent(idx+1,idy+1,  pull(count_d));
	    hratio->SetBinContent(idx+1,idy+1, y(count_d)!=0. ? (y(count_d)-delta(count_d))/y(count_d) : 1.0 );	
	    hratio->SetBinError(idx+1,idy+1, y(count_d)!=0. ? y_err(count_d)/y(count_d) : 0.0 );	
	    hdata->SetBinContent(idx+1,idy+1,  y(count_d));
	    hdata->SetBinError(idx+1,idy+1,  y_err(count_d));
	    hexp->SetBinContent(idx+1,idy+1,  y(count_d)-delta(count_d));
	    //double exp_err = (J.row(count_d)*W*J.row(count_d).transpose())(0,0);
	    //hexp->SetBinError(idx+1,idy+1, exp_err>0. ? TMath::Sqrt(exp_err) : 0.0);
	    count_d++;
	  }
	}

	fout->cd();
	hpull->Write();
	hratio->Write();
	hdelta->Write();
	hdata->Write();
	hexp->Write();    
      }
      fin_syst->Close();
	    
      fout->cd();
      hinfo->Write();
      hstart->Write();    
      fout->Close();
    }

    if(syst_scale){

      TString pr = proc.size()>1 ?  proc[1] : proc[0];
      cout << "Doing syst_scale on " << pr << endl;      

      //TString hname =  "ang_coeff_"+TString(run.c_str())+"_"+TString(xvar.c_str())+"_vs_absy_A_"+TString(pr[1]);
      TString hname =  "ang_coeff_"+TString(run.c_str())+"_2d_"+TString(xvar.c_str())+"_vs_absy_A_"+TString(pr[1]);
      if(pr=="UL")
	hname =  "ul_"+TString(run.c_str())+"_2d_"+TString(xvar.c_str())+"_vs_absy";

      TH2D* h_nom = (TH2D*)fin_nom->Get(hname);    
      if(h_nom==0){
	cout << "Nominal histo not found. Continue." << endl;
	return 0;
      }
      
      // X-axis
      int X_nbins  = h_nom->GetXaxis()->GetNbins();
      double X_max = h_nom->GetXaxis()->GetXmax();
      if(xf_max>0.){
	X_nbins = h_nom->GetXaxis()->FindBin(xf_max);
	X_max = h_nom->GetXaxis()->GetBinUpEdge(X_nbins);
      }
      double X_edges[X_nbins+1];
      for(int ib = 0; ib<X_nbins; ib++) X_edges[ib] = h_nom->GetXaxis()->GetBinLowEdge(ib+1);
      X_edges[X_nbins] = h_nom->GetXaxis()->GetBinUpEdge(X_nbins);
      
      // Y-axis
      int Y_nbins  = h_nom->GetYaxis()->GetNbins();
      double Y_max = h_nom->GetYaxis()->GetXmax();
      if(yf_max>0.){
	Y_nbins = h_nom->GetYaxis()->FindBin(yf_max);
	Y_max = h_nom->GetYaxis()->GetBinUpEdge(Y_nbins);
      }
      double Y_edges[Y_nbins+1];
      for(int ib = 0; ib<Y_nbins; ib++) Y_edges[ib] = h_nom->GetYaxis()->GetBinLowEdge(ib+1);
      Y_edges[Y_nbins] = h_nom->GetYaxis()->GetBinUpEdge(Y_nbins) ;
      
      int nd = X_nbins*Y_nbins;
      cout << "Total bins in acceptance: " << X_nbins << "*" << Y_nbins << " = " << nd << endl;
      
      TH2D* hstart = new TH2D("h_start", h_nom->GetTitle(), X_nbins, X_edges, Y_nbins, Y_edges );
      TH1D* hinfo = new TH1D("h_info", "", 3, 0,3);
      
      TFile *fin_jac = TFile::Open(("fout_"+intag+".root").c_str(), "READ");
      if(fin_jac==0){
	cout << "Cannot find jac file" << endl;
	return 0;
      }
      cout << "Jac file " << fin_jac->GetName() << " opened.";
      TH1D* h_info = (TH1D*)fin_jac->Get(pr+"/h_info_"+pr);
      int nfpx = h_info->GetBinContent(8);
      int nfpy = h_info->GetBinContent(9);
      int np   = nfpx*nfpy;
      cout << " " << np << " parameters found" << endl;
      
      hinfo->SetBinContent(1, nd);
      hinfo->SetBinContent(2, np);
    
      MatrixXd J = MatrixXd::Zero(nd,np);
      MatrixXd V_inv_sqrt = MatrixXd::Zero(nd, nd);
      VectorXd y(nd);
      VectorXd y0(nd);
      VectorXd y_err(nd);
      
      int counter = 0;
      for(unsigned int idx=1; idx<=X_nbins; idx++){
	for(unsigned int idy=1; idy<=Y_nbins; idy++){
	  y0(counter)    = h_nom->GetBinContent(idx,idy);
	  y_err(counter) = h_nom->GetBinError(idx,idy);
	  for(unsigned int isyst=0; isyst<np; isyst++){
	    TH2D* h_jac = (TH2D*)fin_jac->Get(pr+"/h_pdf2data_"+pr+"_jac"+TString(Form("%d", isyst)));
	    J(counter,isyst) = h_jac->GetBinContent(idx,idy);
	  }
	  V_inv_sqrt(counter,counter) = 1.0/y_err(counter);
	  hstart->SetBinContent(idx,idy,  y0(counter));
	  hstart->SetBinError(idx,idy,  y_err(counter));
	  counter++;
	}
      }

      TFile* fin_syst = TFile::Open(("/scratchnvme/tanmay/OutPut_2016/Final_Uses/Plot_root_Files_ang_coeff_"+xvar+"_qcd_vars_2d/root_files_ang_coeff_"+xvar+"_qcd_vars_2d.root").c_str(), "READ");
      if(fin_syst==0){
	cout << "Cannot find syst file" << endl;
	return 0;
      }
      
      vector<TString> qcd_syst_names = {
	"oneone",
	"twotwo",
	"point5point5",
	"onepoint5",
	"point5one",
	"twoone",
	"onetwo"
      };
      
      hinfo->SetBinContent(3, qcd_syst_names.size());

      TFile *fout = TFile::Open("fout_fit_scale_"+pr+"_"+TString(outtag.c_str())+".root", "RECREATE");
      
      for(unsigned int isyst=0; isyst<qcd_syst_names.size(); isyst++){
	
	TString iproc = qcd_syst_names[isyst];
	TH2D* hdummy = new TH2D(Form("hdummy_%d",isyst), iproc+";q_{T}/Q or q_{T};|y|", X_nbins, X_edges, Y_nbins, Y_edges );

	TString hname = "ang_coeff_"+TString(run.c_str())+"_2d_"+TString(xvar.c_str())+"_vs_absy_"+iproc+"_A_"+TString(pr[1]);
	if(pr=="UL")
	  hname = "ul_"+TString(run.c_str())+"_2d_"+TString(xvar.c_str())+"_vs_absy_"+iproc;

	TH2D* h_syst = (TH2D*)fin_syst->Get(hname);
	if(h_syst==0){
	  cout << "Histo not found. Continue." << endl;
	  continue;
	}
	cout << "Histo " << h_syst->GetName() << " found" << endl;
	
	counter = 0;
	for(unsigned int idx=1; idx<=X_nbins; idx++){
	  for(unsigned int idy=1; idy<=Y_nbins; idy++){
	    y(counter)  = h_syst->GetBinContent(idx,idy);
	    counter++;
	  }
	}
	
	MatrixXd A = V_inv_sqrt*J;
	MatrixXd b = V_inv_sqrt*(y-y0);
	VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	VectorXd delta = y-y0-J*x;
	VectorXd pull  = b-A*x;
	if(verbose) cout << delta << endl;
	MatrixXd chi2start = b.transpose()*b;
	MatrixXd chi2 = pull.transpose()*pull;
	double chi2startval = chi2start(0,0);
	double chi2val = chi2(0,0);
	int ndof = y.size() - np;
	cout <<  "Chi2 start = " << chi2startval << " --> chi2/ndof = " << chi2val << " / " << ndof << " = " << chi2val/ndof  << endl;

	TH2D* hdelta = (TH2D*)hdummy->Clone(Form("h_delta_%d",isyst));
	TH2D* hpull  = (TH2D*)hdummy->Clone(Form("h_pull_%d",isyst));
	TH2D* hratio = (TH2D*)hdummy->Clone(Form("h_ratio_%d",isyst));
	TH2D* hdata  = (TH2D*)hdummy->Clone(Form("h_data_%d",isyst));
	TH2D* hexp   = (TH2D*)hdummy->Clone(Form("h_exp_%d",isyst));

	int count_d = 0;
	for(unsigned int idx = 0; idx<X_nbins; idx++){
	  for(unsigned int idy = 0; idy<Y_nbins; idy++){	
	    hdelta->SetBinContent(idx+1,idy+1, delta(count_d));
	    hpull->SetBinContent(idx+1,idy+1,  pull(count_d));
	    hratio->SetBinContent(idx+1,idy+1, y(count_d)!=0. ? (y(count_d)-delta(count_d))/y(count_d) : 1.0 );	
	    hratio->SetBinError(idx+1,idy+1, y(count_d)!=0. ? y_err(count_d)/y(count_d) : 0.0 );	
	    hdata->SetBinContent(idx+1,idy+1,  y(count_d));
	    hdata->SetBinError(idx+1,idy+1,  y_err(count_d));
	    hexp->SetBinContent(idx+1,idy+1,  y(count_d)-delta(count_d));
	    //double exp_err = (J.row(count_d)*W*J.row(count_d).transpose())(0,0);
	    //hexp->SetBinError(idx+1,idy+1, exp_err>0. ? TMath::Sqrt(exp_err) : 0.0);
	    count_d++;
	  }
	}
	fout->cd();
	hpull->Write();
	hratio->Write();
	hdelta->Write();
	hdata->Write();
	hexp->Write();    
      }
      fin_syst->Close();
      
      fout->cd();
      hinfo->Write();
      hstart->Write();    
      fout->Close();
    }

    if(syst_altpdf){

      TString pr = proc.size()>1 ?  proc[1] : proc[0];
      cout << "Doing syst_altpdf on " << pr << endl;      

      TString hname = "ang_coeff_"+TString(run.c_str())+"_2d_"+TString(xvar.c_str())+"_vs_absy_A_"+TString(pr[1]);
      if(pr=="UL")
	hname =  "ul_"+TString(run.c_str())+"_2d_"+TString(xvar.c_str())+"_vs_absy";

      TH2D* h_nom = (TH2D*)fin_nom->Get(hname);    
      if(h_nom==0){
	cout << "Nominal histo not found. Continue." << endl;
	return 0;
      }
      
      // X-axis
      int X_nbins  = h_nom->GetXaxis()->GetNbins();
      double X_max = h_nom->GetXaxis()->GetXmax();
      if(xf_max>0.){
	X_nbins = h_nom->GetXaxis()->FindBin(xf_max);
	X_max = h_nom->GetXaxis()->GetBinUpEdge(X_nbins);
      }
      double X_edges[X_nbins+1];
      for(int ib = 0; ib<X_nbins; ib++) X_edges[ib] = h_nom->GetXaxis()->GetBinLowEdge(ib+1);
      X_edges[X_nbins] = h_nom->GetXaxis()->GetBinUpEdge(X_nbins);
      
      // Y-axis
      int Y_nbins  = h_nom->GetYaxis()->GetNbins();
      double Y_max = h_nom->GetYaxis()->GetXmax();
      if(yf_max>0.){
	Y_nbins = h_nom->GetYaxis()->FindBin(yf_max);
	Y_max = h_nom->GetYaxis()->GetBinUpEdge(Y_nbins);
      }
      double Y_edges[Y_nbins+1];
      for(int ib = 0; ib<Y_nbins; ib++) Y_edges[ib] = h_nom->GetYaxis()->GetBinLowEdge(ib+1);
      Y_edges[Y_nbins] = h_nom->GetYaxis()->GetBinUpEdge(Y_nbins) ;
      
      int nd = X_nbins*Y_nbins;
      cout << "Total bins in acceptance: " << X_nbins << "*" << Y_nbins << " = " << nd << endl;
      
      TH2D* hstart = new TH2D("h_start", h_nom->GetTitle(), X_nbins, X_edges, Y_nbins, Y_edges );
      TH1D* hinfo = new TH1D("h_info", "", 3, 0,3);
      
      TFile *fin_jac = TFile::Open(("fout_"+intag+".root").c_str(), "READ");
      if(fin_jac==0){
	cout << "Cannot find jac file" << endl;
	return 0;
      }
      cout << "Jac file " << fin_jac->GetName() << " opened.";
      TH1D* h_info = (TH1D*)fin_jac->Get(pr+"/h_info_"+pr);
      int nfpx = h_info->GetBinContent(8);
      int nfpy = h_info->GetBinContent(9);
      int np   = nfpx*nfpy;
      cout << " " << np << " parameters found" << endl;
      
      hinfo->SetBinContent(1, nd);
      hinfo->SetBinContent(2, np);
    
      MatrixXd J = MatrixXd::Zero(nd,np);
      MatrixXd V_inv_sqrt = MatrixXd::Zero(nd, nd);
      VectorXd y(nd);
      VectorXd y0(nd);
      VectorXd y_err(nd);
      
      int counter = 0;
      for(unsigned int idx=1; idx<=X_nbins; idx++){
	for(unsigned int idy=1; idy<=Y_nbins; idy++){
	  y0(counter)    = h_nom->GetBinContent(idx,idy);
	  y_err(counter) = h_nom->GetBinError(idx,idy);
	  for(unsigned int isyst=0; isyst<np; isyst++){	    
	    TH2D* h_jac = (TH2D*)fin_jac->Get(pr+"/h_pdf2data_"+pr+"_jac"+TString(Form("%d", isyst)));
	    J(counter,isyst) = h_jac->GetBinContent(idx,idy);
	  }
	  V_inv_sqrt(counter,counter) = 1.0/y_err(counter);
	  hstart->SetBinContent(idx,idy,  y0(counter));
	  hstart->SetBinError(idx,idy,  y_err(counter));
	  counter++;
	}
      }

      vector<TString> pdf_alt_names = {
	"ct18",
	"nnpdf40",
	"mmht",
	"ct18z",
	//"nnpdf30",
	"nnpdf31",
	"atlasWZj20",
	"pdf4lhc21"	
      };

      map<TString, TString> altpdf_nom;
      altpdf_nom.insert( std::make_pair<TString, string >("ct18", "pdf0CT18") );
      altpdf_nom.insert( std::make_pair<TString, string >("nnpdf40", "pdf1NNPDF40") );
      altpdf_nom.insert( std::make_pair<TString, string >("mmht", "pdf0MMHT") );
      altpdf_nom.insert( std::make_pair<TString, string >("ct18z", "pdf0CT18Z") );
      //altpdf_nom.insert( std::make_pair<TString, string >("nnpdf30", "pdf1NNPDF30") );
      altpdf_nom.insert( std::make_pair<TString, string >("nnpdf31", "pdf1NNPDF31") );
      altpdf_nom.insert( std::make_pair<TString, string >("atlasWZj20", "pdf0ATLASWZJ20") );
      altpdf_nom.insert( std::make_pair<TString, string >("pdf4lhc21", "pdf1PDF4LHC21") );
      
      hinfo->SetBinContent(3, pdf_alt_names.size());

      TFile *fout = TFile::Open("fout_fit_altpdf_"+pr+"_"+TString(outtag.c_str())+".root", "RECREATE");

      for(unsigned int isyst=0; isyst<pdf_alt_names.size(); isyst++){

	TString iproc = pdf_alt_names[isyst];	
	TFile* fin_syst = TFile::Open( "/scratchnvme/tanmay/OutPut_2016/Final_Uses/Plot_root_Files_ang_coeff_"+TString(xvar.c_str())+"_pdf_"+iproc+"_vars_2d/root_files_ang_coeff_"+TString(xvar.c_str())+"_pdf_"+iproc+"_vars_2d.root", "READ");
	if(fin_syst==0){
	  cout << "Cannot find syst file" << endl;
	  continue;
	}

	TH2D* hdummy = new TH2D(Form("hdummy_%d",isyst), iproc+";q_{T}/Q or q_{T};|y|", X_nbins, X_edges, Y_nbins, Y_edges );

	TString hname = "ang_coeff_"+TString(run.c_str())+"_2d_"+TString(xvar.c_str())+"_vs_absy_"+altpdf_nom[iproc]+"_A_"+TString(pr[1]);
	if(pr=="UL")
	  hname = "ul_"+TString(run.c_str())+"_2d_"+TString(xvar.c_str())+"_vs_absy_"+altpdf_nom[iproc];

	TH2D* h_syst = (TH2D*)fin_syst->Get( hname );    
	if(h_syst==0){
	  cout << "Histo not found. Continue." << endl;
	  continue;
	}
	cout << "Histo " << h_syst->GetName() << " found" << endl;

	counter = 0;
	for(unsigned int idx=1; idx<=X_nbins; idx++){
	  for(unsigned int idy=1; idy<=Y_nbins; idy++){
	    y(counter)  = h_syst->GetBinContent(idx,idy);
	    counter++;
	  }
	}
	
	MatrixXd A = V_inv_sqrt*J;
	MatrixXd b = V_inv_sqrt*(y-y0);
	VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	VectorXd delta = y-y0-J*x;
	VectorXd pull  = b-A*x;
	if(verbose) cout << delta << endl;
	MatrixXd chi2start = b.transpose()*b;
	MatrixXd chi2 = pull.transpose()*pull;
	double chi2startval = chi2start(0,0);
	double chi2val = chi2(0,0);
	int ndof = y.size() - np;
	cout <<  "Chi2 start = " << chi2startval << " --> chi2/ndof = " << chi2val << " / " << ndof << " = " << chi2val/ndof  << endl;

	TH2D* hdelta = (TH2D*)hdummy->Clone(Form("h_delta_%d",isyst));
	TH2D* hpull  = (TH2D*)hdummy->Clone(Form("h_pull_%d",isyst));
	TH2D* hratio = (TH2D*)hdummy->Clone(Form("h_ratio_%d",isyst));
	TH2D* hdata  = (TH2D*)hdummy->Clone(Form("h_data_%d",isyst));
	TH2D* hexp   = (TH2D*)hdummy->Clone(Form("h_exp_%d",isyst));

	int count_d = 0;
	for(unsigned int idx = 0; idx<X_nbins; idx++){
	  for(unsigned int idy = 0; idy<Y_nbins; idy++){	
	    hdelta->SetBinContent(idx+1,idy+1, delta(count_d));
	    hpull->SetBinContent(idx+1,idy+1,  pull(count_d));
	    hratio->SetBinContent(idx+1,idy+1, y(count_d)!=0. ? (y(count_d)-delta(count_d))/y(count_d) : 1.0 );	
	    hratio->SetBinError(idx+1,idy+1, y(count_d)!=0. ? y_err(count_d)/y(count_d) : 0.0 );	
	    hdata->SetBinContent(idx+1,idy+1,  y(count_d));
	    hdata->SetBinError(idx+1,idy+1,  y_err(count_d));
	    hexp->SetBinContent(idx+1,idy+1,  y(count_d)-delta(count_d));
	    //double exp_err = (J.row(count_d)*W*J.row(count_d).transpose())(0,0);
	    //hexp->SetBinError(idx+1,idy+1, exp_err>0. ? TMath::Sqrt(exp_err) : 0.0);
	    count_d++;
	  }
	}

	fout->cd();
	hpull->Write();
	hratio->Write();
	hdelta->Write();
	hdata->Write();
	hexp->Write();    
	fin_syst->Close();
      }
	    
      fout->cd();
      hinfo->Write();
      hstart->Write();    
      fout->Close();
    }
   
    fin_nom->Close();
    return 0;
  }

  TFile *fout = TFile::Open(("fout_"+outtag+".root").c_str(), "RECREATE");
  
  // dummy file
  if(debug){
    TFile* f = TFile::Open( "root/file_qtbyQ_and_qt_vs_absy_v3_debug.root", "RECREATE");  
    for(auto& pr : proc){
      TString hname = pr=="UL" ?
	run2+(xvar=="qtbyQ" ? "_ptqVgen" : "_ptVgen")+"_2d_absY_vs_"+TString(xvar.c_str())+"_differential_ul" :
	"ang_coeff_"+TString(run.c_str())+"_"+TString(xvar.c_str())+"_vs_absy_A_"+TString(pr[1]);
      
      TH2D* h = new TH2D(hname, "", 20, 0.0, 0.5, 20, 0.0, 2.5 );
      for(int ibx=1; ibx<=h->GetXaxis()->GetNbins() ; ibx++ ){
	for(int iby=1; iby<=h->GetYaxis()->GetNbins() ; iby++ ){
	  if(pr=="A1" || pr=="A3"){
	    h->SetBinContent(ibx,iby, 10000.*iby/h->GetYaxis()->GetNbins()*ibx/h->GetXaxis()->GetNbins());
	  }
	  else if(pr=="UL" || pr=="A0" || pr=="A2"){
	    h->SetBinContent(ibx,iby, 10000.*ibx/h->GetXaxis()->GetNbins());
	    //h->SetBinContent(ibx,iby, 10000.);
	  }
	  else if(pr=="A4"){
	    h->SetBinContent(ibx,iby, 10000.*iby/h->GetYaxis()->GetNbins());
	  }
	  h->SetBinError(ibx,iby, TMath::Sqrt( h->GetBinContent(ibx,iby)));
	}
      }
      h->Write();
    }
    f->Close();
  }

  //TFile* fin = TFile::Open(TString("root/file_qtbyQ_and_qt_vs_absy_v3")+(debug ? "_debug.root" : ".root"), "READ");
  TFile* fin  = TFile::Open(("/scratchnvme/tanmay/OutPut_2016/Final_Uses/Plot_root_Files_ang_coeff_"+xvar+"_2d/root_files_ang_coeff_"+xvar+"_2d.root").c_str(), "READ");
  if(fin==0 || fin==nullptr || fin->IsZombie()){
    cout << "File NOT found" << endl;
    return 0;
  }

  fin->cd();
  
  std::map<TString, std::array<int,2> > deg_map;
  std::map<TString, std::array<int,2> > par_map;
  std::map<TString, std::array<int,2> > ctr_map;
  std::map<TString, std::array<int,2> > degf_map;

  //UL
  deg_map.insert  ( std::make_pair<TString, std::array<int,2> >("UL", {dULx,  dULy}) );
  par_map.insert  ( std::make_pair<TString, std::array<int,2> >("UL", {0, +1}) );
  ctr_map.insert  ( std::make_pair<TString, std::array<int,2> >("UL", {cULx,  cULy}) ); 
  //A0
  deg_map.insert  ( std::make_pair<TString, std::array<int,2> >("A0", {dA0x,  dA0y}) );
  par_map.insert  ( std::make_pair<TString, std::array<int,2> >("A0", {0, +1}) );
  ctr_map.insert  ( std::make_pair<TString, std::array<int,2> >("A0", {cA0x,  cA0y}) ); 
  //A1
  deg_map.insert  ( std::make_pair<TString, std::array<int,2> >("A1", {dA1x,  dA1y}) );
  par_map.insert  ( std::make_pair<TString, std::array<int,2> >("A1", {0, -1}) );
  ctr_map.insert  ( std::make_pair<TString, std::array<int,2> >("A1", {cA1x,  cA1y}) ); 
  //A2
  deg_map.insert  ( std::make_pair<TString, std::array<int,2> >("A2", {dA2x,  dA2y}) );
  par_map.insert  ( std::make_pair<TString, std::array<int,2> >("A2", {0, +1}) );
  ctr_map.insert  ( std::make_pair<TString, std::array<int,2> >("A2", {cA2x,  cA2y}) ); 
  //A3
  deg_map.insert  ( std::make_pair<TString, std::array<int,2> >("A3", {dA3x,  dA3y}) );
  par_map.insert  ( std::make_pair<TString, std::array<int,2> >("A3", {0, +1}) );
  ctr_map.insert  ( std::make_pair<TString, std::array<int,2> >("A3", {cA3x,  cA3y}) ); 
  //A4
  deg_map.insert  ( std::make_pair<TString, std::array<int,2> >("A4", {dA4x,  dA4y}) );
  par_map.insert  ( std::make_pair<TString, std::array<int,2> >("A4", {0, -1}) );
  ctr_map.insert  ( std::make_pair<TString, std::array<int,2> >("A4", {cA4x,  cA4y}) );

  degf_map.insert  ( std::make_pair<TString, std::array<int,2> >("UL", {fULx,  fULy}) );
  degf_map.insert  ( std::make_pair<TString, std::array<int,2> >("A0", {fA0x,  fA0y}) );
  degf_map.insert  ( std::make_pair<TString, std::array<int,2> >("A1", {fA1x,  fA1y}) );
  degf_map.insert  ( std::make_pair<TString, std::array<int,2> >("A2", {fA2x,  fA2y}) );
  degf_map.insert  ( std::make_pair<TString, std::array<int,2> >("A3", {fA3x,  fA3y}) );
  degf_map.insert  ( std::make_pair<TString, std::array<int,2> >("A4", {fA4x,  fA4y}) );

  TH2D* hpdf_UL = 0;
  TH2D* h_interpol_UL = 0;
  TH2D* hpdf2data_UL = 0;
  
  // loop over all histos
  for(unsigned int i=0; i <proc.size(); i++){

    TString iproc = proc[i];
    cout << "Doing proc " << iproc << endl;
    //TString hname = iproc=="UL" ?
    //run2+(xvar=="qtbyQ" ? "_ptqVgen" : "_ptVgen")+"_2d_absY_vs_"+TString(xvar.c_str())+"_differential_ul" :
    //"ang_coeff_"+TString(run.c_str())+"_qtbyQ_vs_absy_A_"+TString(iproc[1]);
    TString hname = iproc=="UL" ?
      "ul_"+TString(run.c_str())+"_2d_"+TString(xvar.c_str())+"_vs_absy" :
      "ang_coeff_"+TString(run.c_str())+"_2d_"+TString(xvar.c_str())+"_vs_absy_A_"+TString(iproc[1]);
    
    TH2D* h = (TH2D*)fin->Get(hname);    
    if(h==0){
      cout << "Histo not found. Continue." << endl;
      continue;
    }
    cout << "Histo " << h->GetName() << " found" << endl;

    if(iproc=="UL"){
      h->Scale(1./ h->Integral());
      if(interpolate){
	h_interpol_UL = (TH2D*)h->Clone("h_interpol_UL");
	for(int ibx=1; ibx<=h->GetXaxis()->GetNbins(); ibx++){
	  for(int iby=1; iby<=h->GetYaxis()->GetNbins(); iby++){
	    h_interpol_UL->SetBinContent(ibx,iby, h->GetBinContent(ibx,iby)/
					 (h->GetXaxis()->GetBinWidth(ibx)*
					  h->GetYaxis()->GetBinWidth(iby)));
	  }
	}
      }
    }
    
    // X-axis
    int X_nbins  = h->GetXaxis()->GetNbins();
    double X_min = h->GetXaxis()->GetXmin();
    double X_max = h->GetXaxis()->GetXmax();
    if(x_max>0.){
      X_nbins = h->GetXaxis()->FindBin(x_max);
      X_max = h->GetXaxis()->GetBinUpEdge(X_nbins);
    }
    double X_scale  = (X_max-X_min)*0.5;
    double X_offset = (X_max+X_min)*0.5;
    double X_edges[X_nbins+1];
    for(int ib = 0; ib<X_nbins; ib++) X_edges[ib] = h->GetXaxis()->GetBinLowEdge(ib+1);
    X_edges[X_nbins] = h->GetXaxis()->GetBinUpEdge(X_nbins);
    if(verbose) cout << "X_edges: ";
    if(verbose) for(int i = 0; i <= X_nbins; i++) cout << X_edges[i] << "," ;
    if(verbose) cout << endl;
    
    // Y-axis
    int Y_nbins  = h->GetYaxis()->GetNbins();
    double Y_max = h->GetYaxis()->GetXmax();
    if(y_max>0.){
      Y_nbins = h->GetYaxis()->FindBin(y_max);
      Y_max = h->GetYaxis()->GetBinUpEdge(Y_nbins);
    }
    double Y_min = -Y_max; // assume |y| is plotted
    double Y_scale  = (Y_max-Y_min)*0.5;
    double Y_offset = (Y_max+Y_min)*0.5;
    double Y_edges[Y_nbins+1];
    for(int ib = 0; ib<Y_nbins; ib++) Y_edges[ib] = h->GetYaxis()->GetBinLowEdge(ib+1);
    Y_edges[Y_nbins] = h->GetYaxis()->GetBinUpEdge(Y_nbins) ;
    if(verbose) cout << "Y_edges: ";
    if(verbose) for(int i = 0; i <= Y_nbins; i++) cout << Y_edges[i] << "," ;
    if(verbose) cout << endl;
    
    assert( ctr_map[iproc].at(1)==0 || (ctr_map[iproc].at(1) && deg_map[iproc].at(1)%2==0));

    TH2D* hdummy = new TH2D("hdummy_"+iproc, iproc+";q_{T}/Q or q_{T};|y|", X_nbins, X_edges, Y_nbins, Y_edges );
    
    int nd = X_nbins*Y_nbins;
    int npx = deg_map[iproc].at(0) + 1 - ctr_map[iproc].at(0);
    int npy = par_map[iproc].at(1)>0 ? (deg_map[iproc].at(1)/2 + 1) : (deg_map[iproc].at(1) + 1)/2;
    if( ctr_map[iproc].at(1) ) npy--;
    int np = npx*npy;
    if(verbose) cout << "np = " << npx << " * " << npy << " = " << np << endl;

    int nfpx = degf_map[iproc].at(0) + 1;
    int nfpy = degf_map[iproc].at(1)/2 + 1;
    //int nfp = nfpx*nfpy;
    
    TF1* cheb_x = new TF1("cheb_x", cheb_fct, X_edges[0], X_edges[X_nbins], 4); 
    cheb_x->SetParNames("n","offset","scale","m");
    
    TF1* cheb_y = new TF1("cheb_y", cheb_fct, Y_edges[0], Y_edges[Y_nbins], 4); 
    cheb_y->SetParNames("n","offset","scale","m");

    cheb_x->SetParameter("offset", X_offset);
    cheb_x->SetParameter("scale",  X_scale);
    cheb_y->SetParameter("offset", Y_offset);
    cheb_y->SetParameter("scale",  Y_scale);

    cheb_x->SetParameter("n", deg_map[iproc].at(0) );
    cheb_y->SetParameter("n", deg_map[iproc].at(1) );
    
    int count_d=0;
    VectorXd y(nd);
    VectorXd y_err(nd);
    MatrixXd J(nd,np);

    std::vector<TString> p_names = {};

    // here we integrate over the bin
    if(iproc=="UL"){

      // loop over data bins
      cout << "Start the loop over X..." << endl;
      MatrixXd integrals_x = MatrixXd::Zero(X_nbins, deg_map[iproc].at(0) + 1);
      for(unsigned int idx = 0; idx<X_nbins; idx++){
	for(unsigned int ipx = 0; ipx< (deg_map[iproc].at(0) + 1); ipx++){
	  // constraint at 0
	  if( ctr_map[iproc].at(0) && ipx==0){
	    continue;
	  }	  
	  cheb_x->SetParameter("m", ipx);
	  double totx = cheb_x->Integral( X_edges[idx], X_edges[idx+1],1.e-12 );
	  if(verbose) cout << "\tmX=" << ipx << ": int_x [" << X_edges[idx] << "," << X_edges[idx+1] << "] = " << totx << endl;
	  integrals_x(idx, ipx) = totx;
	}       
      }
      
      cout << "Start the loop over Y..." << endl;
      MatrixXd integrals_y = MatrixXd::Zero(Y_nbins, deg_map[iproc].at(1) + 1);
      for(unsigned int idy = 0; idy<Y_nbins; idy++){
	for(unsigned int ipy = 0; ipy<deg_map[iproc].at(1) + 1; ipy++){	    
	  cheb_y->SetParameter("m", ipy);
	  double cheby1 = cheb_y->Integral( Y_edges[idy], Y_edges[idy+1],1.e-12 ); 
	  cheb_y->SetParameter("m", deg_map[iproc].at(1) - ipy);
	  double cheby2 = cheb_y->Integral( Y_edges[idy], Y_edges[idy+1],1.e-12 ); 	    
	  // even function:
	  if( par_map[iproc].at(1)>0 && ipy < (deg_map[iproc].at(1)/2 + 1 - ctr_map[iproc].at(1)) ){
	    double toty = cheby1+cheby2 ;
	    if(verbose) cout << "\tmY=[" << ipy << "," << deg_map[iproc].at(1) - ipy << "] : int_y [" << Y_edges[idy] << "," << Y_edges[idy+1] << "] = (" << cheby1  << " + " << cheby2 << ")";
	    if( deg_map[iproc].at(1)%2==0 && ipy==deg_map[iproc].at(1)/2 ){
	      toty *= 0.5;
	      if(verbose) cout << "/2" << endl;
	    }
	    else{
	      if(verbose) cout << endl;
	    }
	    integrals_y(idy,ipy) = toty;
	  }
	  // odd function
	  else if( par_map[iproc].at(1)<0 && ipy < (deg_map[iproc].at(1)+1)/2){
	    double toty = cheby1 - cheby2 ;
	    if(verbose) cout << "\tmY=[" << ipy << "," << deg_map[iproc].at(1) - ipy << "] : int_y [" << Y_edges[idy] << "," << Y_edges[idy+1] << "] = (" << cheby1  << " - " << cheby2 << ")";
	    integrals_y(idy,ipy) = toty;
	  }
	}	
      }
      
      // loop over data bins and fill outer product
      cout << "Start the final loop..." << endl;
      for(unsigned int idx = 0; idx<X_nbins; idx++){
	for(unsigned int idy = 0; idy<Y_nbins; idy++){
	  if(verbose) cout << "Doing bin number... " << count_d << endl;
	  y(count_d) = h->GetBinContent(idx+1, idy+1);
	  y_err(count_d) = h->GetBinError(idx+1, idy+1);
	  // loop over parameters
	  int count_p=0;
	  for(unsigned int ipx = 0; ipx<(deg_map[iproc].at(0) + 1); ipx++){
	    // constraint at 0
	    if( ctr_map[iproc].at(0) && ipx==0){
	      continue;
	    }	  	  
	    for(unsigned int ipy = 0; ipy<deg_map[iproc].at(1) + 1; ipy++){
	      // even function:
	      if( (par_map[iproc].at(1)>0 && ipy < (deg_map[iproc].at(1)/2 + 1 - ctr_map[iproc].at(1))) 
		  || (par_map[iproc].at(1)<0 && ipy < (deg_map[iproc].at(1)+1)/2)){
		J(count_d, count_p) = integrals_x(idx,ipx)*integrals_y(idy,ipy);
		if(verbose) cout << "jac(" << count_d << "," << count_p << ") = " << integrals_x(idx,ipx) << " * " << integrals_y(idy,ipy) << " = " << J(count_d, count_p) << endl;
		p_names.emplace_back( TString(Form("f_%d_%d", ipx, ipy)) ); 
		count_p++;
	      }
	    }
	  }
	  count_d++;
	}
      }
    }

    // here we integrate over the UL pdf
    else{

      // check that UL is there
      assert(hpdf_UL!=0);
      assert(hpdf2data_UL!=0);

      auto func = [hpdf_UL,cheb_x,cheb_y](double *x, double *p)->double{
	double x1=x[0];
	double x2=x[1];
	int m1 = int(p[0]);
	int m2 = int(p[1]);
	int d2 = int(p[2]);
	int p2 = int(p[3]);
	//int norm = int(p[4]);
	double pdf = hpdf_UL->GetBinContent( hpdf_UL->FindBin(x1,x2) );
	//if(norm>0){
	//cout << "x = " << x1 << ", y = " << x2 << ", m=(" << m1 << "," << m2 << ")"", val = " << val << endl; 
	//return pdf;
	///}
	cheb_x->SetParameter("m", m1);
	double totx = cheb_x->Eval(x1);
	cheb_y->SetParameter("m", m2);
	double cheby1 = cheb_y->Eval(x2);
	cheb_y->SetParameter("m", d2 - m2);
	double cheby2 = cheb_y->Eval(x2);
	double toty = 0.0;
	if(p2>0 && m2<(d2/2 + 1)){
	  toty = cheby1+cheby2;
	  if(d2%2==0 && m2==d2/2) toty *= 0.5;
	}
	else if(p2<0 && m2<(d2+1)/2) toty = (cheby1 - cheby2);
	double val = totx*toty*pdf;
	//cout << "x = " << x1 << ", y = " << x2 << ", m=(" << m1 << "," << m2 << "), totx = " << totx << ", toty = " << toty << ", pdf = " << pdf << ", val = " << val << endl; 
	return val;
      };
      TF2* funcAi = new TF2("func_"+iproc, func, X_edges[0], X_edges[X_nbins],
			    Y_edges[0], Y_edges[Y_nbins], 5 );
		
      // loop over data bins and fill outer product
      for(unsigned int idx = 0; idx<X_nbins; idx++){
	double xl = X_edges[idx];
	double xh = X_edges[idx+1];
	for(unsigned int idy = 0; idy<Y_nbins; idy++){
	  double yl = Y_edges[idy];
	  double yh = Y_edges[idy+1];
	  if(verbose) cout << "Doing bin number... " << count_d << endl;
	  y(count_d) = h->GetBinContent(idx+1, idy+1);
	  y_err(count_d) = h->GetBinError(idx+1, idy+1);
	  // loop over parameters
	  int count_p=0;
	  for(unsigned int ipx = 0; ipx<(deg_map[iproc].at(0) + 1); ipx++){
	    // constraint at 0
	    if( ctr_map[iproc].at(0) && ipx==0){
	      continue;
	    }	  	  
	    for(unsigned int ipy = 0; ipy<deg_map[iproc].at(1) + 1; ipy++){
	      // even/odd function:
	      if( (par_map[iproc].at(1)>0 && ipy < (deg_map[iproc].at(1)/2 + 1 - ctr_map[iproc].at(1))) || 
		  (par_map[iproc].at(1)<0 && ipy < (deg_map[iproc].at(1)+1)/2)){
		funcAi->SetParameter(0, ipx);
		funcAi->SetParameter(1, ipy);
		funcAi->SetParameter(2, deg_map[iproc].at(1));
		funcAi->SetParameter(3, par_map[iproc].at(1));
		//funcAi->SetParameter(4, 0);
		double val = funcAi->Integral(xl, xh, yl, yh, 1.e-6); 
		//if(verbose) cout << "val = " << val;
		//funcAi->SetParameter(4, 1);
		//double norm = funcAi->Integral(X_edges[idx], X_edges[idx+1], Y_edges[idy], Y_edges[idy+1], 1.e-6);		
		double norm = hpdf2data_UL->GetBinContent( idx+1, idy+1 );
		//if(verbose) cout << " norm = " << norm;
		J(count_d, count_p) = val/norm;
		if(verbose) cout << "jac(" << count_d << "," << count_p << ") = " << " = " << J(count_d, count_p) << endl;
		p_names.emplace_back( TString(Form("f_%d_%d", ipx, ipy)) ); 
		count_p++;
	      }
	    }
	  }
	  count_d++;
	}
      }      
    }

    MatrixXd V_inv_sqrt = MatrixXd::Zero(nd, nd);
    for(int i=0; i<nd; i++){
      V_inv_sqrt(i,i) = 1.0/y_err(i);
    }

    MatrixXd A = V_inv_sqrt*J;
    MatrixXd b = V_inv_sqrt*y;
    VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
    VectorXd delta = y-J*x;
    VectorXd pull  = b-A*x;
    if(verbose) cout << delta << endl;
    MatrixXd chi2 = pull.transpose()*pull;
    double chi2val = chi2(0,0);
    int ndof = y.size() - np;
    cout << "Chi2/ndof = " << chi2val << " / " << ndof << " = " << chi2val/ndof  << endl;
    MatrixXd W = (A.transpose()*A).inverse();
    
    TH1D* hinfo = new TH1D("h_info_"+iproc, "", 10, 0, 10);
    hinfo->SetBinContent(1, chi2val);
    hinfo->GetXaxis()->SetBinLabel(1, "chi2");
    hinfo->SetBinContent(2, ndof);
    hinfo->GetXaxis()->SetBinLabel(2, "ndof");
    hinfo->SetBinContent(3, TMath::Prob(chi2val,ndof));
    hinfo->GetXaxis()->SetBinLabel(3, "pvalue");
    hinfo->SetBinContent(4, npx);
    hinfo->GetXaxis()->SetBinLabel(4, "npx");
    hinfo->SetBinContent(5, npy);
    hinfo->GetXaxis()->SetBinLabel(5, "npy");
    hinfo->SetBinContent(6, deg_map[iproc].at(0));
    hinfo->GetXaxis()->SetBinLabel(6, "dx");
    hinfo->SetBinContent(7, deg_map[iproc].at(1));
    hinfo->GetXaxis()->SetBinLabel(7, "dy");
    hinfo->SetBinContent(8, nfpx);
    hinfo->GetXaxis()->SetBinLabel(8, "nfpx");
    hinfo->SetBinContent(9, nfpy);
    hinfo->GetXaxis()->SetBinLabel(9, "nfpy");
    hinfo->SetBinContent(10, shift_map[iproc]);
    hinfo->GetXaxis()->SetBinLabel(10, "shift");
    TH2D* hdelta = (TH2D*)hdummy->Clone("h_delta_"+iproc);
    TH2D* hpull  = (TH2D*)hdummy->Clone("h_pull_"+iproc);
    TH2D* hratio = (TH2D*)hdummy->Clone("h_ratio_"+iproc);
    TH2D* hdata  = (TH2D*)hdummy->Clone("h_data_"+iproc);
    TH2D* hexp   = (TH2D*)hdummy->Clone("h_exp_"+iproc);
    TH1D* hpar   = new TH1D("h_par_"+iproc, "", np, 0, np);
    TH2D* hcov   = new TH2D("h_cov_"+iproc, "", np, 0, np, np, 0, np);
    TH2D* hcor   = new TH2D("h_cor_"+iproc, "", np, 0, np, np, 0, np);

    //std::cout << "The matrix J.row is of size "
    //	      << J.row(0).rows() << "x" << J.row(0).cols() << std::endl;

    count_d=0;
    for(unsigned int idx = 0; idx<X_nbins; idx++){
      for(unsigned int idy = 0; idy<Y_nbins; idy++){	
	hdelta->SetBinContent(idx+1,idy+1, delta(count_d));
	hpull->SetBinContent(idx+1,idy+1,  pull(count_d));
	hratio->SetBinContent(idx+1,idy+1, y(count_d)!=0. ? (y(count_d)-delta(count_d))/y(count_d) : 1.0 );	
	hratio->SetBinError(idx+1,idy+1, y(count_d)!=0. ? y_err(count_d)/y(count_d) : 0.0 );	
	hdata->SetBinContent(idx+1,idy+1,  y(count_d));
	hdata->SetBinError(idx+1,idy+1,  y_err(count_d));
	hexp->SetBinContent(idx+1,idy+1,  y(count_d)-delta(count_d));
	double exp_err = (J.row(count_d)*W*J.row(count_d).transpose())(0,0);
	hexp->SetBinError(idx+1,idy+1, exp_err>0. ? TMath::Sqrt(exp_err) : 0.0);
	count_d++;
      }
    }
    
    for(unsigned int ip1 = 0; ip1<np; ip1++){      
      hpar->GetXaxis()->SetBinLabel(ip1+1, p_names[ip1]);
      hpar->SetBinContent(ip1+1, x(ip1));
      hcov->GetXaxis()->SetBinLabel(ip1+1, p_names[ip1]);
      hcor->GetXaxis()->SetBinLabel(ip1+1, p_names[ip1]);
      for(unsigned int ip2 = 0; ip2<np; ip2++){
	hcov->GetYaxis()->SetBinLabel(ip2+1, p_names[ip2]);
	hcor->GetYaxis()->SetBinLabel(ip2+1, p_names[ip2]);
	hcov->SetBinContent(ip1+1, ip2+1, W(ip1,ip2));
	hcor->SetBinContent(ip1+1, ip2+1, W(ip1,ip2)/TMath::Sqrt(W(ip1,ip1)*W(ip2,ip2)));
      }
    }


    cout << "Filling high-stat histo..." << endl;
    double X_edges_pdf[X_nbins*extrabinsX+1];
    for(unsigned int il=0; il<X_nbins; il++){
      for(unsigned int ih=0; ih<extrabinsX; ih++)
	X_edges_pdf[il*extrabinsX + ih] = X_edges[il] + double(ih)/extrabinsX*(X_edges[il+1] - X_edges[il]);
    }
    X_edges_pdf[X_nbins*extrabinsX] = X_edges[X_nbins];
    if(verbose){
      for(unsigned int ih=0; ih<=X_nbins*extrabinsX; ih++) cout << X_edges_pdf[ih] << ", ";
      cout << endl;
    }
    double Y_edges_pdf[Y_nbins*extrabinsY+1];
    for(unsigned int il=0; il<Y_nbins; il++){
      for(unsigned int ih=0; ih<extrabinsY; ih++)
	Y_edges_pdf[il*extrabinsY + ih] = Y_edges[il] + double(ih)/extrabinsY*(Y_edges[il+1] - Y_edges[il]);
    }
    Y_edges_pdf[Y_nbins*extrabinsY] = Y_edges[Y_nbins];
    if(verbose){
      for(unsigned int ih=0; ih<=Y_nbins*extrabinsY; ih++) cout << Y_edges_pdf[ih] << ", ";
      cout << endl;
    }
          
    TH2D* hpdf = new TH2D("h_pdf_"+iproc, "", X_nbins*extrabinsX, X_edges_pdf, Y_nbins*extrabinsY, Y_edges_pdf);
    for(unsigned int idx=1; idx<=hpdf->GetXaxis()->GetNbins(); idx++){
      double x_i = hpdf->GetXaxis()->GetBinCenter(idx); 
      for(unsigned int idy=1; idy<=hpdf->GetYaxis()->GetNbins(); idy++){
	double y_i = hpdf->GetYaxis()->GetBinCenter(idy); 

	// now computing the polynomial approximant
	double val = 0.0;
	int count_p=0;
	for(unsigned int ipx = 0; ipx<(deg_map[iproc].at(0) + 1); ipx++){
	  // constraint at 0
 	  if( ctr_map[iproc].at(0) && ipx==0){
	    continue;
	  }
	  cheb_x->SetParameter("m", ipx);
	  double totx = cheb_x->Eval(x_i); 
	  for(unsigned int ipy = 0; ipy<deg_map[iproc].at(1) + 1; ipy++){
	    cheb_y->SetParameter("m", ipy);
	    double cheby1 = cheb_y->Eval(y_i); 
	    cheb_y->SetParameter("m", deg_map[iproc].at(1) - ipy);
	    double cheby2 = cheb_y->Eval(y_i); 
	    // even function:
	    if( par_map[iproc].at(1)>0 && ipy < (deg_map[iproc].at(1)/2 + 1 - ctr_map[iproc].at(1))){
	      double toty = cheby1 + cheby2 ;
	      if( deg_map[iproc].at(1)%2==0 && ipy==deg_map[iproc].at(1)/2 ) toty *= 0.5;	      
	      val += x(count_p)*totx*toty;
	      count_p++;
	    }
	    // odd function
	    else if( par_map[iproc].at(1)<0 && ipy < (deg_map[iproc].at(1)+1)/2){
	      double toty = cheby1 - cheby2 ;
	      val += x(count_p)*totx*toty;
	      count_p++;
	    }	    
	  }
 	}
	if(interpolate) val = h_interpol_UL->Interpolate(x_i,y_i);
	hpdf->SetBinContent(idx,idy, val);
      }
    }        

    if(iproc=="UL")
      hpdf_UL = (TH2D*)hpdf->Clone("hpdf_UL");
    
    fout->mkdir(iproc);
    fout->cd(iproc+"/");
    
    cheb_x->SetParameter("n", degf_map[iproc].at(0) );
    cheb_y->SetParameter("n", degf_map[iproc].at(1) );

    if(xf_max>0.){
      cheb_x->SetParameter("scale",  xf_max*0.5 );
      cheb_x->SetParameter("offset", xf_max*0.5 );
    }
    else xf_max = X_max;
    if(yf_max>0.){
      cheb_y->SetParameter("scale",  yf_max );
      cheb_y->SetParameter("offset", 0.0 );
    } 
    else yf_max = Y_max; 

    TH2D* hpdf2data_p = (TH2D*)hdummy->Clone("h_pdf2data_"+iproc);        
    for(unsigned int idx = 0; idx<X_nbins; idx++){
      double xl = X_edges[idx];
      double xh = X_edges[idx+1];
      for(unsigned int idy = 0; idy<Y_nbins; idy++){
	double yl = Y_edges[idy];
	double yh = Y_edges[idy+1];
	double val = 0.0;
	for(unsigned int idfx = 0; idfx<hpdf->GetXaxis()->GetNbins(); idfx++){
	  if(!(hpdf->GetXaxis()->GetBinCenter(idfx+1)>=xl && hpdf->GetXaxis()->GetBinCenter(idfx+1)<=xh)) continue;
	  for(unsigned int idfy = 0; idfy<hpdf->GetYaxis()->GetNbins(); idfy++){
	    if(!(hpdf->GetYaxis()->GetBinCenter(idfy+1)>=yl && hpdf->GetYaxis()->GetBinCenter(idfy+1)<=yh)) continue;
	    val += hpdf->GetBinContent(idfx+1,idfy+1)*
	      (iproc=="UL" ? 1.0 : hpdf_UL->GetBinContent(idfx+1,idfy+1))*
	      hpdf->GetXaxis()->GetBinWidth(idfx+1)*
	      hpdf->GetYaxis()->GetBinWidth(idfy+1);
	  }
	}
	if(iproc!="UL"){
	  //cout << idx << "," << idy << " -> " << val << "/" << hpdf2data_UL->GetBinContent(idx+1,idy+1) << endl;
	  val /= hpdf2data_UL->GetBinContent(idx+1,idy+1);
	}
	hpdf2data_p->SetBinContent(idx+1,idy+1, val); 
      }
    }
    if(savePdf2data)
      hpdf2data_p->Write();
    if( iproc=="UL" ) hpdf2data_UL = (TH2D*)hpdf2data_p->Clone("hpdf2data_"+iproc);
    
    
    // now computing the polynomial approximant
    int count_pf = 0;
    for(unsigned int ipx = 0; ipx<(degf_map[iproc].at(0) + 1); ipx++){
      for(unsigned int ipy = 0; ipy<degf_map[iproc].at(1) + 1; ipy++){	  
	// it's a POI
	if( ipy < (degf_map[iproc].at(1)/2 + 1) ){
	  //cout << count_pf << ", " << ipx << ":" << ipy << endl;
	  for(int syst=0; syst<3; syst++){
	    TString syst_name = "";
	    if(syst==0)      syst_name = "";
	    else if(syst==1) syst_name = "Up";
	    else if(syst==2) syst_name = "Down";

	    // save pdf of varied f_syst as fine-grained TH2D
	    TH2D* hpdf_p = (TH2D*)hpdf->Clone( syst==0 ?
					       "h_pdf_"+iproc+Form("_jac%d",count_pf) :
					       "h_pdf_"+iproc+Form("_syst%d",count_pf)+syst_name);
	    for(unsigned int idx=1; idx<=hpdf_p->GetXaxis()->GetNbins(); idx++){
	      double x_i = hpdf_p->GetXaxis()->GetBinCenter(idx); 
	      cheb_x->SetParameter("m", ipx);
	      double totx = cheb_x->Eval(x_i); 
	      for(unsigned int idy=1; idy<=hpdf_p->GetYaxis()->GetNbins(); idy++){
		// restrict to fiducial phase-space
		if( (hpdf_p->GetXaxis()->GetBinLowEdge(idx) >= xf_max) || (hpdf_p->GetYaxis()->GetBinLowEdge(idy) >= yf_max) ){
		  hpdf_p->SetBinContent(idx,idy, syst==0 ? 0.0 : hpdf->GetBinContent(idx,idy));
		  continue;
		}
		double y_i = hpdf_p->GetYaxis()->GetBinCenter(idy); 
		cheb_y->SetParameter("m", ipy);
		double cheby1 = cheb_y->Eval(y_i); 
		cheb_y->SetParameter("m", degf_map[iproc].at(1) - ipy);
		double cheby2 = cheb_y->Eval(y_i);
		double toty = cheby1 + cheby2 ;
		if( degf_map[iproc].at(1)%2==0 && ipy==degf_map[iproc].at(1)/2 ) toty *= 0.5;	      
		double val = 0.0;
		if(syst==0)
		  val = hpdf->GetBinContent(idx,idy)*totx*toty ;
		else{
		  double shift = syst==1 ? +shift_map[iproc] : -shift_map[iproc];
		  val = hpdf->GetBinContent(idx,idy)*(1 + shift*totx*toty) ;
		}
		hpdf_p->SetBinContent(idx,idy, val);
	      }
	    }
	    if((savePdf && (syst==0 && saveJac)) || (savePdf && (syst>0 && saveSyst)) )
	      hpdf_p->Write();

	    // save data expected of varifed f_syst as binned density
	    TH2D* hpdf2data_p_j = (TH2D*)hdummy->Clone( syst==0 ?
							"h_pdf2data_"+iproc+Form("_jac%d",count_pf) :
							"h_pdf2data_"+iproc+Form("_syst%d",count_pf)+syst_name
							);        
	    for(unsigned int idx = 0; idx<X_nbins; idx++){
	      double xl = X_edges[idx];
	      double xh = X_edges[idx+1];
	      for(unsigned int idy = 0; idy<Y_nbins; idy++){
		double yl = Y_edges[idy];
		double yh = Y_edges[idy+1];
		double val = 0.0;
		for(unsigned int idfx = 0; idfx<hpdf_p->GetXaxis()->GetNbins(); idfx++){
		  if(!(hpdf_p->GetXaxis()->GetBinCenter(idfx+1)>=xl && hpdf_p->GetXaxis()->GetBinCenter(idfx+1)<=xh)) continue;
		  for(unsigned int idfy = 0; idfy<hpdf_p->GetYaxis()->GetNbins(); idfy++){
		    if(!(hpdf_p->GetYaxis()->GetBinCenter(idfy+1)>=yl && hpdf_p->GetYaxis()->GetBinCenter(idfy+1)<=yh)) continue;
		    val += hpdf_p->GetBinContent(idfx+1,idfy+1)*
		      (iproc=="UL" ? 1.0 : hpdf_UL->GetBinContent(idfx+1,idfy+1))*
		      hpdf_p->GetXaxis()->GetBinWidth(idfx+1)*
		      hpdf_p->GetYaxis()->GetBinWidth(idfy+1);
		  }
		}
		if(iproc!="UL"){
		  val /= hpdf2data_UL->GetBinContent(idx+1,idy+1);
		}
		hpdf2data_p_j->SetBinContent(idx+1,idy+1, val); 
	      }
	    }
	    if((savePdf2data && (syst==0 && saveJac)) || (savePdf2data && (syst>0 && saveSyst)) )
	      hpdf2data_p_j->Write();
	  }	  
	  count_pf++;
	}	  
      }	
    }   

    hinfo->Write();
    hpull->Write();
    hratio->Write();
    hdelta->Write();
    hdata->Write();
    hexp->Write();
    hcov->Write();
    hcor->Write();
    hpar->Write();
    hpdf->Write();    
  }
  
  fout->Close();
  fin->Close();


  TRandom3* ran = new TRandom3();
  delete ran;
  sw.Stop();
  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;
  return 1;
}
