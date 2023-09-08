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
	("dUL_x",   value<int>()->default_value(2), "max degree in x of corrxy")
	("dUL_y",   value<int>()->default_value(2), "max degree in y of corrxy")
	("dA0_x",   value<int>()->default_value(2), "max degree in x for A0")
	("dA0_y",   value<int>()->default_value(2), "max degree in y for A0")
	("dA1_x",   value<int>()->default_value(2), "max degree in x for A1")
	("dA1_y",   value<int>()->default_value(2), "max degree in y for A1")
	("dA2_x",   value<int>()->default_value(2), "max degree in x for A2")
	("dA2_y",   value<int>()->default_value(2), "max degree in y for A2")
	("dA3_x",   value<int>()->default_value(2), "max degree in x for A3")
	("dA3_y",   value<int>()->default_value(2), "max degree in y for A3")
	("dA4_x",   value<int>()->default_value(2), "max degree in x for A4")
	("dA4_y",   value<int>()->default_value(2), "max degree in y for A4")
	("fUL_x",   value<int>()->default_value(2), "max degree of modifier in x of corrxy")
	("fUL_y",   value<int>()->default_value(2), "max degree of modifier in y of corrxy")
	("fA0_x",   value<int>()->default_value(2), "max degree of modifier in x for A0")
	("fA0_y",   value<int>()->default_value(2), "max degree of modifier in y for A0")
	("fA1_x",   value<int>()->default_value(2), "max degree of modifier in x for A1")
	("fA1_y",   value<int>()->default_value(2), "max degree of modifier in y for A1")
	("fA2_x",   value<int>()->default_value(2), "max degree of modifier in x for A2")
	("fA2_y",   value<int>()->default_value(2), "max degree of modifier in y for A2")
	("fA3_x",   value<int>()->default_value(2), "max degree of modifier in x for A3")
	("fA3_y",   value<int>()->default_value(2), "max degree of modifier in y for A3")
	("fA4_x",   value<int>()->default_value(2), "max degree of modifier in x for A4")
	("fA4_y",   value<int>()->default_value(2), "max degree of modifier in y for A4")
	("x_max",   value<double>()->default_value(-1.0), "max x value for syst")
	("y_max",   value<double>()->default_value(-1.0), "max y value for syst")
	("tag",     value<std::string>()->default_value("default"), "tag name")
      	("doA0",    bool_switch()->default_value(false), "")
	("doA1",    bool_switch()->default_value(false), "")
	("doA2",    bool_switch()->default_value(false), "")
	("doA3",    bool_switch()->default_value(false), "")
	("doA4",    bool_switch()->default_value(false), "")	
      	("verbose",   bool_switch()->default_value(false), "")
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
  std::string tag = vm["tag"].as<std::string>();
  int verbose     = vm["verbose"].as<bool>();
  int debug       = vm["debug"].as<bool>();
  int doA0        = vm["doA0"].as<bool>();
  int doA1        = vm["doA1"].as<bool>();
  int doA2        = vm["doA2"].as<bool>();
  int doA3        = vm["doA3"].as<bool>();
  int doA4        = vm["doA4"].as<bool>();

  int dUL_x = vm["dUL_x"].as<int>();
  int dUL_y = vm["dUL_y"].as<int>();
  int dA0_x   = vm["dA0_x"].as<int>();
  int dA0_y   = vm["dA0_y"].as<int>();
  int dA1_x   = vm["dA1_x"].as<int>();
  int dA1_y   = vm["dA1_y"].as<int>();
  int dA2_x   = vm["dA2_x"].as<int>();
  int dA2_y   = vm["dA2_y"].as<int>();
  int dA3_x   = vm["dA3_x"].as<int>();
  int dA3_y   = vm["dA3_y"].as<int>();
  int dA4_x   = vm["dA4_x"].as<int>();
  int dA4_y   = vm["dA4_y"].as<int>();

  int fUL_x   = vm["fUL_x"].as<int>();
  int fUL_y   = vm["fUL_y"].as<int>();
  int fA0_x   = vm["fA0_x"].as<int>();
  int fA0_y   = vm["fA0_y"].as<int>();
  int fA1_x   = vm["fA1_x"].as<int>();
  int fA1_y   = vm["fA1_y"].as<int>();
  int fA2_x   = vm["fA2_x"].as<int>();
  int fA2_y   = vm["fA2_y"].as<int>();
  int fA3_x   = vm["fA3_x"].as<int>();
  int fA3_y   = vm["fA3_y"].as<int>();
  int fA4_x   = vm["fA4_x"].as<int>();
  int fA4_y   = vm["fA4_y"].as<int>();

  double x_max   = vm["x_max"].as<double>();
  double y_max   = vm["y_max"].as<double>();

  TFile *fout = TFile::Open("fout.root", "RECREATE");

  std::vector<TString> proc = {"UL"};
  if(doA0) proc.emplace_back("A0");
  if(doA1) proc.emplace_back("A1");
  if(doA2) proc.emplace_back("A2");
  if(doA3) proc.emplace_back("A3");
  if(doA4) proc.emplace_back("A4");


  // dummy file
  if(debug){
    TFile* f = TFile::Open( "root/ai_2dmap_qtbyQ_and_qt_vs_absy.root", "RECREATE");  
    for(auto& pr : proc){
      TH2D* h = new TH2D("ang_coeff_wp_qtbyQ_vs_absy_A_"+TString(pr[1]), "", 20, 0.0, 0.5, 20, 0.0, 2.5 );
      for(int ibx=1; ibx<=h->GetXaxis()->GetNbins() ; ibx++ ){
	for(int iby=1; iby<=h->GetYaxis()->GetNbins() ; iby++ ){
	  if(pr=="A1" || pr=="A3"){
	    h->SetBinContent(ibx,iby, 10000.*iby/h->GetYaxis()->GetNbins()*ibx/h->GetXaxis()->GetNbins());
	  }
	  else if(pr=="UL" || pr=="A0" || pr=="A2"){
	    h->SetBinContent(ibx,iby, 10000.*ibx/h->GetXaxis()->GetNbins());
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
  
  TFile* fin = TFile::Open("root/ai_2dmap_qtbyQ_and_qt_vs_absy.root", "READ");
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
  deg_map.insert  ( std::make_pair<TString, std::array<int,2> >("UL", {dUL_x,  dUL_y}) );
  par_map.insert  ( std::make_pair<TString, std::array<int,2> >("UL", {0, +1}) );
  ctr_map.insert  ( std::make_pair<TString, std::array<int,2> >("UL", {1,  0}) ); 
  //A0
  deg_map.insert  ( std::make_pair<TString, std::array<int,2> >("A0", {dA0_x,  dA0_y}) );
  par_map.insert  ( std::make_pair<TString, std::array<int,2> >("A0", {0, +1}) );
  ctr_map.insert  ( std::make_pair<TString, std::array<int,2> >("A0", {1,  0}) ); 
  //A1
  deg_map.insert  ( std::make_pair<TString, std::array<int,2> >("A1", {dA1_x,  dA1_y}) );
  par_map.insert  ( std::make_pair<TString, std::array<int,2> >("A1", {0, -1}) );
  ctr_map.insert  ( std::make_pair<TString, std::array<int,2> >("A1", {1,  0}) ); 
  //A2
  deg_map.insert  ( std::make_pair<TString, std::array<int,2> >("A2", {dA2_x,  dA2_y}) );
  par_map.insert  ( std::make_pair<TString, std::array<int,2> >("A2", {0, +1}) );
  ctr_map.insert  ( std::make_pair<TString, std::array<int,2> >("A2", {1,  0}) ); 
  //A3
  deg_map.insert  ( std::make_pair<TString, std::array<int,2> >("A3", {dA3_x,  dA3_y}) );
  par_map.insert  ( std::make_pair<TString, std::array<int,2> >("A3", {0, +1}) );
  ctr_map.insert  ( std::make_pair<TString, std::array<int,2> >("A3", {1,  1}) ); 
  //A4
  deg_map.insert  ( std::make_pair<TString, std::array<int,2> >("A4", {dA4_x,  dA4_y}) );
  par_map.insert  ( std::make_pair<TString, std::array<int,2> >("A4", {0, -1}) );
  ctr_map.insert  ( std::make_pair<TString, std::array<int,2> >("A4", {0,  0}) );

  degf_map.insert  ( std::make_pair<TString, std::array<int,2> >("UL", {fUL_x,  fUL_y}) );
  degf_map.insert  ( std::make_pair<TString, std::array<int,2> >("A0", {fA0_x,  fA0_y}) );
  degf_map.insert  ( std::make_pair<TString, std::array<int,2> >("A1", {fA1_x,  fA1_y}) );
  degf_map.insert  ( std::make_pair<TString, std::array<int,2> >("A2", {fA2_x,  fA2_y}) );
  degf_map.insert  ( std::make_pair<TString, std::array<int,2> >("A3", {fA3_x,  fA3_y}) );
  degf_map.insert  ( std::make_pair<TString, std::array<int,2> >("A4", {fA4_x,  fA4_y}) );

  TH2D* hdatafine_UL = 0;
  
  // loop over all histos
  for(unsigned int i=0; i <proc.size(); i++){

    TString iproc = proc[i];
    cout << "Doing proc " << iproc << endl;
    TH2D* h = (TH2D*)fin->Get("ang_coeff_wp_qtbyQ_vs_absy_A_"+TString(iproc[1]));    
    if(h==0){
      cout << "Histo not found. Continue." << endl;
      continue;
    }
    cout << "Histo " << h->GetName() << " found" << endl;
    if(iproc=="UL") h->Scale(1./ h->Integral());

    // X-axis
    int X_nbins  = h->GetXaxis()->GetNbins();
    double X_min = h->GetXaxis()->GetXmin();
    double X_max = h->GetXaxis()->GetXmax();
    double X_scale  = (X_max-X_min)*0.5;
    double X_offset = (X_max+X_min)*0.5;
    double X_edges[X_nbins+1];
    //for(int ib = 0; ib<X_nbins; ib++) X_edges[ib] = (h->GetXaxis()->GetBinLowEdge(ib+1)+X_min)/X_range;
    //X_edges[X_nbins] = (h->GetXaxis()->GetBinUpEdge(X_nbins)+X_min)/X_range;
    for(int ib = 0; ib<X_nbins; ib++) X_edges[ib] = h->GetXaxis()->GetBinLowEdge(ib+1);
    X_edges[X_nbins] = h->GetXaxis()->GetBinUpEdge(X_nbins);
    if(verbose) cout << "X_edges: ";
    if(verbose) for(int i = 0; i <= X_nbins; i++) cout << X_edges[i] << "," ;
    if(verbose) cout << endl;
    
    // Y-axis
    int Y_nbins  = h->GetYaxis()->GetNbins();
    double Y_max = h->GetYaxis()->GetXmax();
    double Y_min = -Y_max; // assume |y| is plotted
    double Y_scale  = (Y_max-Y_min)*0.5;
    double Y_offset = (Y_max+Y_min)*0.5;
    double Y_edges[Y_nbins+1];
    //for(int ib = 0; ib<Y_nbins; ib++) Y_edges[ib] = (h->GetYaxis()->GetBinLowEdge(ib+1)+Y_min)/Y_range;
    //Y_edges[Y_nbins] = (h->GetYaxis()->GetBinUpEdge(Y_nbins)+Y_min)/Y_range;
    for(int ib = 0; ib<Y_nbins; ib++) Y_edges[ib] = h->GetYaxis()->GetBinLowEdge(ib+1);
    Y_edges[Y_nbins] = h->GetYaxis()->GetBinUpEdge(Y_nbins) ;
    if(verbose) cout << "Y_edges: ";
    if(verbose) for(int i = 0; i <= Y_nbins; i++) cout << Y_edges[i] << "," ;
    if(verbose) cout << endl;
    
    assert( ctr_map[iproc].at(0) || (ctr_map[iproc].at(1) && deg_map[iproc].at(1)%2==0));

    int nd = X_nbins*Y_nbins;
    int npx = deg_map[iproc].at(0) + 1 - ctr_map[iproc].at(0);
    int npy = par_map[iproc].at(1)>0 ? (deg_map[iproc].at(1)/2 + 1) : (deg_map[iproc].at(1) + 1)/2;
    if( ctr_map[iproc].at(1) ) npy--;
    int np = npx*npy;
    if(verbose) cout << "np = " << npx << " * " << npy << " = " << np << endl;

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
      assert(hdatafine_UL!=0);

      auto func = [hdatafine_UL,cheb_x,cheb_y](double *x, double *p)->double{
	double x1=x[0];
	double x2=x[1];
	int m1 = int(p[0]);
	int m2 = int(p[1]);
	int d2 = int(p[2]);
	int p2 = int(p[3]);
	int norm = int(p[4]);
	double pdf = hdatafine_UL->GetBinContent( hdatafine_UL->FindBin(x1,x2) );
	if(norm>0) return pdf;
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
	//cout << "x = " << x1 << ", y = " << x2 << ", m=(" << m1 << "," << m2 << "), totx = " << totx << ", toty = " << toty << ", val = " << val << endl; 
	return val;
      };
      TF2* funcAi = new TF2("func_"+iproc, func, X_edges[0], X_edges[X_nbins],
			    Y_edges[0], Y_edges[Y_nbins], 5 );
		
      // loop over data bins and fill outer product
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
	      // even/odd function:
	      if( (par_map[iproc].at(1)>0 && ipy < (deg_map[iproc].at(1)/2 + 1 - ctr_map[iproc].at(1))) || 
		  (par_map[iproc].at(1)<0 && ipy < (deg_map[iproc].at(1)+1)/2)){
		funcAi->SetParameter(0, ipx);
		funcAi->SetParameter(1, ipy);
		funcAi->SetParameter(2, deg_map[iproc].at(1));
		funcAi->SetParameter(3, par_map[iproc].at(1));
		funcAi->SetParameter(4, 0);
		double val = funcAi->Integral(X_edges[idx], X_edges[idx+1], Y_edges[idy], Y_edges[idy+1], 1.e-6); 
		//if(verbose) cout << "val = " << val;
		funcAi->SetParameter(4, 1);
		double norm = funcAi->Integral(X_edges[idx], X_edges[idx+1], Y_edges[idy], Y_edges[idy+1], 1.e-6);		
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
    
    TH1D* hinfo = new TH1D("h_info_"+iproc, "", 4, 0,4);
    hinfo->SetBinContent(1, chi2val);
    hinfo->GetXaxis()->SetBinLabel(1, "chi2");
    hinfo->SetBinContent(2, ndof);
    hinfo->GetXaxis()->SetBinLabel(2, "ndof");
    hinfo->SetBinContent(3, npx);
    hinfo->GetXaxis()->SetBinLabel(3, "npx");
    hinfo->SetBinContent(4, npy);
    hinfo->GetXaxis()->SetBinLabel(4, "npy");
    TH2D* hdelta = (TH2D*)h->Clone("h_delta_"+iproc);
    TH2D* hpull  = (TH2D*)h->Clone("h_pull_"+iproc);
    TH2D* hdata  = (TH2D*)h->Clone("h_data_"+iproc);
    TH1D* hpar   = new TH1D("h_par_"+iproc, "", np, 0, np);
    TH2D* hcov   = new TH2D("h_cov_"+iproc, "", np, 0, np, np, 0, np);
    TH2D* hcor   = new TH2D("h_cor_"+iproc, "", np, 0, np, np, 0, np);

    count_d=0;
    for(unsigned int idx = 0; idx<X_nbins; idx++){
      for(unsigned int idy = 0; idy<Y_nbins; idy++){	
	hdelta->SetBinContent(idx+1,idy+1, delta(count_d));
	hpull->SetBinContent(idx+1,idy+1,  pull(count_d));
	hdata->SetBinContent(idx+1,idy+1,  y(count_d));
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
    TH2D* hdatafine = new TH2D("h_pdffine_"+iproc, "", 100, X_edges[0], X_edges[X_nbins], 100, Y_edges[0], Y_edges[Y_nbins]);
    for(unsigned int idx=1; idx<=hdatafine->GetXaxis()->GetNbins(); idx++){
      double x_i = hdatafine->GetXaxis()->GetBinCenter(idx); 
      for(unsigned int idy=1; idy<=hdatafine->GetYaxis()->GetNbins(); idy++){
	double y_i = hdatafine->GetYaxis()->GetBinCenter(idy); 

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
	hdatafine->SetBinContent(idx,idy, val);
      }
    }        

    if(iproc=="UL")
      hdatafine_UL = (TH2D*)hdatafine->Clone("hdatafine_UL");
    
    fout->mkdir(iproc);
    fout->cd(iproc+"/");
    
    cheb_x->SetParameter("n", degf_map[iproc].at(0) );
    cheb_y->SetParameter("n", degf_map[iproc].at(1) );

    if(x_max>0.){
      cheb_x->SetParameter("scale",  x_max*0.5 );
      cheb_x->SetParameter("offset", x_max*0.5 );
    }
    else x_max = X_max;
    if(y_max>0.){
      cheb_y->SetParameter("scale",  y_max );
      cheb_y->SetParameter("offset", 0.0 );
    } 
    else y_max = Y_max; 

    auto funcUL_p = [hdatafine_UL](double* x, double* p)->double{
      double val = 0.0;
      double x1 = x[0];
      double x2 = x[1];
      val = hdatafine_UL->GetBinContent( hdatafine_UL->FindBin(x1,x2) );      
      return val;      
    };
    TF2* funcUL_data_p = new TF2("funcUL_data_"+iproc, funcUL_p, X_edges[0], X_edges[X_nbins], Y_edges[0], Y_edges[Y_nbins], 0 );

    // now computing the polynomial approximant
    int count_pf = 0;
    for(unsigned int ipx = 0; ipx<(degf_map[iproc].at(0) + 1); ipx++){
      for(unsigned int ipy = 0; ipy<degf_map[iproc].at(1) + 1; ipy++){	  
	// it's a POI
	if( ipy < (deg_map[iproc].at(1)/2 + 1)){
	  for(int syst=0; syst<2; syst++){
	    TString syst_name = syst==0 ? "up" : "down";
	    TH2D* hdatafine_p = (TH2D*)hdatafine->Clone("h_pdffine_"+iproc+Form("_jac%d_",count_pf)+syst_name);
	    for(unsigned int idx=1; idx<=hdatafine_p->GetXaxis()->GetNbins(); idx++){
	      // restrict to fiducial phase-space
	      if( hdatafine_p->GetXaxis()->GetBinLowEdge(idx) > x_max ) continue;
	      double x_i = hdatafine_p->GetXaxis()->GetBinCenter(idx); 
	      cheb_x->SetParameter("m", ipx);
	      double totx = cheb_x->Eval(x_i); 
	      for(unsigned int idy=1; idy<=hdatafine_p->GetYaxis()->GetNbins(); idy++){
		// restrict to fiducial phase-space
		if( hdatafine_p->GetYaxis()->GetBinLowEdge(idy) > y_max ) continue;
		double y_i = hdatafine_p->GetYaxis()->GetBinCenter(idy); 
		cheb_y->SetParameter("m", ipy);
		double cheby1 = cheb_y->Eval(y_i); 
		cheb_y->SetParameter("m", degf_map[iproc].at(1) - ipy);
		double cheby2 = cheb_y->Eval(y_i);
		double toty = cheby1 + cheby2 ;
		if( degf_map[iproc].at(1)%2==0 && ipy==degf_map[iproc].at(1)/2 ) toty *= 0.5;	      
		double shift = syst==0 ? +0.1 : -0.1;
		double val = hdatafine->GetBinContent(idx,idy)*(1 + shift*totx*toty) ;
		hdatafine_p->SetBinContent(idx,idy, val);
	      }
	    }
	    hdatafine_p->Write();
	    
	    auto funcAi_p = [hdatafine_p,hdatafine_UL](double* x, double* p)->double{
	      double val = 0.0;
	      double x1 = x[0];
	      double x2 = x[1];
	      //cout << x1 << "," << x2 << " => " << hdatafine_p->GetBinContent( hdatafine_p->FindBin(x1,x2) ) << " * " << hdatafine_UL->GetBinContent( hdatafine_UL->FindBin(x1,x2) ) << endl;
	      val = hdatafine_p->GetBinContent( hdatafine_p->FindBin(x1,x2) )*hdatafine_UL->GetBinContent( hdatafine_UL->FindBin(x1,x2) );
	      return val;      
	    };
	    TF2* funcAi_data_p = new TF2("funcAi_data_"+iproc+Form("_jac%d_",count_pf)+syst_name, funcAi_p, X_edges[0], X_edges[X_nbins], Y_edges[0], Y_edges[Y_nbins], 0 );
	    
	    TH2D* hdata_p = (TH2D*)h->Clone("h_data_"+iproc+Form("_jac%d_",count_pf)+syst_name);
	    for(unsigned int idx = 0; idx<X_nbins; idx++){
	      for(unsigned int idy = 0; idy<Y_nbins; idy++){
		double val = 0.0;
		double den = funcUL_data_p->Integral(X_edges[idx], X_edges[idx+1], Y_edges[idy], Y_edges[idy+1], 1.e-6);
		if(iproc=="UL")
		  val = den;
		else
		  val = den>0. ? funcAi_data_p->Integral(X_edges[idx], X_edges[idx+1], Y_edges[idy], Y_edges[idy+1], 1.e-6)/den : 0.0;
		hdata_p->SetBinContent(idx+1,idy+1, val);
	      }
	    }
	    hdata_p->Write();	    	    	    
	  }
	  count_pf++;
	}	  
      }	
    }   

    hinfo->Write();
    hpull->Write();
    hdelta->Write();
    hdata->Write();
    hcov->Write();
    hcor->Write();
    hpar->Write();
    hdatafine->Write();    
  }
  
  fout->Close();
  fin->Close();


  TRandom3* ran = new TRandom3();
  delete ran;
  sw.Stop();
  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;
  return 1;
}
