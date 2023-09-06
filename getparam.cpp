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
    double xj = (TMath::Cos((par[0]-i)*TMath::Pi()/par[0]) + par[1])*par[2]; // par[1]=offset, par[2]=scale
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
	("nevents",     value<long>()->default_value(1000), "number of events")
	("ntoys",       value<int>()->default_value(-1), "number of toys")
	("dUL_x", value<int>()->default_value(2), "max degree in x of corrxy")
	("dUL_y", value<int>()->default_value(2), "max degree in y of corrxy")
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
	("tag",     value<std::string>()->default_value("default"), "tag name")
	("debug",       bool_switch()->default_value(false), "");

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
  int debug       = vm["debug"].as<bool>();

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

  TFile *fout = TFile::Open("fout.root", "RECREATE");

  std::vector<TString> proc = {"UL",
			       "A0",
			       "A1",
			       "A2",
			       "A3",
			       "A4"
  };

  // dummy file
  TFile* f = TFile::Open( "root/histos_default.root", "RECREATE");  
  for(auto pr : proc){
    TH2D* h = new TH2D("histo_"+pr, "", 10, 0.0, 0.5, 10, 0.0, 2.5 );
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
  
  TFile* fin = TFile::Open(("root/histos_"+tag+".root").c_str(), "READ");
  if(fin==0 || fin==nullptr || fin->IsZombie()){
    cout << "File NOT found" << endl;
    return 0;
  }


  fin->cd();
  
  // [0, 1]
  TF1* cheb_x = new TF1("cheb_x", cheb_fct, 0.0, 1.0, 4); 
  cheb_x->SetParNames("n","offset","scale","m");
  cheb_x->SetParameter("offset", 1.0);
  cheb_x->SetParameter("scale",  0.5);

  // [-1,+1]
  TF1* cheb_y = new TF1("cheb_y", cheb_fct, 0.0, 1.0, 4); 
  cheb_y->SetParNames("n","offset","scale","m");
  cheb_y->SetParameter("offset", 0.0);
  cheb_y->SetParameter("scale",  1.0);

  std::map<TString, std::array<int,2> > deg_map;
  std::map<TString, std::array<int,2> > par_map;
  std::map<TString, std::array<int,2> > ctr_map;

  //UL
  deg_map.insert  ( std::make_pair<TString, std::array<int,2> >("UL", {dUL_x,  dUL_y}) );
  par_map.insert  ( std::make_pair<TString, std::array<int,2> >("UL", {0, +1}) );
  ctr_map.insert  ( std::make_pair<TString, std::array<int,2> >("UL", {0,  0}) ); 
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
  ctr_map.insert  ( std::make_pair<TString, std::array<int,2> >("A3", {1,  0}) ); 
  //A4
  deg_map.insert  ( std::make_pair<TString, std::array<int,2> >("A4", {dA4_x,  dA4_y}) );
  par_map.insert  ( std::make_pair<TString, std::array<int,2> >("A4", {0, -1}) );
  ctr_map.insert  ( std::make_pair<TString, std::array<int,2> >("A4", {0,  0}) ); 
    
  // loop over all histos
  for(unsigned int i=0; i <proc.size(); i++){

    TString iproc = proc[i];
    cout << "Doing proc " << iproc << endl;
    TH2D* h = (TH2D*)fin->Get("histo_"+iproc);    
    cout << "Histo " << h->GetName() << "found" << endl;
    
    // X-axis
    int X_nbins  = h->GetXaxis()->GetNbins();
    double X_min = h->GetXaxis()->GetXmin();
    double X_max = h->GetXaxis()->GetXmax();
    double X_range = X_max-X_min;
    double X_edges[X_nbins+1];
    for(int ib = 0; ib<X_nbins; ib++) X_edges[ib] = (h->GetXaxis()->GetBinLowEdge(ib+1)+X_min)/X_range;
    X_edges[X_nbins] = (h->GetXaxis()->GetBinUpEdge(X_nbins)+X_min)/X_range;
    cout << "X_edges: ";
    for(int i = 0; i <= X_nbins; i++) cout << X_edges[i] << "," ;
    cout << endl;
    
    // Y-axis
    int Y_nbins  = h->GetYaxis()->GetNbins();
    double Y_min = h->GetYaxis()->GetXmin();
    double Y_max = h->GetYaxis()->GetXmax();
    double Y_range = Y_max-Y_min;
    double Y_edges[Y_nbins+1];
    for(int ib = 0; ib<Y_nbins; ib++) Y_edges[ib] = (h->GetYaxis()->GetBinLowEdge(ib+1)+Y_min)/Y_range;
    Y_edges[Y_nbins] = (h->GetYaxis()->GetBinUpEdge(Y_nbins)+Y_min)/Y_range;
    cout << "Y_edges: ";
    for(int i = 0; i <= Y_nbins; i++) cout << Y_edges[i] << "," ;
    cout << endl;

    int nd = X_nbins*Y_nbins;
    int npx = deg_map[iproc].at(0) + 1 - ctr_map[iproc].at(0);
    int npy = par_map[iproc].at(1)>0 ? (deg_map[iproc].at(1)/2 + 1) : (deg_map[iproc].at(1) + 1)/2;
    int np = npx*npy;
    cout << "np = " << npx << " * " << npy << " = " << np << endl;
    
    cheb_x->SetParameter("n", deg_map[iproc].at(0) );
    cheb_y->SetParameter("n", deg_map[iproc].at(1) );
    
    int count_d=0;
    VectorXd y(nd);
    VectorXd y_err(nd);
    MatrixXd J(nd,np);

    std::vector<TString> p_names = {};
      
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
	cout << "\tmX=" << ipx << ": int_x [" << X_edges[idx] << "," << X_edges[idx+1] << "] = " << totx << endl;
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
	if( par_map[iproc].at(1)>0 && ipy < (deg_map[iproc].at(1)/2 + 1)){
	  double toty = cheby1+cheby2 ;
	  cout << "\tmY=[" << ipy << "," << deg_map[iproc].at(1) - ipy << "] : int_y [" << Y_edges[idy] << "," << Y_edges[idy+1] << "] = (" << cheby1  << " + " << cheby2 << ")";
	  if( deg_map[iproc].at(1)%2==0 && ipy==deg_map[iproc].at(1)/2 ){
	    toty *= 0.5;
	    cout << "/2" << endl;
	  }
	  else cout << endl;
	  integrals_y(idy,ipy) = toty;
	}
	// odd function
	else if( par_map[iproc].at(1)<0 && ipy < (deg_map[iproc].at(1)+1)/2){
	  double toty = cheby1 - cheby2 ;
	  cout << "\tmY=[" << ipy << "," << deg_map[iproc].at(1) - ipy << "] : int_y [" << Y_edges[idy] << "," << Y_edges[idy+1] << "] = (" << cheby1  << " - " << cheby2 << ")";
	  integrals_y(idy,ipy) = toty;
	}
      }	
    }

    // loop over data bins and fill outer product
    cout << "Start the final loop..." << endl;
    for(unsigned int idx = 0; idx<X_nbins; idx++){
      for(unsigned int idy = 0; idy<Y_nbins; idy++){
	cout << "Doing bin number... " << count_d << endl;
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
	    if( (par_map[iproc].at(1)>0 && ipy < deg_map[iproc].at(1)/2 + 1) || (par_map[iproc].at(1)<0 && ipy < (deg_map[iproc].at(1)+1)/2)){
	      J(count_d, count_p) = integrals_x(idx,ipx)*integrals_y(idy,ipy);
	      cout << "jac(" << count_d << "," << count_p << ") = " << integrals_x(idx,ipx) << " * " << integrals_y(idy,ipy) << " = " << J(count_d, count_p) << endl;
	      p_names.emplace_back( TString(Form("f_%d_%d", ipx, ipy)) ); 
	      count_p++;
	    }
	  }
	}
	count_d++;
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
    cout << delta << endl;
    MatrixXd chi2 = pull.transpose()*pull;
    double chi2val = chi2(0,0);
    int ndof = y.size() - np;
    cout << "Chi2 = " << chi2val << endl;
    MatrixXd W = (A.transpose()*A).inverse();
    
    count_d=0;
    TH1D* hinfo = new TH1D("hinfo_"+iproc, "", 4, 0,4);
    hinfo->SetBinContent(1, chi2val);
    hinfo->GetXaxis()->SetBinLabel(1, "chi2");
    hinfo->SetBinContent(2, ndof);
    hinfo->GetXaxis()->SetBinLabel(2, "ndof");
    hinfo->SetBinContent(3, npx);
    hinfo->GetXaxis()->SetBinLabel(3, "npx");
    hinfo->SetBinContent(4, npy);
    hinfo->GetXaxis()->SetBinLabel(4, "npy");
    TH2D* hdelta = (TH2D*)h->Clone("hdelta_"+iproc);
    TH2D* hpull  = (TH2D*)h->Clone("hpull_"+iproc);
    TH2D* hdata  = (TH2D*)h->Clone("hdata_"+iproc);
    TH1D* hpar   = new TH1D("hpar_"+iproc, "", np, 0, np);
    TH2D* hcov   = new TH2D("hcov_"+iproc, "", np, 0, np, np, 0, np);
    TH2D* hcor   = new TH2D("hcor_"+iproc, "", np, 0, np, np, 0, np);
    cout << "Start the loop..." << endl;
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
    fout->cd();
    hinfo->Write();
    hpull->Write();
    hdelta->Write();
    hdata->Write();
    hcov->Write();
    hcor->Write();
    hpar->Write();
  }
  

  fout->Close();
  
  fin->Close();

  TRandom3* ran = new TRandom3();
  delete ran;
  sw.Stop();
  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;
  return 1;
}
