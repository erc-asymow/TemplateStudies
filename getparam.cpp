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
	("tag",         value<std::string>()->default_value("default"), "tag name")
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
  
  TFile* f = TFile::Open( "root/histos_default.root", "RECREATE");
  TH2D* h_UL = new TH2D("histo_UL", "", 20, 0, 1.0, 20, 0., 1.0 );
  for(int ibx=1; ibx<=h_UL->GetXaxis()->GetNbins() ; ibx++ ){
    for(int iby=1; iby<=h_UL->GetYaxis()->GetNbins() ; iby++ ){
      h_UL->SetBinContent(ibx,iby, 1.0);
    }
  }
  h_UL->Write();
  f->Close();
  
  TFile* fin = TFile::Open(("root/histos_"+tag+".root").c_str(), "READ");
  if(fin==0 || fin==nullptr || fin->IsZombie()){
    cout << "File NOT found" << endl;
    return 0;
  }

  TFile *fout = TFile::Open("fout.root", "RECREATE");

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

  std::vector<TString> proc = {"UL",
    //"A0", "A1", "A2", "A3", "A4"
  };
  std::map<TString, std::array<int,2> > deg_map;
  std::map<TString, std::array<int,2> > par_map;
  std::map<TString, std::array<int,2> > ctr_map;
  deg_map.insert  ( std::make_pair<TString, std::array<int,2> >("UL", {4,  6}) );
  par_map.insert  ( std::make_pair<TString, std::array<int,2> >("UL", {0, +1}) );
  ctr_map.insert  ( std::make_pair<TString, std::array<int,2> >("UL", {0,  0}) ); 
  
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
    
    cheb_x->SetParameter("n", deg_map[iproc].at(0) + 1);
    cheb_y->SetParameter("n", deg_map[iproc].at(1) + 1);
    
    int count_d=0;
    VectorXd y(nd);
    MatrixXd J(nd,np);

    MatrixXd V_inv_sqrt = MatrixXd::Zero(nd, nd);
    for(int i=0; i<nd; i++){
      V_inv_sqrt(i,i) = 100.;
    }
  
    // loop over data bins
    cout << "Start the loop..." << endl;
    for(unsigned int idx = 0; idx<X_nbins; idx++){
      for(unsigned int idy = 0; idy<Y_nbins; idy++){
	y(count_d) = h->GetBinContent(idx+1, idy+1);
	cout << "Doing bin number... " << count_d << endl;

	// loop over parameters
	int count_p=0;
	for(unsigned int ipx = 0; ipx<deg_map[iproc].at(0) + 1; ipx++){
	  // constraint at 0
	  if( ctr_map[iproc].at(0) && ipx==0){
	    continue;
	  }	  
	  cheb_x->SetParameter("m", ipx);
	  double totx = cheb_x->Integral( X_edges[idx], X_edges[idx+1],1.e-10 );
	  
	  for(unsigned int ipy = 0; ipy<deg_map[iproc].at(1) + 1; ipy++){
	    
	    cheb_y->SetParameter("m", ipy);
	    double cheby1 = cheb_y->Integral( Y_edges[idy], Y_edges[idy+1],1.e-10 ); 
	    cheb_y->SetParameter("m", deg_map[iproc].at(1) - ipy);
	    double cheby2 = cheb_y->Integral( Y_edges[idy], Y_edges[idy+1],1.e-10 ); 
	    
	    // even function:
	    if( par_map[iproc].at(1)>0 && ipy < deg_map[iproc].at(1)/2 + 1){
	      double toty = cheby1+cheby2 ;
	      cout << "\t(mX=" << ipx << ", mY=[" << ipy << "," << deg_map[iproc].at(1) - ipy << "])" << endl;
	      cout << "\t > int_x [" << X_edges[idx] << "," << X_edges[idx+1] << "] = " << totx << endl;
	      cout << "\t > int_y [" << Y_edges[idy] << "," << Y_edges[idy+1] << "] = (" << cheby1  << " + " << cheby2 << ")";
	      if( deg_map[iproc].at(1)%2==0 && ipy==deg_map[iproc].at(1)/2 ){
		toty *= 0.5;
		cout << "/2" << endl;
	      }
	      else cout << endl;
	      J(count_d, count_p) = totx*toty;
	      cout << "\t >>> jac(" << count_d << "," << count_p << ") = " << J(count_d, count_p) << endl;
	      count_p++;
	    }

	    // odd function
	    else if( par_map[iproc].at(1)<0 && ipy < (deg_map[iproc].at(1)+1)/2){
	      double toty = cheby1 - cheby2 ;
	      J(count_d, count_p) = totx*toty;	     	      
	      cout << "\t(mX=" << ipx << ", mY=[" << ipy << "," << deg_map[iproc].at(1) - ipy << "])" << endl;
	      cout << "\t > int_x [" << X_edges[idx] << "," << X_edges[idx+1] << "] = " << totx << endl;
	      cout << "\t > int_y [" << Y_edges[idy] << "," << Y_edges[idy+1] << "] = (" << cheby1  << " - " << cheby2 << ")" << endl;
	      cout << "\t >>> jac(" << count_d << "," << count_p << ") = " << J(count_d, count_p) << endl;
	      count_p++;
	    }
	  }
	}
	count_d++;
      }
    }

    MatrixXd A = V_inv_sqrt*J;
    MatrixXd b = V_inv_sqrt*y;
    VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
    VectorXd delta = y-J*x;

    count_d=0;
    TH2D* hpull = (TH2D*)h_UL->Clone("hpull_"+iproc);
    cout << "Start the loop..." << endl;
    for(unsigned int idx = 0; idx<X_nbins; idx++){
      for(unsigned int idy = 0; idy<Y_nbins; idy++){
	hpull->SetBinContent(idx+1,idy+1, delta(count_d) - h_UL->GetBinContent(idx+1,idy+1));
	count_d++;
      }
    }
    fout->cd();
    hpull->Write();
  }
  



  fout->Close();
  
  fin->Close();

  TRandom3* ran = new TRandom3();
  delete ran;
  sw.Stop();
  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;
  return 1;
}
