#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include "TRandom3.h"
#include "TVector.h"
#include "TVectorT.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include <TMatrixD.h>
#include <TStopwatch.h>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <boost/program_options.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;
using namespace ROOT;
typedef ROOT::VecOps::RVec<double> RVecD;
using ROOT::RDF::RNode; 

using namespace boost::program_options;

constexpr double MW = 0.500;
constexpr double GW = 1.0; 
constexpr int NMASS = 10;
constexpr double DELTAM = 0.2;

auto cheb = [](double x, double scale, double offset, unsigned int n, unsigned int m){
  double den = 0.;
  double num = 0.;
  for(unsigned int i = 0; i <= n ; i++){
    int sign = i%2==0 ? +1 :-1;
    double xj = (TMath::Cos((n-i)*TMath::Pi()/n) + offset)*scale;
    if(x==xj) return 1.0;// protect from nan      
    double val = sign/(x-xj);
    if(i==0 || i==n) val *= 0.5;
    den += val;
    if(i==m) num = val;
  }
  //std::cout << x << "==>" <<  num << "," << den << std::endl;                                             
  return num/den;
};


int main(int argc, char* argv[])
{

  TStopwatch sw;
  sw.Start();

  ROOT::EnableImplicitMT();

  variables_map vm;
  try
    {
      options_description desc{"Options"};
      desc.add_options()
	("help,h", "Help screen")
	("nevents",     value<long>()->default_value(1000), "number of events")
	("lumi",        value<long>()->default_value(1000), "number of events")
	("scalejac",    value<double>()->default_value(10.), "scale")
	("degs_x",      value<int>()->default_value(2), "max degree in x of corrxy")
	("tag",         value<std::string>()->default_value("closure"), "run type")
	("run",         value<std::string>()->default_value("closure"), "run type")
	("do_fit",      bool_switch()->default_value(false), "do_fit")
	("seed",        value<int>()->default_value(4357), "seed");

      store(parse_command_line(argc, argv, desc), vm);
      notify(vm);
      if (vm.count("help")){
	std::cout << desc << '\n';
	return 0;
      }
      if (vm.count("nevents"))    std::cout << "Number of events: " << vm["nevents"].as<long>() << '\n';
      if (vm.count("tag"))        std::cout << "Tag: " << vm["tag"].as<std::string>() << '\n';
      if (vm.count("run"))        std::cout << "Run: " << vm["run"].as<std::string>() << '\n';
      if (vm.count("degs_x"))     std::cout << "Degree in x of corrxy: " << vm["degs_x"].as<int>() << '\n';
    }
  catch (const error &ex)
    {
      std::cerr << ex.what() << '\n';
    }

  long nevents    = vm["nevents"].as<long>();
  long lumi       = vm["lumi"].as<long>();
  double scalejac = vm["scalejac"].as<double>();
  std::string tag = vm["tag"].as<std::string>();
  std::string run = vm["run"].as<std::string>();
  int degs_x      = vm["degs_x"].as<int>();
  int seed        = vm["seed"].as<int>();
  bool do_fit     = vm["do_fit"].as<bool>();

  if(vm.count("degs_x")) tag += std::string(Form("_%d", degs_x));
      
  int x_nbins   = 100; 
  double x_low  = 0.0;
  double x_high = 1.0;

  unsigned int njacs = degs_x + 1;
  
  auto toy_mass = [&](double x, double M, double G){
    return 1./TMath::Pi()/(1 + (x-M)*(x-M)/(G*G/4))*2./G; // non-relativistic
  };

  TFile* fout = TFile::Open(("root/toy_"+tag+"_"+run+".root").c_str(), "RECREATE");

  ROOT::RDataFrame d(nevents);

  unsigned int nslots = d.GetNSlots();
  std::vector<TRandom3*> rans = {};
  for(unsigned int i = 0; i < nslots; i++){
    rans.emplace_back( new TRandom3(seed + i*10) );
  }

  auto dlast = std::make_unique<RNode>(d);

  dlast = std::make_unique<RNode>(dlast->DefineSlot("x", [&](unsigned int nslot)->double{
    double out;
    out = rans[nslot]->Uniform(0.0, 1.0);
    return out;
  } ));
  
  dlast = std::make_unique<RNode>(dlast->Define("weights_mass", 
						[&](double x)->RVecD{
						  RVecD out;
						  double gen = toy_mass(x, MW, GW);
						  out.emplace_back( toy_mass(x, MW, GW)/gen  );
						  for(unsigned int i=0; i<NMASS; i++)
						    out.emplace_back( toy_mass(x, MW - DELTAM*0.5 + DELTAM/(NMASS-1)*i,GW)/gen );
						  return out; 
						}, {"x"}));

  dlast = std::make_unique<RNode>(dlast->Define("weights_jac", 
						[&](double x)->RVecD{
						  RVecD out;						  
						  for(unsigned int k = 0; k<=degs_x; k++){
						    double corrx = cheb(x, 0.5, 1.0, degs_x, k);
						    out.emplace_back( corrx );						    
						  }
						  return out;
						}, {"x"} ));
  
  std::vector<ROOT::RDF::RResultPtr<TH1D> > histos1D;

  dlast = std::make_unique<RNode>(dlast->Define("weight_jacM",[&](double x)->double {
    double out = 4*TMath::Pi()*(1./TMath::Pi()/(1 + (x-MW)*(x-MW)/(GW*GW/4))*2./GW)*2*(x-MW)/GW;
    return out;
  }, {"x"} ));
        
  dlast = std::make_unique<RNode>(dlast->Define("weight", [](RVecD weights_mass){ return weights_mass.at(0);}, {"weights_mass"} ));
  for(unsigned int i=0; i<NMASS; i++){
    dlast = std::make_unique<RNode>(dlast->Define(Form("weight_mass%d", i), [i](RVecD weights_mass){ return weights_mass.at(1+i);}, {"weights_mass"} ));
  }
  for(unsigned int i = 0; i < njacs; i++){
    dlast = std::make_unique<RNode>(dlast->Define(Form("weight_jac%d",i), [i](RVecD weights_jac, RVecD weights_mass){ return weights_jac.at(i)*weights_mass.at(0);},
						  {"weights_jac", "weights_mass"} ));
  }
    
  histos1D.emplace_back(dlast->Histo1D({"h", "nominal", x_nbins, x_low, x_high}, "x", "weight"));      
  for(unsigned int i=0; i<NMASS; i++){
    histos1D.emplace_back(dlast->Histo1D({Form("h_mass%d", i),"", x_nbins, x_low, x_high}, "x", Form("weight_mass%d",i)));
  }    
  for(unsigned int i = 0; i < njacs; i++){
    string hname = std::string(Form("jac_%d", i));	
    histos1D.emplace_back(dlast->Histo1D({ Form("h_jac%d",i), hname.c_str(), x_nbins, x_low, x_high}, "x", Form("weight_jac%d",i)));
  }    
  std::string hname = "jac_mass: d(pdf) / dM";
  histos1D.emplace_back(dlast->Histo1D({Form("h_jac%d", njacs), hname.c_str(), x_nbins, x_low, x_high}, "x", "weight_jacM"));

  auto colNames = dlast->GetColumnNames();
  std::cout << colNames.size() << " columns created" << std::endl;

  double total = *(dlast->Count());  

  fout->cd();
  std::cout << "Writing histos..." << std::endl;
  for(auto h : histos1D) h->Write();
  sw.Stop();

  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;
  std::cout << "Total slots: " << dlast->GetNSlots() << std::endl;

  if(do_fit){
    double sf = double(lumi)/double(nevents);
    unsigned int nbins = x_nbins;
    TH1D* h_nom = (TH1D*)fout->Get("h");
    h_nom->Scale(sf);
    MatrixXd inv_sqrtV = MatrixXd::Zero(nbins, nbins);
    for(unsigned int ib = 0; ib<nbins; ib++)
      inv_sqrtV(ib,ib) = 1./TMath::Sqrt(h_nom->GetBinContent(ib+1));
    MatrixXd jac(nbins, njacs);
    for(unsigned int ij = 0; ij<njacs; ij++){
      TH1D* h_j = (TH1D*)fout->Get(Form("h_jac%d", ij));
      h_j->Scale(sf);
      for(unsigned int ib = 0; ib<nbins; ib++){
	jac(ib, ij) = rans[0]->Gaus(h_j->GetBinContent(ib+1), h_j->GetBinError(ib+1)*scalejac );
      }
    }
    for(unsigned int im = 0; im<NMASS; im++){
      TH1D* h_m = (TH1D*)fout->Get(Form("h_mass%d", im));
      h_m->Scale(sf);
      VectorXd y(nbins);
      for(unsigned int ib=0;ib<nbins; ib++) y(ib) = h_m->GetBinContent(ib+1) - h_nom->GetBinContent(ib+1);
      MatrixXd A = inv_sqrtV*jac;
      VectorXd b = inv_sqrtV*y;
      VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
      MatrixXd chi2old = b.transpose()*b;
      MatrixXd chi2 = ((b - A*x).transpose())*(b-A*x);
      int ndof = nbins-njacs;
      double chi2min = chi2(0,0);
      double chi2norm = chi2(0,0)/ndof;
      cout << "mass " << im << ": chi2_min = " << chi2min << endl;
    }
  }

  
  fout->Close(); 
  for(auto r : rans) delete r;

  return 1;
}
