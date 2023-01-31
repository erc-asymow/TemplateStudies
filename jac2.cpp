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

std::vector<std::string> helicities = {"A0", "A1", "A2", "A3", "A4"};

constexpr double MW = 80.;
constexpr double GW = 2.0; 
constexpr double MASSSHIFT = 0.01; 
constexpr int NMAX  = 200;
//constexpr int NMASS = 20;
constexpr int NMASS = 3;
constexpr double DELTAM = 0.200;

enum pdf_type { pdf_x=0, pdf_y, corr_x, corr_y, A0_x, A0_y, A1_x, A1_y, A2_x, A2_y, A3_x, A3_y, A4_x, A4_y, unknown};

auto get_pdf_type = [](const std::string& name) {
  if     (name=="A0_x") return pdf_type::A0_x; 
  else if(name=="A0_y") return pdf_type::A0_y;
  else if(name=="A1_x") return pdf_type::A1_x; 
  else if(name=="A1_y") return pdf_type::A1_y;
  else if(name=="A2_x") return pdf_type::A2_x; 
  else if(name=="A2_y") return pdf_type::A2_y;
  else if(name=="A3_x") return pdf_type::A3_x; 
  else if(name=="A3_y") return pdf_type::A3_y;
  else if(name=="A4_x") return pdf_type::A4_x; 
  else if(name=="A4_y") return pdf_type::A4_y;
  else return pdf_type::unknown;
};


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
	("degs_corr_x", value<int>()->default_value(2), "max degree in x of corrxy")
	("degs_corr_y", value<int>()->default_value(2), "max degree in y of corrxy")
	("degs_A0_x",   value<int>()->default_value(2), "max degree in x for A0")
	("degs_A0_y",   value<int>()->default_value(2), "max degree in y for A0")
	("degs_A1_x",   value<int>()->default_value(2), "max degree in x for A1")
	("degs_A1_y",   value<int>()->default_value(2), "max degree in y for A1")
	("degs_A2_x",   value<int>()->default_value(2), "max degree in x for A2")
	("degs_A2_y",   value<int>()->default_value(2), "max degree in y for A2")
	("degs_A3_x",   value<int>()->default_value(2), "max degree in x for A3")
	("degs_A3_y",   value<int>()->default_value(2), "max degree in y for A3")
	("degs_A4_x",   value<int>()->default_value(2), "max degree in x for A4")
	("degs_A4_y",   value<int>()->default_value(2), "max degree in y for A4")
	("max_x",       value<double>()->default_value(0.4), "range in x")
	("max_y",       value<double>()->default_value(2.5), "range in y")
	("tag", value<std::string>()->default_value(""), "tag name")
	("run", value<std::string>()->default_value("closure"), "run type")
	("do_jac_vs_mass",    bool_switch()->default_value(false), "compute jacobians for three mass values")
	("smear",  bool_switch()->default_value(false), "smear")
	//("getbin_extTH2_corr",  bool_switch()->default_value(false), "get bin corr(x,y)")
	//("inter_extTH2_corr",  bool_switch()->default_value(false), "interpol corr(x,y)")
	//("toyTF2_corr",  bool_switch()->default_value(false), "toy TF2 corr(x,y)")
	("seed", value<int>()->default_value(4357), "seed")
	("fit_qt_y",  bool_switch()->default_value(false), "fit qt vs y")
        ("relativistic",  bool_switch()->default_value(false), "relativistic");

      store(parse_command_line(argc, argv, desc), vm);
      notify(vm);
      if (vm.count("help")){
	std::cout << desc << '\n';
	return 0;
      }
      if (vm.count("nevents"))    std::cout << "Number of events: " << vm["nevents"].as<long>() << '\n';
      if (vm.count("tag"))        std::cout << "Tag: " << vm["tag"].as<std::string>() << '\n';
      if (vm.count("run"))        std::cout << "Run: " << vm["run"].as<std::string>() << '\n';
      if (vm.count("degs_corr_y")) std::cout << "Degree in x of corrxy: " << vm["degs_corr_x"].as<int>() << '\n';
      if (vm.count("degs_corr_y")) std::cout << "Degree in y of corrxy: " << vm["degs_corr_y"].as<int>() << '\n';
      for(auto hel : helicities){
	if (vm.count("degs_"+hel+"_x")) std::cout << "Degree in x of "+hel+": " << vm["degs_"+hel+"_x"].as<int>() << '\n';
	if (vm.count("degs_"+hel+"_y")) std::cout << "Degree in y of "+hel+": " << vm["degs_"+hel+"_y"].as<int>() << '\n';
      }
    }
  catch (const error &ex)
    {
      std::cerr << ex.what() << '\n';
    }

  long nevents    = vm["nevents"].as<long>();
  std::string tag = vm["tag"].as<std::string>();
  std::string run = vm["run"].as<std::string>();
  int degs_corr_x = vm["degs_corr_x"].as<int>();
  int degs_corr_y = vm["degs_corr_y"].as<int>();
  int degs_A0_x   = vm["degs_A0_x"].as<int>();
  int degs_A0_y   = vm["degs_A0_y"].as<int>();
  int degs_A1_x   = vm["degs_A1_x"].as<int>();
  int degs_A1_y   = vm["degs_A1_y"].as<int>();
  int degs_A2_x   = vm["degs_A2_x"].as<int>();
  int degs_A2_y   = vm["degs_A2_y"].as<int>();
  int degs_A3_x   = vm["degs_A3_x"].as<int>();
  int degs_A3_y   = vm["degs_A3_y"].as<int>();
  int degs_A4_x   = vm["degs_A4_x"].as<int>();
  int degs_A4_y   = vm["degs_A4_y"].as<int>();
  bool do_jac_vs_mass     = vm["do_jac_vs_mass"].as<bool>();
  bool smear      = vm["smear"].as<bool>();
  //bool getbin_extTH2_corr = vm["getbin_extTH2_corr"].as<bool>();
  //bool inter_extTH2_corr  = vm["inter_extTH2_corr"].as<bool>();
  //bool toyTF2_corr        = vm["toyTF2_corr"].as<bool>();
  bool fit_qt_y           = vm["fit_qt_y"].as<bool>();
  int seed = vm["seed"].as<int>();

  if(vm.count("degs_corr_x")) tag += std::string(Form("_UL_%d", degs_corr_x));
  if(vm.count("degs_corr_y")) tag += std::string(Form("_%d", degs_corr_y));
  if(vm.count("degs_A0_x"))   tag += std::string(Form("_A0_%d", degs_A0_x));
  if(vm.count("degs_A0_y"))   tag += std::string(Form("_%d", degs_A0_y));
  if(vm.count("degs_A1_x"))   tag += std::string(Form("_A1_%d", degs_A1_x));
  if(vm.count("degs_A1_y"))   tag += std::string(Form("_%d", degs_A1_y));
  if(vm.count("degs_A2_x"))   tag += std::string(Form("_A2_%d", degs_A2_x));
  if(vm.count("degs_A2_y"))   tag += std::string(Form("_%d", degs_A2_y));
  if(vm.count("degs_A3_x"))   tag += std::string(Form("_A3_%d", degs_A3_x));
  if(vm.count("degs_A3_y"))   tag += std::string(Form("_%d", degs_A3_y));
  if(vm.count("degs_A4_x"))   tag += std::string(Form("_A4_%d", degs_A4_x));
  if(vm.count("degs_A4_y"))   tag += std::string(Form("_%d", degs_A4_y));

  const double deltapOp = 0.015;
  const double deltakOk = 0.0002;
  const double deltah   = 0.001;
    
  double max_x = vm["max_x"].as<double>();
  double max_y = vm["max_y"].as<double>();

  bool relativistic = vm["relativistic"].as<bool>();

  bool do_cheb_as_modifiers = false;
  if(run=="corr") do_cheb_as_modifiers = true;
  
  int nbins_rap   = 25; 
  double rap_low  = 0.0;
  double rap_high = +2.5;
  int nbins_pt    = 35; 
  double pt_low   = 25.;
  double pt_high  = 60.;

  if(fit_qt_y){
    nbins_rap   = 13; 
    rap_low  = -max_y;
    rap_high = +max_y;
    nbins_pt   = 15; 
    pt_low  = 0.;
    pt_high = +max_x;
  }

  auto degs = [degs_corr_x,degs_corr_y,degs_A0_x,degs_A0_y,degs_A1_x,degs_A1_y,degs_A2_x,degs_A2_y,degs_A3_x,degs_A3_y,degs_A4_x,degs_A4_y]
    (const pdf_type& pdf){
    switch(pdf){
    case pdf_type::corr_x:            // corr(x,.) 
    if(degs_corr_x>0) return degs_corr_x;
    else return 2; 
    case pdf_type::corr_y:            // corr(.,y)
    if(degs_corr_y>0) return degs_corr_y;
    else return 2;
    case pdf_type::A0_x:              // A0(x,.)
    if(degs_A0_x>0) return degs_A0_x; 
    else return  2;
    case pdf_type::A0_y:              // A0(.,y)
    if(degs_A0_y>0) return degs_A0_y; 
    else return  2;    
    case pdf_type::A1_x:              // A1(x,.)
    if(degs_A1_x>0) return degs_A1_x; 
    else return  2;
    case pdf_type::A1_y:              // A1(.,y)
    if(degs_A1_y>0) return degs_A1_y; 
    else return  2;    
    case pdf_type::A2_x:              // A2(x,.)
    if(degs_A2_x>0) return degs_A2_x; 
    else return  2;
    case pdf_type::A2_y:              // A2(.,y)
    if(degs_A2_y>0) return degs_A2_y; 
    else return  2;    
    case pdf_type::A3_x:              // A3(x,.)
    if(degs_A3_x>0) return degs_A3_x; 
    else return  2;
    case pdf_type::A3_y:              // A3(.,y)
    if(degs_A3_y>0) return degs_A3_y; 
    else return  2;    
    case pdf_type::A4_x:              // A4(x,.)
    if(degs_A4_x>0) return degs_A4_x; 
    else return  2;
    case pdf_type::A4_y:              // A4(.,y)
    if(degs_A4_y>0) return degs_A4_y; 
    else return  2;    
    default: return 1;
    }
  };

  unsigned int njacs = 0;
  unsigned int first_jac_corrxy = 0;
  unsigned int first_jac_A0xy   = 0;
  unsigned int first_jac_A1xy   = 0;
  unsigned int first_jac_A2xy   = 0;
  unsigned int first_jac_A3xy   = 0;
  unsigned int first_jac_A4xy   = 0;

  unsigned int njacs_corrxy_x = 0;
  unsigned int njacs_corrxy_y = 0;
  unsigned int njacs_A0_x = 0;
  unsigned int njacs_A0_y = 0;
  unsigned int njacs_A1_x = 0;
  unsigned int njacs_A1_y = 0;
  unsigned int njacs_A2_x = 0;
  unsigned int njacs_A2_y = 0;
  unsigned int njacs_A3_x = 0;
  unsigned int njacs_A3_y = 0;
  unsigned int njacs_A4_x = 0;
  unsigned int njacs_A4_y = 0;

  //A1(.,0)=0 & odd
  //A3(.,0)=0 & even
  //A4(.,0)=0 & odd + A4(0,.)!=0
  if(run=="full"){

    assert( degs(pdf_type::A3_y)%2 == 0 );
    
    njacs_corrxy_x = degs(pdf_type::corr_x);
    njacs_corrxy_y = degs(pdf_type::corr_y)/2 + 1;
    njacs_A0_x     = degs(pdf_type::A0_x) ;
    njacs_A0_y     = degs(pdf_type::A0_y)/2 + 1;
    njacs_A1_x     = degs(pdf_type::A1_x);
    njacs_A1_y     = (degs(pdf_type::A1_y)+1)/2;
    njacs_A2_x     = degs(pdf_type::A2_x);
    njacs_A2_y     = degs(pdf_type::A2_y)/2 + 1;
    njacs_A3_x     = degs(pdf_type::A3_x);
    njacs_A3_y     = degs(pdf_type::A3_y)/2;
    njacs_A4_x     = degs(pdf_type::A4_x) + 1;
    njacs_A4_y     = (degs(pdf_type::A4_y)+1)/2;
  }
  else if(run=="corr"){
    njacs_corrxy_x = degs(pdf_type::corr_x) + 1;
    njacs_corrxy_y = degs(pdf_type::corr_y)/2 + 1;
    njacs_A0_x     = degs(pdf_type::A0_x) + 1;
    njacs_A0_y     = degs(pdf_type::A0_y)/2 + 1;
    njacs_A1_x     = degs(pdf_type::A1_x) + 1;
    njacs_A1_y     = degs(pdf_type::A1_y)/2 + 1;
    njacs_A2_x     = degs(pdf_type::A2_x) + 1;
    njacs_A2_y     = degs(pdf_type::A2_y)/2 + 1;
    njacs_A3_x     = degs(pdf_type::A3_x) + 1;
    njacs_A3_y     = degs(pdf_type::A3_y)/2 + 1;
    njacs_A4_x     = degs(pdf_type::A4_x) + 1;
    njacs_A4_y     = degs(pdf_type::A4_y)/2 + 1;
  }
  else if(run=="grid"){
    njacs_corrxy_x = degs(pdf_type::corr_x);
    njacs_corrxy_y = degs(pdf_type::corr_y);
    njacs_A0_x     = njacs_corrxy_x;
    njacs_A0_y     = njacs_corrxy_y;
    njacs_A1_x     = njacs_corrxy_x;
    njacs_A1_y     = njacs_corrxy_y;
    njacs_A2_x     = njacs_corrxy_x;
    njacs_A2_y     = njacs_corrxy_y;
    njacs_A3_x     = njacs_corrxy_x;
    njacs_A3_y     = njacs_corrxy_y;
    njacs_A4_x     = njacs_corrxy_x;
    njacs_A4_y     = njacs_corrxy_y;
  }

  unsigned int njacs_corrxy = njacs_corrxy_x*njacs_corrxy_y;
  unsigned int njacs_A0xy   = njacs_A0_x*njacs_A0_y;
  unsigned int njacs_A1xy   = njacs_A1_x*njacs_A1_y;
  unsigned int njacs_A2xy   = njacs_A2_x*njacs_A2_y;
  unsigned int njacs_A3xy   = njacs_A3_x*njacs_A3_y;
  unsigned int njacs_A4xy   = njacs_A4_x*njacs_A4_y; 
  njacs += njacs_corrxy;
  first_jac_A0xy = njacs;  
  njacs += njacs_A0xy;
  first_jac_A1xy = njacs;  
  njacs += njacs_A1xy;
  first_jac_A2xy = njacs;  
  njacs += njacs_A2xy;
  first_jac_A3xy = njacs;  
  njacs += njacs_A3xy;
  first_jac_A4xy = njacs;  
  njacs += njacs_A4xy;
  
  auto toy_mass = [&](double Q, double M, double G){
    if (relativistic){
      double gamma = TMath::Sqrt(M*M*(M*M+G*G));
      double k = 2*TMath::Sqrt2()*M*G*gamma/TMath::Pi()/TMath::Sqrt(M*M+gamma);
      return k/((Q*Q-M*M)*(Q*Q-M*M)+M*M*G*G);
    }
    return 1./TMath::Pi()/(1 + (Q-M)*(Q-M)/(G*G/4)); // non-relativistic
  };

  TF1* tf1toy_x = new TF1("toy_x", "[0]*x/TMath::Power(x*x+[1], 1.0)", 0.0, max_x);  
  double p0_x = +2.35e-03;
  tf1toy_x->SetParameter(0, 1.0);
  tf1toy_x->SetParameter(1, p0_x);
  double int_toy_x = tf1toy_x->Integral(0.0, max_x);
  tf1toy_x->SetParameter(0, 1.0/int_toy_x);
  auto toy_x = [&](double x)->double{
    return x/TMath::Power(x*x + p0_x, 1.0)/int_toy_x;
  };
  
  TF1* tf1toy_y = new TF1("toy_y", "[2]/TMath::Sqrt(2*TMath::Pi()*[1])*TMath::Exp(-0.5*(x-[0])*(x-[0])/[1])", -max_y, max_y);
  double sigma2_y = 4.0*4.0;
  tf1toy_y->SetParameter(0, 0.0);
  tf1toy_y->SetParameter(1, sigma2_y);
  tf1toy_y->SetParameter(2, 1.0);
  double int_toy_y = tf1toy_y->Integral(-max_y, max_y);
  tf1toy_y->SetParameter(2, 1.0/int_toy_y);
  auto toy_y = [&](double y)->double{
    return 1.0/int_toy_y/TMath::Sqrt(2*TMath::Pi()*sigma2_y)*TMath::Exp(-0.5*y*y/sigma2_y);
  };

  /*
  TFile* fWJets = TFile::Open("WJets.root","READ");
  TH2F* th2_corrxy = 0;
  */
  TF2* tf2toy_corrxy = new TF2("toy_corrxy","0.1*(1.0 - 0.3*y*y)*TMath::Erf(5.0*x) + 1.0", 0., max_x, -max_y, +max_y);
  auto toy_corrxy = [](double x, double y)->double{
    return 0.1*(1.0 - 0.3*y*y)*TMath::Erf(5.0*x) + 1.0;
  };
  //if(fWJets!=nullptr && !fWJets->IsZombie()){
  //  std::cout << "WJets file found! Taking h2 as corrxy" << std::endl;
  //  th2_corrxy = fWJets->Get<TH2F>("h2");
  //}

  TF2* tf2toy_A0 = new TF2("toy_A0", "2*x*x*(1 - 0.01*y*y)",                    0., max_x, -max_y, max_y);  
  TF2* tf2toy_A1 = new TF2("toy_A1", "(0.5*x + 2*x*x)*( 0.05*y + 0.002*y*y*y)", 0., max_x, -max_y, max_y);
  TF2* tf2toy_A2 = new TF2("toy_A2", "2*x*x*(1 - 0.01*y*y)",                    0., max_x, -max_y, max_y);
  //TF2* tf2toy_A3 = new TF2("toy_A3", "(x + x*x + x*x*x)*(1 - 0.01*y*y)",        0., max_x, -max_y, max_y);
  TF2* tf2toy_A3 = new TF2("toy_A3", "0.3*(x + x*x + x*x*x)*(0.1*y*y)",         0., max_x, -max_y, max_y);
  TF2* tf2toy_A4 = new TF2("toy_A4", "(1-x)*(y + y*y*y/10.)/5",                 0., max_x, -max_y, max_y);

  auto toy_A0 = [](double x, double y)->double{ return 2*x*x*(1 - 0.01*y*y); };
  auto toy_A1 = [](double x, double y)->double{ return (0.5*x + 2*x*x)*(0.05*y + 0.002*y*y*y); };
  auto toy_A2 = [](double x, double y)->double{ return 2*x*x*(1 - 0.01*y*y); };
  auto toy_A3 = [](double x, double y)->double{ return 0.3*(x + x*x + x*x*x)*(0.1*y*y); };
  auto toy_A4 = [](double x, double y)->double{ return (1-x)*(y + y*y*y/10.)/5; };
  
  
  // preprare inputs
  if(true){

    TFile* fout = TFile::Open(("root/input_"+tag+".root").c_str(), "RECREATE");
    TTree* tree = new TTree("tree", "tree");
    
    double corr_xy[NMAX*NMAX];
    for(int k = 0; k<=degs(pdf_type::corr_x); k++){
      double x = (TMath::Cos((degs(pdf_type::corr_x)-k)*TMath::Pi()/degs(pdf_type::corr_x))+1.0)*0.5*max_x;
      for(int l = 0; l<=degs(pdf_type::corr_y); l++){      
	int idx = (degs(pdf_type::corr_y)+1)*k + l; 
	tree->Branch(Form("corrxy_%d_%d", k,l), &(corr_xy[idx]), Form("corrxy_%d_%d/D", k,l));
	double y = TMath::Cos((degs(pdf_type::corr_y)-l)*TMath::Pi()/degs(pdf_type::corr_y))*max_y;
	corr_xy[idx] = toy_x(x)*toy_y(y);
	//if(getbin_extTH2_corr)     corr_xy[idx] *= th2_corrxy->GetBinContent( th2_corrxy->FindBin(TMath::Abs(y), x) );
	//else if(inter_extTH2_corr) corr_xy[idx] *= th2_corrxy->Interpolate(TMath::Abs(y),x);
	//else if(toyTF2_corr)
	corr_xy[idx] *= toy_corrxy(x,y);
      }
    }

    // for A0 we sample in [-max_y,max_y]
    double A0_xy[NMAX];
    for(int m = 0; m<=degs(pdf_type::A0_x); m++){
      double x = (TMath::Cos((degs(pdf_type::A0_x)-m)*TMath::Pi()/degs(pdf_type::A0_x))+1.0)*0.5*max_x;
      for(int n = 0; n<=degs(pdf_type::A0_y); n++){
	int idx = (degs(pdf_type::A0_y)+1)*m + n; 
	tree->Branch(Form("A0_xy_%d_%d", m,n), &(A0_xy[idx]), Form("A0_xy_%d_%d/D", m,n));
	double y = TMath::Cos((degs(pdf_type::A0_y)-n)*TMath::Pi()/degs(pdf_type::A0_y))*max_y;
	A0_xy[idx] = toy_A0(x,y);
      }
    }

    // for A1 we sample in [0,max_y]
    double A1_xy[NMAX];
    for(int m = 0; m<=degs(pdf_type::A1_x); m++){
      double x = (TMath::Cos((degs(pdf_type::A1_x)-m)*TMath::Pi()/degs(pdf_type::A1_x))+1.0)*0.5*max_x;
      for(int n = 0; n<=degs(pdf_type::A1_y); n++){
	int idx = (degs(pdf_type::A1_y)+1)*m + n; 
	tree->Branch(Form("A1_xy_%d_%d", m,n), &(A1_xy[idx]), Form("A1_xy_%d_%d/D", m,n));
	double y = TMath::Cos((degs(pdf_type::A1_y)-n)*TMath::Pi()/degs(pdf_type::A1_y))*max_y; 
	A1_xy[idx] = toy_A1(x,y);
      }
    }

    // for A2 we sample in [-max_y,max_y]
    double A2_xy[NMAX];
    for(int m = 0; m<=degs(pdf_type::A2_x); m++){
      double x = (TMath::Cos((degs(pdf_type::A2_x)-m)*TMath::Pi()/degs(pdf_type::A2_x))+1.0)*0.5*max_x;
      for(int n = 0; n<=degs(pdf_type::A2_y); n++){
	int idx = (degs(pdf_type::A2_y)+1)*m + n; 
	tree->Branch(Form("A2_xy_%d_%d", m,n), &(A2_xy[idx]), Form("A2_xy_%d_%d/D", m,n));
	double y = TMath::Cos((degs(pdf_type::A2_y)-n)*TMath::Pi()/degs(pdf_type::A2_y))*max_y;
	A2_xy[idx] = toy_A2(x,y);
      }
    }

    // for A3 we sample in [-max_y,max_y]
    double A3_xy[NMAX];
    for(int m = 0; m<=degs(pdf_type::A3_x); m++){
      double x = (TMath::Cos((degs(pdf_type::A3_x)-m)*TMath::Pi()/degs(pdf_type::A3_x))+1.0)*0.5*max_x;
      for(int n = 0; n<=degs(pdf_type::A3_y); n++){
	int idx = (degs(pdf_type::A3_y)+1)*m + n; 
	tree->Branch(Form("A3_xy_%d_%d", m,n), &(A3_xy[idx]), Form("A3_xy_%d_%d/D", m,n));
	double y = TMath::Cos((degs(pdf_type::A3_y)-n)*TMath::Pi()/degs(pdf_type::A3_y))*max_y;
	A3_xy[idx] = toy_A3(x,y);
      }
    }

    // for A4 we sample in [0,max_y]
    double A4_xy[NMAX];
    for(int m = 0; m<=degs(pdf_type::A4_x); m++){
      double x = (TMath::Cos((degs(pdf_type::A4_x)-m)*TMath::Pi()/degs(pdf_type::A4_x))+1.0)*0.5*max_x;
      for(int n = 0; n<=degs(pdf_type::A4_y); n++){
	int idx = (degs(pdf_type::A4_y)+1)*m + n; 
	tree->Branch(Form("A4_xy_%d_%d", m,n), &(A4_xy[idx]), Form("A4_xy_%d_%d/D", m,n));
	double y = TMath::Cos((degs(pdf_type::A4_y)-n)*TMath::Pi()/degs(pdf_type::A4_y))*max_y;
	A4_xy[idx] = toy_A4(x,y);
      }
    }

    tree->Fill();
    tree->Write();
    fout->Close();
    
    if(nevents<0) return 0;
  }


  // We prepare an input tree to run on
  TFile* fout = TFile::Open(("root/histos_"+tag+"_"+run+".root").c_str(), "RECREATE");

  ROOT::RDataFrame d(nevents);

  unsigned int nslots = d.GetNSlots();
  std::vector<TRandom3*> rans = {};
  for(unsigned int i = 0; i < nslots; i++){
    rans.emplace_back( new TRandom3(seed + i*10) );
  }

  auto dlast = std::make_unique<RNode>(d);
  std::vector<ROOT::RDF::RResultPtr<double> > sums = {};

  dlast = std::make_unique<RNode>(dlast->DefineSlot("PS",[&](unsigned int nslot)->RVecD{
    RVecD out;
    bool accept = false;
    double Q{0.0};
    double cos{0.0};
    double phi{0.0};
    double x{0.0};
    double y{0.0};
    double pt{0.0};
    double eta{0.0};
    while(!accept){
      Q = TMath::Tan(rans[nslot]->Uniform(-TMath::Pi()*0.5, +TMath::Pi()*0.5))*(GW/2) + MW; 
      while(Q < 50.0){
	Q = TMath::Tan(rans[nslot]->Uniform(-TMath::Pi()*0.5, +TMath::Pi()*0.5))*(GW/2) + MW; 
      }
      cos = rans[nslot]->Uniform(-1.0, 1.0);
      phi = rans[nslot]->Uniform(-TMath::Pi(), +TMath::Pi());
      x   = rans[nslot]->Uniform(0.0, max_x);
      y   = rans[nslot]->Uniform(-max_y, max_y);      
      double qT2 = x*x*Q*Q;
      double qT = TMath::Sqrt(qT2);
      double XT = TMath::Sqrt(qT2 + Q*Q);
      double Q0 = XT*TMath::CosH(y);
      double Q3 = XT*TMath::SinH(y);
      double sin = TMath::Sqrt(1-cos*cos);
      Eigen::Matrix4d A;
      A <<      Q0/Q, -qT/Q,   0.0,       -Q3/Q, 
	-qT*Q0/Q/XT,  XT/Q,  0.0,  qT*Q3/Q/XT, 
	0.0,   0.0,  1.0,         0.0,
	-Q3/XT,   0.0,  0.0,       Q0/XT;
      Eigen::Vector4d p4_CS;
      p4_CS << Q/2, Q/2*sin*TMath::Cos(phi), Q/2*sin*TMath::Sin(phi), Q/2*cos ;
      Eigen::Vector4d p4_lab = A.inverse()*p4_CS;
      pt     = TMath::Sqrt(p4_lab(1)*p4_lab(1) + p4_lab(2)*p4_lab(2));
      eta = TMath::ASinH(p4_lab(3)/pt);
      accept |= (pt>=pt_low && pt<=pt_high && TMath::Abs(eta)<=rap_high ) ; 
    }
    out.emplace_back( Q );
    out.emplace_back( cos );
    out.emplace_back( phi );
    out.emplace_back( x );
    out.emplace_back( y );
    out.emplace_back( pt );		
    out.emplace_back( TMath::Abs(eta) );
    double pt_smear  = rans[nslot]->Gaus(pt*(1.0 + deltakOk), pt*deltapOp);
    double eta_smear = eta + rans[nslot]->Gaus(0.0, deltah);
    out.emplace_back( pt_smear );		
    out.emplace_back( TMath::Abs(eta_smear) );
    return out;
  } ));

  dlast = std::make_unique<RNode>(dlast->Define("Q",   [](RVecD PS){return PS.at(0);}, {"PS"}));
  dlast = std::make_unique<RNode>(dlast->Define("cos", [](RVecD PS){return PS.at(1);}, {"PS"}));
  dlast = std::make_unique<RNode>(dlast->Define("phi", [](RVecD PS){return PS.at(2);}, {"PS"}));
  dlast = std::make_unique<RNode>(dlast->Define("x",   [](RVecD PS){return PS.at(3);}, {"PS"}));
  dlast = std::make_unique<RNode>(dlast->Define("y",   [](RVecD PS){return PS.at(4);}, {"PS"}));
  dlast = std::make_unique<RNode>(dlast->Define("pt",  [](RVecD PS){return PS.at(5);}, {"PS"}));
  dlast = std::make_unique<RNode>(dlast->Define("eta", [](RVecD PS){return PS.at(6);}, {"PS"}));
  dlast = std::make_unique<RNode>(dlast->Define("pt_smear",  [](RVecD PS){return PS.at(7);}, {"PS"}));
  dlast = std::make_unique<RNode>(dlast->Define("eta_smear", [](RVecD PS){return PS.at(8);}, {"PS"}));
  
  dlast = std::make_unique<RNode>(dlast->Define("weights_mass", 
						[&](double Q)->RVecD{
						  RVecD out;
						  double gen = 1./TMath::Pi()/(1 + (Q-MW)*(Q-MW)/(GW*GW/4)); //conversion factor
						  out.emplace_back( toy_mass(Q,MW,GW)/gen );
						  out.emplace_back( toy_mass(Q,MW+MASSSHIFT,GW)/gen );
						  out.emplace_back( toy_mass(Q,MW-MASSSHIFT,GW)/gen );
						  for(unsigned int i=0; i<NMASS; i++)
						    out.emplace_back( toy_mass(Q, MW - DELTAM*0.5 + DELTAM/(NMASS-1)*i,GW)/gen );
						  return out; 
						}, {"Q"}));

  dlast = std::make_unique<RNode>(dlast->Define("corrxy_vec", 
						[&](double x, double y)->RVecD{
						  RVecD out;						  
						  for(unsigned int k = (do_cheb_as_modifiers ? 0 : 1); k<=degs(pdf_type::corr_x); k++){
						    double corrx = cheb(x, 0.5*max_x, 1.0, degs(pdf_type::corr_x), k);
						    unsigned int deg = degs(pdf_type::corr_y);
						    unsigned int mid_deg = deg/2;
						    for(unsigned int l = 0; l<=mid_deg; l++){
						      double cheb_l = cheb(y, max_y, 0.0, deg, l) + cheb(y, max_y, 0.0, deg, deg-l) ;
						      double alpha_l = l<mid_deg ? 1.0 : (deg%2==0 ? 0.5 : 1.0);
						      out.emplace_back( corrx*(cheb_l*alpha_l) );
						    }
						  }
						  return out;
						}, {"x","y"} ));
  for(auto hel : helicities){
    dlast = std::make_unique<RNode>(dlast->Define(hel+"xy_vec", 
						  [&,hel](double x, double y)->RVecD{
						    RVecD out;
						    bool is_odd = !do_cheb_as_modifiers && (hel=="A1" || hel=="A4");
						    unsigned int k0 = (do_cheb_as_modifiers || hel=="A4") ? 0 : 1;
						    for(unsigned int k = k0; k<=degs(get_pdf_type(hel+"_x")); k++){
						      double Ax = cheb(x, 0.5*max_x, 1.0, degs(get_pdf_type(hel+"_x")), k);
						      unsigned int deg = degs(get_pdf_type(hel+"_y"));
						      unsigned int up_deg = is_odd ? (deg+1)/2 : deg/2+1 ;
						      bool set_mid_to_zero = run=="full" && hel=="A3";
						      if(set_mid_to_zero) up_deg--;
						      for(unsigned int l = 0; l<up_deg; l++){
							double cheb_l = (cheb(y, max_y, 0.0, deg, l) + (is_odd? -1 : +1)*cheb(y, max_y, 0.0, deg, deg-l)) ;
							double alpha_l = (!is_odd && l==(up_deg-1) && deg%2==0 && !set_mid_to_zero) ? 0.5 : 1.0;;
							out.emplace_back( Ax*(cheb_l*alpha_l) );
						      }						      
						    }
						    return out;
						  }, {"x","y"} ));
  }
  
  std::vector<ROOT::RDF::RResultPtr<TH1D> > histos1D;
  std::vector<ROOT::RDF::RResultPtr<TH2D> > histos2D;
  std::vector<ROOT::RDF::RResultPtr<TH2D> > histosJac;

  double poi_val[NMAX];
  int poi_cat[NMAX];
  unsigned int poi_idx[NMAX];
  unsigned int poi_counter = 0;
  double points_x[ NMAX ];
  double points_y[ NMAX ];
  unsigned int n_pdfx = 0;
  unsigned int n_pdfy = 0;

  if(run!="grid"){
    n_pdfx = degs(pdf_type::corr_x);
    n_pdfy = degs(pdf_type::corr_y)/2 + 1;
    for(unsigned int i = 0; i<n_pdfx; i++){
      points_x[i] = (TMath::Cos((degs(pdf_type::corr_x)-(i+1))*TMath::Pi()/degs(pdf_type::corr_x))+1.0)*0.5*max_x;
    }
    for(unsigned int j = 0; j<n_pdfy; j++){
      points_y[j] = TMath::Cos((degs(pdf_type::corr_y)-j)*TMath::Pi()/degs(pdf_type::corr_y))*max_y;
    }
  }
  else{
    n_pdfx = degs(pdf_type::corr_x);
    n_pdfy = degs(pdf_type::corr_y);
  }
  
  if(true){

    TFile* fin = TFile::Open(("root/input_"+tag+".root").c_str(), "READ");
    if(fin==nullptr) return 0;
    TTree* tree = fin->Get<TTree>("tree");

    double corr_xy[NMAX*NMAX];
    for(int k = 0; k<=degs(pdf_type::corr_x); k++){
      for(int l = 0; l<=degs(pdf_type::corr_y); l++){      
	int idx = (degs(pdf_type::corr_y)+1)*k + l; 
	tree->SetBranchAddress(Form("corrxy_%d_%d", k,l), &(corr_xy[idx]));
      }
    }
    double A0_xy[NMAX];
    for(int m = 0; m<=degs(pdf_type::A0_x); m++){
      for(int n = 0; n<=degs(pdf_type::A0_y); n++){
	int idx = (degs(pdf_type::A0_y)+1)*m + n; 
	tree->SetBranchAddress(Form("A0_xy_%d_%d", m,n), &(A0_xy[idx]));
      }
    }
    double A1_xy[NMAX];
    for(int m = 0; m<=degs(pdf_type::A1_x); m++){
      for(int n = 0; n<=degs(pdf_type::A1_y); n++){
	int idx = (degs(pdf_type::A1_y)+1)*m + n; 
	tree->SetBranchAddress(Form("A1_xy_%d_%d", m,n), &(A1_xy[idx]));
      }
    }
    double A2_xy[NMAX];
    for(int m = 0; m<=degs(pdf_type::A2_x); m++){
      for(int n = 0; n<=degs(pdf_type::A2_y); n++){
	int idx = (degs(pdf_type::A2_y)+1)*m + n; 
	tree->SetBranchAddress(Form("A2_xy_%d_%d", m,n), &(A2_xy[idx]));
      }
    }
    double A3_xy[NMAX];
    for(int m = 0; m<=degs(pdf_type::A3_x); m++){
      for(int n = 0; n<=degs(pdf_type::A3_y); n++){
	int idx = (degs(pdf_type::A3_y)+1)*m + n; 
	tree->SetBranchAddress(Form("A3_xy_%d_%d", m,n), &(A3_xy[idx]));
      }
    }
    double A4_xy[NMAX];
    for(int m = 0; m<=degs(pdf_type::A4_x); m++){
      for(int n = 0; n<=degs(pdf_type::A4_y); n++){
	int idx = (degs(pdf_type::A4_y)+1)*m + n; 
	tree->SetBranchAddress(Form("A4_xy_%d_%d", m,n), &(A4_xy[idx]));
      }
    }

    tree->GetEntry(0);    
    fin->Close();     

    RVecD corrxy_in;
    RVecD A0xy_in;
    RVecD A1xy_in;
    RVecD A2xy_in;
    RVecD A3xy_in;
    RVecD A4xy_in;
    RVecD points_x_in;
    RVecD points_y_in;

    if(run!="grid"){

      int k0 = do_cheb_as_modifiers ? 0 : 1;
      
      for(int k = k0; k<(njacs_corrxy_x+k0); k++){
	for(int l = 0; l<njacs_corrxy_y; l++){      
	  int idx = (degs(pdf_type::corr_y) + 1)*k + l; 
	  corrxy_in.emplace_back( corr_xy[idx] );
	  poi_val[poi_counter] = corr_xy[idx];
	  poi_cat[poi_counter] = -1;
	  poi_idx[poi_counter] = idx;	
	  poi_counter++;
	}
      }
            
      for(int k = k0; k<(njacs_A0_x+k0); k++){
	for(int l = 0; l<njacs_A0_y; l++){      
	int idx = (degs(pdf_type::A0_y)+1)*k + l; 
	A0xy_in.emplace_back( A0_xy[idx] );
	poi_val[poi_counter] = A0_xy[idx];
	poi_cat[poi_counter] = 0;
	poi_idx[poi_counter] = idx;	
	poi_counter++;
	}
      }

      for(int k = k0; k<(njacs_A1_x+k0); k++){
	for(int l = 0; l<njacs_A1_y; l++){      
	  int idx = (degs(pdf_type::A1_y)+1)*k + l; 
	  A1xy_in.emplace_back( A1_xy[idx] );
	  poi_val[poi_counter] = A1_xy[idx];
	  poi_cat[poi_counter] = 1;
	  poi_idx[poi_counter] = idx;	
	  poi_counter++;
	}
      }
      
      for(int k = k0; k<(njacs_A2_x+k0); k++){
	for(int l = 0; l<njacs_A2_y; l++){      
	  int idx = (degs(pdf_type::A2_y)+1)*k + l; 
	  A2xy_in.emplace_back( A2_xy[idx] );
	  poi_val[poi_counter] = A2_xy[idx];
	  poi_cat[poi_counter] = 2;
	  poi_idx[poi_counter] = idx;	
	  poi_counter++;
	}
      }

      for(int k = k0; k<(njacs_A3_x+k0); k++){
	for(int l = 0; l<njacs_A3_y; l++){      
	  int idx = (degs(pdf_type::A3_y)+1)*k + l; 
	  A3xy_in.emplace_back( A3_xy[idx] );
	  poi_val[poi_counter] = A3_xy[idx];
	  poi_cat[poi_counter] = 3;
	  poi_idx[poi_counter] = idx;	
	  poi_counter++;
	}
      }
      
      for(int k = 0; k<njacs_A4_x; k++){
	for(int l = 0; l<njacs_A4_y; l++){      
	  int idx = (degs(pdf_type::A4_y)+1)*k + l; 
	  A4xy_in.emplace_back( A4_xy[idx] );
	  poi_val[poi_counter] = A4_xy[idx];
	  poi_cat[poi_counter] = 4;
	  poi_idx[poi_counter] = idx;	
	  poi_counter++;
	}
      }
    }    
    else if(run=="grid"){      
      points_x_in.emplace_back(0.0);
      points_x[0] = 0.0;
      for(unsigned int i = 1; i<=n_pdfx; i++){
	double val = max_x/n_pdfx*i;
	points_x[i] = val;
	points_x_in.emplace_back( val );
      }
      points_y_in.emplace_back(0.0);
      points_y[0] = 0.0;
      for(unsigned int j = 1; j<=n_pdfy; j++){
	double val = max_y/n_pdfy*j;
	points_y[j] = val;
	points_y_in.emplace_back( val );
      }      
      // unpol
      for(int k = 0; k<n_pdfx; k++){
	for(int l = 0; l<n_pdfy; l++){
	  poi_cat[poi_counter] = -1;
	  poi_idx[poi_counter] = n_pdfy*k + l;
	  poi_counter++;
	}
      }
      // A0
      for(int k = 0; k<n_pdfx; k++){
	for(int l = 0; l<n_pdfy; l++){
	  poi_cat[poi_counter] = 0;
	  poi_idx[poi_counter] = n_pdfy*k + l;
	  poi_counter++;
	}
      }
      // A1
      for(int k = 0; k<n_pdfx; k++){
	for(int l = 0; l<n_pdfy; l++){
	  poi_cat[poi_counter] = 1;
	  poi_idx[poi_counter] = n_pdfy*k + l;	
	  poi_counter++;
	}
      }
      // A2
      for(int k = 0; k<n_pdfx; k++){
	for(int l = 0; l<n_pdfy; l++){
	  poi_cat[poi_counter] = 2;
	  poi_idx[poi_counter] = n_pdfy*k + l;
	  poi_counter++;
	}
      }
      // A3
      for(int k = 0; k<n_pdfx; k++){
	for(int l = 0; l<n_pdfy; l++){
	  poi_cat[poi_counter] = 3;
	  poi_idx[poi_counter] = n_pdfy*k + l;
	  poi_counter++;
	}
      }
      // A4
      for(int k = 0; k<n_pdfx; k++){
	for(int l = 0; l<n_pdfy; l++){
	  poi_cat[poi_counter] = 4;
	  poi_idx[poi_counter] = n_pdfy*k + l;
	  poi_counter++;
	}
      }
    }
        
    // mass
    poi_val[poi_counter] = 0.0;
    poi_cat[poi_counter] = 5;
    poi_idx[poi_counter] = 0;
    poi_counter++;        
    
    dlast = std::make_unique<RNode>(dlast->Define("harmonics", [&](double x, double y, double cos, double phi) -> RVecD{
	  RVecD out;
	  double cosS = cos;//y>0 ? cos : -cos;
	  double phiS = phi;//y>0 ? phi : -phi;
	  double sinS = TMath::Sqrt(1-cosS*cosS);
	  out.emplace_back( 1.0 + cosS*cosS );
	  out.emplace_back( 0.5*(1-3*cosS*cosS) );
	  out.emplace_back( 2*sinS*cosS*TMath::Cos(phiS) );
	  out.emplace_back( 0.5*sinS*sinS*TMath::Cos(2*phiS) );
	  out.emplace_back( sinS*TMath::Cos(phiS) );
	  out.emplace_back( cosS );
	  //out = {1.0, 0., 0., 0., 0., 0., 0.};
	  return out;
	} , {"x", "y", "cos", "phi"} ));

    dlast = std::make_unique<RNode>(dlast->Define("weights_mctruth", [&](double x, double y, RVecD har)->RVecD {
	  RVecD out;
	  double norm{3./16/TMath::Pi()};
	  double UL = toy_x(x)*toy_y(y);
	  //if(getbin_extTH2_corr)     UL *= th2_corrxy->GetBinContent( th2_corrxy->FindBin(TMath::Abs(y),x) );
	  //else if(inter_extTH2_corr) UL *= th2_corrxy->Interpolate(TMath::Abs(y),x);
	  //else if(toyTF2_corr)
	  UL *= toy_corrxy(x,y);	   
	  double A0 = toy_A0(x,y);
	  double A1 = toy_A1(x,y);
	  double A2 = toy_A2(x,y);
	  double A3 = toy_A3(x,y);
	  double A4 = toy_A4(x,y);
	  //cout << "\tULtoy =" << UL << endl;
	  //cout << "\tA0toy =" << A0 << endl;
	  //cout << "\tA1toy =" << A1 << endl;
	  //cout << "\tA2toy =" << A2 << endl;
	  //cout << "\tA3toy =" << A3 << endl;
	  //cout << "\tA4toy =" << A4 << endl;
	  double norm_UL   = norm*UL;
	  double norm_PiAi = norm*(har.at(0)+A0*har.at(1)+A1*har.at(2)+A2*har.at(3)+A3*har.at(4)+A4*har.at(5));
	  out.emplace_back( UL*norm_PiAi);          // full weight 
	  out.emplace_back( UL*norm_PiAi);          // weight for corr
	  out.emplace_back( norm_UL*A0*har.at(1) ); // weight for A0
	  out.emplace_back( norm_UL*A1*har.at(2) ); // weight for A1
	  out.emplace_back( norm_UL*A2*har.at(3) ); // weight for A2
	  out.emplace_back( norm_UL*A3*har.at(4) ); // weight for A3
	  out.emplace_back( norm_UL*A4*har.at(5) ); // weight for A4
	  return out;
    }, {"x", "y", "harmonics"} ));        
    
    dlast = std::make_unique<RNode>(dlast->Define("weights_cheb", [&,corrxy_in,A0xy_in,A1xy_in,A2xy_in,A3xy_in,A4xy_in]
						  (RVecD corrxy_vec,
						   RVecD A0xy_vec,
						   RVecD A1xy_vec,
						   RVecD A2xy_vec,
						   RVecD A3xy_vec,
						   RVecD A4xy_vec,
						   RVecD har )->RVecD {
      /*
      cout << corrxy_vec.size() << ", " << corrxy_in.size() << endl;
      cout << A0xy_vec.size() << ", " << A0xy_in.size() << endl;
      cout << A1xy_vec.size() << ", " << A1xy_in.size() << endl;
      cout << A2xy_vec.size() << ", " << A2xy_in.size() << endl;
      cout << A3xy_vec.size() << ", " << A3xy_in.size() << endl;
      cout << A4xy_vec.size() << ", " << A4xy_in.size() << endl;
      */
      RVecD out;
      double norm{3./16/TMath::Pi()};	  
      double UL = ROOT::VecOps::Dot(corrxy_vec,corrxy_in);
      double A0 = ROOT::VecOps::Dot(A0xy_vec,A0xy_in); 
      double A1 = ROOT::VecOps::Dot(A1xy_vec,A1xy_in); 
      double A2 = ROOT::VecOps::Dot(A2xy_vec,A2xy_in); 
      double A3 = ROOT::VecOps::Dot(A3xy_vec,A3xy_in);
      double A4 = ROOT::VecOps::Dot(A4xy_vec,A4xy_in);  
      //cout << "\tULchb =" << UL << endl;
      //cout << "\tA0chb =" << A0 << endl;
      //cout << "\tA1chb =" << A1 << endl;
      //cout << "\tA2chb =" << A2 << endl;
      //cout << "\tA3chb =" << A3 << endl;
      //cout << "\tA4chb =" << A4 << endl;
      double norm_UL   = norm*UL;      
      double norm_PiAi = norm*(har.at(0)+A0*har.at(1)+A1*har.at(2)+A2*har.at(3)+A3*har.at(4)+A4*har.at(5));
      out.emplace_back( UL*norm_PiAi);       // full weight 
      out.emplace_back( norm_PiAi );         // weight for corr 
      out.emplace_back( norm_UL*har.at(1) ); // weight for A0
      out.emplace_back( norm_UL*har.at(2) ); // weight for A1
      out.emplace_back( norm_UL*har.at(3) ); // weight for A2
      out.emplace_back( norm_UL*har.at(4) ); // weight for A3
      out.emplace_back( norm_UL*har.at(5) ); // weight for A4
      return out;
    }, {"corrxy_vec", 
	"A0xy_vec",
	"A1xy_vec",
	"A2xy_vec",
	"A3xy_vec",
	"A4xy_vec",
	"harmonics"} ));
    
    dlast = std::make_unique<RNode>(dlast->Define("weights_jac", [&]
						  (RVecD corrxy_vec,
						   RVecD A0xy_vec,
						   RVecD A1xy_vec,
						   RVecD A2xy_vec,
						   RVecD A3xy_vec,
						   RVecD A4xy_vec,
						   RVecD weights )->RVecD {
      RVecD out;						
      for(unsigned int i = 0; i<njacs_corrxy; i++) out.emplace_back( weights.at(1)*corrxy_vec.at(i) );
      for(unsigned int i = 0; i<njacs_A0xy; i++)   out.emplace_back( weights.at(2)*A0xy_vec.at(i));
      for(unsigned int i = 0; i<njacs_A1xy; i++)   out.emplace_back( weights.at(3)*A1xy_vec.at(i));
      for(unsigned int i = 0; i<njacs_A2xy; i++)   out.emplace_back( weights.at(4)*A2xy_vec.at(i));
      for(unsigned int i = 0; i<njacs_A3xy; i++)   out.emplace_back( weights.at(5)*A3xy_vec.at(i));
      for(unsigned int i = 0; i<njacs_A4xy; i++)   out.emplace_back( weights.at(6)*A4xy_vec.at(i));      
      return out;
    }, {"corrxy_vec",
	"A0xy_vec",
	"A1xy_vec",
	"A2xy_vec",
	"A3xy_vec",
	"A4xy_vec",
	(do_cheb_as_modifiers ? "weights_mctruth" : "weights_cheb")} ));


    dlast = std::make_unique<RNode>(dlast->Define("weights_jac_grid", [&,points_x_in,points_y_in]
						  (double x, double y, 						   
						   RVecD weights)->RVecD {      
      RVecD out(njacs);						
      double absy = TMath::Abs(y);
      int ijac = -1; 
      if(x==max_x)    x -= 1e-04;
      if(absy==max_y) y -= 1e-04;
      unsigned int counter = 0;
      for(unsigned int i = 0; i<degs(pdf_type::corr_x); i++){
	for(unsigned int j = 0; j<degs(pdf_type::corr_y); j++){
	  bool is_in_bin =  x>=points_x_in.at(i) && x<points_x_in.at(i+1) && absy>=points_y_in.at(j) && absy<points_y_in.at(j+1);
	  if( is_in_bin ) ijac = counter;
	  counter++;
	}
      }
      //cout << ijac << endl;
      if(ijac<0) return out;
      for(unsigned int k=0; k<6; k++){
	double val = 0.0;
	if(k==0)
	  val = (weights.at(0) - weights.at(2) - weights.at(3) - weights.at(4) - weights.at(5) - weights.at(6));
	else
	  val = weights.at(k+1);
	//cout << "out[" << (k+1)*ijac << "]=" << val << endl;
	out[njacs/6*k + ijac] = val;
      }
      return out;    
    }, {"x", "y", "weights_mctruth"} ));

    dlast = std::make_unique<RNode>(dlast->Define("weight_jacM",[&](RVecD weights, double Q)->double {
      double w = weights.at(0);
      w *= TMath::Pi()*toy_mass(Q,MW,GW)*2*(Q-MW)/(GW*GW/4);
      return w;
    }, {(run=="full" ? "weights_cheb" : "weights_mctruth"), "Q"} ));
    
    
    dlast = std::make_unique<RNode>(dlast->Define("w", [](RVecD weights, RVecD weights_mass){ return weights.at(0)*weights_mass.at(0);}, {"weights_cheb", "weights_mass"} ));
    if(do_jac_vs_mass){
      dlast = std::make_unique<RNode>(dlast->Define("w_up",   [](RVecD weights, RVecD weights_mass){ return weights.at(0)*weights_mass.at(1);}, {"weights_cheb", "weights_mass"} ));
      dlast = std::make_unique<RNode>(dlast->Define("w_down", [](RVecD weights, RVecD weights_mass){ return weights.at(0)*weights_mass.at(2);}, {"weights_cheb", "weights_mass"} ));
    }
    dlast = std::make_unique<RNode>(dlast->Define("wMC",      [](RVecD weights, RVecD weights_mass){ return weights.at(0)*weights_mass.at(0);}, {"weights_mctruth", "weights_mass"} ));
    dlast = std::make_unique<RNode>(dlast->Define("wMC_up",   [](RVecD weights, RVecD weights_mass){ return weights.at(0)*weights_mass.at(1);}, {"weights_mctruth", "weights_mass"} ));
    dlast = std::make_unique<RNode>(dlast->Define("wMC_down", [](RVecD weights, RVecD weights_mass){ return weights.at(0)*weights_mass.at(2);}, {"weights_mctruth", "weights_mass"} ));
    for(unsigned int i=0; i<NMASS; i++){
      dlast = std::make_unique<RNode>(dlast->Define(Form("wMC_mass%d", i), [i](RVecD weights, RVecD weights_mass){ return weights.at(0)*weights_mass.at(3+i);}, {"weights_mctruth", "weights_mass"} ));
    }

    for(unsigned int i = 0; i < njacs; i++){
      if(do_jac_vs_mass){
	for(unsigned int j = 0; j < 3; j++){
	  dlast = std::make_unique<RNode>(dlast->Define(Form("jac%d_%d", j, i), [i,j](RVecD weights, RVecD weights_mass){ return weights.at(i)*weights_mass.at(j);},
							{(run=="grid" ? "weights_jac_grid" : "weights_jac"), "weights_mass"} ));
	}
      }
      else{
	dlast = std::make_unique<RNode>(dlast->Define(Form("jac_%d",i), [i](RVecD weights, RVecD weights_mass){ return weights.at(i)*weights_mass.at(0);},
						      {(run=="grid" ? "weights_jac_grid" : "weights_jac"), "weights_mass"} ));
      }
    }
    
    if(run=="full"){
      sums.emplace_back( dlast->Sum<double>("w") );
      sums.emplace_back( dlast->Sum<double>("wMC") );
    }
    
    string varx = fit_qt_y ? "y" : "eta";
    string vary = fit_qt_y ? "x" : "pt";
    if(fit_qt_y==false && smear){
      varx = "eta_smear";
      vary = "pt_smear"; 
    }
    
    if(run=="full"){
      histos2D.emplace_back(dlast->Histo2D({"h",       "", nbins_rap, rap_low, rap_high, nbins_pt, pt_low, pt_high}, varx, vary, "w"));
      if(do_jac_vs_mass){
	histos2D.emplace_back(dlast->Histo2D({"h_up",    "", nbins_rap, rap_low, rap_high, nbins_pt, pt_low, pt_high}, varx, vary, "w_up"));
	histos2D.emplace_back(dlast->Histo2D({"h_down",  "", nbins_rap, rap_low, rap_high, nbins_pt, pt_low, pt_high}, varx, vary, "w_down"));
      }
    }
    histos2D.emplace_back(dlast->Histo2D({"hMC",     "", nbins_rap, rap_low, rap_high, nbins_pt, pt_low, pt_high}, varx, vary, "wMC"));      
    histos2D.emplace_back(dlast->Histo2D({"hMC_up",  "", nbins_rap, rap_low, rap_high, nbins_pt, pt_low, pt_high}, varx, vary, "wMC_up"));      
    histos2D.emplace_back(dlast->Histo2D({"hMC_down","", nbins_rap, rap_low, rap_high, nbins_pt, pt_low, pt_high}, varx, vary, "wMC_down"));
    for(unsigned int i=0; i<NMASS; i++){
      histos2D.emplace_back(dlast->Histo2D({Form("hMC_mass%d", i),"", nbins_rap, rap_low, rap_high, nbins_pt, pt_low, pt_high}, varx, vary, Form("wMC_mass%d",i)));
    }
    
    for(unsigned int i = 0; i < njacs; i++){
      std::string hname = "";
      if(i>=first_jac_corrxy && i<first_jac_A0xy){
	unsigned int idx_x = (i-first_jac_corrxy) / njacs_corrxy_y;
	unsigned int idx_y = (i-first_jac_corrxy) % njacs_corrxy_y;
	if(run=="full") idx_x++;	  
	hname = std::string(Form("jac_%d: d(pdf) / d(corrxy_in[%d][%d])", i, idx_x, idx_y));	
      }
      else if(i>=first_jac_A0xy && i<first_jac_A1xy){
	unsigned int idx_x = (i-first_jac_A0xy) / njacs_A0_y;
	unsigned int idx_y = (i-first_jac_A0xy) % njacs_A0_y;
	if(run=="full") idx_x++;	  
	hname = std::string(Form("jac_%d: d(pdf) / d(A0xy_in[%d][%d])", i, idx_x, idx_y));	
      }
      else if(i>=first_jac_A1xy && i<first_jac_A2xy){
	unsigned int idx_x = (i-first_jac_A1xy) / njacs_A1_y;
	unsigned int idx_y = (i-first_jac_A1xy) % njacs_A1_y;
	if(run=="full") idx_x++;
	hname = std::string(Form("jac_%d: d(pdf) / d(A1xy_in[%d][%d])", i, idx_x, idx_y));	
      }
      else if(i>=first_jac_A2xy && i<first_jac_A3xy){
	unsigned int idx_x = (i-first_jac_A2xy) / njacs_A2_y;
	unsigned int idx_y = (i-first_jac_A2xy) % njacs_A2_y;
	if(run=="full") idx_x++;	  
	hname = std::string(Form("jac_%d: d(pdf) / d(A2xy_in[%d][%d])", i, idx_x, idx_y));	
      }
      else if(i>=first_jac_A3xy && i<first_jac_A4xy){
	unsigned int idx_x = (i-first_jac_A3xy) / njacs_A3_y;
	unsigned int idx_y = (i-first_jac_A3xy) % njacs_A3_y;
	if(run=="full") idx_x++;	  
	hname = std::string(Form("jac_%d: d(pdf) / d(A3xy_in[%d][%d])", i, idx_x, idx_y));	
      }
      else if(i>=first_jac_A4xy){
	unsigned int idx_x = (i-first_jac_A4xy) / njacs_A4_y;
	unsigned int idx_y = (i-first_jac_A4xy) % njacs_A4_y;
	//if(run=="full") idx_y++;
	hname = std::string(Form("jac_%d: d(pdf) / d(A4xy_in[%d][%d])", i, idx_x, idx_y));	
      }
      
      if(do_jac_vs_mass){
	for(unsigned int j = 0; j < 3; j++){
	  histosJac.emplace_back(dlast->Histo2D({ Form("jac%d_%d",j,i), hname.c_str(), nbins_rap, rap_low, rap_high, nbins_pt, pt_low, pt_high}, varx, vary,
						Form("jac%d_%d",j,i)));
	}
      }
      else{
	histosJac.emplace_back(dlast->Histo2D({ Form("jac_%d",i), hname.c_str(), nbins_rap, rap_low, rap_high, nbins_pt, pt_low, pt_high}, varx, vary,
					      Form("jac_%d",i)));
      }            
    }

    // This is the last jacobian
    std::string hname = "jac_mass: d(pdf) / dM";
    histosJac.emplace_back(dlast->Histo2D({Form("jac_%d", njacs), hname.c_str(), nbins_rap, rap_low, rap_high, nbins_pt, pt_low, pt_high}, varx, vary,
					  "weight_jacM"));
    
    /*
    histos1D.emplace_back(dlast->Histo1D({"w_pdfx",    "", 20, 0.0, max_x}, "x", "w"));      
    histos1D.emplace_back(dlast->Histo1D({"w_pdfy",    "", 20, 0.0, max_y}, "y", "w"));      
    histos2D.emplace_back(dlast->Histo2D({"w_corrxy",  "", 20, 0.0, max_y, 20, 0.0, max_x}, "y", "x", "w"));      
    histos1D.emplace_back(dlast->Histo1D({"wMC_pdfx",  "", 20, 0.0, max_x}, "x", "wMC"));      
    histos1D.emplace_back(dlast->Histo1D({"wMC_pdfy",  "", 20, 0.0, max_y}, "y", "wMC"));      
    histos2D.emplace_back(dlast->Histo2D({"wMC_corrxy","", 20, 0.0, max_y, 20, 0.0, max_x}, "y", "x", "wMC"));      
    */

    dlast = std::make_unique<RNode>(dlast->Define("ibin",
						  [&](double pt, double eta){
						    int ibin = histos2D.at(0)->FindBin(eta,pt);
						    //cout << "ibin " << ibin << endl;
						    return ibin;
						  }, {"pt", "eta"}));
    dlast = std::make_unique<RNode>(dlast->Define("ibin_smear",
						  [&](double pt, double eta){
						    int ibin = histos2D.at(0)->FindBin(eta,pt);
						    //cout << "ibin " << ibin << endl;
						    return ibin;
						  }, {"pt_smear", "eta_smear"}));
    int nbins_transfer = histos2D.at(0)->GetNcells(); 
    //histos2D.emplace_back(dlast->Histo2D({"transfer_MC", "", nbins_transfer, 0, double(nbins_transfer), nbins_transfer, 0, double(nbins_transfer)}, "ibin", "ibin_smear", "wMC"));
    
    fin->Close();
  }


  auto colNames = dlast->GetColumnNames();
  std::cout << colNames.size() << " columns created" << std::endl;

  double total = *(dlast->Count());  
  for(auto sum : sums) std::cout << (*sum)*(max_x*2*max_y*4*TMath::Pi())/total << std::endl;

  fout->cd();
  std::cout << "Writing histos..." << std::endl;
  for(auto h : histos1D) h->Write();
  for(auto h : histos2D) h->Write();
  for(auto h : histosJac){
    //continue;
    h->Scale(1./total);
    h->Write();
  }

  //clone existing histos (for consistency)
  if(run!="full"){
    for(auto h : histos2D){
      if( std::string(h->GetName())=="hMC" )      ((TH2D*)h->Clone("h"))->Write();
      if( std::string(h->GetName())=="hMC_up" )   ((TH2D*)h->Clone("h_up"))->Write();
      if( std::string(h->GetName())=="hMC_down" ) ((TH2D*)h->Clone("h_down"))->Write();
    }
  }
  
  sw.Stop();
  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;
  std::cout << "Total slots: " << dlast->GetNSlots() << std::endl;

  TTree* outtree = new TTree("outtree", "tree");

  outtree->Branch("n_pdfx",   &n_pdfx, "n_pdfx/i");
  outtree->Branch("n_pdfy",   &n_pdfy, "n_pdfy/i");
  outtree->Branch("points_x", &points_x, "points_x[n_pdfx]/D");
  outtree->Branch("points_y", &points_y, "points_y[n_pdfy]/D");
  outtree->Branch("poi_counter", &poi_counter, "poi_counter/i");
  outtree->Branch("poi_val", &poi_val, "poi_val[poi_counter]/D");
  outtree->Branch("poi_cat", &poi_cat, "poi_cat[poi_counter]/I");
  outtree->Branch("poi_idx", &poi_idx, "poi_idx[poi_counter]/i");
  outtree->Branch("nevents", &nevents, "nevents/L");

  outtree->Fill();
  outtree->Write();
  
  fout->Close();

  for(auto r : rans) delete r;

  return 1;
}
