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
constexpr double MASSSHIFT = 0.050;
constexpr int NMAX  = 200;
constexpr int NMASS = 20;
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
	("tag", value<std::string>()->default_value(""), "tag name")
	("run", value<std::string>()->default_value("closure"), "run type")
	("do_absy",    bool_switch()->default_value(false), "polycheb in abs(y)")
	("getbin_extTH2_corr",  bool_switch()->default_value(false), "get bin corr(x,y)")
	("inter_extTH2_corr",  bool_switch()->default_value(false), "interpol corr(x,y)")
	("toyTF2_corr",  bool_switch()->default_value(false), "toy TF2 corr(x,y)")
	("seed", value<int>()->default_value(4357), "seed")
	("fit_qt_y",  bool_switch()->default_value(false), "fit qt vs y");

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
      if (vm.count("do_absy"))     std::cout << "Do abs(Y): " << vm["do_absy"].as<bool>() << '\n';
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
  bool do_absy = vm["do_absy"].as<bool>();
  bool getbin_extTH2_corr= vm["getbin_extTH2_corr"].as<bool>();
  bool inter_extTH2_corr= vm["inter_extTH2_corr"].as<bool>();
  bool toyTF2_corr= vm["toyTF2_corr"].as<bool>();
  bool fit_qt_y = vm["fit_qt_y"].as<bool>();
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

  const double max_x = 0.4;
  const double max_y = 2.5;
  
  int nbinsX   = 36; 
  double xLow  = 0.0;
  double xHigh = +2.5;
  int nbinsY   = 60; 
  double yLow  = 25.;
  double yHigh = 55.;

  if(fit_qt_y){
    nbinsX   = 13; 
    xLow  = -max_y;
    xHigh = +max_y;
    nbinsY   = 15; 
    yLow  = 0.;
    yHigh = +max_x;
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
  unsigned int first_jac_corrxy = njacs;  
  njacs += (degs(pdf_type::corr_x))*( (do_absy ? degs(pdf_type::corr_y) : degs(pdf_type::corr_y)/2) + 1 );
  unsigned int first_jac_A0xy = njacs;  
  njacs += (degs(pdf_type::A0_x))*( (do_absy ? degs(pdf_type::A0_y) : degs(pdf_type::A0_y)/2) + 1 );
  unsigned int first_jac_A1xy = njacs;  
  njacs += (degs(pdf_type::A1_x))*( degs(pdf_type::A1_y) );
  unsigned int first_jac_A2xy = njacs;  
  njacs += (degs(pdf_type::A2_x))*( (do_absy ? degs(pdf_type::A2_y) : degs(pdf_type::A2_y)/2) + 1 );
  unsigned int first_jac_A3xy = njacs;  
  njacs += (degs(pdf_type::A3_x))*( (do_absy ? degs(pdf_type::A3_y) : degs(pdf_type::A3_y)/2) + 1 );
  unsigned int first_jac_A4xy = njacs;  
  njacs += (degs(pdf_type::A4_x)+1)*( degs(pdf_type::A4_y) );

  auto toy_mass = [](double Q, double M, double G){
    return 1./TMath::Pi()/(1 + (Q-M)*(Q-M)/G/G);
  };

  TF1* toy_x = new TF1("toy_x", "[0]*x/(x*x+[1])", 0.0, max_x);  
  toy_x->SetParameter(0, 1.0);
  toy_x->SetParameter(1, +2.35e-03);
  double int_toy_x = toy_x->Integral(0.0, max_x);
  toy_x->SetParameter(0, 1.0/int_toy_x);
  //TF1* toy_x = new TF1("toy_x", "[0]", 0.0, max_x);  
  //toy_x->SetParameter(0, 1./max_x);

  TF1* toy_y = new TF1("toy_y", "[2]/TMath::Sqrt(2*TMath::Pi()*[1])*TMath::Exp(-0.5*(x-[0])*(x-[0])/[1]/[1])", -max_y, max_y);
  toy_y->SetParameter(0, 0.0);
  toy_y->SetParameter(1, 5.0);
  toy_y->SetParameter(2, 1.0);
  double int_toy_y = toy_y->Integral(-max_y, max_y);
  toy_y->SetParameter(2, 1.0/int_toy_y);

  TFile* fWJets = TFile::Open("WJets.root","READ");
  TH2F* th2_corrxy = 0;
  TF2* toy_corrxy = new TF2("toy_corrxy","0.1*(1.0 - 0.3*x*x)*TMath::Erf(5.0*y) + 1.0", 0., max_y, 0., max_x);
  if(fWJets!=nullptr && !fWJets->IsZombie()){
    std::cout << "WJets file found! Taking h2 as corrxy" << std::endl;
    th2_corrxy = fWJets->Get<TH2F>("h2");
  }

  /*
  TF2* toy_A0 = new TF2("toy_A0", "2*y*y*(1/(1 + 0.02*x*x))", 0., max_y, 0., max_x);
  TF2* toy_A1 = new TF2("toy_A1", "(0.5*y + 2*y*y)*(x/(1 + x))", 0., max_y, 0., max_x);
  TF2* toy_A2 = new TF2("toy_A2", "2*y*y*(1/(1 + 0.02*x*x))", 0., max_y, 0., max_x);
  TF2* toy_A3 = new TF2("toy_A3", "(y + y*y + y*y*y)*(1/(1 + 0.02*x*x))", 0., max_y, 0., max_x);
  TF2* toy_A4 = new TF2("toy_A4", "(y+1)*(0.5*x + x*x)/10", 0., max_y, 0., max_x);
  */

  TF2* toy_A0 = new TF2("toy_A0", "2*y*y*(1 - 0.01*x*x)", 0., max_y, 0., max_x);
  TF2* toy_A1 = new TF2("toy_A1", "(0.5*y + 2*y*y)*( 0.05*x + 0.01*x*x)", 0., max_y, 0., max_x);
  TF2* toy_A2 = new TF2("toy_A2", "2*y*y*(1 - 0.01*x*x)", 0., max_y, 0., max_x);
  TF2* toy_A3 = new TF2("toy_A3", "(y + y*y + y*y*y)*(1 - 0.01*x*x)", 0., max_y, 0., max_x);
  TF2* toy_A4 = new TF2("toy_A4", "(y+1)*(0.5*x + x*x)/6.0", 0., max_y, 0., max_x);

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
	double y = do_absy ? 
	  (TMath::Cos((degs(pdf_type::corr_y)-l)*TMath::Pi()/degs(pdf_type::corr_y))+1.0)*0.5*max_y : 
	  TMath::Cos((degs(pdf_type::corr_y)-l)*TMath::Pi()/degs(pdf_type::corr_y))*max_y;
	corr_xy[idx] = toy_x->Eval(x)*toy_y->Eval(y);
	if(getbin_extTH2_corr) 
	  corr_xy[idx] *= th2_corrxy->GetBinContent( th2_corrxy->FindBin(TMath::Abs(y), x) );
	else if(inter_extTH2_corr)  
	  corr_xy[idx] *= th2_corrxy->Interpolate(TMath::Abs(y),x);
	else if(toyTF2_corr) 
	  corr_xy[idx] *= toy_corrxy->Eval(TMath::Abs(y),x);
      }
    }
  
    double A0_xy[NMAX];
    for(int m = 0; m<=degs(pdf_type::A0_x); m++){
      double x = (TMath::Cos((degs(pdf_type::A0_x)-m)*TMath::Pi()/degs(pdf_type::A0_x))+1.0)*0.5*max_x;
      for(int n = 0; n<=degs(pdf_type::A0_y); n++){
	int idx = (degs(pdf_type::A0_y)+1)*m + n; 
	tree->Branch(Form("A0_xy_%d_%d", m,n), &(A0_xy[idx]), Form("A0_xy_%d_%d/D", m,n));
	double y = do_absy ? (TMath::Cos((degs(pdf_type::A0_y)-n)*TMath::Pi()/degs(pdf_type::A0_y))+1.0)*0.5*max_y : 
	  TMath::Cos((degs(pdf_type::A0_y)-n)*TMath::Pi()/degs(pdf_type::A0_y))*max_y;
	A0_xy[idx] = toy_A0->Eval(TMath::Abs(y),x);      
      }
    }

    double A1_xy[NMAX];
    for(int m = 0; m<=degs(pdf_type::A1_x); m++){
      double x = (TMath::Cos((degs(pdf_type::A1_x)-m)*TMath::Pi()/degs(pdf_type::A1_x))+1.0)*0.5*max_x;
      for(int n = 0; n<=degs(pdf_type::A1_y); n++){
	int idx = (degs(pdf_type::A1_y)+1)*m + n; 
	tree->Branch(Form("A1_xy_%d_%d", m,n), &(A1_xy[idx]), Form("A1_xy_%d_%d/D", m,n));
	//double y = do_absy ? (TMath::Cos((degs(pdf_type::A1_y)-n)*TMath::Pi()/degs(pdf_type::A1_y))+1.0)*0.5*max_y : 
	//TMath::Cos((degs(pdf_type::A1_y)-n)*TMath::Pi()/degs(pdf_type::A1_y))*max_y;
	double y = (TMath::Cos((degs(pdf_type::A1_y)-n)*TMath::Pi()/degs(pdf_type::A1_y))+1.0)*0.5*max_y; 
	A1_xy[idx] = toy_A1->Eval(TMath::Abs(y),x);      
      }
    }

    double A2_xy[NMAX];
    for(int m = 0; m<=degs(pdf_type::A2_x); m++){
      double x = (TMath::Cos((degs(pdf_type::A2_x)-m)*TMath::Pi()/degs(pdf_type::A2_x))+1.0)*0.5*max_x;
      for(int n = 0; n<=degs(pdf_type::A2_y); n++){
	int idx = (degs(pdf_type::A2_y)+1)*m + n; 
	tree->Branch(Form("A2_xy_%d_%d", m,n), &(A2_xy[idx]), Form("A2_xy_%d_%d/D", m,n));
	double y = do_absy ? (TMath::Cos((degs(pdf_type::A2_y)-n)*TMath::Pi()/degs(pdf_type::A2_y))+1.0)*0.5*max_y : 
	  TMath::Cos((degs(pdf_type::A2_y)-n)*TMath::Pi()/degs(pdf_type::A2_y))*max_y;
	A2_xy[idx] = toy_A2->Eval(TMath::Abs(y),x);      
      }
    }

    double A3_xy[NMAX];
    for(int m = 0; m<=degs(pdf_type::A3_x); m++){
      double x = (TMath::Cos((degs(pdf_type::A3_x)-m)*TMath::Pi()/degs(pdf_type::A3_x))+1.0)*0.5*max_x;
      for(int n = 0; n<=degs(pdf_type::A3_y); n++){
	int idx = (degs(pdf_type::A3_y)+1)*m + n; 
	tree->Branch(Form("A3_xy_%d_%d", m,n), &(A3_xy[idx]), Form("A3_xy_%d_%d/D", m,n));
	double y = do_absy ? (TMath::Cos((degs(pdf_type::A3_y)-n)*TMath::Pi()/degs(pdf_type::A3_y))+1.0)*0.5*max_y : 
	  TMath::Cos((degs(pdf_type::A3_y)-n)*TMath::Pi()/degs(pdf_type::A3_y))*max_y;
	A3_xy[idx] = toy_A3->Eval(TMath::Abs(y),x);      
      }
    }

    double A4_xy[NMAX];
    for(int m = 0; m<=degs(pdf_type::A4_x); m++){
      double x = (TMath::Cos((degs(pdf_type::A4_x)-m)*TMath::Pi()/degs(pdf_type::A4_x))+1.0)*0.5*max_x;
      for(int n = 0; n<=degs(pdf_type::A4_y); n++){
	int idx = (degs(pdf_type::A4_y)+1)*m + n; 
	tree->Branch(Form("A4_xy_%d_%d", m,n), &(A4_xy[idx]), Form("A4_xy_%d_%d/D", m,n));
	//double y = do_absy ? (TMath::Cos((degs(pdf_type::A4_y)-n)*TMath::Pi()/degs(pdf_type::A4_y))+1.0)*0.5*max_y : 
	//TMath::Cos((degs(pdf_type::A4_y)-n)*TMath::Pi()/degs(pdf_type::A4_y))*max_y;
	double y = (TMath::Cos((degs(pdf_type::A4_y)-n)*TMath::Pi()/degs(pdf_type::A4_y))+1.0)*0.5*max_y;
	A4_xy[idx] = toy_A4->Eval(TMath::Abs(y),x);      
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

  dlast = std::make_unique<RNode>(dlast->DefineSlot("Q",  [&](unsigned int nslot){ 
	double Q = TMath::Tan(rans[nslot]->Uniform(-TMath::Pi()*0.5, +TMath::Pi()*0.5))*GW + MW; 
	while(Q < 50.0){
	  Q = TMath::Tan(rans[nslot]->Uniform(-TMath::Pi()*0.5, +TMath::Pi()*0.5))*GW + MW; 
	}
	return Q;
      } ));
  dlast = std::make_unique<RNode>(dlast->DefineSlot("cos",[&](unsigned int nslot){ return rans[nslot]->Uniform(-1.0, 1.0);} ));
  dlast = std::make_unique<RNode>(dlast->DefineSlot("phi",[&](unsigned int nslot){ return rans[nslot]->Uniform(-TMath::Pi(), +TMath::Pi());} ));
  dlast = std::make_unique<RNode>(dlast->DefineSlot("x",  [&](unsigned int nslot){ return rans[nslot]->Uniform(0.0, max_x); }));
  dlast = std::make_unique<RNode>(dlast->DefineSlot("y",  [&](unsigned int nslot){ return rans[nslot]->Uniform(-max_y, max_y);}));
  dlast = std::make_unique<RNode>(dlast->Define("weightsM", 
						[&](double Q)->RVecD{
						  RVecD out;
						  double gen = 1./TMath::Pi()/(1 + (Q-MW)*(Q-MW)/GW/GW);
						  out.emplace_back( toy_mass(Q,MW,GW)/gen );
						  out.emplace_back( toy_mass(Q,MW+MASSSHIFT,GW)/gen );
						  out.emplace_back( toy_mass(Q,MW-MASSSHIFT,GW)/gen );
						  for(unsigned int i=0; i<NMASS; i++)
						    out.emplace_back( toy_mass(Q, MW - DELTAM*0.5 + DELTAM/NMASS*i,GW)/gen );
						  return out; 
						}, {"Q"}));
  dlast = std::make_unique<RNode>(dlast->Define("p4lab",
						[&](double Q, double cos, double phi, double x, double y)->RVecD {
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
						  double pt  = TMath::Sqrt(p4_lab(1)*p4_lab(1) + p4_lab(2)*p4_lab(2));
						  double eta = TMath::ASinH(p4_lab(3)/pt);
						  RVecD p4lab{pt,eta};
						  return p4lab;
						}, {"Q", "cos", "phi", "x", "y"}));
    
  dlast = std::make_unique<RNode>(dlast->Define("corrxy_vec", 
						[&](double x, double y)->RVecD{
						  RVecD out;						  
						  for(unsigned int k = 1; k<=degs(pdf_type::corr_x); k++){
						    double corrx = cheb(x, 0.5*max_x, 1.0, degs(pdf_type::corr_x), k);
						    if(do_absy){
						      for(unsigned int l = 0; l<=degs(pdf_type::corr_y); l++){
							double corry = cheb(TMath::Abs(y), 0.5*max_y, 1.0, degs(pdf_type::corr_y), l);
							out.emplace_back( corrx*corry );
						      }
						    }
						    else{
						      unsigned int deg = degs(pdf_type::corr_y);
						      unsigned int mid_deg = deg/2;
						      for(unsigned int l = 0; l<=mid_deg; l++){
							double cheb_l = cheb(y, max_y, 0.0, deg, l) + cheb(y, max_y, 0.0, deg, deg-l) ;
							double alpha_l = l<mid_deg ? 1.0 : (deg%2==0 ? 0.5 : 1.0);
							out.emplace_back( corrx*(cheb_l*alpha_l) );
						      }
						    }
						  }
						  return out;
						}, {"x","y"} ));
  for(auto hel : helicities){
    dlast = std::make_unique<RNode>(dlast->Define(hel+"xy_vec", 
						  [&,hel](double x, double y)->RVecD{
						    RVecD out;
						    for(unsigned int k = (hel!="A4") ? 1: 0; k<=degs(get_pdf_type(hel+"_x")); k++){
						      double Ax = cheb(x, 0.5*max_x, 1.0, degs(get_pdf_type(hel+"_x")), k);
						      if(do_absy || (hel=="A1" || hel=="A4")){
							for(unsigned int l = (hel=="A1" || hel=="A4") ? 1 : 0; l<=degs(get_pdf_type(hel+"_y")); l++){
							  double Ay = cheb(TMath::Abs(y), 0.5*max_y, 1.0, degs(get_pdf_type(hel+"_y")), l);
							  out.emplace_back( Ax*Ay );
							}
						      }
						      else{
							unsigned int deg = degs(get_pdf_type(hel+"_y"));
							unsigned int mid_deg = deg/2;
							for(unsigned int l = 0; l<=mid_deg; l++){
							  double cheb_l = (cheb(y, max_y, 0.0, deg, l) + cheb(y, max_y, 0.0, deg, deg-l)) ;
							  double alpha_l = l<mid_deg ? 1.0 : (deg%2==0 ? 0.5 : 1.0);
							  out.emplace_back( Ax*(cheb_l*alpha_l) );
							}
						      }
						    }
						    return out;
						  }, {"x","y"} ));
  }
  
  dlast = std::make_unique<RNode>(dlast->Define("pt",  [](RVecD p4lab){ return p4lab[0];}, {"p4lab"}));
  if(xLow>=0.)
    dlast = std::make_unique<RNode>(dlast->Define("eta", [](RVecD p4lab){ return TMath::Abs(p4lab[1]);}, {"p4lab"}));
  else
    dlast = std::make_unique<RNode>(dlast->Define("eta", [](RVecD p4lab){ return p4lab[1];}, {"p4lab"}));

  std::vector<ROOT::RDF::RResultPtr<TH1D> > histos1D;
  
  std::vector<ROOT::RDF::RResultPtr<TH2D> > histos2D;
  std::vector<ROOT::RDF::RResultPtr<TH2D> > histosJac;

  double poi_val[NMAX];
  int poi_cat[NMAX];
  unsigned int poi_idx[NMAX];
  unsigned int poi_counter = 0;
  unsigned int n_pdfx = degs(pdf_type::corr_x);
  unsigned int n_pdfy = (do_absy ? degs(pdf_type::corr_y) : degs(pdf_type::corr_y)/2 ) + 1;
  double points_x[ NMAX ];
  double points_y[ NMAX ];

  for(unsigned int i = 0; i<n_pdfx; i++){
    points_x[i] = (TMath::Cos((degs(pdf_type::corr_x)-(i+1))*TMath::Pi()/degs(pdf_type::corr_x))+1.0)*0.5*max_x;
  }
  for(unsigned int j = 0; j<n_pdfy; j++){
    points_y[j]   = do_absy ? 
      (TMath::Cos((degs(pdf_type::corr_y)-j)*TMath::Pi()/degs(pdf_type::corr_y))+1.0)*0.5*max_y : 
      TMath::Cos((degs(pdf_type::corr_y)-j)*TMath::Pi()/degs(pdf_type::corr_y))*max_y;
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
    for(int k = 1; k<=degs(pdf_type::corr_x); k++){
      for(int l = 0; l<=(do_absy ? degs(pdf_type::corr_y) : degs(pdf_type::corr_y)/2); l++){      
	int idx = (degs(pdf_type::corr_y) + 1)*k + l; 
	corrxy_in.emplace_back( corr_xy[idx] );
	poi_val[poi_counter] = corr_xy[idx];
	poi_cat[poi_counter] = -1;
	poi_idx[poi_counter] = idx;	
	poi_counter++;
      }
    }

    RVecD A0xy_in;
    for(int k = 1; k<=degs(pdf_type::A0_x); k++){
      for(int l = 0; l<=(do_absy ? degs(pdf_type::A0_y) : degs(pdf_type::A0_y)/2); l++){      
	int idx = (degs(pdf_type::A0_y)+1)*k + l; 
	A0xy_in.emplace_back( A0_xy[idx] );
	poi_val[poi_counter] = A0_xy[idx];
	poi_cat[poi_counter] = 0;
	poi_idx[poi_counter] = idx;	
	poi_counter++;
      }
    }
    RVecD A1xy_in;
    for(int k = 1; k<=degs(pdf_type::A1_x); k++){
      for(int l = 1; l<=degs(pdf_type::A1_y); l++){      
	int idx = (degs(pdf_type::A1_y)+1)*k + l; 
	A1xy_in.emplace_back( A1_xy[idx] );
	poi_val[poi_counter] = A1_xy[idx];
	poi_cat[poi_counter] = 1;
	poi_idx[poi_counter] = idx;	
	poi_counter++;
      }
    }
    RVecD A2xy_in;
    for(int k = 1; k<=degs(pdf_type::A2_x); k++){
      for(int l = 0; l<=(do_absy ? degs(pdf_type::A2_y) : degs(pdf_type::A2_y)/2); l++){      
	int idx = (degs(pdf_type::A2_y)+1)*k + l; 
	A2xy_in.emplace_back( A2_xy[idx] );
	poi_val[poi_counter] = A2_xy[idx];
	poi_cat[poi_counter] = 2;
	poi_idx[poi_counter] = idx;	
	poi_counter++;
      }
    }
    RVecD A3xy_in;
    for(int k = 1; k<=degs(pdf_type::A3_x); k++){
      for(int l = 0; l<=(do_absy ? degs(pdf_type::A3_y) : degs(pdf_type::A3_y)/2); l++){      
	int idx = (degs(pdf_type::A3_y)+1)*k + l; 
	A3xy_in.emplace_back( A3_xy[idx] );
	poi_val[poi_counter] = A3_xy[idx];
	poi_cat[poi_counter] = 3;
	poi_idx[poi_counter] = idx;	
	poi_counter++;
      }
    }
    RVecD A4xy_in;
    for(int k = 0; k<=degs(pdf_type::A4_x); k++){
      for(int l = 1; l<=degs(pdf_type::A4_y); l++){      
	int idx = (degs(pdf_type::A4_y)+1)*k + l; 
	A4xy_in.emplace_back( A4_xy[idx] );
	poi_val[poi_counter] = A4_xy[idx];
	poi_cat[poi_counter] = 4;
	poi_idx[poi_counter] = idx;	
	poi_counter++;
      }
    }
		
    dlast = std::make_unique<RNode>(dlast->Define("harmonics", [&](double x, double y, double cos, double phi) -> RVecD{
	  RVecD out;
	  double cosS = y>0 ? cos : -cos;
	  double phiS = y>0 ? phi : -phi;
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

    dlast = std::make_unique<RNode>(dlast->Define("weightsMC", [&](double x, double y, RVecD harmonics, RVecD weightsM  )->RVecD {
	  RVecD out;
	  double wMC{3./16/TMath::Pi()};
	  wMC *= toy_x->Eval(x);
	  wMC *= toy_y->Eval(y);
	  double corrxyMC  = 1.0;
	  if(getbin_extTH2_corr) 
	    corrxyMC = th2_corrxy->GetBinContent( th2_corrxy->FindBin(TMath::Abs(y),x) );
	  else if(inter_extTH2_corr)  
	    corrxyMC = th2_corrxy->Interpolate(TMath::Abs(y),x);
	  else if(toyTF2_corr)        
	    corrxyMC = toy_corrxy->Eval(TMath::Abs(y),x);	   
	  wMC *= corrxyMC;
	  wMC *= (harmonics.at(0) + 
	  	  toy_A0->Eval(TMath::Abs(y),x)*harmonics.at(1) +
	  	  toy_A1->Eval(TMath::Abs(y),x)*harmonics.at(2) + 
	  	  toy_A2->Eval(TMath::Abs(y),x)*harmonics.at(3) +
	  	  toy_A3->Eval(TMath::Abs(y),x)*harmonics.at(4) +
	  	  toy_A4->Eval(TMath::Abs(y),x)*harmonics.at(5)
	  	  );
	  out.emplace_back( wMC*weightsM.at(0) );
	  out.emplace_back( wMC*weightsM.at(1) );
	  out.emplace_back( wMC*weightsM.at(2) );
	  for(unsigned int i=0; i<NMASS; i++) out.emplace_back( wMC*weightsM.at(3+i) );
	  return out;
	}, {"x", "y", "harmonics", "weightsM"} ));
    
    dlast = std::make_unique<RNode>(dlast->Define("weights", [&,corrxy_in,A0xy_in,A1xy_in,A2xy_in,A3xy_in,A4xy_in]
						  (RVecD corrxy_vec,
						   RVecD A0xy_vec,
						   RVecD A1xy_vec,
						   RVecD A2xy_vec,
						   RVecD A3xy_vec,
						   RVecD A4xy_vec,
						   RVecD harmonics,
						   RVecD weightsM  )->RVecD {
						    RVecD out;
						    double wUL{3./16/TMath::Pi()};	  
						    wUL *= ROOT::VecOps::Dot(corrxy_vec,corrxy_in);
						    double w = wUL*(harmonics.at(0) + 
								    ROOT::VecOps::Dot(A0xy_vec,A0xy_in)*harmonics.at(1) +
								    ROOT::VecOps::Dot(A1xy_vec,A1xy_in)*harmonics.at(2) +
								    ROOT::VecOps::Dot(A2xy_vec,A2xy_in)*harmonics.at(3) +
								    ROOT::VecOps::Dot(A3xy_vec,A3xy_in)*harmonics.at(4) +
								    ROOT::VecOps::Dot(A4xy_vec,A4xy_in)*harmonics.at(5)
								    ); 
						    out.emplace_back( w*weightsM.at(0) );
						    out.emplace_back( w*weightsM.at(1) );
						    out.emplace_back( w*weightsM.at(2) );
						    return out;
						  }, {"corrxy_vec", 
						      "A0xy_vec",
						      "A1xy_vec",
						      "A2xy_vec",
						      "A3xy_vec",
						      "A4xy_vec",
						      "harmonics",
						      "weightsM"} ));

    dlast = std::make_unique<RNode>(dlast->Define("weights_jac", [&,corrxy_in,A0xy_in,A1xy_in,A2xy_in,A3xy_in,A4xy_in]
						  (RVecD corrxy_vec,
						   RVecD A0xy_vec,
						   RVecD A1xy_vec,
						   RVecD A2xy_vec,
						   RVecD A3xy_vec,
						   RVecD A4xy_vec,
						   RVecD harmonics,
						   RVecD weightsM  )->RVecD {
						    RVecD out;						
						    double wUL{3./16/TMath::Pi()};	  
						    double A = ROOT::VecOps::Dot(corrxy_vec,corrxy_in);
 						    double B = ROOT::VecOps::Dot(A0xy_vec,A0xy_in); 
 						    double C = ROOT::VecOps::Dot(A1xy_vec,A1xy_in); 
 						    double D = ROOT::VecOps::Dot(A2xy_vec,A2xy_in); 
 						    double E = ROOT::VecOps::Dot(A3xy_vec,A3xy_in);
 						    double F = ROOT::VecOps::Dot(A4xy_vec,A4xy_in);  

						    unsigned int njacs_corrxy = degs(pdf_type::corr_x)*( (do_absy ? degs(pdf_type::corr_y) : degs(pdf_type::corr_y)/2) + 1 );
						    for(unsigned int i = 0; i<njacs_corrxy; i++){
						      RVecD corrxy_in_copy( corrxy_in.size(), 0.0 );
						      corrxy_in_copy[i] = 1.0;
						      out.emplace_back( wUL*
									ROOT::VecOps::Dot(corrxy_vec,corrxy_in_copy)*
									(  
									 harmonics.at(0) + 
									 B*harmonics.at(1) +
									 C*harmonics.at(2) +
									 D*harmonics.at(3) +
									 E*harmonics.at(4) +
									 F*harmonics.at(5)
									   )*
									1.0 /*weightsM.at(0)*/ );
						    }						   
						    
 						    unsigned int njacs_A0xy = degs(pdf_type::A0_x)*( (do_absy ? degs(pdf_type::A0_y) : degs(pdf_type::A0_y)/2) + 1 );
						    for(unsigned int i = 0; i<njacs_A0xy; i++){
 						      RVecD A0xy_in_copy( A0xy_in.size(), 0.0 );
						      A0xy_in_copy[i] = 1.0;
						      out.emplace_back( wUL*
									A*
									(
									 ROOT::VecOps::Dot(A0xy_vec,A0xy_in_copy)*harmonics.at(1)
									 )*
									1.0 /*weightsM.at(0)*/ );
						    }						   

						    unsigned int njacs_A1xy = degs(pdf_type::A1_x)*degs(pdf_type::A1_y);
						    for(unsigned int i = 0; i<njacs_A1xy; i++){
 						      RVecD A1xy_in_copy( A1xy_in.size(), 0.0 );
						      A1xy_in_copy[i] = 1.0;
						      out.emplace_back( wUL*
									A*
									(
									 ROOT::VecOps::Dot(A1xy_vec,A1xy_in_copy)*harmonics.at(2)
									 )* 
									1.0 /*weightsM.at(0)*/ );
						    }						   

						    unsigned int njacs_A2xy = degs(pdf_type::A2_x)*( (do_absy ? degs(pdf_type::A2_y) : degs(pdf_type::A2_y)/2) + 1 );
						    for(unsigned int i = 0; i<njacs_A2xy; i++){
 						      RVecD A2xy_in_copy( A2xy_in.size(), 0.0 );
						      A2xy_in_copy[i] = 1.0;
						      out.emplace_back( wUL*
									A*
									(
									 ROOT::VecOps::Dot(A2xy_vec,A2xy_in_copy)*harmonics.at(3)
									 )* 
									1.0 /*weightsM.at(0)*/ );
						    }						   

						    unsigned int njacs_A3xy = degs(pdf_type::A3_x)*( (do_absy ? degs(pdf_type::A3_y) : degs(pdf_type::A3_y)/2) + 1 );
						    for(unsigned int i = 0; i<njacs_A3xy; i++){
 						      RVecD A3xy_in_copy( A3xy_in.size(), 0.0 );
						      A3xy_in_copy[i] = 1.0;
						      out.emplace_back( wUL*
									A*
									(
									 ROOT::VecOps::Dot(A3xy_vec,A3xy_in_copy)*harmonics.at(4)
									 )* 
									1.0 /*weightsM.at(0)*/ );
						    }

						    unsigned int njacs_A4xy = (degs(pdf_type::A4_x)+1)*degs(pdf_type::A4_y);
						    for(unsigned int i = 0; i<njacs_A4xy; i++){
 						      RVecD A4xy_in_copy( A4xy_in.size(), 0.0 );
						      A4xy_in_copy[i] = 1.0;
						      out.emplace_back( wUL*
									A*
									(
									 ROOT::VecOps::Dot(A4xy_vec,A4xy_in_copy)*harmonics.at(5)
									 )* 
									1.0 /*weightsM.at(0)*/ );
						    }						   


						    return out;
						  }, {"corrxy_vec",
						      "A0xy_vec",
						      "A1xy_vec",
						      "A2xy_vec",
						      "A3xy_vec",
						      "A4xy_vec",
						      "harmonics",
						      "weightsM"} ));

    dlast = std::make_unique<RNode>(dlast->Define("w",        [](RVecD weights){ return weights.at(0);}, {"weights"} ));
    dlast = std::make_unique<RNode>(dlast->Define("w_up",     [](RVecD weights){ return weights.at(1);}, {"weights"} ));
    dlast = std::make_unique<RNode>(dlast->Define("w_down",   [](RVecD weights){ return weights.at(2);}, {"weights"} ));
    dlast = std::make_unique<RNode>(dlast->Define("wMC",      [](RVecD weights){ return weights.at(0);}, {"weightsMC"} ));
    dlast = std::make_unique<RNode>(dlast->Define("wMC_up",   [](RVecD weights){ return weights.at(1);}, {"weightsMC"} ));
    dlast = std::make_unique<RNode>(dlast->Define("wMC_down", [](RVecD weights){ return weights.at(2);}, {"weightsMC"} ));
    for(unsigned int i=0; i<NMASS; i++){
      dlast = std::make_unique<RNode>(dlast->Define(Form("wMC_mass%d", i), [i](RVecD weights){ return weights.at(3+i);}, {"weightsMC"} ));
    }

    for(unsigned int i = 0; i < njacs; i++){
      for(unsigned int j = 0; j < 3; j++){
	dlast = std::make_unique<RNode>(dlast->Define(Form("jac%d_%d", j, i), [i,j](RVecD weights, RVecD weightsM){ return weights.at(i)*weightsM.at(j);}, {"weights_jac", "weightsM"} ));
      }
    }
    
    sums.emplace_back( dlast->Sum<double>("w") );
    sums.emplace_back( dlast->Sum<double>("wMC") );

    if(!fit_qt_y){
      histos2D.emplace_back(dlast->Histo2D({"h",       "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", "w"));
      histos2D.emplace_back(dlast->Histo2D({"h_up",    "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", "w_up"));
      histos2D.emplace_back(dlast->Histo2D({"h_down",  "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", "w_down"));      
      histos2D.emplace_back(dlast->Histo2D({"hMC",     "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", "wMC"));      
      histos2D.emplace_back(dlast->Histo2D({"hMC_up",  "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", "wMC_up"));      
      histos2D.emplace_back(dlast->Histo2D({"hMC_down","", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", "wMC_down"));
      for(unsigned int i=0; i<NMASS; i++){
	histos2D.emplace_back(dlast->Histo2D({Form("hMC_mass%d", i),"", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", Form("wMC_mass%d",i)));
      }
    }
    else{
      histos2D.emplace_back(dlast->Histo2D({"h",       "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "y", "x", "w"));      
      histos2D.emplace_back(dlast->Histo2D({"hMC",     "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "y", "x", "wMC"));      
      histos2D.emplace_back(dlast->Histo2D({"hMC_up",  "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "y", "x", "wMC_up"));      
      histos2D.emplace_back(dlast->Histo2D({"hMC_down","", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "y", "x", "wMC_down"));
    }

    for(unsigned int i = 0; i < njacs; i++){
      std::string hname = "";
      if(i>=first_jac_corrxy && i<first_jac_A0xy){
	unsigned int idx_x = (i-first_jac_corrxy) / ((do_absy ? degs(pdf_type::corr_y) : degs(pdf_type::corr_y)/2) +1) + 1;
	unsigned int idx_y = (i-first_jac_corrxy) % ((do_absy ? degs(pdf_type::corr_y) : degs(pdf_type::corr_y)/2) +1);
	hname = std::string(Form("jac_%d: d(pdf) / d(corrxy_in[%d][%d])", i, idx_x, idx_y));	
      }
      else if(i>=first_jac_A0xy && i<first_jac_A1xy){
	unsigned int idx_x = (i-first_jac_A0xy) / ((do_absy ? degs(pdf_type::A0_y) : degs(pdf_type::A0_y)/2) +1) + 1;
	unsigned int idx_y = (i-first_jac_A0xy) % ((do_absy ? degs(pdf_type::A0_y) : degs(pdf_type::A0_y)/2) +1);
	hname = std::string(Form("jac_%d: d(pdf) / d(A0xy_in[%d][%d])", i, idx_x, idx_y));	
      }
      else if(i>=first_jac_A1xy && i<first_jac_A2xy){
	unsigned int idx_x = (i-first_jac_A1xy) / degs(pdf_type::A1_y) + 1;
	unsigned int idx_y = (i-first_jac_A1xy) % degs(pdf_type::A1_y) + 1;
	hname = std::string(Form("jac_%d: d(pdf) / d(A1xy_in[%d][%d])", i, idx_x, idx_y));	
      }
      else if(i>=first_jac_A2xy && i<first_jac_A3xy){
	unsigned int idx_x = (i-first_jac_A2xy) / ((do_absy ? degs(pdf_type::A2_y) : degs(pdf_type::A2_y)/2) +1) + 1;
	unsigned int idx_y = (i-first_jac_A2xy) % ((do_absy ? degs(pdf_type::A2_y) : degs(pdf_type::A2_y)/2) +1);
	hname = std::string(Form("jac_%d: d(pdf) / d(A2xy_in[%d][%d])", i, idx_x, idx_y));	
      }
      else if(i>=first_jac_A3xy && i<first_jac_A4xy){
	unsigned int idx_x = (i-first_jac_A3xy) / ((do_absy ? degs(pdf_type::A3_y) : degs(pdf_type::A3_y)/2) +1) + 1;
	unsigned int idx_y = (i-first_jac_A3xy) % ((do_absy ? degs(pdf_type::A3_y) : degs(pdf_type::A3_y)/2) +1);
	hname = std::string(Form("jac_%d: d(pdf) / d(A3xy_in[%d][%d])", i, idx_x, idx_y));	
      }
      else if(i>=first_jac_A4xy){
	unsigned int idx_x = (i-first_jac_A4xy) / degs(pdf_type::A4_y);
	unsigned int idx_y = (i-first_jac_A4xy) % degs(pdf_type::A4_y) + 1;
	hname = std::string(Form("jac_%d: d(pdf) / d(A4xy_in[%d][%d])", i, idx_x, idx_y));	
      }
      
      if(!fit_qt_y){
	for(unsigned int j = 0; j < 3; j++){
	  histosJac.emplace_back(dlast->Histo2D({ Form("jac%d_%d",j,i), hname.c_str(), nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", Form("jac%d_%d",j,i)));
	}
      }
      else{
	histosJac.emplace_back(dlast->Histo2D({ Form("jac_%d",i), hname.c_str(), nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "y", "x", Form("jac_%d",i)));
      }
    }

    /*
    histos1D.emplace_back(dlast->Histo1D({"w_pdfx",    "", 20, 0.0, max_x}, "x", "w"));      
    histos1D.emplace_back(dlast->Histo1D({"w_pdfy",    "", 20, 0.0, max_y}, "y", "w"));      
    histos2D.emplace_back(dlast->Histo2D({"w_corrxy",  "", 20, 0.0, max_y, 20, 0.0, max_x}, "y", "x", "w"));      
    histos1D.emplace_back(dlast->Histo1D({"wMC_pdfx",  "", 20, 0.0, max_x}, "x", "wMC"));      
    histos1D.emplace_back(dlast->Histo1D({"wMC_pdfy",  "", 20, 0.0, max_y}, "y", "wMC"));      
    histos2D.emplace_back(dlast->Histo2D({"wMC_corrxy","", 20, 0.0, max_y, 20, 0.0, max_x}, "y", "x", "wMC"));      
    */

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
