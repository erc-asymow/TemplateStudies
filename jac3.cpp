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
constexpr int NMAX  = 200;
constexpr int NMASS = 50;

enum pdf_type { pdf_x=0, pdf_y, corr_x, corr_y, unknown};

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
	("tag", value<std::string>()->default_value(""), "tag name")
	("run", value<std::string>()->default_value("closure"), "run type")
	("do_absy",    bool_switch()->default_value(false), "polycheb in abs(y)")
	("getbin_extTH2_corr",  bool_switch()->default_value(false), "get bin corr(x,y)")
	("inter_extTH2_corr",  bool_switch()->default_value(false), "interpol corr(x,y)")
	("toyTF2_corr",  bool_switch()->default_value(false), "toy TF2 corr(x,y)")
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
      if (vm.count("degs_corr_x")) std::cout << "Degree in x of corrxy: " << vm["degs_corr_x"].as<int>() << '\n';
      if (vm.count("degs_corr_y")) std::cout << "Degree in y of corrxy: " << vm["degs_corr_y"].as<int>() << '\n';
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
  bool do_absy = vm["do_absy"].as<bool>();
  bool getbin_extTH2_corr= vm["getbin_extTH2_corr"].as<bool>();
  bool inter_extTH2_corr= vm["inter_extTH2_corr"].as<bool>();
  bool toyTF2_corr= vm["toyTF2_corr"].as<bool>();
  bool fit_qt_y = vm["fit_qt_y"].as<bool>();

  if(vm.count("degs_corr_x")) tag += std::string(Form("_UL_%d", degs_corr_x));
  if(vm.count("degs_corr_y")) tag += std::string(Form("_%d", degs_corr_y));

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

  auto degs = [degs_corr_x,degs_corr_y]
    (const pdf_type& pdf){
    switch(pdf){
    case pdf_type::corr_x:            // corr(x,.) 
    if(degs_corr_x>0) return degs_corr_x;
    else return 2; 
    case pdf_type::corr_y:            // corr(.,y)
    if(degs_corr_y>0) return degs_corr_y;
    else return 2;
    default: return 1;
    }
  };

  unsigned int njacs = 0;
  unsigned int first_jac_corrxy = njacs;  
  njacs += degs(pdf_type::corr_x) * degs(pdf_type::corr_y) * 6;

  auto toy_mass = [](double Q, double M, double G){
    return 1./TMath::Pi()/(1 + (Q-M)*(Q-M)/G/G);
  };

  TF1* toy_x = new TF1("toy_x", "[0]*x/(x*x+[1])", 0.0, max_x);  
  toy_x->SetParameter(0, 1.0);
  toy_x->SetParameter(1, +2.35e-03);
  double int_toy_x = toy_x->Integral(0.0, max_x);
  toy_x->SetParameter(0, 1.0/int_toy_x);

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

  TF2* toy_A0 = new TF2("toy_A0", "2*y*y*(1 - 0.01*x*x)", 0., max_y, 0., max_x);
  TF2* toy_A1 = new TF2("toy_A1", "(0.5*y + 2*y*y)*( 0.05*x + 0.01*x*x)", 0., max_y, 0., max_x);
  TF2* toy_A2 = new TF2("toy_A2", "2*y*y*(1 - 0.01*x*x)", 0., max_y, 0., max_x);
  TF2* toy_A3 = new TF2("toy_A3", "(y + y*y + y*y*y)*(1 - 0.01*x*x)", 0., max_y, 0., max_x);
  TF2* toy_A4 = new TF2("toy_A4", "(y+1)*(0.5*x + x*x)/6.0", 0., max_y, 0., max_x);

  // We prepare an input tree to run on
  TFile* fout = TFile::Open(("root/histos_"+tag+"_"+run+".root").c_str(), "RECREATE");

  ROOT::RDataFrame d(nevents);

  unsigned int nslots = d.GetNSlots();
  std::vector<TRandom3*> rans = {};
  for(unsigned int i = 0; i < nslots; i++){
    rans.emplace_back( new TRandom3(4357 + i*10) );
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
						  out.emplace_back( toy_mass(Q,MW+0.010,GW)/gen );
						  out.emplace_back( toy_mass(Q,MW-0.010,GW)/gen );
						  for(unsigned int i=0; i<NMASS; i++)
						    out.emplace_back( toy_mass(Q, MW - 0.100 + 0.200/NMASS*i,GW)/gen );
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
    
  dlast = std::make_unique<RNode>(dlast->Define("pt",  [](RVecD p4lab){ return p4lab[0];}, {"p4lab"}));
  if(xLow>=0.)
    dlast = std::make_unique<RNode>(dlast->Define("eta", [](RVecD p4lab){ return TMath::Abs(p4lab[1]);}, {"p4lab"}));
  else
    dlast = std::make_unique<RNode>(dlast->Define("eta", [](RVecD p4lab){ return p4lab[1];}, {"p4lab"}));

  std::vector<ROOT::RDF::RResultPtr<TH1D> > histos1D;
  std::vector<ROOT::RDF::RResultPtr<TH2D> > histos2D;
  std::vector<ROOT::RDF::RResultPtr<TH2D> > histosJac;

  int poi_cat[NMAX];
  unsigned int poi_idx[NMAX];
  unsigned int poi_counter = 0;
  unsigned int n_pdfx = degs(pdf_type::corr_x);
  unsigned int n_pdfy = degs(pdf_type::corr_y);
  double points_x[ NMAX ];
  double points_y[ NMAX ];

  points_x[0] = 0.0;
  for(unsigned int i = 1; i<=n_pdfx; i++){
    points_x[i] = max_x/n_pdfx*i;
    cout << points_x[i] << endl;
  }
  points_y[0] = 0.0;
  for(unsigned int j = 1; j<=n_pdfy; j++){
    points_y[j] = max_y/n_pdfy*j;
    cout << points_y[j] << endl;
  }
 
  // unpol
  for(int k = 0; k<n_pdfx; k++){
    for(int l = 0; l<n_pdfy; l++){
      poi_cat[poi_counter] = -1;
      poi_idx[poi_counter] = n_pdfx*k + l;
      poi_counter++;
    }
  }
  // A0
  for(int k = 0; k<n_pdfx; k++){
    for(int l = 0; l<n_pdfy; l++){
      poi_cat[poi_counter] = 0;
      poi_idx[poi_counter] = n_pdfx*k + l;
      poi_counter++;
    }
  }
  // A1
  for(int k = 0; k<n_pdfx; k++){
    for(int l = 0; l<n_pdfy; l++){
      poi_cat[poi_counter] = 1;
      poi_idx[poi_counter] = n_pdfx*k + l;	
      poi_counter++;
    }
  }
  // A2
  for(int k = 0; k<n_pdfx; k++){
    for(int l = 0; l<n_pdfy; l++){
      poi_cat[poi_counter] = 2;
      poi_idx[poi_counter] = n_pdfx*k + l;
      poi_counter++;
    }
  }
  // A3
  for(int k = 0; k<n_pdfx; k++){
    for(int l = 0; l<n_pdfy; l++){
      poi_cat[poi_counter] = 3;
      poi_idx[poi_counter] = n_pdfx*k + l;
      poi_counter++;
    }
  }
  // A4
  for(int k = 0; k<n_pdfx; k++){
    for(int l = 0; l<n_pdfy; l++){
      poi_cat[poi_counter] = 4;
      poi_idx[poi_counter] = n_pdfx*k + l;
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

  dlast = std::make_unique<RNode>(dlast->Define("weightsMC_angular", [&](double x, double y, RVecD harmonics)->RVecD {
	RVecD out;
	double wMC{1.0};
	wMC *= (harmonics.at(0) + 
		toy_A0->Eval(TMath::Abs(y),x)*harmonics.at(1) +
		toy_A1->Eval(TMath::Abs(y),x)*harmonics.at(2) + 
		toy_A2->Eval(TMath::Abs(y),x)*harmonics.at(3) +
		toy_A3->Eval(TMath::Abs(y),x)*harmonics.at(4) +
		toy_A4->Eval(TMath::Abs(y),x)*harmonics.at(5)
		);
	out.emplace_back( wMC );
	return out;
      }, {"x", "y", "harmonics"} ));

  dlast = std::make_unique<RNode>(dlast->Define("weights_jac", [&,points_x,points_y]
						(double x, double y, 
						 RVecD harmonics,
						 RVecD weightsMC,
						 RVecD weightsMC_angular)->RVecD {
						  RVecD out;						
						  double absy = TMath::Abs(y);
						  // UL
						  for(unsigned int i = 0; i<degs(pdf_type::corr_x); i++){
						    for(unsigned int j = 0; j<degs(pdf_type::corr_y); j++){
						      bool is_in_bin = 
							(x>=points_x[i] && x<points_x[i+1]) && 
							(absy>=points_y[j] && absy<points_y[j+1]);
						      out.emplace_back( is_in_bin ? weightsMC.at(0)*harmonics.at(0)/weightsMC_angular.at(0) : 0.0);
						    }
						  }
						  // A0
						  for(unsigned int i = 0; i<degs(pdf_type::corr_x); i++){
						    for(unsigned int j = 0; j<degs(pdf_type::corr_y); j++){
						      bool is_in_bin = 
							(x>=points_x[i] && x<points_x[i+1]) && 
							(absy>=points_y[j] && absy<points_y[j+1]);
						      out.emplace_back( is_in_bin ? weightsMC.at(0)*harmonics.at(1)/weightsMC_angular.at(0) : 0.0);
						    }
						  }
						  // A1
						  for(unsigned int i = 0; i<degs(pdf_type::corr_x); i++){
						    for(unsigned int j = 0; j<degs(pdf_type::corr_y); j++){
						      bool is_in_bin = 
							(x>=points_x[i] && x<points_x[i+1]) && 
							(absy>=points_y[j] && absy<points_y[j+1]);
						      out.emplace_back( is_in_bin ? weightsMC.at(0)*harmonics.at(2)/weightsMC_angular.at(0) : 0.0);
						    }
						  }
						  // A2
						  for(unsigned int i = 0; i<degs(pdf_type::corr_x); i++){
						    for(unsigned int j = 0; j<degs(pdf_type::corr_y); j++){
						      bool is_in_bin = 
							(x>=points_x[i] && x<points_x[i+1]) && 
							(absy>=points_y[j] && absy<points_y[j+1]);
						      out.emplace_back( is_in_bin ? weightsMC.at(0)*harmonics.at(3)/weightsMC_angular.at(0) : 0.0);
						    }
						  }
						  // A3
						  for(unsigned int i = 0; i<degs(pdf_type::corr_x); i++){
						    for(unsigned int j = 0; j<degs(pdf_type::corr_y); j++){
						      bool is_in_bin = 
							(x>=points_x[i] && x<points_x[i+1]) && 
							(absy>=points_y[j] && absy<points_y[j+1]);
						      out.emplace_back( is_in_bin ? weightsMC.at(0)*harmonics.at(4)/weightsMC_angular.at(0) : 0.0);
						    }
						  }
						  // A4
						  for(unsigned int i = 0; i<degs(pdf_type::corr_x); i++){
						    for(unsigned int j = 0; j<degs(pdf_type::corr_y); j++){
						      bool is_in_bin = 
							(x>=points_x[i] && x<points_x[i+1]) && 
							(absy>=points_y[j] && absy<points_y[j+1]);
						      out.emplace_back( is_in_bin ? weightsMC.at(0)*harmonics.at(5)/weightsMC_angular.at(0) : 0.0);
						    }
						  }
						  cout << out;
						  return out;
						}, {"x", "y",
						    "harmonics",
						    "weightsMC", "weightsMC_angular"} ));

    dlast = std::make_unique<RNode>(dlast->Define("wMC",      [](RVecD weights){ return weights.at(0);}, {"weightsMC"} ));
    dlast = std::make_unique<RNode>(dlast->Define("wMC_up",   [](RVecD weights){ return weights.at(1);}, {"weightsMC"} ));
    dlast = std::make_unique<RNode>(dlast->Define("wMC_down", [](RVecD weights){ return weights.at(2);}, {"weightsMC"} ));
    for(unsigned int i=0; i<NMASS; i++){
      dlast = std::make_unique<RNode>(dlast->Define(Form("wMC_mass%d", i), [i](RVecD weights){ return weights.at(3+i);}, {"weightsMC"} ));
    }
    for(unsigned int i = 0; i < njacs; i++){
      dlast = std::make_unique<RNode>(dlast->Define(Form("jac_%d", i), [i](RVecD weights){ return weights.at(i);}, {"weights_jac"} ));
    }
    
    sums.emplace_back( dlast->Sum<double>("wMC") );

    if(!fit_qt_y){
      histos2D.emplace_back(dlast->Histo2D({"h",       "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", "wMC"));      
      histos2D.emplace_back(dlast->Histo2D({"hMC",     "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", "wMC"));      
      histos2D.emplace_back(dlast->Histo2D({"hMC_up",  "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", "wMC_up"));      
      histos2D.emplace_back(dlast->Histo2D({"hMC_down","", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", "wMC_down"));
      for(unsigned int i=0; i<NMASS; i++){
	histos2D.emplace_back(dlast->Histo2D({Form("hMC_mass%d", i),"", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", Form("wMC_mass%d",i)));
      }
    }
    else{
      histos2D.emplace_back(dlast->Histo2D({"h",       "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "y", "x", "wMC"));      
      histos2D.emplace_back(dlast->Histo2D({"hMC",     "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "y", "x", "wMC"));      
      histos2D.emplace_back(dlast->Histo2D({"hMC_up",  "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "y", "x", "wMC_up"));      
      histos2D.emplace_back(dlast->Histo2D({"hMC_down","", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "y", "x", "wMC_down"));
    }
    
    int njacs_base = njacs/6;
    for(unsigned int i = 0; i < njacs; i++){
      std::string hname = "";
      if(i>=0 && i<njacs_base){
	unsigned int idx_x = (i-0) / degs(pdf_type::corr_y);
	unsigned int idx_y = (i-0) % degs(pdf_type::corr_y);
	hname = std::string(Form("jac_%d: d(pdf) / d(corrxy_in[%d][%d])", i, idx_x, idx_y));	
      }
      else if(i>=njacs_base && i<2*njacs_base){
	unsigned int idx_x = (i-njacs_base) / degs(pdf_type::corr_y);
	unsigned int idx_y = (i-njacs_base) % degs(pdf_type::corr_y);
	hname = std::string(Form("jac_%d: d(pdf) / d(A0xy_in[%d][%d])", i, idx_x, idx_y));	
      }
      else if(i>=2*njacs_base && i<3*njacs_base){
	unsigned int idx_x = (i-2*njacs_base) / degs(pdf_type::corr_y);
	unsigned int idx_y = (i-2*njacs_base) % degs(pdf_type::corr_y);
	hname = std::string(Form("jac_%d: d(pdf) / d(A1xy_in[%d][%d])", i, idx_x, idx_y));	
      }
      else if(i>=3*njacs_base && i<4*njacs_base){
	unsigned int idx_x = (i-3*njacs_base) / degs(pdf_type::corr_y);
	unsigned int idx_y = (i-3*njacs_base) % degs(pdf_type::corr_y);
	hname = std::string(Form("jac_%d: d(pdf) / d(A2xy_in[%d][%d])", i, idx_x, idx_y));	
      }
      else if(i>=4*njacs_base && i<5*njacs_base){
	unsigned int idx_x = (i-4*njacs_base) / degs(pdf_type::corr_y);
	unsigned int idx_y = (i-4*njacs_base) % degs(pdf_type::corr_y);
	hname = std::string(Form("jac_%d: d(pdf) / d(A3xy_in[%d][%d])", i, idx_x, idx_y));	
      }
      else if(i>=5*njacs_base && i<6*njacs_base){
	unsigned int idx_x = (i-5*njacs_base) / degs(pdf_type::corr_y);
	unsigned int idx_y = (i-5*njacs_base) % degs(pdf_type::corr_y);
	hname = std::string(Form("jac_%d: d(pdf) / d(A4xy_in[%d][%d])", i, idx_x, idx_y));	
      }

      if(!fit_qt_y){
	histosJac.emplace_back(dlast->Histo2D({ Form("jac_%d",i), hname.c_str(), nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", Form("jac_%d",i)));
      }
      else{
	histosJac.emplace_back(dlast->Histo2D({ Form("jac_%d",i), hname.c_str(), nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "y", "x", Form("jac_%d",i)));
      }
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
    outtree->Branch("poi_cat", &poi_cat, "poi_cat[poi_counter]/I");
    outtree->Branch("poi_idx", &poi_idx, "poi_idx[poi_counter]/i");
    outtree->Branch("nevents", &nevents, "nevents/L");

    outtree->Fill();
    outtree->Write();
    
    fout->Close();
    
    for(auto r : rans) delete r;
    
    return 1;
}
