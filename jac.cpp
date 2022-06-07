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

//#include <Eigen/Core>
//using Eigen::MatrixXd;

using namespace std;
using namespace ROOT;
typedef ROOT::VecOps::RVec<double> RVecD;
using ROOT::RDF::RNode; 

using namespace boost::program_options;

//std::vector<std::string> helicities = {"UL", "A0", "A1", "A4"};

std::vector<std::string> helicities = {"UL", "A0"};

constexpr double MW = 80.;
constexpr double GW = 2.0;
constexpr int NMAX  = 100;

enum pdf_type { pdf_x=0, pdf_y, corr_x, corr_y, A0_x, A0_y, A1_x, A1_y, A4_x, A4_y, unknown};

auto get_pdf_type = [](const std::string& name) {
  if     (name=="A0_x") return pdf_type::A0_x; 
  else if(name=="A0_y") return pdf_type::A0_y;
  if     (name=="A1_x") return pdf_type::A1_x; 
  else if(name=="A1_y") return pdf_type::A1_y;
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
	("degs_pdf_x",  value<int>()->default_value(5), "max degree of pdf_x")
	("degs_pdf_y" , value<int>()->default_value(5), "max degree of pdf_y")
	("degs_corr_x", value<int>()->default_value(2), "max degree in x of corrxy")
	("degs_corr_y", value<int>()->default_value(2), "max degree in y of corrxy")
	("degs_A0_x",   value<int>()->default_value(2), "max degree in x for A0")
	("degs_A0_y",   value<int>()->default_value(2), "max degree in y for A0")
	("tag", value<std::string>()->default_value(""), "tag name")
	("run", value<std::string>()->default_value("closure"), "run type")
	("do_absy",    bool_switch()->default_value(false), "polycheb in abs(y)")
	("flat_corr",  bool_switch()->default_value(false), "flat corr(x,y)")
	("getbin_extTH2_corr",  bool_switch()->default_value(false), "get bin corr(x,y)")
	("inter_extTH2_corr",  bool_switch()->default_value(false), "interpol corr(x,y)")
	("toyTF2_corr",  bool_switch()->default_value(false), "toy TF2 corr(x,y)")
	("normalize_pdfx",  bool_switch()->default_value(false), "normalize pdfx interpolant")
	("normalize_pdfy",  bool_switch()->default_value(false), "normalize pdfy interpolant");

      store(parse_command_line(argc, argv, desc), vm);
      notify(vm);
      if (vm.count("help")){
	std::cout << desc << '\n';
	return 0;
      }
      if (vm.count("nevents"))    std::cout << "Number of events: " << vm["nevents"].as<long>() << '\n';
      if (vm.count("tag"))        std::cout << "Tag: " << vm["tag"].as<std::string>() << '\n';
      if (vm.count("run"))        std::cout << "Run: " << vm["run"].as<std::string>() << '\n';
      if (vm.count("degs_pdf_x")) std::cout << "Degree of pdf_x: " << vm["degs_pdf_x"].as<int>() << '\n';
      if (vm.count("degs_pdf_y")) std::cout << "Degree of pdf_y: " << vm["degs_pdf_y"].as<int>() << '\n';
      if (vm.count("degs_corr_y")) std::cout << "Degree in x of corrxy: " << vm["degs_corr_x"].as<int>() << '\n';
      if (vm.count("degs_corr_y")) std::cout << "Degree in y of corrxy: " << vm["degs_corr_y"].as<int>() << '\n';
      if (vm.count("degs_A0_x")) std::cout << "Degree in x of A0: " << vm["degs_A0_x"].as<int>() << '\n';
      if (vm.count("degs_A0_y")) std::cout << "Degree in y of A0: " << vm["degs_A0_y"].as<int>() << '\n';
      if (vm.count("do_absy"))     std::cout << "Do abs(Y): " << vm["do_absy"].as<bool>() << '\n';
    }
  catch (const error &ex)
    {
      std::cerr << ex.what() << '\n';
    }

  std::vector<double> norms_cheb4  = {0.0332073, 0.266696, 0.400194, 0.266696, 0.0332073};
  std::vector<double> norms_cheb5  = {0.0199529, 0.180419, 0.299628, 0.299628, 0.180419, 0.0199529};
  std::vector<double> norms_cheb6  = {0.0141919, 0.126953, 0.228635, 0.26044, 0.228635, 0.126953, 0.0141919};
  std::vector<double> norms_cheb7  = {0.0100919, 0.0951067, 0.1761, 0.218702, 0.218702, 0.1761, 0.0951067, 0.0100919};
  std::vector<double> norms_cheb8  = {0.00792723, 0.0731216, 0.139826, 0.180815, 0.196621, 0.180815, 0.139826, 0.0731216, 0.00792723};
  std::vector<double> norms_cheb9  = {0.00613895, 0.0583999, 0.112778, 0.15088, 0.171803, 0.171803, 0.15088, 0.112778, 0.0583999, 0.00613895};
  std::vector<double> norms_cheb10 = {0.00502034, 0.0472209, 0.0927546, 0.127049, 0.149529, 0.156853, 0.149529, 0.127049, 0.0927546, 0.0472209, 0.00502034};
  std::vector<double> norms_cheb11 = {0.00411846, 0.0392313, 0.0775334, 0.10785, 0.129849, 0.141418, 0.141418, 0.129849, 0.10785, 0.0775334, 0.0392313, 0.00411846};
  std::vector<double> norms_cheb12 = {0.00348717, 0.0329594, 0.065827, 0.0923741, 0.113462, 0.126353, 0.131075, 0.126353, 0.113462, 0.0923741, 0.065827, 0.0329594, 0.00348717};
  std::vector<double> norms_cheb13 = {0.00297783, 0.0281563, 0.0564534, 0.0799701, 0.099526, 0.112965, 0.119952, 0.119952, 0.112965, 0.099526, 0.0799701, 0.0564534, 0.0281563, 0.00297783};
  std::vector<double> norms_cheb14 = {0.00252998, 0.02428, 0.0489618, 0.0697433, 0.0878415, 0.100956, 0.109496, 0.112383, 0.109496, 0.100956, 0.0878415, 0.0697433, 0.0489618, 0.02428, 0.00252998};

  auto nu_chebs = [norms_cheb4,norms_cheb5,norms_cheb6,norms_cheb7,norms_cheb8,norms_cheb9,norms_cheb10,norms_cheb11,norms_cheb12,norms_cheb13,norms_cheb14]
    (const int& n, const int& m){
    switch(n){
    case 4: return norms_cheb4[m];
    case 5: return norms_cheb5[m];
    case 6: return norms_cheb6[m];
    case 7: return norms_cheb7[m];
    case 8: return norms_cheb8[m];
    case 9: return norms_cheb9[m];
    case 10: return norms_cheb10[m];
    case 11: return norms_cheb11[m];
    case 12: return norms_cheb12[m];
    case 13: return norms_cheb13[m];
    case 14: return norms_cheb14[m];
    default: return 1.0;
    }
  };

  long nevents    = vm["nevents"].as<long>();
  std::string tag = vm["tag"].as<std::string>();
  std::string run = vm["run"].as<std::string>();
  int degs_pdf_x  = vm["degs_pdf_x"].as<int>();
  int degs_pdf_y  = vm["degs_pdf_y"].as<int>();
  int degs_corr_x = vm["degs_corr_x"].as<int>();
  int degs_corr_y = vm["degs_corr_y"].as<int>();
  int degs_A0_x   = vm["degs_A0_x"].as<int>();
  int degs_A0_y   = vm["degs_A0_y"].as<int>();
  bool do_absy = vm["do_absy"].as<bool>();
  bool flat_corr = vm["flat_corr"].as<bool>();
  bool getbin_extTH2_corr= vm["getbin_extTH2_corr"].as<bool>();
  bool inter_extTH2_corr= vm["inter_extTH2_corr"].as<bool>();
  bool toyTF2_corr= vm["toyTF2_corr"].as<bool>();
  bool normalize_pdfx = vm["normalize_pdfx"].as<bool>();
  bool normalize_pdfy = vm["normalize_pdfy"].as<bool>();

  if(vm.count("degs_pdf_x"))  tag += std::string(Form("_%d", degs_pdf_x));
  if(vm.count("degs_pdf_y"))  tag += std::string(Form("_%d", degs_pdf_y));
  if(vm.count("degs_corr_x")) tag += std::string(Form("_%d", degs_corr_x));
  if(vm.count("degs_corr_y")) tag += std::string(Form("_%d", degs_corr_y));
  if(vm.count("degs_A0_x"))   tag += std::string(Form("_%d", degs_A0_x));
  if(vm.count("degs_A0_y"))   tag += std::string(Form("_%d", degs_A0_y));

  const double max_x = 0.4;
  const double max_y = 3.0;
  const int nbinsX   = 12; 
  const double xLow  = 0.0;
  const double xHigh = +2.5;
  const int nbinsY   = 15; 
  const double yLow  = 25.;
  const double yHigh = 55.;

  auto degs = [degs_pdf_x,degs_pdf_y,degs_corr_x,degs_corr_y,degs_A0_x,degs_A0_y]
    (const pdf_type& pdf){
    switch(pdf){
    case pdf_type::pdf_x:             // f(x)
    if(degs_pdf_x>0) return degs_pdf_x;
    else return 2;
    case pdf_type::pdf_y:             // f(y|0)  
    if(degs_pdf_y>0) return degs_pdf_y;
    else return 2; 
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
    case pdf_type::A1_x:   return  2; // A1(x,.)
    case pdf_type::A1_y:   return  2; // A1(.,y)
    case pdf_type::A4_x:   return  2; // A4(x,.)
    case pdf_type::A4_y:   return  2; // A4(.,y)
    default: return 1;
    }
  };

  unsigned int njacs = 0;
  unsigned int first_jac_pdfx = njacs;
  njacs += (degs(pdf_type::pdf_x) + 1 - 1 - int(normalize_pdfx));
  unsigned int first_jac_pdfy = njacs;
  njacs += (do_absy ? degs(pdf_type::pdf_y) : degs(pdf_type::pdf_y)/2 ) + 1 - int(normalize_pdfy);
  unsigned int first_jac_corrxy = njacs;  
  njacs += (degs(pdf_type::corr_x) + 1 - 1)*( (do_absy ? degs(pdf_type::corr_y) : degs(pdf_type::corr_y)/2) + 1 );

  //njacs += (degs(pdf_type::pdf_y) - int(normalize_pdfy) );
  //njacs += (degs(pdf_type::corr_x) - 1)*degs(pdf_type::corr_y);
  //njacs += degs(pdf_type::A0_x)*degs(pdf_type::A0_y);

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
  //TF1* toy_y = new TF1("toy_y", "[0]+[1]", -max_y, max_y);
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

  TF2* toy_A0 = new TF2("toy_A0", "2*y*y*(1/(1 + 0.02*x*x))", 0., max_y, 0., max_x);
		       
  // preprare inputs
  if(true){

    TFile* fout = TFile::Open(("root/input_"+tag+".root").c_str(), "RECREATE");
    TTree* tree = new TTree("tree", "tree");

    double pdf_x[NMAX];
    for(int i = 0; i<=degs(pdf_type::pdf_x); i++){
      tree->Branch(Form("pdfx_%d", i), &(pdf_x[i]), Form("pdfx_%d/D", i));
      //pdf_x[i] = 1.0/max_x;
      pdf_x[i] = toy_x->Eval( (TMath::Cos((degs(pdf_type::pdf_x)-i)*TMath::Pi()/degs(pdf_type::pdf_x))+1.0)*0.5*max_x );
    }
     
    double pdf_y[NMAX];
    for(int j = 0; j<=degs(pdf_type::pdf_y); j++){
      tree->Branch(Form("pdfy_%d", j), &(pdf_y[j]), Form("pdfy_%d/D", j));
      if(do_absy) 
	pdf_y[j] = toy_y->Eval( (TMath::Cos((degs(pdf_type::pdf_y)-j)*TMath::Pi()/degs(pdf_type::pdf_y))+1.0)*0.5*max_y );
      else
	pdf_y[j] = toy_y->Eval( TMath::Cos((degs(pdf_type::pdf_y)-j)*TMath::Pi()/degs(pdf_type::pdf_y))*max_y );
    }
    
  double corr_xy[NMAX*NMAX];
  for(int k = 0; k<=degs(pdf_type::corr_x); k++){
    for(int l = 0; l<=degs(pdf_type::corr_y); l++){      
      int idx = (degs(pdf_type::corr_y)+1)*k + l; 
      tree->Branch(Form("corrxy_%d_%d", k,l), &(corr_xy[idx]), Form("corrxy_%d_%d/D", k,l));
      corr_xy[idx] = 1.0;

      // Set corr_xy(0, y) to 1.0;
      if(k==0) continue;

      double x = (TMath::Cos((degs(pdf_type::corr_x)-k)*TMath::Pi()/degs(pdf_type::corr_x))+1.0)*0.5*max_x;
      double y = do_absy ? 
	(TMath::Cos((degs(pdf_type::corr_y)-l)*TMath::Pi()/degs(pdf_type::corr_y))+1.0)*0.5*max_y : 
	TMath::Cos((degs(pdf_type::corr_y)-l)*TMath::Pi()/degs(pdf_type::corr_y))*max_y;
      if(flat_corr) corr_xy[idx] = 1.0;
      else if(getbin_extTH2_corr) corr_xy[idx] = th2_corrxy->GetBinContent( th2_corrxy->FindBin(TMath::Abs(y), x) );
      else if(inter_extTH2_corr)  corr_xy[idx] = th2_corrxy->Interpolate(TMath::Abs(y),x);
      else if(toyTF2_corr) corr_xy[idx] = toy_corrxy->Eval(TMath::Abs(y),x);
    }
  }
  
  double A0_xy[NMAX];
  for(int m = 0; m<=degs(pdf_type::A0_x); m++){
    for(int n = 0; n<=degs(pdf_type::A0_y); n++){
      int idx = (degs(pdf_type::A0_y)+1)*m + n; 
      tree->Branch(Form("A0_xy_%d_%d", m,n), &(A0_xy[idx]), Form("A0_xy_%d_%d/D", m,n));
      double x = (TMath::Cos((degs(pdf_type::A0_x)-m)*TMath::Pi()/degs(pdf_type::A0_x))+1.0)*0.5*max_x;
      double y = do_absy ? (TMath::Cos((degs(pdf_type::A0_y)-n)*TMath::Pi()/degs(pdf_type::A0_y))+1.0)*0.5*max_y : 
	TMath::Cos((degs(pdf_type::A0_y)-n)*TMath::Pi()/degs(pdf_type::A0_y))*max_y;
      A0_xy[idx] = toy_A0->Eval(TMath::Abs(y),x);      
    }
  }

  tree->Fill();
  tree->Write();
  fout->Close();

  if(nevents<0) return 0;
  }


  // We prepare an input tree to run on
  TFile* fout = TFile::Open(("root/histos_"+tag+"_"+run+".root").c_str(), "RECREATE");

  TRandom3 R(1);
  ROOT::RDataFrame d(nevents);
  auto dlast = std::make_unique<RNode>(d);
  std::vector<ROOT::RDF::RResultPtr<double> > sums = {};

  dlast = std::make_unique<RNode>(dlast->Define("Q",  [&](){ 
	double Q = TMath::Tan(R.Uniform(-TMath::Pi()*0.5, +TMath::Pi()*0.5))*GW + MW; 
	while(Q < 50.0){
	  Q = TMath::Tan(R.Uniform(-TMath::Pi()*0.5, +TMath::Pi()*0.5))*GW + MW; 
	}
	return Q;
      } ));
  dlast = std::make_unique<RNode>(dlast->Define("cos",[&](){ return R.Uniform(-1.0, 1.0);} ));
  dlast = std::make_unique<RNode>(dlast->Define("phi",[&](){ return R.Uniform(-TMath::Pi(), +TMath::Pi());} ));
  dlast = std::make_unique<RNode>(dlast->Define("x",  [&](){ return R.Uniform(0.0, max_x); }));
  dlast = std::make_unique<RNode>(dlast->Define("y",  [&](){ return R.Uniform(-max_y, max_y); }));
  dlast = std::make_unique<RNode>(dlast->Define("weightsM", 
						[&](double Q)->RVecD{
						  RVecD out;
						  double gen = 1./TMath::Pi()/(1 + (Q-MW)*(Q-MW)/GW/GW);
						  out.emplace_back( toy_mass(Q,MW,GW)/gen );
						  out.emplace_back( toy_mass(Q,MW+0.010,GW)/gen );
						  out.emplace_back( toy_mass(Q,MW-0.010,GW)/gen );
						  return out; 
						}, {"Q"}));
  dlast = std::make_unique<RNode>(dlast->Define("wM",     [](RVecD weights){ return weights.at(0);}, {"weightsM"} ));
  dlast = std::make_unique<RNode>(dlast->Define("wM_up",  [](RVecD weights){ return weights.at(1);}, {"weightsM"} ));
  dlast = std::make_unique<RNode>(dlast->Define("wM_down",[](RVecD weights){ return weights.at(2);}, {"weightsM"} ));
  dlast = std::make_unique<RNode>(dlast->Define("p4lab",
						[&](double Q, double cos, double phi, double x, double y)->RVecD {
						  double qT2 = x*x*Q*Q;
						  double qT = TMath::Sqrt(qT2);
						  double XT = TMath::Sqrt(qT2 + Q*Q);
						  double Q0 = XT*TMath::CosH(y);
						  double Q3 = XT*TMath::SinH(y);
						  TMatrixD A(4,4);
						  double row0[4] = {Q0/Q,       -qT/Q, 0., -Q3/Q};
						  double row1[4] = {-qT*Q0/Q/XT,  XT/Q, 0., qT*Q3/Q/XT};
						  double row2[4] = {0., 0., 1.0, 0.};
						  double row3[4] = {-Q3/XT, 0., 0., Q0/XT};
						  A.InsertRow(0, 0, row0);
						  A.InsertRow(1, 0, row1);
						  A.InsertRow(2, 0, row2);
						  A.InsertRow(3, 0, row3);
						  double cosS = y>0 ? cos : -cos;
						  double phiS = y>0 ? phi : -phi;
						  double sinS = TMath::Sqrt(1-cosS*cosS);
						  //std::cout << Q/2 << "," << sinS << "," << TMath::Cos(phiS) << std::endl;
						  TVectorD p4_CS(0,3, Q/2, Q/2*sinS*TMath::Cos(phi), Q/2*sinS*TMath::Sin(phi), Q/2*cos, "END");
						  //std::cout << "CS" << std::endl;
						  //p4_CS.Print();
						  TVectorD p4_lab = A.Invert()*p4_CS;
						  //std::cout << "LAB" << std::endl;
						  //p4_lab.Print();
						  double pt  = TMath::Sqrt(p4_lab[1]*p4_lab[1] + p4_lab[2]*p4_lab[2]);
						  int sign = TMath::Abs(p4_lab[3])>0. ? TMath::Abs(p4_lab[3])/p4_lab[3] : 1;
						  double eta = TMath::ACosH(p4_lab[0]/pt)*sign;
						  RVecD p4lab{pt,eta};
						  return p4lab;
						}, {"Q", "cos", "phi", "x", "y"}));
    
  dlast = std::make_unique<RNode>(dlast->Define("pdfx_vec", 
						[&](double x)->RVecD{
						  RVecD out;
						  unsigned int deg = degs(pdf_type::pdf_x); 
						  double cheb_last = cheb(x, 0.5*max_x, 1.0, deg, deg);
						  double nu_last   = nu_chebs(deg, deg);
						  for(unsigned int i = 0; i<=deg; i++){
						    double cheb_i = cheb(x, 0.5*max_x, 1.0, deg, i);
						    if(!normalize_pdfx){
						      out.emplace_back( cheb_i );
						      continue;
						    }
						    double nu_i = nu_chebs(deg, i);
						    out.emplace_back(i<deg ? cheb_i - cheb_last*nu_i/nu_last : cheb_last/nu_last/max_x);
						  }
						  return out;
						}, {"x"} ));

  dlast = std::make_unique<RNode>(dlast->Define("pdfy_vec", 
						[&](double y)->RVecD{
						  RVecD out;
						  unsigned int deg = degs(pdf_type::pdf_y); 
						  unsigned int mid_deg = deg/2;
						  if(do_absy){
						    if(!normalize_pdfy){
						      for(unsigned int j = 0; j<=deg; j++){
							double cheb_j = cheb(TMath::Abs(y), 0.5*max_y, 1.0, deg, j);
							out.emplace_back( cheb_j );
						      }
						    }
						    else{
						      double cheb_last = cheb(TMath::Abs(y), 0.5*max_y, 1.0, deg, deg);
						      double nu_last   = nu_chebs(deg, deg);
						      for(unsigned int j = 0; j<=deg; j++){
							double cheb_j = cheb(TMath::Abs(y), 0.5*max_y, 1.0, deg, j);
							double nu_j = nu_chebs(deg, j);
							out.emplace_back(j<deg ? cheb_j - cheb_last*nu_j/nu_last : cheb_last/nu_last/(2*max_y));
						      }
						    }
						  }
						  else{
						    if(!normalize_pdfy){
						      for(unsigned int j = 0; j<=mid_deg; j++){
							double cheb_j  = cheb(y, max_y, 0.0, deg, j) + cheb(y, max_y, 0.0, deg, deg-j) ;
							double alpha_j = j<mid_deg ? 1.0 : (deg%2==0 ? 0.5 : 1.0);
							out.emplace_back( cheb_j*alpha_j );
						      }
						    }
						    else{
						      double cheb_last = cheb(y, max_y, 0.0, deg, 0) + cheb(y, max_y, 0.0, deg, deg);
						      double nu_last   = nu_chebs(deg, 0);
						      double alpha_last = 0<mid_deg ? 1.0 : (deg%2==0 ? 0.5 : 1.0);
						      for(unsigned int j = 0; j<=mid_deg; j++){
							double cheb_j  = cheb(y, max_y, 0.0, deg, j) + cheb(y, max_y, 0.0, deg, deg-j) ;
							double alpha_j = j<mid_deg ? 1.0 : (deg%2==0 ? 0.5 : 1.0);
							double nu_j = nu_chebs(deg, j);
							out.emplace_back( j==0 ? 
									  0.5*cheb_last/nu_last/alpha_last/(2*max_y) : 
									  (alpha_j/alpha_last)*(cheb_j - cheb_last*nu_j/nu_last));
						      } 
						    }
						  }
						  return out;
						}, {"y"} ));

  dlast = std::make_unique<RNode>(dlast->Define("corrxy_vec", 
						[&](double x, double y)->RVecD{
						  RVecD out;
						  for(unsigned int k = 0; k<=degs(pdf_type::corr_x); k++){
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
						    for(unsigned int k = 0; k<=degs(get_pdf_type(hel+"_x")); k++){
						      double Ax = cheb(x, 0.5*max_x, 1.0, degs(get_pdf_type(hel+"_x")), k);
						      if(do_absy){
							for(unsigned int l = 0; l<=degs(get_pdf_type(hel+"_y")); l++){
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
  dlast = std::make_unique<RNode>(dlast->Define("eta", [](RVecD p4lab){ return TMath::Abs(p4lab[1]);}, {"p4lab"}));
  std::vector<ROOT::RDF::RResultPtr<TH1D> > histos1D;
  //histos1D.emplace_back(dlast->Histo1D({"pt", "", 100, 0, 100}, "pt"));
  
  std::vector<ROOT::RDF::RResultPtr<TH2D> > histos2D;
  std::vector<ROOT::RDF::RResultPtr<TH2D> > histosJac;

  double poi_val[NMAX];
  unsigned int poi_cat[NMAX];
  unsigned int poi_idx[NMAX];
  unsigned int poi_counter = 0;
  
  if(run.find("closure")!=string::npos){

    TFile* fin = TFile::Open(("root/input_"+tag+".root").c_str(), "READ");
    if(fin==nullptr) return 0;
    TTree* tree = fin->Get<TTree>("tree");

    double pdf_x[NMAX];
    for(int i = 0; i<=degs(pdf_type::pdf_x); i++){
      tree->SetBranchAddress(Form("pdfx_%d", i), &(pdf_x[i]));
    }
    double pdf_y[NMAX];
    for(int j = 0; j<=degs(pdf_type::pdf_y); j++){
      tree->SetBranchAddress(Form("pdfy_%d", j), &(pdf_y[j]));
    }
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
    tree->GetEntry(0);    
    fin->Close();
     
    RVecD pdfx_in;
    for(int i = 0; i<=degs(pdf_type::pdf_x); i++){
      pdfx_in.emplace_back( (normalize_pdfx && i==degs(pdf_type::pdf_x)) ? 1.0 : pdf_x[i] );
      if(i==0 || (normalize_pdfx && i==degs(pdf_type::pdf_x)) ) continue;
      poi_val[poi_counter] = pdf_x[i];
      poi_cat[poi_counter] = 0;
      poi_idx[poi_counter] = i;
      poi_counter++;
    }

    RVecD pdfy_in;
    if(do_absy){
      for(int j = 0; j<=degs(pdf_type::pdf_y); j++){ 
	pdfy_in.emplace_back( (normalize_pdfy && j==degs(pdf_type::pdf_y)) ? 1.0 : pdf_y[j] );
	if(normalize_pdfy && j==degs(pdf_type::pdf_y) ) continue;
	poi_val[poi_counter] = pdf_y[j];
	poi_cat[poi_counter] = 1;
	poi_idx[poi_counter] = j;	
	poi_counter++;
      }
    }
    else{
      unsigned int mid_deg = degs(pdf_type::pdf_y)/2;
      for(int j = 0; j<=mid_deg; j++){
	pdfy_in.emplace_back( (normalize_pdfy && j==0) ? 1.0 : pdf_y[j] );
	if( normalize_pdfy && j==0 ) continue;
	poi_val[poi_counter] = pdf_y[j];
	poi_cat[poi_counter] = 1;
	poi_idx[poi_counter] = j;	
	poi_counter++;
      }
    }

    RVecD corrxy_in;
    for(int k = 0; k<=degs(pdf_type::corr_x); k++){
      for(int l = 0; l<=(do_absy ? degs(pdf_type::corr_y) : degs(pdf_type::corr_y)/2); l++){      
	// index to input corr_xy[]
	int idx = (degs(pdf_type::corr_y) + 1)*k + l; 
	corrxy_in.emplace_back( corr_xy[idx] );
	if(k==0) continue;
	poi_val[poi_counter] = corr_xy[idx];
	poi_cat[poi_counter] = 2;
	poi_idx[poi_counter] = idx;	
	poi_counter++;
      }
    }

    // FIX
    RVecD A0xy_in;
    for(int k = 0; k<=degs(pdf_type::A0_x); k++){
      if(do_absy){
	for(int l = 0; l<=degs(pdf_type::A0_y); l++){      
	  int idx = (degs(pdf_type::A0_y)+1)*k + l; 
	  A0xy_in.emplace_back( A0_xy[idx] );
	}
      }
      else{
	unsigned int mid_deg = degs(pdf_type::A0_y)/2;
	for(int l = 0; l<=mid_deg; l++){      
	  int idx = (degs(pdf_type::A0_y)+1)*k + l; 
	  A0xy_in.emplace_back( A0_xy[idx] );
	}
      }
    }
		
    dlast = std::make_unique<RNode>(dlast->Define("harmonics", [&](double x, double y, double cos, double phi) -> RVecD{
	  RVecD out;
	  double cosS = y>0 ? cos : -cos;
	  double phiS = y>0 ? phi : -phi;
	  out.emplace_back( 1.0 + cosS*cosS );
	  out.emplace_back( 0.5*(1-3*cosS*cosS) );
	  /* ... */
	  return out;
	} , {"x", "y", "cos", "phi"} ));

    dlast = std::make_unique<RNode>(dlast->Define("weightsMC", [&](double x, double y, RVecD harmonics, RVecD weightsM  )->RVecD {
	  RVecD out;
	  double wMC{3./16/TMath::Pi()};
	  wMC *= toy_x->Eval(x);
	  wMC *= toy_y->Eval(y);
	  double corrxyMC  = 1.0;
	  if(flat_corr) corrxyMC  = 1.0;
	  else if(getbin_extTH2_corr) corrxyMC = th2_corrxy->GetBinContent( th2_corrxy->FindBin(TMath::Abs(y),x) );
	  else if(inter_extTH2_corr)  corrxyMC = th2_corrxy->Interpolate(TMath::Abs(y),x);
	  else if(toyTF2_corr)        corrxyMC = toy_corrxy->Eval(TMath::Abs(y),x);	   
	  wMC *= corrxyMC;
	  wMC *= (harmonics.at(0) + 
		  toy_A0->Eval(TMath::Abs(y),x)*harmonics.at(1) 
		  /* + ... */ 
		  );
	  out.emplace_back( wMC*weightsM.at(0) );
	  out.emplace_back( wMC*weightsM.at(1) );
	  out.emplace_back( wMC*weightsM.at(2) );
	  return out;
	}, {"x", "y", "harmonics", "weightsM"} ));
    
    dlast = std::make_unique<RNode>(dlast->Define("weights", [& ,pdfx_in,pdfy_in,corrxy_in,A0xy_in ]
						  (RVecD pdfx_vec,//   RVecD pdfx_in,
						   RVecD pdfy_vec,//   RVecD pdfy_in,
						   RVecD corrxy_vec,// RVecD corrxy_in,
						   RVecD A0xy_vec,//   RVecD A0xy_in,
						   RVecD harmonics,
						   RVecD weightsM  )->RVecD {
						    RVecD out;
						    double wUL{3./16/TMath::Pi()};	  
						    wUL *= ROOT::VecOps::Dot(pdfx_vec,pdfx_in);
						    wUL *= ROOT::VecOps::Dot(pdfy_vec,pdfy_in);
						    wUL *= ROOT::VecOps::Dot(corrxy_vec,corrxy_in);
						    double w = wUL*(harmonics.at(0) + 
								    ROOT::VecOps::Dot(A0xy_vec,A0xy_in)*harmonics.at(1) 
								    /* + ... */
								    ); 
						    out.emplace_back( w*weightsM.at(0) );
						    return out;
						  }, {"pdfx_vec", //"pdfx_in", 
						      "pdfy_vec", //"pdfy_in",
						      "corrxy_vec", //"corrxy_in", 
						      "A0xy_vec", //"A0xy_in", 
						      "harmonics",
						      "weightsM"} ));

    dlast = std::make_unique<RNode>(dlast->Define("weights_jac", [&,pdfx_in,pdfy_in,corrxy_in,A0xy_in]
						  (RVecD pdfx_vec,   //RVecD pdfx_in,
						   RVecD pdfy_vec,   //RVecD pdfy_in,
						   RVecD corrxy_vec, //RVecD corrxy_in,
						   RVecD A0xy_vec,   //RVecD A0xy_in,
						   RVecD harmonics,
						   RVecD weightsM  )->RVecD {
						    RVecD out;
						
						    double wUL{3./16/TMath::Pi()};	  
						    double A = ROOT::VecOps::Dot(pdfx_vec,pdfx_in);
						    double B = ROOT::VecOps::Dot(pdfy_vec,pdfy_in);
						    double C = ROOT::VecOps::Dot(corrxy_vec,corrxy_in);
						    double D = (harmonics.at(0) + 
								ROOT::VecOps::Dot(A0xy_vec,A0xy_in)*harmonics.at(1) 
								);

						    // Jacobian pdfx
						    unsigned int njacs_pdfx = (degs(pdf_type::pdf_x) - int(normalize_pdfx)); // +1 -1
						    for(unsigned int i = 0; i<njacs_pdfx; i++){
						      RVecD pdfx_in_copy( pdfx_in.size(), 0.0 );
						      //+1 offet accounts for pdfx[0] constrained to 0.0
						      pdfx_in_copy[i+1] = 1.0;
						      out.emplace_back( wUL*
									ROOT::VecOps::Dot(pdfx_vec,pdfx_in_copy)*
									B*C*D*
									weightsM.at(0) );
						    }

						    // Jacobian pdfy
						    unsigned int njacs_pdfy = (do_absy ? degs(pdf_type::pdf_y) : degs(pdf_type::pdf_y)/2) + 1 - int(normalize_pdfy);
						    for(unsigned int j = 0; j<njacs_pdfy; j++){
						      RVecD pdfy_in_copy( pdfy_in.size(), 0.0 );
						      pdfy_in_copy[(!do_absy && normalize_pdfy ? j+1 : j)] = 1.0;
						      out.emplace_back( wUL*
									A*
									ROOT::VecOps::Dot(pdfy_vec,pdfy_in_copy)*
									C*D*
									weightsM.at(0) );
						    }
						    
						    // Jacobian corrxy
						    // +1 offset accounts for corr[0,y]==1
						    //unsigned int njacs_corrxy = degs(pdf_type::corr_x)*(do_absy ? degs(pdf_type::corr_y)/2+1 : degs(pdf_type::corr_y)+1);
						    for(int k = 1; k<=degs(pdf_type::corr_x); k++){
							for(int l = 0; l<=(do_absy ? degs(pdf_type::corr_y) : degs(pdf_type::corr_y)/2); l++){      
							  int idx = ( (do_absy ? degs(pdf_type::corr_y) : degs(pdf_type::corr_y)/2) +1)*k + l;  
							  RVecD corrxy_in_copy( corrxy_in.size(), 0.0 );
							  corrxy_in_copy[idx] = 1.0;
							  out.emplace_back( wUL*
									    A*B*
									    ROOT::VecOps::Dot(corrxy_vec,corrxy_in_copy)*
									    D*
									    weightsM.at(0) );

							}						      						     
						    }
						    
						    /*
						    for(unsigned int k = 0; k<njacs_corrxy; k++){
						    RVecD corrxy_in_copy( corrxy_in.size(), 0.0 );						    
						    }
						    */

						    return out;

						  }, {"pdfx_vec", //"pdfx_in", 
						      "pdfy_vec", //"pdfy_in",
						      "corrxy_vec", //"corrxy_in", 
						      "A0xy_vec", //"A0xy_in", 
						      "harmonics",
						      "weightsM"} ));

    dlast = std::make_unique<RNode>(dlast->Define("w",        [](RVecD weights){ return weights.at(0);}, {"weights"} ));
    dlast = std::make_unique<RNode>(dlast->Define("wMC",      [](RVecD weights){ return weights.at(0);}, {"weightsMC"} ));
    dlast = std::make_unique<RNode>(dlast->Define("wMC_up",   [](RVecD weights){ return weights.at(1);}, {"weightsMC"} ));
    dlast = std::make_unique<RNode>(dlast->Define("wMC_down", [](RVecD weights){ return weights.at(2);}, {"weightsMC"} ));

    for(unsigned int i = 0; i < njacs; i++){
      dlast = std::make_unique<RNode>(dlast->Define(Form("jac_%d", i), [i](RVecD weights){ return weights.at(i);}, {"weights_jac"} ));
    }
    
    sums.emplace_back( dlast->Sum<double>("w") );
    sums.emplace_back( dlast->Sum<double>("wMC") );

    histos2D.emplace_back(dlast->Histo2D({"w",       "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", "w"));      
    histos2D.emplace_back(dlast->Histo2D({"wMC",     "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", "wMC"));      
    histos2D.emplace_back(dlast->Histo2D({"wMC_up",  "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", "wMC_up"));      
    histos2D.emplace_back(dlast->Histo2D({"wMC_down","", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", "wMC_down"));

    for(unsigned int i = 0; i < njacs; i++){
      std::string hname = "";
      if(i>=first_jac_pdfx && i<first_jac_pdfy) 
	hname = std::string(Form("jac_%d: d(pdf) / d(pdfx_in[%d])", i, i-first_jac_pdfx+1));
      else if(i>=first_jac_pdfy && i<first_jac_corrxy)
	hname = std::string(Form("jac_%d: d(pdf) / d(pdfy_in[%d])", i, i-first_jac_pdfy + 1*(!do_absy && normalize_pdfy) ));
      else if(i>=first_jac_corrxy){
	unsigned int idx_x = (i-first_jac_corrxy) / ((do_absy ? degs(pdf_type::corr_y) : degs(pdf_type::corr_y)/2) +1) + 1;
	unsigned int idx_y = (i-first_jac_corrxy) % ((do_absy ? degs(pdf_type::corr_y) : degs(pdf_type::corr_y)/2) +1);
	hname = std::string(Form("jac_%d: d(pdf) / d(corrxy_in[%d][%d])", i, idx_x, idx_y));	
      }
      histosJac.emplace_back(dlast->Histo2D({ Form("jac_%d",i), hname.c_str(), nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", Form("jac_%d",i)));
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


  //histos2D.emplace_back(d1.Histo2D({"pteta", "", 50, -2.5, 2.5, 40, 20, 60}, "eta", "pt", "pdf_x"));
  auto colNames = dlast->GetColumnNames();
  //for (auto &&colName : colNames) std::cout << colName << std::endl;
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

  int n_pdfx = degs(pdf_type::pdf_x) + 1;
  double norms_pdfx[ degs(pdf_type::pdf_x) + 1];
  double points_x[ degs(pdf_type::pdf_x) + 1];
  outtree->Branch("n_pdfx", &n_pdfx, "n_pdfx/I");
  outtree->Branch("norms_pdfx", &norms_pdfx, "norms_pdfx[n_pdfx]/D");
  outtree->Branch("points_x", &points_x, "points_x[n_pdfx]/D");

  int n_pdfy = degs(pdf_type::pdf_y) + 1;
  double norms_pdfy[ degs(pdf_type::pdf_y) + 1];
  double points_y[ degs(pdf_type::pdf_y) + 1];
  outtree->Branch("n_pdfy", &n_pdfy, "n_pdfy/I");
  outtree->Branch("norms_pdfy", &norms_pdfy, "norms_pdfy[n_pdfy]/D");
  outtree->Branch("points_y", &points_y, "points_y[n_pdfy]/D");

  outtree->Branch("poi_counter", &poi_counter, "poi_counter/i");
  outtree->Branch("poi_val", &poi_val, "poi_val[poi_counter]/D");
  outtree->Branch("poi_cat", &poi_cat, "poi_cat[poi_counter]/i");
  outtree->Branch("poi_idx", &poi_idx, "poi_idx[poi_counter]/i");

  outtree->Branch("nevents", &nevents, "nevents/L");

  for(int i = 0; i<=degs(pdf_type::pdf_x); i++){
    //norms_pdfx[i] = *(sums[i])/total;
    points_x[i]   = (TMath::Cos((degs(pdf_type::pdf_x)-i)*TMath::Pi()/degs(pdf_type::pdf_x))+1.0)*0.5*max_x;
  }
  for(int j = 0; j<=degs(pdf_type::pdf_y); j++){
    //norms_pdfy[j] = *(sums[degs(pdf_type::pdf_x) + 1 + j])/total;
    points_y[j]   = do_absy ? (TMath::Cos((degs(pdf_type::pdf_y)-j)*TMath::Pi()/degs(pdf_type::pdf_y))+1.0)*0.5*max_x : 
      TMath::Cos((degs(pdf_type::pdf_y)-j)*TMath::Pi()/degs(pdf_type::pdf_y))*max_x;
  }
  outtree->Fill();
  outtree->Write();
  
  fout->Close();

  return 1;
}
