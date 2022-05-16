#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include "TRandom3.h"
#include "TVector.h"
#include "TVectorT.h"
#include "TMath.h"
#include "TF1.h"
#include <TMatrixD.h>
#include <ROOT/RVec.hxx>
#include <iostream>
//#include <boost/program_options.hpp>

using namespace std;
using namespace ROOT;
typedef ROOT::VecOps::RVec<double> RVecD;
using ROOT::RDF::RNode; 

//std::vector<std::string> helicities = {"UL", "A0", "A1", "A4"};
std::vector<std::string> helicities = {"UL", "A0"};

constexpr double MW = 80.;
constexpr double GW = 2.0;
constexpr int NMAX = 100;

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


auto cheb = [](double x, double scale, double offset, int n, int m){
  double den = 0.;
  double num = 0.;
  for(unsigned int i = 0; i <= n ; i++){
    int sign = i%2==0 ? +1 :-1;
    double xj = (TMath::Cos(i*TMath::Pi()/n) + offset)*scale;
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

  ROOT::EnableImplicitMT();

  /*
  try
    {
      options_description desc{"Options"};
      desc.add_options()
	("help,h", "Help screen")
	("pi", value<float>()->default_value(3.14f), "Pi")
	("age", value<int>()->notifier(on_age), "Age");
      variables_map vm;
      store(parse_command_line(argc, argv, desc), vm);
      notify(vm);
      if (vm.count("help"))
	std::cout << desc << '\n';
      else if (vm.count("age"))
	std::cout << "Age: " << vm["age"].as<int>() << '\n';
      else if (vm.count("pi"))
	std::cout << "Pi: " << vm["pi"].as<float>() << '\n';
    }
  catch (const error &ex)
    {
      std::cerr << ex.what() << '\n';
    }
  */

  long nevents = 1000;
  std::string tag{""};
  std::string run{""};
  int degs_pdf_x = -1;
  int degs_pdf_y = -1;
  if(argc>1) nevents = strtol(argv[1], NULL, 10);
  if(argc>2) tag = std::string(argv[2]);
  if(argc>3) run = std::string(argv[3]);
  if(argc>4){
    degs_pdf_x = strtol(argv[4], NULL, 10); 
    tag += "_"+std::string(argv[4]);
  }
  if(argc>5){
    degs_pdf_y = strtol(argv[5], NULL, 10);
    tag += "_"+std::string(argv[5]); 
  }

  const double max_x = 0.3;
  const double max_y = 2.5;
  const int nbinsX   = 12; 
  const double xLow  = 0.0;
  const double xHigh = +2.5;
  const int nbinsY   = 15; 
  const double yLow  = 25.;
  const double yHigh = 55.;

  auto degs = [degs_pdf_x,degs_pdf_y](const pdf_type& pdf){
    switch(pdf){
    case pdf_type::pdf_x:             // f(x)
    if(degs_pdf_x>0) return degs_pdf_x;
    else return 2;
    case pdf_type::pdf_y:             // f(y|0)  
    if(degs_pdf_y>0) return degs_pdf_y;
    else return 2; 
    case pdf_type::corr_x: return  2; // P(x,.)
    case pdf_type::corr_y: return  2; // P(.,y)
    case pdf_type::A0_x:   return  2; // A0(x,.)
    case pdf_type::A0_y:   return  2; // A0(.,y)
    case pdf_type::A1_x:   return  2; // A1(x,.)
    case pdf_type::A1_y:   return  2; // A1(.,y)
    case pdf_type::A4_x:   return  2; // A4(x,.)
    case pdf_type::A4_y:   return  2; // A4(.,y)
    default: return 1;
    }
  };
  

  TF1* toy_x = new TF1("toy_x", "[0]*x/(x*x+[1])", 0.0, max_x);  
  //TF1* toy_x = new TF1("toy_x", "[0]+[1]", 0.0, max_x);  
  toy_x->SetParameter(0, 1.0);
  toy_x->SetParameter(1, +2.35e-03);
  TF1* toy_y = new TF1("toy_y", "1/TMath::Sqrt(2*TMath::Pi()*[1])*TMath::Exp(-0.5*(x-[0])*(x-[0])/[1]/[1])", -max_y, max_y);
  //TF1* toy_y = new TF1("toy_y", "[0]+[1]", -max_y, max_y);
  toy_y->SetParameter(0, 0.0);
  toy_y->SetParameter(1, 5.0);

  // preprare inputs
  if(nevents<0){

    TFile* fout = TFile::Open(("root/input_"+tag+".root").c_str(), "RECREATE");
    TTree* tree = new TTree("tree", "tree");

    double pdf_x[NMAX];
    for(int i = 0; i<=degs(pdf_type::pdf_x); i++){
      tree->Branch(Form("pdfx_%d", i), &(pdf_x[i]), Form("pdfx_%d/D", i));
      //pdf_x[i] = 1.0/max_x;
      pdf_x[i] = toy_x->Eval( (TMath::Cos(i*TMath::Pi()/degs(pdf_type::pdf_x))+1.0)*0.5*max_x );
    }
    
    double pdf_y[NMAX];
    for(int j = 0; j<=degs(pdf_type::pdf_y); j++){
      tree->Branch(Form("pdfy_%d", j), &(pdf_y[j]), Form("pdfy_%d/D", j));
      //pdf_y[j] = 1.0/(2*max_y);
      pdf_y[j] = toy_y->Eval( (TMath::Cos(j*TMath::Pi()/degs(pdf_type::pdf_y))+1.0)*0.5*max_y );
    }
    
  double corr_xy[NMAX];
  for(int k = 0; k<=degs(pdf_type::corr_x); k++){
    for(int l = 0; l<=degs(pdf_type::corr_y); l++){      
      int idx = (degs(pdf_type::corr_y)+1)*k + l; 
      tree->Branch(Form("corrxy_%d_%d", k,l), &(corr_xy[idx]), Form("corrxy_%d_%d/D", k,l));
      corr_xy[idx] = 1.0;
    }
  }
  
  double A0_xy[NMAX];
  for(int m = 0; m<=degs(pdf_type::A0_x); m++){
    for(int n = 0; n<=degs(pdf_type::A0_y); n++){
      int idx = (degs(pdf_type::A0_y)+1)*m + n; 
      tree->Branch(Form("A0_xy_%d_%d", m,n), &(A0_xy[idx]), Form("A0_xy_%d_%d/D", m,n));
      A0_xy[idx] = 0.0;      
    }
  }

  tree->Fill();
  tree->Write();
  fout->Close();

  return 0;
  }


  // We prepare an input tree to run on
  TFile* fout = TFile::Open(("root/histos_"+tag+"_"+run+".root").c_str(), "RECREATE");

  //auto fileName = "rdf.root";
  //auto treeName = "tree";
  //auto d = fill_tree(treeName, fileName, nevents, false);
  //std::cout << "Tree filled" << std::endl;
  //ROOT::RDataFrame d(treeName, fileName);
  //auto d1 = d.Define("pt", [](RVecD p4lab){ return p4lab[0];}, {"p4lab"})
  //  .Define("eta", [](RVecD p4lab){ return p4lab[1];}, {"p4lab"});

  TRandom3 R(1);
  ROOT::RDataFrame d(nevents);
  auto dlast = std::make_unique<RNode>(d);

  dlast = std::make_unique<RNode>(dlast->Define("Q",  [&](){ return TMath::Tan(R.Uniform(-TMath::Pi()*0.5, +TMath::Pi()*0.5))*GW + MW; } ));
  dlast = std::make_unique<RNode>(dlast->Define("cos",[&](){ return R.Uniform(-1.0, 1.0);} ));
  dlast = std::make_unique<RNode>(dlast->Define("phi",[&](){ return R.Uniform(-TMath::Pi(), +TMath::Pi());} ));
  dlast = std::make_unique<RNode>(dlast->Define("x",  [&](){ return R.Uniform(0.0, max_x); }));
  dlast = std::make_unique<RNode>(dlast->Define("y",  [&](){ return R.Uniform(-max_y, max_y); }));
  dlast = std::make_unique<RNode>(dlast->Define("wM", 
						[&](double x){
						  double gen = 1./TMath::Pi()/(1 + (x-MW)*(x-MW)/GW/GW);
						  double target = gen;
						  return target/gen; 
						}, {"Q"}));
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
  
  
  std::vector<ROOT::RDF::RResultPtr<double> > sums = {};
  for(int i = 0; i<=degs(pdf_type::pdf_x); i++){
    dlast = std::make_unique<RNode>(dlast->Define(Form("pdfx_%d", i), [i,max_x,degs](double x){return cheb(x, 0.5*max_x, 1.0, degs(pdf_type::pdf_x), i);}, {"x"} ));
    sums.emplace_back( dlast->Sum<double>( Form("pdfx_%d", i)) );
  }

  for(int j = 0; j<=degs(pdf_type::pdf_y); j++){
    dlast = std::make_unique<RNode>(dlast->Define(Form("pdfy_%d", j), [j,max_y,degs](double y){return cheb(TMath::Abs(y), 0.5*max_y, 1.0, degs(pdf_type::pdf_y), j);} , {"y"} ));
    //sums.emplace_back( dlast->Sum( Form("pdfy_%d", j)) );
  }
 
  for(int k = 0; k<=degs(pdf_type::corr_x); k++){
    dlast = std::make_unique<RNode>(dlast->Define(Form("corrxy_%d", k), [k,max_x,degs](double x){return cheb(x, 0.5*max_x, 1.0, degs(pdf_type::corr_x), k); } , {"x"} ));
  }
  for(int l = 0; l<=degs(pdf_type::corr_y); l++){
    dlast = std::make_unique<RNode>(dlast->Define(Form("corryx_%d", l), [l,max_y,degs](double y){return cheb(TMath::Abs(y), 0.5*max_y, 1.0, degs(pdf_type::corr_y), l); }, {"y"} ));
  }

  for(auto name : helicities){

    dlast = std::make_unique<RNode>(dlast->Define(name, [name](double cos, double phi, double y ){
	  double val  = 1.0;
	  double cosS = y>0 ? cos : -cos;
	  double phiS = y>0 ? phi : -phi;
	  if     (name=="UL") val = 1+cosS*cosS;
	  else if(name=="A0") val = 0.5*(1-3*cosS*cosS);
	  else if(name=="A1") val = 2*cosS*TMath::Sqrt(1-cosS*cosS)*TMath::Cos(phiS);
	  else if(name=="A4") val = cosS;
	  return val;
	}, {"cos", "phi", "y"})
      );
    
    if(name=="UL") continue;
    
    for(int m = 0; m<=degs(get_pdf_type(name+"_x")); m++){
      std::string wname = name+"x_"+std::string(Form("%d",m));
      dlast = std::make_unique<RNode>(dlast->Define(wname.c_str(), [m,max_x,name,degs](double x){return cheb(x, 0.5*max_x, 1.0, degs(get_pdf_type(name+"_x")), m); } , {"x"} ));
    }
    for(int n = 0; n<=degs(get_pdf_type(name+"_y")); n++){
      std::string wname = name+"y_"+std::string(Form("%d",n));
      dlast = std::make_unique<RNode>(dlast->Define(wname.c_str(), [n,max_y,name,degs](double y){return cheb(TMath::Abs(y), 0.5*max_y, 1.0, degs(get_pdf_type(name+"_y")), n); }, {"y"} ));
    }
    
  }

  /*
  auto h1 = dlast->Histo1D({"cos","", 10, -1,1}, "cos" );
  TFile* fout = TFile::Open("out.root","RECREATE");
  fout->cd();
  h1->Write();
  fout->Close();
  if(snapshot){
    std::cout << "making a snapshot to disk" << std::endl;
    dlast->Snapshot(treeName, fileName);
  }
  */

  dlast = std::make_unique<RNode>(dlast->Define("pt",  [](RVecD p4lab){ return p4lab[0];}, {"p4lab"}));
  dlast = std::make_unique<RNode>(dlast->Define("eta", [](RVecD p4lab){ return TMath::Abs(p4lab[1]);}, {"p4lab"}));
  //dlast = std::make_unique<RNode>(dlast->Define("eta", [](RVecD p4lab){ return p4lab[1];}, {"p4lab"}));
  std::vector<ROOT::RDF::RResultPtr<TH1D> > histos1D;
  //histos1D.emplace_back(dlast->Histo1D({"pt", "", 100, 0, 100}, "pt"));
  
  std::vector<ROOT::RDF::RResultPtr<TH2D> > histos2D;
  auto prod2 = [](double a, double b){ return a*b;};
  auto prod6 = [](double a, double b, double c, double d, double e, double f){ return a*b*c*d*e*f;};
  auto prod8 = [](double a, double b, double c, double d, double e, double f, double g, double h){ return a*b*c*d*e*f*g*h;};

  if(run=="templates"){
    for(int i = 0; i<=degs(pdf_type::pdf_x); i++){
      for(int j = 0; j<=degs(pdf_type::pdf_y); j++){
	for(int k = 0; k<=degs(pdf_type::corr_x); k++){
	  for(int l = 0; l<=degs(pdf_type::corr_y); l++){
	    for(auto name : helicities){
	      if(name=="UL"){	    
		std::string whname(Form("w_%d_%d_%d_%d",i,j,k,l));
		whname += ("_"+name); 
		dlast = std::make_unique<RNode>(dlast->Define(whname, prod6, {"wM", Form("pdfx_%d",i), Form("pdfy_%d",j), Form("corrxy_%d",k), Form("corryx_%d",l), name } ));	    
		histos2D.emplace_back(dlast->Histo2D({whname.c_str(), "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", whname));
		continue;
	      }
	      for(int m = 0; m<=degs(get_pdf_type(name+"_x")); m++){
		std::string wx = name+"x_"+std::string(Form("%d",m)); 
		for(int n = 0; n<=degs(get_pdf_type(name+"_y")); n++){
		  std::string wy = name+"y_"+std::string(Form("%d",n)); 
		  std::string whname(Form("w_%d_%d_%d_%d_%d_%d",i,j,k,l,m,n));
		  whname += ("_"+name);
		  dlast = std::make_unique<RNode>(dlast->Define(whname, prod8, {"wM", Form("pdfx_%d",i), Form("pdfy_%d",j), Form("corrxy_%d",k), Form("corryx_%d",l), wx, wy, name } ));	    
		  histos2D.emplace_back(dlast->Histo2D({whname.c_str(), "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", whname));
		}
	      }
	      
	    }
	  }
	}
      }
    }
  }

  else if(run.find("closure")!=string::npos){

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
    double corr_xy[NMAX];
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

    dlast = std::make_unique<RNode>(dlast->Define("weights", [&](double x, double y, double cos, double phi, double mW )->RVecD {
	  RVecD out;
	  double   w{0.0};
	  double wMC{0.0};
	  for(auto name : helicities){	    

	    double whel  = 1.0;
	    double cosS = y>0 ? cos : -cos;
	    double phiS = y>0 ? phi : -phi;
	    if     (name=="UL") whel = (1+cosS*cosS);
	    else if(name=="A0") whel = (0.5*(1-3*cosS*cosS));
	    else if(name=="A1") whel = 2*cosS*TMath::Sqrt(1-cosS*cosS)*TMath::Cos(phiS);
	    else if(name=="A4") whel = cosS;
	    double whelMC  = whel;

	    double pdfx_ = 0.0;
	    for(int i = 0; i<=degs(pdf_type::pdf_x); i++){
	      pdfx_ += cheb(x, 0.5*max_x, 1.0, degs(pdf_type::pdf_x), i)*pdf_x[i];
	    }
	    //////
	    double pdfxMC_ = toy_x->Eval(x);
	    //std::cout << "pdfx " << pdfx_ << std::endl;

	    double pdfy_ = 0.0;
	    for(int j = 0; j<=degs(pdf_type::pdf_y); j++){
	      pdfy_ += cheb(TMath::Abs(y), 0.5*max_y, 1.0, degs(pdf_type::pdf_y), j)*pdf_y[j];
	    }
	    //////
	    double  pdfyMC_ = toy_y->Eval(y);
	    //std::cout << "pdfy " << pdfy_ << std::endl;
	    //std::cout << "pdfyMC " << pdfyMC_ << std::endl;

	    double corrxy_  = 0.0;
	    for(int k = 0; k<=degs(pdf_type::corr_x); k++){
	      for(int l = 0; l<=degs(pdf_type::corr_y); l++){       
		corrxy_ += 
		  cheb(x, 0.5*max_x, 1.0, degs(pdf_type::corr_x), k)*
		  cheb(TMath::Abs(y), 0.5*max_y, 1.0, degs(pdf_type::corr_y), l)*
		  corr_xy[ (degs(pdf_type::corr_y)+1)*k + l];
	      }
	    }
	    double corrxyMC_  = corrxy_;

	    whel   *= (pdfx_*pdfy_*corrxy_);
	    whelMC *= (pdfxMC_*pdfyMC_*corrxyMC_);
	    //std::cout << "corrxy " << corrxy_ << std::endl;

	    if(name=="A0"){
	      double pdfA0_ = 0.0;
	      for(int m = 0; m<=degs(pdf_type::A0_x); m++){
		for(int n = 0; n<=degs(pdf_type::A0_y); n++){
		 pdfA0_ += 
		   cheb(x, 0.5*max_x, 1.0, degs(pdf_type::A0_x), m)*
		   cheb(TMath::Abs(y), 0.5*max_y, 1.0, degs(pdf_type::A0_y), n)*
		   A0_xy[ (degs(pdf_type::A0_y)+1)*m + n];
		}
	      }
	      double pdfA0MC_ = pdfA0_;
	      whel   *= pdfA0_;
	      whelMC *= pdfA0MC_;
	      //std::cout << "pdfA0 " << pdfA0_ << std::endl;
	    }

	    // others Ai go here 
	    // ...

	    // sum it up
	    w   += whel;  
	    wMC += whelMC;  
	  }
	  out.emplace_back(w*mW);
	  out.emplace_back(wMC*mW);
	  return out;
	}, {"x","y","cos","phi","wM"}));

    dlast = std::make_unique<RNode>(dlast->Define("w",   [](RVecD weights){ return weights.at(0);}, {"weights"} ));
    dlast = std::make_unique<RNode>(dlast->Define("wMC", [](RVecD weights){ return weights.at(1);}, {"weights"} ));
    
    histos2D.emplace_back(dlast->Histo2D({"w",   "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", "w"));      
    histos2D.emplace_back(dlast->Histo2D({"wMC", "", nbinsX, xLow, xHigh, nbinsY, yLow, yHigh}, "eta", "pt", "wMC"));      
    fin->Close();
  }


  //histos2D.emplace_back(d1.Histo2D({"pteta", "", 50, -2.5, 2.5, 40, 20, 60}, "eta", "pt", "pdf_x"));
  auto colNames = dlast->GetColumnNames();
  //for (auto &&colName : colNames) std::cout << colName << std::endl;
  std::cout << colNames.size() << " columns created" << std::endl;

  fout->cd();
  std::cout << "Writing histos..." << std::endl;
  for(auto h : histos1D) h->Write();
  for(auto h : histos2D) h->Write();
  double total = *(dlast->Count());
  //for(auto sum : sums) std::cout << *sum/total << std::endl;
  fout->Close();
  return 1;
}
