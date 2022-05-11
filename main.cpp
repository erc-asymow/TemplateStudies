#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include "TRandom3.h"
#include "TVector.h"
#include "TVectorT.h"
#include "TMath.h"
#include <TMatrixD.h>
#include <ROOT/RVec.hxx>
#include <stdlib.h>

using namespace std;
using namespace ROOT;
typedef ROOT::VecOps::RVec<double> RVecD;
using ROOT::RDF::RNode; 

//std::vector<std::string> helicities = {"UL", "A0", "A1", "A4"};
std::vector<std::string> helicities = {"A4"};

constexpr double MW = 80.;
constexpr double GW = 2.0;

enum pdf_type { pdf_x=0, pdf_y, corr_x, corr_y, A0_x, A0_y, A1_x, A1_y, A4_x, A4_y, unknown};

auto degs = [](const pdf_type& pdf){
  switch(pdf){
  case pdf_type::pdf_x:  return  2; // f(x)
  case pdf_type::pdf_y:  return  2; // f(y|0)
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
  
  long nevents = 1000;
  std::string tag = "";
  if(argc>1) nevents = strtol(argv[1], NULL, 10);
  if(argc>2) tag = std::string(argv[2]);

  const double max_x = 0.15;
  const double max_y = 2.5;
  const int nbinsX   = 12; 
  const double xLow  = 0.0;
  const double xHigh = +2.5;
  const int nbinsY   = 15; 
  const double yLow  = 25.;
  const double yHigh = 55.;

  // We prepare an input tree to run on
  TFile* fout = TFile::Open(("histos_"+tag+".root").c_str(), "RECREATE");

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
						  double qT2 = x*Q*Q;
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
    dlast = std::make_unique<RNode>(dlast->Define(Form("pdfx_%d", i), [i,max_x](double x){return cheb(x, 0.5*max_x, 1.0, degs(pdf_type::pdf_x), i);}, {"x"} ));
    sums.emplace_back( dlast->Sum<double>( Form("pdfx_%d", i)) );
  }

  for(int j = 0; j<=degs(pdf_type::pdf_y); j++){
    dlast = std::make_unique<RNode>(dlast->Define(Form("pdfy_%d", j), [j,max_y](double y){return cheb(TMath::Abs(y), 0.5*max_y, 1.0, degs(pdf_type::pdf_y), j);} , {"y"} ));
    //sums.emplace_back( dlast->Sum( Form("pdfy_%d", j)) );
  }
 
  for(int k = 0; k<=degs(pdf_type::corr_x); k++){
    dlast = std::make_unique<RNode>(dlast->Define(Form("corrxy_%d", k), [k,max_x](double x){return cheb(x, 0.5*max_x, 1.0, degs(pdf_type::corr_x), k); } , {"x"} ));
  }
  for(int l = 0; l<=degs(pdf_type::corr_y); l++){
    dlast = std::make_unique<RNode>(dlast->Define(Form("corryx_%d", l), [l,max_y](double y){return cheb(TMath::Abs(y), 0.5*max_y, 1.0, degs(pdf_type::corr_y), l); }, {"y"} ));
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
      dlast = std::make_unique<RNode>(dlast->Define(wname.c_str(), [m,max_x,name](double x){return cheb(x, 0.5*max_x, 1.0, degs(get_pdf_type(name+"_x")), m); } , {"x"} ));
    }
    for(int n = 0; n<=degs(get_pdf_type(name+"_y")); n++){
      std::string wname = name+"y_"+std::string(Form("%d",n));
      dlast = std::make_unique<RNode>(dlast->Define(wname.c_str(), [n,max_y,name](double y){return cheb(TMath::Abs(y), 0.5*max_y, 1.0, degs(get_pdf_type(name+"_y")), n); }, {"y"} ));
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

  //histos2D.emplace_back(d1.Histo2D({"pteta", "", 50, -2.5, 2.5, 40, 20, 60}, "eta", "pt", "pdf_x"));
  auto colNames = dlast->GetColumnNames();
  //for (auto &&colName : colNames) std::cout << colName << std::endl;
  std::cout << colNames.size() << " columns created" << std::endl;

  fout->cd();
  std::cout << "Writing histos..." << std::endl;
  for(auto h : histos1D) h->Write();
  for(auto h : histos2D) h->Write();
  double total = *(dlast->Count());
  for(auto sum : sums) std::cout << *sum/total << std::endl;
  fout->Close();
  return 1;
}
