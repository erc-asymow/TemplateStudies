#include<ROOT/RDataFrame.hxx>
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

std::vector<std::string> helicities = {"UL", "P0"};

constexpr double MW = 80.;
constexpr double GW = 2.0;

enum pdf_type { pdf_x=0, pdf_y, corr_x, corr_y, A0_x, A0_y};

auto degs = [](const pdf_type& pdf){
  switch(pdf){
  case pdf_type::pdf_x:  return  5; // f(x)
  case pdf_type::pdf_y:  return  5; // f(y|0)
  case pdf_type::corr_x: return  2; // P(x,.)
  case pdf_type::corr_y: return  2; // P(.,y)
  case pdf_type::A0_x:   return  3; // A0(x,.)
  case pdf_type::A0_y:   return  3; // A0(.,y)
  default: return 1;
  }
};

auto cheb = [](double x, double scale, double offset, int n, int m){
  double den = 0.;
  double num = 0.;
  for(unsigned int i = 0; i <= n ; i++){
    double xj = (TMath::Cos(i*TMath::Pi()/n) + offset)*scale;
    if(x==xj) return 1.0;                                                          
    double val = 0.;
    if(i==0) val = 0.5/(x-xj);
    else if(i==n) val = 0.5*(n%2==0?1:-1)/(x-xj);
    else val = (i%2==0?1:-1)/(x-xj);
    if(i==m) num=val;
    //cout << n << "," << m << "," << i << "," << x << "," << xj << "," << 
    //TMath::Cos(i*TMath::Pi()/n) << "," << offset << "," << scale << "==>" << val << endl;
    den += val;
  }
  //std::cout << x << "==>" <<  num << "," << den << std::endl;                                             
  return num/den;
};


/*
auto w_ijklmno = [](double x, double xscale, double xoffset, double y, double yscale, double yoffset, int i, int j, int k, int l, int m, int n, int o){
  double pdf0 = cheb(x, xscale, xoffset, n, degs(m));
  double pdf1 = cheb(y, yscale, yoffset, o, degs(m+1));
  return w_ijkl(x,xscale,xoffset,y,yscale,yoffset,i,j,k,l)*pdf0*pdf1;
};
*/

int main(int argc, char* argv[])
{

  ROOT::EnableImplicitMT();

  long nevents = 1000;
  std::string tag = "";
  if(argc>1) nevents = strtol(argv[1], NULL, 10);
  if(argc>2) tag = std::string(argv[2]);

  const double max_x = 0.15;
  const double max_y = 2.5;

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

  
  for(auto name : helicities){
    dlast = std::make_unique<RNode>(dlast->Define(name, [&](double cos, double phi, double y ){
	  double val  = 1.0;
	  double cosS = y>0 ? cos : -cos;
	  double phiS = y>0 ? phi : -phi;
	  if(name=="UL")      val = 1+cosS*cosS;
	  else if(name=="P0") val = 0.5*(1-3*cosS*cosS);
	  return val;
	}, {"cos", "phi", "y"})
      );
  }

  for(int i = 0; i<=degs(pdf_type::pdf_x); i++){
    for(int j = 0; j<=degs(pdf_type::pdf_y); j++){
      for(int k = 0; k<=degs(pdf_type::corr_x); k++){
	for(int l = 0; l<=degs(pdf_type::corr_y); l++){
	  auto w_ijkl = [i,j,k,l,max_x, max_y](double x, double y){	    
	    double pdf0 = cheb(x, 0.5*max_x, 1.0, degs(pdf_type::pdf_x), i);
	    double pdf1 = cheb(TMath::Abs(y), 0.5*max_y, 1.0, degs(pdf_type::pdf_y), j);
	    double pdf2 = cheb(x, 0.5*max_x, 1.0, degs(pdf_type::corr_x), k);
	    double pdf3 = cheb(TMath::Abs(y), 0.5*max_y, 1.0, degs(pdf_type::corr_y), l);
	    //cout << pdf0 << "," << pdf1 << "," << pdf2 << "," << pdf3 << endl;
	    return pdf0*pdf1*pdf2*pdf3;
	  };
	  dlast = std::make_unique<RNode>(dlast->Define(Form("w_%d%d%d%d", i,j,k,l), w_ijkl, {"x","y"} ));

	  for(int n = 0; n<=degs(pdf_type::A0_x); n++){
	    for(int o = 0; o<=degs(pdf_type::A0_y); o++){
	      auto w_ijkl_A0_no = [n,o,max_x, max_y, w_ijkl](double x, double y){	    
		double pdf0 = cheb(x, 0.5*max_x, 1.0, degs(pdf_type::A0_x), n);
		double pdf1 = cheb(TMath::Abs(y), 0.5*max_y, 1.0, degs(pdf_type::A0_y), o);
		return w_ijkl(x,y)*pdf0*pdf1;
	      };	      
	      dlast = std::make_unique<RNode>(dlast->Define(Form("w_%d%d%d%d_A0_%d%d", i,j,k,l,n,o), w_ijkl_A0_no, {"x","y"} ));
	    }
	  }	  

	}
      }
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
  std::vector<ROOT::RDF::RResultPtr<TH1D> > histos1D;
  histos1D.emplace_back(dlast->Histo1D({"pt", "", 100, 0, 100}, "pt"));
  
  std::vector<ROOT::RDF::RResultPtr<TH2D> > histos2D;
  for(int i = 0; i<=degs(pdf_type::pdf_x); i++){
    for(int j = 0; j<=degs(pdf_type::pdf_y); j++){
      for(int k = 0; k<=degs(pdf_type::corr_x); k++){
	for(int l = 0; l<=degs(pdf_type::corr_y); l++){
	  auto prod = [](double a, double b){ return a*b;};
	  std::string wname(Form("w_%d%d%d%d_UL",i,j,k,l));
	  dlast = std::make_unique<RNode>(dlast->Define(wname, prod, {Form("w_%d%d%d%d", i,j,k,l), "UL"} ));
	  histos2D.emplace_back(dlast->Histo2D({Form("pteta_%d%d%d%d_UL", i,j,k,l), "", 25, 0., 2.5, 40, 20, 60}, "eta", "pt", wname));
	}
      }
    }
  }
  //histos2D.emplace_back(d1.Histo2D({"pteta", "", 50, -2.5, 2.5, 40, 20, 60}, "eta", "pt", "pdf_x"));

  //auto colNames = dlast->GetColumnNames();
  //for (auto &&colName : colNames) std::cout << colName << std::endl;

  fout->cd();
  std::cout << "Writing histos..." << std::endl;
  for(auto h : histos1D) h->Write();
  for(auto h : histos2D) h->Write();
  fout->Close();
  return 1;
}
