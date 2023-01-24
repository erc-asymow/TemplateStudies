#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include "TRandom3.h"
#include "TVector.h"
#include "TVectorT.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include <TMatrixD.h>
#include <TStopwatch.h>
#include <ROOT/RVec.hxx>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;
using namespace ROOT;
typedef ROOT::VecOps::RVec<double> RVecD;
using ROOT::RDF::RNode; 

int main()
{  
  double x, y;
  double max_x = 0.8;
  double max_y = 3.5;
  int n_bins_x = 10;
  int n_bins_y = 6; 
  
  //TF2* wxy = new TF2("w_xy", "[0]*x/TMath::Power(x*x+[1], [2])*[3]/TMath::Sqrt(2*TMath::Pi()*[4])*TMath::Exp(-0.5*(y-[5])*(y-[5])/2/[4]/[4])*(0.1*(1.0 - 0.3*y*y)*TMath::Erf(5.0*x) + 1.0)", 0.0, max_x, 0.0, max_y);  
  TF2* wxy = new TF2("w_xy", "[0]*x/TMath::Power(x*x+[1], [2])*[3]/TMath::Sqrt(2*TMath::Pi()*[4])*TMath::Exp(-0.5*(y-[5])*(y-[5])/2/[4]/[4])", 0.0, max_x, -max_y, max_y);
  double sigma2_y = 4.0*4.0;
  wxy->SetParameter(0, 1.0);
  wxy->SetParameter(1, 0.00235);
  wxy->SetParameter(2, 1.0);
  wxy->SetParameter(3, 1.0);
  wxy->SetParameter(4, sigma2_y);
  wxy->SetParameter(5, 0.0);  
  
  TH2D* hist = new TH2D("hist","weights", n_bins_x, 0., max_x, n_bins_y, -max_y, max_y);
  for(int i=0; i<=n_bins_x; i++){ //<= or < 
    for(int k=0; k<=n_bins_y; k++){
      double x_bin_width = max_x/n_bins_x;
      double y_bin_width = 2*max_y/n_bins_y;
      double weight=wxy->Eval(i*x_bin_width,-max_y + k*y_bin_width); //at what value should the fct be eval?
      hist->Fill(i*x_bin_width, -max_y + k*y_bin_width, weight);
    }    
  }

  double integral = hist->Integral("width");
  std::cout<<"integral with width is: "<<integral;
  
  TCanvas *c1 = new TCanvas("c1", "weights", 600, 600);
  c1->cd();
  hist->Draw();
  c1->SaveAs("hist.pdf");
}


