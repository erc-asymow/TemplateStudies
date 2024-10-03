#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVector.h"
#include "TVectorT.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
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

using namespace boost::program_options;


int main(int argc, char* argv[])
{

  TStopwatch sw;
  sw.Start();

  variables_map vm;
  try
    {
      options_description desc{"Options"};
      desc.add_options()
	("help,h", "Help screen")
	("nevents",     value<long>()->default_value(2000000), "number of events")
	("tag",         value<std::string>()->default_value("closure"), "run type")
	("do_fit",      bool_switch()->default_value(false), "do_fit")
	("nbins",       value<int>()->default_value(200), "nbins")
	("frac",       value<float>()->default_value(0.5), "frac")
	("asym",       value<float>()->default_value(0.015), "asym")
	("lumiscale",       value<float>()->default_value(1.0), "lumiscale")
	("seed",        value<int>()->default_value(4357), "seed");

      store(parse_command_line(argc, argv, desc), vm);
      notify(vm);
      if (vm.count("help")){
	std::cout << desc << '\n';
	return 0;
      }
      if (vm.count("nevents"))    std::cout << "Number of events: " << vm["nevents"].as<long>() << '\n';
      if (vm.count("tag"))        std::cout << "Tag: " << vm["tag"].as<std::string>() << '\n';
    }
  catch (const error &ex)
    {
      std::cerr << ex.what() << '\n';
    }

  long nevents    = vm["nevents"].as<long>();
  std::string tag = vm["tag"].as<std::string>();
  int seed        = vm["seed"].as<int>();
  int nbins       = vm["nbins"].as<int>();
  float frac      = vm["frac"].as<float>();
  float asym      = vm["asym"].as<float>();
  float lumiscale    = vm["lumiscale"].as<float>();
  bool do_fit     = vm["do_fit"].as<bool>();

  std::vector<TRandom3*> rans = {};
  for(unsigned int i = 0; i < 10; i++){
    rans.emplace_back( new TRandom3(seed + i*10) );
  }

  TFile* fout = TFile::Open(("root/mcstat_"+tag+".root").c_str(), "RECREATE");
  TTree* tree = new TTree("tree", "");
  double mu0, mu1, err0, err1;
  tree->Branch("mu0", &mu0, "mu0/D");
  tree->Branch("mu1", &mu1, "mu1/D");
  tree->Branch("err0", &err0, "err0/D");
  tree->Branch("err1", &err1, "err1/D");
  
  TH1D* h_true_0 = new TH1D("h_true_0", "", nbins, 0.0, 1.0);
  TH1D* h_true_1 = new TH1D("h_true_1", "", nbins, 0.0, 1.0);
  for(int ib = 1; ib<=nbins; ib++){
    double ival_0 = nevents/nbins;
    ival_0  *= (1-frac)*(1-asym);
    double ival_1 = ib<=nbins/2 ? nevents/nbins*(1+asym) : nevents/nbins*(1-asym);
    ival_1  *= (frac);
    h_true_0->SetBinContent(ib, ival_0);
    h_true_1->SetBinContent(ib, ival_1);
  }

  double rescale = nevents/(h_true_0->Integral() + h_true_1->Integral());
  h_true_0->Scale(rescale);
  h_true_1->Scale(rescale);
  
  TH1D* h_true_tot = (TH1D*)h_true_0->Clone("h_true_tot");
  h_true_tot->Add(h_true_1, 1.0);
  
  cout << "h_true_0: integral = " << h_true_0->Integral() << endl;
  cout << "h_true_1: integral = " << h_true_1->Integral() << endl;
  cout << "h_true_tot: integral = " << h_true_tot->Integral() << endl;

  MatrixXd A_true = MatrixXd::Zero(nbins, 2);
  MatrixXd invV_true = MatrixXd::Zero(nbins, nbins);

  for(unsigned int ir=0; ir<nbins; ir++){
    A_true(ir, 0) = h_true_0->GetBinContent(ir+1);
    A_true(ir, 1) = h_true_1->GetBinContent(ir+1);
  }
  for(unsigned int ir=0; ir<nbins; ir++){
    invV_true(ir,ir) = 1./h_true_tot->GetBinContent(ir+1);
  }

  MatrixXd C_true = ( A_true.transpose()*invV_true*A_true ).inverse();
  cout << "True errors on mu_[0,1] = [" << TMath::Sqrt(C_true(0,0)) << "," << TMath::Sqrt(C_true(1,1)) << "]" << endl;

  /*
  MatrixXd C_inv = MatrixXd::Zero(2, 2);
  for(unsigned int ib=0; ib<nbins; ib++){
    C_inv += (A_true.row(ib).transpose()*A_true.row(ib))/h_true_tot->GetBinContent(ib+1);
  }
  MatrixXd C_alt = C_inv.inverse();
  cout << "Errors on mu_[0,1] = [" << TMath::Sqrt(C_alt(0,0)) << "," << TMath::Sqrt(C_alt(1,1)) << "]" << endl;
  */

  MatrixXd A_itoy = MatrixXd::Zero(nbins, 2);
  MatrixXd invV_itoy = MatrixXd::Zero(nbins, nbins);
  VectorXd y_itoy(nbins);
  VectorXd y0_itoy(nbins);
  
  unsigned int ntoys = 3000;
  for(unsigned int itoy=0; itoy<ntoys; itoy++){

    MatrixXd invC_itoy = MatrixXd::Zero(2, 2);
    for(unsigned int ir=0; ir<nbins; ir++){      
      A_itoy(ir,0) = rans[0]->Poisson( A_true(ir, 0)*lumiscale ) / lumiscale;
      A_itoy(ir,1) = rans[0]->Poisson( A_true(ir, 1)*lumiscale ) / lumiscale;
      double norm_rowi = A_itoy.row(ir).sum();
      y_itoy(ir) = rans[0]->Poisson( A_true.row(ir).sum() );
      y0_itoy(ir) = norm_rowi;
      invV_itoy(ir,ir) = 1./norm_rowi;
      invC_itoy += (A_itoy.row(ir).transpose()*A_itoy.row(ir))/norm_rowi;
    }

    MatrixXd C_itoy = invC_itoy.inverse();
    VectorXd x_itoy = C_itoy*A_itoy.transpose()*invV_itoy*(y_itoy - y0_itoy);
    mu0 = x_itoy(0);
    mu1 = x_itoy(1);
    err0 = TMath::Sqrt(C_itoy(0,0));
    err1 = TMath::Sqrt(C_itoy(1,1));
    tree->Fill();
  }

    
  fout->cd();
  h_true_0->Write();
  h_true_1->Write();
  tree->Write();
  
  sw.Stop();

  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;

  
  fout->Close(); 
  for(auto r : rans) delete r;

  return 1;

}
