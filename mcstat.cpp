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
#include "TString.h"
#include "TGraphErrors.h"
#include <TMatrixD.h>
#include <TStopwatch.h>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <boost/program_options.hpp>
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/FCNGradientBase.h"

#include <Eigen/Core>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;
using namespace ROOT;
using namespace ROOT::Minuit2;

using namespace boost::program_options;


class Likelihood : public FCNBase {

public:
  
  Likelihood(const int& debug, const int& nbins, const bool& BBlite, const bool& doPoisson, const bool& profileMCNP, const bool& decorrelate )
    : errorDef_(1.0), debug_(debug), nbins_(nbins),  bblite_(BBlite), doPoisson_(doPoisson), profileMCNP_(profileMCNP), decorrelate_(decorrelate) {
    //npars_ = BBlite ? (2 + nbins) : (2 + 2*nbins);
    removeBB_ = false;
    npars_ = 2;
    ndof_ = nbins - 2;
    data_.reserve( nbins );
    dataErr2_.reserve( nbins );
    mc0_.reserve( nbins );
    mc1_.reserve( nbins );
    mc0err_.reserve( nbins );
    mc1err_.reserve( nbins );
    for(unsigned int ibin = 0; ibin<nbins; ibin++) data_.push_back(0.);
    for(unsigned int ibin = 0; ibin<nbins; ibin++) dataErr2_.push_back(0.);
    for(unsigned int ibin = 0; ibin<nbins; ibin++) mc0_.push_back(0.);
    for(unsigned int ibin = 0; ibin<nbins; ibin++) mc1_.push_back(0.);
    for(unsigned int ibin = 0; ibin<nbins; ibin++) mc0err_.push_back(0.);
    for(unsigned int ibin = 0; ibin<nbins; ibin++) mc1err_.push_back(0.);
  }

  unsigned int get_n_dof(){ return ndof_;}

  void set_removeBB(const bool& set){ removeBB_ = set; }

  void restore_dataErr2(){
    for(unsigned int ibin = 0; ibin<nbins_; ibin++){
      dataErr2_[ibin] = data_[ibin];
    }
  }
  
  void set_data(const unsigned int& ibin, const double& val){
    data_[ibin] = val;
    dataErr2_[ibin] = val;
  }
  void set_dataErr2(const unsigned int& ibin, const double& val){
    dataErr2_[ibin] = val; 
  }
  void set_mc0(const unsigned int& ibin, const double& val){
    mc0_[ibin] = val; 
  }
  void set_mc1(const unsigned int& ibin, const double& val){
    mc1_[ibin] = val; 
  }
  void set_mc0err(const unsigned int& ibin, const double& val){
    mc0err_[ibin] = val; 
  }
  void set_mc1err(const unsigned int& ibin, const double& val){
    mc1err_[ibin] = val; 
  }
  
  virtual double Up() const {return errorDef_;}
  virtual void SetErrorDef(double def) {errorDef_ = def;}

  virtual double operator()(const vector<double>&) const;
  //virtual vector<double> Gradient(const vector<double>& ) const;
  //virtual bool CheckGradient() const {return true;} 

private:

  vector<double> data_;
  vector<double> dataErr2_;
  vector<double> mc0_;
  vector<double> mc1_;
  vector<double> mc0err_;
  vector<double> mc1err_;
  unsigned int nbins_;
  unsigned int npars_;
  unsigned int ndof_;
  bool bblite_;
  int debug_;
  bool profileMCNP_;
  bool doPoisson_;
  bool decorrelate_;
  double errorDef_;
  bool removeBB_;
};

// par[0] = mu0, par[1] = mu1, par[2,..] = BB 
double Likelihood::operator()(const vector<double>& par) const {  
  double val  = 0.0;
  double par0 = par[0];
  double par1 = par[1];
  if(decorrelate_){
    par0 =  par[0] + par[1];
    par1 = -par[0] + par[1];
  }
  for(unsigned int ibin=0; ibin<nbins_; ibin++){
    double data_i = data_[ibin];
    double dataErr2_i = dataErr2_[ibin];
    double mc_tot = (mc0_[ibin]+mc1_[ibin]);
    double mc_tot_err = TMath::Sqrt( mc0err_[ibin]*mc0err_[ibin] + mc1err_[ibin]*mc1err_[ibin] );
    double rel_err_mc = mc_tot_err/mc_tot;
    double rel_err_mc0 = mc0err_[ibin]/mc0_[ibin];
    double rel_err_mc1 = mc1err_[ibin]/mc1_[ibin];

    if(removeBB_){
      double exp_data_i = (1+par0)*mc0_[ibin] + (1+par1)*mc1_[ibin];
      double res  = data_i - exp_data_i;
      double res2 = res*res;
      double chi2 = 0.5*res2/dataErr2_i;
      if(doPoisson_){
	chi2 *= (1.0 - 2./3.*(res/dataErr2_i) + 0.5*res2/dataErr2_i/dataErr2_i /* + ... */ );
      }
      val += chi2;
      continue;
    }

    if(profileMCNP_){
      double exp_data_i = bblite_ ?
	((1+par0)*mc0_[ibin] + (1+par1)*mc1_[ibin])*(1+par[2+ibin]) :
	(1.0 + par0)*(1.0 + par[2+ibin])*mc0_[ibin] + (1.0 + par1)*(1.0 + par[2+nbins_+ibin])*mc1_[ibin];
      double res = data_i - exp_data_i;
      double res2 = res*res;
      double prior  = bblite_ ?  0.5*par[2+ibin]*par[2+ibin]/(rel_err_mc*rel_err_mc) : 0.0;
      if(!bblite_){
	prior += 0.5*par[2+ibin]*par[2+ibin]/(rel_err_mc0*rel_err_mc0);
	prior += 0.5*par[2+nbins_+ibin]*par[2+nbins_+ibin]/(rel_err_mc1*rel_err_mc1);
      }
      double chi2 = 0.5*res2/dataErr2_i;
      if(doPoisson_){
	chi2 *= (1.0 - 2./3.*(res/dataErr2_i) + 0.5*res2/dataErr2_i/dataErr2_i /* + ... */ );
      }
      val  += chi2;
      val  += prior;
    }
    else{
      double exp_data_i = (1+par0)*mc0_[ibin] + (1+par1)*mc1_[ibin] ;
      double res  = data_i - exp_data_i;
      double res2 = res*res;
      double chi2   = bblite_ ?
	0.5*res2/( dataErr2_i + TMath::Power(rel_err_mc*exp_data_i, 2.0 ) ) :
	0.5*res2/( dataErr2_i + TMath::Power(mc0err_[ibin]*(1.0 + par0), 2.0) + TMath::Power(mc1err_[ibin]*(1.0 + par1), 2.0) ) ;
      val  += chi2;
    }    
  }    
  return val;
}

double get_quantile(TH1D* h){
  double q[1];
  vector<double> probsum = {0.5};
  if(h->Integral()>0)
    h->GetQuantiles(1, q, probsum.data());
  return q[0];
}
  
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
	("ntoys",       value<long>()->default_value(3000), "number of toys")
	("ntoysFC",       value<long>()->default_value(3000), "number of toys")
	("tag",         value<std::string>()->default_value("closure"), "run type")
	("doFC",      bool_switch()->default_value(false), "doFC")
	("doFCcheat",      bool_switch()->default_value(false), "doFCcheat")
	("doBarlett",      bool_switch()->default_value(false), "doBarlett")
	("saveHistos",       bool_switch()->default_value(false), "saveHistos")
	("FCfixToTrue",      bool_switch()->default_value(false), "FCfixToTrue")
	("verbose",      bool_switch()->default_value(false), "verbose")
	("profileMCNP",      bool_switch()->default_value(false), "profileMCNP")
	("doPoisson",      bool_switch()->default_value(false), "doPoisson")
	("decorrelate",    bool_switch()->default_value(false), "decorrelate")
	("bExact",    bool_switch()->default_value(false), "bExact")
	("jExact",    bool_switch()->default_value(false), "jExact")
	("computeJtildeError", bool_switch()->default_value(false), "computeJtildeError")
	("computePropError", bool_switch()->default_value(false), "computePropError")
	("computePropErrorTrue", bool_switch()->default_value(false), "computePropErrorTrue")
	("computePropErrorBB", bool_switch()->default_value(false), "computePropErrorBB")
	("nbins",       value<int>()->default_value(200), "nbins")
	("nsigmas",     value<int>()->default_value(5), "nsigmas")
	("frac",       value<float>()->default_value(0.5), "frac")
	("asym",       value<float>()->default_value(0.015), "asym")
	("alpha",       value<float>()->default_value(45.), "alpha")
	("lumiscale",   value<float>()->default_value(1.0), "lumiscale")
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
  long ntoys    = vm["ntoys"].as<long>();
  long ntoysFC    = vm["ntoysFC"].as<long>();
  std::string tag = vm["tag"].as<std::string>();
  int seed        = vm["seed"].as<int>();
  int nbins       = vm["nbins"].as<int>();
  int nsigmas     = vm["nsigmas"].as<int>();
  float frac      = vm["frac"].as<float>();
  float asym      = vm["asym"].as<float>();
  float alpha      = vm["alpha"].as<float>();
  float lumiscale    = vm["lumiscale"].as<float>();
  bool doFC     = vm["doFC"].as<bool>();
  bool doFCcheat = vm["doFCcheat"].as<bool>();
  bool doBarlett = vm["doBarlett"].as<bool>();
  bool saveHistos  = vm["saveHistos"].as<bool>();
  bool FCfixToTrue = vm["FCfixToTrue"].as<bool>();
  bool verbose     = vm["verbose"].as<bool>();
  bool profileMCNP = vm["profileMCNP"].as<bool>();
  bool decorrelate = vm["decorrelate"].as<bool>();
  bool bExact = vm["bExact"].as<bool>();
  bool jExact = vm["jExact"].as<bool>();
  bool doPoisson = vm["doPoisson"].as<bool>();
  bool computeJtildeError = vm["computeJtildeError"].as<bool>();
  bool computePropError = vm["computePropError"].as<bool>();
  bool computePropErrorTrue = vm["computePropErrorTrue"].as<bool>();
  bool computePropErrorBB = vm["computePropErrorBB"].as<bool>();

  alpha *= (TMath::Pi()/180.);
  
  std::vector<TRandom3*> rans = {};
  for(unsigned int i = 0; i < 10; i++){
    rans.emplace_back( new TRandom3(seed + i*10) );
  }

  TFile* fout = TFile::Open(("root/mcstat_"+tag+".root").c_str(), "RECREATE");
  TTree* tree = new TTree("tree", "");
  TTree* treeFC = new TTree("treeFC", "");
  double mu_data, err_data;
  double err_adhocdata, err_adhocdata5s;
  double err_propdata, err_propdata5s;
  double mu_poisdata, err_poisdata;
  double mu_poisdata_BB, err_poisdata_BB;
  double mu_poisdata_BBfull, err_poisdata_BBfull;
  double errPLRLow_poisdata_BBfull,errPLRUp_poisdata_BBfull;
  double mu_poisdata5s, err_poisdata5s;
  double mu_poisdata5s_BB, err_poisdata5s_BB;
  double mu_poisdata5s_BBfull, err_poisdata5s_BBfull;
  double errPLRLow_poisdata5s_BBfull, errPLRUp_poisdata5s_BBfull;
  double errFCLow_poisdata_BBfull, errFCUp_poisdata_BBfull;
  double errFCLow_poisdata5s_BBfull, errFCUp_poisdata5s_BBfull;
  double mu_data5s, err_data5s;
  double mu_mc, err_mc;
  double mu_data_BB, err_data_BB;
  double mu_data5s_BB, err_data5s_BB;
  double mu_true, err_true;
  double mu_true5s, err_true5s;
  double tstat, tstat5s;
  tree->Branch("mu_data",  &mu_data, "mu_data/D");
  tree->Branch("err_data", &err_data, "err_data/D");
  tree->Branch("err_adhocdata", &err_adhocdata, "err_adhocdata/D");
  tree->Branch("err_adhocdata5s", &err_adhocdata5s, "err_adhocdata5s/D");
  tree->Branch("err_propdata", &err_propdata, "err_propdata/D");
  tree->Branch("err_propdata5s", &err_propdata5s, "err_propdata5s/D");
  tree->Branch("mu_poisdata",  &mu_poisdata, "mu_poisdata/D");
  tree->Branch("err_poisdata", &err_poisdata, "err_poisdata/D");
  tree->Branch("mu_poisdata_BB",  &mu_poisdata_BB, "mu_poisdata_BB/D");
  tree->Branch("err_poisdata_BB", &err_poisdata_BB, "err_poisdata_BB/D");
  tree->Branch("mu_poisdata_BBfull",  &mu_poisdata_BBfull, "mu_poisdata_BBfull/D");
  tree->Branch("err_poisdata_BBfull", &err_poisdata_BBfull, "err_poisdata_BBfull/D");
  tree->Branch("errPLRLow_poisdata_BBfull", &errPLRLow_poisdata_BBfull, "errPLRLow_poisdata_BBfull/D");
  tree->Branch("errPLRUp_poisdata_BBfull", &errPLRUp_poisdata_BBfull, "errPLRUp_poisdata_BBfull/D");
  treeFC->Branch("errFCLow_poisdata_BBfull", &errFCLow_poisdata_BBfull, "errFCLow_poisdata_BBfull/D");
  treeFC->Branch("errFCUp_poisdata_BBfull", &errFCUp_poisdata_BBfull, "errFCUp_poisdata_BBfull/D");
  tree->Branch("mu_poisdata5s",  &mu_poisdata5s, "mu_poisdata5s/D");
  tree->Branch("err_poisdata5s", &err_poisdata5s, "err_poisdata5s/D");
  tree->Branch("mu_poisdata5s_BB",  &mu_poisdata5s_BB, "mu_poisdata5s_BB/D");
  tree->Branch("err_poisdata5s_BB", &err_poisdata5s_BB, "err_poisdata5s_BB/D");
  tree->Branch("mu_poisdata5s_BBfull",  &mu_poisdata5s_BBfull, "mu_poisdata5s_BBfull/D");
  tree->Branch("err_poisdata5s_BBfull", &err_poisdata5s_BBfull, "err_poisdata5s_BBfull/D");
  tree->Branch("errPLRLow_poisdata5s_BBfull", &errPLRLow_poisdata5s_BBfull, "errPLRLow_poisdata5s_BBfull/D");
  tree->Branch("errPLRUp_poisdata5s_BBfull", &errPLRUp_poisdata5s_BBfull, "errPLRUp_poisdata5s_BBfull/D");
  treeFC->Branch("errFCLow_poisdata5s_BBfull", &errFCLow_poisdata5s_BBfull, "errFCLow_poisdata5s_BBfull/D");
  treeFC->Branch("errFCUp_poisdata5s_BBfull", &errFCUp_poisdata5s_BBfull, "errFCUp_poisdata5s_BBfull/D");
  tree->Branch("mu_data5s",  &mu_data5s, "mu_data5s/D");
  tree->Branch("err_data5s", &err_data5s, "err_data5s/D");
  tree->Branch("mu_data_BB",  &mu_data_BB, "mu_data_BB/D");
  tree->Branch("err_data_BB", &err_data_BB, "err_data_BB/D");
  tree->Branch("mu_data5s_BB",  &mu_data5s_BB, "mu_data5s_BB/D");
  tree->Branch("err_data5s_BB", &err_data5s_BB, "err_data5s_BB/D");
  tree->Branch("mu_mc",  &mu_mc, "mu_mc/D");
  tree->Branch("err_mc", &err_mc, "err_mc/D");
  tree->Branch("mu_true",  &mu_true, "mu_true/D");
  tree->Branch("err_true", &err_true, "err_true/D");
  tree->Branch("mu_true5s",  &mu_true5s, "mu_true5s/D");
  tree->Branch("err_true5s", &err_true5s, "err_true5s/D");
  treeFC->Branch("tstat",  &tstat, "tstat/D");
  treeFC->Branch("tstat5s",  &tstat5s, "tstat5s/D");
   
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
  MatrixXd J_true = MatrixXd::Zero(nbins, 2);

  for(unsigned int ir=0; ir<nbins; ir++){
    A_true(ir, 0) = h_true_0->GetBinContent(ir+1);
    A_true(ir, 1) = h_true_1->GetBinContent(ir+1);
    //cout << A_true(ir, 0) << " : " << A_true(ir, 1) << endl;
    J_true(ir, 0) = A_true(ir, 0);
    J_true(ir, 1) = A_true(ir, 1);
    if(decorrelate){
      J_true(ir, 0) = TMath::Sqrt(2.)*(TMath::Cos(alpha)*A_true(ir, 0) - TMath::Sin(alpha)*A_true(ir, 1));
      J_true(ir, 1) = TMath::Sqrt(2.)*(TMath::Sin(alpha)*A_true(ir, 0) + TMath::Cos(alpha)*A_true(ir, 1));
    }
  }
  for(unsigned int ir=0; ir<nbins; ir++){
    invV_true(ir,ir) = 1./( A_true(ir, 0) + A_true(ir, 1));
  }

  // A_true --> J_true
  MatrixXd C_true = ( J_true.transpose()*invV_true*J_true ).inverse();
  //cout << "True errors on mu_[0,1] = [" << TMath::Sqrt(C_true(0,0)) << "," << TMath::Sqrt(C_true(1,1)) << "]" << endl;
  double rho_true = C_true(0,1)/TMath::Sqrt(C_true(0,0)*C_true(1,1));
  double condition_true = 0.;
  cout << "Correlation: " << rho_true << ", err: " << TMath::Sqrt(C_true(0,0)) << endl;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(C_true);
  if (eigensolver.info() != Eigen::Success){
    cout << "Could not eigendecompose U" << endl;
  }
  else{    
    auto eigenvals = eigensolver.eigenvalues();
    auto eigenvecs = eigensolver.eigenvectors();
    condition_true = eigenvals(1)/eigenvals(0);
    cout << "Condition number of V: " << eigenvals(1) << "/" << eigenvals(0) << " = " << condition_true << endl;
  }

  
  mu_true = 0.0;
  err_true = TMath::Sqrt(C_true(0,0));

  mu_true5s = mu_true + err_true*nsigmas;
  
  MatrixXd A_true5s = MatrixXd::Zero(nbins, 2);
  MatrixXd invV_true5s = MatrixXd::Zero(nbins, nbins);
  MatrixXd J_true5s = MatrixXd::Zero(nbins, 2);
  for(unsigned int ir=0; ir<nbins; ir++){
    A_true5s(ir, 0) = (1.0 + mu_true5s)*h_true_0->GetBinContent(ir+1);
    A_true5s(ir, 1) = decorrelate ? (1.0 - mu_true5s)*h_true_1->GetBinContent(ir+1) : h_true_1->GetBinContent(ir+1);
    //J_true5s(ir, 0) = A_true5s(ir, 0);
    //J_true5s(ir, 1) = A_true5s(ir, 1);
    J_true5s(ir, 0) = h_true_0->GetBinContent(ir+1);
    J_true5s(ir, 1) = h_true_1->GetBinContent(ir+1);
    if(decorrelate){
      //J_true5s(ir, 0) = A_true5s(ir, 0) - A_true5s(ir, 1);
      //J_true5s(ir, 1) = A_true5s(ir, 0) + A_true5s(ir, 1);
      J_true5s(ir, 0) = h_true_0->GetBinContent(ir+1) - h_true_1->GetBinContent(ir+1);
      J_true5s(ir, 1) = h_true_0->GetBinContent(ir+1) + h_true_1->GetBinContent(ir+1);
    }
  }
  for(unsigned int ir=0; ir<nbins; ir++){
    invV_true5s(ir,ir) = 1./( A_true5s(ir, 0) + A_true5s(ir, 1));
  }
  MatrixXd C_true5s = ( J_true5s.transpose()*invV_true5s*J_true5s ).inverse();
  err_true5s = TMath::Sqrt(C_true5s(0,0));

  MatrixXd A_itoy = MatrixXd::Zero(nbins, 2);
  MatrixXd J_itoy = MatrixXd::Zero(nbins, 2);
  MatrixXd invV_itoy      = MatrixXd::Zero(nbins, nbins);
  MatrixXd invsqrtV_itoy  = MatrixXd::Zero(nbins, nbins);
  MatrixXd invsqrtV5s_itoy = MatrixXd::Zero(nbins, nbins);
  MatrixXd invV5s_itoy    = MatrixXd::Zero(nbins, nbins);
  MatrixXd invV_BB_itoy   = MatrixXd::Zero(nbins, nbins);
  MatrixXd invV5s_BB_itoy = MatrixXd::Zero(nbins, nbins);
  MatrixXd invsqrtV_BB_itoy  = MatrixXd::Zero(nbins, nbins);
  MatrixXd invsqrtV5s_BB_itoy = MatrixXd::Zero(nbins, nbins);
  
  MatrixXd invVMC_itoy    = MatrixXd::Zero(nbins, nbins);
  VectorXd y_itoy(nbins);
  VectorXd yFC_itoy(nbins);
  VectorXd y5s_itoy(nbins);
  VectorXd yMC_itoy(nbins);
  VectorXd ynom_itoy(nbins);
  
  double prob_data = 0.;
  double prob_data_BB = 0.;
  double prob_poisdata = 0.;
  double prob_poisdata_BB = 0.;
  double prob_poisdata_BBfull = 0.;
  double prob_poisdata_BBfullPLR = 0.;
  double prob_poisdata_BBfullFC = 0.;

  double prob_data5s = 0.;  
  double prob_data5s_BB = 0.;
  double prob_poisdata5s = 0.;
  double prob_poisdata5s_BB = 0.;
  double prob_poisdata5s_BBfull = 0.;
  double prob_poisdata5s_BBfullPLR = 0.;
  double prob_poisdata5s_BBfullFC = 0.;
  double prob_mc = 0.;

  double prob_adhocdata = 0.;
  double prob_adhocdata5s = 0.;
  double prob_propdata = 0.;
  double prob_propdata5s = 0.;
  
  unsigned int maxfcn(numeric_limits<unsigned int>::max());
  double tolerance(0.001);
  int verbosity = int(nevents<2);
  int npars_lite = 2 ;
  if(profileMCNP) npars_lite += nbins;
  int npars_full = 2 ;
  if(profileMCNP) npars_full +=	2*nbins;

  ROOT::Minuit2::MnPrint::SetGlobalLevel(verbosity);

  MatrixXd Vmin     = MatrixXd::Zero(npars_lite,npars_lite);    
  VectorXd xmin     = VectorXd::Zero(npars_lite);
  VectorXd xmin_Err = VectorXd::Zero(npars_lite);

  MatrixXd Vmin_lite    = MatrixXd::Zero(npars_lite,npars_lite);    
  VectorXd xmin_lite    = VectorXd::Zero(npars_lite);
  VectorXd xmin_liteErr = VectorXd::Zero(npars_lite);

  MatrixXd Vmin_full    = MatrixXd::Zero(npars_full,npars_full);    
  VectorXd xmin_full    = VectorXd::Zero(npars_full);
  VectorXd xmin_fullErr = VectorXd::Zero(npars_full);

  MatrixXd Vmin_5s        = MatrixXd::Zero(npars_lite,npars_lite);    
  VectorXd xmin_5s        = VectorXd::Zero(npars_lite);
  VectorXd xmin_5sErr     = VectorXd::Zero(npars_lite);

  MatrixXd Vmin_lite5s    = MatrixXd::Zero(npars_lite,npars_lite);    
  VectorXd xmin_lite5s    = VectorXd::Zero(npars_lite);
  VectorXd xmin_lite5sErr = VectorXd::Zero(npars_lite);

  MatrixXd Vmin_full5s    = MatrixXd::Zero(npars_full,npars_full);    
  VectorXd xmin_full5s    = VectorXd::Zero(npars_full);
  VectorXd xmin_full5sErr = VectorXd::Zero(npars_full);

  //int ntoysFC = ntoys;    
  if(doFC){

    for(unsigned int itoy=0; itoy<ntoysFC; itoy++){

      if(itoy%100==0) cout << "Doing FC toy " << itoy << " / " << ntoysFC << endl;

      // idata=0: data nominal, idata=1: data 5s
      for(int idata=0; idata<2; idata++){       

 	Likelihood* likelihoodFC_full = new Likelihood(0, nbins, false, doPoisson, false, decorrelate);
	likelihoodFC_full->SetErrorDef(0.5);

	for(unsigned int ir=0; ir<nbins; ir++){      
	  A_itoy(ir,0) = rans[1+idata]->Poisson( A_true(ir, 0)*lumiscale ) / lumiscale;
	  A_itoy(ir,1) = rans[1+idata]->Poisson( A_true(ir, 1)*lumiscale ) / lumiscale;
	  //ynom_itoy(ir) = A_itoy.row(ir).sum();
	  if(idata==0)
	    y_itoy(ir)   = rans[1+idata]->Poisson( A_true(ir,0) + A_true(ir,1) );
	  else{
	    y_itoy(ir)   = rans[1+idata]->Poisson( A_true(ir,0)*(1.0 + mu_true5s) + A_true(ir,1) );
	    if(decorrelate)
	      y_itoy(ir)   = rans[1+idata]->Poisson( A_true(ir,0)*(1.0 + mu_true5s) + A_true(ir,1)*(1.0 - mu_true5s) );
	  }
	  likelihoodFC_full->set_mc0(ir, A_itoy(ir,0));
	  likelihoodFC_full->set_mc1(ir, A_itoy(ir,1));
	  likelihoodFC_full->set_mc0err(ir, TMath::Sqrt(A_itoy(ir,0)/lumiscale));
	  likelihoodFC_full->set_mc1err(ir, TMath::Sqrt(A_itoy(ir,1)/lumiscale));
	  likelihoodFC_full->set_data(ir, y_itoy(ir) );
	}
	
	MnUserParameters uparFC_full;
	uparFC_full.Add("mu0", 0.0, 0.01);
	uparFC_full.Add("mu1", 0.0, 0.01);      
	
	MnMigrad migradFC_full(*likelihoodFC_full, uparFC_full, 1);
	FunctionMinimum minFC_full = migradFC_full(maxfcn, tolerance);
	double mu0 = minFC_full.UserState().Value(0);
	double mu1 = minFC_full.UserState().Value(1);
	double fval_full = double(minFC_full.Fval());      
	double mu0Hat = minFC_full.UserState().Value(0);
	double mu1Hat = minFC_full.UserState().Value(1);

	uparFC_full.SetValue("mu0", idata==0 ? 0.0 : mu_true5s );
	uparFC_full.Fix("mu0");
	MnMigrad migradFC_fullFix(*likelihoodFC_full, uparFC_full, 1);
	FunctionMinimum minFC_fullFix = migradFC_fullFix(maxfcn, tolerance);	
	double mu0Fix = FCfixToTrue ? minFC_fullFix.UserState().Value(0) : minFC_full.UserState().Value(0);
	double mu1Fix = FCfixToTrue ? minFC_fullFix.UserState().Value(1) : minFC_full.UserState().Value(1);

	//cout << "Toy: mu0  " << minFC_full.UserState().Value(0) << endl;
	//cout << "Toy: mu1  " << minFC_full.UserState().Value(1) << endl;
	//cout << "Toy: mu0 fix " << mu0Fix << endl;
	//cout << "Toy: mu1 fix " << mu1Fix << endl;
	
	double fval_fullFix = double(minFC_fullFix.Fval());
	double test_stat = 2*(fval_fullFix-fval_full);
	if(idata==0)
	  tstat = test_stat;
	else
	  tstat5s = test_stat;
	uparFC_full.Release("mu0");

	//cout << "Toy: fval " << fval_full <<  " --> " << fval_fullFix << endl;

	VectorXd jtilde = VectorXd::Zero(2*nbins);
	VectorXd jhat   = VectorXd::Zero(2*nbins);
	VectorXd error_on_sum = VectorXd::Zero(nbins);
	
	MatrixXd invsqrtV  = MatrixXd::Zero(1,1);
	MatrixXd invVj = MatrixXd::Zero(2,2);
	VectorXd jhati = VectorXd::Zero(2);
	VectorXd y = VectorXd::Zero(1);
	MatrixXd X = MatrixXd::Zero(1, 2);
	for(unsigned int i=0; i<nbins; i++){
	  invsqrtV(0,0) = 1./TMath::Sqrt(y_itoy(i));
	  invVj(0,0) = lumiscale/A_itoy(i,0);
	  invVj(1,1) = lumiscale/A_itoy(i,1);       
	  jhati(0) = A_itoy(i,0);
	  jhati(1) = A_itoy(i,1);
	  jhat(2*i) = jhati(0);
	  jhat(2*i+1) = jhati(1);
	  y(0) = y_itoy(i);
	  VectorXd b = invsqrtV*y;
	  X(0,0) = 1.0 + mu0Fix;
	  X(0,1) = 1.0 + mu1Fix;
	  if(decorrelate){
	    X(0,0) = 1.0 + mu0Fix + mu1Fix;
	    X(0,1) = 1.0 - mu0Fix + mu1Fix;
	  }
	  MatrixXd D = invsqrtV*X;
	  MatrixXd B = D.transpose()*D + invVj;
	  VectorXd g = -(D.transpose()*b + invVj*jhati );
	  VectorXd jtildei = B.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(-g);
	  jtilde(2*i)   = jtildei(0);
	  jtilde(2*i+1) = jtildei(1);

	  if(computeJtildeError){
	    MatrixXd Vjtildei = ( X.transpose()*invsqrtV*invsqrtV*X + invVj ).inverse();
	    VectorXd ei = VectorXd::Zero(2);
	    ei << X(0,0), X(0,1);
	    MatrixXd error_on_sumi = ei.transpose()*Vjtildei*ei;
	    //cout << "Error on sum: " << TMath::Sqrt( 1./invVj(0,0)*X(0,0)*X(0,0)  + 1./invVj(1,1)*X(0,1)*X(0,1) ) << " --> " << TMath::Sqrt(error_on_sumi(0,0)) << endl;
	    error_on_sum(i) = TMath::Sqrt(error_on_sumi(0,0));
	  }
	  //cout << "Bin " << i << endl;
	  //cout << "\tj0: " << jtildei(0) << " --> " << jhati(0) << " (" << (jtildei(0) - jhati(0))/TMath::Sqrt(A_itoy(i,0)*lumiscale) << ")" << endl;
	  //cout << "\tj1: " << jtildei(1) << " --> " << jhati(1) << " (" << (jtildei(1) - jhati(1))/TMath::Sqrt(A_itoy(i,1)*lumiscale) << ")" << endl;
	}

	// check that solution for jtilde is correct
	if(false){
	  for(unsigned int ir=0; ir<nbins; ir++){
	    likelihoodFC_full->set_mc0(ir, jtilde(2*ir));
	    likelihoodFC_full->set_mc1(ir, jtilde(2*ir+1));
	  }
	  likelihoodFC_full->set_removeBB(true);
	  MnMigrad migrad_debug(*likelihoodFC_full, uparFC_full, 1);
	  FunctionMinimum min_debug = migrad_debug(maxfcn, tolerance);
	  double mu0_debug = min_debug.UserState().Value(0);
	  double mu1_debug = min_debug.UserState().Value(1);
	  double fval_debug = double(min_debug.Fval());
	  cout << "Mu0 from full fit: " << mu0Fix << ", mu0 after fixing jtilde: " << mu0_debug << endl;
	  likelihoodFC_full->set_removeBB(false);
	  for(unsigned int ir=0; ir<nbins; ir++){
	    likelihoodFC_full->set_mc0(ir, A_itoy(ir,0));
	    likelihoodFC_full->set_mc1(ir, A_itoy(ir,1));
	  }
	}

	// debug histograms containing the pre and post-fit distributions of templates
	TH1D* h_jtilde0_itoy = (TH1D*)h_true_0->Clone(Form("h_data%d_jtilde0_%d",idata,itoy));
	TH1D* h_jtilde1_itoy = (TH1D*)h_true_1->Clone(Form("h_data%d_jtilde1_%d",idata,itoy));
	TH1D* h_j0_itoy      = (TH1D*)h_true_0->Clone(Form("h_data%d_j0_%d",idata,itoy));
	TH1D* h_j1_itoy      = (TH1D*)h_true_1->Clone(Form("h_data%d_j1_%d",idata,itoy));
	TH1D* h_obs_itoy     = (TH1D*)h_true_0->Clone(Form("h_data%d_obs_%d",idata,itoy));
	h_jtilde0_itoy->Reset();
	h_jtilde1_itoy->Reset();
	h_j0_itoy->Reset();
	h_j1_itoy->Reset();
	h_obs_itoy->Reset();
	for(unsigned int ir=0; ir<nbins; ir++){
	  h_j0_itoy->SetBinContent(ir+1, jhat(2*ir)  );
	  h_j1_itoy->SetBinContent(ir+1, jhat(2*ir+1));
	  h_jtilde0_itoy->SetBinContent(ir+1, (1.0 + mu0Fix)*jtilde(2*ir) );
	  h_jtilde1_itoy->SetBinContent(ir+1, (1.0 + mu1Fix)*jtilde(2*ir+1) );	  
	  if(decorrelate){
	    h_jtilde0_itoy->SetBinContent(ir+1, (1.0 + mu0Fix + mu1Fix)*jtilde(2*ir) );
	    h_jtilde1_itoy->SetBinContent(ir+1, (1.0 - mu0Fix + mu1Fix)*jtilde(2*ir+1) );
	  }
	  h_obs_itoy->SetBinContent(ir+1, y_itoy(ir));
	}

	// histogram containing the sampling distribution of the test-statistic
	TString hname = "";
	if(idata==0)
	  hname = Form("h_data_teststat_%d", itoy );
	else
	  hname = Form("h_data5s_teststat_%d", itoy );	
	TH1D* h_teststat_itoy = new TH1D( hname, "", 200, 0., 20);

	for(unsigned int itoyFC=0; itoyFC<TMath::Max(ntoysFC, long(1000)); itoyFC++){

	  // fill the data for the FC toys
	  for(unsigned int ir=0; ir<nbins; ir++){      	  

	    // cheat means that we draw the data from the ideal case
	    if(doFCcheat){
	      yFC_itoy(ir) = rans[1+idata]->Poisson(  ( 1.0 + (idata==1)*mu_true5s)*A_true(ir,0) + A_true(ir,1) );
	      if(decorrelate)
		yFC_itoy(ir) = rans[1+idata]->Poisson(  ( 1.0 + (idata==1)*mu_true5s  )*A_true(ir,0) + (1.0 - (idata==1)*mu_true5s )*A_true(ir,1) );
	    }
	    // otherwise, we draw from the best fit to the data
	    else{
	      yFC_itoy(ir) = rans[1+idata]->Poisson( (1.0 + mu0Fix)*jtilde(2*ir) + (1.0 + mu1Fix)*jtilde(2*ir+1) );
	      if(decorrelate)
		yFC_itoy(ir) = rans[1+idata]->Poisson( (1.0 + mu0Fix + mu1Fix)*jtilde(2*ir) + (1.0 - mu0Fix + mu1Fix)*jtilde(2*ir+1) );

	      // Optionally, add extra noise to the data corresponding to the error on the postfit MC
	      if(computeJtildeError){
		yFC_itoy(ir) += rans[3+idata]->Gaus(0., error_on_sum(ir));
	      }
	      
	    }
	    likelihoodFC_full->set_data(ir, yFC_itoy(ir) );
	  }

	  // fit the FC toy
	  uparFC_full.Release("mu0");
	  uparFC_full.SetValue("mu0", 0.0);
	  uparFC_full.SetValue("mu1", 0.0);
	  MnMigrad migradFC_full_itoy(*likelihoodFC_full, uparFC_full, 1);
	  FunctionMinimum minFC_full_itoy = migradFC_full_itoy(maxfcn, tolerance);
	  double fval_full_itoy = double(minFC_full_itoy.Fval());
	  uparFC_full.Fix("mu0");
	  if(!doFCcheat){
	    uparFC_full.SetValue("mu0", mu0Fix );
	  }
	  else{
	    uparFC_full.SetValue("mu0", (idata==1)*mu_true5s );
	  }
	  uparFC_full.SetValue("mu1", 0.0);
	  MnMigrad migradFC_full_itoyFix(*likelihoodFC_full, uparFC_full, 1);
	  FunctionMinimum minFC_full_itoyFix = migradFC_full_itoyFix(maxfcn, tolerance);
	  double fval_full_itoyFix = double(minFC_full_itoyFix.Fval());
	  //cout << "\t FC Tot fval " << fval_full_itoy <<  " --> " << fval_full_itoyFix << endl;
	  //cout << 2*(fval_full_itoyFix - fval_full_itoy) << endl;
	  h_teststat_itoy->Fill( 2*(fval_full_itoyFix - fval_full_itoy) );
	  uparFC_full.Release("mu0");
	}

	// restore initial data
	//if(computeJtildeError) likelihoodFC_full->restore_dataErr2();
	for(unsigned int ir=0; ir<nbins; ir++){
	  likelihoodFC_full->set_data(ir, y_itoy(ir) );  
	}
	
	double q[1];
	double scale = 1.0;
	vector<double> probsum = {1.0 - TMath::Prob(1.0, 1)};
	if(h_teststat_itoy->Integral()>0){
	  h_teststat_itoy->GetQuantiles(1, q, probsum.data());
	  if(doBarlett){
	    scale = h_teststat_itoy->GetMean();
	    q[0] = 1.0*scale;
	  }
	  likelihoodFC_full->SetErrorDef(0.5*q[0]);
	  uparFC_full.Release("mu0");
	  uparFC_full.SetValue("mu0", 0.0);
	  uparFC_full.SetValue("mu1", 0.0);
	  MnMigrad migradFC_full_minos(*likelihoodFC_full, uparFC_full, 1);
	  FunctionMinimum minFC_full_minos = migradFC_full_minos(maxfcn, tolerance);
	  MnMinos minos_FCfull(*likelihoodFC_full, minFC_full_minos, 1);	
	  std::pair< double, double >  minos_err_FCfull = minos_FCfull(0);
	  if(idata==0){
	    errFCLow_poisdata_BBfull =  minos_err_FCfull.first;
	    errFCUp_poisdata_BBfull  =  minos_err_FCfull.second;
	  }
	  else{
	    errFCLow_poisdata5s_BBfull =  minos_err_FCfull.first;
	    errFCUp_poisdata5s_BBfull  =  minos_err_FCfull.second;
	  }
	}
	else{
	  q[0] = numeric_limits<double>::max();
	  errFCLow_poisdata_BBfull = -99.;
	  errFCUp_poisdata_BBfull = -99.;
	  errFCLow_poisdata5s_BBfull = -99.;
	  errFCUp_poisdata5s_BBfull = -99.;
	}

	if( idata==0 ){
	  if(ntoysFC<=100)
	    cout << "t0: " << test_stat << " < " << q[0] << " ? " << "t<q0 vs eL<0<eT: " << (test_stat<=q[0]) << ":" <<  ( 0.0 >= mu0Hat+errFCLow_poisdata_BBfull && 0.0 <= mu0Hat+errFCUp_poisdata_BBfull  )  << endl;
	  if(  mu_true >= mu0Hat+errFCLow_poisdata_BBfull && mu_true <= mu0Hat+errFCUp_poisdata_BBfull  )
	    prob_poisdata_BBfullFC += 1./ntoysFC;
	}
	else{
	  if( mu_true5s >= mu0Hat+errFCLow_poisdata5s_BBfull && mu_true5s <= mu0Hat+errFCUp_poisdata5s_BBfull)
	    prob_poisdata5s_BBfullFC += 1./ntoysFC;
	}

	/*
	if( test_stat<=q[0] ){	  
	  if(idata==0)
	    prob_poisdata_BBfullFC += 1./ntoysFC;
	  else
	    prob_poisdata5s_BBfullFC += 1./ntoysFC;
	}
	*/

	
	fout->cd();
	if(saveHistos && itoy<10){
	  h_teststat_itoy->Write();
	  h_jtilde0_itoy->Write();
	  h_jtilde1_itoy->Write();
	  h_j0_itoy->Write();
	  h_j1_itoy->Write();
	  h_obs_itoy->Write();
	}
	
	delete likelihoodFC_full;
	delete h_teststat_itoy;
	delete h_jtilde0_itoy;
	delete h_jtilde1_itoy;
	delete h_j0_itoy;
	delete h_j1_itoy;
      }      
      treeFC->Fill();
    }
  }

  
  for(unsigned int itoy=0; itoy<ntoys; itoy++){
    if(itoy%10000==0) cout << "Doing toy " << itoy << " / " << ntoys << endl;
      
    MatrixXd invC_itoy      = MatrixXd::Zero(2, 2);
    MatrixXd invC5s_itoy    = MatrixXd::Zero(2, 2);
    MatrixXd invCMC_itoy    = MatrixXd::Zero(2, 2);
    MatrixXd invC_BB_itoy   = MatrixXd::Zero(2, 2);
    MatrixXd invC5s_BB_itoy = MatrixXd::Zero(2, 2);
    MatrixXd invCMC_BB_itoy = MatrixXd::Zero(2, 2);

    Likelihood* likelihood = new Likelihood(0, nbins, true, doPoisson, profileMCNP, decorrelate);
    likelihood->SetErrorDef(0.5);
    likelihood->set_removeBB(true);
    Likelihood* likelihood_lite = new Likelihood(0, nbins, true, doPoisson, profileMCNP, decorrelate);
    likelihood_lite->SetErrorDef(0.5);
    Likelihood* likelihood_full = new Likelihood(0, nbins, false, doPoisson, profileMCNP, decorrelate);
    likelihood_full->SetErrorDef(0.5);
    Likelihood* likelihood_5s = new Likelihood(0, nbins, true, doPoisson, profileMCNP, decorrelate);
    likelihood_5s->SetErrorDef(0.5);    
    likelihood_5s->set_removeBB(true);
    Likelihood* likelihood_lite5s = new Likelihood(0, nbins, true, doPoisson, profileMCNP, decorrelate);
    likelihood_lite5s->SetErrorDef(0.5);
    Likelihood* likelihood_full5s = new Likelihood(0, nbins, false, doPoisson, profileMCNP, decorrelate);
    likelihood_full5s->SetErrorDef(0.5);
    
    for(unsigned int ir=0; ir<nbins; ir++){      
      A_itoy(ir,0) = rans[0]->Poisson( A_true(ir, 0)*lumiscale ) / lumiscale;
      A_itoy(ir,1) = rans[0]->Poisson( A_true(ir, 1)*lumiscale ) / lumiscale;
      J_itoy(ir,0) = A_itoy(ir,0);
      J_itoy(ir,1) = A_itoy(ir,1);
      if(decorrelate){
	J_itoy(ir,0) = TMath::Sqrt(2.)*(TMath::Cos(alpha)*A_itoy(ir,0) - TMath::Sin(alpha)*A_itoy(ir,1));
	J_itoy(ir,1) = TMath::Sqrt(2.)*(TMath::Sin(alpha)*A_itoy(ir,0) + TMath::Cos(alpha)*A_itoy(ir,1));
      }
      
      if(bExact) J_itoy(ir,0) = J_true(ir,0);
      if(jExact) J_itoy(ir,1) = J_true(ir,1);
      
      // "MC" template drawn from true nominal, including lumi scale
      ynom_itoy(ir) = A_itoy.row(ir).sum();
      
      // toy data drawn from true nominal 
      y_itoy(ir)   = rans[0]->Poisson( A_true.row(ir).sum() );
      
      // toy data drawn from true 5sigma 
      y5s_itoy(ir) = rans[0]->Poisson( A_true(ir,0)*(1.0 + mu_true5s) + A_true(ir,1)  );
      if(decorrelate)
	y5s_itoy(ir) = rans[0]->Poisson( A_true(ir,0)*(1.0 + mu_true5s) + A_true(ir,1)*(1.0 - mu_true5s)  );
      
      // toy data drawn from "MC" 
      yMC_itoy(ir) = rans[0]->Poisson( A_itoy.row(ir).sum() );

      invV_itoy(ir,ir)      = 1./y_itoy(ir);
      invsqrtV_itoy(ir,ir)  = 1./TMath::Sqrt(y_itoy(ir));
      invsqrtV5s_itoy(ir,ir)  = 1./TMath::Sqrt(y5s_itoy(ir));
      invsqrtV_BB_itoy(ir,ir)   = 1./TMath::Sqrt(y_itoy(ir)   + ynom_itoy(ir)/lumiscale);
      invsqrtV5s_BB_itoy(ir,ir) = 1./TMath::Sqrt(y5s_itoy(ir) + ynom_itoy(ir)/lumiscale);
      invV5s_itoy(ir,ir)    = 1./y5s_itoy(ir);
      invV_BB_itoy(ir,ir)   = 1./(y_itoy(ir)   + ynom_itoy(ir)/lumiscale);
      invV5s_BB_itoy(ir,ir) = 1./(y5s_itoy(ir) + ynom_itoy(ir)/lumiscale);
      invVMC_itoy(ir,ir)    = 1./yMC_itoy(ir);

      // this is a common piece
      MatrixXd K_ir_itoy = J_itoy.row(ir).transpose()*J_itoy.row(ir);      
      invC_itoy      += K_ir_itoy/y_itoy(ir);
      invC5s_itoy    += K_ir_itoy/y5s_itoy(ir);
      invCMC_itoy    += K_ir_itoy/yMC_itoy(ir);
      invC_BB_itoy   += K_ir_itoy/(y_itoy(ir)   + ynom_itoy(ir)/lumiscale );
      invC5s_BB_itoy += K_ir_itoy/(y5s_itoy(ir) + ynom_itoy(ir)/lumiscale );

      // numerical
      likelihood->set_mc0(ir, A_itoy(ir,0));
      likelihood->set_mc1(ir, A_itoy(ir,1));
      likelihood->set_data(ir, y_itoy(ir) );
      likelihood_lite->set_mc0(ir, A_itoy(ir,0));
      likelihood_lite->set_mc1(ir, A_itoy(ir,1));
      likelihood_lite->set_mc0err(ir, TMath::Sqrt(A_itoy(ir,0)/lumiscale));
      likelihood_lite->set_mc1err(ir, TMath::Sqrt(A_itoy(ir,1)/lumiscale));
      likelihood_lite->set_data(ir, y_itoy(ir) );
      likelihood_full->set_mc0(ir, A_itoy(ir,0));
      likelihood_full->set_mc1(ir, A_itoy(ir,1));
      likelihood_full->set_mc0err(ir, TMath::Sqrt(A_itoy(ir,0)/lumiscale));
      likelihood_full->set_mc1err(ir, TMath::Sqrt(A_itoy(ir,1)/lumiscale));
      likelihood_full->set_data(ir, y_itoy(ir) );
      likelihood_5s->set_mc0(ir, A_itoy(ir,0));
      likelihood_5s->set_mc1(ir, A_itoy(ir,1));
      likelihood_5s->set_data(ir, y5s_itoy(ir) );
      likelihood_lite5s->set_mc0(ir, A_itoy(ir,0));
      likelihood_lite5s->set_mc1(ir, A_itoy(ir,1));
      likelihood_lite5s->set_mc0err(ir, TMath::Sqrt(A_itoy(ir,0)/lumiscale));
      likelihood_lite5s->set_mc1err(ir, TMath::Sqrt(A_itoy(ir,1)/lumiscale));
      likelihood_lite5s->set_data(ir, y5s_itoy(ir) );
      likelihood_full5s->set_mc0(ir, A_itoy(ir,0));
      likelihood_full5s->set_mc1(ir, A_itoy(ir,1));
      likelihood_full5s->set_mc0err(ir, TMath::Sqrt(A_itoy(ir,0)/lumiscale));
      likelihood_full5s->set_mc1err(ir, TMath::Sqrt(A_itoy(ir,1)/lumiscale));
      likelihood_full5s->set_data(ir, y5s_itoy(ir) );
    }

    MatrixXd C_itoy      = invC_itoy.inverse();
    MatrixXd C5s_itoy    = invC5s_itoy.inverse();
    MatrixXd C_BB_itoy   = invC_BB_itoy.inverse();
    MatrixXd C5s_BB_itoy = invC5s_BB_itoy.inverse();
    MatrixXd CMC_itoy    = invCMC_itoy.inverse();
    VectorXd x_itoy      = C_itoy*J_itoy.transpose()*invV_itoy*(y_itoy - ynom_itoy);
    VectorXd x5s_itoy    = C5s_itoy*J_itoy.transpose()*invV5s_itoy*(y5s_itoy - ynom_itoy);
    VectorXd xMC_itoy    = CMC_itoy*J_itoy.transpose()*invVMC_itoy*(yMC_itoy - ynom_itoy);
    VectorXd x_BB_itoy   = C_BB_itoy*J_itoy.transpose()*invV_BB_itoy*(y_itoy - ynom_itoy);
    VectorXd x5s_BB_itoy = C5s_BB_itoy*J_itoy.transpose()*invV5s_BB_itoy*(y5s_itoy - ynom_itoy);

    if(computePropError){
      
      MatrixXd J_templ   = computePropErrorTrue ? J_true : J_itoy;
      MatrixXd W_templ   = computePropErrorTrue ? C_true   : (computePropErrorBB ?  C_BB_itoy   : C_itoy );
      MatrixXd W5s_templ = computePropErrorTrue ? C_true5s : (computePropErrorBB ?  C5s_BB_itoy : C5s_itoy);
      MatrixXd invsqrtV_templ   = computePropErrorBB ? invsqrtV_BB_itoy : invsqrtV_itoy;
      MatrixXd invsqrtV5s_templ = computePropErrorBB ? invsqrtV5s_BB_itoy : invsqrtV5s_itoy;
      MatrixXd invV_templ       = computePropErrorBB ? invV_BB_itoy : invV_itoy;
      MatrixXd invV5s_templ     = computePropErrorBB ? invV5s_BB_itoy : invV5s_itoy;
      
      VectorXd x_templ   = W_templ*J_templ.transpose()*invV_templ*y_itoy;
      VectorXd x5s_templ = W5s_templ*J_templ.transpose()*invV5s_templ*y5s_itoy;

      MatrixXd A_templ = invsqrtV_templ*J_templ;
      //VectorXd r_templ = invsqrtV_templ*(y_itoy - ynom_itoy);
      VectorXd r_templ = invsqrtV_templ*y_itoy;
      VectorXd rstar_templ = r_templ - A_templ*x_templ;

      MatrixXd A5s_templ = invsqrtV5s_templ*J_templ;
      //VectorXd r5s_templ = invsqrtV5s_templ*(y5s_itoy - ynom_itoy);
      VectorXd r5s_templ = invsqrtV5s_templ*y5s_itoy;
      VectorXd rstar5s_templ = r5s_templ - A5s_templ*x5s_templ;

      double var_templ = 0.;
      double var5s_templ = 0.;
      for( unsigned int it_i = 0; it_i<y_itoy.size(); it_i++ ){
	for( unsigned int it_j = 0; it_j<x_templ.size(); it_j++ ){
	  VectorXd vij = VectorXd::Zero( x_templ.size() );
	  VectorXd v5sij = VectorXd::Zero( x_templ.size() );
	  for(unsigned int it_k = 0; it_k < x_templ.size(); it_k++){
	    vij(it_k) = - x_templ(it_j)*A_templ(it_i, it_k);
	    if( it_k == it_j) vij(it_k) += rstar_templ(it_i);  
	    v5sij(it_k) = - x5s_templ(it_j)*A5s_templ(it_i, it_k);
	    if( it_k == it_j) v5sij(it_k) += rstar5s_templ(it_i);  
	  }
	  double gij = (W_templ*vij)(0);
	  double varij = decorrelate ? A_itoy(it_i,0) + A_itoy(it_i,1) :  A_itoy(it_i,it_j);
	  varij /= (lumiscale*(computePropErrorBB ? y_itoy(it_i) + ynom_itoy(it_i)/lumiscale : y_itoy(it_i) ) );
	  var_templ += gij*gij*varij;

	  double g5sij = (W_templ*v5sij)(0);
	  double var5sij = decorrelate ? A_itoy(it_i,0) + A_itoy(it_i,1) :  A_itoy(it_i,it_j);
	  var5sij /= (lumiscale*(computePropErrorBB ? y5s_itoy(it_i) + ynom_itoy(it_i)/lumiscale : y5s_itoy(it_i)));
	  var5s_templ += g5sij*g5sij*var5sij;
	}
      }
      //cout << "Stat.: " << TMath::Sqrt(W_templ(0,0)) << ", from templates: " << TMath::Sqrt(var_templ) << endl;
      err_propdata   = TMath::Sqrt( W_templ(0,0) + var_templ );
      //cout << "\t" << err_propdata << endl;
      err_propdata5s = TMath::Sqrt( W5s_templ(0,0) + var5s_templ );
      //err_propdata = TMath::Sqrt(var_templ);
      //err_propdata5s = TMath::Sqrt(var5s_templ);
    }

    mu_data      = x_itoy(0);
    mu_data5s    = x5s_itoy(0);
    mu_data_BB   = x_BB_itoy(0);
    mu_data5s_BB = x5s_BB_itoy(0);
    mu_mc        = xMC_itoy(0); 
    err_data      = TMath::Sqrt(C_itoy(0,0));
    err_data5s    = TMath::Sqrt(C5s_itoy(0,0));
    err_data_BB   = TMath::Sqrt(C_BB_itoy(0,0));
    err_data5s_BB = TMath::Sqrt(C5s_BB_itoy(0,0));
    err_mc        = TMath::Sqrt(CMC_itoy(0,0));

    //cout << x5s_itoy(0) << " +/- " << err_data5s << " : " << x5s_itoy(1) << " +/- " << TMath::Sqrt(C5s_itoy(1,1)) << endl;
    
    MnUserParameters upar_lite;
    //upar_lite.Add("mu0", x_BB_itoy(0), TMath::Sqrt(C_BB_itoy(0,0)), x_BB_itoy(0) - 5.0*TMath::Sqrt(C_BB_itoy(0,0)), x_BB_itoy(0) + 5.0*TMath::Sqrt(C_BB_itoy(0,0)) );
    //upar_lite.Add("mu1", x_BB_itoy(1), TMath::Sqrt(C_BB_itoy(1,1)), x_BB_itoy(1) - 5.0*TMath::Sqrt(C_BB_itoy(1,1)), x_BB_itoy(1) + 5.0*TMath::Sqrt(C_BB_itoy(1,1)) );
    upar_lite.Add("mu0", 0., 0.01);
    upar_lite.Add("mu1", 0., 0.01);
    if(profileMCNP){
      for (int i=0; i<nbins; i++){
	upar_lite.Add(Form("BBLite%d",i), 0.0, 0.001);
      }
    }
   

    // Gaussian
    if(verbose) cout << "\tMigrad..." << endl;
    MnMigrad migrad(*likelihood, upar_lite, 1);
    FunctionMinimum min = migrad(maxfcn, tolerance);
    double edm = double(min.Edm());
    double fmin = double(min.Fval());
    if(verbose) cout << "\tHesse..." << endl;
    MnHesse hesse(1);
    hesse(*likelihood, min);       
    for(unsigned int i = 0 ; i<npars_lite; i++){    
      for(unsigned int j = 0 ; j<npars_lite; j++){
	Vmin(i,j) = i>j ?
	  min.UserState().Covariance().Data()[j+ i*(i+1)/2] :
	  min.UserState().Covariance().Data()[i+ j*(j+1)/2];;
      }
    }
    for(unsigned int i = 0 ; i<npars_lite; i++){
      xmin(i)     = min.UserState().Value(i) ;
      xmin_Err(i) = min.UserState().Error(i) ;
    }

    // Gaussian 5s
    if(verbose) cout << "\tMigrad..." << endl;
    MnMigrad migrad_5s(*likelihood_5s, upar_lite, 1);
    FunctionMinimum min_5s = migrad_5s(maxfcn, tolerance);
    double edm_5s = double(min_5s.Edm());
    double fmin_5s = double(min_5s.Fval());
    if(verbose) cout << "\tHesse..." << endl;
    MnHesse hesse_5s(1);
    hesse(*likelihood_5s, min_5s);       
    for(unsigned int i = 0 ; i<npars_lite; i++){    
      for(unsigned int j = 0 ; j<npars_lite; j++){
	Vmin_5s(i,j) = i>j ?
	  min_5s.UserState().Covariance().Data()[j+ i*(i+1)/2] :
	  min_5s.UserState().Covariance().Data()[i+ j*(j+1)/2];;
      }
    }
    for(unsigned int i = 0 ; i<npars_lite; i++){
      xmin_5s(i)     = min_5s.UserState().Value(i) ;
      xmin_5sErr(i)  = min_5s.UserState().Error(i) ;
    }

    // Gaussian + BBlite
    if(verbose) cout << "\tMigrad..." << endl;
    MnMigrad migrad_lite(*likelihood_lite, upar_lite, 1);
    FunctionMinimum min_lite = migrad_lite(maxfcn, tolerance);
    double edm_lite = double(min_lite.Edm());
    double fmin_lite = double(min_lite.Fval());
    if(verbose) cout << "\tHesse..." << endl;
    MnHesse hesse_lite(1);
    hesse_lite(*likelihood_lite, min_lite);       
    for(unsigned int i = 0 ; i<npars_lite; i++){    
      for(unsigned int j = 0 ; j<npars_lite; j++){
	Vmin_lite(i,j) = i>j ?
	  min_lite.UserState().Covariance().Data()[j+ i*(i+1)/2] :
	  min_lite.UserState().Covariance().Data()[i+ j*(j+1)/2];;
      }
    }
    for(unsigned int i = 0 ; i<npars_lite; i++){
      xmin_lite(i)    = min_lite.UserState().Value(i) ;
      xmin_liteErr(i) = min_lite.UserState().Error(i) ;
      //xmin_liteErr(i) = TMath::Sqrt(Vmin_lite(i,i)) ;
    }

    // Gaussian + BBlite 5s
    if(verbose) cout << "\tMigrad..." << endl;
    MnMigrad migrad_lite5s(*likelihood_lite5s, upar_lite, 1);
    FunctionMinimum min_lite5s = migrad_lite5s(maxfcn, tolerance);
    double edm_lite5s = double(min_lite5s.Edm());
    double fmin_lite5s = double(min_lite5s.Fval());
    if(verbose) cout << "\tHesse..." << endl;
    MnHesse hesse_lite5s(1);
    hesse_lite5s(*likelihood_lite5s, min_lite5s);       
    for(unsigned int i = 0 ; i<npars_lite; i++){    
      for(unsigned int j = 0 ; j<npars_lite; j++){
	Vmin_lite5s(i,j) = i>j ?
	  min_lite5s.UserState().Covariance().Data()[j+ i*(i+1)/2] :
	  min_lite5s.UserState().Covariance().Data()[i+ j*(j+1)/2];;
      }
    }
    for(unsigned int i = 0 ; i<npars_lite; i++){
      xmin_lite5s(i)    = min_lite5s.UserState().Value(i) ;
      xmin_lite5sErr(i) = min_lite5s.UserState().Error(i) ;
      //xmin_lite5sErr(i) = TMath::Sqrt(Vmin_lite5s(i,i)) ;
    }    
    
    if(verbose){
      cout << ">>>>> Lite:" << endl;
      cout << "Edm: " << min_lite.Edm() << std::endl;
      cout << "Val: " << min_lite.Fval() << std::endl;
      cout << "min is valid: " << min_lite.IsValid() << std::endl;
      cout << "HesseFailed: " << min_lite.HesseFailed() << std::endl;
      cout << "HasCovariance: " << min_lite.HasCovariance() << std::endl;
      cout << "HasValidCovariance: " << min_lite.HasValidCovariance() << std::endl;
      cout << "HasValidParameters: " << min_lite.HasValidParameters() << std::endl;
      cout << "IsAboveMaxEdm: " << min_lite.IsAboveMaxEdm() << std::endl;
      cout << "HasReachedCallLimit: " << min_lite.HasReachedCallLimit() << std::endl;
      cout << "HasAccurateCovar: " << min_lite.HasAccurateCovar() << std::endl;
      cout << "HasPosDefCovar : " << min_lite.HasPosDefCovar() << std::endl;
      cout << "HasMadePosDefCovar : " << min_lite.HasMadePosDefCovar() << std::endl;
    }

    MnUserParameters upar_full;
    //upar_full.Add("mu0", x_BB_itoy(0), TMath::Sqrt(C_BB_itoy(0,0)), x_BB_itoy(0) - 5.0*TMath::Sqrt(C_BB_itoy(0,0)), x_BB_itoy(0) + 5.0*TMath::Sqrt(C_BB_itoy(0,0)) );
    //upar_full.Add("mu1", x_BB_itoy(1), TMath::Sqrt(C_BB_itoy(1,1)), x_BB_itoy(1) - 5.0*TMath::Sqrt(C_BB_itoy(1,1)), x_BB_itoy(1) + 5.0*TMath::Sqrt(C_BB_itoy(1,1)) );
    upar_full.Add("mu0", 0.0, 0.01);
    upar_full.Add("mu1", 0.0, 0.01);
    if(profileMCNP){
      for (int i=0; i<nbins; i++){
	double exp_error = TMath::Sqrt( 1.0/(lumiscale*A_itoy(i,0) ) );
	upar_full.Add(Form("BB0Full%d",i), 0.0, 0.001 );
      }
      for (int i=0; i<nbins; i++){
	double exp_error = TMath::Sqrt( 1.0/(lumiscale*A_itoy(i,1) ) );
	upar_full.Add(Form("BB1Full%d",i), 0.0, 0.001);
      }
    }
    
    MnMigrad migrad_full(*likelihood_full, upar_full, 1);
    MnMigrad migrad_full5s(*likelihood_full5s, upar_full, 1);

    if(verbose) cout << "\tMigrad..." << endl;
    FunctionMinimum min_full = migrad_full(maxfcn, tolerance);
    double edm_full = double(min_full.Edm());
    double fmin_full = double(min_full.Fval());
    if(verbose) cout << "\tHesse..." << endl;
    MnHesse hesse_full(1);
    hesse_full(*likelihood_full, min_full);       
    for(unsigned int i = 0 ; i<npars_full; i++){    
      for(unsigned int j = 0 ; j<npars_full; j++){
	Vmin_full(i,j) = i>j ?
	  min_full.UserState().Covariance().Data()[j+ i*(i+1)/2] :
	  min_full.UserState().Covariance().Data()[i+ j*(j+1)/2];;
      }
    }
    for(unsigned int i = 0 ; i<npars_full; i++){
      xmin_full(i)    = min_full.UserState().Value(i) ;
      xmin_fullErr(i) = min_full.UserState().Error(i) ;
    }

    MnMinos minos_full(*likelihood_full, min_full, 1);
    std::pair< double, double >  minos_err_full = minos_full(0);
    errPLRLow_poisdata_BBfull =  minos_err_full.first;
    errPLRUp_poisdata_BBfull  =  minos_err_full.second;

    if(verbose) cout << "\tMigrad..." << endl;
    FunctionMinimum min_full5s = migrad_full5s(maxfcn, tolerance);
    double edm_full5s = double(min_full5s.Edm());
    double fmin_full5s = double(min_full5s.Fval());
    if(verbose) cout << "\tHesse..." << endl;
    MnHesse hesse_full5s(1);
    hesse_full5s(*likelihood_full5s, min_full5s);       
    for(unsigned int i = 0 ; i<npars_full; i++){    
      for(unsigned int j = 0 ; j<npars_full; j++){
	Vmin_full5s(i,j) = i>j ?
	  min_full5s.UserState().Covariance().Data()[j+ i*(i+1)/2] :
	  min_full5s.UserState().Covariance().Data()[i+ j*(j+1)/2];;
      }
    }
    for(unsigned int i = 0 ; i<npars_full; i++){
      xmin_full5s(i)    = min_full5s.UserState().Value(i) ;
      xmin_full5sErr(i) = min_full5s.UserState().Error(i) ;
    }

    MnMinos minos_full5s(*likelihood_full5s, min_full5s, 1);
    std::pair< double, double >  minos_err_full5s = minos_full5s(0);
    errPLRLow_poisdata5s_BBfull =  minos_err_full5s.first;
    errPLRUp_poisdata5s_BBfull  =  minos_err_full5s.second;

    
    if(verbose){
      cout << ">>>>> Full:" << endl;
      cout << "Edm: " << min_full.Edm() << std::endl;
      cout << "Val: " << min_full.Fval() << std::endl;
      cout << "min is valid: " << min_full.IsValid() << std::endl;
      cout << "HesseFailed: " << min_full.HesseFailed() << std::endl;
      cout << "HasCovariance: " << min_full.HasCovariance() << std::endl;
      cout << "HasValidCovariance: " << min_full.HasValidCovariance() << std::endl;
      cout << "HasValidParameters: " << min_full.HasValidParameters() << std::endl;
      cout << "IsAboveMaxEdm: " << min_full.IsAboveMaxEdm() << std::endl;
      cout << "HasReachedCallLimit: " << min_full.HasReachedCallLimit() << std::endl;
      cout << "HasAccurateCovar: " << min_full.HasAccurateCovar() << std::endl;
      cout << "HasPosDefCovar : " << min_full.HasPosDefCovar() << std::endl;
      cout << "HasMadePosDefCovar : " << min_full.HasMadePosDefCovar() << std::endl;    
    }

    mu_poisdata           = xmin(0);
    err_poisdata          = xmin_Err(0);
    mu_poisdata_BB        = xmin_lite(0);
    err_poisdata_BB       = xmin_liteErr(0);
    mu_poisdata_BBfull    = xmin_full(0);
    err_poisdata_BBfull   = xmin_fullErr(0);
    mu_poisdata5s         = xmin_5s(0);
    err_poisdata5s        = xmin_5sErr(0);
    mu_poisdata5s_BB      = xmin_lite5s(0);
    err_poisdata5s_BB     = xmin_lite5sErr(0);
    mu_poisdata5s_BBfull  = xmin_full5s(0);
    err_poisdata5s_BBfull = xmin_full5sErr(0);

    err_adhocdata   = err_true * TMath::Sqrt( 1 + 1./lumiscale ) * (err_poisdata_BBfull / err_poisdata_BB);
    err_adhocdata5s = err_true5s * TMath::Sqrt( 1 + 1./lumiscale ) * (err_poisdata5s_BBfull / err_poisdata5s_BB);

    if( TMath::Abs(mu_poisdata_BBfull-mu_true)/err_adhocdata <= 1.0 )
      prob_adhocdata += 1./ntoys;
    if( TMath::Abs(mu_poisdata5s_BBfull-mu_true5s)/err_adhocdata5s <= 1.0 )
      prob_adhocdata5s += 1./ntoys;
    if( TMath::Abs(mu_data-mu_true)/err_propdata <= 1.0 )
      prob_propdata += 1./ntoys;
    if( TMath::Abs(mu_data5s-mu_true5s)/err_propdata5s <= 1.0 )
      prob_propdata5s += 1./ntoys;
    
    
    if( TMath::Abs(mu_data-mu_true)/err_data <= 1.0 )
      prob_data += 1./ntoys;
    if( TMath::Abs(mu_data5s-mu_true5s)/err_data5s <= 1.0 )
      prob_data5s += 1./ntoys;
    if( TMath::Abs(mu_data_BB-mu_true)/err_data_BB <= 1.0 )
      prob_data_BB += 1./ntoys;
    if( TMath::Abs(mu_data5s_BB- mu_true5s)/err_data5s_BB <= 1.0 )
      prob_data5s_BB += 1./ntoys;
    if( TMath::Abs(mu_mc-mu_true)/err_mc <= 1.0 )
      prob_mc += 1./ntoys;
    if( TMath::Abs(mu_poisdata-mu_true)/err_poisdata <= 1.0 )
      prob_poisdata += 1./ntoys;
    if( TMath::Abs(mu_poisdata_BB-mu_true)/err_poisdata_BB <= 1.0 )
      prob_poisdata_BB += 1./ntoys;
    if( TMath::Abs(mu_poisdata_BBfull-mu_true)/err_poisdata_BBfull <= 1.0 )
      prob_poisdata_BBfull += 1./ntoys;
    if( TMath::Abs(mu_poisdata5s- mu_true5s)/err_poisdata5s <= 1.0 )
      prob_poisdata5s += 1./ntoys;
    if( TMath::Abs(mu_poisdata5s_BB- mu_true5s)/err_poisdata5s_BB <= 1.0 )
      prob_poisdata5s_BB += 1./ntoys;
    if( TMath::Abs(mu_poisdata5s_BBfull- mu_true5s)/err_poisdata5s_BBfull <= 1.0 )
      prob_poisdata5s_BBfull += 1./ntoys;
    if( mu_poisdata_BBfull >= (mu_true + errPLRLow_poisdata_BBfull) && mu_poisdata_BBfull <= (mu_true + errPLRUp_poisdata_BBfull)  )
      prob_poisdata_BBfullPLR += 1./ntoys;
    if( mu_poisdata5s_BBfull >= (mu_true5s + errPLRLow_poisdata5s_BBfull) && mu_poisdata5s_BBfull <= (mu_true5s + errPLRUp_poisdata5s_BBfull) )
      prob_poisdata5s_BBfullPLR += 1./ntoys;
    
    if(verbose){
      cout << "Gaus:          " << mu_data << " +/- " << err_data << endl;
      cout << "Gaus + bblite: " << mu_data_BB << " +/- " << err_data_BB << endl;
      cout << "Pois + bblite: " << xmin_lite(0) << " +/- " << xmin_liteErr(0) << endl;
      cout << "Pois + mcstat: " << xmin_full(0) << " +/- " << xmin_fullErr(0) << endl;
    }
    
    delete likelihood_lite;
    delete likelihood_full;
    
    tree->Fill();
  }

  TTree* treesum = new TTree("treesum","");

  int n, nFC;
  double asym_corr, asym_cond;
  double asym_err, asym_derr, asym_med, asym_cov;
  double asym5s_err;
  double data_err, data_derr, data_med, data_cov;
  double data5s_err, data5s_derr, data5s_med, data5s_cov;
  double dataBB_err, dataBB_derr, dataBB_med, dataBB_cov;
  double data5sBB_err, data5sBB_derr, data5sBB_med, data5sBB_cov;
  double dataPois_err, dataPois_derr, dataPois_med, dataPois_cov;
  double dataPoisBB_err, dataPoisBB_derr, dataPoisBB_med, dataPoisBB_cov;
  double data5sPois_err, data5sPois_derr, data5sPois_med, data5sPois_cov;
  double data5sPoisBB_err, data5sPoisBB_derr, data5sPoisBB_med, data5sPoisBB_cov;
  double dataPoisBBfull_err, dataPoisBBfull_derr, dataPoisBBfull_med, dataPoisBBfull_cov;
  double data5sPoisBBfull_err, data5sPoisBBfull_derr, data5sPoisBBfull_med, data5sPoisBBfull_cov;
  double dataPoisBBfullPLR_err, dataPoisBBfullPLR_derr, dataPoisBBfullPLR_med, dataPoisBBfullPLR_cov;
  double data5sPoisBBfullPLR_err, data5sPoisBBfullPLR_derr, data5sPoisBBfullPLR_med, data5sPoisBBfullPLR_cov;
  double dataPoisBBfullFC_err, dataPoisBBfullFC_derr, dataPoisBBfullFC_med, dataPoisBBfullFC_cov;
  double data5sPoisBBfullFC_err, data5sPoisBBfullFC_derr, data5sPoisBBfullFC_med, data5sPoisBBfullFC_cov;
  double adhocdata_err, adhocdata_derr, adhocdata_med, adhocdata_cov;
  double adhocdata5s_err, adhocdata5s_derr, adhocdata5s_med, adhocdata5s_cov;
  double propdata_err, propdata_derr, propdata_med, propdata_cov;
  double propdata5s_err, propdata5s_derr, propdata5s_med, propdata5s_cov;
  
  treesum->Branch("n",  &n,  "n/I");
  treesum->Branch("nFC",&nFC,  "nFC/I");
  treesum->Branch("asym_corr",  &asym_corr,  "asym_corr/D");
  treesum->Branch("asym_cond",  &asym_cond,  "asym_cond/D");
  treesum->Branch("asym_err",  &asym_err,  "asym_err/D");
  treesum->Branch("asym5s_err", &asym5s_err,  "asym5s_err/D");
  treesum->Branch("asym_derr", &asym_derr, "asym_derr/D");
  treesum->Branch("asym_med",  &asym_med,  "asym_med/D");
  treesum->Branch("asym_cov",  &asym_cov,  "asym_cov/D");
  treesum->Branch("data_err",  &data_err,  "data_err/D");
  treesum->Branch("data_derr", &data_derr, "data_derr/D");
  treesum->Branch("data_med",  &data_med,  "data_med/D");
  treesum->Branch("data_cov",  &data_cov,  "data_cov/D");
  treesum->Branch("data5s_err",  &data5s_err,  "data5s_err/D");
  treesum->Branch("data5s_derr", &data5s_derr, "data5s_derr/D");
  treesum->Branch("data5s_med",  &data5s_med,  "data5s_med/D");
  treesum->Branch("data5s_cov",  &data5s_cov,  "data5s_cov/D");
  treesum->Branch("dataBB_err",  &dataBB_err,  "dataBB_err/D");
  treesum->Branch("dataBB_derr", &dataBB_derr, "dataBB_derr/D");
  treesum->Branch("dataBB_med",  &dataBB_med,  "dataBB_med/D");
  treesum->Branch("dataBB_cov",  &dataBB_cov,  "dataBB_cov/D");
  treesum->Branch("data5sBB_err",  &data5sBB_err,  "data5sBB_err/D");
  treesum->Branch("data5sBB_derr", &data5sBB_derr, "data5sBB_derr/D");
  treesum->Branch("data5sBB_med",  &data5sBB_med,  "data5sBB_med/D");
  treesum->Branch("data5sBB_cov",  &data5sBB_cov,  "data5sBB_cov/D");
  treesum->Branch("dataPois_err",  &dataPois_err,  "dataPois_err/D");
  treesum->Branch("dataPois_derr", &dataPois_derr, "dataPois_derr/D");
  treesum->Branch("dataPois_med",  &dataPois_med,  "dataPois_med/D");
  treesum->Branch("dataPois_cov",  &dataPois_cov,  "dataPois_cov/D");
  treesum->Branch("dataPoisBB_err",  &dataPoisBB_err,  "dataPoisBB_err/D");
  treesum->Branch("dataPoisBB_derr", &dataPoisBB_derr, "dataPoisBB_derr/D");
  treesum->Branch("dataPoisBB_med",  &dataPoisBB_med,  "dataPoisBB_med/D");
  treesum->Branch("dataPoisBB_cov",  &dataPoisBB_cov,  "dataPoisBB_cov/D");
  treesum->Branch("data5sPois_err",  &data5sPois_err,  "data5sPois_err/D");
  treesum->Branch("data5sPois_derr", &data5sPois_derr, "data5sPois_derr/D");
  treesum->Branch("data5sPois_med",  &data5sPois_med,  "data5sPois_med/D");
  treesum->Branch("data5sPois_cov",  &data5sPois_cov,  "data5sPois_cov/D");  
  treesum->Branch("data5sPoisBB_err",  &data5sPoisBB_err,  "data5sPoisBB_err/D");
  treesum->Branch("data5sPoisBB_derr", &data5sPoisBB_derr, "data5sPoisBB_derr/D");
  treesum->Branch("data5sPoisBB_med",  &data5sPoisBB_med,  "data5sPoisBB_med/D");
  treesum->Branch("data5sPoisBB_cov",  &data5sPoisBB_cov,  "data5sPoisBB_cov/D");
  treesum->Branch("dataPoisBBfull_err",  &dataPoisBBfull_err,  "dataPoisBBfull_err/D");
  treesum->Branch("dataPoisBBfull_derr", &dataPoisBBfull_derr, "dataPoisBBfull_derr/D");
  treesum->Branch("dataPoisBBfull_med",  &dataPoisBBfull_med,  "dataPoisBBfull_med/D");
  treesum->Branch("dataPoisBBfull_cov",  &dataPoisBBfull_cov,  "dataPoisBBfull_cov/D");
  treesum->Branch("data5sPoisBBfull_err",  &data5sPoisBBfull_err,  "data5sPoisBBfull_err/D");
  treesum->Branch("data5sPoisBBfull_derr", &data5sPoisBBfull_derr, "data5sPoisBBfull_derr/D");
  treesum->Branch("data5sPoisBBfull_med",  &data5sPoisBBfull_med,  "data5sPoisBBfull_med/D");
  treesum->Branch("data5sPoisBBfull_cov",  &data5sPoisBBfull_cov,  "data5sPoisBBfull_cov/D");
  treesum->Branch("dataPoisBBfullPLR_err",  &dataPoisBBfullPLR_err,  "dataPoisBBfullPLR_err/D");
  treesum->Branch("dataPoisBBfullPLR_derr", &dataPoisBBfullPLR_derr, "dataPoisBBfullPLR_derr/D");
  treesum->Branch("dataPoisBBfullPLR_med",  &dataPoisBBfullPLR_med,  "dataPoisBBfullPLR_med/D");
  treesum->Branch("dataPoisBBfullPLR_cov",  &dataPoisBBfullPLR_cov,  "dataPoisBBfullPLR_cov/D");
  treesum->Branch("data5sPoisBBfullPLR_err",  &data5sPoisBBfullPLR_err,  "data5sPoisBBfullPLR_err/D");
  treesum->Branch("data5sPoisBBfullPLR_derr", &data5sPoisBBfullPLR_derr, "data5sPoisBBfullPLR_derr/D");
  treesum->Branch("data5sPoisBBfullPLR_med",  &data5sPoisBBfullPLR_med,  "data5sPoisBBfullPLR_med/D");
  treesum->Branch("data5sPoisBBfullPLR_cov",  &data5sPoisBBfullPLR_cov,  "data5sPoisBBfullPLR_cov/D");
  treesum->Branch("dataPoisBBfullFC_err",  &dataPoisBBfullFC_err,  "dataPoisBBfullFC_err/D");
  treesum->Branch("dataPoisBBfullFC_derr", &dataPoisBBfullFC_derr, "dataPoisBBfullFC_derr/D");
  treesum->Branch("dataPoisBBfullFC_med",  &dataPoisBBfullFC_med,  "dataPoisBBfullFC_med/D");
  treesum->Branch("dataPoisBBfullFC_cov",  &dataPoisBBfullFC_cov,  "dataPoisBBfullFC_cov/D");
  treesum->Branch("data5sPoisBBfullFC_err",  &data5sPoisBBfullFC_err,  "data5sPoisBBfullFC_err/D");
  treesum->Branch("data5sPoisBBfullFC_derr", &data5sPoisBBfullFC_derr, "data5sPoisBBfullFC_derr/D");
  treesum->Branch("data5sPoisBBfullFC_med",  &data5sPoisBBfullFC_med,  "data5sPoisBBfullFC_med/D");
  treesum->Branch("data5sPoisBBfullFC_cov",  &data5sPoisBBfullFC_cov,  "data5sPoisBBfullFC_cov/D");
  treesum->Branch("adhocdata_err",  &adhocdata_err,  "adhocdata_err/D");
  treesum->Branch("adhocdata_derr", &adhocdata_derr, "adhocdata_derr/D");
  treesum->Branch("adhocdata_med",  &adhocdata_med,  "adhocdata_med/D");
  treesum->Branch("adhocdata_cov",  &adhocdata_cov,  "adhocdata_cov/D");
  treesum->Branch("adhocdata5s_err",  &adhocdata5s_err,  "adhocdata5s_err/D");
  treesum->Branch("adhocdata5s_derr", &adhocdata5s_derr, "adhocdata5s_derr/D");
  treesum->Branch("adhocdata5s_med",  &adhocdata5s_med,  "adhocdata5s_med/D");
  treesum->Branch("adhocdata5s_cov",  &adhocdata5s_cov,  "adhocdata5s_cov/D");
  treesum->Branch("propdata_err",  &propdata_err,  "propdata_err/D");
  treesum->Branch("propdata_derr", &propdata_derr, "propdata_derr/D");
  treesum->Branch("propdata_med",  &propdata_med,  "propdata_med/D");
  treesum->Branch("propdata_cov",  &propdata_cov,  "propdata_cov/D");
  treesum->Branch("propdata5s_err",  &propdata5s_err,  "propdata5s_err/D");
  treesum->Branch("propdata5s_derr", &propdata5s_derr, "propdata5s_derr/D");
  treesum->Branch("propdata5s_med",  &propdata5s_med,  "propdata5s_med/D");
  treesum->Branch("propdata5s_cov",  &propdata5s_cov,  "propdata5s_cov/D");

  n = ntoys;
  nFC = ntoysFC;
  asym_corr = rho_true;
  asym_cond = condition_true;
  asym_err = TMath::Sqrt(C_true(0,0));
  asym5s_err = TMath::Sqrt(C_true5s(0,0));
  asym_derr = 0.0;
  asym_med = asym_err;
  asym_cov = 1.0 - TMath::Prob(1.0, 1);
  
  TH1D* haux = new TH1D("haux","", 500, 0., 0.5);

  tree->Draw("err_mc>>haux");
  cout << "Asympt. err:            " << TMath::Sqrt(C_true(0,0)) << ", coverage = " << 1.0 - TMath::Prob(1.0, 1)  << endl;
  cout << "Asympt. err 5s:         " << TMath::Sqrt(C_true5s(0,0)) << ", coverage = " << 1.0 - TMath::Prob(1.0, 1)  << endl;
  cout << "MC:                     " << prob_mc << " +/- " << TMath::Sqrt(prob_mc*(1-prob_mc)/ntoys) <<  " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();

  tree->Draw("err_data>>haux");  
  data_err  = haux->GetMean();
  data_derr = haux->GetMeanError();
  data_med  = get_quantile(haux);
  data_cov  = prob_data;
  cout << "Data:                   " << prob_data << " +/- " << TMath::Sqrt(prob_data*(1-prob_data)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();

  tree->Draw("err_data5s>>haux");  
  data5s_err  = haux->GetMean();
  data5s_derr = haux->GetMeanError();
  data5s_med  = get_quantile(haux);
  data5s_cov  = prob_data5s;
  cout << "Data 5s:                " << prob_data5s << " +/- " << TMath::Sqrt(prob_data5s*(1-prob_data5s)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();

  tree->Draw("err_data_BB>>haux");  
  dataBB_err  = haux->GetMean();
  dataBB_derr = haux->GetMeanError();
  dataBB_med  = get_quantile(haux); 
  dataBB_cov  = prob_data_BB;
  cout << "Data BB:                " << prob_data_BB << " +/- " << TMath::Sqrt(prob_data_BB*(1-prob_data_BB)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();

  tree->Draw("err_data5s_BB>>haux");  
  data5sBB_err  = haux->GetMean();
  data5sBB_derr = haux->GetMeanError();
  data5sBB_med  = get_quantile(haux); 
  data5sBB_cov  = prob_data5s_BB;
  cout << "Data 5s BB:             " << prob_data5s_BB << " +/- " << TMath::Sqrt(prob_data5s_BB*(1-prob_data5s_BB)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();

  tree->Draw("err_poisdata>>haux");  
  dataPois_err  = haux->GetMean();
  dataPois_derr = haux->GetMeanError();
  dataPois_med  = get_quantile(haux); 
  dataPois_cov  = prob_poisdata;
  cout << "Data Pois               " << prob_poisdata << " +/- " << TMath::Sqrt(prob_poisdata*(1-prob_poisdata)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();

  tree->Draw("err_poisdata5s>>haux");  
  data5sPois_err  = haux->GetMean();
  data5sPois_derr = haux->GetMeanError();
  data5sPois_med  = get_quantile(haux); 
  data5sPois_cov  = prob_poisdata5s;
  cout << "Data5s Pois             " << prob_poisdata5s << " +/- " << TMath::Sqrt(prob_poisdata5s*(1-prob_poisdata5s)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();

  tree->Draw("err_poisdata_BB>>haux");  
  dataPoisBB_err  = haux->GetMean();
  dataPoisBB_derr = haux->GetMeanError();
  dataPoisBB_med  = get_quantile(haux); 
  dataPoisBB_cov  = prob_poisdata_BB;
  cout << "Data Pois BB:           " << prob_poisdata_BB << " +/- " << TMath::Sqrt(prob_poisdata_BB*(1-prob_poisdata_BB)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();

  tree->Draw("err_poisdata5s_BB>>haux");  
  data5sPoisBB_err  = haux->GetMean();
  data5sPoisBB_derr = haux->GetMeanError();
  data5sPoisBB_med  = get_quantile(haux); 
  data5sPoisBB_cov  = prob_poisdata5s_BB;
  cout << "Data5s Pois BB:         " << prob_poisdata5s_BB << " +/- " << TMath::Sqrt(prob_poisdata5s_BB*(1-prob_poisdata5s_BB)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();

  tree->Draw("err_poisdata_BBfull>>haux");  
  dataPoisBBfull_err  = haux->GetMean();
  dataPoisBBfull_derr = haux->GetMeanError();
  dataPoisBBfull_med  = get_quantile(haux); 
  dataPoisBBfull_cov  = prob_poisdata_BBfull;
  cout << "Data Pois BBFull:       " << prob_poisdata_BBfull << " +/- " << TMath::Sqrt(prob_poisdata_BBfull*(1-prob_poisdata_BBfull)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();
  
  tree->Draw("err_poisdata5s_BBfull>>haux");  
  data5sPoisBBfull_err  = haux->GetMean();
  data5sPoisBBfull_derr = haux->GetMeanError();
  data5sPoisBBfull_med  = get_quantile(haux); 
  data5sPoisBBfull_cov  = prob_poisdata5s_BBfull;
  cout << "Data5s Pois BBFull:     " << prob_poisdata5s_BBfull << " +/- " << TMath::Sqrt(prob_poisdata5s_BBfull*(1-prob_poisdata5s_BBfull)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();

  tree->Draw("errPLRUp_poisdata_BBfull>>haux");  
  dataPoisBBfullPLR_err  = haux->GetMean();
  dataPoisBBfullPLR_derr = haux->GetMeanError();
  dataPoisBBfullPLR_med  = get_quantile(haux); 
  dataPoisBBfullPLR_cov  = prob_poisdata_BBfullPLR;
  cout << "Data PLR BBFull:        " << prob_poisdata_BBfullPLR << " +/- " << TMath::Sqrt(prob_poisdata_BBfullPLR*(1-prob_poisdata_BBfullPLR)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();
  
  tree->Draw("errPLRUp_poisdata5s_BBfull>>haux");  
  data5sPoisBBfullPLR_err  = haux->GetMean();
  data5sPoisBBfullPLR_derr = haux->GetMeanError();
  data5sPoisBBfullPLR_med  = get_quantile(haux); 
  data5sPoisBBfullPLR_cov  = prob_poisdata5s_BBfullPLR;
  cout << "Data5s PLR BBFull:      " << prob_poisdata5s_BBfullPLR << " +/- " << TMath::Sqrt(prob_poisdata5s_BBfullPLR*(1-prob_poisdata5s_BBfullPLR)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();

  treeFC->Draw("errFCUp_poisdata_BBfull>>haux");  
  dataPoisBBfullFC_err  = haux->GetMean();
  dataPoisBBfullFC_derr = haux->GetMeanError();
  dataPoisBBfullFC_med  = get_quantile(haux); 
  dataPoisBBfullFC_cov  = prob_poisdata_BBfullFC;
  cout << "Data FC BBFull" << (doFCcheat ? "cheat:    " : " :        ")  << prob_poisdata_BBfullFC << " +/- " << TMath::Sqrt(prob_poisdata_BBfullFC*(1-prob_poisdata_BBfullFC)/ntoysFC) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();

  treeFC->Draw("errFCUp_poisdata5s_BBfull>>haux");  
  data5sPoisBBfullFC_err  = haux->GetMean();
  data5sPoisBBfullFC_derr = haux->GetMeanError();
  data5sPoisBBfullFC_med  = get_quantile(haux); 
  data5sPoisBBfullFC_cov  = prob_poisdata5s_BBfullFC;
  cout << "Data5s FC BBFull" << (doFCcheat ? "cheat:  " : " :      ") << prob_poisdata5s_BBfullFC << " +/- " << TMath::Sqrt(prob_poisdata5s_BBfullFC*(1-prob_poisdata5s_BBfullFC)/ntoysFC) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();

  tree->Draw("err_adhocdata>>haux");  
  adhocdata_err  = haux->GetMean();
  adhocdata_derr = haux->GetMeanError();
  adhocdata_med  = get_quantile(haux); 
  adhocdata_cov  = prob_adhocdata;
  cout << "Ad hoc                  " << prob_adhocdata << " +/- " << TMath::Sqrt(prob_adhocdata*(1-prob_adhocdata)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();

  tree->Draw("err_adhocdata5s>>haux");  
  adhocdata5s_err  = haux->GetMean();
  adhocdata5s_derr = haux->GetMeanError();
  adhocdata5s_med  = get_quantile(haux); 
  adhocdata5s_cov  = prob_adhocdata5s;
  cout << "Ad hoc                  " << prob_adhocdata5s << " +/- " << TMath::Sqrt(prob_adhocdata5s*(1-prob_adhocdata5s)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();

  tree->Draw("err_propdata>>haux");  
  propdata_err  = haux->GetMean();
  propdata_derr = haux->GetMeanError();
  propdata_med  = get_quantile(haux); 
  propdata_cov  = prob_propdata;
  cout << "Prop                    " << prob_propdata << " +/- " << TMath::Sqrt(prob_propdata*(1-prob_propdata)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")";
  haux->Reset();

  tree->Draw("TMath::Sqrt(err_propdata*err_propdata - err_data*err_data)>>haux");  
  double jvar_propdata_err  = haux->GetMean();
  double jvar_propdata_derr = haux->GetMeanError();
  cout << ":  jac var =  " << haux->GetMean() << " +/- " << haux->GetMeanError() << " => r = " << haux->GetMean()/(data_err/TMath::Sqrt(lumiscale)) << endl;
  haux->Reset();
  
  tree->Draw("err_propdata5s>>haux");  
  propdata5s_err  = haux->GetMean();
  propdata5s_derr = haux->GetMeanError();
  propdata5s_med  = get_quantile(haux); 
  propdata5s_cov  = prob_propdata5s;
  cout << "Prop                    " << prob_propdata5s << " +/- " << TMath::Sqrt(prob_propdata5s*(1-prob_propdata5s)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ", median: " << get_quantile(haux) <<  ")" << endl;
  haux->Reset();

  treesum->Fill();
  
  delete haux;
  
  fout->cd();
  h_true_0->Write();
  h_true_1->Write();
  tree->Write();
  treeFC->Write();
  treesum->Write();
  
  sw.Stop();

  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;

  
  fout->Close(); 
  for(auto r : rans) delete r;

  return 1;

}
