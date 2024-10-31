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
    npars_ = 2;
    ndof_ = nbins - 2;
    data_.reserve( nbins );
    mc0_.reserve( nbins );
    mc1_.reserve( nbins );
    mc0err_.reserve( nbins );
    mc1err_.reserve( nbins );    
    for(unsigned int ibin = 0; ibin<nbins; ibin++) data_.push_back(0.);
    for(unsigned int ibin = 0; ibin<nbins; ibin++) mc0_.push_back(0.);
    for(unsigned int ibin = 0; ibin<nbins; ibin++) mc1_.push_back(0.);
    for(unsigned int ibin = 0; ibin<nbins; ibin++) mc0err_.push_back(0.);
    for(unsigned int ibin = 0; ibin<nbins; ibin++) mc1err_.push_back(0.);
  }

  unsigned int get_n_dof(){ return ndof_;}

  void set_data(const unsigned int& ibin, const double& val){
    data_[ibin] = val; 
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
    double mc_tot = (mc0_[ibin]+mc1_[ibin]);
    double mc_tot_err = TMath::Sqrt( mc0err_[ibin]*mc0err_[ibin] + mc1err_[ibin]*mc1err_[ibin] );
    double rel_err_mc = mc_tot_err/mc_tot;
    double rel_err_mc0 = mc0err_[ibin]/mc0_[ibin];
    double rel_err_mc1 = mc1err_[ibin]/mc1_[ibin];
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
      double chi2 = 0.5*res2/data_i;
      if(doPoisson_){
	chi2 *= (1.0 - 2./3.*(res/data_i) + 0.5*res2/data_i/data_i /* + ... */ );
      }
      val  += chi2;
      val  += prior;
    }
    else{
      double exp_data_i = (1+par0)*mc0_[ibin] + (1+par1)*mc1_[ibin] ;
      double res  = data_i - exp_data_i;
      double res2 = res*res;
      double chi2   = bblite_ ?
	0.5*res2/( data_i + TMath::Power(rel_err_mc*exp_data_i, 2.0 ) ) :
	0.5*res2/( data_i + TMath::Power(mc0err_[ibin]*(1.0 + par0), 2.0) + TMath::Power(mc1err_[ibin]*(1.0 + par1), 2.0) ) ;
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
	("nbins",       value<int>()->default_value(200), "nbins")
	("nsigmas",     value<int>()->default_value(5), "nsigmas")
	("frac",       value<float>()->default_value(0.5), "frac")
	("asym",       value<float>()->default_value(0.015), "asym")
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
  float lumiscale    = vm["lumiscale"].as<float>();
  bool doFC     = vm["doFC"].as<bool>();
  bool doFCcheat = vm["doFCcheat"].as<bool>();
  bool doBarlett = vm["doBarlett"].as<bool>();
  bool saveHistos  = vm["saveHistos"].as<bool>();
  bool FCfixToTrue = vm["FCfixToTrue"].as<bool>();
  bool verbose     = vm["verbose"].as<bool>();
  bool profileMCNP = vm["profileMCNP"].as<bool>();
  bool decorrelate = vm["decorrelate"].as<bool>();
  bool doPoisson = vm["doPoisson"].as<bool>();
  
  std::vector<TRandom3*> rans = {};
  for(unsigned int i = 0; i < 10; i++){
    rans.emplace_back( new TRandom3(seed + i*10) );
  }

  TFile* fout = TFile::Open(("root/mcstat_"+tag+".root").c_str(), "RECREATE");
  TTree* tree = new TTree("tree", "");
  TTree* treeFC = new TTree("treeFC", "");
  double mu_data, err_data;
  double mu_poisdata_BB, err_poisdata_BB;
  double mu_poisdata_BBfull, err_poisdata_BBfull;
  double errPLRLow_poisdata_BBfull,errPLRUp_poisdata_BBfull;
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
  tree->Branch("mu_data",  &mu_data, "mu_data/D");
  tree->Branch("err_data", &err_data, "err_data/D");
  tree->Branch("mu_poisdata_BB",  &mu_poisdata_BB, "mu_poisdata_BB/D");
  tree->Branch("err_poisdata_BB", &err_poisdata_BB, "err_poisdata_BB/D");
  tree->Branch("mu_poisdata_BBfull",  &mu_poisdata_BBfull, "mu_poisdata_BBfull/D");
  tree->Branch("err_poisdata_BBfull", &err_poisdata_BBfull, "err_poisdata_BBfull/D");
  tree->Branch("errPLRLow_poisdata_BBfull", &errPLRLow_poisdata_BBfull, "errPLRLow_poisdata_BBfull/D");
  tree->Branch("errPLRUp_poisdata_BBfull", &errPLRUp_poisdata_BBfull, "errPLRUp_poisdata_BBfull/D");
  treeFC->Branch("errFCLow_poisdata_BBfull", &errFCLow_poisdata_BBfull, "errFCLow_poisdata_BBfull/D");
  treeFC->Branch("errFCUp_poisdata_BBfull", &errFCUp_poisdata_BBfull, "errFCUp_poisdata_BBfull/D");
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
  double proc_ratio = h_true_0->Integral()/h_true_1->Integral();
  
  cout << "h_true_0: integral = " << h_true_0->Integral() << endl;
  cout << "h_true_1: integral = " << h_true_1->Integral() << endl;
  cout << "h_true_tot: integral = " << h_true_tot->Integral() << endl;

  MatrixXd A_true = MatrixXd::Zero(nbins, 2);
  MatrixXd invV_true = MatrixXd::Zero(nbins, nbins);
  MatrixXd J_true = MatrixXd::Zero(nbins, 2);

  for(unsigned int ir=0; ir<nbins; ir++){
    A_true(ir, 0) = h_true_0->GetBinContent(ir+1);
    A_true(ir, 1) = h_true_1->GetBinContent(ir+1);
    J_true(ir, 0) = A_true(ir, 0);
    J_true(ir, 1) = A_true(ir, 1);
    if(decorrelate){
      J_true(ir, 0) = A_true(ir, 0) - A_true(ir, 1);
      J_true(ir, 1) = A_true(ir, 0) + A_true(ir, 1);
    }
  }
  for(unsigned int ir=0; ir<nbins; ir++){
    invV_true(ir,ir) = 1./h_true_tot->GetBinContent(ir+1);
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
  
  /*
  MatrixXd C_inv = MatrixXd::Zero(2, 2);
  for(unsigned int ib=0; ib<nbins; ib++){
    C_inv += (A_true.row(ib).transpose()*A_true.row(ib))/h_true_tot->GetBinContent(ib+1);
  }
  MatrixXd C_alt = C_inv.inverse();
  cout << "Errors on mu_[0,1] = [" << TMath::Sqrt(C_alt(0,0)) << "," << TMath::Sqrt(C_alt(1,1)) << "]" << endl;
  */

  MatrixXd A_itoy = MatrixXd::Zero(nbins, 2);
  MatrixXd J_itoy = MatrixXd::Zero(nbins, 2);
  MatrixXd invV_itoy      = MatrixXd::Zero(nbins, nbins);
  MatrixXd invV5s_itoy    = MatrixXd::Zero(nbins, nbins);
  MatrixXd invV_BB_itoy   = MatrixXd::Zero(nbins, nbins);
  MatrixXd invV5s_BB_itoy = MatrixXd::Zero(nbins, nbins);
  MatrixXd invVMC_itoy    = MatrixXd::Zero(nbins, nbins);
  VectorXd y_itoy(nbins);
  VectorXd yFC_itoy(nbins);
  VectorXd y5s_itoy(nbins);
  VectorXd yMC_itoy(nbins);
  VectorXd ynom_itoy(nbins);
  
  double prob_data = 0.;
  double prob_data_BB = 0.;
  double prob_poisdata_BB = 0.;
  double prob_poisdata_BBfull = 0.;
  double prob_poisdata_BBfullPLR = 0.;
  double prob_poisdata_BBfullFC = 0.;
  double prob_data5s = 0.;
  double prob_data5s_BB = 0.;
  double prob_poisdata5s_BB = 0.;
  double prob_poisdata5s_BBfull = 0.;
  double prob_poisdata5s_BBfullPLR = 0.;
  double prob_poisdata5s_BBfullFC = 0.;
  double prob_mc = 0.;

  unsigned int maxfcn(numeric_limits<unsigned int>::max());
  double tolerance(0.001);
  int verbosity = int(nevents<2);
  int npars_lite = 2 ;
  if(profileMCNP) npars_lite += nbins;
  int npars_full = 2 ;
  if(profileMCNP) npars_full +=	2*nbins;

  ROOT::Minuit2::MnPrint::SetGlobalLevel(verbosity);
  MatrixXd Vmin_lite    = MatrixXd::Zero(npars_lite,npars_lite);    
  VectorXd xmin_lite    = VectorXd::Zero(npars_lite);
  VectorXd xmin_liteErr = VectorXd::Zero(npars_lite);
  MatrixXd Vmin_full    = MatrixXd::Zero(npars_full,npars_full);    
  VectorXd xmin_full    = VectorXd::Zero(npars_full);
  VectorXd xmin_fullErr = VectorXd::Zero(npars_full);
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
	  ynom_itoy(ir) = A_itoy.row(ir).sum();
	  if(idata==0)
	    y_itoy(ir)   = rans[1+idata]->Poisson( A_true(ir,0) + A_true(ir,1) );
	  else{
	    y_itoy(ir)   = rans[1+idata]->Poisson( A_true(ir,0)*(1.0 + err_true*nsigmas) + A_true(ir,1) );
	    if(decorrelate)
	      y_itoy(ir)   = rans[1+idata]->Poisson( A_true(ir,0)*(1.0 + err_true*nsigmas) + A_true(ir,1)*(1.0 - err_true*nsigmas) );
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

	uparFC_full.SetValue("mu0", idata==0 ? 0.0 : err_true*nsigmas );
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
	uparFC_full.Release("mu0");

	//cout << "Toy: fval " << fval_full <<  " --> " << fval_fullFix << endl;

	VectorXd jtilde = VectorXd::Zero(2*nbins);
	VectorXd jhat   = VectorXd::Zero(2*nbins);
	
	MatrixXd invsqrtV  = MatrixXd::Zero(1,1);
	MatrixXd invVj = MatrixXd::Zero(2,2);
	VectorXd jhati = VectorXd::Zero(2);
	VectorXd y = VectorXd::Zero(1);
	MatrixXd X = MatrixXd::Zero(1, 2);
	for(unsigned int i=0; i<nbins; i++){
	  invsqrtV(0,0) = 1./TMath::Sqrt(y_itoy(i));
	  invVj(0,0) = 1./(A_itoy(i,0)*lumiscale);
	  invVj(1,1) = 1./(A_itoy(i,1)*lumiscale);       
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
	  //cout << "Bin " << i << endl;
	  //cout << "\tj0: " << jtildei(0) << " --> " << jhati(0) << " (" << (jtildei(0) - jhati(0))/TMath::Sqrt(A_itoy(i,0)*lumiscale) << ")" << endl;
	  //cout << "\tj1: " << jtildei(1) << " --> " << jhati(1) << " (" << (jtildei(1) - jhati(1))/TMath::Sqrt(A_itoy(i,1)*lumiscale) << ")" << endl;
	}
	
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
	  // <--
	  h_obs_itoy->SetBinContent(ir+1, y_itoy(ir));
	}

	TString hname = "";
	if(idata==0)
	  hname = Form("h_data_teststat_%d", itoy );
	else
	  hname = Form("h_data5s_teststat_%d", itoy );
	TH1D* h_teststat_itoy = new TH1D( hname, "", 200, 0., 20);

	//int count = 0;
	for(unsigned int itoyFC=0; itoyFC<1000; itoyFC++){
	  for(unsigned int ir=0; ir<nbins; ir++){      	  
	    if(doFCcheat){
	      yFC_itoy(ir) = rans[1+idata]->Poisson(  ( 1.0 + (idata==1)*err_true*nsigmas)*A_true(ir,0) + A_true(ir,1) );
	      if(decorrelate)
		yFC_itoy(ir) = rans[1+idata]->Poisson(  ( 1.0 + (idata==1)*err_true*nsigmas  )*A_true(ir,0) + (1.0 - (idata==1)*err_true*nsigmas )*A_true(ir,1) );
	    }
	    else{
	      yFC_itoy(ir) = rans[1+idata]->Poisson( (1.0 + mu0Fix)*jtilde(2*ir) + (1.0 + mu1Fix)*jtilde(2*ir+1) );
	      if(decorrelate)
		yFC_itoy(ir) = rans[1+idata]->Poisson( (1.0 + mu0Fix + mu1Fix)*jtilde(2*ir) + (1.0 - mu0Fix + mu1Fix)*jtilde(2*ir+1) );
	    }
	    likelihoodFC_full->set_data(ir, yFC_itoy(ir) );
	  }
	  uparFC_full.Release("mu0");
	  uparFC_full.SetValue("mu0", 0.0);
	  uparFC_full.SetValue("mu1", 0.0);
	  MnMigrad migradFC_full_itoy(*likelihoodFC_full, uparFC_full, 1);
	  FunctionMinimum minFC_full_itoy = migradFC_full_itoy(maxfcn, tolerance);
	  double fval_full_itoy = double(minFC_full_itoy.Fval());
	  uparFC_full.Fix("mu0");
	  uparFC_full.SetValue("mu0", mu0Fix );
	  uparFC_full.SetValue("mu1", 0.0);
	  MnMigrad migradFC_full_itoyFix(*likelihoodFC_full, uparFC_full, 1);
	  FunctionMinimum minFC_full_itoyFix = migradFC_full_itoyFix(maxfcn, tolerance);
	  double fval_full_itoyFix = double(minFC_full_itoyFix.Fval());
	  //cout << "\t FC Tot fval " << fval_full_itoy <<  " --> " << fval_full_itoyFix << endl;
	  //cout << 2*(fval_full_itoyFix - fval_full_itoy) << endl;
	  h_teststat_itoy->Fill( 2*(fval_full_itoyFix - fval_full_itoy) );
	  //if( 2*(fval_full_itoyFix - fval_full_itoy) < 1.0 ) count++;
	  uparFC_full.Release("mu0");
	}

	// restore initial data
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
	  if(  0.0 >= mu0Hat+errFCLow_poisdata_BBfull && 0.0 <= mu0Hat+errFCUp_poisdata_BBfull  )
	    prob_poisdata_BBfullFC += 1./ntoysFC;
	}
	else{
	  if( err_true*nsigmas >= mu0Hat+errFCLow_poisdata5s_BBfull && err_true*nsigmas <= mu0Hat+errFCUp_poisdata5s_BBfull)
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

    Likelihood* likelihood_lite = new Likelihood(0, nbins, true, doPoisson, profileMCNP, decorrelate);
    likelihood_lite->SetErrorDef(0.5);
    Likelihood* likelihood_full = new Likelihood(0, nbins, false, doPoisson, profileMCNP, decorrelate);
    likelihood_full->SetErrorDef(0.5);
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
	J_itoy(ir,0) = A_itoy(ir,0) - A_itoy(ir,1);
	J_itoy(ir,1) = A_itoy(ir,0) + A_itoy(ir,1);
      }
      
      // "MC" template drawn from true nominal, including lumi scale
      ynom_itoy(ir) = A_itoy.row(ir).sum();
      
      // toy data drawn from true nominal 
      y_itoy(ir)   = rans[0]->Poisson( A_true.row(ir).sum() );
      
      // toy data drawn from true 5sigma 
      y5s_itoy(ir) = rans[0]->Poisson( A_true(ir,0)*(1.0 + err_true*nsigmas) + A_true(ir,1)  );
      if(decorrelate)
	y5s_itoy(ir) = rans[0]->Poisson( A_true(ir,0)*(1.0 + err_true*nsigmas) + A_true(ir,1)*(1.0 - err_true*nsigmas)  );
      
      // toy data drawn from "MC" 
      yMC_itoy(ir) = rans[0]->Poisson( A_itoy.row(ir).sum() );

      invV_itoy(ir,ir)      = 1./y_itoy(ir);
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
    
    MnMigrad migrad_lite(*likelihood_lite, upar_lite, 1);
    MnMigrad migrad_lite5s(*likelihood_lite5s, upar_lite, 1);

    if(verbose) cout << "\tMigrad..." << endl;
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
    
    if(verbose) cout << "\tMigrad..." << endl;
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

    mu_poisdata_BB      = xmin_lite(0);
    err_poisdata_BB     = xmin_liteErr(0);
    mu_poisdata_BBfull  = xmin_full(0);
    err_poisdata_BBfull = xmin_fullErr(0);
    mu_poisdata5s_BB      = xmin_lite5s(0);
    err_poisdata5s_BB     = xmin_lite5sErr(0);
    mu_poisdata5s_BBfull  = xmin_full5s(0);
    err_poisdata5s_BBfull = xmin_full5sErr(0);
    
    if( TMath::Abs(mu_data-mu_true)/err_data <= 1.0 )
      prob_data += 1./ntoys;
    if( TMath::Abs(mu_data5s-(mu_true + err_true*nsigmas))/err_data5s <= 1.0 )
      prob_data5s += 1./ntoys;
    if( TMath::Abs(mu_data_BB-mu_true)/err_data_BB <= 1.0 )
      prob_data_BB += 1./ntoys;
    if( TMath::Abs(mu_data5s_BB-(mu_true + err_true*nsigmas))/err_data5s_BB <= 1.0 )
      prob_data5s_BB += 1./ntoys;
    if( TMath::Abs(mu_mc-mu_true)/err_mc <= 1.0 )
      prob_mc += 1./ntoys;
    if( TMath::Abs(mu_poisdata_BB-mu_true)/err_poisdata_BB <= 1.0 )
      prob_poisdata_BB += 1./ntoys;
    if( TMath::Abs(mu_poisdata_BBfull-mu_true)/err_poisdata_BBfull <= 1.0 )
      prob_poisdata_BBfull += 1./ntoys;
    if( TMath::Abs(mu_poisdata5s_BB-(mu_true + err_true*nsigmas))/err_poisdata5s_BB <= 1.0 )
      prob_poisdata5s_BB += 1./ntoys;
    if( TMath::Abs(mu_poisdata5s_BBfull-(mu_true + err_true*nsigmas))/err_poisdata5s_BBfull <= 1.0 )
      prob_poisdata5s_BBfull += 1./ntoys;
    if( mu_poisdata_BBfull >= (mu_true + errPLRLow_poisdata_BBfull) && mu_poisdata_BBfull <= (mu_true + errPLRUp_poisdata_BBfull)  )
      prob_poisdata_BBfullPLR += 1./ntoys;
    if( mu_poisdata5s_BBfull >= (mu_true + err_true*nsigmas + errPLRLow_poisdata5s_BBfull) && mu_poisdata5s_BBfull <= (mu_true + err_true*nsigmas + errPLRUp_poisdata5s_BBfull) )
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
  double data_err, data_derr, data_med, data_cov;
  double data5s_err, data5s_derr, data5s_med, data5s_cov;
  double dataBB_err, dataBB_derr, dataBB_med, dataBB_cov;
  double data5sBB_err, data5sBB_derr, data5sBB_med, data5sBB_cov;
  double dataPoisBB_err, dataPoisBB_derr, dataPoisBB_med, dataPoisBB_cov;
  double data5sPoisBB_err, data5sPoisBB_derr, data5sPoisBB_med, data5sPoisBB_cov;
  double dataPoisBBfull_err, dataPoisBBfull_derr, dataPoisBBfull_med, dataPoisBBfull_cov;
  double data5sPoisBBfull_err, data5sPoisBBfull_derr, data5sPoisBBfull_med, data5sPoisBBfull_cov;
  double dataPoisBBfullPLR_err, dataPoisBBfullPLR_derr, dataPoisBBfullPLR_med, dataPoisBBfullPLR_cov;
  double data5sPoisBBfullPLR_err, data5sPoisBBfullPLR_derr, data5sPoisBBfullPLR_med, data5sPoisBBfullPLR_cov;
  double dataPoisBBfullFC_err, dataPoisBBfullFC_derr, dataPoisBBfullFC_med, dataPoisBBfullFC_cov;
  double data5sPoisBBfullFC_err, data5sPoisBBfullFC_derr, data5sPoisBBfullFC_med, data5sPoisBBfullFC_cov;

  treesum->Branch("n",  &n,  "n/I");
  treesum->Branch("nFC",&nFC,  "nFC/I");
  treesum->Branch("asym_corr",  &asym_corr,  "asym_corr/D");
  treesum->Branch("asym_cond",  &asym_cond,  "asym_cond/D");
  treesum->Branch("asym_err",  &asym_err,  "asym_err/D");
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
  treesum->Branch("dataPoisBB_err",  &dataPoisBB_err,  "dataPoisBB_err/D");
  treesum->Branch("dataPoisBB_derr", &dataPoisBB_derr, "dataPoisBB_derr/D");
  treesum->Branch("dataPoisBB_med",  &dataPoisBB_med,  "dataPoisBB_med/D");
  treesum->Branch("dataPoisBB_cov",  &dataPoisBB_cov,  "dataPoisBB_cov/D");
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

  n = ntoys;
  nFC = ntoysFC;
  asym_corr = rho_true;
  asym_cond = condition_true;
  asym_err = TMath::Sqrt(C_true(0,0));
  asym_derr = 0.0;
  asym_med = asym_err;
  asym_cov = 1.0 - TMath::Prob(1.0, 1);
  
  TH1D* haux = new TH1D("haux","", 500, 0., 0.5);

  tree->Draw("err_mc>>haux");
  cout << "Asympt. err:            " << TMath::Sqrt(C_true(0,0)) << ", coverage = " << 1.0 - TMath::Prob(1.0, 1)  << endl;
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
