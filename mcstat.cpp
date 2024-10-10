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
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnHesse.h"
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
  
  Likelihood(const int& debug, const int& nbins, const double& lumiscale, const bool& BBlite )
    : errorDef_(1.0), debug_(debug), nbins_(nbins), lumiscale_(lumiscale), bblite_(BBlite) {
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
  double lumiscale_;
  unsigned int nbins_;
  unsigned int npars_;
  unsigned int ndof_;
  bool bblite_;
  int debug_;
  double errorDef_;
};

// par[0] = mu0, par[1] = mu1, par[2,..] = BB 
double Likelihood::operator()(const vector<double>& par) const {  
  double val  = 0.0;
  for(unsigned int ibin=0; ibin<nbins_; ibin++){
    double mc_tot = (mc0_[ibin]+mc1_[ibin]);
    double mc_tot_err = TMath::Sqrt( mc0err_[ibin]*mc0err_[ibin] + mc1err_[ibin]*mc1err_[ibin] );
    double rel_err_mc = mc_tot_err/mc_tot;
    double rel_err_mc0 = mc0err_[ibin]/mc0_[ibin];
    double rel_err_mc1 = mc1err_[ibin]/mc1_[ibin];
    double exp_data_i = //bblite_ ?
      //((1+par[0])*mc0_[ibin] + (1+par[1])*mc1_[ibin])*(1+par[2+ibin]) :
      //(1.0 + par[0])*(1.0 + par[2+ibin])*mc0_[ibin] + (1.0 + par[1])*(1.0 + par[2+nbins_+ibin])*mc1_[ibin];
      ((1+par[0])*mc0_[ibin] + (1+par[1])*mc1_[ibin]) ;
    double prior  = bblite_ ?  0.5*par[2+ibin]*par[2+ibin]/(rel_err_mc*rel_err_mc) : 0.0;
    if(false && !bblite_){
      prior += 0.5*par[2+ibin]*par[2+ibin]/(rel_err_mc0*rel_err_mc0);
      prior += 0.5*par[2+nbins_+ibin]*par[2+nbins_+ibin]/(rel_err_mc1*rel_err_mc1);
    }
    double chi2   = bblite_ ?
      //0.5*(data_[ibin] - exp_data_i)*(data_[ibin] - exp_data_i)/exp_data_i :
      0.5*(data_[ibin] - exp_data_i)*(data_[ibin] - exp_data_i)/( data_[ibin] + rel_err_mc*rel_err_mc*exp_data_i*exp_data_i ) :
      //0.5*(data_[ibin] - exp_data_i)*(data_[ibin] - exp_data_i)/( data_[ibin] + mc0err_[ibin]*mc0err_[ibin]*(1.0 + par[0]*par[0]) + mc1err_[ibin]*mc1err_[ibin]*(1.0 + par[1]*par[1]) ) ;
      0.5*(data_[ibin] - exp_data_i)*(data_[ibin] - exp_data_i)/( data_[ibin] + mc0err_[ibin]*mc0err_[ibin]*(1.0 + par[0])*(1.0 + par[0]) + mc1err_[ibin]*mc1err_[ibin]*(1.0 + par[1])*(1.0 +par[1]) ) ;
      //0.5*(data_[ibin] - exp_data_i)*(data_[ibin] - exp_data_i)/( mc0err_[ibin]*mc0err_[ibin] + mc1err_[ibin]*mc1err_[ibin] ) ;
    val  += chi2;
    //val  += prior;
  }
  
  /*
  if(bblite_){
    //cout << par[0] << ", " << par[1] << endl;
    for(unsigned int ibin=0; ibin<nbins_; ibin++){
      double mc_tot = (mc0_[ibin]+mc1_[ibin]);
      double mc_tot_err = TMath::Sqrt( mc0err_[ibin]*mc0err_[ibin] + mc1err_[ibin]*mc1err_[ibin] );
      double rel_err_mc = mc_tot_err/mc_tot;
      double exp_data_i = ((1+par[0])*mc0_[ibin] + (1+par[1])*mc1_[ibin])*(1+par[2+ibin]);
      double prior = 0.5*par[2+ibin]*par[2+ibin]/(rel_err_mc*rel_err_mc);
      double chi2  = 0.5*(data_[ibin] - exp_data_i)*(data_[ibin] - exp_data_i)/exp_data_i;
      val  += chi2;
      val  += prior;      
      //double exp_data_i = (1+par[0])*mc0_[ibin] + (1+par[1])*mc1_[ibin];
      //double pois  = exp_data_i  - data_[ibin]*TMath::Log(exp_data_i);
      //double pois0 = mc_tot - data_[ibin]*TMath::Log(mc_tot) ;
      //val0  += pois0;
      //val  += pois;
      //double chi2 = 0.5*(data_[ibin] - exp_data_i)*(data_[ibin] - exp_data_i)/data_[ibin];
      //cout << "-lnL += " << -1.0*pois << ",  chi2 += " << chi2 << endl; 
      //val  += 0.5*(data_[ibin] - exp_data_i)*(data_[ibin] - exp_data_i)/exp_data_i;
      //cout << "Bin" << ibin << ": " << -2.0*(pois - pois0) <<  " + " << prior << endl;
      //rel_err_mc  <<  ", prior: " << prior << " --> " << val << endl;
    }
  }
  else{
    for(unsigned int ibin=0; ibin<nbins_; ibin++){
      double mc_tot = (mc0_[ibin]+mc1_[ibin]);
      double rel_err_mc0 = mc0err_[ibin]/mc0_[ibin];
      double rel_err_mc1 = mc1err_[ibin]/mc1_[ibin];
      double exp_data_i = (1.0+par[0])*(1.0+par[2+ibin])*mc0_[ibin] + (1.0+par[1])*(1.0+par[2+nbins_+ibin])*mc1_[ibin];
      //double pois   = -exp_data_i + data_[ibin]*TMath::Log(exp_data_i);
      //double pois0  = -mc_tot     + data_[ibin]*TMath::Log(mc_tot) ;
      double chi2  = 0.5*(data_[ibin] - exp_data_i)*(data_[ibin] - exp_data_i)/exp_data_i;
      double prior0 = par[2+ibin]*par[2+ibin]/(rel_err_mc0*rel_err_mc0);
      double prior1 = par[2+nbins_+ibin]*par[2+nbins_+ibin]/(rel_err_mc1*rel_err_mc1);
      //cout << -2.0*(pois - pois0) << " + " << prior0 << " + " << prior1 << endl;
      //val += -2.0*(pois - pois0) ;
      val += chi2;
      val += prior0;
      val += prior1;
    }
  }
  //cout << val << endl;
  //val0 *= -2.0;
  */
  
  return val;
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
	("tag",         value<std::string>()->default_value("closure"), "run type")
	("do_fit",      bool_switch()->default_value(false), "do_fit")
	("verbose",      bool_switch()->default_value(false), "verbose")
	("nbins",       value<int>()->default_value(200), "nbins")
	("nsigmas",     value<int>()->default_value(5), "nsigmas")
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
  long ntoys    = vm["ntoys"].as<long>();
  std::string tag = vm["tag"].as<std::string>();
  int seed        = vm["seed"].as<int>();
  int nbins       = vm["nbins"].as<int>();
  int nsigmas     = vm["nsigmas"].as<int>();
  float frac      = vm["frac"].as<float>();
  float asym      = vm["asym"].as<float>();
  float lumiscale    = vm["lumiscale"].as<float>();
  bool do_fit     = vm["do_fit"].as<bool>();
  bool verbose     = vm["verbose"].as<bool>();

  std::vector<TRandom3*> rans = {};
  for(unsigned int i = 0; i < 10; i++){
    rans.emplace_back( new TRandom3(seed + i*10) );
  }

  TFile* fout = TFile::Open(("root/mcstat_"+tag+".root").c_str(), "RECREATE");
  TTree* tree = new TTree("tree", "");
  double mu_data, err_data;
  double mu_poisdata_BB, err_poisdata_BB;
  double mu_poisdata_BBfull, err_poisdata_BBfull;
  double mu_poisdata5s_BB, err_poisdata5s_BB;
  double mu_poisdata5s_BBfull, err_poisdata5s_BBfull;
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
  tree->Branch("mu_poisdata5s_BB",  &mu_poisdata5s_BB, "mu_poisdata5s_BB/D");
  tree->Branch("err_poisdata5s_BB", &err_poisdata5s_BB, "err_poisdata5s_BB/D");
  tree->Branch("mu_poisdata5s_BBfull",  &mu_poisdata5s_BBfull, "mu_poisdata5s_BBfull/D");
  tree->Branch("err_poisdata5s_BBfull", &err_poisdata5s_BBfull, "err_poisdata5s_BBfull/D");
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

  for(unsigned int ir=0; ir<nbins; ir++){
    A_true(ir, 0) = h_true_0->GetBinContent(ir+1);
    A_true(ir, 1) = h_true_1->GetBinContent(ir+1);
  }
  for(unsigned int ir=0; ir<nbins; ir++){
    invV_true(ir,ir) = 1./h_true_tot->GetBinContent(ir+1);
  }

  MatrixXd C_true = ( A_true.transpose()*invV_true*A_true ).inverse();
  //cout << "True errors on mu_[0,1] = [" << TMath::Sqrt(C_true(0,0)) << "," << TMath::Sqrt(C_true(1,1)) << "]" << endl;

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
  MatrixXd invV_itoy      = MatrixXd::Zero(nbins, nbins);
  MatrixXd invV5s_itoy    = MatrixXd::Zero(nbins, nbins);
  MatrixXd invV_BB_itoy   = MatrixXd::Zero(nbins, nbins);
  MatrixXd invV5s_BB_itoy = MatrixXd::Zero(nbins, nbins);
  MatrixXd invVMC_itoy    = MatrixXd::Zero(nbins, nbins);
  VectorXd y_itoy(nbins);
  VectorXd y5s_itoy(nbins);
  VectorXd yMC_itoy(nbins);
  VectorXd ynom_itoy(nbins);
  
  double prob_data = 0.;
  double prob_data_BB = 0.;
  double prob_poisdata_BB = 0.;
  double prob_poisdata_BBfull = 0.;
  double prob_data5s = 0.;
  double prob_data5s_BB = 0.;
  double prob_poisdata5s_BB = 0.;
  double prob_poisdata5s_BBfull = 0.;
  double prob_mc = 0.;

  unsigned int maxfcn(numeric_limits<unsigned int>::max());
  double tolerance(0.001);
  int verbosity = int(nevents<2);
  //int npars_lite = 2 + nbins;
  int npars_lite = 2 ;
  //int npars_full = 2 + 2*nbins;
  int npars_full = 2 ; //+ 2*nbins;
  //ROOT::Minuit2::MnPrint::SetGlobalLevel(verbosity);
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
    
  for(unsigned int itoy=0; itoy<ntoys; itoy++){
    if(itoy%10000==0) cout << "Doing toy " << itoy << " / " << ntoys << endl;
      
    MatrixXd invC_itoy      = MatrixXd::Zero(2, 2);
    MatrixXd invC5s_itoy    = MatrixXd::Zero(2, 2);
    MatrixXd invCMC_itoy    = MatrixXd::Zero(2, 2);
    MatrixXd invC_BB_itoy   = MatrixXd::Zero(2, 2);
    MatrixXd invC5s_BB_itoy = MatrixXd::Zero(2, 2);
    MatrixXd invCMC_BB_itoy = MatrixXd::Zero(2, 2);

    Likelihood* likelihood_lite = new Likelihood(0, nbins, lumiscale, true);
    likelihood_lite->SetErrorDef(0.5);
    Likelihood* likelihood_full = new Likelihood(0, nbins, lumiscale, false);
    likelihood_full->SetErrorDef(0.5);
    Likelihood* likelihood_lite5s = new Likelihood(0, nbins, lumiscale, true);
    likelihood_lite5s->SetErrorDef(0.5);
    Likelihood* likelihood_full5s = new Likelihood(0, nbins, lumiscale, false);
    likelihood_full5s->SetErrorDef(0.5);
    
    for(unsigned int ir=0; ir<nbins; ir++){      
      A_itoy(ir,0) = rans[0]->Poisson( A_true(ir, 0)*lumiscale ) / lumiscale;
      A_itoy(ir,1) = rans[0]->Poisson( A_true(ir, 1)*lumiscale ) / lumiscale;
      
      // "MC" template drawn from true nominal, including lumi scale
      ynom_itoy(ir) = A_itoy.row(ir).sum();
      
      // toy data drawn from true nominal 
      y_itoy(ir)   = rans[0]->Poisson( A_true.row(ir).sum() );
      
      // toy data drawn from true 5sigma 
      y5s_itoy(ir) = rans[0]->Poisson( A_true(ir,0)*(1+err_true*nsigmas) + A_true(ir,1)  );

      // toy data drawn from "MC" 
      yMC_itoy(ir) = rans[0]->Poisson( A_itoy.row(ir).sum() );

      invV_itoy(ir,ir)      = 1./y_itoy(ir);
      invV5s_itoy(ir,ir)    = 1./y5s_itoy(ir);
      invV_BB_itoy(ir,ir)   = 1./(y_itoy(ir)   + ynom_itoy(ir)/lumiscale);
      invV5s_BB_itoy(ir,ir) = 1./(y5s_itoy(ir) + ynom_itoy(ir)/lumiscale);
      invVMC_itoy(ir,ir)    = 1./yMC_itoy(ir);

      // this is a common piece
      MatrixXd K_ir_itoy = A_itoy.row(ir).transpose()*A_itoy.row(ir);      
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
    VectorXd x_itoy      = C_itoy*A_itoy.transpose()*invV_itoy*(y_itoy - ynom_itoy);
    VectorXd x5s_itoy    = C5s_itoy*A_itoy.transpose()*invV5s_itoy*(y5s_itoy - ynom_itoy);
    VectorXd xMC_itoy    = CMC_itoy*A_itoy.transpose()*invVMC_itoy*(yMC_itoy - ynom_itoy);
    VectorXd x_BB_itoy   = C_BB_itoy*A_itoy.transpose()*invV_BB_itoy*(y_itoy - ynom_itoy);
    VectorXd x5s_BB_itoy = C5s_BB_itoy*A_itoy.transpose()*invV5s_BB_itoy*(y5s_itoy - ynom_itoy);
    mu_data      = x_itoy(0);
    mu_data5s    = x5s_itoy(0);
    mu_data_BB   = x_BB_itoy(0);
    mu_data5s_BB = x5s_BB_itoy(0);
    mu_mc = xMC_itoy(0); 
    err_data      = TMath::Sqrt(C_itoy(0,0));
    err_data5s    = TMath::Sqrt(C5s_itoy(0,0));
    err_data_BB   = TMath::Sqrt(C_BB_itoy(0,0));
    err_data5s_BB = TMath::Sqrt(C5s_BB_itoy(0,0));
    err_mc        = TMath::Sqrt(CMC_itoy(0,0));

    MnUserParameters upar_lite;
    //upar_lite.Add("mu0", x_BB_itoy(0), TMath::Sqrt(C_BB_itoy(0,0)), x_BB_itoy(0) - 5.0*TMath::Sqrt(C_BB_itoy(0,0)), x_BB_itoy(0) + 5.0*TMath::Sqrt(C_BB_itoy(0,0)) );
    //upar_lite.Add("mu1", x_BB_itoy(1), TMath::Sqrt(C_BB_itoy(1,1)), x_BB_itoy(1) - 5.0*TMath::Sqrt(C_BB_itoy(1,1)), x_BB_itoy(1) + 5.0*TMath::Sqrt(C_BB_itoy(1,1)) );
    upar_lite.Add("mu0", 0., 0.01);
    upar_lite.Add("mu1", 0., 0.01);
    for (int i=0; i<nbins; i++){
      double exp_error = TMath::Sqrt( 1./(lumiscale*(A_itoy(i,0)+A_itoy(i,1))) );
      //upar_lite.Add(Form("BBLite%d",i), 0.0, exp_error*0.001, -10*exp_error, 10*exp_error );
      //upar_lite.Add(Form("BBLite%d",i), 0.0, 0.001);
    }
    
    MnMigrad migrad_lite(*likelihood_lite, upar_lite, 1);
    MnMigrad migrad_lite5s(*likelihood_lite5s, upar_lite, 1);

    if(verbose) cout << "\tMigrad..." << endl;
    FunctionMinimum min_lite = migrad_lite(maxfcn, tolerance);
    double edm_lite = double(min_lite.Edm());
    double fmin_lite = double(min_lite.Fval());
    if(verbose) cout << "\tHesse POST..." << endl;
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
    if(verbose) cout << "\tHesse POST..." << endl;
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
    for (int i=0; i<nbins; i++){
      double exp_error = TMath::Sqrt( 1.0/(lumiscale*A_itoy(i,0) ) );
      //upar_full.Add(Form("BB0Full%d",i), 0.0, exp_error, -5*exp_error, 5*exp_error );
      //upar_full.Add(Form("BB0Full%d",i), 0.0, 0.001 );
    }
    for (int i=0; i<nbins; i++){
      double exp_error = TMath::Sqrt( 1.0/(lumiscale*A_itoy(i,1) ) );
      //upar_full.Add(Form("BB1Full%d",i), 0.0, exp_error, -5*exp_error, 5*exp_error );
      //upar_full.Add(Form("BB1Full%d",i), 0.0, 0.001);
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
    
    if( TMath::Abs(mu_data-mu_true)/err_data <= 1.0 ) prob_data += 1./ntoys;
    if( TMath::Abs(mu_data5s-(mu_true + err_true*nsigmas))/err_data5s <= 1.0 ) prob_data5s += 1./ntoys;
    if( TMath::Abs(mu_data_BB-mu_true)/err_data_BB <= 1.0 ) prob_data_BB += 1./ntoys;
    if( TMath::Abs(mu_data5s_BB-(mu_true + err_true*nsigmas))/err_data5s_BB <= 1.0 ) prob_data5s_BB += 1./ntoys;
    if( TMath::Abs(mu_mc-mu_true)/err_mc <= 1.0 ) prob_mc += 1./ntoys;
    if( TMath::Abs(mu_poisdata_BB-mu_true)/err_poisdata_BB <= 1.0 ) prob_poisdata_BB += 1./ntoys;
    if( TMath::Abs(mu_poisdata_BBfull-mu_true)/err_poisdata_BBfull <= 1.0 ) prob_poisdata_BBfull += 1./ntoys;
    if( TMath::Abs(mu_poisdata5s_BB-(mu_true + err_true*nsigmas))/err_poisdata5s_BB <= 1.0 ) prob_poisdata5s_BB += 1./ntoys;
    if( TMath::Abs(mu_poisdata5s_BBfull-(mu_true + err_true*nsigmas))/err_poisdata5s_BBfull <= 1.0 ) prob_poisdata5s_BBfull += 1./ntoys;

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
  
  TH1D* haux = new TH1D("haux","", 500, 0., 0.5);

  tree->Draw("err_mc>>haux");  
  cout << "Asympt. err:        " << TMath::Sqrt(C_true(0,0)) << ", coverage = " << 1.0 - TMath::Prob(1.0, 1)  << endl;
  cout << "MC:                 " << prob_mc << " +/- " << TMath::Sqrt(prob_mc*(1-prob_mc)/ntoys) <<  " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ")" << endl;
  haux->Reset();

  tree->Draw("err_data>>haux");  
  cout << "Data:               " << prob_data << " +/- " << TMath::Sqrt(prob_data*(1-prob_data)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ")" << endl;
  haux->Reset();

  tree->Draw("err_data5s>>haux");  
  cout << "Data 5s:            " << prob_data5s << " +/- " << TMath::Sqrt(prob_data5s*(1-prob_data5s)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ")" << endl;
  haux->Reset();

  tree->Draw("err_data_BB>>haux");  
  cout << "Data BB:            " << prob_data_BB << " +/- " << TMath::Sqrt(prob_data_BB*(1-prob_data_BB)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ")" << endl;
  haux->Reset();

  tree->Draw("err_data5s_BB>>haux");  
  cout << "Data 5s BB:         " << prob_data5s_BB << " +/- " << TMath::Sqrt(prob_data5s_BB*(1-prob_data5s_BB)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ")" << endl;
  haux->Reset();

  tree->Draw("err_poisdata_BB>>haux");  
  cout << "Data Pois BB:       " << prob_poisdata_BB << " +/- " << TMath::Sqrt(prob_poisdata_BB*(1-prob_poisdata_BB)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ")" << endl;
  haux->Reset();

  tree->Draw("err_poisdata5s_BB>>haux");  
  cout << "Data5s Pois BB:     " << prob_poisdata5s_BB << " +/- " << TMath::Sqrt(prob_poisdata5s_BB*(1-prob_poisdata5s_BB)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ")" << endl;
  haux->Reset();

  tree->Draw("err_poisdata_BBfull>>haux");  
  cout << "Data Pois BBFull:   " << prob_poisdata_BBfull << " +/- " << TMath::Sqrt(prob_poisdata_BBfull*(1-prob_poisdata_BBfull)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ")" << endl;
  haux->Reset();
  
  tree->Draw("err_poisdata5s_BBfull>>haux");  
  cout << "Data5s Pois BBFull: " << prob_poisdata5s_BBfull << " +/- " << TMath::Sqrt(prob_poisdata5s_BBfull*(1-prob_poisdata5s_BBfull)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ")" << endl;
  haux->Reset();

  delete haux;
  
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
