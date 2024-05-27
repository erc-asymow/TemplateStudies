#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include "TRandom3.h"
#include "TVector.h"
#include "TVectorT.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include <TMatrixD.h>
#include <TMatrixDSymfwd.h>
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
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

//#include <Eigen/Core>
//#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;
using namespace ROOT;
using namespace ROOT::Minuit2;

typedef ROOT::VecOps::RVec<double> RVecD;
using ROOT::RDF::RNode; 

using namespace boost::program_options;

constexpr double MZ = 91.;
constexpr double GW = 2.5;


class TheoryFcn : public FCNGradientBase {
  //class TheoryFcn : public FCNBase {

public:
  TheoryFcn(const int& debug, const int& seed) : errorDef_(1.0), debug_(debug), seed_(seed)
  {

    ran_ = new TRandom3(seed);
    pt_edges_  = {25, 30, 35, 40, 45, 50, 55}; 
    eta_edges_ = {-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};

    // pt_edges_  = {25, 30, 35}; 
    //eta_edges_ = {-3.0, -2.5, -2.0};

    n_pt_bins_  = pt_edges_.size()-1;
    n_eta_bins_ = eta_edges_.size()-1;

    n_pars_ = 3*n_eta_bins_;
    
    for(unsigned int i = 0; i < pt_edges_.size(); i++){
      k_edges_.emplace_back(1./pt_edges_[i]);
    }
    for(unsigned int i = 0; i < n_pt_bins_; i++){
      kmean_vals_.emplace_back( 0.5*(k_edges_[i]+k_edges_[i+1]) );
    }

    n_data_ = n_eta_bins_*n_eta_bins_*n_pt_bins_*n_pt_bins_;
    scales2_.reserve(n_data_);
    scales2Err_.reserve(n_data_);
    for(unsigned int idata = 0; idata<n_data_; idata++){
      scales2_.push_back( 0.0 );
      scales2Err_.push_back( 0.0 );
    }
        
    // generate initial set of data points;
    //generate_data();

    n_data_ = scales2_.size();
    n_dof_ = n_data_ - n_pars_;
  }
  
  ~TheoryFcn() { delete ran_;}

  void generate_data();

  void set_seed(const int& seed){ ran_->SetSeed(seed);}
  
  unsigned int get_n_params(){ return n_pars_;}
  unsigned int get_n_data(){ return n_data_;}
  unsigned int get_n_dof(){ return n_dof_;} 
  
  virtual double Up() const {return errorDef_;}
  virtual void SetErrorDef(double def) {errorDef_ = def;}

  virtual double operator()(const vector<double>&) const;
  virtual vector<double> Gradient(const vector<double>& ) const;
  virtual bool CheckGradient() const {return true;} 

  static constexpr double scaling_A = 0.001;
  static constexpr double scaling_e = 0.001 * 40.0;
  static constexpr double scaling_M = 0.001 / 40.0;
  static constexpr double scaling_e_prime = 0.001 / 0.01; 

private:

  vector<double> scales2_;
  vector<double> scales2Err_;
  vector<double> pt_edges_;
  vector<double> k_edges_;
  vector<double> kmean_vals_;
  vector<double> eta_edges_;
  unsigned int n_pt_bins_;
  unsigned int n_eta_bins_;
  unsigned int n_data_;
  unsigned int n_pars_;
  unsigned int n_dof_;
  int debug_;
  int seed_;
  double errorDef_;
  TRandom3* ran_;
};

void TheoryFcn::generate_data(){
  //ran_->SetSeed(seed_);
  double chi2_start = 0.;
  unsigned int ibin = 0;
  for(unsigned int ieta_p = 0; ieta_p<n_eta_bins_; ieta_p++){
    for(unsigned int ipt_p = 0; ipt_p<n_pt_bins_; ipt_p++){
      for(unsigned int ieta_m = 0; ieta_m<n_eta_bins_; ieta_m++){
	for(unsigned int ipt_m = 0; ipt_m<n_pt_bins_; ipt_m++){
	  double ierr2_nom = 0.001*(1+double(ieta_p)/n_eta_bins_)*(1+double(ieta_m)/n_eta_bins_);
	  //*(2-0.1*double(ipt_p)/n_pt_bins_)*(2-0.1*double(ipt_m)/n_pt_bins_);
	  double ierr2 = ran_->Gaus(ierr2_nom,  ierr2_nom*0.1);
	  while(ierr2<=0.){
	    ierr2 = ran_->Gaus(ierr2_nom,  ierr2_nom*0.1);
	  }
	  double iscale2 = ran_->Gaus(1.00, ierr2);
	  //if(ibin<3) cout << iscale2 << endl;
	  scales2_[ibin]    = iscale2 ;
	  scales2Err_[ibin] =  ierr2 ;
	  double dchi2 = (scales2_[ibin]-1.0)/scales2Err_[ibin];
	  //cout << dchi2*dchi2 << endl;
	  chi2_start += dchi2*dchi2 ;
	  ibin++;
	}
      }
    }
  }
  cout << "Inistial chi2/ndata = " << chi2_start/(n_data_-1) << " has prob " << TMath::Prob(chi2_start, n_data_-1 ) <<  endl;
  return;
}

double TheoryFcn::operator()(const vector<double>& par) const {

  double val = 0.0;
  const unsigned int npars = par.size();

  unsigned int ibin = 0;
  for(unsigned int ieta_p = 0; ieta_p < n_eta_bins_; ieta_p++){
    double A_p = par[ieta_p];
    double e_p = par[ieta_p+n_eta_bins_];
    double M_p = par[ieta_p+2*n_eta_bins_]; 
    for(unsigned int ipt_p = 0; ipt_p < n_pt_bins_; ipt_p++){   
      double k_p = kmean_vals_[ipt_p];
      double p_term = (1.0 + A_p + e_p*k_p - M_p/k_p );
      for(unsigned int ieta_m = 0; ieta_m < n_eta_bins_; ieta_m++){
	double A_m = par[ieta_m];
	double e_m = par[ieta_m+n_eta_bins_];
	double M_m = par[ieta_m+2*n_eta_bins_];
	for(unsigned int ipt_m = 0; ipt_m < n_pt_bins_; ipt_m++){	  
	  double k_m = kmean_vals_[ipt_m];
	  double m_term = (1.0 + A_m + e_m*k_m + M_m/k_m);
	  double ival = (scales2_[ibin] - p_term*m_term)/scales2Err_[ibin];
	  double ival2 = ival*ival;
	  val += ival2;
	  ibin++;
	}
      }
    }
  }
  val /= n_dof_;

  val -= 1.0;
  
  //val = 0.0;
  //for(unsigned int ipar=0; ipar<par.size(); ipar++) val += (par[ipar]-0.001)*(par[ipar]-0.001);
  //cout << val << endl;
  
  return val;
}

vector<double> TheoryFcn::Gradient(const vector<double> &par ) const {

  //cout << "Using gradient" << endl; 
  vector<double> grad(par.size(), 0.0);

  for(unsigned int ipar = 0; ipar < par.size(); ipar++){
    unsigned int ieta     = ipar % n_eta_bins_;
    unsigned int par_type = ipar / n_eta_bins_;    
    //cout << "ipar " << ipar << ": " << ieta << ", " << par_type << endl;
    double grad_i = 0.0;    
    unsigned int ibin = 0;
    for(unsigned int ieta_p = 0; ieta_p < n_eta_bins_; ieta_p++){
      double A_p = par[ieta_p];
      double e_p = par[ieta_p+n_eta_bins_];
      double M_p = par[ieta_p+2*n_eta_bins_]; 
      for(unsigned int ipt_p = 0; ipt_p < n_pt_bins_; ipt_p++){   
	double k_p = kmean_vals_[ipt_p];
	double p_term = 0.;
	if(ieta_p != ieta) p_term = (1.0 + A_p + e_p*k_p - M_p/k_p);
	else{
	  if(par_type==0)       p_term = 1.0;
	  else if( par_type==1) p_term = k_p;
	  else                  p_term = -1./k_p;
	}
	for(unsigned int ieta_m = 0; ieta_m < n_eta_bins_; ieta_m++){
	  double A_m = par[ieta_m];
	  double e_m = par[ieta_m+n_eta_bins_];
	  double M_m = par[ieta_m+2*n_eta_bins_];
	  for(unsigned int ipt_m = 0; ipt_m < n_pt_bins_; ipt_m++){	  
	    double k_m = kmean_vals_[ipt_m];
	    double m_term = 0.;
	    if(ieta_m != ieta) m_term = (1.0 + A_m + e_m*k_m + M_m/k_m);
	    else{
	      if(par_type==0)       m_term = 1.0;
	      else if( par_type==1) m_term = k_m;
	      else                  m_term = +1./k_m;
	    }
	    double ival = -2*(scales2_[ibin] - (1.0 + A_p + e_p*k_p - M_p/k_p)*(1.0 + A_m + e_m*k_m + M_m/k_m))/scales2Err_[ibin]/scales2Err_[ibin];
	    double term = 0.0;
	    if(ieta_p==ieta || ieta_m==ieta){
	      if(ieta_p!=ieta_m) term = p_term*m_term;
	      else{
		if(par_type==0)      term = 1.0*(1.0 + A_m + e_m*k_m + M_m/k_m) + (1.0 + A_p + e_p*k_p - M_p/k_p)*1.0;
		else if(par_type==1) term = k_p*(1.0 + A_m + e_m*k_m + M_m/k_m) + (1.0 + A_p + e_p*k_p - M_p/k_p)*k_m;
		else                 term = -1/k_p*(1.0 + A_m + e_m*k_m + M_m/k_m) + (1.0 + A_p + e_p*k_p - M_p/k_p)/k_m;
	      }
	    }
	    //cout << "ival " << ival << "," << term << endl; 
	    double ig = ival*term;
	    ig /= n_dof_;
	    //cout << "ibin " << ibin << " += " << ig << endl;
	    grad_i += ig;
	    ibin++;
	  }
	}
      }
    }
    //cout << "\t" << ipar << ": " << grad_i << endl;
    //grad_i = 2*(par[ipar]-0.001);
    grad[ipar] = grad_i;
  }

  return grad; 
}


  
int main(int argc, char* argv[])
{

  TStopwatch sw;
  sw.Start();

  //ROOT::EnableImplicitMT();

  variables_map vm;
  try
    {
      options_description desc{"Options"};
      desc.add_options()
	("help,h", "Help screen")
	("nevents",     value<long>()->default_value(1000), "number of events")
	("lumi",        value<long>()->default_value(1000), "number of events")
	("tag",         value<std::string>()->default_value("closure"), "run type")
	("run",         value<std::string>()->default_value("closure"), "run type")
	("seed",        value<int>()->default_value(4357), "seed");

      store(parse_command_line(argc, argv, desc), vm);
      notify(vm);
      if (vm.count("help")){
	std::cout << desc << '\n';
	return 0;
      }
      if (vm.count("nevents"))    std::cout << "Number of events: " << vm["nevents"].as<long>() << '\n';
      if (vm.count("tag"))        std::cout << "Tag: " << vm["tag"].as<std::string>() << '\n';
      if (vm.count("run"))        std::cout << "Run: " << vm["run"].as<std::string>() << '\n';
    }
  catch (const error &ex)
    {
      std::cerr << ex.what() << '\n';
    }

  long nevents    = vm["nevents"].as<long>();
  long lumi       = vm["lumi"].as<long>();
  std::string tag = vm["tag"].as<std::string>();
  std::string run = vm["run"].as<std::string>();
  int seed        = vm["seed"].as<int>();

  TFile* fout = TFile::Open(("./massfit_"+tag+"_"+run+".root").c_str(), "RECREATE");
  TTree* tree = new TTree("tree", "tree");

  double edm, fmin, prob;
  int isvalid, hasAccurateCovar, hasPosDefCovar;
  tree->Branch("edm", &edm, "edm/D");
  tree->Branch("fmin", &fmin, "fmin/D");
  tree->Branch("prob", &prob, "prob/D");
  tree->Branch("isvalid", &isvalid, "isvalid/I");
  tree->Branch("hasAccurateCovar", &hasAccurateCovar, "hasAccurateCovar/I");
  tree->Branch("hasPosDefCovar", &hasPosDefCovar, "hasPosDefCovar/I");
  //fFCN->set_seed(seed);

  int debug = 0;
  TheoryFcn* fFCN_0 = new TheoryFcn(debug, seed);   
  unsigned int n_parameters = fFCN_0->get_n_params();
  
  vector<double> tpar(n_parameters);
  vector<double> tparErr(n_parameters);

   for (int i=0; i<n_parameters/3; i++){
    //tpar.emplace_back(0.0);
    //tparErr.emplace_back(0.0);
    tree->Branch(Form("A%d",i), &tpar[i], Form("A%d/D",i));
    tree->Branch(Form("A%derr",i), &tparErr[i], Form("A%derr/D",i));
  }
  for (int i=0; i<n_parameters/3; i++){
    //tpar.emplace_back(0.0);
    //tparErr.emplace_back(0.0);
    tree->Branch(Form("e%d",i), &tpar[i+n_parameters/3], Form("e%d/D",i));
    tree->Branch(Form("e%derr",i), &tparErr[i+n_parameters/3], Form("e%derr/D",i));
  }
  for (int i=0; i<n_parameters/3; i++){
    //tpar.emplace_back(0.0);
    //tparErr.emplace_back(0.0);
    tree->Branch(Form("M%d",i), &tpar[i+2*n_parameters/3], Form("M%d/D",i));
    tree->Branch(Form("M%derr",i), &tparErr[i+2*n_parameters/3], Form("M%derr/D",i));
  }
      

  unsigned int maxfcn(numeric_limits<unsigned int>::max());
  double tolerance(0.001);
  int verbosity = 1; 
  ROOT::Minuit2::MnPrint::SetGlobalLevel(verbosity);
  
  for(unsigned int itoy=0; itoy<nevents; itoy++){

    int debug = 0;
    TheoryFcn* fFCN = new TheoryFcn(debug, seed+itoy);  
    fFCN->generate_data();
    fFCN->SetErrorDef(1.0 / fFCN->get_n_dof());
    
    MnUserParameters upar;
    double start=0.0, par_error=0.01;
    unsigned int n_parameters = fFCN->get_n_params();
    for (int i=0; i<n_parameters/3; i++){
      upar.Add(Form("A%d",i), start, par_error);
    }
    for (int i=0; i<n_parameters/3; i++){
      upar.Add(Form("e%d",i), start, par_error);
    }
    for (int i=0; i<n_parameters/3; i++){
      upar.Add(Form("M%d",i), start, par_error);      
    }

    MnMigrad migrad(*fFCN, upar, 1);    

    //fFCN->set_seed(seed);

    cout << "Migrad" << endl;
    FunctionMinimum min = migrad(maxfcn, tolerance);

    edm = double(min.Edm());
    fmin = double(min.Fval());
    prob = TMath::Prob((min.Fval()+1)*fFCN->get_n_dof(), fFCN->get_n_dof() );
    isvalid = int(min.IsValid());
    hasAccurateCovar = int(min.HasAccurateCovar());
    hasPosDefCovar = int(min.HasPosDefCovar());

    cout << "Hesse" << endl;
    MnHesse hesse(1);
    hesse(*fFCN, min);
    
    for(unsigned int ipar = 0 ; ipar<n_parameters; ipar++){
      //cout << upar.GetName(ipar) << ":" << min.UserState().Value(ipar) << " +/- " << min.UserState().Error(ipar) << endl;
      tpar[ipar]    = min.UserState().Value(ipar) ;
      tparErr[ipar] = min.UserState().Error(ipar) ;
    }
    tree->Fill();
    
    TH2D* hcov = new TH2D(Form("hcov_%d", itoy), "", n_parameters, 0, n_parameters, n_parameters, 0, n_parameters);
    TH2D* hcor = new TH2D(Form("hcor_%d", itoy), "", n_parameters, 0, n_parameters, n_parameters, 0, n_parameters);  
    for(unsigned int iparX = 0 ; iparX<n_parameters; iparX++){    
      hcov->GetXaxis()->SetBinLabel(iparX+1, TString(upar.GetName(iparX).c_str()) );
      hcor->GetXaxis()->SetBinLabel(iparX+1, TString(upar.GetName(iparX).c_str()) );
      for(unsigned int iparY = 0 ; iparY<n_parameters; iparY++){
	//cout << min.UserState().Covariance().Data()[iparX*n_parameters+iparY] << endl;
	double cov_XY = iparX>iparY ?
	  min.UserState().Covariance().Data()[iparY+ iparX*(iparX+1)/2] :
	  min.UserState().Covariance().Data()[iparX+ iparY*(iparY+1)/2];
	double cor_XY = cov_XY / (TMath::Sqrt(min.UserState().Covariance().Data()[iparX+ iparX*(iparX+1)/2])*
				  TMath::Sqrt(min.UserState().Covariance().Data()[iparY+ iparY*(iparY+1)/2]));
	hcov->SetBinContent(iparX+1, iparY+1, cov_XY);
	hcov->GetYaxis()->SetBinLabel(iparY+1, TString(upar.GetName(iparY).c_str()) );
	hcor->SetBinContent(iparX+1, iparY+1, cor_XY);
	hcor->GetYaxis()->SetBinLabel(iparY+1, TString(upar.GetName(iparY).c_str()) );
      }
    }
    hcor->SetMinimum(-1.0);
    hcor->SetMaximum(+1.0);
    fout->cd();
    if(itoy<20){
      hcor->Write();
      //hcov->Write();
    }
    
    if(verbosity>1){
      cout << "Data points: " << fFCN->get_n_data() << endl;
      cout << "Number of parameters: " << fFCN->get_n_params() << endl;
      cout << "chi2/ndf: " << min.Fval()+1 << " (prob: " << TMath::Prob((min.Fval()+1)*fFCN->get_n_dof(), fFCN->get_n_dof() ) << ")" <<  endl;;
      cout << "min is valid: " << min.IsValid() << std::endl;
      cout << "HesseFailed: " << min.HesseFailed() << std::endl;
      cout << "HasCovariance: " << min.HasCovariance() << std::endl;
      cout << "HasValidCovariance: " << min.HasValidCovariance() << std::endl;
      cout << "HasValidParameters: " << min.HasValidParameters() << std::endl;
      cout << "IsAboveMaxEdm: " << min.IsAboveMaxEdm() << std::endl;
      cout << "HasReachedCallLimit: " << min.HasReachedCallLimit() << std::endl;
      cout << "HasAccurateCovar: " << min.HasAccurateCovar() << std::endl;
      cout << "HasPosDefCovar : " << min.HasPosDefCovar() << std::endl;
      cout << "HasMadePosDefCovar : " << min.HasMadePosDefCovar() << std::endl;
    }
    delete fFCN;
  }

  fout->cd();
  tree->Write();
  
  /*
  int x_nbins   = 100; 
  double x_low  = 0.0;
  double x_high = 1.0;

  auto toy_mass = [&](double x, double M, double G){
    return 1./TMath::Pi()/(1 + (x-M)*(x-M)/(G*G/4))*2./G; // non-relativistic
  };


  TFile* fout = TFile::Open(("./massfit_"+tag+"_"+run+".root").c_str(), "RECREATE");

  ROOT::RDataFrame d(nevents);

  unsigned int nslots = d.GetNSlots();
  std::vector<TRandom3*> rans = {};
  for(unsigned int i = 0; i < nslots; i++){
    rans.emplace_back( new TRandom3(seed + i*10) );
  }

  auto dlast = std::make_unique<RNode>(d);

  dlast = std::make_unique<RNode>(dlast->DefineSlot("x", [&](unsigned int nslot)->double{
    double out;
    out = rans[nslot]->Uniform(0.0, 1.0);
    return out;
  } ));
          
  dlast = std::make_unique<RNode>(dlast->Define("weight", [](RVecD weights_mass, double x){
    double out = weights_mass.at(0);
    // select even-numbered events for templates
    if( int(x*1e+07)%2==1) out *= 0.;
    return out;
  }, {"weights_mass", "x"} ));


  std::vector<ROOT::RDF::RResultPtr<TH1D> > histos1D;
  histos1D.emplace_back(dlast->Histo1D({"h", "nominal", x_nbins, x_low, x_high}, "x", "weight"));      

  auto colNames = dlast->GetColumnNames();
  std::cout << colNames.size() << " columns created" << std::endl;

  double total = *(dlast->Count());  

  fout->cd();
  std::cout << "Writing histos..." << std::endl;
  double sf = double(lumi)/double(nevents);
  for(auto h : histos1D){
    h->Scale(sf);
    string h_name = std::string(h->GetName());
    h->Write();
  }
  //std::cout << "Total slots: " << dlast->GetNSlots() << std::endl;
  */
  sw.Stop();

  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;

  fout->Close(); 
  //for(auto r : rans) delete r;

  return 1;
}
