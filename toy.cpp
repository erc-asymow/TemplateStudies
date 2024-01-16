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
typedef ROOT::VecOps::RVec<double> RVecD;
using ROOT::RDF::RNode; 

using namespace boost::program_options;

constexpr double MW = 0.500;
constexpr double GW = 1.0; 
constexpr int NMASS = 9;
constexpr double DELTAM = 0.2;

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
	("lumi",        value<long>()->default_value(1000), "number of events")
	("scalejac",    value<double>()->default_value(10.), "scale")
	("degs_x",      value<int>()->default_value(2), "max degree in x of corrxy")
	("tag",         value<std::string>()->default_value("closure"), "run type")
	("run",         value<std::string>()->default_value("closure"), "run type")
	("do_fit",      bool_switch()->default_value(false), "do_fit")
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
      if (vm.count("degs_x"))     std::cout << "Degree in x of corrxy: " << vm["degs_x"].as<int>() << '\n';
    }
  catch (const error &ex)
    {
      std::cerr << ex.what() << '\n';
    }

  long nevents    = vm["nevents"].as<long>();
  long lumi       = vm["lumi"].as<long>();
  double scalejac = vm["scalejac"].as<double>();
  std::string tag = vm["tag"].as<std::string>();
  std::string run = vm["run"].as<std::string>();
  int degs_x      = vm["degs_x"].as<int>();
  int seed        = vm["seed"].as<int>();
  bool do_fit     = vm["do_fit"].as<bool>();

  if(vm.count("degs_x")) tag += std::string(Form("_%d", degs_x));
      
  int x_nbins   = 100; 
  double x_low  = 0.0;
  double x_high = 1.0;

  unsigned int njacs = degs_x + 1;
  
  auto toy_mass = [&](double x, double M, double G){
    return 1./TMath::Pi()/(1 + (x-M)*(x-M)/(G*G/4))*2./G; // non-relativistic
  };

  TFile* fout = TFile::Open(("root/toy_"+tag+"_"+run+".root").c_str(), "RECREATE");

  ROOT::RDataFrame d(nevents*2);

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
  
  dlast = std::make_unique<RNode>(dlast->Define("weights_mass", 
						[&](double x)->RVecD{
						  RVecD out;
 						  double gen = toy_mass(x, MW, GW);
						  out.emplace_back( toy_mass(x, MW, GW)/gen  );
						  for(unsigned int i=0; i<NMASS; i++)
						    out.emplace_back( toy_mass(x, MW - DELTAM*0.5 + DELTAM/(NMASS-1)*i,GW)/gen );
						  return out; 
						}, {"x"}));

  dlast = std::make_unique<RNode>(dlast->Define("weights_jac", 
						[&](double x , ULong64_t rdfentry)->RVecD{
						  RVecD out;						  
						  for(unsigned int k = 0; k<=degs_x; k++){
						    double corrx = cheb(x, 0.5, 1.0, degs_x, k);
						    //if(rdfentry%njacs!=k) corrx *= 0.;
						    // divide events for jac_k based on pseudo-random number
						    if( int(x*1e+05)%njacs!=k) corrx *= 0.;
						    // select odd-numbered events for jacs 
						    if( int(x*1e+07)%2==0) corrx *= 0.;
						    out.emplace_back( corrx );						    
						  }
						  return out;
						}, {"x", "rdfentry_"} ));
  
  std::vector<ROOT::RDF::RResultPtr<TH1D> > histos1D;

  dlast = std::make_unique<RNode>(dlast->Define("weight_jacM",[&](double x)->double {
    double out = 4*TMath::Pi()*(1./TMath::Pi()/(1 + (x-MW)*(x-MW)/(GW*GW/4))*2./GW)*2*(x-MW)/GW;
    return out;
  }, {"x"} ));
        
  dlast = std::make_unique<RNode>(dlast->Define("weight", [](RVecD weights_mass, double x){
    double out = weights_mass.at(0);
    // select even-numbered events for templates
    if( int(x*1e+07)%2==1) out *= 0.;
    return out;
  }, {"weights_mass", "x"} ));

  for(unsigned int i=0; i<NMASS; i++){
    dlast = std::make_unique<RNode>(dlast->Define(Form("weight_mass%d", i), [i](RVecD weights_mass, double x){
      double out = weights_mass.at(1+i); 
      if( int(x*1e+07)%2==1) out *= 0.;
      return out;
    }, {"weights_mass", "x"} ));
  }

  for(unsigned int i = 0; i < njacs; i++){
    dlast = std::make_unique<RNode>(dlast->Define(Form("weight_jac%d",i), [i](RVecD weights_jac, RVecD weights_mass){ return weights_jac.at(i)*weights_mass.at(0);},
						  {"weights_jac", "weights_mass"} ));
  }
    
  histos1D.emplace_back(dlast->Histo1D({"h", "nominal", x_nbins, x_low, x_high}, "x", "weight"));      
  for(unsigned int i=0; i<NMASS; i++){
    histos1D.emplace_back(dlast->Histo1D({Form("h_mass%d", i),"", x_nbins, x_low, x_high}, "x", Form("weight_mass%d",i)));
  }    
  for(unsigned int i = 0; i < njacs; i++){
    string hname = std::string(Form("jac_%d", i));	
    histos1D.emplace_back(dlast->Histo1D({ Form("h_jac%d",i), hname.c_str(), x_nbins, x_low, x_high}, "x", Form("weight_jac%d",i)));
  }    
  std::string hname = "jac_mass: d(pdf) / dM";
  histos1D.emplace_back(dlast->Histo1D({Form("h_jac%d", njacs), hname.c_str(), x_nbins, x_low, x_high}, "x", "weight_jacM"));

  auto colNames = dlast->GetColumnNames();
  std::cout << colNames.size() << " columns created" << std::endl;

  double total = *(dlast->Count());  

  fout->cd();
  std::cout << "Writing histos..." << std::endl;
  double sf = double(lumi)/double(nevents);
  for(auto h : histos1D){
    h->Scale(sf);
    string h_name = std::string(h->GetName());
    if(h_name.find("jac")!=string::npos){
      cout << "Scaling " << h_name << " by " << njacs << endl;
      h->Scale(double(njacs));
    }
    h->Write();
  }

  TH1D* h_nom = (TH1D*)fout->Get("h");
  TH1D* h_asy = (TH1D*)h_nom->Clone("h_asy");
  h_asy->Reset();
  for(unsigned int ib = 0; ib<x_nbins; ib++){
    h_asy->SetBinContent(ib+1, double(lumi)/x_nbins);
    h_asy->SetBinError(ib+1, 0.0);
  }
  h_asy->Write();
  for(unsigned int ij = 0; ij<njacs; ij++){
    TH1D* h_j = (TH1D*)fout->Get(Form("h_jac%d", ij));
    TH1D* h_j_asy = (TH1D*)h_j->Clone(Form("h_jac%d_asy", ij));
    h_j_asy->Reset();
    for(unsigned int ib = 0; ib<x_nbins; ib++){
      double x_i = h_j_asy->GetBinCenter(ib+1);
      double ratio = cheb(x_i, 0.5, 1.0, degs_x, ij);
      h_j_asy->SetBinContent(ib+1, h_asy->GetBinContent(ib+1)*ratio  );
      h_j_asy->SetBinError(ib+1, 0.0 );
    }
    h_j_asy->Write();
  }
  for(unsigned int im=0; im<NMASS; im++){
    TH1D* h_m = (TH1D*)fout->Get(Form("h_mass%d", im));
    TH1D* h_m_asy = (TH1D*)h_m->Clone(Form("h_mass%d_asy", im));
    h_m_asy->Reset();
    for(unsigned int ib = 0; ib<x_nbins; ib++){
      double x_i = h_m_asy->GetBinCenter(ib+1);
      double ratio = toy_mass(x_i, MW - DELTAM*0.5 + DELTAM/(NMASS-1)*im,GW)/toy_mass(x_i, MW, GW);
      h_m_asy->SetBinContent(ib+1, h_asy->GetBinContent(ib+1)*ratio  );
      h_m_asy->SetBinError(ib+1, 0.0 );
    }
    h_m_asy->Write();
  }

  
  sw.Stop();

  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;
  std::cout << "Total slots: " << dlast->GetNSlots() << std::endl;

  if(do_fit){
    unsigned int nbins = x_nbins;
    TH1D* h_nom = (TH1D*)fout->Get("h");
    TH1D* h_nom_asy = (TH1D*)fout->Get("h_asy");
    TH1D* h_sum = (TH1D*)h_nom->Clone("h_sum");
    h_sum->Reset();
    MatrixXd inv_sqrtV = MatrixXd::Zero(nbins, nbins);
    MatrixXd inv_sqrtV_asy = MatrixXd::Zero(nbins, nbins);
    MatrixXd jac(nbins, njacs);
    MatrixXd jac_asy(nbins, njacs);
    for(unsigned int ij = 0; ij<njacs; ij++){
      TH1D* h_j = (TH1D*)fout->Get(Form("h_jac%d", ij));
      TH1D* h_j_asy = (TH1D*)fout->Get(Form("h_jac%d_asy", ij));
      for(unsigned int ib = 0; ib<nbins; ib++){
	double val = h_sum->GetBinContent(ib+1);
	double err = h_sum->GetBinError(ib+1);
	double err2 = err*err;
	double valj = h_j->GetBinContent(ib+1);
	//valj *= njacs;
	double errj = h_j->GetBinError(ib+1);
	//errj *= njacs;
	double errj2 = errj*errj;
	h_sum->SetBinContent(ib+1, val + valj);
	h_sum->SetBinError(ib+1, TMath::Sqrt(err2 + errj2));
	jac(ib, ij) = valj;
	jac_asy(ib, ij) = h_j_asy->GetBinContent(ib+1);
      }
    }
    for(unsigned int ib = 0; ib<nbins; ib++){
      double var_poisson = h_sum->GetBinContent(ib+1);
      double var_poisson_asy = h_asy->GetBinContent(ib+1);
      double var_mc = h_sum->GetBinError(ib+1)*h_sum->GetBinError(ib+1)*scalejac;
      double var_tot = var_poisson + var_mc;
      inv_sqrtV(ib,ib) = 1./TMath::Sqrt( var_tot );
      inv_sqrtV_asy(ib,ib) = 1./TMath::Sqrt( var_poisson_asy );
    }

    double xx_mass[NMASS], yy_chi2[NMASS], exx_mass[NMASS], eyy_chi2[NMASS];
    double yy_jacasy_chi2[NMASS];
    double yy_basy_chi2[NMASS];
    double yy_asy_chi2[NMASS];	
    for(unsigned int im = 0; im<NMASS; im++){
      TH1D* h_m = (TH1D*)fout->Get(Form("h_mass%d", im));
      TH1D* h_m_asy = (TH1D*)fout->Get(Form("h_mass%d_asy", im));
      VectorXd y(nbins);
      VectorXd y_asy(nbins);
      for(unsigned int ib=0;ib<nbins; ib++){
	y(ib)     = h_m->GetBinContent(ib+1) - h_nom->GetBinContent(ib+1);
	y_asy(ib) = h_m_asy->GetBinContent(ib+1) - h_nom_asy->GetBinContent(ib+1);
      }
      MatrixXd A = inv_sqrtV*jac;
      VectorXd b = inv_sqrtV*y;
      MatrixXd A_jacasy = inv_sqrtV*jac_asy;
      MatrixXd A_asy = inv_sqrtV_asy*jac_asy;
      VectorXd b_asy  = inv_sqrtV_asy*y_asy;
      VectorXd b_basy = inv_sqrtV*y_asy;
      VectorXd x        = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

      if(im==0 && true){
	MatrixXd U = MatrixXd::Identity(y.size(),y.size()) - A*(A.transpose()*A).inverse()*A.transpose();
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(U);
	if (eigensolver.info() != Eigen::Success){
	  cout << "Could not eigendecompose U" << endl;
	}
	else{
 	  auto eigenvals = eigensolver.eigenvalues();
	  auto eigenvecs = eigensolver.eigenvectors();
	  if(true){
	    cout << "Number of parameters: " << x.size() << ", points: " << y.size() << ", d.o.f.: " << y.size()-x.size() << endl;
	    cout << "Sum of eigenvalues:" << eigenvals.sum() << " (num of 0/1 eigenvals = " << eigenvals.size() << ")" << endl;
	    for(unsigned int ei = 0; ei<y.size(); ei++)
	      std::cout << "Eigen(" << ei << ") = " << eigenvals(ei) << ": bT*uj = " << eigenvecs.col(ei).transpose()*b << endl;//"\n" << eigenvecs.col(ei) << std::endl;
	    
	  }
	}
      }



      VectorXd x_asy    = A_asy.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_asy);
      VectorXd x_jacasy = A_jacasy.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
      VectorXd x_basy   = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_basy);
      MatrixXd chi2old     = b.transpose()*b;
      MatrixXd chi2        = ((b - A*x).transpose())*(b-A*x);
      MatrixXd chi2_asy    = ((b_asy - A_asy*x_asy).transpose())*(b_asy-A_asy*x_asy);
      MatrixXd chi2_jacasy = ((b - A_jacasy*x_jacasy).transpose())*(b-A_jacasy*x_jacasy);
      MatrixXd chi2_basy   = ((b_basy - A*x_basy).transpose())*(b_basy-A*x_basy);
      int ndof = nbins-njacs;
      double chi2min = chi2(0,0);
      double chi2min_asy = chi2_asy(0,0);
      double chi2min_jacasy = chi2_jacasy(0,0);
      double chi2min_basy = chi2_basy(0,0);
      double chi2norm = chi2(0,0)/ndof;
      xx_mass[im] = MW - DELTAM*0.5 + DELTAM/(NMASS-1)*im; 
      exx_mass[im] = 0.0;
      yy_chi2[im] = chi2min;
      yy_jacasy_chi2[im] = chi2min_jacasy;
      yy_asy_chi2[im] = chi2min_asy;
      yy_basy_chi2[im] = chi2min_basy;
      eyy_chi2[im] = 0.0;
      cout << "POI[" << im << "]" << ": chi2_min = " << chi2min << " (jac-asy: " << chi2min_jacasy << ", b-asy: " << chi2min_basy << ", full-asy: " << chi2min_asy << ") at " << xx_mass[im] << endl;
    }

    TGraphErrors* chi2_fit = new TGraphErrors(NMASS,xx_mass,yy_chi2,exx_mass,eyy_chi2);
    chi2_fit->Fit("pol2", "Q");
    chi2_fit->Write("chi2_vs_POI");      
    TF1* parabola = chi2_fit->GetFunction("pol2");
    float param0 = parabola->GetParameter(0); 
    float param1 = parabola->GetParameter(1); 
    float param2 = parabola->GetParameter(2); 
    float deltaM = 1./TMath::Sqrt(param2);
    float biasM = -param1/param2*0.5;
    float pullM = (biasM-MW)/deltaM; 
    cout << "DeltaM        = " << deltaM << " " 
	 << " -- bias: " << (biasM-MW)  << " "
	 << ", pull: " << pullM 
	 << endl;

    TGraphErrors* chi2_fit_basy = new TGraphErrors(NMASS,xx_mass,yy_basy_chi2,exx_mass,eyy_chi2);
    chi2_fit_basy->Fit("pol2", "Q");
    chi2_fit_basy->Write("chi2_basy_vs_POI");      
    TF1* parabola_basy = chi2_fit_basy->GetFunction("pol2");
    float param0_basy = parabola_basy->GetParameter(0); 
    float param1_basy = parabola_basy->GetParameter(1); 
    float param2_basy = parabola_basy->GetParameter(2); 
    float deltaM_basy = 1./TMath::Sqrt(param2_basy);
    float biasM_basy = -param1_basy/param2_basy*0.5;
    float pullM_basy = (biasM_basy-MW)/deltaM_basy; 
    cout << "DeltaM basy   = " << deltaM_basy << " " 
	 << " -- bias: " << (biasM_basy-MW)  << " "
	 << ", pull: " << pullM_basy 
	 << endl;

    TGraphErrors* chi2_fit_jacasy = new TGraphErrors(NMASS,xx_mass,yy_jacasy_chi2,exx_mass,eyy_chi2);
    chi2_fit_jacasy->Fit("pol2", "Q");
    chi2_fit_jacasy->Write("chi2_jacasy_vs_POI");      
    TF1* parabola_jacasy = chi2_fit_jacasy->GetFunction("pol2");
    float param0_jacasy = parabola_jacasy->GetParameter(0); 
    float param1_jacasy = parabola_jacasy->GetParameter(1); 
    float param2_jacasy = parabola_jacasy->GetParameter(2); 
    float deltaM_jacasy = 1./TMath::Sqrt(param2_jacasy);
    float biasM_jacasy = -param1_jacasy/param2_jacasy*0.5;
    float pullM_jacasy = (biasM_jacasy-MW)/deltaM_jacasy; 
    cout << "DeltaM jacasy = " << deltaM_jacasy << " " 
	 << " -- bias: " << (biasM_jacasy-MW)  << " "
	 << ", pull: " << pullM_jacasy 
	 << endl;

    TGraphErrors* chi2_fit_asy = new TGraphErrors(NMASS,xx_mass,yy_asy_chi2,exx_mass,eyy_chi2);
    chi2_fit_asy->Fit("pol2", "Q");
    chi2_fit_asy->Write("chi2_asy_vs_POI");      
    TF1* parabola_asy = chi2_fit_asy->GetFunction("pol2");
    float param0_asy = parabola_asy->GetParameter(0); 
    float param1_asy = parabola_asy->GetParameter(1); 
    float param2_asy = parabola_asy->GetParameter(2); 
    float deltaM_asy = 1./TMath::Sqrt(param2_asy);
    float biasM_asy = -param1_asy/param2_asy*0.5;
    float pullM_asy = (biasM_asy-MW)/deltaM_asy; 
    cout << "DeltaM asy    = " << deltaM_asy << " " 
	 << " -- bias: " << (biasM_asy-MW)  << " "
	 << ", pull: " << pullM_asy 
	 << endl;

    h_sum->Write();

    TTree* tree = new TTree("tree", "tree");
    float n_data;
    float n_effmc;
    float n_genmc;
    int n_par;
    tree->Branch("sigma",        &deltaM,     "sigma/F");
    tree->Branch("sigma_asy",    &deltaM_asy, "sigma_asy/F");
    tree->Branch("sigma_jacasy", &deltaM_jacasy, "sigma_jacasy/F");
    tree->Branch("sigma_basy",   &deltaM_basy, "sigma_basy/F");
    tree->Branch("n_data",   &n_data,  "n_data/F");
    tree->Branch("n_effmc",  &n_effmc, "n_effmc/F");
    tree->Branch("n_genmc",  &n_genmc, "n_genmc/F");
    tree->Branch("n_par",    &n_par,   "n_par/I");
    n_effmc = float( h_sum->GetEffectiveEntries() );
    n_genmc = float(nevents);
    n_data = float(lumi);
    n_par = njacs;
    tree->Fill();
    tree->Write();
  }

  
  fout->Close(); 
  for(auto r : rans) delete r;

  return 1;
}
