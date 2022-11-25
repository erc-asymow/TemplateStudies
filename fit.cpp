#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TH2D.h"
#include "TF1.h"
#include "TF2.h"
#include "TFitResultPtr.h"
#include <TStopwatch.h>
#include <TGraphErrors.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>
#include <boost/program_options.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

using namespace boost::program_options;

constexpr double MW = 80.;
constexpr int NMAX  = 200;
//constexpr int NMASS = 50;
constexpr int NMASS = 3;
constexpr double DELTAM = 0.200;
constexpr double MASSSHIFT = 0.050;


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
	("nevents",     value<long>()->default_value(1000), "number of events")
	("ntoys",       value<int>()->default_value(-1), "number of toys")
	("degs_pdf_x",  value<int>()->default_value(5), "max degree of pdf_x")
	("degs_pdf_y" , value<int>()->default_value(5), "max degree of pdf_y")
	("degs_corr_x", value<int>()->default_value(2), "max degree in x of corrxy")
	("degs_corr_y", value<int>()->default_value(2), "max degree in y of corrxy")
	("degs_A0_x",   value<int>()->default_value(2), "max degree in x for A0")
	("degs_A0_y",   value<int>()->default_value(2), "max degree in y for A0")
	("degs_A1_x",   value<int>()->default_value(2), "max degree in x for A1")
	("degs_A1_y",   value<int>()->default_value(2), "max degree in y for A1")
	("degs_A2_x",   value<int>()->default_value(2), "max degree in x for A2")
	("degs_A2_y",   value<int>()->default_value(2), "max degree in y for A2")
	("degs_A3_x",   value<int>()->default_value(2), "max degree in x for A3")
	("degs_A3_y",   value<int>()->default_value(2), "max degree in y for A3")
	("degs_A4_x",   value<int>()->default_value(2), "max degree in x for A4")
	("degs_A4_y",   value<int>()->default_value(2), "max degree in y for A4")
	("rebinX",   value<int>()->default_value(-1), "rebin X axis")
	("rebinY",   value<int>()->default_value(-1), "rebin Y axis")
	("jUL",   bool_switch()->default_value(false), "")
	("j0",    bool_switch()->default_value(false), "")
	("j1",    bool_switch()->default_value(false), "")
	("j2",    bool_switch()->default_value(false), "")
	("j3",    bool_switch()->default_value(false), "")
	("j4",    bool_switch()->default_value(false), "")
	("jM",    bool_switch()->default_value(false), "")
	("scale0",bool_switch()->default_value(false), "")
	("scale1",bool_switch()->default_value(false), "")
	("scale2",bool_switch()->default_value(false), "")
	("scale3",bool_switch()->default_value(false), "")
	("scale4",bool_switch()->default_value(false), "")
	("jacmass", value<int>()->default_value(-1), "")
	("add_MC_uncert", bool_switch()->default_value(false), "")
	("X_max", value<float>()->default_value(99.), "")
	("X_min", value<float>()->default_value(-99.), "")
	("Y_max", value<float>()->default_value(99.), "")
	("Y_min", value<float>()->default_value(-99.), "")
	("verbose", bool_switch()->default_value(false), "")
	("debug", bool_switch()->default_value(false), "")
	("tag", value<std::string>()->default_value(""), "tag name")
	("post_tag", value<std::string>()->default_value(""), "post tag name")
	("run", value<std::string>()->default_value("closure"), "run type");

      store(parse_command_line(argc, argv, desc), vm);
      notify(vm);
      if (vm.count("help")){
	std::cout << desc << '\n';
	return 0;
      }
    }
  catch (const error &ex)
    {
      std::cerr << ex.what() << '\n';
    }

  long nevents    = vm["nevents"].as<long>();
  int ntoys       = vm["ntoys"].as<int>();
  std::string tag = vm["tag"].as<std::string>();
  std::string post_tag = vm["post_tag"].as<std::string>();
  std::string run = vm["run"].as<std::string>();
  int degs_pdf_x  = vm["degs_pdf_x"].as<int>();
  int degs_pdf_y  = vm["degs_pdf_y"].as<int>();
  int degs_corr_x = vm["degs_corr_x"].as<int>();
  int degs_corr_y = vm["degs_corr_y"].as<int>();
  int degs_A0_x   = vm["degs_A0_x"].as<int>();
  int degs_A0_y   = vm["degs_A0_y"].as<int>();
  int degs_A1_x   = vm["degs_A1_x"].as<int>();
  int degs_A1_y   = vm["degs_A1_y"].as<int>();
  int degs_A2_x   = vm["degs_A2_x"].as<int>();
  int degs_A2_y   = vm["degs_A2_y"].as<int>();
  int degs_A3_x   = vm["degs_A3_x"].as<int>();
  int degs_A3_y   = vm["degs_A3_y"].as<int>();
  int degs_A4_x   = vm["degs_A4_x"].as<int>();
  int degs_A4_y   = vm["degs_A4_y"].as<int>();
  int rebinX      = vm["rebinX"].as<int>();
  int rebinY      = vm["rebinY"].as<int>();
  int jacmass     = vm["jacmass"].as<int>();
  float X_max     = vm["X_max"].as<float>();
  float X_min     = vm["X_min"].as<float>();
  float Y_max     = vm["Y_max"].as<float>();
  float Y_min     = vm["Y_min"].as<float>();
  int jUL  = vm["jUL"].as<bool>();
  int j0   = vm["j0"].as<bool>();
  int j1   = vm["j1"].as<bool>();
  int j2   = vm["j2"].as<bool>();
  int j3   = vm["j3"].as<bool>();
  int j4   = vm["j4"].as<bool>();
  int jM   = vm["jM"].as<bool>();
  bool add_MC_uncert = vm["add_MC_uncert"].as<bool>();
  int verbose = vm["verbose"].as<bool>();
  int debug = vm["debug"].as<bool>();

  if(vm.count("degs_corr_x")) tag += std::string(Form("_UL_%d", degs_corr_x));
  if(vm.count("degs_corr_y")) tag += std::string(Form("_%d", degs_corr_y));
  if(vm.count("degs_A0_x"))   tag += std::string(Form("_A0_%d", degs_A0_x));
  if(vm.count("degs_A0_y"))   tag += std::string(Form("_%d", degs_A0_y));
  if(vm.count("degs_A1_x"))   tag += std::string(Form("_A1_%d", degs_A1_x));
  if(vm.count("degs_A1_y"))   tag += std::string(Form("_%d", degs_A1_y));
  if(vm.count("degs_A2_x"))   tag += std::string(Form("_A2_%d", degs_A2_x));
  if(vm.count("degs_A2_y"))   tag += std::string(Form("_%d", degs_A2_y));
  if(vm.count("degs_A3_x"))   tag += std::string(Form("_A3_%d", degs_A3_x));
  if(vm.count("degs_A3_y"))   tag += std::string(Form("_%d", degs_A3_y));
  if(vm.count("degs_A4_x"))   tag += std::string(Form("_A4_%d", degs_A4_x));
  if(vm.count("degs_A4_y"))   tag += std::string(Form("_%d", degs_A4_y));

  TFile* fin = TFile::Open(("root/histos_"+tag+"_"+run+".root").c_str(), "READ");
  if(fin==0 || fin==nullptr || fin->IsZombie()){
    cout << "File NOT found" << endl;
    return 0;
  }

  auto get_nbins_XY = [X_max,X_min,Y_max,Y_min](TH2D* h)-> std::pair<int,int> {
    int nx = 0;
    int ny = 0;
    for(unsigned int ix = 1; ix<=h->GetXaxis()->GetNbins(); ix++ ){	
      if( h->GetXaxis()->GetBinLowEdge(ix)>=X_min &&
	  (h->GetXaxis()->GetBinLowEdge(ix)+h->GetXaxis()->GetBinWidth(ix))<=X_max ) nx++;
    }
    for(unsigned int iy = 1; iy<=h->GetYaxis()->GetNbins(); iy++ ){	
      if( h->GetYaxis()->GetBinLowEdge(iy)>=Y_min &&
	  (h->GetYaxis()->GetBinLowEdge(iy)+h->GetYaxis()->GetBinWidth(iy))<=Y_max ) ny++;
    }
    return std::make_pair(nx,ny);
  };

  auto accept_bin = [&](TH2D* h, int ix, int iy)->bool{
    if( h->GetXaxis()->GetBinLowEdge(ix)>=X_min && (h->GetXaxis()->GetBinLowEdge(ix)+h->GetXaxis()->GetBinWidth(ix))<=X_max &&
	h->GetYaxis()->GetBinLowEdge(iy)>=Y_min && (h->GetYaxis()->GetBinLowEdge(iy)+h->GetYaxis()->GetBinWidth(iy))<=Y_max)
      return true;
    return false;
  };

  auto rescale_aftercut = [&](TH2D* h)->double{
    double all = h->Integral();
    double cut{0.0};
    for(unsigned int ix = 1; ix<=h->GetXaxis()->GetNbins(); ix++ ){
      for(unsigned int iy = 1; iy<=h->GetYaxis()->GetNbins(); iy++ ){
	if(accept_bin(h,ix,iy)) cut += h->GetBinContent(ix,iy); 
      }
    }
    cout << "After range restriction, " << h->GetName() << " scaled by "
	 << cut << "/" << all << "=" << cut/all << endl;  
    return cut/all;
  };

  
  TRandom3* ran = new TRandom3();
  
  Long64_t N;
  double poi_val[NMAX];
  int poi_cat[NMAX];
  unsigned int poi_idx[NMAX];
  unsigned int poi_counter;
  unsigned int n_pdfx;
  unsigned int n_pdfy;
  double points_x[ NMAX ];
  double points_y[ NMAX ];

  TTree* outtree = fin->Get<TTree>("outtree");
  outtree->SetBranchAddress("nevents", &N);
  outtree->SetBranchAddress("poi_counter", &poi_counter);
  outtree->SetBranchAddress("poi_val", &poi_val);
  outtree->SetBranchAddress("poi_cat", &poi_cat);
  outtree->SetBranchAddress("poi_idx", &poi_idx);
  outtree->SetBranchAddress("n_pdfx",   &n_pdfx);
  outtree->SetBranchAddress("n_pdfy",   &n_pdfy);
  outtree->SetBranchAddress("points_x", &points_x);
  outtree->SetBranchAddress("points_y", &points_y);
  outtree->GetEntry(0);    

  if(verbose) cout << "Listing all POIs..." << endl;
  std::vector<unsigned int> active_pois = {};
  for(unsigned int p = 0; p < poi_counter; p++){
    if(verbose) cout << "(" << p << ", " << poi_cat[p] << "), ";
    if(poi_cat[p]==-1 && jUL) active_pois.emplace_back(p);
    if(poi_cat[p]==0  && j0)  active_pois.emplace_back(p);
    if(poi_cat[p]==1  && j1)  active_pois.emplace_back(p);
    if(poi_cat[p]==2  && j2)  active_pois.emplace_back(p);      
    if(poi_cat[p]==3  && j3)  active_pois.emplace_back(p);      
    if(poi_cat[p]==4  && j4)  active_pois.emplace_back(p);
    if(poi_cat[p]==5  && jM)  active_pois.emplace_back(p);      
  }  
  poi_counter = active_pois.size();
  if(verbose){
    cout << endl;
    cout << poi_counter << " selected POIs: ";
    for(auto p : active_pois) cout << p << ", ";
    cout << endl;
  }

  std::vector<unsigned int> cleaned_active_pois = {};
  std::vector<int> helicities = {-1, 0, 1, 2, 3, 4, 5};
  for(auto hel : helicities){
    TH2D* hjac_first = 0;
    for(auto p : active_pois){
      if((hel<0 || hel==5) && poi_cat[p]==hel){
	cleaned_active_pois.emplace_back(p);
	continue;
      }
      else if( poi_cat[p]==hel && vm[std::string(Form("scale%d", hel))].as<bool>()  ){
	if(hjac_first==0){
	  hjac_first = fin->Get<TH2D>(Form("jac_%d", p));
	  hjac_first->Scale( poi_val[p] );
	  if(verbose) cout << "Found first jacobian of A" << hel << ", scaled by initial value " << poi_val[p] << endl;
	  cleaned_active_pois.emplace_back(p);
	}
	else{
	  TH2D* hjac_next = fin->Get<TH2D>(Form("jac_%d", p));
	  hjac_next->Scale( poi_val[p]);
	  hjac_first->Add( hjac_next );
	  if(verbose) cout << "Found another jacobian of A" << hel << ", scaled by initial value " << poi_val[p] << ": adding up..." << endl;
	}
      }
      else if(poi_cat[p]==hel && !vm[std::string(Form("scale%d", hel))].as<bool>()){
	cleaned_active_pois.emplace_back(p);
      }
    }
  }  
    
  active_pois.clear();
  active_pois = cleaned_active_pois;
  if(verbose){
    cout << active_pois.size() << " POIs after cleaning: " << endl;
    for(auto p : active_pois) cout << p << ", ";
    cout << endl;
  }
  poi_counter = active_pois.size();

  // hw is the "approximated" template
  TString hw_name = "h";
  if(jacmass==0)      hw_name += "";
  else if(jacmass==1) hw_name += "_up";
  else if(jacmass==2) hw_name += "_down";

  cout << "Starting template is histo=" << hw_name << endl;

  TH2D* hwMC = fin->Get<TH2D>("hMC");
  TH2D* hw   = fin->Get<TH2D>( hw_name );
  vector<TH2D*> all_histos = { hw, hwMC };
  for(unsigned int i=0; i<NMASS; i++) all_histos.push_back( fin->Get<TH2D>(Form("hMC_mass%d",i)) );
  if(rebinX>0 || rebinY>0){
    for(auto h : all_histos) h->Rebin2D(rebinX,rebinY);
  }

  int nx = hw->GetXaxis()->GetNbins(); 
  int ny = hw->GetYaxis()->GetNbins(); 
  int nx_cut = get_nbins_XY(hw).first;
  int ny_cut = get_nbins_XY(hw).second;
  double extra_scale = rescale_aftercut(hwMC);
  int nbins = nx_cut*ny_cut;

  double lumi_fact = nevents>0 ? nevents/hwMC->Integral() : 1.0;
  lumi_fact /= extra_scale;
  cout << "Scaling input histos by " << lumi_fact << endl;  
  N *= lumi_fact;
  for(auto h : all_histos) h->Scale(lumi_fact);

  cout << "Number of data points: " << nx_cut << "*" << ny_cut << " = " << nbins << endl;
  cout << "Data integral: " << hwMC->Integral() << endl;
  
  unsigned int bin_counter = 0;

  MatrixXd jac(nbins, poi_counter);
  for(unsigned int j = 0; j<poi_counter; j++){
    unsigned int idx = active_pois[j];
    TString jac_name(Form("jac_%d", idx));
    if(jacmass==0)      jac_name = TString(Form("jac0_%d", idx));
    else if(jacmass==1) jac_name = TString(Form("jac1_%d", idx));
    else if(jacmass==2) jac_name = TString(Form("jac2_%d", idx));
    TH2D* hjac = fin->Get<TH2D>( jac_name );
    // fall-back for consistency
    if( jacmass<0 && (hjac==0 || hjac==nullptr) ){
      hjac = fin->Get<TH2D>( TString(Form("jac0_%d", idx)) );
    }
    if(verbose) cout << "\tUsing jacobian histo=" << hjac->GetName() << endl;
    
    if(rebinX>0 || rebinY>0){
      hjac->Rebin2D(rebinX,rebinY);
    }
    bin_counter = 0;
    for(unsigned int ix = 1; ix<=nx; ix++ ){
      for(unsigned int iy = 1; iy<=ny; iy++ ){
	if(!accept_bin(hjac,ix,iy)) continue;
	jac(bin_counter,j) = N*hjac->GetBinContent(ix,iy); 
	bin_counter++;
      }
    }
  }

  //cout << "Jacobians filled." << endl;

  MatrixXd inv_sqrtV(nbins, nbins);
  MatrixXd inv_V(nbins, nbins);
  for(unsigned int ib = 0; ib<nbins; ib++ ){
    for(unsigned int jb = 0; jb<nbins; jb++ ){
      inv_sqrtV(ib,jb) = 0.;
      inv_V(ib,jb) = 0.;
    }
  }
  bin_counter = 0;
  for(unsigned int ix = 1; ix<=nx; ix++ ){
    for(unsigned int iy = 1; iy<=ny; iy++ ){
      if(!accept_bin(hwMC,ix,iy)) continue;
      double var = hwMC->GetBinContent(ix,iy);
      if(add_MC_uncert) var += (hwMC->GetBinError(ix,iy)*hwMC->GetBinError(ix,iy));
      inv_sqrtV(bin_counter,bin_counter) = 1./TMath::Sqrt(var);
      inv_V(bin_counter,bin_counter) = 1./var;
      bin_counter++;
    }
  }
  //cout << "V matrix filled" << endl;

  VectorXd norms_chebDummy(20);
  for(unsigned int i=0; i<20; i++) norms_chebDummy(i) = 1.0; 
    
  VectorXd norms_cheb4(5);
  norms_cheb4 << 0.0332073, 0.266696, 0.400194, 0.266696, 0.0332073;
  VectorXd norms_cheb6(7);
  norms_cheb6 << 0.0141919, 0.126953, 0.228635, 0.26044, 0.228635, 0.126953, 0.0141919; 
  VectorXd norms_cheb8(9);
  norms_cheb8 << 0.00792723, 0.0731216, 0.139826, 0.180815, 0.196621, 0.180815, 0.139826, 0.0731216, 0.00792723;

  //assert(degs_corr_y==4 || degs_corr_y==6 || degs_corr_y==8);

  // plotting
  vector<string> fit_tags = {"UL", "0", "1", "2", "3", "4" };
  string fit_tag = "";
  for(auto t : fit_tags){
    if(vm["j"+t].as<bool>()) fit_tag += "_j"+t;
  }
  TFile* fout = TFile::Open(("root/fit_"+tag+"_"+run+fit_tag+"_"+post_tag+".root").c_str(), "RECREATE");
  TTree* tree = new TTree("tree", "tree");
  double best_mass;
  double best_mass_err;
  tree->Branch("best_mass",     &best_mass,     "best_mass/D");
  tree->Branch("best_mass_err", &best_mass_err, "best_mass_err/D");
  //cout << "Commons done. Output file opened." << endl;

  TTree* tree2 = new TTree("tree2", "tree2");
  double chi2_start;
  double chi2_min;
  double mass_test;
  tree2->Branch("chi2_start", &chi2_start, "chi2_start/D");
  tree2->Branch("chi2_min", &chi2_min,     "chi2_min/D");
  tree2->Branch("mass_test", &mass_test,   "mass_test/D");

  bool do_toys = ntoys>0;
  if(!do_toys) ntoys = 1;

  double xx_mass[NMASS], yy_chi2[NMASS], exx_mass[NMASS], eyy_chi2[NMASS];
  
  for(unsigned int itoy=0; itoy<ntoys; itoy++ ){

    TString toy_tag(!do_toys ? "" : Form("_toy%d", itoy));

    VectorXd mu_ran(nbins);
    // randomize
    if(do_toys){
      //cout << "Toy " << itoy << endl; 
      bin_counter = 0;
      for(unsigned int ix = 1; ix<=nx; ix++ ){
	for(unsigned int iy = 1; iy<=ny; iy++ ){
	  if(!accept_bin(all_histos[0],ix,iy)) continue;
	  double mu = all_histos[0]->GetBinContent(ix,iy);
	  mu_ran(bin_counter) = ran->PoissonD(mu);  
	  bin_counter++;
	}
      }	
    }

    // mass specific
    for(int m = (!do_toys ? -1 : 0); m<NMASS; m++){

      if(jacmass>=0 && m>=0) continue;
      
      //cout << "Now doing mass idx " << m << endl;
      //TH2D* hMC = m<0 ? all_histos[1] : fin->Get<TH2D>(Form("hMC_mass%d",m));
      TH2D* hMC = all_histos[m+2];
      if(hMC==0 || hMC==nullptr){
	cout << "Null pointer" << endl;
	continue;
      }
      cout << "Pseudo-data is histo=" << hMC->GetName() << endl;
      
      // account for slight change in V for different mass hypos
      bin_counter = 0;
      for(unsigned int ix = 1; ix<=nx; ix++ ){
	for(unsigned int iy = 1; iy<=ny; iy++ ){
	  if(!accept_bin(hMC,ix,iy)) continue;
	  double val = hMC->GetBinContent(ix,iy);
	  if(do_toys){
	    val = -1.0;
	    while(val<0) val = mu_ran(bin_counter);
	  }
	  if(add_MC_uncert) val += (hMC->GetBinError(ix,iy)*hMC->GetBinError(ix,iy));
	  inv_sqrtV(bin_counter,bin_counter) = 1./TMath::Sqrt(val);
	  inv_V(bin_counter,bin_counter) = 1./val;
	  bin_counter++;
	}
      }
      
      VectorXd y(nbins);
      bin_counter = 0;
      for(unsigned int ix = 1; ix<=nx; ix++ ){
	for(unsigned int iy = 1; iy<=ny; iy++ ){
	  if(!accept_bin(hMC,ix,iy)) continue;
	  if(!do_toys)
	    y(bin_counter) = -all_histos[0]->GetBinContent(ix,iy)+hMC->GetBinContent(ix,iy);
	  else
	    y(bin_counter) = -mu_ran(bin_counter)+hMC->GetBinContent(ix,iy);
	  bin_counter++;
	}
      }

      //cout << inv_V << endl;
      MatrixXd A = inv_sqrtV*jac;
      MatrixXd b = inv_sqrtV*y;
      VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
      
      if(debug && false){
	double chi2_BDCSVD       = (b-A*A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b)).squaredNorm();
	double chi2_HouseholderQR = (b-A*A.householderQr().solve(b)).squaredNorm();
	double chi2_ColPivHouseholderQR = (b-A*A.colPivHouseholderQr().solve(b)).squaredNorm();
	double chi2_FullPivHouseholderQR = (b-A*A.fullPivHouseholderQr().solve(b)).squaredNorm();
	double chi2_LLT  = (b-A*((A.transpose()*A).llt().solve(A.transpose()*b))).squaredNorm();
	double chi2_LDLT = (b-A*((A.transpose()*A).ldlt().solve(A.transpose()*b))).squaredNorm();
	double chi2_JacobiSVD = (b-A*A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b)).squaredNorm();
	//double chi2_PartialPivLU = (b-A*A.partialPivLu().solve(b)).squaredNorm();
	//double chi2_FullPivLU = (b-A*A.fullPivLu().solve(b)).squaredNorm();
	//double chi2_CompleteOrthogonalDecomposition = (b-A*A.completeOrthogonalDecomposition().solve(b)).squaredNorm();
	cout << "DEBUG linear system for m_idx = " << m << endl;
	cout << "\tBDCSVD               = " << chi2_BDCSVD << endl;
	cout << "\tHouseholderQR        = " << chi2_HouseholderQR << endl;
	cout << "\tColPivHouseholderQR  = " << chi2_ColPivHouseholderQR << endl;
	cout << "\tFullPivHouseholderQR = " << chi2_FullPivHouseholderQR << endl;
	cout << "\tLLT                  = " << chi2_LLT << endl;
	cout << "\tLDLT                 = " << chi2_LDLT << endl;
	cout << "\tJacobiSVD            = " << chi2_JacobiSVD << endl;
	//cout << "\tPartialPivLU = " << chi2_PartialPivLU << endl;
	//cout << "\tFullPivLU = " << chi2_FullPivLU << endl;
	//cout << "\tCompleteOrthogonalDecomposition = " << chi2_CompleteOrthogonalDecomposition << endl;
      }
      
      MatrixXd C = (jac.transpose()*inv_V*jac).inverse();
      MatrixXd rho( C.rows(), C.rows() ) ;
      for(unsigned int ir = 0; ir<C.rows(); ir++){
	for(unsigned int ic = 0; ic<C.rows(); ic++){
	  rho(ir,ic) = C(ir,ic)/TMath::Sqrt(C(ir,ir)*C(ic,ic));
	}
      }
      //cout << rho << endl;
      TH2D* rho_th2 = 0;
      TH2D* A_th2 = 0;
      TH1D* b_th1 = 0;
      if(m<0){	
	rho_th2 = new TH2D("rho_th2"+toy_tag, ";POI;POI", rho.rows(), 0, rho.rows(), rho.cols(), 0, rho.cols());
	for(unsigned int i=0; i<rho.rows(); i++){
	  for(unsigned int j=0; j<rho.rows(); j++){
	    rho_th2->SetBinContent(i+1,j+1, rho(i,j));
	  }
	}
	A_th2 = new TH2D("A_th2", ";bin;POI", A.rows(), 0, A.rows(), A.cols(), 0, A.cols());
	for(unsigned int i=0; i<A.rows(); i++){
	  for(unsigned int j=0; j<A.cols(); j++){
	    A_th2->SetBinContent(i+1,j+1, A(i,j));
	  }
	}
	b_th1 = new TH1D("b_th1", ";bin;value", b.size(), 0, b.size());
	for(unsigned int i=0; i<b.size(); i++){
	  b_th1->SetBinContent(i+1, b(i));
	}      	
      }
      
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(C);
      if (eigensolver.info() != Eigen::Success){
	cout << "Could not eigendecompose C" << endl;
	abort();
      }
      if(verbose && m<0) std::cout << "The eigenvalues of C are:\n" << eigensolver.eigenvalues() << std::endl;
      //eigensolver.eigenvectors();

      TH2D* rhox_int_th2 = 0;
      VectorXd x_int;
      VectorXd xMC_int;
      MatrixXd Cx_int;
      if(run=="full"){
	int degs_xy = n_pdfx*n_pdfy;
	MatrixXd C_xy = C.block(0,0,degs_xy,degs_xy);
	//cout << C_xy << endl;
	VectorXd x_xy = x(Eigen::seqN(0,degs_xy));
	VectorXd xMC_xy(degs_xy);
	for(unsigned int i=0; i<degs_xy; i++ ) xMC_xy(i) = poi_val[ active_pois[i] ];
	MatrixXd D(n_pdfx, degs_xy);
	for(unsigned int i = 0 ; i<D.rows(); i++){
	  unsigned int j_first = i*n_pdfy; 
	  unsigned int j_last  = (i+1)*n_pdfy; 
	  for(unsigned int j = 0 ; j < D.cols() ; j++){
	    if(j>=j_first && j<j_last){ 
	      if(degs_corr_y==4)
		D(i,j) = norms_cheb4(j-j_first)*( j==(j_last-1) ? 1.0 : 2.0);
	      else if(degs_corr_y==6)
		D(i,j) = norms_cheb6(j-j_first)*( j==(j_last-1) ? 1.0 : 2.0);
	      else if(degs_corr_y==8)
		D(i,j) = norms_cheb8(j-j_first)*( j==(j_last-1) ? 1.0 : 2.0);
	      else
		D(i,j) = norms_chebDummy(j-j_first)*( j==(j_last-1) ? 1.0 : 2.0);
	    }
	    else 
	      D(i,j) = 0.0;
	  }
	}
	x_int   = D*x_xy;
	xMC_int = D*xMC_xy;
	Cx_int  = D*C_xy*D.transpose();
	
	MatrixXd rhox_int(Cx_int.rows(), Cx_int.cols());
	for(unsigned int i = 0 ; i<rhox_int.rows(); i++){
	  for(unsigned int j = 0 ; j < rhox_int.cols() ; j++){
	    rhox_int(i,j) = Cx_int(i,j)/TMath::Sqrt(Cx_int(i,i)*Cx_int(j,j));
	  }
	}
	if(m<0 && verbose){
	  cout << "Correlation matrix rhox_int:" << endl;
	  cout << rhox_int << endl;
	}
	
	if(m<0){
	  rhox_int_th2 = new TH2D("rhox_int_th2"+toy_tag, ";POI;POI", rhox_int.rows(), 0, rhox_int.rows(), rhox_int.cols(), 0, rhox_int.cols());
	  for(unsigned int i=0; i<rhox_int.rows(); i++){
	    for(unsigned int j=0; j<rhox_int.rows(); j++){
	      rhox_int_th2->SetBinContent(i+1,j+1, rhox_int(i,j));
	    }
	  }
	}
      }
      
      //cout << D << endl;
      //cout << xMC_int << endl;
      //cout << Cx_int << endl;
      
      MatrixXd chi2old = b.transpose()*b;
      MatrixXd chi2 = ((b - A*x).transpose())*(b-A*x);
      int ndof = nbins-poi_counter;
      double chi2norm = chi2(0,0)/ndof;

      chi2_start = chi2old(0,0);
      chi2_min = chi2(0,0);
      if(m>=0){
	mass_test =  (MW - DELTAM*0.5 + DELTAM/NMASS*m);
      }
      else if(jacmass==1) mass_test = MW + MASSSHIFT;
      else if(jacmass==2) mass_test = MW - MASSSHIFT;
      else mass_test = MW;
      tree2->Fill();
      
      if(m>=0){
	xx_mass[m] = MW - DELTAM*0.5 + DELTAM/NMASS*m; 
	exx_mass[m] = 0.0;
	yy_chi2[m] = chi2(0,0);
	eyy_chi2[m] = 0.0;
      }
      
      VectorXd pulls(x.size());
      for(unsigned int ip = 0; ip<pulls.size(); ip++){
	pulls(ip) = x(ip) / TMath::Sqrt(C(ip,ip));
      }
      //cout << pulls << endl;
      
      if(m<0 && verbose){
	for(unsigned int j = 0; j<poi_counter; j++){
	  unsigned int idx = active_pois[j];
	  cout << "POI " << idx << ": " << poi_val[idx] << " +/- " << TMath::Sqrt(C(j,j)) 
	       << " (pull = " << pulls(j) << ")" <<  endl;
	}
      }
      
      if(m<0) cout << "chi2     : " << chi2old(0,0) << " --> " << chi2 << "; ndof = " << ndof << " => chi2/ndof = " << chi2norm << endl; 
      
      fout->cd();
      if(m<0){
	rho_th2->Write();
	if(run=="full") rhox_int_th2->Write();
	A_th2->Write();
	b_th1->Write();
      }
      
      if(debug && m>=0){
	TH2D* hw_postfit = (TH2D*)hw->Clone(TString(Form("hw_postfit%d", m))+toy_tag);
	bin_counter = 0;
	for(unsigned int ix = 1; ix<=nx; ix++ ){
	  for(unsigned int iy = 1; iy<=ny; iy++ ){
	    if(!accept_bin(hw,ix,iy)) continue;
	    double val = hw->GetBinContent(ix,iy);
	    val += (jac*x)(bin_counter);
	    hw_postfit->SetBinContent(ix,iy,val);
	    bin_counter++;
	  }
	}
	hMC->Write(TString(Form("hw_prefit%d", m))+toy_tag);
	hw_postfit->Write(TString(Form("hw_postfit%d", m))+toy_tag);
      }

      if(verbose && m<0) cout << "Saving results to file..." << endl;
      std::vector<int> helicities = {-1, 0, 1, 2, 3, 4};
      for(auto hel : helicities) {
	vector< std::pair<unsigned int, unsigned int> > active;
	for(unsigned int i = 0; i < active_pois.size(); i++){
	  if(poi_cat[ active_pois[i] ]==hel){
	    active.emplace_back( std::make_pair(i,active_pois[i]) );
	  }
	}  
	unsigned int n = active.size();
	if(verbose && m<0){
	  cout << "\t" << n << " active POIs for hel=" << hel << endl;
	}
	double xx[n], yy[n], yyMC[n], exx[n], eyy[n];
	for(unsigned int i = 0; i < active.size(); i++){
	  auto p = active[i];
	  xx[i]  = p.second;
	  yy[i]  = poi_val[p.second] + x(p.first);
	  yyMC[i]= poi_val[p.second];
	  exx[i] = 0.0;
	  eyy[i] = TMath::Sqrt(C(p.first,p.first));
	}
	fout->cd();
	TGraphErrors* fit   = new TGraphErrors(n,xx,yy,exx,eyy);
	TGraphErrors* fitMC = new TGraphErrors(n,xx,yyMC,exx,exx);
	string name = hel==-1 ? "corrxy" : std::string(Form("A%d", hel));
	
	if(m>=0) name += std::string(Form("_mass%d", m));
	
	if(m<0) fit->Write(("fit_"+name).c_str());
	if(m<0) fitMC->Write(("fitMC_"+name).c_str());
	
	if(run=="full" && hel==-1){
	  n = n_pdfx;
	  double xx[n], yy[n], yyMC[n], exx[n], eyy[n];
	  for(unsigned int i = 0; i < n_pdfx; i++){
	    xx[i]  = points_x[i];
	    yy[i]  = xMC_int(i)+x_int(i); 
	    yyMC[i]= xMC_int(i);
	    exx[i] = 0.0;
	    eyy[i] = TMath::Sqrt( Cx_int(i,i) );
	  }
	  TGraphErrors* fit   = new TGraphErrors(n,xx,yy,exx,eyy);
	  TGraphErrors* fitMC = new TGraphErrors(n,xx,yyMC,exx,exx);
	  string name = "x_inty";
	  if(m>=0) name += std::string(Form("_mass%d", m));
	  if(m<0) fit->Write(("fit_"+name).c_str());
	  if(m<0) fitMC->Write(("fitMC_"+name).c_str());      
	}
      }
    }
    
    //hwMC->Write("hw_nominal");
    if(NMASS>1){
      TGraphErrors* chi2_fit = new TGraphErrors(NMASS,xx_mass,yy_chi2,exx_mass,eyy_chi2);
      chi2_fit->Write("chi2_vs_mass"+toy_tag);
      
      chi2_fit->Fit("pol2", "Q");
      TF1* parabola = chi2_fit->GetFunction("pol2");
      float param0 = parabola->GetParameter(0); 
      float param1 = parabola->GetParameter(1); 
      float param2 = parabola->GetParameter(2); 
      float deltaM = 1./TMath::Sqrt(param2);
      float biasM = -param1/param2*0.5;
      float pullM = (biasM-MW)/deltaM; 
      if(do_toys) cout << itoy << ": ";
      cout << "DeltaM = " << deltaM*1e+03 << " MeV" 
	   << " -- bias: " << (biasM-MW)*1e+03  << " MeV"
	   << ", pull: " << pullM 
	 << endl;
      
      int idx = 0;
      for(int m=0; m<NMASS; m++){
	float mass_m = MW - DELTAM*0.5 + DELTAM/NMASS*m;
	if(mass_m > (biasM - deltaM)){
	  idx = m;
	  break;
	}
      }
      cout << "Intersection at mass index [" << idx-1 << "," << idx << "]" << endl; 
      best_mass     = (biasM-MW)*1e+03;
      best_mass_err = deltaM*1e+03;
      tree->Fill();
    }
    
  }

  tree->Write();
  tree2->Write();
  fout->Close();
  fin->Close();

  sw.Stop();
  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;
  return 1;
}
