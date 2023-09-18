#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2D.h"
#include "TF1.h"
#include "TF2.h"
#include "TFitResultPtr.h"
#include "TRandom3.h"
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
constexpr int NMAX  = 1000;
constexpr int NMASS = 3;
constexpr double DELTAM = 0.200;

int main(int argc, char* argv[])
{

  TRandom3* ran = new TRandom3();
  
  TStopwatch sw;
  sw.Start();

  variables_map vm;
  try
    {
      options_description desc{"Options"};
      desc.add_options()
	("help,h", "Help screen")
	("nevents",     value<long>()->default_value(1000), "number of events")
	("degs_corr_x", value<int>()->default_value(2), "max degree in x of corrxy")
	("degs_corr_y", value<int>()->default_value(2), "max degree in y of corrxy")
	("rebinX",   value<int>()->default_value(-1), "rebin X axis")
	("rebinY",   value<int>()->default_value(-1), "rebin Y axis")
	("prior",   value<float>()->default_value(0.5), "prior")
	("jUL",   bool_switch()->default_value(false), "")
	("j0",    bool_switch()->default_value(false), "")
	("j1",    bool_switch()->default_value(false), "")
	("j2",    bool_switch()->default_value(false), "")
	("j3",    bool_switch()->default_value(false), "")
	("j4",    bool_switch()->default_value(false), "")
	("scale0",bool_switch()->default_value(false), "")
	("scale1",bool_switch()->default_value(false), "")
	("scale2",bool_switch()->default_value(false), "")
	("scale3",bool_switch()->default_value(false), "")
	("scale4",bool_switch()->default_value(false), "")
	("verbose", bool_switch()->default_value(false), "")
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
  std::string tag = vm["tag"].as<std::string>();
  std::string post_tag = vm["post_tag"].as<std::string>();
  std::string run = vm["run"].as<std::string>();
  int degs_corr_x = vm["degs_corr_x"].as<int>();
  float prior_sigma = vm["prior"].as<float>();
  int degs_corr_y = vm["degs_corr_y"].as<int>();
  int rebinX      = vm["rebinX"].as<int>();
  int rebinY      = vm["rebinY"].as<int>();
  int jUL  = vm["jUL"].as<bool>();
  int j0   = vm["j0"].as<bool>();
  int j1   = vm["j1"].as<bool>();
  int j2   = vm["j2"].as<bool>();
  int j3   = vm["j3"].as<bool>();
  int j4   = vm["j4"].as<bool>();
  int verbose = vm["verbose"].as<bool>();

  if(vm.count("degs_corr_x")) tag += std::string(Form("_UL_%d", degs_corr_x));
  if(vm.count("degs_corr_y")) tag += std::string(Form("_%d", degs_corr_y));

  TFile* fin = TFile::Open(("root/histos_"+tag+"_"+run+".root").c_str(), "READ");
  if(fin==0 || fin==nullptr || fin->IsZombie()){
    cout << "File NOT found" << endl;
    return 0;
  }

  Long64_t N;
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
  outtree->SetBranchAddress("poi_cat", &poi_cat);
  outtree->SetBranchAddress("poi_idx", &poi_idx);
  outtree->SetBranchAddress("n_pdfx",   &n_pdfx);
  outtree->SetBranchAddress("n_pdfy",   &n_pdfy);
  outtree->SetBranchAddress("points_x", &points_x);
  outtree->SetBranchAddress("points_y", &points_y);
  outtree->GetEntry(0);    
  
  std::vector<unsigned int> active_pois = {};
  for(unsigned int p = 0; p < poi_counter; p++){
    if(verbose) cout << p << ", " << poi_cat[p] << endl;
    if(poi_cat[p]==-1 && jUL) active_pois.emplace_back(p);
    if(poi_cat[p]==0  && j0)  active_pois.emplace_back(p);
    if(poi_cat[p]==1  && j1)  active_pois.emplace_back(p);
    if(poi_cat[p]==2  && j2)  active_pois.emplace_back(p);      
    if(poi_cat[p]==3  && j3)  active_pois.emplace_back(p);      
    if(poi_cat[p]==4  && j4)  active_pois.emplace_back(p);      
  }
  poi_counter = active_pois.size();
  if(verbose){
    for(auto p : active_pois) cout << p << ", ";
    cout << endl;
  }

  std::vector<unsigned int> cleaned_active_pois = {};
  std::vector<int> helicities = {-1, 0, 1, 2, 3, 4};
  for(auto hel : helicities){
    TH2D* hjac_first = 0;
    for(auto p : active_pois){
      if(hel<0 && poi_cat[p]==hel){
	cleaned_active_pois.emplace_back(p);
	continue;
      }
      else if( poi_cat[p]==hel && vm[std::string(Form("scale%d", hel))].as<bool>()  ){
	if(hjac_first==0){
	  if(verbose) cout << "Found first jacobian of A" << hel << endl;
	  hjac_first = fin->Get<TH2D>(Form("jac_%d", p));
	  cleaned_active_pois.emplace_back(p);
	}
	else{
	  if(verbose) cout << "Found >1 jacobian of A" << hel << ": adding up..." << endl;
	  hjac_first->Add(fin->Get<TH2D>(Form("jac_%d", p)));
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
    cout << "After cleaning:" << endl;
    for(auto p : active_pois) cout << p << ", ";
    cout << endl;
  }
  poi_counter = active_pois.size();


  
  // common
  TH2D* hw   = fin->Get<TH2D>("h");
  TH2D* hwMC = fin->Get<TH2D>("hMC");
  vector<TH2D*> all_histos = { hw, hwMC };
  for(unsigned int i=0; i<NMASS; i++) all_histos.push_back( fin->Get<TH2D>(Form("hMC_mass%d",i)) );
  if(rebinX>0 || rebinY>0){
    for(auto h : all_histos) h->Rebin2D(rebinX,rebinY);
  }
  double lumi_fact = nevents>0 ? nevents/hwMC->Integral() : 1.0; 
  cout << "Scaling input histos by " << lumi_fact << endl;  
  N *= lumi_fact;
  for(auto h : all_histos) h->Scale(lumi_fact);
  int nx = hw->GetXaxis()->GetNbins(); 
  int ny = hw->GetYaxis()->GetNbins(); 
  int nbins = nx*ny;
  cout << "Number of data points: " << nx << "*" << ny << " = " << nbins << endl;
  cout << "Data integral: " << hwMC->Integral() << endl;
  
  unsigned int bin_counter = 0;

  MatrixXd jac(nbins, poi_counter);
  for(unsigned int j = 0; j<poi_counter; j++){
    unsigned int idx = active_pois[j];
    TH2D* hjac = fin->Get<TH2D>(Form("jac_%d", idx));
    if(rebinX>0 || rebinY>0){
      hjac->Rebin2D(rebinX,rebinY);
    }
    bin_counter = 0;
    for(unsigned int ix = 1; ix<=nx; ix++ ){
      for(unsigned int iy = 1; iy<=ny; iy++ ){
	jac(bin_counter,j) = N*hjac->GetBinContent(ix,iy); 
	bin_counter++;
      }
    }
  }

  //cout << "Jacobians filled." << endl;

  VectorXd y0(nbins);
  VectorXd x0(poi_counter);
  double x0norm = 1.0;
  MatrixXd invVp0 =  MatrixXd::Zero(poi_counter,poi_counter);
  MatrixXd Vp0 =  MatrixXd::Zero(poi_counter,poi_counter);
  MatrixXd inv_sqrtV(nbins, nbins);
  MatrixXd inv_V(nbins, nbins);
  for(unsigned int ix = 0; ix<nbins; ix++ ){
    for(unsigned int iy = 0; iy<nbins; iy++ ){
      inv_sqrtV(ix,iy) = 0.;
      inv_V(ix,iy) = 0.;
    }
  }
  bin_counter = 0;
  for(unsigned int ix = 1; ix<=nx; ix++ ){
    for(unsigned int iy = 1; iy<=ny; iy++ ){
      inv_sqrtV(bin_counter,bin_counter) = 1./TMath::Sqrt(hwMC->GetBinContent(ix,iy));
      inv_V(bin_counter,bin_counter) = 1./hwMC->GetBinContent(ix,iy);
      y0(bin_counter) = all_histos[2]->GetBinContent(ix,iy) - all_histos[0]->GetBinContent(ix,iy);
      bin_counter++;
    }
  }
  //cout << "V matrix filled" << endl;
  
  
  // plotting
  vector<string> fit_tags = {"UL", "0", "1", "2", "3", "4" };
  string fit_tag = "";
  for(auto t : fit_tags){
    if(vm["j"+t].as<bool>()) fit_tag += "_j"+t;
  }
  TFile* fout = TFile::Open(("root/fit_"+tag+"_"+run+fit_tag+"_"+post_tag+".root").c_str(), "RECREATE");
  //cout << "Commons done. Output file opened." << endl;

  TH1D* hx0  = new TH1D("x0", "", poi_counter, 0, poi_counter);
  TH2D* hVp0 = new TH2D("Vp0", "", poi_counter, 0, poi_counter, poi_counter, 0, poi_counter);
  TH2D* hCp0 = new TH2D("Cp0", "", poi_counter, 0, poi_counter, poi_counter, 0, poi_counter);
  
  double xx_mass[NMASS], yy_chi2[NMASS], yy_chi22[NMASS], exx_mass[NMASS], eyy_chi2[NMASS];
  // mass specific
  for(int m=-1; m<NMASS; m++){

    TH2D* hMC = all_histos[m+2];
    if(hMC==0 || hMC==nullptr){
      cout << "Null pointer" << endl;
      continue;
    }

    VectorXd y(nbins);
    bin_counter = 0;
    for(unsigned int ix = 1; ix<=nx; ix++ ){
      for(unsigned int iy = 1; iy<=ny; iy++ ){
	y(bin_counter) = all_histos[0]->GetBinContent(ix,iy)-hMC->GetBinContent(ix,iy);
	//y(bin_counter) = hMC->GetBinContent(ix,iy);
	bin_counter++;
      }
    }
    
    MatrixXd invVp = MatrixXd::Zero(poi_counter,poi_counter);
    MatrixXd Vp = MatrixXd::Zero(poi_counter,poi_counter);
    for(unsigned int ix = 0; ix<poi_counter ; ix++){
      for(unsigned int iy = 0; iy<poi_counter ; iy++){
	if(ix==iy){
	  invVp(ix,iy) = 1./(prior_sigma*prior_sigma);
	  Vp(ix,iy) = (prior_sigma*prior_sigma);
	}
      }
    }   

    //cout << inv_V << endl;
    MatrixXd A = inv_sqrtV*jac;
    MatrixXd b = inv_sqrtV*y;

    MatrixXd B = A.transpose()*A + invVp;
    VectorXd g = -A.transpose()*b;

    //VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
    VectorXd x = B.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(-g);        
    MatrixXd prior = x.transpose()*invVp*x;

    if(m==-1){
      cout << "Doing x0" << endl;
      x0 = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(inv_sqrtV*y0);
      for(unsigned int ix=0; ix<poi_counter;ix++){
	hx0->SetBinContent(ix+1,x0(ix));
      }

      //x0 = Eigen::VectorXd::Zero(poi_counter);
      //x0(0) = 1.;
      //x0(1) = -1.;
      
      x0norm = x0.norm();
      cout << "x0norm = " << x0norm << endl;
      // to be checked
      x0 /= x0norm;
      std::vector<VectorXd> new_eigenvecs;
      new_eigenvecs.emplace_back(x0);
      //for(unsigned int iv=0; iv<poi_counter-1; iv++){
      for(unsigned int iv=1; iv<poi_counter; iv++){
	//cout << "Doing eigenvector " << iv << endl;
	VectorXd vi = Eigen::VectorXd::Zero(poi_counter);
	vi(iv) = 1.0;
	VectorXd ui = vi;
	for(unsigned int jv=0; jv<new_eigenvecs.size(); jv++){
	  VectorXd update = vi.dot(new_eigenvecs[jv])*new_eigenvecs[jv];
	  ui -= update;
	}
	if(ui.norm()>0.){
	  ui /= ui.norm();
	}
	cout << iv << " norm is " << ui.norm() << endl;	  
	new_eigenvecs.emplace_back(ui);
      }
      //cout << "Check orthogonality:" << endl;
      for(unsigned int jv=0; jv<new_eigenvecs.size(); jv++){
	cout << "eigen " << jv << ": " << new_eigenvecs[jv] << endl;
      }
      //for(unsigned int jv=0; jv<new_eigenvecs.size(); jv++){
      //for(unsigned int kv=0; kv<new_eigenvecs.size(); kv++){
      //  cout << "(" << jv << "," << kv << "): " << new_eigenvecs[jv].dot(new_eigenvecs[kv]) << endl;
      //}      
      //}
      //x0norm = 1.0;
      for(unsigned int jv=0; jv<new_eigenvecs.size(); jv++){
	invVp0 += (jv==0 ? 1./prior_sigma/prior_sigma/x0norm/x0norm :
		   1.0/prior_sigma/prior_sigma)*new_eigenvecs[jv]*new_eigenvecs[jv].transpose();
	Vp0 += (jv==0 ? prior_sigma*prior_sigma*x0norm*x0norm :
		prior_sigma*prior_sigma )*new_eigenvecs[jv]*new_eigenvecs[jv].transpose();
      }

      //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver0(invVp0);
      //std::cout << "The eigenvalues of invVp0 are:\n" << eigensolver0.eigenvalues() << std::endl;
      //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver1(Vp0);
      //std::cout << "The eigenvalues of Vp0 are:\n" << eigensolver1.eigenvalues() << std::endl;
      //std::cout << "The eigenvector 0 of invVp0 are:\n" << std::endl;
      //cout << eigensolver0.eigenvectors().col(0) << endl;
      //std::cout << "The eigenvector n-1 of Vp0 are:\n" << std::endl;
      //cout << eigensolver1.eigenvectors().col(poi_counter-1) << endl;
      
      //MatrixXd Vp0 = invVp0.inverse();
      MatrixXd Vp0scaled = MatrixXd::Zero(poi_counter,poi_counter); 
      for(unsigned int xp = 0 ; xp<poi_counter; xp++){
	for(unsigned int yp = 0 ; yp<poi_counter; yp++){
	  //Vp0scaled(xp,yp) = Vp0(xp,yp)/TMath::Sqrt(Vp0(xp,xp)*Vp0(yp,yp))*prior_sigma*prior_sigma;
	  Vp0scaled(xp,yp) = xp==yp ? Vp0(xp,yp) : 0.0;
	}
      }
      Vp0 = Vp0scaled;
      invVp0 = Vp0.inverse();


      //Vp0 = invVp0.inverse();
      //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver2(Vp0);
      //std::cout << "The eigenvalues of invVp0^-1 are:\n" << eigensolver2.eigenvalues() << std::endl;
      
      for(unsigned int xp = 0 ; xp<poi_counter; xp++){
	for(unsigned int yp = 0 ; yp<poi_counter; yp++){
	  hVp0->SetBinContent(xp+1,yp+1,Vp0(xp,yp));
	  hCp0->SetBinContent(xp+1,yp+1, Vp0(xp,yp)/TMath::Sqrt(Vp0(xp,xp))/TMath::Sqrt(Vp0(yp,yp)) );
	  if(xp==yp) cout << "(" << xp << "," << yp << ") => " << TMath::Sqrt(Vp0(xp,yp)) << endl;
	}
      }
    }
    
    MatrixXd invCp = (A.transpose()*A + invVp);
    MatrixXd Cp = invCp.inverse();
    MatrixXd Vp2 = Cp;
    for(unsigned int xp = 0 ; xp<poi_counter; xp++){
      for(unsigned int yp = 0 ; yp<poi_counter; yp++){
	if(xp==yp)
	  Vp2(xp,yp) = Vp(xp,yp);
	else{
	  //Vp2(xp,yp) = Cp(xp,yp)/TMath::Sqrt(Cp(xp,xp)*Cp(yp,yp))*TMath::Sqrt(Vp(xp,xp)*Vp(yp,yp));
	  Vp2(xp,yp) = ran->Uniform(-0.99,0.99)*TMath::Sqrt(Vp(xp,xp)*Vp(yp,yp));
	}
      }
    }
    //MatrixXd invVp2 = Vp2.inverse();
    MatrixXd invVp2 = invVp0;

    MatrixXd B2 = A.transpose()*A + invVp2;
    VectorXd x2 = B2.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(-g);        
    MatrixXd prior2 = x2.transpose()*invVp2*x2;


    MatrixXd C = (jac.transpose()*inv_V*jac).inverse();
    MatrixXd rho( C.rows(), C.rows() ) ;
    for(unsigned int ir = 0; ir<C.rows(); ir++){
      for(unsigned int ic = 0; ic<C.rows(); ic++){
	rho(ir,ic) = C(ir,ic)/TMath::Sqrt(C(ir,ir)*C(ic,ic));
      }
    }
    if(verbose && m<0) cout << rho << endl;
    
    TH2D* rho_th2 = 0;
    TH2D* A_th2 = 0;
    TH1D* b_th1 = 0;
    if(m<0){
      rho_th2 = new TH2D("rho_th2", ";POI;POI", rho.rows(), 0, rho.rows(), rho.cols(), 0, rho.cols());
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
   
    MatrixXd chi2old = b.transpose()*b;
    MatrixXd chi2 = ((b - A*x).transpose())*(b-A*x) + prior;
    MatrixXd chi22 = ((b - A*x2).transpose())*(b-A*x2) + prior2;
    int ndof = nbins-poi_counter;
    double chi2norm = chi2(0,0)/ndof;

    if(m>=0){
      xx_mass[m] = MW - DELTAM*0.5 + DELTAM/NMASS*m; 
      exx_mass[m] = 0.0;
      yy_chi2[m] = chi2(0,0);
      yy_chi22[m] = chi22(0,0);
      eyy_chi2[m] = 0.0;
    }

    VectorXd pulls(x.size());
    for(unsigned int ip = 0; ip<pulls.size(); ip++){
      pulls(ip) = x(ip) / TMath::Sqrt(C(ip,ip));
      //cout << x(ip) << " +/- " << TMath::Sqrt(C(ip,ip)) << endl;
    }
    //cout << pulls << endl;

    if(m<0){
      for(unsigned int ip = 0; ip<pulls.size(); ip++){
	cout << ip << ": " <<  TMath::Sqrt(1./invVp(ip,ip)) << " --> " << TMath::Sqrt(Cp(ip,ip)) << endl;
      }
    }
    if(m<0 || true){
      cout << "chi2     : " << chi2old(0,0) << " --> " << chi2 << "; ndof = " << ndof << " => chi2/ndof = " << chi2norm << endl; 
      cout << "From prior += " << prior(0,0) << endl;
    }
    {
      cout << "chi22     : " << chi2old(0,0) << " --> " << chi22 << endl;
      cout << "From prior += " << prior2(0,0) << endl;
    }

    fout->cd();
    if(m<0){
      rho_th2->Write();
      A_th2->Write();
      b_th1->Write();
    }

    //std::vector<int> helicities = {-1, 0, 1, 2, 3, 4};
    for(auto hel : helicities) {
      if(verbose && m<0) cout << "helicity=" << hel << endl;
      vector< std::pair<unsigned int, unsigned int> > active;
      for(unsigned int i = 0; i < active_pois.size(); i++){
	if(poi_cat[ active_pois[i] ]==hel){
	  if(verbose && m<0) cout << "(" << i << "," << active_pois[i] << ")" << endl; 
	  active.emplace_back( std::make_pair(i,active_pois[i]) );
	}
      }
      unsigned int n = active.size();
      double xx[n], yy[n], yyMC[n], exx[n], eyy[n];
      for(unsigned int i = 0; i < active.size(); i++){
	auto p = active[i];
	xx[i]  = p.second;
	yy[i]  = x(p.first);
	yyMC[i]= 0.0;
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
    }
    if(m<0){
      hx0->Write();
      hVp0->Write();
      hCp0->Write();
    }
      
  }

  TGraphErrors* chi2_fit = new TGraphErrors(NMASS,xx_mass,yy_chi2,exx_mass,eyy_chi2);
  chi2_fit->Write("chi2_vs_mass");
  TGraphErrors* chi22_fit = new TGraphErrors(NMASS,xx_mass,yy_chi22,exx_mass,eyy_chi2);
  chi22_fit->Write("chi22_vs_mass");
  fout->Close();
  fin->Close();

  chi2_fit->Fit("pol2", "Q");
  TF1* parabola = chi2_fit->GetFunction("pol2");
  float param0 = parabola->GetParameter(0); 
  float param1 = parabola->GetParameter(1); 
  float param2 = parabola->GetParameter(2); 
  float deltaM = 1./TMath::Sqrt(param2);
  float biasM = -param1/param2*0.5;
  float pullM = (MW-biasM)/deltaM; 
  cout << "from 1: dM = " << deltaM*1e+03 << " MeV" 
       << " -- bias: " << (MW-biasM)*1e+03  << " MeV"
       << ", pull: " << pullM 
       << endl;

  chi22_fit->Fit("pol2", "Q");
  TF1* parabola2 = chi22_fit->GetFunction("pol2");
  float param02 = parabola2->GetParameter(0); 
  float param12 = parabola2->GetParameter(1); 
  float param22 = parabola2->GetParameter(2); 
  float deltaM2 = 1./TMath::Sqrt(param22);
  float biasM2 = -param12/param22*0.5;
  float pullM2 = (MW-biasM2)/deltaM2; 
  cout << "from 2: dM = " << deltaM2*1e+03 << " MeV" 
       << " -- bias: " << (MW-biasM2)*1e+03  << " MeV"
       << ", pull: " << pullM2 
       << endl;

  sw.Stop();
  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;
  return 1;
}
