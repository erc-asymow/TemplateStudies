#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
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
constexpr int NMAX  = 100;
constexpr int NMASS = 50;

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
  int jUL  = vm["jUL"].as<bool>();
  int j0   = vm["j0"].as<bool>();
  int j1   = vm["j1"].as<bool>();
  int j2   = vm["j2"].as<bool>();
  int j3   = vm["j3"].as<bool>();
  int j4   = vm["j4"].as<bool>();
  int verbose = vm["verbose"].as<bool>();

  //if(vm.count("degs_pdf_x"))  tag += std::string(Form("_%d", degs_pdf_x));
  //if(vm.count("degs_pdf_y"))  tag += std::string(Form("_%d", degs_pdf_y));
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
  
  std::vector<unsigned int> active_pois = {};
  for(unsigned int p = 0; p < poi_counter; p++){
    //cout << p << ", " << poi_cat[p] << endl;
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
      bin_counter++;
    }
  }
  //cout << "V matrix filled" << endl;

  VectorXd norms_cheb4(5);
  norms_cheb4 << 0.0332073, 0.266696, 0.400194, 0.266696, 0.0332073;
  VectorXd norms_cheb6(7);
  norms_cheb6 << 0.0141919, 0.126953, 0.228635, 0.26044, 0.228635, 0.126953, 0.0141919; 
  assert(degs_corr_y==4 || degs_corr_y==6);

  // plotting
  vector<string> fit_tags = {"UL", "0", "1", "2", "3", "4" };
  string fit_tag = "";
  for(auto t : fit_tags){
    if(vm["j"+t].as<bool>()) fit_tag += "_j"+t;
  }
  TFile* fout = TFile::Open(("root/fit_"+tag+"_"+run+fit_tag+"_"+post_tag+".root").c_str(), "RECREATE");
  //cout << "Commons done. Output file opened." << endl;


  double xx_mass[NMASS], yy_chi2[NMASS], exx_mass[NMASS], eyy_chi2[NMASS];
  // mass specific
  for(int m=-1; m<NMASS; m++){

    //cout << "Now doing mass idx " << m << endl;
    //TH2D* hMC = m<0 ? all_histos[1] : fin->Get<TH2D>(Form("hMC_mass%d",m));
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
	bin_counter++;
      }
    }
    
    //cout << inv_V << endl;
    MatrixXd A = inv_sqrtV*jac;
    MatrixXd b = inv_sqrtV*y;
    VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
    MatrixXd C = (jac.transpose()*inv_V*jac).inverse();
    MatrixXd rho( C.rows(), C.rows() ) ;
    for(unsigned int ir = 0; ir<C.rows(); ir++){
      for(unsigned int ic = 0; ic<C.rows(); ic++){
	rho(ir,ic) = C(ir,ic)/TMath::Sqrt(C(ir,ir)*C(ic,ic));
      }
    }
    //cout << rho << endl;
    TH2D* rho_th2 = 0;
    if(m<0){
      rho_th2 = new TH2D("rho_th2", ";POI;POI", rho.rows(), 0, rho.rows(), rho.cols(), 0, rho.cols());
      for(unsigned int i=0; i<rho.rows(); i++){
	for(unsigned int j=0; j<rho.rows(); j++){
	  rho_th2->SetBinContent(i+1,j+1, rho(i,j));
	}
      }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(C);
    if (eigensolver.info() != Eigen::Success){
      cout << "Could not eigendecompose C" << endl;
      abort();
    }
    if(verbose && m<0) std::cout << "The eigenvalues of C are:\n" << eigensolver.eigenvalues() << std::endl;
    //eigensolver.eigenvectors();

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
	if(j>=j_first && j<j_last) 
	  D(i,j) = degs_corr_y==4 ? norms_cheb4(j-j_first)*( j==(j_last-1) ? 1.0 : 2.0) : norms_cheb6(j-j_first)*( j==(j_last-1) ? 1.0 : 2.0);
	else 
	  D(i,j) = 0.0;
      }
    }
    VectorXd x_int   = D*x_xy;
    VectorXd xMC_int = D*xMC_xy;
    MatrixXd Cx_int  = D*C_xy*D.transpose();

    MatrixXd rhox_int(Cx_int.rows(), Cx_int.cols());
    for(unsigned int i = 0 ; i<rhox_int.rows(); i++){
      for(unsigned int j = 0 ; j < rhox_int.cols() ; j++){
	rhox_int(i,j) = Cx_int(i,j)/TMath::Sqrt(Cx_int(i,i)*Cx_int(j,j));
      }
    }
    if(m<0 && verbose) cout << rhox_int << endl;

    TH2D* rhox_int_th2 = 0;
    if(m<0){
      rhox_int_th2 = new TH2D("rhox_int_th2", ";POI;POI", rhox_int.rows(), 0, rhox_int.rows(), rhox_int.cols(), 0, rhox_int.cols());
      for(unsigned int i=0; i<rhox_int.rows(); i++){
	for(unsigned int j=0; j<rhox_int.rows(); j++){
	  rhox_int_th2->SetBinContent(i+1,j+1, rhox_int(i,j));
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

    if(m>=0){
      xx_mass[m] = MW - 0.100 + 0.200/NMASS*m; 
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
	cout << "poi " << idx << ": " << poi_val[idx] << " +/- " << TMath::Sqrt(C(j,j)) 
	     << ". Pull = " << pulls(j) << endl;
      }
    }

    if(m<0) cout << "chi2     : " << chi2old(0,0) << " --> " << chi2 << "; ndof = " << ndof << " => chi2/ndof = " << chi2norm << endl; 

    fout->cd();
    if(m<0){
      rho_th2->Write();
      rhox_int_th2->Write();
    }

    std::vector<int> helicities = {-1, 0, 1, 2, 3, 4};
    for(auto hel : helicities) {
      vector< std::pair<unsigned int, unsigned int> > active;
      for(unsigned int i = 0; i < active_pois.size(); i++){
	if(poi_cat[ active_pois[i] ]==hel){
	  active.emplace_back( std::make_pair(i,active_pois[i]) );
	}
      }  
      unsigned int n = active.size();
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
    
      if(hel==-1){
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

  TGraphErrors* chi2_fit = new TGraphErrors(NMASS,xx_mass,yy_chi2,exx_mass,eyy_chi2);
  chi2_fit->Write("chi2_vs_mass");
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
  cout << "dM = " << deltaM*1e+03 << " MeV" 
       << " -- bias: " << (MW-biasM)*1e+03  << " MeV"
       << ", pull: " << pullM 
       << endl;
  
  sw.Stop();
  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;
  return 1;
}
