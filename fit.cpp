#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2D.h"
#include "TF1.h"
#include "TF2.h"
#include <TStopwatch.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <boost/program_options.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

using namespace boost::program_options;

constexpr int NMAX  = 100;

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
	("j0",    bool_switch()->default_value(false), "")
	("j1",    bool_switch()->default_value(false), "")
	("j2",    bool_switch()->default_value(false), "")
	("tag", value<std::string>()->default_value(""), "tag name")
	("run", value<std::string>()->default_value("closure"), "run type");

      store(parse_command_line(argc, argv, desc), vm);
      notify(vm);
      if (vm.count("help")){
	std::cout << desc << '\n';
	return 0;
      }
      if (vm.count("nevents"))    std::cout << "Number of events: " << vm["nevents"].as<long>() << '\n';
      if (vm.count("tag"))        std::cout << "Tag: " << vm["tag"].as<std::string>() << '\n';
      if (vm.count("run"))        std::cout << "Run: " << vm["run"].as<std::string>() << '\n';
      if (vm.count("degs_pdf_x")) std::cout << "Degree of pdf_x: " << vm["degs_pdf_x"].as<int>() << '\n';
      if (vm.count("degs_pdf_y")) std::cout << "Degree of pdf_y: " << vm["degs_pdf_y"].as<int>() << '\n';
      if (vm.count("degs_corr_y")) std::cout << "Degree in x of corrxy: " << vm["degs_corr_x"].as<int>() << '\n';
      if (vm.count("degs_corr_y")) std::cout << "Degree in y of corrxy: " << vm["degs_corr_y"].as<int>() << '\n';
      if (vm.count("degs_A0_x")) std::cout << "Degree in x of A0: " << vm["degs_A0_x"].as<int>() << '\n';
      if (vm.count("degs_A0_y")) std::cout << "Degree in y of A0: " << vm["degs_A0_y"].as<int>() << '\n';
    }
  catch (const error &ex)
    {
      std::cerr << ex.what() << '\n';
    }

  long nevents    = vm["nevents"].as<long>();
  std::string tag = vm["tag"].as<std::string>();
  std::string run = vm["run"].as<std::string>();
  int degs_pdf_x  = vm["degs_pdf_x"].as<int>();
  int degs_pdf_y  = vm["degs_pdf_y"].as<int>();
  int degs_corr_x = vm["degs_corr_x"].as<int>();
  int degs_corr_y = vm["degs_corr_y"].as<int>();
  int degs_A0_x   = vm["degs_A0_x"].as<int>();
  int degs_A0_y   = vm["degs_A0_y"].as<int>();
  int j0   = vm["j0"].as<bool>();
  int j1   = vm["j1"].as<bool>();
  int j2   = vm["j2"].as<bool>();

  if(vm.count("degs_pdf_x"))  tag += std::string(Form("_%d", degs_pdf_x));
  if(vm.count("degs_pdf_y"))  tag += std::string(Form("_%d", degs_pdf_y));
  if(vm.count("degs_corr_x")) tag += std::string(Form("_%d", degs_corr_x));
  if(vm.count("degs_corr_y")) tag += std::string(Form("_%d", degs_corr_y));
  if(vm.count("degs_A0_x"))   tag += std::string(Form("_%d", degs_A0_x));
  if(vm.count("degs_A0_y"))   tag += std::string(Form("_%d", degs_A0_y));

  TFile* fout = TFile::Open(("root/histos_"+tag+"_"+run+".root").c_str(), "READ");
  if(fout==0 || fout==nullptr || fout->IsZombie()){
    cout << "File NOT found" << endl;
    return 0;
  }

  Long64_t N;
  double poi_val[NMAX];
  unsigned int poi_cat[NMAX];
  unsigned int poi_idx[NMAX];
  unsigned int poi_counter;
  TTree* outtree = fout->Get<TTree>("outtree");
  outtree->SetBranchAddress("nevents", &N);
  outtree->SetBranchAddress("poi_counter", &poi_counter);
  outtree->SetBranchAddress("poi_val", &poi_val);
  outtree->SetBranchAddress("poi_cat", &poi_cat);
  outtree->SetBranchAddress("poi_idx", &poi_idx);
  outtree->GetEntry(0);    
  
  std::vector<unsigned int> active_pois = {};
  for(unsigned int p = 0; p < poi_counter; p++){
    if(poi_cat[p]==0 && j0) active_pois.emplace_back(p);
    if(poi_cat[p]==1 && j1) active_pois.emplace_back(p);
    if(poi_cat[p]==2 && j2){
      active_pois.emplace_back(p);
    }
  }
  poi_counter = active_pois.size();
  for(auto p : active_pois) cout << p << ", " << endl;

  TH2D* hw = fout->Get<TH2D>("w");
  TH2D* hwMC = fout->Get<TH2D>("wMC");
  
  double lumi_fact = nevents>0 ? nevents/hwMC->Integral() : 1.0; 
  cout << "Scaling input histos by " << lumi_fact << endl;
  N *= lumi_fact;
  hw->Scale(lumi_fact);
  hwMC->Scale(lumi_fact);

  int nx = hw->GetXaxis()->GetNbins(); 
  int ny = hw->GetYaxis()->GetNbins(); 
  int nbins = nx*ny;

  //MatrixXd jac(nbins, poi_counter+1);
  MatrixXd jac(nbins, poi_counter);
  VectorXd y(nbins);
  MatrixXd inv_sqrtV(nbins, nbins);
  MatrixXd inv_V(nbins, nbins);
  inv_sqrtV *= 0.;
  inv_V *= 0.;

  unsigned int bin_counter = 0;
  for(unsigned int ix = 1; ix<=nx; ix++ ){
    for(unsigned int iy = 1; iy<=ny; iy++ ){
      //jac(bin_counter,0) = hw->GetBinContent(ix,iy);
      y(bin_counter) = hw->GetBinContent(ix,iy)-hwMC->GetBinContent(ix,iy);
      inv_sqrtV(bin_counter,bin_counter) = 1./TMath::Sqrt(hwMC->GetBinContent(ix,iy));
      inv_V(bin_counter,bin_counter) = 1./hwMC->GetBinContent(ix,iy);
      bin_counter++;
    }
  }

  for(unsigned int j = 0; j<poi_counter; j++){
    //if( std::find(active_pois.begin(), active_pois.end(), j) == active_pois.end() ) continue;    
    unsigned int idx = active_pois[j];
    TH2D* hjac = fout->Get<TH2D>(Form("jac_%d", idx));
    bin_counter = 0;
    for(unsigned int ix = 1; ix<=nx; ix++ ){
      for(unsigned int iy = 1; iy<=ny; iy++ ){
	//jac(bin_counter,j+1) = N*hjac->GetBinContent(ix,iy); 
	jac(bin_counter,j) = N*hjac->GetBinContent(ix,iy); 
	bin_counter++;
      }
    }
    //cout << jac << endl;
  }

  fout->Close();
  
  MatrixXd A = inv_sqrtV*jac;
  MatrixXd b = inv_sqrtV*y;
  //cout << b << endl;
  VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
  MatrixXd C = (jac.transpose()*inv_V*jac).inverse();
  MatrixXd rho( C.rows(), C.rows() ) ;
  for(unsigned int ir = 0; ir<C.rows(); ir++){
    for(unsigned int ic = 0; ic<C.rows(); ic++){
      rho(ir,ic) = C(ir,ic)/TMath::Sqrt(C(ir,ir)*C(ic,ic));
    }
  }

  cout << rho << endl;

  VectorXd pulls(x.size());
  for(unsigned int ip = 0; ip<pulls.size(); ip++){
    pulls(ip) = x(ip) / TMath::Sqrt(C(ip,ip));
  }
  cout << pulls << endl;

  for(unsigned int j = 0; j<poi_counter; j++){
    unsigned int idx = active_pois[j];
    cout << "poi " << idx << ": " << poi_val[idx] << " +/- " << TMath::Sqrt(C(j,j)) << endl;
  }

  sw.Stop();
  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;
  return 1;
}
