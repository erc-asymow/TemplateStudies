#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2D.h"
#include "TF1.h"
#include "TF2.h"
#include <TStopwatch.h>
#include <iostream>
#include <boost/program_options.hpp>

#include <Eigen/Core>

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

  long n;
  double poi_val[NMAX];
  unsigned int poi_cat[NMAX];
  unsigned int poi_idx[NMAX];
  unsigned int poi_counter;
  TTree* outtree = fout->Get<TTree>("outtree");
  outtree->SetBranchAddress("events", &n);
  outtree->SetBranchAddress("poi_counter", &poi_counter);
  outtree->SetBranchAddress("poi_val", &poi_val);
  outtree->SetBranchAddress("poi_cat", &poi_cat);
  outtree->SetBranchAddress("poi_idx", &poi_idx);
  outtree->GetEntry(0);    

  TH2D* hw = fout->Get<TH2D>("w");
  TH2D* hwMC = fout->Get<TH2D>("wMC");
  int nx = hw->GetXaxis()->GetNbins(); 
  int ny = hw->GetYaxis()->GetNbins(); 

  MatrixXd jac((nx*ny), poi_counter+1);

  // FIX jacobian wrt N

  unsigned int counter = 0;
  for(unsigned int ix = 1; ix<=nx; ix++ ){
    for(unsigned int iy = 1; iy<=ny; iy++ ){
      jac(counter,0) = (hw->GetBinContent(ix,iy) - hwMC->GetBinContent(ix,iy)); 
      counter++;
    }
  }
    
  for(unsigned int j = 0; j<poi_counter; j++){
    TH2D* hjac = fout->Get<TH2D>(Form("jac_%d", j));
    //VectorXd jac_v(nx*ny);
    unsigned int counter = 0;
    for(unsigned int ix = 1; ix<=nx; ix++ ){
      for(unsigned int iy = 1; iy<=ny; iy++ ){
	jac(counter,j+1) = n*hjac->GetBinContent(ix,iy); 
	counter++;
      }
    }
    cout << jac << endl;
  }


  fout->Close();
  /*
  for(unsigned int p = 0 ; p<poi_counter; p++){
    cout << "poi_val[" << p << "]" << " = " << poi_val[p] << endl; 
  }
  */
  
  
  sw.Stop();
  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;
  return 1;
}
