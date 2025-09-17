#include <Eigen/Core>
#include <Eigen/Dense>
#include "TRandom3.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TFile.h"
#include <boost/program_options.hpp>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;
using namespace ROOT;

using namespace boost::program_options;

int main(int argc, char* argv[])
{

    variables_map vm;
  try
    {
      options_description desc{"Options"};
      desc.add_options()
	("help,h", "Help screen")
	("ndimX",       value<int>()->default_value(20), "")
	("ndimY",       value<int>()->default_value(20), "")
	("V",           value<double>()->default_value(100), "")
	("ntoys",       value<int>()->default_value(10), "");
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

  TRandom3* ran = new TRandom3();
  TFile* fout = TFile::Open("relaxation.root", "RECREATE");
  
  int ndimX  = vm["ndimX"].as<int>();
  int ndimY  = vm["ndimY"].as<int>();
  double V   = vm["V"].as<double>();
  int ntoys  = vm["ntoys"].as<int>();

  TH2D* h2 = new TH2D("grid","", ndimX, 0, ndimX, ndimY, 0, ndimY);
  h2->SetStats(0);
  h2->SetTitle(Form("V_{0}=%.0f V, %d iterations", V, ntoys));
  TCanvas* c = new TCanvas("c", "canvas", 800, int(800*float(ndimY)/ndimX));
  
  unsigned int iyG = (unsigned int)(ndimY*0.3);
  unsigned int iyV = (unsigned int)(ndimY*0.7);
  unsigned int ixL = (unsigned int)(ndimX*0.1);
  unsigned int ixH = (unsigned int)(ndimX*0.9);
  
  MatrixXd grid = MatrixXd::Zero(ndimX,ndimY);

  for(unsigned int iter=0; iter<ntoys; iter++){
    for(unsigned int iX=0; iX<ndimX; iX++){
      for(unsigned int iY=0; iY<ndimY; iY++){
	bool is_G  = 
	  (iX==0 || iX==ndimX-1) ||
	  (iY==0 || iY==ndimY-1) ||
	  (iY==iyG && iX>ixL && iX<=ixH );
	bool is_V  = 
	  (iY==iyV && iX>ixL && iX<=ixH );
	if(is_G)
	  grid(iX,iY) = 0.;
	else if(is_V)
	  grid(iX,iY) = V;
	else
	  grid(iX,iY) = 0.25*( grid(iX-1,iY) + grid(iX+1,iY) +
			       grid(iX,iY-1) + grid(iX,iY+1) );
	h2->SetBinContent(iX+1,iY+1, grid(iX,iY));
      }
    }
  }

  fout->cd();
  c->cd();
  h2->Draw("colz");
  h2->Write();
  c->Write();
  fout->Close();
  //cout << grid << endl;
  
  delete ran;
  return 0;
}
