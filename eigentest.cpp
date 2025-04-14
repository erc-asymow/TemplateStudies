#include <Eigen/Core>
#include <Eigen/Dense>
#include "TRandom3.h"
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
	("ndim",       value<int>()->default_value(100), "")
	("noise",       value<double>()->default_value(0.01), "")
	("ntoys",       value<int>()->default_value(10000), "");
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
  
  int ndim  = vm["ndim"].as<int>();
  int ntoys  = vm["ntoys"].as<int>();
  double noise  = vm["noise"].as<double>();

  VectorXd u0 = VectorXd::Zero(ndim);
  VectorXd u = VectorXd::Zero(ndim);
  MatrixXd nunuT = MatrixXd::Zero(ndim,ndim);
  MatrixXd identity = MatrixXd::Zero(ndim,ndim);
  for(unsigned int i=0; i<ndim; i++){
    u0(i) = ran->Gaus(0., 1.);
    identity(i,i) = 1;
  }
  u0 /= u0.norm();

  double nu2 = 0.;
  double nu_mean = 0.;
  for(unsigned int itoy=0; itoy < ntoys; itoy++){
    for(unsigned int i=0; i<ndim; i++){
      u(i) = u0(i) + ran->Gaus(0., noise);
    }
    u /= u.norm();
    VectorXd nu = u - u0;
    nu2 += nu.squaredNorm();
    nu_mean += nu.norm();
    nunuT += nu*nu.transpose();
  }
  nunuT /= ntoys;
  nu2 /= ntoys;
  nu_mean /= ntoys;

  nunuT /= (nu2);
  nunuT *= ndim-1;
  
  cout << "nu2 = " << nu2 << endl;
  cout << "nu_mean2 = " << nu_mean*nu_mean << endl;
  cout << "nu*nu.T: " << endl;
  cout << nunuT << endl;
  cout <<  "(1 - nu0*nu0.T)/<nu2>*(n-1): " << endl;
  cout << identity-(u0*u0.transpose()) << endl;
  cout << "Difference: " << endl;
  cout << nunuT - (identity-(u0*u0.transpose())) << endl;
  
  delete ran;
  return 0;
}
