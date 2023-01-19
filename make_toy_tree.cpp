#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
//#include <boost/program_options.hpp>

using namespace std;
using namespace ROOT;
//using namespace boost::program_options;

enum pdf_type { pdf_x=0, pdf_y, corr_x, corr_y, A0_x, A0_y, A1_x, A1_y, A4_x, A4_y, unknown};

auto degs = [](const pdf_type& pdf){
  switch(pdf){
  case pdf_type::pdf_x:  return  30; // f(x)
  case pdf_type::pdf_y:  return  20; // f(y|0)
  case pdf_type::corr_x: return  2; // P(x,.)
  case pdf_type::corr_y: return  2; // P(.,y)
  case pdf_type::A0_x:   return  2; // A0(x,.)
  case pdf_type::A0_y:   return  2; // A0(.,y)
  case pdf_type::A1_x:   return  2; // A1(x,.)
  case pdf_type::A1_y:   return  2; // A1(.,y)
  case pdf_type::A4_x:   return  2; // A4(x,.)
  case pdf_type::A4_y:   return  2; // A4(.,y)
  default: return 1;
  }
};

constexpr int NMAX = 100;

int main(int argc, char* argv[]){

  /*
  try
    {
      options_description desc{"Options"};
      desc.add_options()
	("help,h", "Help screen")
	("pi", value<float>()->default_value(3.14f), "Pi");
      variables_map vm;
      store(parse_command_line(argc, argv, desc), vm);
      notify(vm);
      if (vm.count("help"))
	std::cout << desc << '\n';
      else if (vm.count("pi"))
	std::cout << "Pi: " << vm["pi"].as<float>() << '\n';
    }
  catch (const error &ex)
    {
      std::cerr << ex.what() << '\n';
    }
  */

  std::string tag = "";
  if(argc>1) tag = std::string(argv[1]);

  TFile* fout = TFile::Open(("input_"+tag+".root").c_str(), "RECREATE");
  TTree* tree = new TTree("tree", "tree");

  const double max_x = 0.15;
  const double max_y = 2.5;

  TF1* toy_x = new TF1("toy_x", "[0]/(x-[1])", 0.0, max_x);  
  toy_x->SetParameter(0, 1.0);
  toy_x->SetParameter(1, -2.35e-03);
  double pdf_x[NMAX];
  for(int i = 0; i<=degs(pdf_type::pdf_x); i++){
    tree->Branch(Form("pdfx_%d", i), &(pdf_x[i]), Form("pdfx_%d/D", i));
    //pdf_x[i] = 1.0/max_x;
    pdf_x[i] = toy_x->Eval( (TMath::Cos(i*TMath::Pi()/degs(pdf_type::pdf_x))+1.0)*0.5*max_x );
  }

  TF1* toy_y = new TF1("toy_y", "1/TMath::Sqrt(2*TMath::Pi()*[1])*TMath::Exp(-0.5*(x-[0])*(x-[0])/[1]/[1])", -max_y, max_y);
  toy_y->SetParameter(0, 0.0);
  toy_y->SetParameter(1, 2.0);
  double pdf_y[NMAX];
  for(int j = 0; j<=degs(pdf_type::pdf_y); j++){
    tree->Branch(Form("pdfy_%d", j), &(pdf_y[j]), Form("pdfy_%d/D", j));
    //pdf_y[j] = 1.0/(2*max_y);
    pdf_y[j] = toy_y->Eval( (TMath::Cos(j*TMath::Pi()/degs(pdf_type::pdf_y))+1.0)*0.5*max_y );
  }

  double corr_xy[NMAX];
  for(int k = 0; k<=degs(pdf_type::corr_x); k++){
    for(int l = 0; l<=degs(pdf_type::corr_y); l++){      
      int idx = (degs(pdf_type::corr_y)+1)*k + l; 
      tree->Branch(Form("corrxy_%d_%d", k,l), &(corr_xy[idx]), Form("corrxy_%d_%d/D", k,l));
      corr_xy[idx] = 1.0;
    }
  }

  double A0_xy[NMAX];
  for(int m = 0; m<=degs(pdf_type::A0_x); m++){
    for(int n = 0; n<=degs(pdf_type::A0_y); n++){
      int idx = (degs(pdf_type::A0_y)+1)*m + n; 
      tree->Branch(Form("A0_xy_%d_%d", m,n), &(A0_xy[idx]), Form("A0_xy_%d_%d/D", m,n));
      A0_xy[idx] = 0.0;      
    }
  }

  tree->Fill();
  tree->Write();
  fout->Close();

  return 0;
}
