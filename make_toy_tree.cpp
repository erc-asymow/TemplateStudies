#include "TFile.h"
#include "TTree.h"

using namespace std;
using namespace ROOT;

enum pdf_type { pdf_x=0, pdf_y, corr_x, corr_y, A0_x, A0_y, A1_x, A1_y, A4_x, A4_y, unknown};

auto degs = [](const pdf_type& pdf){
  switch(pdf){
  case pdf_type::pdf_x:  return  5; // f(x)
  case pdf_type::pdf_y:  return  5; // f(y|0)
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

  std::string tag = "";
  if(argc>1) tag = std::string(argv[1]);

  TFile* fout = TFile::Open(("input_"+tag+".root").c_str(), "RECREATE");
  TTree* tree = new TTree("tree", "tree");

  double pdf_x[NMAX];
  for(int i = 0; i<=degs(pdf_type::pdf_x); i++){
    tree->Branch(Form("pdfx_%d", i), &(pdf_x[i]), Form("pdfx_%d/D", i));
    pdf_x[i] = 1.0;
  }

  tree->Fill();
  tree->Write();
  fout->Close();

  return 0;
}
