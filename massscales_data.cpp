#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include "TRandom3.h"
#include "TVector.h"
#include "TVectorT.h"
#include "TMath.h"
#include "TError.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include <TMatrixD.h>
#include <TMatrixDSymfwd.h>
#include <TStopwatch.h>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <boost/program_options.hpp>
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/FCNGradientBase.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include "RooRealVar.h"
#include "RooDerivative.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooCrystalBall.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooPlot.h"

using namespace RooFit;

//#include <Eigen/Core>
//#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;
using namespace ROOT;
using namespace ROOT::Minuit2;

typedef ROOT::VecOps::RVec<double> RVecD;
typedef ROOT::VecOps::RVec<unsigned int> RVecUI;
typedef ROOT::VecOps::RVec<int> RVecI;
typedef ROOT::VecOps::RVec<float> RVecF;
typedef ROOT::VecOps::RVec<bool> RVecB;
using ROOT::RDF::RNode;


using namespace boost::program_options;

constexpr double lumiMC2016 = 3.33369e+08/2001.9e+03;
constexpr double lumiMC2017 = 4.9803e+07/2001.9e+03;
constexpr double lumiMC2018 = 6.84093e+07/2001.9e+03;

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
	("minNumEvents",     value<int>()->default_value(100), "number of events")
	("minNumEventsPerBin",   value<int>()->default_value(10), "bias")
	("lumi",        value<float>()->default_value(16.1), "number of events")
	("tag",         value<std::string>()->default_value("closure"), "run type")
	("run",         value<std::string>()->default_value("closure"), "run type")
	("saveHistos",  bool_switch()->default_value(false), "")
	("firstIter",     value<int>()->default_value(-1), "firstIter")
	("lastIter",     value<int>()->default_value(-1), "lastIter")
	("nRMSforGausFit",     value<float>()->default_value(-1.), "number of events")
	("minNumMassBins",     value<int>()->default_value(4), "number of events")
	("rebin",       value<int>()->default_value(2), "rebin")
	("fitWidth",    bool_switch()->default_value(false), "")
	("fitNorm",     bool_switch()->default_value(false), "")
	("usePrevFit",  bool_switch()->default_value(false), "")	
	("tagPrevFit",         value<std::string>()->default_value("closure"), "run type")
	("runPrevFit",         value<std::string>()->default_value("closure"), "run type")
	("useSmearFit",         bool_switch()->default_value(false), "")
	("tagSmearFit",         value<std::string>()->default_value("closure"), "run type")
	("runSmearFit",         value<std::string>()->default_value("closure"), "run type")
	("useKf",             bool_switch()->default_value(false), "")
	("scaleToData",       bool_switch()->default_value(false), "")
	("useCBpdf",          bool_switch()->default_value(false), "")
	("minFrac",           value<float>()->default_value(0.995), "")
	("y2016",             bool_switch()->default_value(false), "")
	("y2017",             bool_switch()->default_value(false), "")
	("y2018",             bool_switch()->default_value(false), "")
	("seed",        value<int>()->default_value(4357), "seed");

      store(parse_command_line(argc, argv, desc), vm);
      notify(vm);
      if (vm.count("help")){
	std::cout << desc << '\n';
	return 0;
      }
      if (vm.count("tag"))        std::cout << "Tag: " << vm["tag"].as<std::string>() << '\n';
      if (vm.count("run"))        std::cout << "Run: " << vm["run"].as<std::string>() << '\n';
    }
  catch (const error &ex)
    {
      std::cerr << ex.what() << '\n';
    }
  
  int minNumEvents     = vm["minNumEvents"].as<int>();
  float lumi       = vm["lumi"].as<float>();
  float nRMSforGausFit = vm["nRMSforGausFit"].as<float>();
  std::string tag = vm["tag"].as<std::string>();
  std::string run = vm["run"].as<std::string>();
  int seed        = vm["seed"].as<int>();
  int minNumEventsPerBin   = vm["minNumEventsPerBin"].as<int>();
  int minNumMassBins = vm["minNumMassBins"].as<int>();
  int rebin       = vm["rebin"].as<int>();
  int firstIter   = vm["firstIter"].as<int>();
  int lastIter    = vm["lastIter"].as<int>();
  bool saveHistos     = vm["saveHistos"].as<bool>();
  bool fitWidth       = vm["fitWidth"].as<bool>();
  bool fitNorm        = vm["fitNorm"].as<bool>();
  bool usePrevFit        = vm["usePrevFit"].as<bool>();
  bool useSmearFit       = vm["useSmearFit"].as<bool>();
  bool useKf             = vm["useKf"].as<bool>();
  bool y2016             = vm["y2016"].as<bool>();
  bool y2017             = vm["y2017"].as<bool>();
  bool y2018             = vm["y2018"].as<bool>();
  std::string tagPrevFit = vm["tagPrevFit"].as<std::string>();
  std::string runPrevFit = vm["runPrevFit"].as<std::string>();
  std::string tagSmearFit = vm["tagSmearFit"].as<std::string>();
  std::string runSmearFit = vm["runSmearFit"].as<std::string>();
  bool scaleToData        = vm["scaleToData"].as<bool>();
  bool useCBpdf           = vm["useCBpdf"].as<bool>();
  float minFrac           = vm["minFrac"].as<float>();
    
  assert( y2016 || y2017 || y2018 );
  
  TRandom3* ran0 = new TRandom3(seed);

  vector<float> pt_edges  = {25, 30, 35, 40, 45, 55}; 
  vector<float> eta_edges = {-2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0,
			      0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};

  TH1F* h_pt_edges  = new TH1F("h_pt_edges", "",  pt_edges.size()-1, pt_edges.data());
  TH1F* h_eta_edges = new TH1F("h_eta_edges", "", eta_edges.size()-1, eta_edges.data());
  
  unsigned int n_pt_bins  = pt_edges.size()-1;
  unsigned int n_eta_bins = eta_edges.size()-1;
  int n_bins = n_pt_bins*n_pt_bins*n_eta_bins*n_eta_bins;

  TH1F* h_A_vals_nom = new TH1F("h_A_vals_nom", "", n_eta_bins, 0, n_eta_bins );
  TH1F* h_e_vals_nom = new TH1F("h_e_vals_nom", "", n_eta_bins, 0, n_eta_bins );
  TH1F* h_M_vals_nom = new TH1F("h_M_vals_nom", "", n_eta_bins, 0, n_eta_bins );
  TH1F* h_A_vals_prevfit = new TH1F("h_A_vals_prevfit", "", n_eta_bins, 0, n_eta_bins );
  TH1F* h_e_vals_prevfit = new TH1F("h_e_vals_prevfit", "", n_eta_bins, 0, n_eta_bins );
  TH1F* h_M_vals_prevfit = new TH1F("h_M_vals_prevfit", "", n_eta_bins, 0, n_eta_bins );
  TH1F* h_c_vals_prevfit = new TH1F("h_c_vals_prevfit", "", n_eta_bins, 0, n_eta_bins );
  TH1F* h_d_vals_prevfit = new TH1F("h_d_vals_prevfit", "", n_eta_bins, 0, n_eta_bins );
  
  float kmean_val = 0.5*( 1./pt_edges[0] + 1./pt_edges[ pt_edges.size()-1] );
  VectorXd A_vals_fit( n_eta_bins );
  VectorXd e_vals_fit( n_eta_bins );
  VectorXd M_vals_fit( n_eta_bins );
  VectorXd c_vals_fit( n_eta_bins );
  VectorXd d_vals_fit( n_eta_bins );

  // bias for A out
  for(unsigned int i=0; i<n_eta_bins; i++){
    h_A_vals_nom->SetBinContent(i+1, 0.0);
    h_A_vals_prevfit->SetBinContent(i+1, 0.0);
    h_c_vals_prevfit->SetBinContent(i+1, 0.0);
    h_d_vals_prevfit->SetBinContent(i+1, 0.0);
    A_vals_fit(i) = 0.0;
  }
  // bias for e out
  for(unsigned int i=0; i<n_eta_bins; i++){
    h_e_vals_nom->SetBinContent(i+1, 0.0);
    h_e_vals_prevfit->SetBinContent(i+1, 0.0);
    e_vals_fit(i) = 0.0;
  }
  // bias for M out
  for(unsigned int i=0; i<n_eta_bins; i++){
    h_M_vals_nom->SetBinContent(i+1, 0.0);
    h_M_vals_prevfit->SetBinContent(i+1, 0.0);
    M_vals_fit(i) = 0.0;
  }            
  
  std::vector<string> recos = {"reco", "smear0"}; 

  std::map<string, TH1D*> h_map;
  for(unsigned int r = 0; r<recos.size(); r++){
    h_map.insert( std::make_pair<string, TH1D* >("mean_"+recos[r], 0 ) );
    h_map.insert( std::make_pair<string, TH1D* >("rms_"+recos[r],  0 ) );
    h_map.insert( std::make_pair<string, TH1D* >("mask_"+recos[r],  0 ) );
  }

  std::map<string, TH2D*> hCB_map;
  for(unsigned int r = 0; r<recos.size(); r++){    
    hCB_map.insert( std::make_pair<string, TH2D* >("lnder_"+recos[r], 0 ) );
  }

  // map of positions in RVecF "masses"
  std::map<string, unsigned int> idx_map;
  idx_map.insert( std::make_pair<string, unsigned int >("reco",   1 ) );
  idx_map.insert( std::make_pair<string, unsigned int >("smear0", 2 ) );
  
  if(usePrevFit){
    TFile* ffit = TFile::Open(("./massfit_"+tagPrevFit+"_"+runPrevFit+".root").c_str(), "READ");
    if(ffit!=0){    
      cout << "Using fit results from " <<  std::string(ffit->GetName()) << " as new nominal for smear0" << endl;
      TH1D* h_A_vals_prevfit_in = (TH1D*)ffit->Get("h_A_vals_prevfit");
      TH1D* h_e_vals_prevfit_in = (TH1D*)ffit->Get("h_e_vals_prevfit");
      TH1D* h_M_vals_prevfit_in = (TH1D*)ffit->Get("h_M_vals_prevfit");
      for(unsigned int i=0; i<n_eta_bins; i++){
	A_vals_fit(i) = -h_A_vals_prevfit_in->GetBinContent(i+1);
	e_vals_fit(i) = -h_e_vals_prevfit_in->GetBinContent(i+1);
	M_vals_fit(i) = -h_M_vals_prevfit_in->GetBinContent(i+1);
      }
      h_A_vals_prevfit->Add(h_A_vals_prevfit_in, -1.0);
      h_e_vals_prevfit->Add(h_e_vals_prevfit_in, -1.0);
      h_M_vals_prevfit->Add(h_M_vals_prevfit_in, -1.0);
      ffit->Close();
    }
    else{
      cout << "No mass fit file!" << endl;
    }
  }

  if(useSmearFit){
    TFile* ffit = TFile::Open(("./resolfit_"+tagSmearFit+"_"+runSmearFit+".root").c_str(), "READ");
    if(ffit!=0){    
      cout << "Using fit results from " <<  std::string(ffit->GetName()) << " as MC smear" << endl;
      TH1D* h_c_vals_prevfit_in = (TH1D*)ffit->Get("h_c_vals_prevfit");
      TH1D* h_d_vals_prevfit_in = (TH1D*)ffit->Get("h_d_vals_prevfit");
      for(unsigned int i=0; i<n_eta_bins; i++){
	c_vals_fit(i) = h_c_vals_prevfit_in->GetBinContent(i+1);
	d_vals_fit(i) = h_d_vals_prevfit_in->GetBinContent(i+1);
      }
      h_c_vals_prevfit->Add(h_c_vals_prevfit_in, +1.0);
      h_d_vals_prevfit->Add(h_d_vals_prevfit_in, +1.0);
      ffit->Close();
    }
    else{
      cout << "No smear fit file!" << endl;
    }    
  }
    
  TFile* fout = TFile::Open(("./massscales_"+tag+"_"+run+".root").c_str(), firstIter<2 ? "RECREATE" : "UPDATE");
  
  for(int iter=-1; iter<3; iter++){

    if( !(iter>=firstIter && iter<=lastIter) ) continue;

    cout << "Doing iter " << iter << endl;
    //TTree* tree = new TTree("tree", "tree");

    vector<string> in_files = {};
    if(iter>=0){

      if(y2016){
	in_files = {
	  "/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_040854/0000/NanoV9MCPostVFP_*.root",
	  "/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_041233/0000/NanoV9MCPostVFP_*.root",
	  "/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_041233/0001/NanoV9MCPostVFP_*.root",
	  "/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_041233/0002/NanoV9MCPostVFP_*.root"
	};
      }
      else if(y2017){
	in_files = {
	  "/scratch/wmass/y2017/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2017_TrackFitV722_NanoProdv3/NanoV9MC2017_*.root"
	};
      }
      else if(y2018){
	in_files = {
	  "/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0000/NanoV9MC2018_*.root",
	  "/scratch/wmass/y2018/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MC2018_TrackFitV722_NanoProdv3/240124_121800/0001/NanoV9MC2018_*.root"
	};
      }      
    }
    else{
      if(y2016){
	in_files = {
	  "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016FDataPostVFP_TrackFitV722_NanoProdv6/240509_051502/0000/NanoV9DataPostVFP_*.root",
	  "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016GDataPostVFP_TrackFitV722_NanoProdv6/240509_051653/0000/NanoV9DataPostVFP_*.root",
	  "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016GDataPostVFP_TrackFitV722_NanoProdv6/240509_051653/0001/NanoV9DataPostVFP_*.root",
	  "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016GDataPostVFP_TrackFitV722_NanoProdv6/240509_051653/0002/NanoV9DataPostVFP_*.root",
	  "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016GDataPostVFP_TrackFitV722_NanoProdv6/240509_051653/0003/NanoV9DataPostVFP_*.root",
	  "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016GDataPostVFP_TrackFitV722_NanoProdv6/240509_051653/0004/NanoV9DataPostVFP_*.root",
	  "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV722_NanoProdv6/240509_051807/0000/NanoV9DataPostVFP_*.root",
	  "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV722_NanoProdv6/240509_051807/0001/NanoV9DataPostVFP_*.root",
	  "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV722_NanoProdv6/240509_051807/0002/NanoV9DataPostVFP_*.root",
	  "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV722_NanoProdv6/240509_051807/0003/NanoV9DataPostVFP_*.root",
	  "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV722_NanoProdv6/240509_051807/0004/NanoV9DataPostVFP_*.root",
	  "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV722_NanoProdv6/240509_051807/0005/NanoV9DataPostVFP_*.root",
	  "/scratch/wmass/y2016/SingleMuon/NanoV9Run2016HDataPostVFP_TrackFitV722_NanoProdv6/240509_051807/0006/NanoV9DataPostVFP_*.root"
	};
      }
      else if(y2017){
	in_files = {
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017B_TrackFitV722_NanoProdv3/240127_110915/0000/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017B_TrackFitV722_NanoProdv3/240127_110915/0001/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017B_TrackFitV722_NanoProdv3/240127_110915/0002/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017B_TrackFitV722_NanoProdv3/240127_110915/0003/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017C_TrackFitV722_NanoProdv3/240127_115941/0000/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017C_TrackFitV722_NanoProdv3/240127_115941/0001/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017C_TrackFitV722_NanoProdv3/240127_115941/0002/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017C_TrackFitV722_NanoProdv3/240127_115941/0003/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017C_TrackFitV722_NanoProdv3/240127_115941/0004/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017C_TrackFitV722_NanoProdv3/240127_115941/0005/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017D_TrackFitV722_NanoProdv3/240127_120137/0000/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017D_TrackFitV722_NanoProdv3/240127_120137/0001/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017D_TrackFitV722_NanoProdv3/240127_120137/0002/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017E_TrackFitV722_NanoProdv3/240127_121346/0000/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017E_TrackFitV722_NanoProdv3/240127_121346/0001/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017E_TrackFitV722_NanoProdv3/240127_121346/0002/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017E_TrackFitV722_NanoProdv3/240127_121346/0003/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017E_TrackFitV722_NanoProdv3/240127_121346/0004/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017F_TrackFitV722_NanoProdv3/240127_122701/0000/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017F_TrackFitV722_NanoProdv3/240127_122701/0001/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017F_TrackFitV722_NanoProdv3/240127_122701/0002/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017F_TrackFitV722_NanoProdv3/240127_122701/0003/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017F_TrackFitV722_NanoProdv3/240127_122701/0004/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017F_TrackFitV722_NanoProdv3/240127_122701/0005/NanoV9Data2017_*.root",
	  "/scratch/wmass/y2017/SingleMuon/NanoV9Run2017F_TrackFitV722_NanoProdv3/240127_122701/0006/NanoV9Data2017_*.root"
	};
      }
      else if(y2018){
	in_files = {
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018A_TrackFitV722_NanoProdv3/231102_185937/0000/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018A_TrackFitV722_NanoProdv3/231102_185937/0001/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018A_TrackFitV722_NanoProdv3/231102_185937/0002/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018A_TrackFitV722_NanoProdv3/231102_185937/0003/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018A_TrackFitV722_NanoProdv3/231102_185937/0004/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018A_TrackFitV722_NanoProdv3/231102_185937/0005/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018A_TrackFitV722_NanoProdv3/231102_185937/0006/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018B_TrackFitV722_NanoProdv3/231103_093816/0000/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018B_TrackFitV722_NanoProdv3/231103_093816/0001/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018B_TrackFitV722_NanoProdv3/231103_093816/0002/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018C_TrackFitV722_NanoProdv3/231103_101410/0000/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018C_TrackFitV722_NanoProdv3/231103_101410/0001/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018C_TrackFitV722_NanoProdv3/231103_101410/0002/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018D_TrackFitV722_NanoProdv3/231107_134901/0000/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018D_TrackFitV722_NanoProdv3/231107_134901/0001/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018D_TrackFitV722_NanoProdv3/231107_134901/0002/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018D_TrackFitV722_NanoProdv3/231107_134901/0003/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018D_TrackFitV722_NanoProdv3/231107_134901/0004/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018D_TrackFitV722_NanoProdv3/231107_134901/0005/NanoV9Data2018_*.root",
	  "/scratch/wmass/y2018/SingleMuon/NanoV9Run2018D_TrackFitV722_NanoProdv3/231107_134901/0006/NanoV9Data2018_*.root"
	};
      }
    }
    
    ROOT::RDataFrame d( "Events", in_files );
    unsigned int nslots = d.GetNSlots();
    std::vector<TRandom3*> rans = {};
    for(unsigned int i = 0; i < nslots; i++){
      rans.emplace_back( new TRandom3(seed + i*10 + iter) );
    }

    auto dlast = std::make_unique<RNode>(d);
        
    // MC
    if(iter>=0){

      dlast = std::make_unique<RNode>(dlast->Define("idxs", [&](UInt_t nMuon, RVecB Muon_looseId, RVecF Muon_dxybs, RVecB Muon_isGlobal,
								RVecB Muon_highPurity, RVecB Muon_mediumId, RVecF Muon_pfRelIso04_all,
								RVecF Muon_pt, RVecF Muon_eta)->RVecUI {
	RVecUI out;
	for(unsigned int i = 0; i < nMuon; i++){
	  if( Muon_looseId[i] && TMath::Abs(Muon_dxybs[i]) < 0.05 && Muon_isGlobal[i] && Muon_highPurity[i] && Muon_mediumId[i] && Muon_pfRelIso04_all[i]<0.15 &&
	      Muon_pt[i] >= pt_edges[0] && Muon_pt[i] < pt_edges[ n_pt_bins ]  && Muon_eta[i]>=eta_edges[0] && Muon_eta[i]<=eta_edges[ n_eta_bins ]  ){
	    out.emplace_back(i);
	  }
	}
	return out;
      }, {"nMuon", "Muon_looseId", "Muon_dxybs", "Muon_isGlobal",
	  "Muon_highPurity","Muon_mediumId", "Muon_pfRelIso04_all",
	  useKf ? "Muon_pt" : "Muon_cvhidealPt", useKf ? "Muon_eta" : "Muon_cvhidealEta" } ));

      dlast = std::make_unique<RNode>(dlast->Filter( [](RVecUI idxs, RVecI Muon_charge, bool HLT_IsoMu24 ){
	if( idxs.size()!=2 || !HLT_IsoMu24) return false;
	if( Muon_charge[idxs[0]]*Muon_charge[idxs[1]] > 0 ) return false;
	return true;
      }, {"idxs", "Muon_charge", "HLT_IsoMu24"} ));
      
      dlast = std::make_unique<RNode>(dlast->Define("weight", [](float weight)->float{
	return std::copysign(1.0, weight);
      }, {"Generator_weight"} ));          
      
      dlast = std::make_unique<RNode>(dlast->DefineSlot("Muon_ksmear", [&](unsigned int nslot, RVecUI idxs,
									   RVecF Muon_pt, RVecF Muon_eta, RVecF Muon_phi, RVecF Muon_mass, RVecI Muon_charge,
									   UInt_t nGenPart, RVecI GenPart_status, RVecI GenPart_statusFlags, RVecI GenPart_pdgId,
									   RVecF GenPart_pt, RVecF GenPart_eta, RVecF GenPart_phi, RVecF GenPart_mass)->RVecF {
	RVecF out;
	unsigned int idxP = Muon_charge[idxs[0]]>0 ? idxs[0] : idxs[1];
	unsigned int idxM = Muon_charge[idxs[0]]>0 ? idxs[1] : idxs[0];
	ROOT::Math::PtEtaPhiMVector muP( Muon_pt[ idxP ], Muon_eta[ idxP ], Muon_phi[ idxP ], Muon_mass[ idxP ] );
	ROOT::Math::PtEtaPhiMVector muM( Muon_pt[ idxM ], Muon_eta[ idxM ], Muon_phi[ idxM ], Muon_mass[ idxM ] );
	ROOT::Math::PtEtaPhiMVector gmuP;
	ROOT::Math::PtEtaPhiMVector gmuM;
	for(unsigned int i = 0; i < nGenPart; i++){
	  bool isGoodGenPart = (GenPart_status[i]==1 && (GenPart_statusFlags[i] & 1 || (GenPart_statusFlags[i] & (1<<5))) && TMath::Abs(GenPart_pdgId[i])==13);
	  if(!isGoodGenPart) continue;
	  ROOT::Math::PtEtaPhiMVector gen(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
	  if( ROOT::Math::VectorUtil::DeltaR(gen, muP) < 0.1 && ROOT::Math::VectorUtil::DeltaR(gen, muM) > 0.1) gmuP = gen;
	  else if( ROOT::Math::VectorUtil::DeltaR(gen, muP) > 0.1 && ROOT::Math::VectorUtil::DeltaR(gen, muM) < 0.1) gmuM = gen;
	}
	
	if( gmuP.Pt()>10. && gmuM.Pt()>10.){
	  
	  float kmuP = 1./muP.Pt();
	  float kmuM = 1./muM.Pt();
	  float kgmuP = 1./gmuP.Pt();
	  float kgmuM = 1./gmuM.Pt();
	  
	  float scale_smear0P = 1.0;
	  float scale_smear0M = 1.0;
	  float extraSmear0P = 0.0;
	  float extraSmear0M = 0.0;
	  
	  unsigned int ietaP = n_eta_bins;
	  for(unsigned int ieta_p = 0; ieta_p<n_eta_bins; ieta_p++){
	    float eta_p_low = eta_edges[ieta_p];
	    float eta_p_up  = eta_edges[ieta_p+1];
	    if( muP.Eta()>=eta_p_low && muP.Eta()<eta_p_up) ietaP = ieta_p;	  
	  }
	  unsigned int ietaM = n_eta_bins;
	  for(unsigned int ieta_m = 0; ieta_m<n_eta_bins; ieta_m++){
	    float eta_m_low = eta_edges[ieta_m];
	    float eta_m_up  = eta_edges[ieta_m+1];
	    if( muM.Eta()>=eta_m_low && muM.Eta()<eta_m_up) ietaM = ieta_m;	  
	  }		  
	  if(ietaP<n_eta_bins && ietaM<n_eta_bins){
	    scale_smear0P = (1. + A_vals_fit(ietaP) + e_vals_fit(ietaP)*kmuP - M_vals_fit(ietaP)/kmuP);
	    scale_smear0M = (1. + A_vals_fit(ietaM) + e_vals_fit(ietaM)*kmuM + M_vals_fit(ietaM)/kmuM);
	    //cout << "smear0:" << scale_smear0P << ": " << 1 << " + " << A_vals_fit(ietaP) << " + " << e_vals_fit(ietaP)*kmuP << " - " << M_vals_fit(ietaP)/kmuP << endl;
	    //cout << "smear1:" << scale_smear1P << ": " << 1 << " + " << A_vals_nom(ietaP) << " + " << e_vals_nom(ietaP)*kmuP << " - " << M_vals_nom(ietaP)/kmuP << endl;
	    if(useSmearFit){
	      extraSmear0P = TMath::Sqrt( TMath::Max( 1.0 + c_vals_fit(ietaP) + d_vals_fit(ietaP)*kmuP, 0.0)  ) - 1.0;
	      extraSmear0M = TMath::Sqrt( TMath::Max( 1.0 + c_vals_fit(ietaM) + d_vals_fit(ietaM)*kmuM, 0.0)  ) - 1.0;
	      //cout << "extraSmear0P: sqrt( max(1.0 + " << c_vals_fit(ietaP)  << " + " << d_vals_fit(ietaP)*kmuP << ")) - 1.0 = " << extraSmear0P << endl;
	      //cout << "extraSmear0M: sqrt( max(1.0 + " << c_vals_fit(ietaM)  << " + " << d_vals_fit(ietaM)*kmuM << ")) - 1.0 = " << extraSmear0M << endl;  
	    }
	  }	  
	  float kmuPsmear0 = (kgmuP + (kmuP - kgmuP)*(1.0 + extraSmear0P))*scale_smear0P;
	  float kmuMsmear0 = (kgmuM + (kmuM - kgmuM)*(1.0 + extraSmear0M))*scale_smear0M;	  
	  out.emplace_back( kmuPsmear0 );
	  out.emplace_back( kmuMsmear0 );
	}
	else{
	  out.emplace_back(0.0);
	  out.emplace_back(0.0);
	}
	return out;
      }, {"idxs",
	  useKf ? "Muon_pt" : "Muon_cvhidealPt", useKf ? "Muon_eta" : "Muon_cvhidealEta", useKf ? "Muon_phi" : "Muon_cvhidealPhi", "Muon_mass", "Muon_charge",
	  "nGenPart", "GenPart_status", "GenPart_statusFlags", "GenPart_pdgId",
	  "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass"} ));
      
      dlast = std::make_unique<RNode>(dlast->Define("indexes", [&](RVecUI idxs, RVecF Muon_pt, RVecF Muon_eta, RVecI Muon_charge, RVecF Muon_ksmear)-> RVecUI {
	unsigned int idxP = Muon_charge[idxs[0]]>0 ? idxs[0] : idxs[1];
	unsigned int idxM = Muon_charge[idxs[0]]>0 ? idxs[1] : idxs[0];
	float ptP  = Muon_pt[idxP];
	float ptM  = Muon_pt[idxM];
	float ksmear0P = Muon_ksmear[0]>0. ? Muon_ksmear[0] : 1./(pt_edges[0]-0.01);
	float ksmear0M = Muon_ksmear[1]>0. ? Muon_ksmear[1] : 1./(pt_edges[0]-0.01);
	float etaP = Muon_eta[idxP];
	float etaM = Muon_eta[idxM];
	RVecUI out;
	out.emplace_back(n_bins);
	out.emplace_back(n_bins);
	
	unsigned int ibin = 0;
	for(unsigned int ieta_p = 0; ieta_p<n_eta_bins; ieta_p++){
	  float eta_p_low = eta_edges[ieta_p];
	  float eta_p_up  = eta_edges[ieta_p+1];      
	  for(unsigned int ipt_p = 0; ipt_p<n_pt_bins; ipt_p++){
	    float pt_p_low = pt_edges[ipt_p];
	    float pt_p_up  = pt_edges[ipt_p+1];
	    for(unsigned int ieta_m = 0; ieta_m<n_eta_bins; ieta_m++){
	      float eta_m_low = eta_edges[ieta_m];
	      float eta_m_up  = eta_edges[ieta_m+1];      
	      for(unsigned int ipt_m = 0; ipt_m<n_pt_bins; ipt_m++){
		float pt_m_low = pt_edges[ipt_m];
		float pt_m_up  = pt_edges[ipt_m+1];
	      if( etaP>=eta_p_low && etaP<eta_p_up &&
		  etaM>=eta_m_low && etaM<eta_m_up &&
		  ptP>=pt_p_low   && ptP<pt_p_up &&
		  ptM>=pt_m_low   && ptM<pt_m_up 
		  ) out[0] = ibin;
	      if( etaP>=eta_p_low && etaP<eta_p_up &&
		  etaM>=eta_m_low && etaM<eta_m_up &&
		  1./ksmear0P>=pt_p_low   && 1./ksmear0P<pt_p_up &&
		  1./ksmear0M>=pt_m_low   && 1./ksmear0M<pt_m_up 
		  ) out[1] = ibin;
	      ibin++;
	      }
	    }
	  }
	}
	return out;
      }, {"idxs", useKf ? "Muon_pt" : "Muon_cvhidealPt", useKf ? "Muon_eta" : "Muon_cvhidealEta", "Muon_charge", "Muon_ksmear"} ));
      
      for(unsigned int r = 0 ; r<recos.size(); r++){
	dlast = std::make_unique<RNode>(dlast->Define( TString(("index_"+recos[r]).c_str()), [r](RVecUI indexes){
	  return indexes.at(r);
	}, {"indexes"} ));
      }
      
      dlast = std::make_unique<RNode>(dlast->Define("masses", [&](RVecUI idxs,
								  RVecF Muon_pt, RVecF Muon_eta, RVecF Muon_phi, RVecF Muon_mass, RVecI Muon_charge,
								  UInt_t nGenPart, RVecI GenPart_status, RVecI GenPart_statusFlags, RVecI GenPart_pdgId,
								  RVecF GenPart_pt, RVecF GenPart_eta, RVecF GenPart_phi, RVecF GenPart_mass,
								  RVecF Muon_ksmear)->RVecF {
	RVecF out;
	unsigned int idxP = Muon_charge[idxs[0]]>0 ? idxs[0] : idxs[1];
	unsigned int idxM = Muon_charge[idxs[0]]>0 ? idxs[1] : idxs[0];
	ROOT::Math::PtEtaPhiMVector muP( Muon_pt[ idxP ], Muon_eta[ idxP ], Muon_phi[ idxP ], Muon_mass[ idxP ] );
	ROOT::Math::PtEtaPhiMVector muM( Muon_pt[ idxM ], Muon_eta[ idxM ], Muon_phi[ idxM ], Muon_mass[ idxM ] );
	ROOT::Math::PtEtaPhiMVector gmuP;
	ROOT::Math::PtEtaPhiMVector gmuM;
	for(unsigned int i = 0; i < nGenPart; i++){
	  bool isGoodGenPart = (GenPart_status[i]==1 && (GenPart_statusFlags[i] & 1 || (GenPart_statusFlags[i] & (1<<5))) && TMath::Abs(GenPart_pdgId[i])==13);
	  if(!isGoodGenPart) continue;
	  ROOT::Math::PtEtaPhiMVector gen(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
	  if( ROOT::Math::VectorUtil::DeltaR(gen, muP) < 0.1 && ROOT::Math::VectorUtil::DeltaR(gen, muM) > 0.1) gmuP = gen;
	  else if( ROOT::Math::VectorUtil::DeltaR(gen, muP) > 0.1 && ROOT::Math::VectorUtil::DeltaR(gen, muM) < 0.1) gmuM = gen;
	}
	
	if( gmuP.Pt()>10. && gmuM.Pt()>10.){
	  out.emplace_back( (gmuP + gmuM).M() );
	  out.emplace_back( (muP + muM).M() );	  
	  float ksmear0P = Muon_ksmear[0]>0. ? Muon_ksmear[0] : 1./(pt_edges[0]-0.01);
	  float ksmear0M = Muon_ksmear[1]>0. ? Muon_ksmear[1] : 1./(pt_edges[0]-0.01);	  
	  ROOT::Math::PtEtaPhiMVector muP_smear0( 1./ksmear0P,
						  Muon_eta[ idxP ], Muon_phi[ idxP ], Muon_mass[ idxP ] );
	  ROOT::Math::PtEtaPhiMVector muM_smear0( 1./ksmear0M,
						  Muon_eta[ idxM ], Muon_phi[ idxM ], Muon_mass[ idxM ] );      
	  out.emplace_back( (muP_smear0 + muM_smear0).M() );	
	}
	
	return out;
      }, {"idxs",
	  useKf ? "Muon_pt" : "Muon_cvhidealPt", useKf ? "Muon_eta" : "Muon_cvhidealEta", useKf ? "Muon_phi" : "Muon_cvhidealPhi", "Muon_mass", "Muon_charge",
	  "nGenPart", "GenPart_status", "GenPart_statusFlags", "GenPart_pdgId",
	  "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass",
	  "Muon_ksmear"} ));
      
      for(unsigned int r = 0 ; r<recos.size(); r++){
	unsigned int mpos = idx_map.at(recos[r]);
	dlast = std::make_unique<RNode>(dlast->Define(TString( (recos[r]+"_m").c_str() ), [mpos](RVecF masses){
	  return masses.size()>0 ? masses.at( mpos ) : -99.;
	}, {"masses"} ));      
	dlast = std::make_unique<RNode>(dlast->Define(TString( (recos[r]+"_dm").c_str() ), [mpos](RVecF masses){
	  return masses.size()>0 ? masses.at( mpos ) - masses.at(0) : -99.;
	}, {"masses"} ));            
      }

      for(unsigned int r = 0 ; r<recos.size(); r++){
	string rname = recos[r];
	unsigned int rpos = idx_map.at(rname);

	dlast = std::make_unique<RNode>(dlast->Define(TString( (recos[r]+"_weights_jac").c_str() ), [n_bins,r,rpos,rname,&h_map,&hCB_map,useCBpdf](RVecF masses, RVecUI indexes)->RVecF{
	  
	  RVecF out;
	  if(masses.size()==0){
	    out.emplace_back(0.0);
	    out.emplace_back(0.0);
	    return out;
	  }
	
	  float gen_m  = masses.at(0);

	  if(!useCBpdf){
	    TH1D* h_mean = h_map.at("mean_"+rname);
	    TH1D* h_rms  = h_map.at("rms_"+rname);
	    float reco_m = masses.at( rpos );	
	    float reco_delta = 0.;
	    float reco_sigma = 0.;
	    if(indexes[r]<n_bins){
	      reco_delta = h_mean->GetBinContent(indexes[r]+1);
	      reco_sigma = h_rms->GetBinContent(indexes[r]+1);
	    }
	    float reco_jscale = reco_sigma>0. ? +(reco_m - (gen_m+reco_delta) )*(gen_m+reco_delta)/reco_sigma/reco_sigma : 0.0;
	    float reco_jwidth = reco_sigma>0. ? +(reco_m - (gen_m+reco_delta) )*(reco_m - (gen_m+reco_delta) )/reco_sigma/reco_sigma - 1.0 : 0.0;
	    out.emplace_back(reco_jscale);
	    out.emplace_back(reco_jwidth);
	  }
	  else{
	    TH2D* h_lnder = hCB_map.at("lnder_"+rname);
	    float reco_m = masses.at( rpos );
	    float dm = reco_m-gen_m;
	    int bin_lnder = h_lnder->GetYaxis()->FindBin( dm );
	    if( bin_lnder==0 )
	      bin_lnder = 1;
	    else if( bin_lnder==h_lnder->GetYaxis()->GetNbins()+1 )
	      bin_lnder = h_lnder->GetYaxis()->GetNbins();
	    float lnder = indexes[r]<n_bins ? h_lnder->GetBinContent(indexes[r]+1, bin_lnder) : 0.;
	    float reco_jscale = -lnder*reco_m;
	    float reco_jwidth = -(1.0 + lnder*dm);
	    out.emplace_back(reco_jscale);
	    out.emplace_back(reco_jwidth);
	  }
	  
	  return out;
	}, {"masses", "indexes"} ));
      }
      
      for(unsigned int r = 0 ; r<recos.size(); r++){
	dlast = std::make_unique<RNode>(dlast->Define( TString((recos[r]+"_jscale_weight").c_str()), [](RVecF weights_jac, float weight)->float{
	  return weights_jac.at( 0 )*weight;
	}, { recos[r]+"_weights_jac", "weight" } ));
	dlast = std::make_unique<RNode>(dlast->Define( TString((recos[r]+"_jwidth_weight").c_str()), [](RVecF weights_jac, float weight)->float{
	  return weights_jac.at( 1 )*weight;
	}, { recos[r]+"_weights_jac", "weight" } ));
      }
      
    }
    
    // data
    else{

      dlast = std::make_unique<RNode>(dlast->Define("idxs", [&](UInt_t nMuon, RVecB Muon_looseId, RVecF Muon_dxybs, RVecB Muon_isGlobal,
								RVecB Muon_highPurity, RVecB Muon_mediumId, RVecF Muon_pfRelIso04_all,
								RVecF Muon_pt, RVecF Muon_eta)->RVecUI {
	RVecUI out;
	for(unsigned int i = 0; i < nMuon; i++){
	  if( Muon_looseId[i] && TMath::Abs(Muon_dxybs[i]) < 0.05 && Muon_isGlobal[i] && Muon_highPurity[i] && Muon_mediumId[i] && Muon_pfRelIso04_all[i]<0.15 &&
	      Muon_pt[i] >= pt_edges[0] && Muon_pt[i] < pt_edges[ n_pt_bins ]  && Muon_eta[i]>=eta_edges[0] && Muon_eta[i]<=eta_edges[ n_eta_bins ]  ){
	    out.emplace_back(i);
	  }
	}
	return out;
      }, {"nMuon", "Muon_looseId", "Muon_dxybs", "Muon_isGlobal",
	  "Muon_highPurity","Muon_mediumId", "Muon_pfRelIso04_all",
	  useKf ? "Muon_pt" : "Muon_cvhPt", useKf ? "Muon_eta" : "Muon_cvhEta" } ));
      
      dlast = std::make_unique<RNode>(dlast->Filter( [](RVecUI idxs, RVecI Muon_charge, bool HLT_IsoMu24 ){
	if( idxs.size()!=2 || !HLT_IsoMu24) return false;
	if( Muon_charge[idxs[0]]*Muon_charge[idxs[1]] > 0 ) return false;
	return true;
      }, {"idxs", "Muon_charge", "HLT_IsoMu24"} ));      
	  
      dlast = std::make_unique<RNode>(dlast->Define("weight", []()->float{ return 1.0; }, {} ));          
            
      dlast = std::make_unique<RNode>(dlast->Define("index_data", [&](RVecUI idxs, RVecF Muon_pt, RVecF Muon_eta, RVecI Muon_charge)-> unsigned int {
	unsigned int idxP = Muon_charge[idxs[0]]>0 ? idxs[0] : idxs[1];
	unsigned int idxM = Muon_charge[idxs[0]]>0 ? idxs[1] : idxs[0];
	float ptP  = Muon_pt[idxP];
	float ptM  = Muon_pt[idxM];
	float etaP = Muon_eta[idxP];
	float etaM = Muon_eta[idxM];
	unsigned int out = n_bins;	
	unsigned int ibin = 0;
	for(unsigned int ieta_p = 0; ieta_p<n_eta_bins; ieta_p++){
	  float eta_p_low = eta_edges[ieta_p];
	  float eta_p_up  = eta_edges[ieta_p+1];      
	  for(unsigned int ipt_p = 0; ipt_p<n_pt_bins; ipt_p++){
	    float pt_p_low = pt_edges[ipt_p];
	    float pt_p_up  = pt_edges[ipt_p+1];
	    for(unsigned int ieta_m = 0; ieta_m<n_eta_bins; ieta_m++){
	      float eta_m_low = eta_edges[ieta_m];
	      float eta_m_up  = eta_edges[ieta_m+1];      
	      for(unsigned int ipt_m = 0; ipt_m<n_pt_bins; ipt_m++){
		float pt_m_low = pt_edges[ipt_m];
		float pt_m_up  = pt_edges[ipt_m+1];
	      if( etaP>=eta_p_low && etaP<eta_p_up &&
		  etaM>=eta_m_low && etaM<eta_m_up &&
		  ptP>=pt_p_low   && ptP<pt_p_up &&
		  ptM>=pt_m_low   && ptM<pt_m_up 
		  ) out = ibin;
	      ibin++;
	      }
	    }
	  }
	}
	return out;
      }, {"idxs", useKf ? "Muon_pt" : "Muon_cvhPt", useKf ? "Muon_eta" : "Muon_cvhEta", "Muon_charge"} ));
            
      dlast = std::make_unique<RNode>(dlast->Define("data_m", [&](RVecUI idxs,
								  RVecF Muon_pt, RVecF Muon_eta, RVecF Muon_phi, RVecF Muon_mass, RVecI Muon_charge)->float {
	float out = 0.0;
	unsigned int idxP = Muon_charge[idxs[0]]>0 ? idxs[0] : idxs[1];
	unsigned int idxM = Muon_charge[idxs[0]]>0 ? idxs[1] : idxs[0];
	ROOT::Math::PtEtaPhiMVector muP( Muon_pt[ idxP ], Muon_eta[ idxP ], Muon_phi[ idxP ], Muon_mass[ idxP ] );
	ROOT::Math::PtEtaPhiMVector muM( Muon_pt[ idxM ], Muon_eta[ idxM ], Muon_phi[ idxM ], Muon_mass[ idxM ] );
	out = (muP + muM).M();	  
	return out;
      }, {"idxs",
	  useKf ? "Muon_pt" : "Muon_cvhPt", useKf ? "Muon_eta" : "Muon_cvhEta", useKf ? "Muon_phi" : "Muon_cvhPhi", "Muon_mass", "Muon_charge"} ));           
    }
    
    std::vector<ROOT::RDF::RResultPtr<TH1D> > histos1D;
    std::vector<ROOT::RDF::RResultPtr<TH2D> > histos2D;
    
    const int x_nbins   = 40;
    const double x_low  = 70.0;
    const double x_high = 110.0;
  
    //histos1D.emplace_back(dlast->Histo1D({"h_gen_m", "nominal", x_nbins, x_low, x_high}, "gen_m", "weight"));      
    //histos1D.emplace_back(dlast->Histo1D({"h_reco_m", "nominal", x_nbins, x_low, x_high}, "reco_m", "weight"));
    //histos1D.emplace_back(dlast->Histo1D({"h_smear_m", "nominal", x_nbins, x_low, x_high}, "smear_m", "weight"));

    if(iter==-1){
      histos2D.emplace_back(dlast->Histo2D({ "h_data_bin_m", "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index_data", "data_m", "weight" ));
      auto colNames = dlast->GetColumnNames();
      double total = *(dlast->Count());  
      std::cout << colNames.size() << " columns created. Total event count is " << total  << std::endl;
    }
    else if(iter==0){
      for(unsigned int r = 0 ; r<recos.size(); r++){
	histos2D.emplace_back(dlast->Histo2D({ "h_"+TString(recos[r].c_str())+"_bin_m",    "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_m", "weight" ));
	histos2D.emplace_back(dlast->Histo2D({ "h_"+TString(recos[r].c_str())+"_bin_dm",   "nominal", n_bins, 0, double(n_bins), 40, -10.0, 10.0},          "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_dm", "weight"));
      }
      auto colNames = dlast->GetColumnNames();
      double total = *(dlast->Count());  
      std::cout << colNames.size() << " columns created. Total event count is " << total  << std::endl;
    }
    else if(iter==1){
      for(unsigned int r = 0 ; r<recos.size(); r++){
	if(recos[r]!="smear0") continue;
	histos2D.emplace_back(dlast->Histo2D({"h_"+TString(recos[r].c_str())+"_bin_jac_scale", "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_m", TString(recos[r].c_str())+"_jscale_weight"));
	histos2D.emplace_back(dlast->Histo2D({"h_"+TString(recos[r].c_str())+"_bin_jac_width", "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_m", TString(recos[r].c_str())+"_jwidth_weight"));
      }
    }

    if(iter<2){
      fout->cd();
      std::cout << "Writing histos..." << std::endl;
      double lumiMC = lumiMC2016;
      if(y2017)      lumiMC = lumiMC2017;
      else if(y2018) lumiMC = lumiMC2018;
      double sf = lumi>0. ? lumi/lumiMC : 1.0; //double(lumi)/double(minNumEvents);
      for(auto h : histos1D){
	if(iter>=0) h->Scale(sf);
	string h_name = std::string(h->GetName());
	h->Write();
      }
      for(auto h : histos2D){
	if(iter>=0) h->Scale(sf);
	string h_name = std::string(h->GetName());
	std::cout << "Total number of events in 2D histo " << h_name << ": " << h->GetEntries() << std::endl;
	h->Write();
      }    
      std::cout << "Total slots: " << dlast->GetNSlots() << std::endl;
    }
    
    // fill histos
    if(iter==0){

      cout << "Writing aux files" << endl;
      h_pt_edges->Write();
      h_eta_edges->Write();
      h_A_vals_nom->Write();
      h_e_vals_nom->Write();
      h_M_vals_nom->Write();
      h_A_vals_prevfit->Write();
      h_e_vals_prevfit->Write();
      h_M_vals_prevfit->Write();
      h_c_vals_prevfit->Write();
      h_d_vals_prevfit->Write();
      
      for(unsigned int r = 0 ; r<recos.size(); r++){

	if(recos[r]!="smear0") continue;
	
	TH2D* h_reco_dm = (TH2D*)fout->Get(TString( ("h_"+recos[r]+"_bin_dm").c_str()) );
	TH2D* h_reco_m  = (TH2D*)fout->Get(TString( ("h_"+recos[r]+"_bin_m").c_str()) );

	if( h_reco_dm==0 || h_reco_m==0 ){
	  cout << "h_reco_dm/h_reco_m NOT FOUND" << endl;
	  continue;
	}

	if(useCBpdf){
	  hCB_map["lnder_"+recos[r]] = new TH2D( TString( ("h_lnder_"+recos[r]+"_bin_dm").c_str() ),"", n_bins, 0, double(n_bins),
						 h_reco_dm->GetYaxis()->GetNbins()*2, h_reco_dm->GetYaxis()->GetXmin(), h_reco_dm->GetYaxis()->GetXmax());
	}
	else{
	  h_map["mean_"+recos[r]] = new TH1D( TString( ("h_mean_"+recos[r]+"_bin_dm").c_str() ),"", n_bins, 0, double(n_bins));
	  h_map["rms_"+recos[r]]  = new TH1D( TString( ("h_rms_"+recos[r]+"_bin_dm").c_str() ), "", n_bins, 0, double(n_bins));
	}
	h_map["mask_"+recos[r]] = new TH1D( TString( ("h_mask_"+recos[r]+"_bin_dm").c_str() ),"", n_bins, 0, double(n_bins));

	for(unsigned int i = 0; i<n_bins; i++ ){
	  if(i%1000==0){
	    if(!useCBpdf)
	      cout << "Doing gaus fit for bin " << i << " / " << n_bins << " of " << recos[r] << endl;
	    else
	      cout << "Doing CB fit for bin " << i << " / " << n_bins << " of " << recos[r] << endl;
	  }
	  TString projname(Form("bin_%d_", i));
	  projname += TString( recos[r].c_str() );
	  TH1D* hi   = (TH1D*)h_reco_dm->ProjectionY( projname+"_dm", i+1, i+1 );
	  TH1D* hi_m = (TH1D*)h_reco_m->ProjectionY( projname+"_m", i+1, i+1 );

	  // gaus pdf
	  if(!useCBpdf){
	    double mean_i = 0.0;
	    double meanerr_i = 0.0;
	    double rms_i = 0.0;
	    double rmserr_i = 0.0;
	    //cout << hi_m->Integral() << ", " << hi->Integral() << ", " << hi_m->GetMean() << endl;
	    if( hi_m->Integral() > minNumEvents && hi->Integral() > minNumEvents  &&  hi_m->GetMean()>75. && hi_m->GetMean()<105. ){
	      h_map.at("mask_"+recos[r])->SetBinContent(i+1, 1);
	      TF1* gf = new TF1("gf","[0]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp( -0.5*(x-[1])*(x-[1])/[2]/[2] )",
				hi->GetXaxis()->GetBinLowEdge(1), hi->GetXaxis()->GetBinUpEdge( hi->GetXaxis()->GetNbins() ));      
	      gf->SetParameter(0, hi->Integral());
	      gf->SetParameter(1, hi->GetMean());
	      gf->SetParameter(2, hi->GetRMS() );
	      float m_min = nRMSforGausFit>0. ? TMath::Max(-nRMSforGausFit*hi->GetRMS(), -6.0) : -6.0;
	      float m_max = nRMSforGausFit>0. ? TMath::Min(+nRMSforGausFit*hi->GetRMS(), +6.0) : +6.0;
	      hi->Fit("gf", "QR", "", m_min, m_max );
	      mean_i    = gf->GetParameter(1);
	      meanerr_i = gf->GetParError(1);
	      rms_i     = TMath::Abs(gf->GetParameter(2));
	      rmserr_i  = gf->GetParError(2);
	      //cout << "Fit " << mean_i << endl;
	      delete gf;
	    }
	    else{
	      h_map.at("mask_"+recos[r])->SetBinContent(i+1, 0);
	    }
	    h_map.at("mean_"+recos[r])->SetBinContent(i+1, mean_i);
	    h_map.at("mean_"+recos[r])->SetBinError(i+1, meanerr_i);
	    h_map.at("rms_"+recos[r])->SetBinContent(i+1, rms_i);
	    h_map.at("rms_"+recos[r])->SetBinError(i+1, rmserr_i);
	  }
	  // CB pdf
	  else{

	    if( !(hi_m->Integral() > minNumEvents && hi->Integral() > minNumEvents  &&  hi_m->GetMean()>75. && hi_m->GetMean()<105.) ){
	      h_map.at("mask_"+recos[r])->SetBinContent(i+1, 0);
	      continue;
	    }
	    h_map.at("mask_"+recos[r])->SetBinContent(i+1, 1);

	    //continue;

	    bool verbosity = false;
	    int printlevel = -1;
	    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
	    gErrorIgnoreLevel = 6001;

	    int status = -99;
	    int flag = -99;
	    float fr = -99.;
	    
	    RooRealVar mass("mass", Form("mass for bin %d", i), hi->GetXaxis()->GetXmin(), hi->GetXaxis()->GetXmax());
	    mass.setRange("r1", -10.0, 10.0);
	    
	    RooDataHist data("data", "", RooArgList(mass), hi );
	    
	    RooRealVar x0("x0", "", hi->GetMean(), mass.getMin(), mass.getMax() );
	    RooRealVar sigmaL("sigmaL", "", hi->GetRMS(), hi->GetRMS()*0.5, hi->GetRMS()*2 );
	    RooRealVar sigmaR("sigmaR", "", hi->GetRMS(), hi->GetRMS()*0.5, hi->GetRMS()*2 );
	    RooRealVar alphaL("alphaL", "", 1.0, 0.2, +10 );
	    RooRealVar alphaR("alphaR", "", 1.0, 0.2, +10 );
	    RooRealVar nL("nL", "", 2, 1, 100 );
	    RooRealVar nR("nR", "", 2, 1, 100 );
	    RooRealVar tau("tau", "", 0., -10, 10);    
	    RooRealVar frac("frac", "", 0.9, 0., 1.);
	    
	    RooCrystalBall pdf("pdf", "", mass, x0, sigmaL, sigmaR, alphaL, nL, alphaR, nR);
	    RooGaussian gaus("pdf", "", mass, x0, sigmaL);	    
	    RooExponential bkg("bkg", "", mass, tau);
	    RooAddPdf pdfTot("pdfTot", "", {pdf, bkg}, frac);

	    TString rname = "r1";
	    pdfTot.fixCoefRange( rname.Data() );
	    pdfTot.fixCoefNormalization(mass);

	    std::shared_ptr<RooFitResult> rfit{pdfTot.fitTo(data,
							    InitialHesse(true),
							    Minimizer("Minuit2"),
							    Range( rname.Data() ),
							    Save(), SumW2Error(true),
							    PrintLevel(printlevel),
							    Verbose(verbosity) )};
	    status = rfit->status();
	    fr = frac.getVal();
	    flag = 0;

	    if(fr>minFrac){
	      std::shared_ptr<RooFitResult> rn{pdf.fitTo(data,
							 InitialHesse(true),
							 Minimizer("Minuit2"),
							 Range( rname.Data() ),
							 Save(), SumW2Error(true),
							 PrintLevel(printlevel),
							 Verbose(verbosity) )}; 
	      rfit = rn;
	      flag = 1;
	    }
	    status = rfit->status();

	    if(status!=0){
	      rname = "r2";
	      mass.setRange( rname.Data(), -8.0, 8.0);
	      std::shared_ptr<RooFitResult> rn{pdf.fitTo(data,
							 InitialHesse(true),
							 Minimizer("Minuit2"),
							 Range( rname.Data() ),
							 Save(), SumW2Error(true),
							 PrintLevel(printlevel),
							 Verbose(verbosity) )};      
	      rfit = rn;
	      flag = 2;
	    }
	    status = rfit->status();

	    if(status!=0){
	      rname = "r3";
	      mass.setRange( rname.Data() , -6.0, 6.0);
	      std::shared_ptr<RooFitResult> rn{pdf.fitTo(data,
							 InitialHesse(true),
							 Minimizer("Minuit2"),
							 Range( rname.Data() ),
							 Save(), SumW2Error(true),
							 PrintLevel(printlevel),
							 Verbose(verbosity) )};      
	      rfit = rn;
	      flag = 3;
	    }
	    status = rfit->status();
    
	    if(status!=0){
	      rname = "r4";
	      mass.setRange( rname.Data() , -3.0, 3.0);
	      std::shared_ptr<RooFitResult> rn{gaus.fitTo(data,
							  InitialHesse(true),
							  Minimizer("Minuit2"),
							  Range( rname.Data() ),
							  Save(), SumW2Error(true),
							  PrintLevel(printlevel),
							  Verbose(verbosity) )};      
	      rfit = rn;
	      flag = 4;
	    }
	    status = rfit->status();

	    RooDerivative* der = 0;
	    if(flag==0){
	      der = pdfTot.derivative(mass, 1, 0.00005 );
	    }
	    else if(flag==1 || flag==2 || flag==3){
	      der = pdf.derivative(mass, 1, 0.00005 );
	    }
	    else if(flag==4){
	      der = gaus.derivative(mass, 1, 0.00005 );
	    }

	    TH2D* h_der = hCB_map.at("lnder_"+recos[r]);
	    for(int ib=0; ib<h_der->GetYaxis()->GetNbins();ib++){
	      double x = h_der->GetYaxis()->GetBinCenter(ib+1);
	      mass.setVal( x );
	      double fprime = der->getVal();
	      double f = 0.;
	      if(flag==0)
		f = pdfTot.getVal();
	      else if(flag==1 || flag==2 || flag==3)
		f = pdf.getVal();
	      else if(flag==4)
		f = gaus.getVal();	
	      h_der->SetBinContent(i+1,ib+1, fprime/f);
	    }
 	    //if(i==139){
	    //TH1D* proj = (TH1D*)h_der->ProjectionY("py", i+1, i+1);
	    //proj->Print("all");
	    //}
	  }
	  
	  delete hi;
	  delete hi_m;
	} 
      }
      
      fout->cd();
    
      for(unsigned int r = 0 ; r<recos.size(); r++){	  
	if(h_map["mask_"+recos[r]]!=0)
	  h_map["mask_"+recos[r]]->Write();
	if(!useCBpdf){
	  if(h_map["mean_"+recos[r]]!=0)
	    h_map["mean_"+recos[r]]->Write();
	  if(h_map["rms_"+recos[r]]!=0)
	    h_map["rms_"+recos[r]]->Write();
	}
	//else{
	//cout << "Writing " << hCB_map["lnder_"+recos[r]]->GetName() << endl;
	//hCB_map["lnder_"+recos[r]]->Write();
	//}
      }
    }

    else if(iter==2){
      if(saveHistos && fout->GetDirectory("postfit")==0){
	fout->mkdir("postfit");
      }
      TTree* treescales = new TTree("treescales","treescales");
      int inmassbins, indof, ibinIdx  ;
      float inevents, ibeta, ibetaErr, ialpha, ialphaErr, inu, inuErr, iprob, ichi2old, ichi2new; 
      treescales->Branch("nevents",&inevents,"nevents/F");
      treescales->Branch("beta",&ibeta,"beta/F");
      treescales->Branch("betaErr",&ibetaErr,"betaErr/F");
      treescales->Branch("alpha",&ialpha,"alpha/F");
      treescales->Branch("alphaErr",&ialphaErr,"alphaErr/F");
      treescales->Branch("nu",&inu,"nu/F");
      treescales->Branch("nuErr",&inuErr,"nuErr/F");
      treescales->Branch("prob",&iprob,"prob/F");
      treescales->Branch("chi2old",&ichi2old,"chi2old/F");
      treescales->Branch("chi2new",&ichi2new,"chi2new/F");
      treescales->Branch("nmassbins",&inmassbins,"nmassbins/I");
      treescales->Branch("ndof",&indof,"ndof/I");
      treescales->Branch("binIdx",&ibinIdx,"binIdx/I");
      
      TH1D* h_scales  = new TH1D("h_scales", "", n_bins, 0, double(n_bins));
      TH1D* h_widths  = new TH1D("h_widths", "", n_bins, 0, double(n_bins));
      TH1D* h_norms   = new TH1D("h_norms", "", n_bins, 0, double(n_bins));
      TH1D* h_probs   = new TH1D("h_probs", "", n_bins, 0, double(n_bins));
      TH1D* h_masks   = new TH1D("h_masks", "", n_bins, 0, double(n_bins));

      TH2D* h_data_2D   = (TH2D*)fout->Get("h_data_bin_m");
      TH2D* h_nom_2D    = (TH2D*)fout->Get("h_smear0_bin_m");
      TH1D* h_nom_mask  = (TH1D*)fout->Get("h_mask_smear0_bin_dm");
      TH2D* h_jscale_2D = (TH2D*)fout->Get("h_smear0_bin_jac_scale");
      TH2D* h_jwidth_2D = (TH2D*)fout->Get("h_smear0_bin_jac_width");

      for(unsigned int ibin=0; ibin<n_bins; ibin++){

	// skip empty bins
	if( h_nom_mask->GetBinContent(ibin+1)<0.5 ) continue;

	ibinIdx = ibin;
	
	TH1D* h_data_i    = (TH1D*)h_data_2D->ProjectionY( Form("h_data_i_%d",ibin ),     ibin+1, ibin+1 );
	TH1D* h_nom_i     = (TH1D*)h_nom_2D->ProjectionY( Form("h_nom_i_%d", ibin),       ibin+1, ibin+1 );
	TH1D* h_jscale_i  = (TH1D*)h_jscale_2D->ProjectionY( Form("h_jscale_i_%d", ibin), ibin+1, ibin+1 );
	TH1D* h_jwidth_i  = (TH1D*)h_jwidth_2D->ProjectionY( Form("h_jwidth_i_%d", ibin), ibin+1, ibin+1 );	

	if(scaleToData){
	  double data_norm_i = h_data_i->Integral();
	  double mc_norm_i = h_nom_i->Integral();
	  if(data_norm_i>0. && mc_norm_i>0.){
	    double sf_i = data_norm_i/mc_norm_i;
	    h_nom_i->Scale(sf_i);
	    h_jscale_i->Scale(sf_i);
	    h_jwidth_i->Scale(sf_i);
	  }
	}
	
	if(rebin>1){
	  h_data_i->Rebin(rebin);
	  h_nom_i->Rebin(rebin);
	  h_jscale_i->Rebin(rebin);
	  h_jwidth_i->Rebin(rebin);
	}
	
	//const float minNumEventsPerBin = 10.; 
	unsigned int n_mass_bins = 0;

	// skip bins with less than 4 high-stat (>10) mass bins
	for(int im = 1 ; im<=h_data_i->GetXaxis()->GetNbins(); im++){
	  if( h_data_i->GetBinContent(im)>minNumEventsPerBin ) n_mass_bins++;
	}
	if( n_mass_bins<minNumMassBins ){
	  h_scales->SetBinContent(ibin+1, 0.0);
	  h_widths->SetBinContent(ibin+1, 0.0);
	  h_norms->SetBinContent(ibin+1, 0.0);
	  h_probs->SetBinContent(ibin+1, 0.0);
	  h_masks->SetBinContent(ibin+1, 0.0);
	  continue;
	}

	inevents = h_data_i->Integral();
	inmassbins = n_mass_bins;
	
	// the data
	MatrixXd inv_sqrtV(n_mass_bins,n_mass_bins);
	MatrixXd inv_V(n_mass_bins,n_mass_bins);
	for(unsigned int ibm = 0; ibm<n_mass_bins; ibm++ ){
	  for(unsigned int jbm = 0; jbm<n_mass_bins; jbm++ ){
	    inv_sqrtV(ibm,jbm) = 0.;
	    inv_V(ibm,jbm) = 0.;
	  }
	}
	VectorXd y(n_mass_bins);
	VectorXd y0(n_mass_bins);
	VectorXd jscale(n_mass_bins);
	VectorXd jwidth(n_mass_bins);
	unsigned int bin_counter = 0;
	for(int im = 0 ; im<h_data_i->GetXaxis()->GetNbins(); im++){
	  if( h_data_i->GetBinContent(im+1)>minNumEventsPerBin ){
	    y(bin_counter)  = h_data_i->GetBinContent(im+1);
	    y0(bin_counter) = h_nom_i->GetBinContent(im+1);	    
	    jscale(bin_counter) = h_jscale_i->GetBinContent(im+1);
	    jwidth(bin_counter) = h_jwidth_i->GetBinContent(im+1);  
	    double mcErr_im = h_nom_i->GetBinError(im+1);
	    inv_V(bin_counter,bin_counter) = lumi>0. ?
	      1./(y(bin_counter)  + mcErr_im*mcErr_im ) :
	      1/(2*mcErr_im*mcErr_im);
	    //cout << TMath::Sqrt(y(bin_counter)) << " (+) " << h_nom_i->GetBinError(im+1) << endl;
	    inv_sqrtV(bin_counter,bin_counter) = TMath::Sqrt( inv_V(bin_counter,bin_counter) );
	    bin_counter++;
	  }
	}

	unsigned int n_fit_params = 3;
	if(!fitWidth) n_fit_params--;
	if(!fitNorm)  n_fit_params--;
	MatrixXd jac(n_mass_bins, n_fit_params);
	for(unsigned int ib=0; ib<n_mass_bins;ib++){
	  jac(ib, 0) = jscale(ib);
	  if(fitWidth)
	    jac(ib, 1) = jwidth(ib);
	  if(fitNorm)
	    jac(ib, 2) = y0(ib);
	}

	MatrixXd A = inv_sqrtV*jac;
	VectorXd b = inv_sqrtV*(y-y0);
	VectorXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
	MatrixXd C = (jac.transpose()*inv_V*jac).inverse();
	MatrixXd rho( C.rows(), C.rows() ) ;
	for(unsigned int ir = 0; ir<C.rows(); ir++){
	  for(unsigned int ic = 0; ic<C.rows(); ic++){
	    rho(ir,ic) = C(ir,ic)/TMath::Sqrt(C(ir,ir)*C(ic,ic));
	  }
	}
	MatrixXd chi2old = b.transpose()*b;
	MatrixXd chi2new = ((b - A*x).transpose())*(b-A*x);
	int ndof = n_mass_bins-n_fit_params;
	double chi2norm_old = chi2old(0,0)/(n_mass_bins);
	double chi2norm_new = chi2new(0,0)/ndof;
	cout << "Bin: " << ibin << ": mass fits with " << n_mass_bins << " mass bins: " << chi2norm_old << " ---> " << chi2norm_new << " (prob = " << TMath::Prob(chi2norm_new*ndof, ndof ) << ")" << endl;

	indof = ndof;
	ibeta = x(0);
	ibetaErr = TMath::Sqrt(C(0,0));
	ialpha = fitWidth ? x(1) : 0.;
	ialphaErr = fitWidth ? TMath::Sqrt(C(1,1)) : 0.;
	inu = fitNorm ? x(2) : 0.;
	inuErr = fitNorm ? TMath::Sqrt(C(2,2)) : 0.;
 	ichi2old = chi2norm_old;
	ichi2new = chi2norm_new;
	iprob =  TMath::Prob(chi2norm_new*ndof, ndof );	
	treescales->Fill();
	//cout << "Filling tree" << endl;
	
	h_scales->SetBinContent(ibin+1, ibeta+1.0);
	h_scales->SetBinError(ibin+1, ibetaErr);
	h_norms->SetBinContent(ibin+1, inu+1.0);
	h_norms->SetBinError(ibin+1, inuErr);
	h_widths->SetBinContent(ibin+1, ialpha+1.0);
	h_widths->SetBinError(ibin+1, ialphaErr);
	h_probs->SetBinContent(ibin+1, TMath::Prob(chi2norm_new*ndof, ndof ));
	h_probs->SetBinError(ibin+1, 0.);
	h_masks->SetBinContent(ibin+1, 1.0);

	if(saveHistos){
	  TH1D* h_pre_i   = (TH1D*)h_nom_i->Clone(Form("h_prefit_%d", ibin));
	  TH1D* h_post_i  = (TH1D*)h_nom_i->Clone(Form("h_postfit_%d", ibin));
	  unsigned int bin_counter = 0;
	  for(int im = 0 ; im<h_post_i->GetXaxis()->GetNbins(); im++){	  
	    if( h_data_i->GetBinContent(im+1)>minNumEventsPerBin ){
	      h_post_i->SetBinContent( im+1, y0(bin_counter)+(jac*x)(bin_counter) );
	      bin_counter++;
	    }
	  }
	  fout->cd("postfit/");
	  h_data_i->Write(Form("h_data_%d", ibin) ,TObject::kOverwrite);
	  h_pre_i->Write(TString(h_pre_i->GetName()) ,TObject::kOverwrite);
	  h_post_i->Write(TString(h_post_i->GetName()),TObject::kOverwrite);
	}	
      }
      
      fout->cd();
      h_scales->Write(0,TObject::kOverwrite);
      h_norms->Write(0,TObject::kOverwrite);
      h_widths->Write(0,TObject::kOverwrite);
      h_probs->Write(0,TObject::kOverwrite);
      h_masks->Write(0,TObject::kOverwrite);
      treescales->Write(0,TObject::kOverwrite);
	
      cout << h_masks->Integral() << " scales have been computed" << endl;
    }

    for(auto r : rans) delete r;
  }
  
  sw.Stop();

  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;

  fout->Close(); 

  delete ran0;
  
  return 1;
}
