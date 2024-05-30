#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include "TRandom3.h"
#include "TVector.h"
#include "TVectorT.h"
#include "TMath.h"
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

constexpr double mZ = 91.1;
constexpr double GZ = 2.5;

  
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
	("nevents",     value<int>()->default_value(100), "number of events")
	("lumi",        value<long>()->default_value(1000), "number of events")
	("tag",         value<std::string>()->default_value("closure"), "run type")
	("run",         value<std::string>()->default_value("closure"), "run type")
	("dofits",      bool_switch()->default_value(false), "")
	("savehistos",  bool_switch()->default_value(false), "")
	("bias",        value<int>()->default_value(0), "bias")
	("event_cut",   value<int>()->default_value(10), "bias")
	("rebin",       value<int>()->default_value(1), "rebin")
	("seed",        value<int>()->default_value(4357), "seed");

      store(parse_command_line(argc, argv, desc), vm);
      notify(vm);
      if (vm.count("help")){
	std::cout << desc << '\n';
	return 0;
      }
      if (vm.count("nevents"))    std::cout << "Number of events: " << vm["nevents"].as<int>() << '\n';
      if (vm.count("tag"))        std::cout << "Tag: " << vm["tag"].as<std::string>() << '\n';
      if (vm.count("run"))        std::cout << "Run: " << vm["run"].as<std::string>() << '\n';
    }
  catch (const error &ex)
    {
      std::cerr << ex.what() << '\n';
    }

  int nevents     = vm["nevents"].as<int>();
  long lumi       = vm["lumi"].as<long>();
  std::string tag = vm["tag"].as<std::string>();
  std::string run = vm["run"].as<std::string>();
  int bias        = vm["bias"].as<int>();
  int seed        = vm["seed"].as<int>();
  int event_cut   = vm["event_cut"].as<int>();
  int rebin       = vm["rebin"].as<int>();
  bool dofits     = vm["dofits"].as<bool>();
  bool savehistos     = vm["savehistos"].as<bool>();

  vector<float> pt_edges  = {25, 30, 35, 40, 45, 55}; 
  vector<float> eta_edges = {-2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0,
			      0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};
		  
  unsigned int n_pt_bins  = pt_edges.size()-1;
  unsigned int n_eta_bins = eta_edges.size()-1;
  int n_bins = n_pt_bins*n_pt_bins*n_eta_bins*n_eta_bins;

  std::vector<string> recos = {"reco", "smear0", "smear1"}; 

  std::map<string, TH1D*> h_map;
  for(unsigned int r = 0; r<recos.size(); r++){
    h_map.insert( std::make_pair<string, TH1D* >("mean_"+recos[r], 0 ) );
    h_map.insert( std::make_pair<string, TH1D* >("rms_"+recos[r],  0 ) );
    h_map.insert( std::make_pair<string, TH1D* >("mask_"+recos[r],  0 ) );
  }

  // map of positions in RVecF "masses"
  std::map<string, unsigned int> idx_map;
  idx_map.insert( std::make_pair<string, unsigned int >("reco",   1 ) );
  idx_map.insert( std::make_pair<string, unsigned int >("smear0", 2 ) );
  idx_map.insert( std::make_pair<string, unsigned int >("smear1", 3 ) );

  TFile* fout = TFile::Open(("./massscales_"+tag+"_"+run+".root").c_str(), !dofits ? "RECREATE" : "UPDATE");
  
  for(unsigned int iter=0; iter<3; iter++){

    if(dofits && iter!=2) continue;    
    cout << "Iter " << iter << endl;
    //TTree* tree = new TTree("tree", "tree");

    ROOT::RDataFrame d( "Events",
			{//"/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_040854/0000/NanoV9MCPostVFP_1.root",
			  "/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_040854/0000/NanoV9MCPostVFP_*.root",
			  "/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_041233/0000/NanoV9MCPostVFP_*.root",
			  "/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_041233/0001/NanoV9MCPostVFP_*.root",
			  "/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_PDFExt_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_041233/0002/NanoV9MCPostVFP_*.root"
			} );
    unsigned int nslots = d.GetNSlots();
    std::vector<TRandom3*> rans = {};
    for(unsigned int i = 0; i < nslots; i++){
      rans.emplace_back( new TRandom3(seed + i*10) );
    }

    auto dlast = std::make_unique<RNode>(d);
    dlast = std::make_unique<RNode>(dlast->DefineSlot("x", [&](unsigned int nslot)->double{
      double out;
      out = rans[nslot]->Uniform(0.0, 1.0);
      return out;
    } ));
    
    dlast = std::make_unique<RNode>(dlast->Define("weight", [](float weight)->float{
      return std::copysign(1.0, weight);
    }, {"Generator_weight"} ));
    
    dlast = std::make_unique<RNode>(dlast->Define("idxs", [&](UInt_t nMuon, RVecB Muon_looseId, RVecF Muon_dxybs, RVecB Muon_isGlobal,
							      RVecB Muon_highPurity, RVecB Muon_mediumId, RVecF Muon_pfRelIso04_all,
							      RVecF Muon_pt, RVecF Muon_eta)->RVecUI {
      RVecUI out;
      for(unsigned int i = 0; i < nMuon; i++){
	if( Muon_looseId[i] && TMath::Abs(Muon_dxybs[i]) < 0.05 && Muon_isGlobal[i] && Muon_highPurity[i] && Muon_mediumId[i] && Muon_pfRelIso04_all[i]<0.15 &&
	    Muon_pt[i] >= pt_edges[0] && Muon_pt[i] < pt_edges[ n_pt_bins ]  && Muon_eta[i]>=eta_edges[0] && Muon_eta[i]<=eta_edges[ n_eta_bins ]  ){
	  out.emplace_back(i);
	}
	else{
	  continue;
	  /*
	    cout << i << ":  " <<  Muon_looseId[i] << ", "
	    << TMath::Abs(Muon_dxybs[i]) << ", "
	    << Muon_isGlobal[i] << ", "
	    << Muon_highPurity[i] << ", "
	    << Muon_mediumId[i] << ", "
	    << Muon_pt[i]<< ", "
	    << TMath::Abs(Muon_eta[i])<< endl;
	  */
	}
      }
      //cout << ">>>>> " << out.size() << endl;
      return out;
    }, {"nMuon", "Muon_looseId", "Muon_dxybs", "Muon_isGlobal",
	"Muon_highPurity","Muon_mediumId", "Muon_pfRelIso04_all",
	"Muon_pt", "Muon_eta" } ));

    dlast = std::make_unique<RNode>(dlast->Filter( [](RVecUI idxs, RVecI Muon_charge, bool HLT_IsoMu24 ){
      if( idxs.size()!=2 || !HLT_IsoMu24) return false;
      if( Muon_charge[idxs[0]]*Muon_charge[idxs[1]] > 0 ) return false;
      //cout << "pass" << endl;
      return true;
    }, {"idxs", "Muon_charge", "HLT_IsoMu24"} ));

    dlast = std::make_unique<RNode>(dlast->DefineSlot("Muon_ptsmear", [&](unsigned int nslot, RVecUI idxs,
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
	float scale_smear0 = 1.0;
	out.emplace_back( rans[nslot]->Gaus(gmuP.Pt()*scale_smear0, gmuP.Pt()*0.02) );
	out.emplace_back( rans[nslot]->Gaus(gmuM.Pt()*scale_smear0, gmuM.Pt()*0.02) );
	float scale_smear1 = 1.001;
	out.emplace_back( rans[nslot]->Gaus(gmuP.Pt()*scale_smear1, gmuP.Pt()*0.02) );
	out.emplace_back( rans[nslot]->Gaus(gmuM.Pt()*scale_smear1, gmuM.Pt()*0.02) );
      }
      else{
	out.emplace_back(0.0);
	out.emplace_back(0.0);
	out.emplace_back(0.0);
	out.emplace_back(0.0);
      }
      return out;
	}, {"idxs",
	"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge",
	"nGenPart", "GenPart_status", "GenPart_statusFlags", "GenPart_pdgId",
	"GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass"} ));

    dlast = std::make_unique<RNode>(dlast->Define("indexes", [&](RVecUI idxs, RVecF Muon_pt, RVecF Muon_eta, RVecI Muon_charge, RVecF Muon_ptsmear)-> RVecUI {
      unsigned int idxP = Muon_charge[idxs[0]]>0 ? idxs[0] : idxs[1];
      unsigned int idxM = Muon_charge[idxs[0]]>0 ? idxs[1] : idxs[0];
      float ptP  = Muon_pt[idxP];
      float ptM  = Muon_pt[idxM];
      float ptsmear0P = Muon_ptsmear[0];
      float ptsmear0M = Muon_ptsmear[1];
      float ptsmear1P = Muon_ptsmear[2];
      float ptsmear1M = Muon_ptsmear[3];
      float etaP = Muon_eta[idxP];
      float etaM = Muon_eta[idxM];
      RVecUI out;
      out.emplace_back(n_bins);
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
		  ptsmear0P>=pt_p_low   && ptsmear0P<pt_p_up &&
		  ptsmear0M>=pt_m_low   && ptsmear0M<pt_m_up 
		  ) out[1] = ibin;
	      if( etaP>=eta_p_low && etaP<eta_p_up &&
		  etaM>=eta_m_low && etaM<eta_m_up &&
		  ptsmear1P>=pt_p_low   && ptsmear1P<pt_p_up &&
		  ptsmear1M>=pt_m_low   && ptsmear1M<pt_m_up 
		  ) out[2] = ibin;	      
	      ibin++;
	    }
	  }
	}
      }
      return out;
    }, {"idxs", "Muon_pt", "Muon_eta", "Muon_charge", "Muon_ptsmear"} ));

    for(unsigned int r = 0 ; r<recos.size(); r++){
      dlast = std::make_unique<RNode>(dlast->Define( TString(("index_"+recos[r]).c_str()), [r](RVecUI indexes){
	return indexes.at(r);
      }, {"indexes"} ));
    }
    
    dlast = std::make_unique<RNode>(dlast->Define("masses", [&](RVecUI idxs,
								RVecF Muon_pt, RVecF Muon_eta, RVecF Muon_phi, RVecF Muon_mass, RVecI Muon_charge,
								UInt_t nGenPart, RVecI GenPart_status, RVecI GenPart_statusFlags, RVecI GenPart_pdgId,
								RVecF GenPart_pt, RVecF GenPart_eta, RVecF GenPart_phi, RVecF GenPart_mass,
								RVecF Muon_ptsmear)->RVecF {
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

	float scale_smear0 = 1.0;
	ROOT::Math::PtEtaPhiMVector muP_smear0( Muon_ptsmear.at(0),
					       Muon_eta[ idxP ], Muon_phi[ idxP ], Muon_mass[ idxP ] );
	ROOT::Math::PtEtaPhiMVector muM_smear0( Muon_ptsmear.at(1),
					       Muon_eta[ idxM ], Muon_phi[ idxM ], Muon_mass[ idxM ] );      
	out.emplace_back( (muP_smear0 + muM_smear0).M() );	
	ROOT::Math::PtEtaPhiMVector muP_smear1( Muon_ptsmear.at(2),
					       Muon_eta[ idxP ], Muon_phi[ idxP ], Muon_mass[ idxP ] );
	ROOT::Math::PtEtaPhiMVector muM_smear1( Muon_ptsmear.at(3),
					       Muon_eta[ idxM ], Muon_phi[ idxM ], Muon_mass[ idxM ] );      
	out.emplace_back( (muP_smear1 + muM_smear1).M() );
      }
      
      return out;
    }, {"idxs",
	"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge",
	"nGenPart", "GenPart_status", "GenPart_statusFlags", "GenPart_pdgId",
	"GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass",
	"Muon_ptsmear"} ));

    // single column
    dlast = std::make_unique<RNode>(dlast->Define("gen_m", [](RVecF masses){
      return masses.size()>0 ? masses.at(0) : -99.;
    }, {"masses"} ));

    for(unsigned int r = 0 ; r<recos.size(); r++){
      unsigned int mpos = idx_map.at(recos[r]);
      dlast = std::make_unique<RNode>(dlast->Define(TString( (recos[r]+"_m").c_str() ), [mpos](RVecF masses){
	return masses.size()>0 ? masses.at( mpos ) : -99.;
      }, {"masses"} ));      
      dlast = std::make_unique<RNode>(dlast->Define(TString( (recos[r]+"_dm").c_str() ), [mpos](RVecF masses){
	return masses.size()>0 ? masses.at( mpos ) - masses.at(0) : -99.;
      }, {"masses"} ));            
    }
    
    dlast = std::make_unique<RNode>(dlast->Define("weights_jac", [n_bins,recos,h_map,idx_map](RVecF masses, RVecUI indexes)->RVecF{
      RVecF out;
      if(masses.size()==0){
	for(unsigned int r = 0 ; r<recos.size(); r++){
	  out.emplace_back(0.0);
	  out.emplace_back(0.0);
	}
	return out;
      }
      
      float gen_m  = masses.at(0);
      for(unsigned int r = 0 ; r<recos.size(); r++){
	unsigned int rpos = idx_map.at(recos[r]);
	TH1D* h_mean = h_map.at("mean_"+recos[r]);
	TH1D* h_rms  = h_map.at("rms_"+recos[r]);
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
      return out;
    }, {"masses", "indexes"} ));

    for(unsigned int r = 0 ; r<recos.size(); r++){
      unsigned int jpos = (idx_map.at(recos[r])-1)*2;
      dlast = std::make_unique<RNode>(dlast->Define( TString((recos[r]+"_jscale_weight").c_str()), [jpos](RVecF weights_jac, float weight)->float{
	//unisgned int jpos = (idx_map[recos[r]]-1)*2;
	return weights_jac.at( jpos )*weight;
      }, {"weights_jac", "weight"} ));
      dlast = std::make_unique<RNode>(dlast->Define( TString((recos[r]+"_jwidth_weight").c_str()), [jpos](RVecF weights_jac, float weight)->float{
	//unisgned int jpos = (idx_map[recos[r]]-1)*2 + 1;
	return weights_jac.at( jpos+1 )*weight;
      }, {"weights_jac", "weight"} ));
    }
    
    std::vector<ROOT::RDF::RResultPtr<TH1D> > histos1D;
    std::vector<ROOT::RDF::RResultPtr<TH2D> > histos2D;
    
    const int x_nbins   = 40;
    const double x_low  = 70.0;
    const double x_high = 110.0;
  
    //histos1D.emplace_back(dlast->Histo1D({"h_gen_m", "nominal", x_nbins, x_low, x_high}, "gen_m", "weight"));      
    //histos1D.emplace_back(dlast->Histo1D({"h_reco_m", "nominal", x_nbins, x_low, x_high}, "reco_m", "weight"));
    //histos1D.emplace_back(dlast->Histo1D({"h_smear_m", "nominal", x_nbins, x_low, x_high}, "smear_m", "weight"));

    if(iter==0){
      histos2D.emplace_back(dlast->Histo2D({"h_gen_bin_m",    "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index_reco", "gen_m", "weight"));
      for(unsigned int r = 0 ; r<recos.size(); r++){
	histos2D.emplace_back(dlast->Histo2D({ "h_"+TString(recos[r].c_str())+"_bin_m",    "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_m", "weight"));
	histos2D.emplace_back(dlast->Histo2D({ "h_"+TString(recos[r].c_str())+"_bin_dm",   "nominal", n_bins, 0, double(n_bins), 24, -6.0, 6.0},          "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_dm", "weight"));
      }
      auto colNames = dlast->GetColumnNames();
      double total = *(dlast->Count());  
      std::cout << colNames.size() << " columns created. Total event count is " << total  << std::endl;
    }
    else if(iter==1){
      for(unsigned int r = 0 ; r<recos.size(); r++){
	histos2D.emplace_back(dlast->Histo2D({"h_"+TString(recos[r].c_str())+"_bin_jac_scale", "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_m", TString(recos[r].c_str())+"_jscale_weight"));
	histos2D.emplace_back(dlast->Histo2D({"h_"+TString(recos[r].c_str())+"_bin_jac_width", "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index_"+TString(recos[r].c_str()), TString(recos[r].c_str())+"_m", TString(recos[r].c_str())+"_jwidth_weight"));
      }
    }
    
    fout->cd();
    std::cout << "Writing histos..." << std::endl;
    double sf = 1.0; //double(lumi)/double(nevents);
    for(auto h : histos1D){
      h->Scale(sf);
      string h_name = std::string(h->GetName());
      h->Write();
    }
    for(auto h : histos2D){
      h->Scale(sf);
      string h_name = std::string(h->GetName());
      std::cout << "Total number of events in 2D histo " << h_name << ": " << h->GetEntries() << std::endl;
      h->Write();
    }    
    std::cout << "Total slots: " << dlast->GetNSlots() << std::endl;

    // fill histos
    if(iter==0){
      for(unsigned int r = 0 ; r<recos.size(); r++){
	TH2D* h_reco_dm = (TH2D*)fout->Get(TString( ("h_"+recos[r]+"_bin_dm").c_str()) );
	TH2D* h_reco_m  = (TH2D*)fout->Get(TString( ("h_"+recos[r]+"_bin_m").c_str()) );
	if( h_reco_dm==0 || h_reco_m==0 ){
	  cout << "h_reco_dm/h_reco_m NOT FOUND" << endl;
	  continue;
	}
	h_map["mean_"+recos[r]] = new TH1D( TString( ("h_mean_"+recos[r]+"_bin_dm").c_str() ),"", n_bins, 0, double(n_bins));
	h_map["rms_"+recos[r]]  = new TH1D( TString( ("h_rms_"+recos[r]+"_bin_dm").c_str() ),"", n_bins, 0, double(n_bins));
	h_map["mask_"+recos[r]] = new TH1D( TString( ("h_mask_"+recos[r]+"_bin_dm").c_str() ),"", n_bins, 0, double(n_bins));
	for(unsigned int i = 0; i<n_bins; i++ ){
	  TString projname(Form("bin_%d_", i));
	  projname += TString( recos[r].c_str() );
	  TH1D* hi   = (TH1D*)h_reco_dm->ProjectionY( projname+"_dm", i+1, i+1 );
	  TH1D* hi_m = (TH1D*)h_reco_m->ProjectionY( projname+"_m", i+1, i+1 );
	  double mean_i = 0.0;
	  double meanerr_i = 0.0;
	  double rms_i = 0.0;
	  double rmserr_i = 0.0;
	  //cout << hi_m->Integral() << ", " << hi->Integral() << ", " << hi_m->GetMean() << endl;
	  if( hi_m->Integral() > nevents && hi->Integral() > nevents  &&  hi_m->GetMean()>75. && hi_m->GetMean()<105. ){
	    h_map.at("mask_"+recos[r])->SetBinContent(i+1, 1);
	    TF1* gf = new TF1("gf","[0]/TMath::Sqrt(2)/[2]*TMath::Exp( -0.5*(x-[1])*(x-[1])/[2]/[2] )",
			      hi->GetXaxis()->GetBinLowEdge(1), hi->GetXaxis()->GetBinUpEdge( hi->GetXaxis()->GetNbins() ));      
	    gf->SetParameter(0, hi->Integral());
	    gf->SetParameter(1, hi->GetMean());
	    gf->SetParameter(2, hi->GetRMS() );
	    hi->Fit("gf", "QR", "", -6.0, 6.0);
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
      }
      
      fout->cd();
    
      for(unsigned int r = 0 ; r<recos.size(); r++){	  
	h_map["mean_"+recos[r]]->Write();
	h_map["rms_"+recos[r]]->Write();
	h_map["mask_"+recos[r]]->Write();
      }
    }

    if(iter==2){

      if(savehistos)
	fout->mkdir("postfit");
      
      TH1D* h_scales  = new TH1D("h_scales", "", n_bins, 0, double(n_bins));
      TH1D* h_widths  = new TH1D("h_widths", "", n_bins, 0, double(n_bins));
      TH1D* h_norms   = new TH1D("h_norms", "", n_bins, 0, double(n_bins));
      TH1D* h_probs   = new TH1D("h_probs", "", n_bins, 0, double(n_bins));
      TH1D* h_masks   = new TH1D("h_masks", "", n_bins, 0, double(n_bins));

      TH2D* h_data_2D   = (TH2D*)fout->Get("h_smear1_bin_m");
      TH1D* h_nom_mask  = (TH1D*)fout->Get("h_mask_smear0_bin_dm");
      TH2D* h_nom_2D    = (TH2D*)fout->Get("h_smear0_bin_m");
      TH2D* h_jscale_2D = (TH2D*)fout->Get("h_smear0_bin_jac_scale");
      TH2D* h_jwidth_2D = (TH2D*)fout->Get("h_smear0_bin_jac_width");

      for(unsigned int ibin=0; ibin<n_bins; ibin++){

	if( h_nom_mask->GetBinContent(ibin+1)<0.5 ) continue;

	TH1D* h_data_i    = (TH1D*)h_data_2D->ProjectionY( Form("h_data_i_%d",ibin ),   ibin+1, ibin+1 );
	TH1D* h_nom_i     = (TH1D*)h_nom_2D->ProjectionY( Form("h_nom_i_%d", ibin),       ibin+1, ibin+1 );
	TH1D* h_jscale_i  = (TH1D*)h_jscale_2D->ProjectionY( Form("h_jscale_i_%d", ibin), ibin+1, ibin+1 );
	TH1D* h_jwidth_i  = (TH1D*)h_jwidth_2D->ProjectionY( Form("h_jwidth_i_%d", ibin), ibin+1, ibin+1 );	

	if(rebin>1){
	  h_data_i->Rebin(rebin);
	  h_nom_i->Rebin(rebin);
	  h_jscale_i->Rebin(rebin);
	  h_jwidth_i->Rebin(rebin);
	}
	
	//const float event_cut = 10.; 
	unsigned int n_mass_bins = 0;
	// skip bins with less than 4 high-stat (>10) mass bins
	for(int im = 1 ; im<=h_data_i->GetXaxis()->GetNbins(); im++){
	  if( h_data_i->GetBinContent(im)>event_cut ) n_mass_bins++;
	}
	if( n_mass_bins<4 ){
	  h_scales->SetBinContent(ibin+1, 0.0);
	  h_widths->SetBinContent(ibin+1, 0.0);
	  h_norms->SetBinContent(ibin+1, 0.0);
	  h_probs->SetBinContent(ibin+1, 0.0);
	  h_masks->SetBinContent(ibin+1, 0.0);
	  continue;
	}

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
	  if( h_data_i->GetBinContent(im+1)>event_cut ){
	    y0(bin_counter) = h_nom_i->GetBinContent(im+1);	    
	    y(bin_counter)  = h_data_i->GetBinContent(im+1) - y0(bin_counter);
	    jscale(bin_counter) = h_jscale_i->GetBinContent(im+1);
	    jwidth(bin_counter) = h_jwidth_i->GetBinContent(im+1);  
	    inv_V(bin_counter,bin_counter) = 1./(y(bin_counter) + h_nom_i->GetBinError(im+1)*h_nom_i->GetBinError(im+1) );
	    inv_sqrtV(bin_counter,bin_counter) = TMath::Sqrt( inv_V(bin_counter,bin_counter) );
	    bin_counter++;
	  }
	}
	
	MatrixXd jac(n_mass_bins, 3);
	for(unsigned int ib=0; ib<n_mass_bins;ib++){
	  jac(ib, 0) = y0(ib);
	  jac(ib, 1) = jscale(ib);
	  jac(ib, 2) = jwidth(ib);
	}

	MatrixXd A = inv_sqrtV*jac;
	VectorXd b = inv_sqrtV*y;
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
	int ndof = n_mass_bins-3;
	double chi2norm_old = chi2old(0,0)/ndof;
	double chi2norm_new = chi2new(0,0)/ndof;
	cout << "Bin: " << ibin << ": mass fits with " << n_mass_bins << " mass bins: " << chi2norm_old << " ---> " << chi2norm_new << " (prob = " << TMath::Prob(chi2norm_new*ndof, ndof ) << endl;      	
	h_scales->SetBinContent(ibin+1, x(1)+1.0);
	h_scales->SetBinError(ibin+1, TMath::Sqrt(C(1,1)));
	h_norms->SetBinContent(ibin+1, x(0)+1.0);
	h_norms->SetBinError(ibin+1, TMath::Sqrt(C(0,0)));
	h_widths->SetBinContent(ibin+1, x(2)+1.0);
	h_widths->SetBinError(ibin+1, TMath::Sqrt(C(2,2)));
	h_probs->SetBinContent(ibin+1, TMath::Prob(chi2norm_new*ndof, ndof ));
	h_probs->SetBinError(ibin+1, 0.);
	h_masks->SetBinContent(ibin+1, 1.0);

	if(savehistos){
	  TH1D* h_pre_i   = (TH1D*)h_nom_i->Clone(Form("h_prefit_%d", ibin));
	  TH1D* h_post_i  = (TH1D*)h_nom_i->Clone(Form("h_postfit_%d", ibin));
	  unsigned int bin_counter = 0;
	  for(int im = 0 ; im<h_post_i->GetXaxis()->GetNbins(); im++){	  
	    if( h_data_i->GetBinContent(im+1)>event_cut ){
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

      cout << h_masks->Integral() << " scales have been computed" << endl;
    }

  }
  
  sw.Stop();

  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;

  fout->Close(); 
  //for(auto r : rans) delete r;

  return 1;
}
