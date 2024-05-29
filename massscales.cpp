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
	("nevents",     value<long>()->default_value(1000), "number of events")
	("lumi",        value<long>()->default_value(1000), "number of events")
	("tag",         value<std::string>()->default_value("closure"), "run type")
	("run",         value<std::string>()->default_value("closure"), "run type")
	("bias",        value<int>()->default_value(0), "bias")
	("seed",        value<int>()->default_value(4357), "seed");

      store(parse_command_line(argc, argv, desc), vm);
      notify(vm);
      if (vm.count("help")){
	std::cout << desc << '\n';
	return 0;
      }
      if (vm.count("nevents"))    std::cout << "Number of events: " << vm["nevents"].as<long>() << '\n';
      if (vm.count("tag"))        std::cout << "Tag: " << vm["tag"].as<std::string>() << '\n';
      if (vm.count("run"))        std::cout << "Run: " << vm["run"].as<std::string>() << '\n';
    }
  catch (const error &ex)
    {
      std::cerr << ex.what() << '\n';
    }

  long nevents    = vm["nevents"].as<long>();
  long lumi       = vm["lumi"].as<long>();
  std::string tag = vm["tag"].as<std::string>();
  std::string run = vm["run"].as<std::string>();
  int bias        = vm["bias"].as<int>();
  int seed        = vm["seed"].as<int>();

  vector<float> pt_edges  = {25, 30, 35, 40, 45, 55}; 
  vector<float> eta_edges = {-2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0,
			      0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};
		  
  unsigned int n_pt_bins  = pt_edges.size()-1;
  unsigned int n_eta_bins = eta_edges.size()-1;
  int n_bins = n_pt_bins*n_pt_bins*n_eta_bins*n_eta_bins;

  TH1D* h_mean_reco_bin_dm = 0;
  TH1D* h_rms_reco_bin_dm = 0;

  TFile* fout = TFile::Open(("./massscales_"+tag+"_"+run+".root").c_str(), "RECREATE");
  
  for(unsigned int iter=0; iter<2; iter++){

    cout << "Iter " << iter << endl;
  //TTree* tree = new TTree("tree", "tree");

    ROOT::RDataFrame d( "Events",
			{"/scratch/wmass/y2016/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_040854/0000/NanoV9MCPostVFP_*.root",
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
  
    dlast = std::make_unique<RNode>(dlast->Define("index", [&](RVecUI idxs, RVecF Muon_pt, RVecF Muon_eta, RVecI Muon_charge)-> unsigned int {
      unsigned int idxP = Muon_charge[idxs[0]]>0 ? idxs[0] : idxs[1];
      unsigned int idxM = Muon_charge[idxs[0]]>0 ? idxs[1] : idxs[0];
      float ptP  = Muon_pt[idxP];
      float ptM  = Muon_pt[idxM];
      float etaP = Muon_eta[idxP];
      float etaM = Muon_eta[idxM];
      unsigned int out = 0;    
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
		  ) return out;
	      out++;
	    }
	  }
	}
      }
      return 0;
    }, {"idxs", "Muon_pt", "Muon_eta", "Muon_charge"} ));
    
  
    dlast = std::make_unique<RNode>(dlast->DefineSlot("masses", [&](unsigned int nslot, RVecUI idxs,
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
	ROOT::Math::PtEtaPhiMVector muP_smear( rans[nslot]->Gaus( gmuP.Pt(), gmuP.Pt()*0.01 ) ,
					       Muon_eta[ idxP ], Muon_phi[ idxP ], Muon_mass[ idxP ] );
	ROOT::Math::PtEtaPhiMVector muM_smear( rans[nslot]->Gaus( gmuM.Pt(), gmuM.Pt()*0.01 ) ,
					       Muon_eta[ idxM ], Muon_phi[ idxM ], Muon_mass[ idxM ] );      
	out.emplace_back( (gmuP + gmuM).M() );
	out.emplace_back( (muP + muM).M() );
	out.emplace_back( (muP_smear + muM_smear).M() );
      }
      
      return out;
    }, {"idxs",
	"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge",
	"nGenPart", "GenPart_status", "GenPart_statusFlags", "GenPart_pdgId",
	"GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass"} ));
    
    dlast = std::make_unique<RNode>(dlast->Define("gen_m", [](RVecF masses){
      return masses.size()>1 ? masses.at(0) : -99.;
    }, {"masses"} ));
    dlast = std::make_unique<RNode>(dlast->Define("reco_m", [](RVecF masses){
      return masses.size()>1 ? masses.at(1) : -99.;
    }, {"masses"} ));
    dlast = std::make_unique<RNode>(dlast->Define("smear_m", [](RVecF masses){
      return masses.size()>2 ? masses.at(2) : -99.;
    }, {"masses"} ));
    dlast = std::make_unique<RNode>(dlast->Define("smear_dm", [](RVecF masses){
      return masses.size()>2 ? masses.at(2)-masses.at(0) : -99.;
    }, {"masses"} ));
    dlast = std::make_unique<RNode>(dlast->Define("reco_dm", [](RVecF masses){
      return masses.size()>2 ? masses.at(1)-masses.at(0) : -99.;
    }, {"masses"} ));

    dlast = std::make_unique<RNode>(dlast->Define("weights_jac", [h_mean_reco_bin_dm, h_rms_reco_bin_dm](RVecF masses, unsigned int index)->RVecF{
      RVecF out;
      if(masses.size()<2){
	out.emplace_back(0.0);
	out.emplace_back(0.0);
	return out;
      }
      float reco_m = masses.at(1);
      float gen_m  = masses.at(0);
      float reco_delta = h_mean_reco_bin_dm->GetBinContent(index+1);
      float reco_sigma = h_rms_reco_bin_dm->GetBinContent(index+1);
      float reco_jscale = reco_sigma>0. ? +(reco_m - (gen_m+reco_delta) )*(gen_m+reco_delta)/reco_sigma/reco_sigma : 0.0;
      float reco_jwidth = reco_sigma>0. ? +(reco_m - (gen_m+reco_delta) )*(reco_m - (gen_m+reco_delta) )/reco_sigma/reco_sigma - 1.0 : 0.0;
      out.emplace_back(reco_jscale);
      out.emplace_back(reco_jwidth);
      return out;
    }, {"masses", "index"} ));

    dlast = std::make_unique<RNode>(dlast->Define("reco_jscale_weight", [](RVecF weights_jac, float weight)->float{
      return weights_jac.at(0)*weight;
    }, {"weights_jac", "weight"} ));
    dlast = std::make_unique<RNode>(dlast->Define("reco_jwidth_weight", [](RVecF weights_jac, float weight)->float{
      return weights_jac.at(1)*weight;
    }, {"weights_jac", "weight"} ));

    
    std::vector<ROOT::RDF::RResultPtr<TH1D> > histos1D;
    std::vector<ROOT::RDF::RResultPtr<TH2D> > histos2D;
    
    const int x_nbins   = 40;
    const double x_low  = 70.0;
    const double x_high = 110.0;
  
    //histos1D.emplace_back(dlast->Histo1D({"h_gen_m", "nominal", x_nbins, x_low, x_high}, "gen_m", "weight"));      
    //histos1D.emplace_back(dlast->Histo1D({"h_reco_m", "nominal", x_nbins, x_low, x_high}, "reco_m", "weight"));
    //histos1D.emplace_back(dlast->Histo1D({"h_smear_m", "nominal", x_nbins, x_low, x_high}, "smear_m", "weight"));

    if(iter==0){
      histos1D.emplace_back(dlast->Histo1D({"h_smear_m_iter0", "nominal", x_nbins, x_low, x_high}, "smear_m", "weight"));
      histos2D.emplace_back(dlast->Histo2D({"h_gen_bin_m",    "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index", "gen_m", "weight"));
      histos2D.emplace_back(dlast->Histo2D({"h_reco_bin_m",   "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index", "reco_m", "weight"));
      histos2D.emplace_back(dlast->Histo2D({"h_reco_bin_dm",  "nominal", n_bins, 0, double(n_bins), 24, -6.0, 6.0}, "index", "reco_dm", "weight"));
      //histos2D.emplace_back(dlast->Histo2D({"h_smear_bin_m",  "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index", "smear_m", "weight"));
      //histos2D.emplace_back(dlast->Histo2D({"h_smear_bin_dm", "nominal", n_bins, 0, double(n_bins), 24, -6.0, 6.0}, "index", "smear_dm", "weight"));
    }
    else if(iter==1){
      histos1D.emplace_back(dlast->Histo1D({"h_smear_m_iter1", "nominal", x_nbins, x_low, x_high}, "smear_m", "weight"));
      histos2D.emplace_back(dlast->Histo2D({"h_reco_bin_jac_scale", "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index", "reco_m", "reco_jscale_weight"));
      histos2D.emplace_back(dlast->Histo2D({"h_reco_bin_jac_width", "nominal", n_bins, 0, double(n_bins), x_nbins, x_low, x_high}, "index", "reco_m", "reco_jwidth_weight"));
    }
    
    auto colNames = dlast->GetColumnNames();
    std::cout << colNames.size() << " columns created" << std::endl;
    double total = *(dlast->Count());  

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

    if(iter>0) continue;

    TH2D* h_reco_bin_dm = (TH2D*)fout->Get("h_reco_bin_dm");
    TH2D* h_gen_bin_m = (TH2D*)fout->Get("h_gen_bin_m");
    if( h_reco_bin_dm==0 || h_gen_bin_m==0 ){
      cout << "h_reco_bin_dm/h_gen_bin_m NOT FOUND" << endl;
      continue;
    }
    h_mean_reco_bin_dm = new TH1D("h_mean_reco_bin_dm","", n_bins, 0, double(n_bins));
    h_rms_reco_bin_dm  = new TH1D("h_rms_reco_bin_dm", "", n_bins, 0, double(n_bins));
    for(unsigned int i = 0; i<n_bins; i++ ){
      TH1D* hi = (TH1D*)h_reco_bin_dm->ProjectionY( Form("h_reco_bin_dm_%d", i), i+1, i+1 );
      TH1D* hi_m = (TH1D*)h_gen_bin_m->ProjectionY( Form("h_gen_bin_m_%d", i), i+1, i+1 );
      double mean_i = 0.0;
      double meanerr_i = 0.0;
      double rms_i = 0.0;
      double rmserr_i = 0.0;
      if( hi->Integral() > 100. &&  hi_m->GetMean()>75. && hi_m->GetMean()<105. ){
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
	delete gf;
      }
      h_mean_reco_bin_dm->SetBinContent(i+1, mean_i);
      h_mean_reco_bin_dm->SetBinError(i+1, meanerr_i);
      h_rms_reco_bin_dm->SetBinContent(i+1, rms_i);
      h_rms_reco_bin_dm->SetBinError(i+1, rmserr_i);
    }    
    fout->cd();
    h_mean_reco_bin_dm->Write();
    h_rms_reco_bin_dm->Write();
  }
  
  sw.Stop();

  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;

  fout->Close(); 
  //for(auto r : rans) delete r;

  return 1;
}
