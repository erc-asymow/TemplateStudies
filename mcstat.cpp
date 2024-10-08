#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TVector.h"
#include "TVectorT.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include <TMatrixD.h>
#include <TStopwatch.h>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <boost/program_options.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;
using namespace ROOT;

using namespace boost::program_options;


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
	("nevents",     value<long>()->default_value(2000000), "number of events")
	("ntoys",       value<long>()->default_value(3000), "number of toys")
	("tag",         value<std::string>()->default_value("closure"), "run type")
	("do_fit",      bool_switch()->default_value(false), "do_fit")
	("nbins",       value<int>()->default_value(200), "nbins")
	("nsigmas",     value<int>()->default_value(5), "nsigmas")
	("frac",       value<float>()->default_value(0.5), "frac")
	("asym",       value<float>()->default_value(0.015), "asym")
	("lumiscale",       value<float>()->default_value(1.0), "lumiscale")
	("seed",        value<int>()->default_value(4357), "seed");

      store(parse_command_line(argc, argv, desc), vm);
      notify(vm);
      if (vm.count("help")){
	std::cout << desc << '\n';
	return 0;
      }
      if (vm.count("nevents"))    std::cout << "Number of events: " << vm["nevents"].as<long>() << '\n';
      if (vm.count("tag"))        std::cout << "Tag: " << vm["tag"].as<std::string>() << '\n';
    }
  catch (const error &ex)
    {
      std::cerr << ex.what() << '\n';
    }

  long nevents    = vm["nevents"].as<long>();
  long ntoys    = vm["ntoys"].as<long>();
  std::string tag = vm["tag"].as<std::string>();
  int seed        = vm["seed"].as<int>();
  int nbins       = vm["nbins"].as<int>();
  int nsigmas     = vm["nsigmas"].as<int>();
  float frac      = vm["frac"].as<float>();
  float asym      = vm["asym"].as<float>();
  float lumiscale    = vm["lumiscale"].as<float>();
  bool do_fit     = vm["do_fit"].as<bool>();

  std::vector<TRandom3*> rans = {};
  for(unsigned int i = 0; i < 10; i++){
    rans.emplace_back( new TRandom3(seed + i*10) );
  }

  TFile* fout = TFile::Open(("root/mcstat_"+tag+".root").c_str(), "RECREATE");
  TTree* tree = new TTree("tree", "");
  double mu_data, err_data;
  double mu_data5s, err_data5s;
  double mu_mc, err_mc;
  double mu_data_BB, err_data_BB;
  double mu_data5s_BB, err_data5s_BB;
  double mu_true, err_true;
  tree->Branch("mu_data",  &mu_data, "mu_data/D");
  tree->Branch("err_data", &err_data, "err_data/D");
  tree->Branch("mu_data5s",  &mu_data5s, "mu_data5s/D");
  tree->Branch("err_data5s", &err_data5s, "err_data5s/D");
  tree->Branch("mu_data_BB",  &mu_data_BB, "mu_data_BB/D");
  tree->Branch("err_data_BB", &err_data_BB, "err_data_BB/D");
  tree->Branch("mu_data5s_BB",  &mu_data5s_BB, "mu_data5s_BB/D");
  tree->Branch("err_data5s_BB", &err_data5s_BB, "err_data5s_BB/D");
  tree->Branch("mu_mc",  &mu_mc, "mu_mc/D");
  tree->Branch("err_mc", &err_mc, "err_mc/D");
  tree->Branch("mu_true",  &mu_true, "mu_true/D");
  tree->Branch("err_true", &err_true, "err_true/D");

   
  TH1D* h_true_0 = new TH1D("h_true_0", "", nbins, 0.0, 1.0);
  TH1D* h_true_1 = new TH1D("h_true_1", "", nbins, 0.0, 1.0);
  for(int ib = 1; ib<=nbins; ib++){
    double ival_0 = nevents/nbins;
    ival_0  *= (1-frac)*(1-asym);
    double ival_1 = ib<=nbins/2 ? nevents/nbins*(1+asym) : nevents/nbins*(1-asym);
    ival_1  *= (frac);
    h_true_0->SetBinContent(ib, ival_0);
    h_true_1->SetBinContent(ib, ival_1);
  }

  double rescale = nevents/(h_true_0->Integral() + h_true_1->Integral());
  h_true_0->Scale(rescale);
  h_true_1->Scale(rescale);
  
  TH1D* h_true_tot = (TH1D*)h_true_0->Clone("h_true_tot");
  h_true_tot->Add(h_true_1, 1.0);
  double proc_ratio = h_true_0->Integral()/h_true_1->Integral();
  
  cout << "h_true_0: integral = " << h_true_0->Integral() << endl;
  cout << "h_true_1: integral = " << h_true_1->Integral() << endl;
  cout << "h_true_tot: integral = " << h_true_tot->Integral() << endl;

  MatrixXd A_true = MatrixXd::Zero(nbins, 2);
  MatrixXd invV_true = MatrixXd::Zero(nbins, nbins);

  for(unsigned int ir=0; ir<nbins; ir++){
    A_true(ir, 0) = h_true_0->GetBinContent(ir+1);
    A_true(ir, 1) = h_true_1->GetBinContent(ir+1);
  }
  for(unsigned int ir=0; ir<nbins; ir++){
    invV_true(ir,ir) = 1./h_true_tot->GetBinContent(ir+1);
  }

  MatrixXd C_true = ( A_true.transpose()*invV_true*A_true ).inverse();
  //cout << "True errors on mu_[0,1] = [" << TMath::Sqrt(C_true(0,0)) << "," << TMath::Sqrt(C_true(1,1)) << "]" << endl;

  mu_true = 0.0;
  err_true = TMath::Sqrt(C_true(0,0));
  
  /*
  MatrixXd C_inv = MatrixXd::Zero(2, 2);
  for(unsigned int ib=0; ib<nbins; ib++){
    C_inv += (A_true.row(ib).transpose()*A_true.row(ib))/h_true_tot->GetBinContent(ib+1);
  }
  MatrixXd C_alt = C_inv.inverse();
  cout << "Errors on mu_[0,1] = [" << TMath::Sqrt(C_alt(0,0)) << "," << TMath::Sqrt(C_alt(1,1)) << "]" << endl;
  */

  MatrixXd A_itoy = MatrixXd::Zero(nbins, 2);
  MatrixXd invV_itoy      = MatrixXd::Zero(nbins, nbins);
  MatrixXd invV5s_itoy    = MatrixXd::Zero(nbins, nbins);
  MatrixXd invV_BB_itoy   = MatrixXd::Zero(nbins, nbins);
  MatrixXd invV5s_BB_itoy = MatrixXd::Zero(nbins, nbins);
  MatrixXd invVMC_itoy    = MatrixXd::Zero(nbins, nbins);
  VectorXd y_itoy(nbins);
  VectorXd y5s_itoy(nbins);
  VectorXd yMC_itoy(nbins);
  VectorXd ynom_itoy(nbins);
  
  double prob_data = 0.;
  double prob_data5s = 0.;
  double prob_data_BB = 0.;
  double prob_data5s_BB = 0.;
  double prob_mc = 0.;

  for(unsigned int itoy=0; itoy<ntoys; itoy++){
    if(itoy%10000==0) cout << "Doing toy " << itoy << " / " << ntoys << endl;
      
    MatrixXd invC_itoy      = MatrixXd::Zero(2, 2);
    MatrixXd invC5s_itoy    = MatrixXd::Zero(2, 2);
    MatrixXd invCMC_itoy    = MatrixXd::Zero(2, 2);
    MatrixXd invC_BB_itoy   = MatrixXd::Zero(2, 2);
    MatrixXd invC5s_BB_itoy = MatrixXd::Zero(2, 2);
    MatrixXd invCMC_BB_itoy = MatrixXd::Zero(2, 2);
    
    for(unsigned int ir=0; ir<nbins; ir++){      
      A_itoy(ir,0) = rans[0]->Poisson( A_true(ir, 0)*lumiscale ) / lumiscale;
      A_itoy(ir,1) = rans[0]->Poisson( A_true(ir, 1)*lumiscale ) / lumiscale;

      // "MC" template drawn from true nominal, including lumi scale
      ynom_itoy(ir) = A_itoy.row(ir).sum();
      
      // toy data drawn from true nominal 
      y_itoy(ir)   = rans[0]->Poisson( A_true.row(ir).sum() );

      // toy data drawn from true 5sigma 
      y5s_itoy(ir) = rans[0]->Poisson( A_true(ir,0)*(1+err_true*nsigmas) + A_true(ir,1)  );

      // toy data drawn from "MC" 
      yMC_itoy(ir) = rans[0]->Poisson( A_itoy.row(ir).sum() );

      invV_itoy(ir,ir)      = 1./y_itoy(ir);
      invV5s_itoy(ir,ir)    = 1./y5s_itoy(ir);
      invV_BB_itoy(ir,ir)   = 1./(y_itoy(ir)   + ynom_itoy(ir)/lumiscale);
      invV5s_BB_itoy(ir,ir) = 1./(y5s_itoy(ir) + ynom_itoy(ir)/lumiscale);
      invVMC_itoy(ir,ir)    = 1./yMC_itoy(ir);

      // this is a common piece
      MatrixXd K_ir_itoy = A_itoy.row(ir).transpose()*A_itoy.row(ir);      
      invC_itoy      += K_ir_itoy/y_itoy(ir);
      invC5s_itoy    += K_ir_itoy/y5s_itoy(ir);
      invCMC_itoy    += K_ir_itoy/yMC_itoy(ir);
      invC_BB_itoy   += K_ir_itoy/(y_itoy(ir)   + ynom_itoy(ir)/lumiscale );
      invC5s_BB_itoy += K_ir_itoy/(y5s_itoy(ir) + ynom_itoy(ir)/lumiscale );
    }

    MatrixXd C_itoy      = invC_itoy.inverse();
    MatrixXd C5s_itoy    = invC5s_itoy.inverse();
    MatrixXd C_BB_itoy   = invC_BB_itoy.inverse();
    MatrixXd C5s_BB_itoy = invC5s_BB_itoy.inverse();
    MatrixXd CMC_itoy    = invCMC_itoy.inverse();
    VectorXd x_itoy      = C_itoy*A_itoy.transpose()*invV_itoy*(y_itoy - ynom_itoy);
    VectorXd x5s_itoy    = C5s_itoy*A_itoy.transpose()*invV5s_itoy*(y5s_itoy - ynom_itoy);
    VectorXd xMC_itoy    = CMC_itoy*A_itoy.transpose()*invVMC_itoy*(yMC_itoy - ynom_itoy);
    VectorXd x_BB_itoy   = C_BB_itoy*A_itoy.transpose()*invV_BB_itoy*(y_itoy - ynom_itoy);
    VectorXd x5s_BB_itoy = C5s_BB_itoy*A_itoy.transpose()*invV5s_BB_itoy*(y5s_itoy - ynom_itoy);
    mu_data      = x_itoy(0);
    mu_data5s    = x5s_itoy(0);
    mu_data_BB   = x_BB_itoy(0);
    mu_data5s_BB = x5s_BB_itoy(0);
    mu_mc = xMC_itoy(0); 
    err_data      = TMath::Sqrt(C_itoy(0,0));
    err_data5s    = TMath::Sqrt(C5s_itoy(0,0));
    err_data_BB   = TMath::Sqrt(C_BB_itoy(0,0));
    err_data5s_BB = TMath::Sqrt(C5s_BB_itoy(0,0));
    err_mc        = TMath::Sqrt(CMC_itoy(0,0));

    if( TMath::Abs(mu_data-mu_true)/err_data <= 1.0 ) prob_data += 1./ntoys;
    if( TMath::Abs(mu_data5s-(mu_true + err_true*nsigmas))/err_data5s <= 1.0 ) prob_data5s += 1./ntoys;
    if( TMath::Abs(mu_data_BB-mu_true)/err_data_BB <= 1.0 ) prob_data_BB += 1./ntoys;
    if( TMath::Abs(mu_data5s_BB-(mu_true + err_true*nsigmas))/err_data5s_BB <= 1.0 ) prob_data5s_BB += 1./ntoys;
    if( TMath::Abs(mu_mc-mu_true)/err_mc <= 1.0 ) prob_mc += 1./ntoys;
    
    tree->Fill();
  }

  cout << "Asympt. err:" << 1.0 - TMath::Prob(1.0, 1) << " (err=" << TMath::Sqrt(C_true(0,0)) << ")" << endl;
  
  TH1D* haux = new TH1D("haux","", 500, 0., 0.5);

  tree->Draw("err_mc>>haux");  
  cout << "MC:         " << prob_mc << " +/- " << TMath::Sqrt(prob_mc*(1-prob_mc)/ntoys) <<  " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ")" << endl;
  haux->Reset();

  tree->Draw("err_data>>haux");  
  cout << "Data:       " << prob_data << " +/- " << TMath::Sqrt(prob_data*(1-prob_data)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ")" << endl;
  haux->Reset();

  tree->Draw("err_data_BB>>haux");  
  cout << "Data BB:    " << prob_data_BB << " +/- " << TMath::Sqrt(prob_data_BB*(1-prob_data_BB)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ")" << endl;
  haux->Reset();

  tree->Draw("err_data5s>>haux");  
  cout << "Data 5s:    " << prob_data5s << " +/- " << TMath::Sqrt(prob_data5s*(1-prob_data5s)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ")" << endl;
  haux->Reset();

  tree->Draw("err_data5s_BB>>haux");  
  cout << "Data 5s BB: " << prob_data5s_BB << " +/- " << TMath::Sqrt(prob_data5s_BB*(1-prob_data5s_BB)/ntoys) << " (err=" << haux->GetMean() << " +/- " << haux->GetMeanError() << ")" << endl;
  haux->Reset();
    
  delete haux;
  
  fout->cd();
  h_true_0->Write();
  h_true_1->Write();
  tree->Write();
  
  sw.Stop();

  std::cout << "Real time: " << sw.RealTime() << " seconds " << "(CPU time:  " << sw.CpuTime() << " seconds)" << std::endl;

  
  fout->Close(); 
  for(auto r : rans) delete r;

  return 1;

}
