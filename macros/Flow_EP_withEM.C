#include "TAxis.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include <TCanvas.h>
#include <TChain.h>
#include <TClass.h>
#include <TCollection.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphMultiErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <THashList.h>
#include <THnSparse.h>
#include <TKey.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TMatrixD.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TParameter.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
namespace fs = std::filesystem;

#include "FlowAnalysis_Fitting.h"
#include "FlowAnalysis_Helper.h"
#include "/home/dhananjaya/alice/O2/Framework/Logger/include/Framework/Logger.h"


// Predefined binnings
vector<double> Bin_pt_mass_3 = {0, 3, 5, 15};
vector<double> Bin_pt_mass_3_bis = {1, 3, 5, 15};
vector<double> Bin_pt_mass_4 = {1, 2, 3, 5, 15};
vector<double> Bin_pt_mass_5 = {0, 2, 4, 6, 8, 15};
vector<double> Bin_pt_mass_8 = {0, 2, 3, 4, 5, 6, 8, 10, 15};
vector<double> Bin_pt_mass_10 = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 20};
vector<double> Bin_pt_mass_11 = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20};
vector<double> Bin_pt_mass_13 = {0., 0.5, 1., 1.5, 2.,  2.5, 3.,
                                 4., 5.,  6., 8.,  10., 12., 15.};
vector<double> Bin_pt_mass_15 = {0., 0.3, 1., 2.,  3.,  4.,  5.,  6.,
                                 7., 8.,  9., 10., 11., 12., 15., 20.};
vector<double> Bin_pt_mass_16 = {0.,  0.3, 1., 1.5, 2., 2.5, 3.,  3.5, 4.,
                                 4.5, 5.,  6., 7.,  8., 10., 12., 20.};
vector<double> Bin_pt_mass_18 = {0,   1, 1.5, 2, 2.5, 3,  3.5, 4,  4.5, 5,
                                 5.5, 6, 7,   8, 9,   10, 12,  15, 20};
vector<double> Bin_PoolEM = {0.,  5.,  10., 20., 30., 40.,
                             50., 60., 70., 80., 90.};

//root -l Event_Plane_withEM.C(2,6,0,0)
// NA60, Exp2, EventMixing(betafix), no run2 data

//string FlowAnalysis_Fitting::model_string[8] = {
//    "CB2(data)",   "CB2(MC)", "NA60", "Cheby",
///    "EventMixing", "VWG",     "Exp2", "PolExp"};
//string FlowAnalysis_Fitting::v2bkg_string[2] = {"EventMixing(beta fix)",
 //                                               "EventMixing(beta free)"};

 using namespace std;

void SetLegend(TLegend *);

void Event_Plane_withEM(int flag_binning, int flag_sig, int flag_bkg, int flag_v2, int flag_run2,
    int flag_run2yield, std::string FileName = "AnalysisResults.root",
    double mass_min = 2.3, double mass_max = 4.3, double cent_min = 10.,
    double cent_max = 50., double chi2max_mass = 2., double chi2max_v2 = 2.,
    double fmin = 1., double fmax = 5., bool sys = false, bool meanPt = true,
    bool SaveSys = false, std::string inputFlag = "goodmedium",
    std::string muonCut = "muonLowPt510SigmaPDCA", std::string dimuonCut = "") 

{
    gSystem->CompileMacro("FlowAnalysis_Fitting.cxx");
    gSystem->Load("FlowAnalysis_Fitting_cxx.so"); // Linux

    gSystem->CompileMacro("FlowAnalysis_Helper.cxx");
    gSystem->Load("FlowAnalysis_Helper_cxx.so");

  LOG(info) << "Start flow analysis...";
  // Init Helper class
  FlowAnalysis_Helper *helper = new FlowAnalysis_Helper();

  // Load input data for analysis
  TProfile3D *tp_V2SEPM;
  TProfile3D *tp_V2SEPP;
  TProfile3D *tp_V2SEMM;
  TProfile3D *tp_V2MEPM;
  TProfile3D *tp_V2MEPP;
  TProfile3D *tp_V2MEMM;

  helper->LoadDataMEProfile(FileName, tp_V2SEPM, tp_V2SEPP, tp_V2SEMM,tp_V2MEPM, tp_V2MEPP, tp_V2MEMM, muonCut,dimuonCut);
 std::cout<<" Get entries :::::::::: "<<tp_V2MEPM->GetEntries()<<endl;
  // Get general binning information
  TAxis *massAxis = new TAxis();
  TAxis *ptAxis = new TAxis();
  TAxis *centAxis = new TAxis();
  massAxis = tp_V2SEPM->GetXaxis();
  ptAxis = tp_V2SEPM->GetYaxis();
  centAxis = tp_V2SEPM->GetZaxis();
  int Nbins_mass = massAxis->GetNbins();
  int Nbins_pt = ptAxis->GetNbins();
  int Nbins_cent = centAxis->GetNbins();
  
  // Initialize fitter
  FlowAnalysis_Fitting *fitter = new FlowAnalysis_Fitting();
  fitter->init();
  fitter->setChi2MaxMass(chi2max_mass);
  fitter->setChi2MaxV2(chi2max_v2);
  fitter->setCentRange(cent_min, cent_max);

  // Define variables' range for analysis
  vector<double> Bin_pt_mass;
  switch (flag_binning) {
  case 3:
    Bin_pt_mass = Bin_pt_mass_3;
    break;
  case 30:
    Bin_pt_mass = Bin_pt_mass_3_bis;
    break;
  case 4:
    Bin_pt_mass = Bin_pt_mass_4;
    break;
  case 5:
    Bin_pt_mass = Bin_pt_mass_5;
    break;
  case 8:
    Bin_pt_mass = Bin_pt_mass_8;
    break;
  case 10:
    Bin_pt_mass = Bin_pt_mass_10;
    break;
  case 11:
    Bin_pt_mass = Bin_pt_mass_11;
    break;
  case 13:
    Bin_pt_mass = Bin_pt_mass_13;
    break;
  case 15:
    Bin_pt_mass = Bin_pt_mass_15;
    break;
  case 16:
    Bin_pt_mass = Bin_pt_mass_16;
    break;
  case 18:
    Bin_pt_mass = Bin_pt_mass_11;
    break;
  default:
    Bin_pt_mass = Bin_pt_mass_10;
    break;
  }

  // Define the pool for systematics:
  // combinations
  double mass_min_sys[3] = {2.2, 2.3, 2.4};
  double mass_max_sys[3] = {4.2, 4.3, 4.4};
  string sig_enum[5] = {"CB2(data)", "CB2(MC)", "NA60", "Chebychev",
                        "EventMixing"};
  string bkg_v2_enum[2] = {"EventMixing(beta fix)", "EventMixing(beta free)"};
  int sig_mass[3] = {0, 1, 2}; // CB2(data,MC),NA60
  int bkg_mass[2] = {3, 4};    // Chebychev, Event-Mixing
  int bkg_v2[2] = {0, 1};      // Event-Mixing, beta fix or free
  int nb_trials = int(size(sig_mass)) * int(size(bkg_mass)) *
                  int(size(bkg_v2)) * int(size(mass_min_sys));

  // Create output file
  TFile *f = new TFile(
      sys ? Form("FlowAnalysisResults_%s_"
                 "EventMixing%d_%s_%g_%"
                 "g_%dBinPt_%s_withSys.root",
                 inputFlag.c_str(), nb_trials, muonCut.c_str(), cent_min,
                 cent_max, int(Bin_pt_mass.size()) - 1,
                 meanPt ? "MeanPt" : "NoMeanPt")
          : Form("FlowAnalysisResults_%s_"
                 "EventMixing_%s_%g_%"
                 "g_%dBinPt_%s.root",
                 inputFlag.c_str(), muonCut.c_str(), cent_min, cent_max,
                 int(Bin_pt_mass.size()) - 1, meanPt ? "MeanPt" : "NoMeanPt"),
      "RECREATE");

///////////////////////////////////////////////////
  ///                                             ///
  ///   Analysis for Differential Flow of J/Psi   ///
  ///                                             ///
  ///////////////////////////////////////////////////

  LOG(info) << "Processing analysis for differential flow ...";

  // Get index of interested centrality bin
  int itmin = std::find(Bin_PoolEM.begin(), Bin_PoolEM.end(), cent_min) -
              Bin_PoolEM.begin();
  int itmax = std::find(Bin_PoolEM.begin(), Bin_PoolEM.end(), cent_max) -
              Bin_PoolEM.begin();

  // Calculate pT-integrated R factors and F factors
  LOG(info) << "Processing analysis for pT-integrated R and F factors ...";
  cout<<" int(Bin_PoolEM.size()) "<<int(Bin_PoolEM.size())<<endl;
  vector<double> ffactor;
  for (int i = 0; i < int(Bin_PoolEM.size()) - 1; i++) {
    TH1D *hist_rfactor = helper->GetRfactorProfile(
        0., 20., 0., 5., Bin_PoolEM[i], Bin_PoolEM[i + 1], tp_V2MEPM, tp_V2MEPP,
        tp_V2MEMM);
    double F_value = helper->GetFfactorProfile(
        0., 20., fmin, fmax, Bin_PoolEM[i], Bin_PoolEM[i + 1], tp_V2SEPP,
        tp_V2SEMM, tp_V2MEPM, hist_rfactor);
    LOG(info) << Form("F factor for [%g - %g] (%%): %g", Bin_PoolEM[i],
                      Bin_PoolEM[i + 1], F_value);
    ffactor.emplace_back(F_value);
    TList *l_rfactor = new TList();
    l_rfactor->Add(hist_rfactor);
    f->cd();
    l_rfactor->Write(
        Form("RNormalizationFactor_%g_%g", Bin_PoolEM[i], Bin_PoolEM[i + 1]),
        TObject::kSingleKey);

    delete hist_rfactor;
    delete l_rfactor;
  }

 LOG(info) << "Processing analysis for pT-differential R and F factors ...";
  vector<double *> ffactor_ptDiff;
  for (int i = 0; i < int(Bin_pt_mass.size()) - 1; i++) {
    ffactor_ptDiff.emplace_back(new double[itmax - itmin]);
  }
  for (int i = 0; i < int(Bin_pt_mass.size()) - 1; i++) {
    for (int j = itmin; j < itmax; j++) {
      TH1D *hist_rfactor = helper->GetRfactorProfile(
          Bin_pt_mass[i], Bin_pt_mass[i + 1], 0., 5., Bin_PoolEM[j],
          Bin_PoolEM[j + 1], tp_V2MEPM, tp_V2MEPP, tp_V2MEMM);
      double F_value = helper->GetFfactorProfile(
          Bin_pt_mass[i], Bin_pt_mass[i + 1], fmin, fmax, Bin_PoolEM[j],
          Bin_PoolEM[j + 1], tp_V2SEPP, tp_V2SEMM, tp_V2MEPM, hist_rfactor);
      LOG(info) << Form("F factor for [%g - %g] (%%) [%g - %g] (GeV/c): %g",
                        Bin_PoolEM[j], Bin_PoolEM[j + 1], Bin_pt_mass[i],
                        Bin_pt_mass[i + 1], F_value);
      ffactor_ptDiff[i][j - itmin] = F_value;
      delete hist_rfactor;
    }
  }

   // Create histogram for pt-differential v2
  double *SNR = new double[int(Bin_pt_mass.size()) - 1];
  
  for (int i = 0; i < int(Bin_pt_mass.size()) - 1; i++) {

    // Normalisation for mixed-event spectra
    TH1D *hs_mass_mepm_proj = new TH1D();
    TH1D *hs_mass_mepp_proj = new TH1D();
    TH1D *hs_mass_memm_proj = new TH1D();
    for (int j = itmin; j < itmax; j++) {
      if (j == itmin) {
        hs_mass_mepm_proj = helper->GetMassProfile(
            Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
            Bin_PoolEM[j], Bin_PoolEM[j + 1], tp_V2MEPM, "MEPM");
        hs_mass_mepm_proj->Scale(ffactor_ptDiff[i][j - itmin]);

        hs_mass_mepp_proj = helper->GetMassProfile(
            Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
            Bin_PoolEM[j], Bin_PoolEM[j + 1], tp_V2MEPP, "MEPP");
        hs_mass_mepp_proj->Scale(ffactor_ptDiff[i][j - itmin]);

        hs_mass_memm_proj = helper->GetMassProfile(
            Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max,
            Bin_PoolEM[j], Bin_PoolEM[j + 1], tp_V2MEMM, "MEMM");
        hs_mass_memm_proj->Scale(ffactor_ptDiff[i][j - itmin]);

      } else {
        hs_mass_mepm_proj->Add(
            helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                                   mass_max, Bin_PoolEM[j], Bin_PoolEM[j + 1],
                                   tp_V2MEPM, "MEPM"),
            ffactor_ptDiff[i][j - itmin]);

        hs_mass_mepp_proj->Add(
            helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                                   mass_max, Bin_PoolEM[j], Bin_PoolEM[j + 1],
                                   tp_V2MEPP, "MEPP"),
            ffactor_ptDiff[i][j - itmin]);

        hs_mass_memm_proj->Add(
            helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                                   mass_max, Bin_PoolEM[j], Bin_PoolEM[j + 1],
                                   tp_V2MEMM, "MEMM"),
            ffactor_ptDiff[i][j - itmin]);
      }
    }

    // Same-event profiles: mass
    TH1D *hs_mass_sepm_proj =
        helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                               mass_max, cent_min, cent_max, tp_V2SEPM, "SEPM");
    TH1D *hs_mass_sepp_proj =
        helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                               mass_max, cent_min, cent_max, tp_V2SEPP, "SEPP");
    TH1D *hs_mass_semm_proj =
        helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                               mass_max, cent_min, cent_max, tp_V2SEMM, "SEMM");

    // Same-event profiles: v2
    TH1D *hs_v2_sepm_proj =
        helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                             mass_max, cent_min, cent_max, tp_V2SEPM, "SEPM");
    TH1D *hs_v2_sepp_proj =
        helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                             mass_max, cent_min, cent_max, tp_V2SEPP, "SEPP");
    TH1D *hs_v2_semm_proj =
        helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                             mass_max, cent_min, cent_max, tp_V2SEMM, "SEMM");

    // Mixed-event profiles: v2
    TH1D *hs_v2_mepm_proj =
        helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                             mass_max, cent_min, cent_max, tp_V2MEPM, "MEPM");
    TH1D *hs_v2_mepp_proj =
        helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                             mass_max, cent_min, cent_max, tp_V2MEPP, "MEPP");
    TH1D *hs_v2_memm_proj =
        helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                             mass_max, cent_min, cent_max, tp_V2MEMM, "MEMM");

    // Mean pT profile
    TH1D *hs_mass_sepm_proj_meanPt = new TH1D();
    if (meanPt) {
      hs_mass_sepm_proj_meanPt =
          helper->GetMeanPt(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                            mass_max, cent_min, cent_max, tp_V2SEPM, "SEPM");
    }

    // Save plots for invariant mass
    TList *l_SE_ME = new TList();
    helper->PlotSEME("PM", "Mass", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                     mass_max, cent_min, cent_max, hs_mass_sepm_proj,
                     hs_mass_mepm_proj, l_SE_ME);
    helper->PlotSEME("PP", "Mass", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                     mass_max, cent_min, cent_max, hs_mass_sepp_proj,
                     hs_mass_mepp_proj, l_SE_ME);
    helper->PlotSEME("MM", "Mass", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                     mass_max, cent_min, cent_max, hs_mass_semm_proj,
                     hs_mass_memm_proj, l_SE_ME);

    TList *l_SE_ME_V2 = new TList();
    helper->PlotSEME("PM", "V2", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                     mass_max, cent_min, cent_max, hs_v2_sepm_proj,
                     hs_v2_mepm_proj, l_SE_ME_V2);
    helper->PlotSEME("PP", "V2", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                     mass_max, cent_min, cent_max, hs_v2_sepp_proj,
                     hs_v2_mepp_proj, l_SE_ME_V2);
    helper->PlotSEME("MM", "V2", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min,
                     mass_max, cent_min, cent_max, hs_v2_semm_proj,
                     hs_v2_memm_proj, l_SE_ME_V2);

    /// Do fitting
    // Configuration for fitting
    fitter->setModel(flag_sig, flag_bkg);
    fitter->setModelV2(flag_v2);
    fitter->setMassRange(mass_min, mass_max);
    fitter->setPtRange(Bin_pt_mass[i], Bin_pt_mass[i + 1]);
    fitter->setOrder(2);
    //fitter->setMode(0); // standard mode, no systematics

    // Fit invariant mass + v2
    TList *l_diff_fit = new TList();
    vector<double> results_v2 =
        meanPt ? fitter->runFittingEM(hs_mass_sepm_proj, hs_mass_mepm_proj,
                                      hs_v2_sepm_proj, hs_v2_mepm_proj,
                                      hs_mass_sepm_proj_meanPt, l_diff_fit)
               : fitter->runFittingEM(hs_mass_sepm_proj, hs_mass_mepm_proj,
                                      hs_v2_sepm_proj, hs_v2_mepm_proj, nullptr,
                                      l_diff_fit);

    

   SNR[i] = results_v2[4];

    f->cd();
    l_SE_ME->SetOwner();
    l_SE_ME->Write(Form("Mass_SEME_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
                   TObject::kSingleKey);
    delete l_SE_ME;

    f->cd();
    l_SE_ME_V2->SetOwner();
    l_SE_ME_V2->Write(Form("V2_SEME_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
                      TObject::kSingleKey);
    delete l_SE_ME_V2;

    f->cd();
    l_diff_fit->SetOwner();
    l_diff_fit->Write(
        Form("DifferentialFlow_Fit_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
        TObject::kSingleKey);
    delete l_diff_fit;

    delete hs_mass_sepm_proj;
    delete hs_mass_sepp_proj;
    delete hs_mass_semm_proj;
    delete hs_mass_mepm_proj;
    delete hs_mass_mepp_proj;
    delete hs_mass_memm_proj;
    delete hs_v2_sepm_proj;
    delete hs_v2_sepp_proj;
    delete hs_v2_semm_proj;
    if (hs_v2_mepm_proj) {
      delete hs_v2_mepm_proj;
    }
    if (hs_v2_mepp_proj) {
      delete hs_v2_mepp_proj;
    }
    if (hs_v2_memm_proj) {
      delete hs_v2_memm_proj;
    }
    if (hs_mass_sepm_proj_meanPt) {
      delete hs_mass_sepm_proj_meanPt;
    }

  }


}
 
