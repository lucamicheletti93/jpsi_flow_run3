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

vector<double> Bin_pt_mass = {0, 2, 4, 6, 8};
vector<double> Bin_PoolEM = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};

using namespace std;

//string FlowAnalysis_Fitting::v2bkg_string[2] = {"EventMixing(beta fix)", "EventMixing(beta free)"};

void flow_fitter(int flag_sig, int flag_bkg, int flag_v2, 
    std::string FileName = "AnalysisResults_OO_merged.root",
    double mass_min = 2.4, double mass_max = 5., 
    double cent_min = 20., double cent_max = 60., 
    double chi2max_mass = 2., double chi2max_v2 = 2.,
    double fmin = 1., double fmax = 5., bool meanPt = false,
    std::string muonCut = "muonQualityCutsMUONStandalone", std::string dimuonCut = "") {
    // other cuts: muonLowPt210SigmaPDCA, muonQualityCutsMUONStandalone
    
    //gSystem->CompileMacro("FlowAnalysis_Fitting.cxx");
    gSystem->Load("FlowAnalysis_Fitting_cxx.so"); // Linux

    //gSystem->CompileMacro("FlowAnalysis_Helper.cxx");
    gSystem->Load("FlowAnalysis_Helper_cxx.so");

    // Init Helper class
    FlowAnalysis_Helper *helper = new FlowAnalysis_Helper();

    // Load input data for analysis
    TProfile3D *tp_V2SEPM;
    TProfile3D *tp_V2SEPP;
    TProfile3D *tp_V2SEMM;
    TProfile3D *tp_V2MEPM;
    TProfile3D *tp_V2MEPP;
    TProfile3D *tp_V2MEMM;
    helper->LoadDataMEProfile(FileName, tp_V2SEPM, tp_V2SEPP, tp_V2SEMM, tp_V2MEPM, tp_V2MEPP, tp_V2MEMM, muonCut, dimuonCut);

    /*TH3D *tp_V2_cp_projxyz = tp_V2SEPM->ProjectionXYZ("tp_V2_cp_projxyz","B");
    TH3D* rebinned = dynamic_cast<TH3D*>(tp_V2_cp_projxyz->RebinX(5, "Mass_Pt_CentFT0C_V2_rebinned"));
    TH1D *tp_V2_cp_proj = dynamic_cast<TH1D *>(rebinned->Project3D("x"));
    tp_V2_cp_proj->Draw();*/

    /*tp_V2SEPM->RebinX(2);
    tp_V2SEPP->RebinX(2);
    tp_V2SEMM->RebinX(2);
    tp_V2MEPM->RebinX(2);
    tp_V2MEPP->RebinX(2);
    tp_V2MEMM->RebinX(2);*/

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

    // Define the pool for systematics:
    // combinations
    double mass_min_sys[3] = {2.2, 2.3, 2.4};
    double mass_max_sys[3] = {4.2, 4.3, 4.4};
    string sig_enum[5] = {"CB2(data)", "CB2(MC)", "NA60", "Chebychev", "EventMixing"};
    string bkg_v2_enum[2] = {"EventMixing(beta fix)", "EventMixing(beta free)"};
    int sig_mass[3] = {0, 1, 2}; // CB2(MC,data) NA60
    int bkg_mass[2] = {3, 4};    // Chebychev Event-Mixing
    int bkg_v2[2] = {0, 1};      // Event-Mixing, beta fix or free
    int nb_trials = int(size(sig_mass)) * int(size(bkg_mass)) *int(size(bkg_v2)) * int(size(mass_min_sys));

    // Create output file
    TFile *fOut = new TFile("outout_test.root", "RECREATE");

    ///////////////////////////////////////////////////
    ///                                             ///
    ///   Analysis for Differential Flow of J/Psi   ///
    ///                                             ///
    ///////////////////////////////////////////////////

    // Get index of interested centrality bin
    int itmin = std::find(Bin_PoolEM.begin(), Bin_PoolEM.end(), cent_min) - Bin_PoolEM.begin();
    int itmax = std::find(Bin_PoolEM.begin(), Bin_PoolEM.end(), cent_max) - Bin_PoolEM.begin();

    // Calculate pT-integrated R factors and F factors
    vector<double> ffactor;
    for (int i = 0; i < int(Bin_PoolEM.size()) - 1; i++) {
        TH1D *hist_rfactor = helper->GetRfactorProfile(0., 20., 0., 5., Bin_PoolEM[i], Bin_PoolEM[i + 1], tp_V2MEPM, tp_V2MEPP, tp_V2MEMM);
        double F_value = helper->GetFfactorProfile(0., 20., fmin, fmax, Bin_PoolEM[i], Bin_PoolEM[i + 1], tp_V2SEPP, tp_V2SEMM, tp_V2MEPM, hist_rfactor);
        ffactor.emplace_back(F_value);
        TList *l_rfactor = new TList();
        l_rfactor->Add(hist_rfactor);
        fOut->cd();
        l_rfactor->Write(Form("RNormalizationFactor_%g_%g", Bin_PoolEM[i], Bin_PoolEM[i + 1]), TObject::kSingleKey);
        //delete hist_rfactor;
        //delete l_rfactor;
    }

    vector<double *> ffactor_ptDiff;
    for (int i = 0; i < int(Bin_pt_mass.size()) - 1; i++) {
        ffactor_ptDiff.emplace_back(new double[itmax - itmin]);
    }
    for (int i = 0; i < int(Bin_pt_mass.size()) - 1; i++) {
        for (int j = itmin; j < itmax; j++) {
            TH1D *hist_rfactor = helper->GetRfactorProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], 0., 5., Bin_PoolEM[j], Bin_PoolEM[j + 1], tp_V2MEPM, tp_V2MEPP, tp_V2MEMM);
            double F_value = helper->GetFfactorProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], fmin, fmax, Bin_PoolEM[j], Bin_PoolEM[j + 1], tp_V2SEPP, tp_V2SEMM, tp_V2MEPM, hist_rfactor);
            ffactor_ptDiff[i][j - itmin] = F_value;
            cout << "F-value: " << F_value << endl;
            delete hist_rfactor;
        }
    }

    // Create histogram for pt-differential v2
    double *x_yield = new double[int(Bin_pt_mass.size()) - 1];
    double *y_yield = new double[int(Bin_pt_mass.size()) - 1];
    double *ex_yield = new double[int(Bin_pt_mass.size()) - 1];
    double *ey_yield = new double[int(Bin_pt_mass.size()) - 1];
    double *eysys_yield = new double[int(Bin_pt_mass.size()) - 1];
    double *SNR = new double[int(Bin_pt_mass.size()) - 1];
    double *x_v2pt = new double[int(Bin_pt_mass.size()) - 1];
    double *y_v2pt = new double[int(Bin_pt_mass.size()) - 1];
    double *ex_v2pt = new double[int(Bin_pt_mass.size()) - 1];
    double *ey_v2pt = new double[int(Bin_pt_mass.size()) - 1];
    double *eysys_v2pt = new double[int(Bin_pt_mass.size()) - 1];

    // Initialize arrays for each trial of systematic study
    const int nbCombo_v2 = int(size(sig_mass)) * int(size(bkg_mass)) * int(size(bkg_v2)) * int(size(mass_max_sys));
    const int nbCombo_yield = int(size(sig_mass)) * int(size(bkg_mass)) * int(size(mass_max_sys));

    vector<double *> x_sys_pt, ex_sys_pt, y_sys_yield, ey_sys_yield, y_sys_v2, ey_sys_v2, chi2_yield, chi2_v2, chi2_meanPt;
    double *bins_sys_v2 = new double[nbCombo_v2 + 1];
    double *bins_sys_yield = new double[nbCombo_yield + 1];
    for (int i = 0; i < nbCombo_v2 + 1; i++) {bins_sys_v2[i] = i;}
    for (int i = 0; i < nbCombo_yield + 1; i++) {bins_sys_yield[i] = i;}

    vector<TH1D *> hist_sys_yield, hist_sys_v2, hist_sys_meanPt;
    for (int i = 0; i < int(Bin_pt_mass.size()) - 1; i++) {
        y_sys_yield.emplace_back(new double[nbCombo_yield]);
        ey_sys_yield.emplace_back(new double[nbCombo_yield]);

        y_sys_v2.emplace_back(new double[nbCombo_v2]);
        ey_sys_v2.emplace_back(new double[nbCombo_v2]);

        chi2_yield.emplace_back(new double[nbCombo_yield]);
        chi2_v2.emplace_back(new double[nbCombo_v2]);

        hist_sys_yield.emplace_back(new TH1D(Form("hist_sys_yield_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]), Form("hist_sys_yield_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]), nbCombo_yield, bins_sys_yield));
        hist_sys_v2.emplace_back( new TH1D(Form("hist_sys_v2_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]), Form("hist_sys_v2_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]), nbCombo_v2, bins_sys_v2));

        if (meanPt) {
            x_sys_pt.emplace_back(new double[nbCombo_yield]);
            ex_sys_pt.emplace_back(new double[nbCombo_yield]);
            chi2_meanPt.emplace_back(new double[nbCombo_yield]);
            hist_sys_meanPt.emplace_back(new TH1D(Form("hist_sys_meanPt_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]), Form("hist_sys_meanPt_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]), nbCombo_yield, bins_sys_yield));
        }
    }

    int rebin = 2;
    for (int i = 0; i < int(Bin_pt_mass.size()) - 1; i++) {
        // Normalisation for mixed-event spectra
        TH1D *hs_mass_mepm_proj = new TH1D();
        TH1D *hs_mass_mepp_proj = new TH1D();
        TH1D *hs_mass_memm_proj = new TH1D();
        for (int j = itmin; j < itmax; j++) {
            if (j == itmin) {
                hs_mass_mepm_proj = helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, Bin_PoolEM[j], Bin_PoolEM[j + 1], tp_V2MEPM, "MEPM", rebin);
                hs_mass_mepm_proj->Scale(ffactor_ptDiff[i][j - itmin]);

                hs_mass_mepp_proj = helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, Bin_PoolEM[j], Bin_PoolEM[j + 1], tp_V2MEPP, "MEPP", rebin);
                hs_mass_mepp_proj->Scale(ffactor_ptDiff[i][j - itmin]);

                hs_mass_memm_proj = helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, Bin_PoolEM[j], Bin_PoolEM[j + 1], tp_V2MEMM, "MEMM", rebin);
                hs_mass_memm_proj->Scale(ffactor_ptDiff[i][j - itmin]);
            } else {
                hs_mass_mepm_proj->Add(helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, Bin_PoolEM[j], Bin_PoolEM[j + 1], tp_V2MEPM, "MEPM", rebin), ffactor_ptDiff[i][j - itmin]);
                hs_mass_mepp_proj->Add(helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, Bin_PoolEM[j], Bin_PoolEM[j + 1], tp_V2MEPP, "MEPP", rebin), ffactor_ptDiff[i][j - itmin]);
                hs_mass_memm_proj->Add(helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, Bin_PoolEM[j], Bin_PoolEM[j + 1], tp_V2MEMM, "MEMM", rebin), ffactor_ptDiff[i][j - itmin]);
            }
        }

        // Same-event profiles: mass
        TH1D *hs_mass_sepm_proj = helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min, cent_max, tp_V2SEPM, "SEPM", rebin);
        TH1D *hs_mass_sepp_proj = helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min, cent_max, tp_V2SEPP, "SEPP", rebin);
        TH1D *hs_mass_semm_proj = helper->GetMassProfile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min, cent_max, tp_V2SEMM, "SEMM", rebin);

        // Same-event profiles: v2
        TH1D *hs_v2_sepm_proj = helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min, cent_max, tp_V2SEPM, "SEPM", rebin);
        TH1D *hs_v2_sepp_proj = helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min, cent_max, tp_V2SEPP, "SEPP", rebin);
        TH1D *hs_v2_semm_proj = helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min, cent_max, tp_V2SEMM, "SEMM", rebin);

        // Mixed-event profiles: v2
        TH1D *hs_v2_mepm_proj = helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min, cent_max, tp_V2MEPM, "MEPM", rebin);
        TH1D *hs_v2_mepp_proj = helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min, cent_max, tp_V2MEPP, "MEPP", rebin);
        TH1D *hs_v2_memm_proj = helper->GetV2Profile(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min, cent_max, tp_V2MEMM, "MEMM", rebin);

        // Mean pT profile
        TH1D *hs_mass_sepm_proj_meanPt = new TH1D();
        if (meanPt) {hs_mass_sepm_proj_meanPt = helper->GetMeanPt(Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min, cent_max, tp_V2SEPM, "SEPM");}

        // Save plots for invariant mass
        TList *l_SE_ME = new TList();
        helper->PlotSEME("PM", "Mass", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min, cent_max, hs_mass_sepm_proj, hs_mass_mepm_proj, l_SE_ME);
        helper->PlotSEME("PP", "Mass", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min, cent_max, hs_mass_sepp_proj, hs_mass_mepp_proj, l_SE_ME);
        helper->PlotSEME("MM", "Mass", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min, cent_max, hs_mass_semm_proj, hs_mass_memm_proj, l_SE_ME);

        TList *l_SE_ME_V2 = new TList();
        helper->PlotSEME("PM", "V2", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min, cent_max, hs_v2_sepm_proj, hs_v2_mepm_proj, l_SE_ME_V2);
        helper->PlotSEME("PP", "V2", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min, cent_max, hs_v2_sepp_proj, hs_v2_mepp_proj, l_SE_ME_V2);
        helper->PlotSEME("MM", "V2", Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, cent_min, cent_max, hs_v2_semm_proj, hs_v2_memm_proj, l_SE_ME_V2);

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
            meanPt ? fitter->runFittingEM(hs_mass_sepm_proj, hs_mass_mepm_proj, hs_v2_sepm_proj, hs_v2_mepm_proj, hs_mass_sepm_proj_meanPt, l_diff_fit, true)
                : fitter->runFittingEM(hs_mass_sepm_proj, hs_mass_mepm_proj, hs_v2_sepm_proj, hs_v2_mepm_proj, nullptr, l_diff_fit, true);

        SNR[i] = results_v2[4];
        fOut->cd();
        l_SE_ME->SetOwner();
        l_SE_ME->Write(Form("Mass_SEME_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]), TObject::kSingleKey);
        delete l_SE_ME;

        fOut->cd();
        l_SE_ME_V2->SetOwner();
        l_SE_ME_V2->Write(Form("V2_SEME_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]), TObject::kSingleKey);
        delete l_SE_ME_V2;

        fOut->cd();
        l_diff_fit->SetOwner();
        l_diff_fit->Write( Form("DifferentialFlow_Fit_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]), TObject::kSingleKey);
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
        if (hs_v2_mepm_proj) {delete hs_v2_mepm_proj;}
        if (hs_v2_mepp_proj) {delete hs_v2_mepp_proj;}
        if (hs_v2_memm_proj) {delete hs_v2_memm_proj;}
        if (hs_mass_sepm_proj_meanPt) {delete hs_mass_sepm_proj_meanPt;}
    }
}