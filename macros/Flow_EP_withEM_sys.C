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

double* Event_Plane_withEM(int, int, int, int, double, double, TFile*);

void load_Library()
{
    gSystem->CompileMacro("FlowAnalysis_Fitting.cxx");
    gSystem->Load("FlowAnalysis_Fitting_cxx.so"); // Linux

    gSystem->CompileMacro("FlowAnalysis_Helper.cxx");
    gSystem->Load("FlowAnalysis_Helper_cxx.so");
}


vector<double> Bin_pt_mass_bins[] = {
    {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, 
    {5, 6}, {6, 8}, {8, 10}, {10, 12}, {12, 15}, {15, 20}
};
double ptBins[] = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15,20};

vector<double> Bin_PoolEM = {0.,  5.,  10., 20., 30., 40.,
                             50., 60., 70., 80., 90.};
//Set Centrality
double cent_min = 10.;
double cent_max = 80.;
//
double chi2max_mass = 2.;
double chi2max_v2 = 2.;
//
bool meanPt = false;
bool rnorm_hist = false;
bool mass_hist = false;
bool yield_sys = false;
bool v2_sys = true;

// Define the pool for systematics:
// combinations
  double mass_min_sys[] = {2.2, 2.3, 2.4};
  double mass_max_sys[] = {4.2, 4.3, 4.4};
  string sig_enum[] = {"CB2(data)", "CB2(MC)", "NA60", "Chebychev",
                        "EventMixing"};
  string bkg_v2_enum[] = {"EventMixing(beta fix)", "EventMixing(beta free)"};
  int sig_mass[] = {0,1,2}; // CB2(data,MC),NA60
  int bkg_mass[] = {3,4};    // Chebychev, Event-Mixing
  int bkg_v2[] = {0,1};      // Event-Mixing, beta fix or free
  int nb_trials = int(size(sig_mass)) * int(size(bkg_mass)) *
                  int(size(bkg_v2)) * int(size(mass_min_sys));

int iSigFunc, iBkgFunc, iRange, iBkgV2Func, iPt;
int nPtBins = 2;
int nTrials = 0;

string dirPath = "/home/dhananjaya/Desktop/v2_Flow_Mixing_substractiion";

void SetLegend(TLegend *);
using namespace std;


void Flow_EP_withEM_sys()
{
    
    load_Library();

  // Histograms for plotting the results
  TH1D *histStatJpsiV2 = new TH1D("histStatJpsiV2", "", 10, ptBins);
  histStatJpsiV2 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
  histStatJpsiV2 -> GetYaxis() -> SetRangeUser(-0.05, 0.4);
  histStatJpsiV2 -> GetYaxis() -> SetTitle("#it{v}_{2}");
  histStatJpsiV2 -> SetMarkerStyle(20);
  histStatJpsiV2 -> SetMarkerColor(kRed+1);
  histStatJpsiV2 -> SetLineColor(kRed+1);

  TH1D *histSystJpsiV2 = new TH1D("histSystJpsiV2", "", 10, ptBins);
  histSystJpsiV2 -> SetMarkerStyle(20);
  histSystJpsiV2 -> SetMarkerColor(kRed+1);
  histSystJpsiV2 -> SetLineColor(kRed+1);
  histSystJpsiV2 -> SetFillStyle(0);

    for (iPt = 0;iPt < nPtBins;iPt++)
    {
    vector<double> jpsiV2s;
    vector<double> errJpsiV2s;
    vector<double> chi2Ndfs;
    vector<string> trialSetNames;


    const vector<double>& current_bin = Bin_pt_mass_bins[iPt];
    double minPtBin = current_bin.front();
    double maxPtBin = current_bin.back();

   TFile *fOut = new TFile(Form("systematics/fitResults_Pt_%0.1f_%0.1f.root", minPtBin, maxPtBin), "RECREATE");
   nTrials = 0;

    for (iSigFunc = 0;iSigFunc < int(size(sig_mass));iSigFunc++) { // loop on mass signal funcs
        
        for (iBkgFunc = 0; iBkgFunc < int(size(bkg_mass)) ;iBkgFunc++) { // loop on mass background funcs

            for (iRange = 0;iRange < int(size(mass_min_sys));iRange++) {   //Mass Range

            for (iBkgV2Func = 0;iBkgV2Func < int(size(bkg_v2));iBkgV2Func++) { // loop on v2 background types

        double *results; 
        results = Event_Plane_withEM(iPt,sig_mass[iSigFunc],bkg_mass[iBkgFunc], bkg_v2[iBkgV2Func],mass_min_sys[iRange], mass_max_sys[iRange],fOut);
        //gROOT->ProcessLine("Event_Plane_withEM(iPt+1,sig_mass[iSigFunc],bkg_mass[iBkgFunc], bkg_v2[iBkgV2Func],mass_min_sys[iRange], mass_max_sys[iRange],fOut)");
        
            
              //v2 from the fit 
              jpsiV2s.push_back(results[0]);
              errJpsiV2s.push_back(results[1]);
              chi2Ndfs.push_back(results[2]);
              trialSetNames.push_back(Form("%s + %s, %2.1f - %2.1f, %s[v2 bkg]", sig_enum[sig_mass[iSigFunc]].c_str(), sig_enum[bkg_mass[iBkgFunc]].c_str(), mass_min_sys[iRange], mass_max_sys[iRange], bkg_v2_enum[iBkgV2Func].c_str()));
              // Update the number of trials
              nTrials++;
          
    
          } // loop on v2 background types
        } //Mass Range
    } // loop on mass background funcs
}  //loop on mass signal funcs



    double sumJpsiV2s = 0; 
    double statErrJpsiV2s = 0;

    TH1D *histSyst = new TH1D("histSyst", "", nTrials, 0, nTrials);
    histSyst -> SetMarkerStyle(20);
    histSyst -> SetMarkerColor(kBlack);
    histSyst -> SetLineColor(kBlack);

    TH1D *histChi2Ndf = new TH1D("histChi2Ndf", "", nTrials, 0, nTrials);
    histChi2Ndf -> SetMarkerStyle(20);
    histChi2Ndf -> SetMarkerColor(kBlack);
    histChi2Ndf -> SetLineColor(kBlack);

    histSyst -> GetYaxis() -> SetRangeUser(0, 0.2);
    for (int iBin = 0;iBin < nTrials;iBin++) {
      //cout << iBin << ") " << jpsiV2s[iBin] << " +/- " << errJpsiV2s[iBin] << " -> " << trialSetNames[iBin].c_str() << endl;
      histSyst -> SetBinContent(iBin+1, jpsiV2s[iBin]);
      histSyst -> SetBinError(iBin+1, errJpsiV2s[iBin]);
      histSyst -> GetXaxis() -> SetBinLabel(iBin+1, trialSetNames[iBin].c_str());

      sumJpsiV2s = sumJpsiV2s + jpsiV2s[iBin];
	    statErrJpsiV2s = statErrJpsiV2s + errJpsiV2s[iBin];

      histChi2Ndf -> SetBinContent(iBin+1, chi2Ndfs[iBin]);
      histChi2Ndf -> GetXaxis() -> SetBinLabel(iBin+1, trialSetNames[iBin].c_str());
    }
    double meanJpsiV2 = sumJpsiV2s / nTrials;
    double meanStatErrJpsiV2 = statErrJpsiV2s / nTrials;

    double sumSyst = 0;
    for (int iBin = 0;iBin < nTrials;iBin++) {
      double dev = (jpsiV2s[iBin] - meanJpsiV2);
	    sumSyst = sumSyst + dev * dev;
    }
    double meanSystErrJpsiV2 = TMath::Sqrt(sumSyst / nTrials);

    histStatJpsiV2 -> SetBinContent(iPt+1, meanJpsiV2);
    histStatJpsiV2 -> SetBinError(iPt+1, meanStatErrJpsiV2);
    histSystJpsiV2 -> SetBinContent(iPt+1, meanJpsiV2);
    histSystJpsiV2 -> SetBinError(iPt+1, meanSystErrJpsiV2);

    double meanStatErrJpsiV2Perc = (meanStatErrJpsiV2 / meanJpsiV2) * 100;
    double meanSystErrJpsiV2Perc = (meanSystErrJpsiV2 / meanJpsiV2) * 100;

    TLine *lineMean = new TLine(0, meanJpsiV2, nTrials, meanJpsiV2);
    lineMean -> SetLineStyle(1);
    lineMean -> SetLineColor(kBlue);
    lineMean -> SetLineWidth(2);

    TLine *lineSystUp = new TLine(0, meanJpsiV2 + meanSystErrJpsiV2, nTrials, meanJpsiV2 + meanSystErrJpsiV2);
    lineSystUp-> SetLineStyle(2);
    lineSystUp-> SetLineColor(kBlue);
    lineSystUp-> SetLineWidth(2);

    TLine *lineSystLw = new TLine(0, meanJpsiV2 - meanSystErrJpsiV2, nTrials, meanJpsiV2 - meanSystErrJpsiV2);
    lineSystLw -> SetLineStyle(2);
    lineSystLw -> SetLineColor(kBlue);
    lineSystLw -> SetLineWidth(2);

    TLine *lineStatUp = new TLine(0, meanJpsiV2 + meanStatErrJpsiV2, nTrials, meanJpsiV2 + meanStatErrJpsiV2);
    lineStatUp-> SetLineStyle(2);
    lineStatUp-> SetLineColor(kGray+1);
    lineStatUp-> SetLineWidth(2);

    TLine *lineStatLw = new TLine(0, meanJpsiV2 - meanStatErrJpsiV2, nTrials, meanJpsiV2 - meanStatErrJpsiV2);
    lineStatLw -> SetLineStyle(2);
    lineStatLw -> SetLineColor(kGray+1);
    lineStatLw -> SetLineWidth(2);

    TCanvas *canvasSyst = new TCanvas("canvasSyst", "", 800, 900);
    canvasSyst -> SetTopMargin(0.05);
    canvasSyst -> SetBottomMargin(0.5);
    histSyst -> SetStats(0);
    histSyst -> GetXaxis() -> LabelsOption("v");
    if (meanJpsiV2 > 0) {
      histSyst -> GetYaxis() -> SetRangeUser(meanJpsiV2 -  (meanJpsiV2 / 2.), meanJpsiV2 + (meanJpsiV2 / 2.));
    } else {
      histSyst -> GetYaxis() -> SetRangeUser(-0.05, 0.05);
    }
    
    histSyst -> Draw();
    lineMean -> Draw();
    lineSystUp-> Draw();
    lineSystLw -> Draw();
    lineStatUp-> Draw();
    lineStatLw -> Draw();

    //TPaveText *display1 = new TPaveText(0.22, 0.80, 0.80, 0.85, "blNDC");
    TPaveText *display1 = new TPaveText(0.20, 0.80, 0.80, 0.85, "blNDC");
    display1 -> SetTextFont(42);
    display1 -> SetTextSize(0.036);
    display1 -> SetTextColor(kBlack);
    display1 -> SetBorderSize(0);
    display1 -> SetFillColor(0);

    //TText *text1 = display1 -> AddText(Form(" v_{2}^{J/#psi} [%1.0f < #it{p}_{T} < %1.0f GeV/#it{c}] = %5.4f #pm %5.4f (%4.3f %%) #pm %5.4f (%4.3f %%)", minPtBin, maxPtBin, meanJpsiV2, meanStatErrJpsiV2, meanStatErrJpsiV2Perc, meanSystErrJpsiV2, meanSystErrJpsiV2Perc));
    TText *text1 = display1 -> AddText(Form(" v_{2}^{J/#psi} [%1.0f < #it{p}_{T} < %1.0f GeV/#it{c}] = %5.4f #pm %5.4f #pm %5.4f", minPtBin, maxPtBin, meanJpsiV2, meanStatErrJpsiV2, meanSystErrJpsiV2));
    display1 -> Draw("same");

    canvasSyst -> Write();
    histChi2Ndf -> Write();
    histSyst  -> Write();
    canvasSyst -> SaveAs(Form("%s/systematic_plot_Pt_%1.0f_%1.0f.pdf", dirPath.c_str(), minPtBin, maxPtBin));

    fOut -> Close();
  } //Loop on Pt-binbs

  TCanvas *canvasJpsiV2 = new TCanvas("canvasJpsiV2", "", 800, 600);
  histStatJpsiV2 -> Draw("EP");
  histSystJpsiV2 -> Draw("E2P SAME");

  cout << "-------------------------" << endl;
  cout << "x_min x_max val stat syst" << endl;
  for (int iPt = 0;iPt < nPtBins;iPt++) {
    Printf("%d %6.5f %6.5f %6.5f ", iPt, histStatJpsiV2 -> GetBinContent(iPt+1), histStatJpsiV2 -> GetBinError(iPt+1), histSystJpsiV2 -> GetBinError(iPt+1));
  }




}


double* Event_Plane_withEM(int flag_binning, int flag_sig, int flag_bkg, int flag_v2, double mass_min, double mass_max, TFile *f)
{
  LOG(info) << "Start flow analysis...";
  double* results = new double[3];
  std::string FileName = "AnalysisResults.root";
  std::string inputFlag = "goodmedium";
  std::string muonCut = "muonLowPt510SigmaPDCA";
  std::string dimuonCut = "";
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
    case 0: case 1: case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 9:
    case 10: case 11: case 12: case 13: case 14: case 15: case 16: case 17:
        Bin_pt_mass = Bin_pt_mass_bins[flag_binning];
        break;
    default:
        Bin_pt_mass = Bin_pt_mass_bins[0]; // Default case
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
          
  /*                
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
   */

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
        0., 20.,mass_min, mass_max, Bin_PoolEM[i], Bin_PoolEM[i + 1], tp_V2SEPP,
        tp_V2SEMM, tp_V2MEPM, hist_rfactor);
    LOG(info) << Form("F factor for [%g - %g] (%%): %g", Bin_PoolEM[i],
                      Bin_PoolEM[i + 1], F_value);
    ffactor.emplace_back(F_value);
   
    
    TList *l_rfactor = new TList();
    l_rfactor->Add(hist_rfactor);
    
    if(rnorm_hist) 
    {
    f->cd();
    l_rfactor->Write(
    Form("RNormalizationFactor_%g_%g", Bin_PoolEM[i], Bin_PoolEM[i + 1]),
    TObject::kSingleKey);
    }

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
          Bin_pt_mass[i], Bin_pt_mass[i + 1], mass_min, mass_max, Bin_PoolEM[j],
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

if(v2_sys)
{
results[0] = results_v2[0];
results[1] = results_v2[1];
results[2] = results_v2[6];
}

if(yield_sys)
{
results[0] = results_v2[2];
results[1] = results_v2[3];
results[2] = results_v2[5];
}


  if(mass_hist)
   {
    f->cd();
    l_SE_ME->SetOwner();
    l_SE_ME->Write(Form("Mass_SEME_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
                   TObject::kSingleKey);

    f->cd();
    l_SE_ME_V2->SetOwner();
    l_SE_ME_V2->Write(Form("V2_SEME_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1]),
                      TObject::kSingleKey);
   }
   delete l_SE_ME;
   delete l_SE_ME_V2;


    f->cd();
    l_diff_fit->SetOwner();
    l_diff_fit->Write(
        Form("DifferentialFlow_Fit_%g_%g__%s_%s_%s__R_%g_%g", Bin_pt_mass[i], Bin_pt_mass[i + 1], sig_enum[flag_sig].data(), sig_enum[flag_bkg].data(), bkg_v2_enum[flag_v2].data(),mass_min,mass_max),
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

    

return results;

}
 
