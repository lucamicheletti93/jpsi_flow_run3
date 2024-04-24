#include "TProfile.h"
#include <fstream>
#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLatex.h"
#include "TStyle.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TGraphErrors.h"

void LoadStyle();
void SetLegend(TLegend *);

TH1D* projectHistogram(THnSparseF *, double , double , double , double );
TProfile* projectProfile(THnSparseF *, double , double , double , double , int);

void plot_track_distributions() {
    gStyle -> SetPalette(kRainBow);

    int runList[] = {544510, 544508, 544490,544477, 544474, 544454, 544451, 544492, 544491, 544476, 
                     544475, 544392, 544391, 544390, 544389, 544185, 544184, 544180, 544124, 544123, 
                     544122, 544121, 544116, 544098, 544095, 544091, 544032, 544028, 544013};

    string pathToFiles = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/data/pass2/LHC23_golden/train_195781";
    string centFT0C[] = {"0 - 10%", "10 - 20%", "20 - 30%", "30 - 40%", "40 - 50%", "50 - 60%", "60 - 70%", "70 - 80%", "80 - 90%"};
    TH1D *histProjPsi2A[29][9];
    TH1D *histProjPsi2B[29][9];
    TH1D *histProjPsi2C[29][9];

    TH1F *histColCounterAll = new TH1F("histColCounterAll", "", 29, 0, 29);
    TH1F *histColCounterAcc = new TH1F("histColCounterAcc", "", 29, 0, 29);
    TH1F *histCounterTVX = new TH1F("histCounterTVX", "", 29, 0, 29);

    int index = 0;
    double counterColAll = 0;
    double counterColAcc = 0;
    double counterTVX = 0;

    for (auto const& run : runList) {
        std::cout << run << std::endl;
        histColCounterAll -> GetXaxis() -> SetBinLabel(index+1, Form("%i", run));
        histColCounterAcc -> GetXaxis() -> SetBinLabel(index+1, Form("%i", run));
        histCounterTVX -> GetXaxis() -> SetBinLabel(index+1, Form("%i", run));

        TFile *fIn = new TFile(Form("%s/%i/AnalysisResults.root", pathToFiles.c_str(), run));
        TList *list1 = (TList*) fIn -> Get("d-q-event-qvector/outputQA");
        TList *list2 = (TList*) list1 -> FindObject("Event_AfterCuts");
        TH2D *histPsi2ACentFT0C = (TH2D*)list2 -> FindObject("Psi2A_CentFT0C");
        TH2D *histPsi2BCentFT0C = (TH2D*)list2 -> FindObject("Psi2B_CentFT0C");
        TH2D *histPsi2CCentFT0C = (TH2D*)list2 -> FindObject("Psi2C_CentFT0C");
        histPsi2ACentFT0C -> RebinX(2);
        histPsi2BCentFT0C -> RebinX(2);
        histPsi2CCentFT0C -> RebinX(2);

        TH1F *histTmpColCounterAll = (TH1F*) fIn -> Get("event-selection-task/hColCounterAll");
        TH1F *histTmpColCounterAcc = (TH1F*) fIn -> Get("event-selection-task/hColCounterAcc");
        TH1F *histTmpCounterTVX = (TH1F*) fIn -> Get("bc-selection-task/hCounterTVX");

        histColCounterAll -> SetBinContent(index+1, histTmpColCounterAll -> GetBinContent(1));
        histColCounterAcc -> SetBinContent(index+1, histTmpColCounterAcc -> GetBinContent(1));
        histCounterTVX -> SetBinContent(index+1, histTmpCounterTVX -> GetBinContent(1));

        counterColAll += histTmpColCounterAll -> GetBinContent(1);
        counterColAcc += histTmpColCounterAcc -> GetBinContent(1);
        counterTVX += histTmpCounterTVX -> GetBinContent(1);
        
        for (int iCent = 0;iCent < 9;iCent++) {
            histProjPsi2A[index][iCent] = (TH1D*) histPsi2ACentFT0C -> ProjectionY(Form("Psi2ACentFT0C_%i_%i", iCent, run), iCent+1);
            histProjPsi2B[index][iCent] = (TH1D*) histPsi2BCentFT0C -> ProjectionY(Form("Psi2BCentFT0C_%i_%i", iCent, run), iCent+1);
            histProjPsi2C[index][iCent] = (TH1D*) histPsi2CCentFT0C -> ProjectionY(Form("Psi2CCentFT0C_%i_%i", iCent, run), iCent+1);
        }
        index++;
    }

    TLatex *latexTitle = new TLatex();
    latexTitle -> SetTextSize(0.065);
    latexTitle -> SetNDC();
    latexTitle -> SetTextFont(42);

    TCanvas *canvasPsi2A = new TCanvas("canvasPsi2A", "", 1800, 1800);
    canvasPsi2A -> Divide(3, 3);
    for (int iCent = 0;iCent < 9;iCent++) {
        canvasPsi2A -> cd(iCent+1);
        for (int i = 0;i < 29;i++) {
            histProjPsi2A[i][iCent] -> Scale(1. / histProjPsi2A[i][iCent] -> Integral());
            histProjPsi2A[i][iCent] -> GetYaxis() -> SetRangeUser(0, 0.03);
            histProjPsi2A[i][iCent] -> Draw("SAME PLC PMC");
        }
        gPad -> BuildLegend(0.7, 0.4, 0.980, 0.935, "", "L");
        latexTitle -> DrawLatex(0.20, 0.80, Form("FT0C %s", centFT0C[iCent].c_str()));
    }
    canvasPsi2A -> SaveAs("LHC23_golden/Psi2A_vs_CentFT0.pdf");

    TCanvas *canvasPsi2B = new TCanvas("canvasPsi2B", "", 1800, 1800);
    canvasPsi2B -> Divide(3, 3);
    for (int iCent = 0;iCent < 9;iCent++) {
        canvasPsi2B -> cd(iCent+1);
        for (int i = 0;i < 29;i++) {
            histProjPsi2B[i][iCent] -> Scale(1. / histProjPsi2B[i][iCent] -> Integral());
            histProjPsi2B[i][iCent] -> GetYaxis() -> SetRangeUser(0, 0.03);
            histProjPsi2B[i][iCent] -> Draw("SAME PLC PMC");
        }
        gPad -> BuildLegend(0.6, 0.6, 0.980, 0.935, "", "L");
        latexTitle -> DrawLatex(0.20, 0.80, Form("FT0C %s", centFT0C[iCent].c_str()));
    }
    canvasPsi2B -> SaveAs("LHC23_golden/Psi2B_vs_CentFT0.pdf");

    TCanvas *canvasPsi2C = new TCanvas("canvasPsi2C", "", 1800, 1800);
    canvasPsi2C -> Divide(3, 3);
    for (int iCent = 0;iCent < 9;iCent++) {
        canvasPsi2C -> cd(iCent+1);
        for (int i = 0;i < 29;i++) {
            histProjPsi2C[i][iCent] -> Scale(1. / histProjPsi2C[i][iCent] -> Integral());
            histProjPsi2C[i][iCent] -> GetYaxis() -> SetRangeUser(0, 0.03);
            histProjPsi2C[i][iCent] -> Draw("SAME PLC PMC");
        }
        gPad -> BuildLegend(0.6, 0.6, 0.980, 0.935, "", "L");
        latexTitle -> DrawLatex(0.20, 0.80, Form("FT0C %s", centFT0C[iCent].c_str()));
    }
    canvasPsi2C -> SaveAs("LHC23_golden/Psi2C_vs_CentFT0.pdf");

    TCanvas *canvasColCounterAcc = new TCanvas("canvasColCounterAcc", "", 1000, 600);
    gStyle -> SetOptStat(0);
    gPad -> SetLogy(1);
    histColCounterAll -> SetLineColor(kBlack);
    histColCounterAll -> GetYaxis() -> SetRangeUser(1e4, 1e11);
    histColCounterAll -> Draw("H");
    histColCounterAcc -> Draw("H SAME");
    latexTitle -> DrawLatex(0.20, 0.80, Form("counter collisions accepted = %1.0f", counterColAcc));

    TCanvas *canvasCounterTVX = new TCanvas("canvasCounterTVX", "", 1000, 600);
    gStyle -> SetOptStat(0);
    gPad -> SetLogy();
    histCounterTVX -> GetYaxis() -> SetRangeUser(1e4, 1e11);
    histCounterTVX -> Draw("H");
    latexTitle -> DrawLatex(0.20, 0.80, Form("counter TVX = %1.0f", counterTVX));


}

void plot_dimuon_distributions(bool integrated = false) {
    //string pathToFiles = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/data/pass2/LHC23_golden/dimuon_train_198686_with_cuts";
    string pathToFiles = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/data/pass2/LHC23_golden/dimuon_train_201824_no_cuts";
    //string muonCut = "muonLowPt10SigmaPDCA";
    string muonCut = "matchedMchMid";
    string suffix = "";
    //string suffix = "_new_EP_calculation";
    string methodName = "EP";
    int method = 0;

    if (methodName == "EP") {
        method = 4;
    } else {
        method = 3;
    }

    TFile *fIn  = new TFile(Form("%s%s/AnalysisResults.root", pathToFiles.c_str(), suffix.c_str()), "READ");
    // Same event pairing
    TList *list1 = (TList*) fIn -> Get("analysis-same-event-pairing/output");
    TList *listSE = (TList*) list1 -> FindObject(Form("PairsMuonSEPM_%s", muonCut.c_str()));
    THnSparseF* histSEPM = (THnSparseF*) listSE -> FindObject("Mass_Pt_centrFT0C_V2"); // 0: mass; 1: pT ; 2:FT0C_Centrality; 3:u2q2_A; 4:cos2(phi-psi)
    // Mixed event pairing
    TList *list2 = (TList*) fIn -> Get("analysis-event-mixing/output");
    TList *listME = (TList*) list2 -> FindObject(Form("PairsMuonMEPM_%s", muonCut.c_str()));
    THnSparseF* histMEPM = (THnSparseF*) listME -> FindObject("Mass_Pt_centrFT0C_V2");


    // v2 and pT distributions
    if (!integrated) {
        const int nPtBins = 5;
        double minCentrBins[] = {10};
        double maxCentrBins[] = {50};
        //double minPtBins[] = {0, 1, 2, 4, 6};
        //double maxPtBins[] = {1, 2, 4, 6, 10};
        double minPtBins[] = {0, 1, 2, 3, 5};
        double maxPtBins[] = {1, 2, 3, 5, 10};

        TFile *fOut = new TFile(Form("%s/Histograms_%s_centr_%1.f_%1.f.root", pathToFiles.c_str(), muonCut.c_str(), minCentrBins[0], maxCentrBins[0]), "RECREATE");

        for (int iPt = 0;iPt < nPtBins;iPt++) {
            TH1D *histMassSEPM = (TH1D*) projectHistogram(histSEPM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]);
            histMassSEPM -> SetName(Form("histMassSEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
            TH1D *histMassMEPM = (TH1D*) projectHistogram(histMEPM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]);
            histMassMEPM -> SetName(Form("histMassMEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
            TProfile *histV2SEPM = (TProfile*) projectProfile(histSEPM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0], method);
            histV2SEPM -> SetName(Form("histV2SEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
            TProfile *histV2MEPM = (TProfile*) projectProfile(histMEPM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0], method);
            histV2MEPM -> SetName(Form("histV2MEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
            TProfile *histPtSEPM = (TProfile*) projectProfile(histSEPM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0], 1);
            histPtSEPM -> SetName(Form("histPtSEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
            TProfile *histPtMEPM = (TProfile*) projectProfile(histMEPM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0], 1);
            histPtMEPM -> SetName(Form("histPtMEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));

            histMassSEPM -> SetTitle(Form("SE - %1.0f < #it{p}_{T} < %1.0f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
            histMassMEPM -> SetTitle(Form("ME - %1.0f < #it{p}_{T} < %1.0f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
            histV2SEPM -> SetTitle(Form("V2 SE - %1.0f < #it{p}_{T} < %1.0f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
            histV2MEPM -> SetTitle(Form("V2 ME - %1.0f < #it{p}_{T} < %1.0f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
            histPtSEPM -> SetTitle(Form("<pT> SE - %1.0f < #it{p}_{T} < %1.0f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
            histPtMEPM -> SetTitle(Form("<pT> ME - %1.0f < #it{p}_{T} < %1.0f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));

            histMassSEPM -> SetLineColor(kBlack);
            histMassSEPM -> SetMarkerStyle(20);
            histMassSEPM -> SetMarkerSize(0.8);
            histMassSEPM -> SetMarkerColor(kBlack);

            histV2SEPM -> SetLineColor(kBlack);
            histV2SEPM -> SetMarkerStyle(20);
            histV2SEPM -> SetMarkerSize(0.8);
            histV2SEPM -> SetMarkerColor(kBlack);

            histPtSEPM -> SetLineColor(kBlack);
            histPtSEPM -> SetMarkerStyle(20);
            histPtSEPM -> SetMarkerSize(0.8);
            histPtSEPM -> SetMarkerColor(kBlack);

            histMassMEPM -> SetLineColor(kAzure+4);
            histV2MEPM -> SetLineColor(kAzure+4);
            histPtMEPM -> SetLineColor(kAzure+4);

            Printf("Compute normalization...\n");
            double minMassBin1 = histMassSEPM -> GetXaxis() -> FindBin(2.00);
            double maxMassBin1 = histMassSEPM -> GetXaxis() -> FindBin(2.50);
            double minMassBin2 = histMassSEPM -> GetXaxis() -> FindBin(3.72);
            double maxMassBin2 = histMassSEPM -> GetXaxis() -> FindBin(5.00);

            double normMassSEPM = histMassSEPM -> Integral(minMassBin1, maxMassBin1) + histMassSEPM->Integral(minMassBin2, maxMassBin2);
            double normMassMEPM = histMassMEPM -> Integral(minMassBin1, maxMassBin1) + histMassMEPM->Integral(minMassBin2, maxMassBin2);

            Printf("Mass normalization factor = %f\n", normMassSEPM / normMassMEPM);


            TCanvas *canvasMassV2Pt = new TCanvas("canvasMassV2Pt", "", 600, 1800);
            canvasMassV2Pt -> Divide(1, 3);

            canvasMassV2Pt -> cd(1);
            histMassSEPM -> GetXaxis() -> SetRangeUser(2, 5);
            histMassMEPM -> GetXaxis() -> SetRangeUser(2, 5);
            histMassSEPM -> Draw("EP");
            histMassMEPM -> Scale(normMassSEPM / normMassMEPM);
            histMassMEPM -> Draw("EP SAME");

            TH1D *histBkgSubtrSEMP = (TH1D*) histMassSEPM -> Clone(Form("histMassBkgSubtrSEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
            histBkgSubtrSEMP -> Add(histMassMEPM, -1.);
            histBkgSubtrSEMP -> SetLineColor(kRed+1);
            histBkgSubtrSEMP -> SetMarkerStyle(24);
            histBkgSubtrSEMP -> SetMarkerColor(kRed+1);
            histBkgSubtrSEMP -> Draw("EP SAME");

            canvasMassV2Pt -> cd(2);
            histV2SEPM -> GetXaxis() -> SetRangeUser(2, 5);
            histV2SEPM -> Draw("EP");
            histV2MEPM -> Scale(histV2SEPM -> GetBinContent(histV2SEPM -> GetXaxis() -> FindBin(2)) / histV2MEPM -> GetBinContent(histV2MEPM -> GetXaxis() -> FindBin(2)));
            histV2MEPM -> Draw("H SAME");

            canvasMassV2Pt -> cd(3);
            histPtSEPM -> GetXaxis() -> SetRangeUser(2, 5);
            histPtSEPM -> Draw("EP");
            histPtMEPM -> Draw("H SAME");

            canvasMassV2Pt -> SaveAs(Form("LHC23_golden/flow/%s_distrib_centr_%s_%1.0f_%1.0f__%1.0f_%1.0f%s.pdf", muonCut.c_str(), methodName.c_str(), minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0], suffix.c_str()));
        
            fOut -> cd();
            histMassSEPM -> Write();
            histMassMEPM -> Write();
            histV2SEPM -> Write();
            histV2MEPM -> Write();
            histPtSEPM -> Write();
            histPtMEPM -> Write();
        }
        fOut -> Close();
    }

    if (integrated) {
        TFile *fOut = new TFile(Form("%s/Histograms_%s.root", pathToFiles.c_str(), muonCut.c_str()), "RECREATE");
        // pT scan
        //const int nPtBins = 9;
        //double minPtBins[] = {0, 1, 2, 3, 5, 7, 10, 15, 20};
        //double maxPtBins[] = {1, 2, 3, 5, 7, 10, 15, 20, 30};
        const int nPtBins = 4;
        double minPtBins[] = {0, 2, 4, 6};
        double maxPtBins[] = {2, 4, 6, 12};

        TCanvas *canvasMassPtScan = new TCanvas("canvasMassPtScan", "", 1800, 1800);
        canvasMassPtScan -> Divide(3, 3);

        for (int iPt = 0;iPt < nPtBins;iPt++) {
            TH1D *histMassSEPM = (TH1D*) projectHistogram(histSEPM, minPtBins[iPt], maxPtBins[iPt], 0., 90.);
            TH1D *histMassMEPM = (TH1D*) projectHistogram(histMEPM, minPtBins[iPt], maxPtBins[iPt], 0., 90.);

            histMassSEPM -> SetName(Form("histMassSEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], 0., 90.));
            histMassMEPM -> SetName(Form("histMassMEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], 0., 90.));
            histMassSEPM -> SetTitle(Form("SE - %1.0f < #it{p}_{T} < %1.0f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], 0., 90.));
            histMassMEPM -> SetTitle(Form("ME - %1.0f < #it{p}_{T} < %1.0f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], 0., 90.));

            histMassSEPM -> SetLineColor(kBlack);
            histMassSEPM -> SetMarkerStyle(20);
            histMassSEPM -> SetMarkerSize(0.8);
            histMassSEPM -> SetMarkerColor(kBlack);
            histMassMEPM -> SetLineColor(kAzure+4);

            Printf("Compute normalization...\n");
            double minMassBin1 = histMassSEPM -> GetXaxis() -> FindBin(1.00);
            double maxMassBin1 = histMassSEPM -> GetXaxis() -> FindBin(2.50);
            double minMassBin2 = histMassSEPM -> GetXaxis() -> FindBin(3.72);
            double maxMassBin2 = histMassSEPM -> GetXaxis() -> FindBin(5.00);

            double normMassSEPM = histMassSEPM -> Integral(minMassBin1, maxMassBin1) + histMassSEPM->Integral(minMassBin2, maxMassBin2);
            double normMassMEPM = histMassMEPM -> Integral(minMassBin1, maxMassBin1) + histMassMEPM->Integral(minMassBin2, maxMassBin2);

            Printf("Normalization factor = %f\n", normMassSEPM / normMassMEPM);

            canvasMassPtScan -> cd(iPt+1);
            gPad -> SetLogy(1);
            histMassSEPM -> GetXaxis() -> SetRangeUser(2, 5);
            histMassSEPM -> GetYaxis() -> SetRangeUser(1, 1e7);
            histMassMEPM -> GetXaxis() -> SetRangeUser(2, 5);
            histMassSEPM -> Draw("EP");
            histMassMEPM -> Scale(normMassSEPM / normMassMEPM);
            histMassMEPM -> Draw("EP SAME");

            TH1D *histBkgSubtrSEMP = (TH1D*) histMassSEPM -> Clone(Form("histMassBkgSubtrSEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], 0., 90.));
            histBkgSubtrSEMP -> Add(histMassMEPM, -1.);
            histBkgSubtrSEMP -> SetLineColor(kRed+1);
            histBkgSubtrSEMP -> SetMarkerStyle(24);
            histBkgSubtrSEMP -> SetMarkerSize(0.8);
            histBkgSubtrSEMP -> SetMarkerColor(kRed+1);
            histBkgSubtrSEMP -> Draw("EP SAME");

            fOut -> cd();
            histMassSEPM -> Write();
            histMassMEPM -> Write();
            histBkgSubtrSEMP -> Write();
        }

        // Centrality scan
        const int nCentrBins = 9;
        double minCentrBins[] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
        double maxCentrBins[] = {10, 20, 30, 40, 50, 60, 70, 80, 90};

        TCanvas *canvasCentrScan = new TCanvas("canvasCentrScan", "", 1800, 1800);
        canvasCentrScan -> Divide(3, 3);

        TCanvas *canvasV2CentrScan = new TCanvas("canvasV2CentrScan", "", 1800, 1800);
        canvasV2CentrScan -> Divide(3, 3);


        for (int iCentr = 0;iCentr < nCentrBins;iCentr++) {
            TH1D *histMassSEPM = (TH1D*) projectHistogram(histSEPM, 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]);
            histMassSEPM -> SetName(Form("histMassSEPM_%1.0f_%1.0f__%1.0f_%1.0f", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));
            TH1D *histMassMEPM = (TH1D*) projectHistogram(histMEPM, 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]);
            histMassMEPM -> SetName(Form("histMassMEPM_%1.0f_%1.0f__%1.0f_%1.0f", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));
            TProfile *histV2SEPM = (TProfile*) projectProfile(histSEPM, 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr], method);
            histV2SEPM -> SetName(Form("histV2SEPM_%1.0f_%1.0f__%1.0f_%1.0f", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));
            TProfile *histV2MEPM = (TProfile*) projectProfile(histMEPM, 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr], method);
            histV2MEPM -> SetName(Form("histV2MEPM_%1.0f_%1.0f__%1.0f_%1.0f", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));

            histMassSEPM -> SetTitle(Form("SE - %1.0f < #it{p}_{T} < %1.0f GeV/#it{c}, %1.0f #minus %1.0f %%", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));
            histMassMEPM -> SetTitle(Form("ME - %1.0f < #it{p}_{T} < %1.0f GeV/#it{c}, %1.0f #minus %1.0f %%", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));
            histV2SEPM -> SetTitle(Form("V2 SE - %1.0f < #it{p}_{T} < %1.0f GeV/#it{c}, %1.0f #minus %1.0f %%", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));
            histV2MEPM -> SetTitle(Form("V2 ME - %1.0f < #it{p}_{T} < %1.0f GeV/#it{c}, %1.0f #minus %1.0f %%", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));

            histMassSEPM -> SetLineColor(kBlack);
            histMassSEPM -> SetMarkerStyle(20);
            histMassSEPM -> SetMarkerSize(0.8);
            histMassSEPM -> SetMarkerColor(kBlack);
            histMassMEPM -> SetLineColor(kAzure+4);

            histV2SEPM -> SetLineColor(kBlack);
            histV2SEPM -> SetMarkerStyle(20);
            histV2SEPM -> SetMarkerSize(0.8);
            histV2SEPM -> SetMarkerColor(kBlack);
            histV2MEPM -> SetLineColor(kAzure+4);

            Printf("Compute normalization...\n");
            double minMassBin1 = histMassSEPM -> GetXaxis() -> FindBin(1.00);
            double maxMassBin1 = histMassSEPM -> GetXaxis() -> FindBin(2.50);
            double minMassBin2 = histMassSEPM -> GetXaxis() -> FindBin(3.72);
            double maxMassBin2 = histMassSEPM -> GetXaxis() -> FindBin(5.00);

            double minV2Bin1 = histV2SEPM -> GetXaxis() -> FindBin(1.00);
            double maxV2Bin1 = histV2SEPM -> GetXaxis() -> FindBin(2.50);
            double minV2Bin2 = histV2SEPM -> GetXaxis() -> FindBin(3.72);
            double maxV2Bin2 = histV2SEPM -> GetXaxis() -> FindBin(5.00);

            double normMassSEPM = histMassSEPM -> Integral(minMassBin1, maxMassBin1) + histMassSEPM->Integral(minMassBin2, maxMassBin2);
            double normMassMEPM = histMassMEPM -> Integral(minMassBin1, maxMassBin1) + histMassMEPM->Integral(minMassBin2, maxMassBin2);
            double normV2SEPM = histV2SEPM -> Integral(minV2Bin1, maxV2Bin1) + histV2SEPM->Integral(minV2Bin2, maxV2Bin2);
            double normV2MEPM = histV2MEPM -> Integral(minV2Bin1, maxV2Bin1) + histV2MEPM->Integral(minV2Bin2, maxV2Bin2);

            Printf("Mass normalization factor = %f\n", normMassSEPM / normMassMEPM);
            Printf("V2 normalization factor = %f (%f - %f)\n", normV2SEPM / normV2MEPM, normV2SEPM, normV2MEPM);

            canvasCentrScan -> cd(iCentr+1);
            gPad -> SetLogy(1);
            histMassSEPM -> GetXaxis() -> SetRangeUser(2, 5);
            histMassSEPM -> GetYaxis() -> SetRangeUser(1, 1e7);
            histMassMEPM -> GetXaxis() -> SetRangeUser(2, 5);
            histMassSEPM -> Draw("EP");
            histMassMEPM -> Scale(normMassSEPM / normMassMEPM);
            histMassMEPM -> Draw("EP SAME");

            TH1D *histBkgSubtrSEMP = (TH1D*) histMassSEPM -> Clone(Form("histMassBkgSubtrSEPM_%1.0f_%1.0f__%1.0f_%1.0f", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));
            histBkgSubtrSEMP -> Add(histMassMEPM, -1.);
            histBkgSubtrSEMP -> SetLineColor(kRed+1);
            histBkgSubtrSEMP -> SetMarkerStyle(24);
            histBkgSubtrSEMP -> SetMarkerSize(0.8);
            histBkgSubtrSEMP -> SetMarkerColor(kRed+1);
            histBkgSubtrSEMP -> Draw("EP SAME");

            canvasV2CentrScan -> cd(iCentr+1);
            histV2SEPM -> GetXaxis() -> SetRangeUser(1, 5);
            histV2SEPM -> Draw("EP");
            //histV2MEPM -> Scale(histV2SEPM -> GetBinContent(histV2SEPM -> GetXaxis() -> FindBin(2)) / histV2MEPM -> GetBinContent(histV2MEPM -> GetXaxis() -> FindBin(2)));
            histV2MEPM -> Scale(normV2SEPM / normV2MEPM);
            histV2MEPM -> Draw("H SAME");

            fOut -> cd();
            histMassSEPM -> Write();
            histMassMEPM -> Write();
            histBkgSubtrSEMP -> Write();
        }

        canvasMassPtScan -> SaveAs(Form("LHC23_golden/flow/plot_mass_pt_scan_%s.pdf", muonCut.c_str()));
        canvasCentrScan -> SaveAs(Form("LHC23_golden/flow/plot_mass_centrality_scan_%s.pdf", muonCut.c_str()));
        canvasV2CentrScan -> SaveAs(Form("LHC23_golden/flow/plot_v2_centrality_scan_%s.pdf", muonCut.c_str()));

        fOut -> ls();
        fOut -> Close();
    }
}

void plot_results () {
    LoadStyle();

    // Run 2
    double ptBinsRun2[] = {0.0, 2.0, 4.0, 6.0, 12.0};
    double centrBinsRun2[] = {0, 20, 40, 60, 90};
    
    double widthPtRun2[] = {0.0685, 0.0681, 0.0692, 0.0698};
    double errWidthPtRun2[] = {0.0005, 0.0005, 0.0007, 0.0009};

    double widthCentrRun2[] = {0.0667, 0.0675, 0.0677, 0.0664};
    double errWidthCentrRun2[] = {0.0004, 0.0005, 0.0006, 0.0007};

    double sigJpsiPtRun2[] = {465742, 315216, 98463, 43169};
    double errSigJpsiPtRun2[] = {3214, 1901, 861, 425};

    TH1D *histPtWidthRun2 = new TH1D("histPtWidthRun2", "", 4, ptBinsRun2);
    histPtWidthRun2 -> SetLineColor(kBlack);
    histPtWidthRun2 -> SetMarkerStyle(20);
    histPtWidthRun2 -> SetMarkerColor(kBlack);

    TH1D *histCentrWidthRun2 = new TH1D("histCentrWidthRun2", "", 4, centrBinsRun2);
    histCentrWidthRun2 -> SetLineColor(kBlack);
    histCentrWidthRun2 -> SetMarkerStyle(20);
    histCentrWidthRun2 -> SetMarkerColor(kBlack);

    TH1D *histPtSigJpsiRun2 = new TH1D("histPtSigJpsiRun2", "", 4, ptBinsRun2);
    histPtSigJpsiRun2 -> SetLineColor(kBlack);
    histPtSigJpsiRun2 -> SetMarkerStyle(20);
    histPtSigJpsiRun2 -> SetMarkerColor(kBlack);


    for (int idx = 0;idx < 4;idx++) {
        histPtWidthRun2 -> SetBinContent(idx+1, widthPtRun2[idx]);
        histPtWidthRun2 -> SetBinError(idx+1, errWidthPtRun2[idx]);

        histCentrWidthRun2 -> SetBinContent(idx+1, widthCentrRun2[idx]);
        histCentrWidthRun2 -> SetBinError(idx+1, errWidthCentrRun2[idx]);

        histPtSigJpsiRun2 -> SetBinContent(idx+1, sigJpsiPtRun2[idx]);
        histPtSigJpsiRun2 -> SetBinError(idx+1, errSigJpsiPtRun2[idx]);
    }

    double ptBinsRun3[] = {0.0, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0};
    double centrBinsRun3[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};

    double widthPtRun3[] = {0.067, 0.070, 0.069, 0.075, 0.081, 0.088, 0.096};
    double errWidthPtRun3[] = {0.002, 0.001, 0.001, 0.001, 0.003, 0.004, 0.013};

    double widthCentrRun3[] = {0.067, 0.069, 0.073, 0.069, 0.076, 0.075, 0.089, 0.075, 0.075};
    double errWidthCentrRun3[] = {0.001, 0.001, 0.002, 0.001, 0.001, 0.002, 0.003, 0.003, 0.004};

    double sigJpsiPtRun3[] = {48352, 86537, 62338, 26232, 16646, 5339, 1289};
    double errSigJpsiPtRun3[] = {1412, 1478, 1073, 559, 410, 262, 129};

    double sigJpsiPtRebinRun3[] = {250811, 170352, 47394, 18019};
    double errSigJpsiPtRebinRun3[] = {16664, 2144, 912, 469};

    TH1D *histPtWidthRun3 = new TH1D("histPtWidthRun3", "", 7, ptBinsRun3);
    histPtWidthRun3 -> SetLineColor(kRed+1);
    histPtWidthRun3 -> SetMarkerStyle(20);
    histPtWidthRun3 -> SetMarkerColor(kRed+1);

    TH1D *histCentrWidthRun3 = new TH1D("histCentrWidthRun3", "", 9, centrBinsRun3);
    histCentrWidthRun3 -> SetLineColor(kRed+1);
    histCentrWidthRun3 -> SetMarkerStyle(20);
    histCentrWidthRun3 -> SetMarkerColor(kRed+1);

    TH1D *histPtSigJpsiRun3 = new TH1D("histPtSigJpsiRun3", "", 7, ptBinsRun3);
    histPtSigJpsiRun3 -> SetLineColor(kRed+1);
    histPtSigJpsiRun3 -> SetMarkerStyle(20);
    histPtSigJpsiRun3 -> SetMarkerColor(kRed+1);

    TH1D *histPtRebinSigJpsiRun3 = new TH1D("histPtRebinSigJpsiRun3", "", 4, ptBinsRun2);
    histPtRebinSigJpsiRun3 -> SetLineColor(kRed+1);
    histPtRebinSigJpsiRun3 -> SetMarkerStyle(20);
    histPtRebinSigJpsiRun3 -> SetMarkerColor(kRed+1);


    for (int idx = 0;idx < 9;idx++) {
        histCentrWidthRun3 -> SetBinContent(idx+1, widthCentrRun3[idx]);
        histCentrWidthRun3 -> SetBinError(idx+1, errWidthCentrRun3[idx]);

        if (idx < 7) {
            histPtWidthRun3 -> SetBinContent(idx+1, widthPtRun3[idx]);
            histPtWidthRun3 -> SetBinError(idx+1, errWidthPtRun3[idx]);

            histPtSigJpsiRun3 -> SetBinContent(idx+1, sigJpsiPtRun3[idx]);
            histPtSigJpsiRun3 -> SetBinError(idx+1, errSigJpsiPtRun3[idx]);
        }

        if (idx < 4) {
            histPtRebinSigJpsiRun3 -> SetBinContent(idx+1, sigJpsiPtRebinRun3[idx]);
            histPtRebinSigJpsiRun3 -> SetBinError(idx+1, errSigJpsiPtRebinRun3[idx]);
        }
    }


    TCanvas *canvasPtWidth = new TCanvas("canvasPtWidth", "", 800, 600);
    histPtWidthRun3 -> SetTitle("");
    histPtWidthRun3 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histPtWidthRun3 -> GetYaxis() -> SetTitle("#sigma_{J/#psi} (GeV/#it{c}^{2})");
    histPtWidthRun3 -> GetYaxis() -> SetRangeUser(0, 0.15);
    histPtWidthRun3 -> Draw("EP");
    histPtWidthRun2 -> Draw("EP SAME");

    TCanvas *canvasCentrWidth = new TCanvas("canvasCentrWidth", "", 800, 600);
    histCentrWidthRun3 -> SetTitle("");
    histCentrWidthRun3 -> GetXaxis() -> SetTitle("FT0 centrality (%)");
    histCentrWidthRun3 -> GetYaxis() -> SetTitle("#sigma_{J/#psi} (GeV/#it{c}^{2})");
    histCentrWidthRun3 -> GetYaxis() -> SetRangeUser(0, 0.15);
    histCentrWidthRun3 -> Draw("EP");
    histCentrWidthRun2 -> Draw("EP SAME");

    TCanvas *canvasPtSigJpsi = new TCanvas("canvasPtSigJpsi", "", 800, 600);
    gPad -> SetLogy(1);
    //histPtSigJpsiRun3 -> Scale(1. / histPtSigJpsiRun3 -> Integral(), "WIDTH");
    //histPtRebinSigJpsiRun3 -> Scale(1. / histPtRebinSigJpsiRun3 -> Integral(), "WIDTH");
    //histPtSigJpsiRun2 -> Scale(1. / histPtSigJpsiRun2 -> Integral(), "WIDTH");

    histPtSigJpsiRun3 -> Scale(1. / 0.30, "WIDTH");
    histPtRebinSigJpsiRun3 -> Scale(1. / 0.30, "WIDTH");
    histPtSigJpsiRun2 -> Scale(1. / 0.75, "WIDTH");

    histPtSigJpsiRun2 -> SetTitle("");
    histPtSigJpsiRun2 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histPtSigJpsiRun2 -> GetYaxis() -> SetTitle("Normalized N_{J/#psi}");
    //histPtSigJpsiRun3 -> Draw("EP");
    histPtSigJpsiRun2 -> Draw("EP");
    histPtRebinSigJpsiRun3 -> Draw("EP SAME");


    // Resolution extrapolation
    double massesRun2[] = {3.098, 9.455};
    double errMassesRun2[] = {0.001, 0.005};
    double widthsRun2[] = {0.069, 0.126};
    double errWidthRun2[] = {0.002, 0.006};

    double massesRun3[] = {3.081, 9.385};
    double errMassesRun3[] = {0.001, 0.029};
    double widthsRun3[] = {0.073, 0.187};
    double errWidthRun3[] = {0.001, 0.024};

    TGraphErrors *graWidthVsMassRun2 = new TGraphErrors(2, massesRun2, widthsRun2, errMassesRun2, errWidthRun2);
    graWidthVsMassRun2 -> SetLineColor(kBlack);
    graWidthVsMassRun2 -> SetMarkerStyle(20);
    graWidthVsMassRun2 -> SetMarkerColor(kBlack);

    TGraphErrors *graWidthVsMassRun3 = new TGraphErrors(2, massesRun3, widthsRun3, errMassesRun3, errWidthRun3);
    graWidthVsMassRun3 -> SetLineColor(kRed+1);
    graWidthVsMassRun3 -> SetMarkerStyle(20);
    graWidthVsMassRun3 -> SetMarkerColor(kRed+1);


    TCanvas *canvasWidthVsMass = new TCanvas("canvasWidthVsMass", "", 800, 600);
    gPad -> SetLogy(1);
    TH2D *histGridWidthVsMass = new TH2D("histGridWidthVsMass", "", 100, 0, 120, 100, 0.01, 1.);
    histGridWidthVsMass -> GetXaxis() -> SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    histGridWidthVsMass -> GetYaxis() -> SetTitle("#sigma (GeV/#it{c}^{2})");
    histGridWidthVsMass -> Draw();
    graWidthVsMassRun2 -> Draw("EP SAME");
    graWidthVsMassRun3 -> Draw("EP SAME");



    canvasPtWidth -> SaveAs("LHC23_golden/flow/width_vs_pT.pdf");
    canvasCentrWidth -> SaveAs("LHC23_golden/flow/width_vs_centrality.pdf");
    canvasPtSigJpsi -> SaveAs("LHC23_golden/flow/signal_vs_pT.pdf");
}
////////////////////////////////////////////////////////////////////////////////
void LoadStyle(){
    int font = 42;
    gStyle -> SetFrameBorderMode(0);
    gStyle -> SetFrameFillColor(0);
    gStyle -> SetCanvasBorderMode(0);
    gStyle -> SetPadBorderMode(0);
    gStyle -> SetPadColor(10);
    gStyle -> SetCanvasColor(10);
    gStyle -> SetTitleFillColor(10);
    gStyle -> SetTitleBorderSize(1);
    gStyle -> SetStatColor(10);
    gStyle -> SetStatBorderSize(1);
    gStyle -> SetLegendBorderSize(1);
    gStyle -> SetDrawBorder(0);
    gStyle -> SetTextFont(font);
    gStyle -> SetStatFontSize(0.05);
    gStyle -> SetStatX(0.97);
    gStyle -> SetStatY(0.98);
    gStyle -> SetStatH(0.03);
    gStyle -> SetStatW(0.3);
    gStyle -> SetTickLength(0.02,"y");
    gStyle -> SetEndErrorSize(3);
    gStyle -> SetLabelSize(0.05,"xyz");
    gStyle -> SetLabelFont(font,"xyz");
    gStyle -> SetLabelOffset(0.01,"xyz");
    gStyle -> SetTitleFont(font,"xyz");
    gStyle -> SetTitleOffset(0.9,"x");
    gStyle -> SetTitleOffset(1.02,"y");
    gStyle -> SetTitleSize(0.05,"xyz");
    gStyle -> SetMarkerSize(1.3);
    gStyle -> SetOptStat(0);
    gStyle -> SetEndErrorSize(0);
    gStyle -> SetCanvasPreferGL(kTRUE);
    gStyle -> SetHatchesSpacing(0.5);
    gStyle -> SetPadLeftMargin(0.15);
    gStyle -> SetPadBottomMargin(0.15);
    gStyle -> SetPadTopMargin(0.05);
    gStyle -> SetPadRightMargin(0.05);
    gStyle -> SetEndErrorSize(0.0);
    gStyle -> SetTitleSize(0.05,"X");
    gStyle -> SetTitleSize(0.045,"Y");
    gStyle -> SetLabelSize(0.045,"X");
    gStyle -> SetLabelSize(0.045,"Y");
    gStyle -> SetTitleOffset(1.2,"X");
    gStyle -> SetTitleOffset(1.35,"Y");
}
////////////////////////////////////////////////////////////////////////////////
void SetLegend(TLegend *legend){
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.04);
}

TH1D* projectHistogram(THnSparseF *histSparse, double minPtRange, double maxPtRange, double minCentrRange, double maxCentrRange) {
    double minPtBin = histSparse -> GetAxis(1) -> FindBin(minPtRange);
    double maxPtBin = histSparse -> GetAxis(1) -> FindBin(maxPtRange - 0.01);
    double minCentrBin = histSparse -> GetAxis(2) -> FindBin(minCentrRange);
    double maxCentrBin = histSparse -> GetAxis(2) -> FindBin(maxCentrRange - 0.01);

    Printf("minPtBin = %1.0f, maxPtBin = %1.0f, minCentrBin = %1.0f, maxCentrBin = %1.0f", minPtBin, maxPtBin, minCentrBin, maxCentrBin);
    histSparse -> GetAxis(1) -> SetRange(minPtBin, maxPtBin);
    histSparse -> GetAxis(2) -> SetRange(minCentrBin, maxCentrBin);

    TH1D *histProj = (TH1D*) histSparse -> Projection(0, Form("histProj_%1.0f_%1.0f__%1.0f_%1.0f", minPtRange, maxPtRange, minCentrRange, maxCentrRange));
    return histProj;
}

TProfile* projectProfile(THnSparseF *histSparse, double minPtRange, double maxPtRange, double minCentrRange, double maxCentrRange, int var) {
    double minPtBin = histSparse -> GetAxis(1) -> FindBin(minPtRange);
    double maxPtBin = histSparse -> GetAxis(1) -> FindBin(maxPtRange - 0.01);
    double minCentrBin = histSparse -> GetAxis(2) -> FindBin(minCentrRange);
    double maxCentrBin = histSparse -> GetAxis(2) -> FindBin(maxCentrRange - 0.01);

    Printf("minPtBin = %1.0f, maxPtBin = %1.0f, minCentrBin = %1.0f, maxCentrBin = %1.0f", minPtBin, maxPtBin, minCentrBin, maxCentrBin);
    histSparse -> GetAxis(1) -> SetRange(minPtBin, maxPtBin);
    histSparse -> GetAxis(2) -> SetRange(minCentrBin, maxCentrBin);

    TH2D *hist2DProj = (TH2D*) histSparse -> Projection(var, 0, "Projection");
    hist2DProj -> Rebin2D(2);
    TProfile *histProj = (TProfile*) hist2DProj -> ProfileX(Form("histProj_%i_%1.0f_%1.0f__%1.0f_%1.0f", var, minPtRange, maxPtRange, minCentrRange, maxCentrRange));
    return histProj;
}