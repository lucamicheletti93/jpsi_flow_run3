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

void subtract_background() {
    //string pathToFiles = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/data/pass3/LHC23_golden/train_231992"; // 20% golden Pb-Pb, muonLowPt 0.7 GeV/c
    //string cuts[] = {"muonLowPt210SigmaPDCA"};
    //string pathToFiles = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/data/pass3/LHC23_full/train_235036"; // full Pb-Pb, muonLowPt 0.7 GeV/c
    //string cuts[] = {"muonLowPt210SigmaPDCA"};
    string pathToFiles = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/data/pass3/LHC23_full/train_235531"; // full Pb-Pb, matchedcedMchMid
    string cuts[] = {"matchedMchMid"};

    const int nPtBins = 12;
    double minPtBins[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20};
    double maxPtBins[] = {1, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20, 30};

    const int nCentrBins = 9;
    double minCentrBins[] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
    double maxCentrBins[] = {10, 20, 30, 40, 50, 60, 70, 80, 90};


    TFile *fIn = new TFile(Form("%s/Histograms_%s_merged.root", pathToFiles.c_str(), cuts[0].c_str()));
    TH1D *histMassIntSEPM = (TH1D*) fIn -> Get("Mass_Int_SEPM");
    TH1D *histMassIntMEPM = (TH1D*) fIn -> Get("Mass_Int_MEPM");

    double minMassBin1 = histMassIntSEPM -> GetXaxis() -> FindBin(1.00);
    double maxMassBin1 = histMassIntSEPM -> GetXaxis() -> FindBin(2.50);
    double minMassBin2 = histMassIntSEPM -> GetXaxis() -> FindBin(3.92);
    double maxMassBin2 = histMassIntSEPM -> GetXaxis() -> FindBin(5.00);

    double normMassSEPM = histMassIntSEPM -> Integral(minMassBin1, maxMassBin1) + histMassIntSEPM -> Integral(minMassBin2, maxMassBin2);
    double normMassMEPM = histMassIntMEPM -> Integral(minMassBin1, maxMassBin1) + histMassIntMEPM -> Integral(minMassBin2, maxMassBin2);

    TCanvas *canvasIntegrated = new TCanvas("canvasIntegrated", "", 800, 600);
    gPad -> SetLogy(1);
    histMassIntSEPM -> SetTitle("0 < #it{p}_{T} < 30 GeV/#it{c}, 0 - 90 %");
    histMassIntSEPM -> GetXaxis() -> SetRangeUser(2, 5);
    histMassIntSEPM -> GetYaxis() -> SetRangeUser(1e2, 1e7);
    histMassIntSEPM -> SetMarkerStyle(20);
    histMassIntSEPM -> SetMarkerSize(0.5);
    histMassIntMEPM -> GetXaxis() -> SetRangeUser(2, 5);
    histMassIntSEPM -> Draw("EP");
    histMassIntMEPM -> Scale(normMassSEPM / normMassMEPM);
    histMassIntMEPM -> Draw("H SAME");

    TH1D *histBkgSubtrIntSEMP = (TH1D*) histMassIntSEPM -> Clone("histMassBkgSubtrSEPM_integrated");
    histBkgSubtrIntSEMP -> Add(histMassIntMEPM, -1.);
    histBkgSubtrIntSEMP -> SetLineColor(kRed+1);
    histBkgSubtrIntSEMP -> SetMarkerStyle(24);
    histBkgSubtrIntSEMP -> SetMarkerSize(0.5);
    histBkgSubtrIntSEMP -> SetMarkerColor(kRed+1);
    histBkgSubtrIntSEMP -> Draw("EP SAME");

    TFile *fOut = new TFile(Form("histograms_%s.root", cuts[0].c_str()), "RECREATE");
    histMassIntSEPM -> Write();
    histMassIntMEPM -> Write();
    histBkgSubtrIntSEMP -> Write();
    fOut -> Close();

    TCanvas *canvasPtScan = new TCanvas("canvasPtScan", "", 1800, 1800);
    canvasPtScan -> Divide(4, 3);
    for (int iPt = 0;iPt < nPtBins;iPt++) {
        TH1D *histMassSEPM = (TH1D*) fIn -> Get(Form("Mass_Pt_%2.1f_%2.1f_SEPM", minPtBins[iPt], maxPtBins[iPt]));
        TH1D *histMassMEPM = (TH1D*) fIn -> Get(Form("Mass_Pt_%2.1f_%2.1f_MEPM", minPtBins[iPt], maxPtBins[iPt]));

        minMassBin1 = histMassSEPM -> GetXaxis() -> FindBin(1.00);
        maxMassBin1 = histMassSEPM -> GetXaxis() -> FindBin(2.50);
        minMassBin2 = histMassSEPM -> GetXaxis() -> FindBin(3.92);
        maxMassBin2 = histMassSEPM -> GetXaxis() -> FindBin(5.00);

        normMassSEPM = histMassSEPM -> Integral(minMassBin1, maxMassBin1) + histMassSEPM -> Integral(minMassBin2, maxMassBin2);
        normMassMEPM = histMassMEPM -> Integral(minMassBin1, maxMassBin1) + histMassMEPM -> Integral(minMassBin2, maxMassBin2);

        canvasPtScan -> cd(iPt+1);
        gPad -> SetLogy(1);
        histMassSEPM -> SetTitle(Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}, 0 - 90 %%", minPtBins[iPt], maxPtBins[iPt]));
        histMassSEPM -> GetXaxis() -> SetRangeUser(2, 5);
        histMassSEPM -> GetYaxis() -> SetRangeUser(1, 1e7);
        histMassSEPM -> SetMarkerStyle(20);
        histMassSEPM -> SetMarkerSize(0.5);
        histMassMEPM -> GetXaxis() -> SetRangeUser(2, 5);
        histMassSEPM -> Draw("EP");
        histMassMEPM -> Scale(normMassSEPM / normMassMEPM);
        histMassMEPM -> Draw("H SAME");

        TH1D *histBkgSubtrSEMP = (TH1D*) histMassSEPM -> Clone(Form("histMassBkgSubtrSEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], 0., 90.));
        histBkgSubtrSEMP -> Add(histMassMEPM, -1.);
        histBkgSubtrSEMP -> SetLineColor(kRed+1);
        histBkgSubtrSEMP -> SetMarkerStyle(24);
        histBkgSubtrSEMP -> SetMarkerSize(0.5);
        histBkgSubtrSEMP -> SetMarkerColor(kRed+1);
        histBkgSubtrSEMP -> Draw("EP SAME");
    }




    TCanvas *canvasCentrScan = new TCanvas("canvasCentrScan", "", 1800, 1800);
    canvasCentrScan -> Divide(3, 3);
    for (int iCentr = 0;iCentr < nCentrBins;iCentr++) {
        TH1D *histMassSEPM = (TH1D*) fIn -> Get(Form("Mass_CentrFT0C_%2.1f_%2.1f_SEPM", minCentrBins[iCentr], maxCentrBins[iCentr]));
        TH1D *histMassMEPM = (TH1D*) fIn -> Get(Form("Mass_CentrFT0C_%2.1f_%2.1f_MEPM", minCentrBins[iCentr], maxCentrBins[iCentr]));

        minMassBin1 = histMassSEPM -> GetXaxis() -> FindBin(1.00);
        maxMassBin1 = histMassSEPM -> GetXaxis() -> FindBin(2.50);
        minMassBin2 = histMassSEPM -> GetXaxis() -> FindBin(3.92);
        maxMassBin2 = histMassSEPM -> GetXaxis() -> FindBin(5.00);

        normMassSEPM = histMassSEPM -> Integral(minMassBin1, maxMassBin1) + histMassSEPM -> Integral(minMassBin2, maxMassBin2);
        normMassMEPM = histMassMEPM -> Integral(minMassBin1, maxMassBin1) + histMassMEPM -> Integral(minMassBin2, maxMassBin2);

        canvasCentrScan -> cd(iCentr+1);
        gPad -> SetLogy(1);
        histMassSEPM -> SetTitle(Form("0 < #it{p}_{T} < 30 GeV/#it{c}, %1.0f - %1.0f %%", minCentrBins[iCentr], maxCentrBins[iCentr]));
        histMassSEPM -> GetXaxis() -> SetRangeUser(2, 5);
        histMassSEPM -> GetYaxis() -> SetRangeUser(1, 1e7);
        histMassSEPM -> SetMarkerStyle(20);
        histMassSEPM -> SetMarkerSize(0.5);
        histMassMEPM -> GetXaxis() -> SetRangeUser(2, 5);
        histMassSEPM -> Draw("EP");
        histMassMEPM -> Scale(normMassSEPM / normMassMEPM);
        histMassMEPM -> Draw("H SAME");

        TH1D *histBkgSubtrSEMP = (TH1D*) histMassSEPM -> Clone(Form("histMassBkgSubtrSEPM_%1.0f_%1.0f__%1.0f_%1.0f", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));
        histBkgSubtrSEMP -> Add(histMassMEPM, -1.);
        histBkgSubtrSEMP -> SetLineColor(kRed+1);
        histBkgSubtrSEMP -> SetMarkerStyle(24);
        histBkgSubtrSEMP -> SetMarkerSize(0.5);
        histBkgSubtrSEMP -> SetMarkerColor(kRed+1);
        histBkgSubtrSEMP -> Draw("EP SAME");
    }

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