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

void project_histograms() {
    bool integrated = false;
    string pathToFin = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/LHC25ae";
    string muonCut = "matchedMchMid";
    string suffix = "";
    string methodName = "SP";
    int method = methodName == "SP" ? 4 : 5;

    TFile *fInSE  = new TFile("/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/LHC25ae/AnalysisResults_654169.root", "READ");
    TList *listSE = (TList*) fInSE -> Get("analysis-same-event-pairing/output");

    TList *listSEPM = (TList*) listSE -> FindObject(Form("PairsMuonSEPM_%s", muonCut.c_str()));
    THnSparseF* histSEPM = (THnSparseF*) listSEPM->FindObject("Mass_Pt_centrFT0C_V2");

    TList *listSEPP = (TList*) listSE -> FindObject(Form("PairsMuonSEPP_%s", muonCut.c_str()));
    THnSparseF* histSEPP = (THnSparseF*) listSEPP -> FindObject("Mass_Pt_centrFT0C_V2");

    TList *listSEMM = (TList*) listSE -> FindObject(Form("PairsMuonSEMM_%s", muonCut.c_str()));
    THnSparseF* histSEMM = (THnSparseF*) listSEMM -> FindObject("Mass_Pt_centrFT0C_V2");

    int nCentrBins = 1;
    double minCentrBins[] = {0};
    double maxCentrBins[] = {10};

    const int nPtBins = 8;
    double minPtBins[] = {0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0};
    double maxPtBins[] = {1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 15.0};

    TFile *fOut = new TFile(Form("%s/Histograms_OO_%s_Centrality_%1.f_%1.f_%s.root", pathToFin.c_str(), muonCut.c_str(), minCentrBins[0], maxCentrBins[0], methodName.c_str()), "RECREATE");
    for (int iCentr = 0;iCentr < nCentrBins;iCentr++) {
        for (int iPt = 0;iPt < nPtBins;iPt++) {
            TH1D *histMassSEPM = (TH1D*) projectHistogram(histSEPM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr]);
            histMassSEPM -> SetName(Form("histMassSEPM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr]));
            TH1D *histMassSEPP = (TH1D*) projectHistogram(histSEPP, minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr]);
            histMassSEPP -> SetName(Form("histMassSEPP_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr]));
            TH1D *histMassSEMM = (TH1D*) projectHistogram(histSEMM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr]);
            histMassSEMM -> SetName(Form("histMassSEMM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr]));

            TProfile *histV2SEPM = (TProfile*) projectProfile(histSEPM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr], method);
            histV2SEPM -> SetName(Form("histV2SEPM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr]));
	        TProfile *histV2SEPP = (TProfile*) projectProfile(histSEPP, minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr], method);
            histV2SEPP -> SetName(Form("histV2SEPP_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr]));
	        TProfile *histV2SEMM = (TProfile*) projectProfile(histSEMM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr], method);
            histV2SEMM -> SetName(Form("histV2SEMM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr]));
	    
            TProfile *histPtSEPM = (TProfile*) projectProfile(histSEPM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr], 1);
            histPtSEPM -> SetName(Form("histPtSEPM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr]));

            histMassSEPM -> SetTitle(Form("SEPM - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr]));
            histMassSEPP -> SetTitle(Form("SEPP - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr]));
            histMassSEMM -> SetTitle(Form("SEMM - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr]));
            histV2SEPM -> SetTitle(Form("V2 SEPM - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr]));
	        histV2SEPP -> SetTitle(Form("V2 SEPP - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr]));
	        histV2SEMM -> SetTitle(Form("V2 SEMM - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr]));
            histPtSEPM -> SetTitle(Form("<pT> SEPM - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[iCentr], maxCentrBins[iCentr]));

            fOut -> cd();
            histMassSEPM -> Write();
	        histMassSEPP -> Write();
	        histMassSEMM -> Write();   

	        histV2SEPM -> Write();
            histV2SEPP -> Write();
            histV2SEMM -> Write();

            histPtSEPM -> Write();
        }
    }
    fOut -> Close();
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
    double minCentrBin = histSparse -> GetAxis(3) -> FindBin(minCentrRange);
    double maxCentrBin = histSparse -> GetAxis(3) -> FindBin(maxCentrRange - 0.01);

    Printf("minPtBin = %0.1f, maxPtBin = %0.1f, minCentrBin = %0.1f, maxCentrBin = %0.1f", minPtBin, maxPtBin, minCentrBin, maxCentrBin);
    histSparse -> GetAxis(1) -> SetRange(minPtBin, maxPtBin);
    histSparse -> GetAxis(3) -> SetRange(minCentrBin, maxCentrBin);

    TH1D *histProj = (TH1D*) histSparse -> Projection(0, Form("histProj_%0.1f_%0.1f__%1.0f_%1.0f", minPtRange, maxPtRange, minCentrRange, maxCentrRange));
    return histProj;
}

TProfile* projectProfile(THnSparseF *histSparse, double minPtRange, double maxPtRange, double minCentrRange, double maxCentrRange, int var) {
    double minPtBin = histSparse -> GetAxis(1) -> FindBin(minPtRange);
    double maxPtBin = histSparse -> GetAxis(1) -> FindBin(maxPtRange - 0.01);
    double minCentrBin = histSparse -> GetAxis(3) -> FindBin(minCentrRange);
    double maxCentrBin = histSparse -> GetAxis(3) -> FindBin(maxCentrRange - 0.01);

    Printf("minPtBin = %0.1f, maxPtBin = %0.1f, minCentrBin = %0.1f, maxCentrBin = %0.1f", minPtBin, maxPtBin, minCentrBin, maxCentrBin);
    histSparse -> GetAxis(1) -> SetRange(minPtBin, maxPtBin);
    histSparse -> GetAxis(3) -> SetRange(minCentrBin, maxCentrBin);

    TH2D *hist2DProj = (TH2D*) histSparse -> Projection(var, 0, "Projection");
    hist2DProj -> Rebin2D(2);
    TProfile *histProj = (TProfile*) hist2DProj -> ProfileX(Form("histProj_%i_%0.1f_%0.1f__%1.0f_%1.0f", var, minPtRange, maxPtRange, minCentrRange, maxCentrRange));
    return histProj;
}
