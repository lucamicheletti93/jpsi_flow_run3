#if !defined(CLING) || defined(ROOTCLING)

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>

#include "TSystemDirectory.h"
#include <TLorentzVector.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TKey.h"
#include "THashList.h"
#include "TProfile.h"
#include "TTreeReader.h"
#include "TLine.h"
#include "TGaxis.h"
#include <ROOT/RDataFrame.hxx>

#endif

void LoadStyle();
void SetLegend(TLegend *);

void read_tree_fast() {
    //LoadStyle();
    float fMass, fPt, fEta, fTauz, fTauxy, fU2Q2, fCos2DeltaPhi, fR2EP, fR2SP, fCentFT0C = -99999;
    int fSign = -99999;

    const int nMassBins = 30;
    const double minMassRange = 2;
    const double maxMassRange = 5;

    TLine *lineJpsi = new TLine(3.096, 1, 3.096, 1e5);
    lineJpsi -> SetLineStyle(kDashed);
    lineJpsi -> SetLineColor(kRed+1);
    lineJpsi -> SetLineWidth(2);

    TH1F *histMass = new TH1F("histMass", "", nMassBins, minMassRange, maxMassRange);
    TH1F *histU2Q2 = new TH1F("histU2Q2", "", 1000, -1, 1);
    TH1F *histR2SP = new TH1F("histR2SP", "", 1000, -1, 1);
    TH1F *histC2DP = new TH1F("histC2DP", "", 1000, -1, 1);
    TH1F *histR2EP = new TH1F("histR2EP", "", 1000, -1, 1);

    TH2F *histMassU2Q2 = new TH2F("histMassU2Q2", "histMassU2Q2", nMassBins, minMassRange, maxMassRange, 200, -1, 1);
    TH2F *histMassR2EP = new TH2F("histMassR2EP", "histMassR2EP", nMassBins, minMassRange, maxMassRange, 200, -1, 1);
    TH2F *histMassC2DP = new TH2F("histMassC2DP", "histMassC2DP", nMassBins, minMassRange, maxMassRange, 200, -1, 1);
    TH2F *histMassR2SP = new TH2F("histMassR2SP", "histMassR2SP", nMassBins, minMassRange, maxMassRange, 200, -1, 1);
    TH2F *histMassPt = new TH2F("histMassPt", "", nMassBins, minMassRange, maxMassRange, 100, 0, 5);

    auto hprofMassU2Q2  = new TProfile("hprofMassU2Q2", "hprofMassU2Q2", 100, minMassRange, maxMassRange, -1, 1, "");
    auto hprofMassR2SP  = new TProfile("hprofMassR2SP", "hprofMassR2SP", 100, minMassRange, maxMassRange, -1, 1, "");

    auto hprofCentU2Q2  = new TProfile("hprofCentU2Q2", "hprofCentU2Q2", 100, 0, 100, -1, 1, "");
    auto hprofCentR2SP  = new TProfile("hprofCentR2SP", "hprofCentR2SP", 100, 0, 100, -1, 1, "");

    hprofMassU2Q2 -> SetLineColor(kRed+1);
    hprofMassR2SP -> SetLineColor(kRed+1);

    string pathToFile = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/data/pass2/LHC22zzh/AO2D_reduced_with_flow.root";
    TFile *fIn = new TFile(pathToFile.c_str(), "READ");
    TIter next(fIn -> GetListOfKeys()); 
    TKey *key; 
    while ((key = (TKey*) next())) { 
        TString dirName = key -> GetName();
        if (!dirName.Contains("DF_")) {
            continue;
        }

        TTree *tree = (TTree*) fIn -> Get(Form("%s/O2rtdimuonall", dirName.Data()));
        tree -> SetBranchAddress("fMass", &fMass);
        tree -> SetBranchAddress("fPt", &fPt);
        tree -> SetBranchAddress("fEta", &fEta);
        tree -> SetBranchAddress("fSign", &fSign);
        tree -> SetBranchAddress("fTauz", &fTauz);
        tree -> SetBranchAddress("fTauxy", &fTauxy);
        tree -> SetBranchAddress("fU2Q2", &fU2Q2);
        tree -> SetBranchAddress("fCos2DeltaPhi", &fCos2DeltaPhi);
        tree -> SetBranchAddress("fR2EP", &fR2EP);
        tree -> SetBranchAddress("fR2SP", &fR2SP);
        tree -> SetBranchAddress("fCentFT0C", &fCentFT0C);


        for (int iEntry = 0;iEntry < tree -> GetEntries();iEntry++) {
            tree -> GetEntry(iEntry);
            histU2Q2 -> Fill(fU2Q2);
            histR2SP -> Fill(fR2SP);
            histC2DP -> Fill(fCos2DeltaPhi);
            histR2EP -> Fill(fR2EP);

            if (fSign == 0 && TMath::Abs(fEta) > 2.5 && TMath::Abs(fEta) < 4) {
                if (fCentFT0C > 0 && fCentFT0C < 100 && fPt > 0 && fPt < 10) {
                    hprofCentU2Q2 -> Fill(fCentFT0C, fU2Q2);
                    hprofCentR2SP -> Fill(fCentFT0C, fR2SP);
                }
                if (fCentFT0C > 20 && fCentFT0C < 50 && fPt > 0 && fPt < 10) {
                    histMassU2Q2 -> Fill(fMass, fU2Q2);
                    histMassR2SP -> Fill(fMass, fR2SP);
                    histMassC2DP -> Fill(fMass, fCos2DeltaPhi);
                    histMassR2EP -> Fill(fMass, fR2EP);
                    histMassPt -> Fill(fMass, fPt);
                    histR2SP -> Fill(fR2SP);

                    hprofMassU2Q2 -> Fill(fMass, fU2Q2);
                    hprofMassR2SP -> Fill(fMass, fR2SP);
                }
            }
        }
    }
    fIn -> Close();

    TCanvas *canvasProfCent = new TCanvas("canvasProfCent", "", 1200, 1200);
    canvasProfCent -> Divide(2, 2);
    canvasProfCent -> cd(1);
    hprofCentU2Q2 -> Draw("EP");
    canvasProfCent -> cd(2);
    hprofCentR2SP -> Draw("EP");

    TCanvas *canvasDistr = new TCanvas("canvasDistr", "canvasDistr", 1200, 1200);
    canvasDistr -> Divide(2, 2);
    canvasDistr -> cd(1);
    histMassU2Q2 -> Draw("COLZ");
    canvasDistr -> cd(2);
    histMassR2EP -> Draw("COLZ");
    canvasDistr -> cd(3);
    histMassC2DP -> Draw("COLZ");
    canvasDistr -> cd(4);
    histMassR2SP -> Draw("COLZ");

    TH1F *histProjMass  = (TH1F*) histMassU2Q2 -> ProjectionX("histProjMass");

    TProfile *histProjU2Q2  = (TProfile*) histMassU2Q2 -> ProfileX("histProjU2Q2");
    TProfile *histProjR2SP  = (TProfile*) histMassR2SP -> ProfileX("histProjR2SP");

    TProfile *histProjC2DP  = (TProfile*) histMassC2DP -> ProfileX("histProjC2DP");
    TProfile *histProjR2EP  = (TProfile*) histMassR2EP -> ProfileX("histProjR2EP");

    TProfile *histProjPt  = (TProfile*) histMassPt -> ProfileX("histProjPt");

    // Compute v2 with SP
    TH1F *histV2SP = new TH1F("histV2SP", "", nMassBins, minMassRange, maxMassRange);
    TH1F *histV2EP = new TH1F("histV2EP", "", nMassBins, minMassRange, maxMassRange);

    for (int i = 0;i < nMassBins;i++) {
        // SP
        float u2q2 = histProjU2Q2 -> GetBinContent(i+1);
        float errU2q2 = histProjU2Q2 -> GetBinError(i+1);
        float r2sp = histProjR2SP -> GetBinContent(i+1);
        float errR2sp = histProjR2SP -> GetBinError(i+1);
        // EP
        float c2dp = histProjC2DP -> GetBinContent(i+1);
        float errC2dp = histProjC2DP -> GetBinError(i+1);
        float r2ep = histProjR2EP -> GetBinContent(i+1);
        float errR2ep = histProjR2EP -> GetBinError(i+1);

        if (r2sp == 0) {
            histV2SP -> SetBinContent(i+1, -999.);
            histV2SP -> SetBinError(i+1, -999.);
        } else {
            histV2SP -> SetBinContent(i+1, u2q2/r2sp);
            histV2SP -> SetBinError(i+1, (u2q2/r2sp) * TMath::Sqrt((errU2q2/u2q2)*(errU2q2/u2q2) + (errR2sp/r2sp)*(errR2sp/r2sp)));
        }

        if (r2ep == 0) {
            histV2EP -> SetBinContent(i+1, -999.);
            histV2EP -> SetBinError(i+1, -999.);
        } else {
            histV2EP -> SetBinContent(i+1, c2dp/r2ep);
            histV2EP -> SetBinError(i+1, (c2dp/r2ep) * TMath::Sqrt((errC2dp/c2dp)*(errC2dp/c2dp) + (errR2ep/r2ep)*(errR2ep/r2ep)));
        }  
    }

    TCanvas *canvasProfMassV2SP = new TCanvas("canvasProfMassV2SP", "", 1200, 1200);
    canvasProfMassV2SP -> Divide(2, 2);
    canvasProfMassV2SP -> cd(1);
    gPad -> SetLogy(1);
    histProjMass -> Draw("EP");
    lineJpsi -> Draw("SAME");
    canvasProfMassV2SP -> cd(2);
    histProjU2Q2 -> Draw("EP");
    hprofMassU2Q2 -> Draw("EP SAME");
    canvasProfMassV2SP -> cd(3);
    histV2SP -> GetYaxis() -> SetRangeUser(-1, 2);
    histV2SP -> Draw("EP");
    lineJpsi -> Draw("SAME");
    canvasProfMassV2SP -> cd(4);
    histProjR2SP -> Draw("EP");
    hprofMassR2SP -> Draw("EP SAME");

    TCanvas *canvasProfMassV2EP = new TCanvas("canvasProfMassV2EP", "", 1200, 1200);
    canvasProfMassV2EP -> Divide(2, 2);
    canvasProfMassV2EP -> cd(1);
    gPad -> SetLogy(1);
    histProjMass -> Draw("EP");
    lineJpsi -> Draw("SAME");
    canvasProfMassV2EP -> cd(2);
    histProjC2DP -> Draw("EP");
    canvasProfMassV2EP -> cd(3);
    histV2EP -> GetYaxis() -> SetRangeUser(-0.2, 0.2);
    histV2EP -> Draw("EP");
    lineJpsi -> Draw("SAME");
    canvasProfMassV2EP -> cd(4);
    histProjR2EP -> Draw("EP");

    TCanvas *canvasPr = new TCanvas("canvasPr", "", 600, 600);
    histProjPt -> Draw("EP");


    TCanvas *canvasU2Q2 = new TCanvas("canvasU2Q2", "", 600, 600);
    histU2Q2 -> Draw("EP");

    TCanvas *canvasR2SP = new TCanvas("canvasR2SP", "", 600, 600);
    histR2SP -> Draw("EP");

    TCanvas *canvasC2DP = new TCanvas("canvasC2DP", "", 600, 600);
    histC2DP -> Draw("EP");

    TCanvas *canvasR2EP = new TCanvas("canvasR2EP", "", 600, 600);
    histR2EP -> Draw("EP");

    TFile *fOut = new TFile("v2_results.root", "RECREATE");
    histProjMass -> Write();
    histProjPt -> Write();
    histV2SP -> Write();
    histV2EP -> Write();
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
    legend -> SetTextSize(0.045);
}