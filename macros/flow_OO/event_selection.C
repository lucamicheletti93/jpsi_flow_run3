#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <string>

#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <THnSparse.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TGrid.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TRandom3.h>

void SetLegend(TLegend *);
void SetHist(TH1D *hist, Color_t color, int size) {
    hist->SetMarkerStyle(size);
    hist->SetMarkerColor(color);
    hist->SetMarkerSize(0.5);
    hist->SetLineColor(color);
    hist->Scale(1./hist->Integral());
}
TCanvas* DoRatioPlot(TH1D *histSum, std::vector<TH1D*> hist, string canvasName, string varName) {
    const int nRuns = int(hist.size());
    gStyle->SetPalette(kRainBow);
    gStyle->SetOptStat(false);

    int firstBin = histSum->FindFirstBinAbove(0.0);
    int lastBin  = histSum->FindLastBinAbove(0.0);
    double xRangeMin = histSum->GetBinLowEdge(firstBin);
    double xRangeMax = histSum->GetBinCenter(lastBin) + histSum->GetBinWidth(lastBin);

    TCanvas *canvas = new TCanvas(canvasName.c_str(), "", 800, 600);

    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);

    pad1->SetBottomMargin(0);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);

    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    //gPad->SetLogy(true);
    SetHist(histSum, kBlack, 20);
    histSum->GetXaxis()->SetRangeUser(xRangeMin, xRangeMax);
    histSum->GetYaxis()->SetRangeUser(0.005, 0.020);
    histSum->SetTitle("");
    histSum->SetName("Sum");
    histSum->Draw("EP");
    for (int iRun = 0;iRun < nRuns;iRun++) {
        hist[iRun]->Scale(1./hist[iRun]->Integral());
        hist[iRun]->Draw("SAME PLC PMC");
    }
    
    TLegend *leg = (TLegend*)gPad->BuildLegend();
    leg->SetTextSize(0.04);
    leg->SetX1(0.70);
    leg->SetY1(0.50);
    leg->SetX2(0.90);
    leg->SetY2(0.90);

    gPad->Modified();
    gPad->Update();

    pad2->cd();
    std::vector<TH1D*> histRatio(nRuns);
    for (int iRun = 0;iRun < nRuns;iRun++) {
        histRatio[iRun] = (TH1D*) hist[iRun]->Clone(Form("histRatio_%d", iRun));

        if (iRun == 0) {
            histRatio[iRun]->SetTitle("");
            histRatio[iRun]->GetYaxis()->SetTitle("Data/MC");
            histRatio[iRun]->GetXaxis()->SetTitle(varName.c_str());
            histRatio[iRun]->GetXaxis()->SetLabelSize(0.08);
            histRatio[iRun]->GetXaxis()->SetTitleSize(0.1);
            histRatio[iRun]->GetYaxis()->SetNdivisions(505);
            histRatio[iRun]->GetYaxis()->SetTitleSize(0.1);
            histRatio[iRun]->GetYaxis()->SetLabelSize(0.08);
        }

        histRatio[iRun]->SetTitle("");
        histRatio[iRun]->GetXaxis()->SetRangeUser(xRangeMin, xRangeMax);
        histRatio[iRun]->Divide(histSum);
        histRatio[iRun]->GetYaxis()->SetRangeUser(0.5, 1.5);
        histRatio[iRun]->Draw("SAME PLC PMC");
    }
    TLine *line = new TLine(xRangeMin, 1, xRangeMax, 1);
    line->SetLineStyle(2);
    line->Draw("SAME");

    canvas->cd();

    return canvas;
}

void event_selection() {
    string pathIn = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/LHC25ae/train_687739";
    string muonCutTm = "muonMinimalCuts10SigmaPDCA";
    string muonCutTr = "muonQualityCutsMUONStandalone";
    string histBlock = "Mass_Pt_centrFT0C_V2";

    TH1D *hSumVtxZ, *hSumCentFT0C, *hSumCentFT0M;
    TH2D *hSumPsi2A_CentFT0C, *hSumPsi2B_CentFT0C, *hSumPsi2C_CentFT0C;
    TH1D *hSumPt, *hSumEta, *hSumPhi, *hSumpdca, *hSumRAtAbsorberEnd, *hSumChi2MCHMID, *hSumChi2MCHMFT;

    std::vector<int> runList;
    int entry;
    std::ifstream fInRunList(Form("%s/run_list.txt", pathIn.c_str()));
    while (fInRunList >> entry) {runList.push_back(entry);}
    const int nRuns = int(runList.size());
    std::cout << nRuns << std::endl;

    int iRun = 0;

    std::vector<TH1D*> hVtxZ(nRuns);
    std::vector<TH1D*> hCentFT0C(nRuns);
    std::vector<TH1D*> hCentFT0M(nRuns);
    std::vector<TH2D*> hPsi2A_CentFT0C(nRuns);
    std::vector<TH2D*> hPsi2B_CentFT0C(nRuns);
    std::vector<TH2D*> hPsi2C_CentFT0C(nRuns);
    std::vector<TH1D*> hPt(nRuns);
    std::vector<TH1D*> hEta(nRuns);
    std::vector<TH1D*> hPhi(nRuns);
    std::vector<TH1D*> hpdca(nRuns);
    std::vector<TH1D*> hRAtAbsorberEnd(nRuns);
    std::vector<TH1D*> hChi2MCHMID(nRuns);
    std::vector<TH1D*> hChi2MCHMFT(nRuns);

    for (auto run : runList) {
        std::cout << "-------- Processing run " << run << " --------" << std::endl;
        TFile *fIn = new TFile(Form("%s/%d/AnalysisResults.root", pathIn.c_str(), run));

        TList *hlistTrEvents = (TList*) fIn->Get("analysis-event-selection/output");
        TList *listTrEventsAfterCuts = (TList*) hlistTrEvents->FindObject("Event_AfterCuts");
        hVtxZ[iRun] = (TH1D*) listTrEventsAfterCuts->FindObject("VtxZ");
        hCentFT0C[iRun] = (TH1D*) listTrEventsAfterCuts->FindObject("CentFT0C");
        hCentFT0M[iRun] = (TH1D*) listTrEventsAfterCuts->FindObject("CentFT0M");
        hPsi2A_CentFT0C[iRun] = (TH2D*) listTrEventsAfterCuts->FindObject("Psi2A_CentFT0C");
        hPsi2B_CentFT0C[iRun] = (TH2D*) listTrEventsAfterCuts->FindObject("Psi2B_CentFT0C");
        hPsi2C_CentFT0C[iRun] = (TH2D*) listTrEventsAfterCuts->FindObject("Psi2C_CentFT0C");

        if (iRun == 0) {
            hSumVtxZ = (TH1D*) hVtxZ[iRun]->Clone("hSumVtxZ");
            hSumCentFT0C = (TH1D*) hCentFT0C[iRun]->Clone("hSumCentFT0C");
            hSumCentFT0M = (TH1D*) hCentFT0M[iRun]->Clone("hSumCentFT0M");
            hSumPsi2A_CentFT0C = (TH2D*) hPsi2A_CentFT0C[iRun]->Clone("hSumPsi2A_CentFT0C");
            hSumPsi2B_CentFT0C = (TH2D*) hPsi2B_CentFT0C[iRun]->Clone("hSumPsi2B_CentFT0C");
            hSumPsi2C_CentFT0C = (TH2D*) hPsi2C_CentFT0C[iRun]->Clone("hSumPsi2C_CentFT0C");
            /*hSumPt = (TH1D*) hPt[iRun]->Clone("hSumPt");
            hSumEta = (TH1D*) hEta[iRun]->Clone("hSumEta");
            hSumPhi = (TH1D*) hPhi[iRun]->Clone("hSumPhi");
            hSumpdca = (TH1D*) hpdca[iRun]->Clone("hSumpdca");
            hSumRAtAbsorberEnd = (TH1D*) hRAtAbsorberEnd[iRun]->Clone("hSumRAtAbsorberEnd");
            hSumChi2MCHMID = (TH1D*) hChi2MCHMID[iRun]->Clone("hSumChi2MCHMID");
            hSumChi2MCHMFT = (TH1D*) hChi2MCHMFT[iRun]->Clone("hSumChi2MCHMFT");*/
        } else {
            hSumVtxZ->Add(hVtxZ[iRun]);
            hSumCentFT0C->Add(hCentFT0C[iRun]);
            hSumCentFT0M->Add(hCentFT0M[iRun]);
            hSumPsi2A_CentFT0C->Add(hPsi2A_CentFT0C[iRun]);
            hSumPsi2B_CentFT0C->Add(hPsi2B_CentFT0C[iRun]);
            hSumPsi2C_CentFT0C->Add(hPsi2C_CentFT0C[iRun]);
            /*hSumPt->Add(hPt[iRun]);
            hSumEta->Add(hEta[iRun]);
            hSumPhi->Add(hPhi[iRun]);
            hSumpdca->Add(hpdca[iRun]);
            hSumRAtAbsorberEnd->Add(hRAtAbsorberEnd[iRun]);
            hSumChi2MCHMID->Add(hChi2MCHMID[iRun]);
            hSumChi2MCHMFT->Add(hChi2MCHMFT[iRun]);*/
        }

        hVtxZ[iRun]->SetTitle(Form("%d", run));
        hCentFT0C[iRun]->SetTitle(Form("%d", run));
        hCentFT0M[iRun]->SetTitle(Form("%d", run));
        hPsi2A_CentFT0C[iRun]->SetTitle(Form("%d", run));
        hPsi2B_CentFT0C[iRun]->SetTitle(Form("%d", run));
        hPsi2C_CentFT0C[iRun]->SetTitle(Form("%d", run));
    
        iRun++;
    }

    
    TLatex *latexTitle = new TLatex();
    latexTitle -> SetTextSize(0.040);
    latexTitle -> SetNDC();
    latexTitle -> SetTextFont(42);

    TCanvas *canvasVtxZ = DoRatioPlot(hSumVtxZ, hVtxZ, "canvasVtxZ", "vtx_{Z} (cm)");
    TCanvas *canvasCentFT0C = DoRatioPlot(hSumCentFT0C, hCentFT0C, "canvasCentFT0C", "Centr FT0C (%)");
    canvasCentFT0C->cd();
    latexTitle->DrawLatex(0.15, 0.87, "OO #sqrt{#it{s}_{NN}} = 5.36 TeV, LHC25ae pass2");
    canvasCentFT0C-> Update();

    canvasVtxZ->SaveAs(Form("%s/figures/VtxZ.pdf", pathIn.c_str()));
    canvasCentFT0C->SaveAs(Form("%s/figures/CentFT0C.pdf", pathIn.c_str()));
}