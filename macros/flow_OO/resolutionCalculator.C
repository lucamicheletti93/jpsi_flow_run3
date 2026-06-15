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
void SetHist(TH1D *hist, Color_t color, int size, bool scale = false) {
    hist->SetMarkerStyle(size);
    hist->SetMarkerColor(color);
    hist->SetMarkerSize(0.5);
    hist->SetLineColor(color);
    if (scale) {hist->Scale(1./hist->Integral());}
}
TH1D* ComputeResolution(TH2F *hR2SPAB, TH2F *hR2SPAC, TH2F *hR2SPBC, string histName) {
    TAxis *axis = (TAxis*) hR2SPAB->GetXaxis();
    int nCentrBins = axis->GetNbins();
    double *centrBins = new double[nCentrBins + 1];
    axis->GetLowEdge(centrBins);
    centrBins[nCentrBins] = axis->GetBinUpEdge(nCentrBins);

    TH1D *hR2SP = new TH1D(Form("%s", histName.c_str()), Form("%s;Centrality FT0C(%%);#it{R}_{2}^{SP}", histName.c_str()), nCentrBins, centrBins);

    for (int iBin = 0; iBin < nCentrBins; iBin++) {
        hR2SPAB->GetXaxis()->SetRangeUser(centrBins[iBin], centrBins[iBin+1]);
        hR2SPAC->GetXaxis()->SetRangeUser(centrBins[iBin], centrBins[iBin+1]);
        hR2SPBC->GetXaxis()->SetRangeUser(centrBins[iBin], centrBins[iBin+1]);

        double R2SPAB = hR2SPAB->GetMean(2);
        double R2SPAC = hR2SPAC->GetMean(2);
        double R2SPBC = hR2SPBC->GetMean(2);
        double R2SPABe = hR2SPAB->GetMean(12);
        double R2SPACe = hR2SPAC->GetMean(12);
        double R2SPBCe = hR2SPBC->GetMean(12);

        double R22SP = R2SPBC != 0 ? R2SPAB * R2SPAC / R2SPBC : 0.0;
        double R2SP = R22SP > 0 ? TMath::Sqrt(R22SP) : 0.0;
        double R2SPe =
            R2SPAB * R2SPAC * R2SPBC == 0
                ? 0.0
                : TMath::Sqrt(
                    1. / 4 * (R2SPAC / (R2SPAB * R2SPBC)) * pow(R2SPABe, 2.) +
                    1. / 4 * (R2SPAB / (R2SPAC * R2SPBC)) * pow(R2SPACe, 2.) +
                    1. / 4 * R2SPAC * R2SPAB / pow(R2SPBC, 3.) *
                        pow(R2SPBCe, 2.));
        R2SPe = isnan(R2SPe) || isinf(R2SPe) ? 0. : R2SPe;

        hR2SP->SetBinContent(iBin+1, R2SP);
        hR2SP->SetBinError(iBin+1, R2SPe);
    }
    return hR2SP;
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
    SetHist(histSum, kBlack, 20);
    histSum->GetXaxis()->SetRangeUser(xRangeMin, xRangeMax);
    histSum->GetXaxis()->SetTitleSize(0.045);
    histSum->GetXaxis()->SetLabelSize(0.045);
    histSum->GetYaxis()->SetRangeUser(0.02, 0.5);
    histSum->GetYaxis()->SetTitleSize(0.045);
    histSum->GetYaxis()->SetLabelSize(0.045);

    histSum->SetTitle("");
    histSum->SetName("Sum");
    histSum->Draw("EP");
    for (int iRun = 0;iRun < nRuns;iRun++) {
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
            histRatio[iRun]->GetXaxis()->SetLabelSize(0.1);
            histRatio[iRun]->GetXaxis()->SetTitleSize(0.1);
            histRatio[iRun]->GetYaxis()->SetNdivisions(505);
            histRatio[iRun]->GetYaxis()->SetTitleSize(0.1);
            histRatio[iRun]->GetYaxis()->SetLabelSize(0.1);
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


void resolutionCalculator(string train = "687739", string detector = "TPC", bool runByRun = false) {
    
    string estimators[3];
    if (detector == "TPCPOS") {
        estimators[0] = "FT0ATPCPOS";
        estimators[1] = "FT0CTPCPOS";
        estimators[2] = "FT0AFT0C";
    } else if (detector == "TPCNEG") {
        estimators[0] = "FT0ATPCNEG";
        estimators[1] = "FT0CTPCNEG";
        estimators[2] = "FT0AFT0C";
    } else if (detector == "TPC") {
        estimators[0] = "TPCFT0A";
        estimators[1] = "TPCFT0C";
        estimators[2] = "FT0AFT0C";
    } else if (detector == "FT0A") {
        estimators[0] = "TPCFT0A";
        estimators[1] = "FT0AFT0C";
        estimators[2] = "TPCFT0C";
    } else if (detector == "FT0C") {
        estimators[0] = "TPCFT0C";
        estimators[1] = "FT0AFT0C";
        estimators[2] = "TPCFT0A";
    } else {
        std::cout << "No detector selected. Exiting..." << std::endl;
        return;
    }

    if (runByRun) {
        // Run by run results
        string pathIn = Form("/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/LHC25ae/train_%s", train.c_str());
        std::vector<int> runList;
        int entry;
        std::ifstream fInRunList(Form("%s/run_list.txt", pathIn.c_str()));
        while (fInRunList >> entry) {runList.push_back(entry);}
        const int nRuns = int(runList.size());
        int iRun = 0;

        std::vector<TH1D*> hR2SP(nRuns);
        TH1D *hSumR2SP;
        TH2F *hSumR2SPAB, *hSumR2SPAC, *hSumR2SPBC;

        for (auto run : runList) {
            TFile *fIn = new TFile(Form("%s/%d/AnalysisResults.root", pathIn.c_str(), run));

            THashList *hList = (THashList*) fIn->Get("analysis-event-selection/output");
            TList *subList = (TList*) hList->FindObject("Event_AfterCuts");
            TH2F *hR2SPAB = (TH2F*) subList->FindObject(Form("R2SP_%s_CentFT0C", estimators[0].c_str()));
            TH2F *hR2SPAC = (TH2F*) subList->FindObject(Form("R2SP_%s_CentFT0C", estimators[1].c_str()));
            TH2F *hR2SPBC = (TH2F*) subList->FindObject(Form("R2SP_%s_CentFT0C", estimators[2].c_str()));

            if (iRun == 0) {
                hSumR2SPAB = (TH2F*) hR2SPAB->Clone("hSumR2SPAB");
                hSumR2SPAC = (TH2F*) hR2SPAC->Clone("hSumR2SPAC");
                hSumR2SPBC = (TH2F*) hR2SPBC->Clone("hSumR2SPBC");
            } else {
                hSumR2SPAB->Add(hR2SPAB);
                hSumR2SPAC->Add(hR2SPAC);
                hSumR2SPBC->Add(hR2SPBC);
            }

            hR2SP[iRun] = (TH1D*) ComputeResolution(hR2SPAB, hR2SPAC, hR2SPBC, std::to_string(run));
            iRun++;
        }
        hSumR2SP = (TH1D*) ComputeResolution(hSumR2SPAB, hSumR2SPAC, hSumR2SPBC, "sum");

        TCanvas *canvasR2SP = DoRatioPlot(hSumR2SP, hR2SP, "canvasR2SP", "Centrality FT0C (%)");

        TFile *fOut = new TFile(Form("output/runbyrun_resolution_%s.root", detector.c_str()), "RECREATE");
        hSumR2SP->Write();
        fOut->Close();

        canvasR2SP->SaveAs(Form("figures/runbyrun_resolution_%s.pdf", detector.c_str()));
    } else {
        string pathIn = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/LHC25ae";
        TFile *fIn = TFile::Open(Form("%s/AnalysisResults_%s.root", pathIn.c_str(), train.c_str()));

        THashList *hList = (THashList*) fIn->Get("analysis-event-selection/output");
        TList *subList = (TList*) hList->FindObject("Event_AfterCuts");
        TH2F *hR2SPAB = (TH2F*) subList->FindObject(Form("R2SP_%s_CentFT0C", estimators[0].c_str()));
        TH2F *hR2SPAC = (TH2F*) subList->FindObject(Form("R2SP_%s_CentFT0C", estimators[1].c_str()));
        TH2F *hR2SPBC = (TH2F*) subList->FindObject(Form("R2SP_%s_CentFT0C", estimators[2].c_str()));

        TH1D *hSumR2SP = (TH1D*) ComputeResolution(hR2SPAB, hR2SPAC, hR2SPBC, "sum");

        TCanvas *canvasR2SP = new TCanvas("canvasR2SP", "", 800, 600);
        hSumR2SP->Draw("EP");

        TFile *fOut = new TFile(Form("output/resolution_%s_train_%s.root", detector.c_str(), train.c_str()), "RECREATE");
        hSumR2SP->Write();
        fOut->Close();

        canvasR2SP->SaveAs(Form("figures/resolution_%s_train_%s.pdf", detector.c_str(), train.c_str()));

    }













    /*TFile *fIn = TFile::Open(Form("/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/LHC25ae/AnalysisResults_%s.root", train.c_str()));

    THashList *hList = (THashList*) fIn->Get("analysis-event-selection/output");
    TList *subList = (TList*) hList->FindObject("Event_AfterCuts");
    TH2F *hs_R2SPAB = (TH2F*) subList->FindObject(Form("R2SP_%s_CentFT0C", estimators[0].c_str()));
    TH2F *hs_R2SPAC = (TH2F*) subList->FindObject(Form("R2SP_%s_CentFT0C", estimators[1].c_str()));
    TH2F *hs_R2SPBC = (TH2F*) subList->FindObject(Form("R2SP_%s_CentFT0C", estimators[2].c_str()));

    TAxis *axis = (TAxis*) hs_R2SPAB->GetXaxis();
    int nCentrBins = axis->GetNbins();
    double *centrBins = new double[nCentrBins + 1];
    axis->GetLowEdge(centrBins);
    centrBins[nCentrBins] = axis->GetBinUpEdge(nCentrBins);

    TH1D *hist_r2sp = new TH1D("R2SP_Cent", "R_{2}^{SP}", nCentrBins, centrBins);
    hist_r2sp->GetXaxis()->SetTitle("Centrality FT0C(%)");
    hist_r2sp->GetYaxis()->SetTitle("R_{2}{SP}");

    for (int iBin = 0; iBin < nCentrBins; iBin++) {
        std::cout << centrBins[iBin] << " " << centrBins[iBin+1] << std::endl;
        hs_R2SPAB->GetXaxis()->SetRangeUser(centrBins[iBin], centrBins[iBin+1]);
        hs_R2SPAC->GetXaxis()->SetRangeUser(centrBins[iBin], centrBins[iBin+1]);
        hs_R2SPBC->GetXaxis()->SetRangeUser(centrBins[iBin], centrBins[iBin+1]);

        double R2SPAB = hs_R2SPAB->GetMean(2);
        double R2SPAC = hs_R2SPAC->GetMean(2);
        double R2SPBC = hs_R2SPBC->GetMean(2);
        double R2SPABe = hs_R2SPAB->GetMean(12);
        double R2SPACe = hs_R2SPAC->GetMean(12);
        double R2SPBCe = hs_R2SPBC->GetMean(12);

        double R22SP = R2SPBC != 0 ? R2SPAB * R2SPAC / R2SPBC : 0.0;
        //double R22SP = R2SPBC != 0 ? R2SPBC * R2SPAC / R2SPAB : 0.0;
        double R2SP = R22SP > 0 ? TMath::Sqrt(R22SP) : 0.0;
        double R2SPe =
            R2SPAB * R2SPAC * R2SPBC == 0
                ? 0.0
                : TMath::Sqrt(
                      1. / 4 * (R2SPAC / (R2SPAB * R2SPBC)) * pow(R2SPABe, 2.) +
                      1. / 4 * (R2SPAB / (R2SPAC * R2SPBC)) * pow(R2SPACe, 2.) +
                      1. / 4 * R2SPAC * R2SPAB / pow(R2SPBC, 3.) *
                          pow(R2SPBCe, 2.));
        R2SPe = isnan(R2SPe) || isinf(R2SPe) ? 0. : R2SPe;

        std::cout << R2SP << " +/- " << R2SPe << std::endl;

        hist_r2sp->SetBinContent(iBin+1, R2SP);
        hist_r2sp->SetBinError(iBin+1, R2SPe);
    }
    hist_r2sp->Draw();

    TFile *fOut = new TFile(Form("output/resolution_%s.root", detector.c_str()), "RECREATE");
    hist_r2sp->Write();
    fOut->Close();*/
}