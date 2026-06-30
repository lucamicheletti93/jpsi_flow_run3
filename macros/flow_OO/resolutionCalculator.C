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

void LoadStyle();
void SetLegend(TLegend *);
void SetHist(TH1D *, Color_t , int , bool );
TH1D* ComputeResolution(TH2F *, TH2F *, TH2F *, string );
TCanvas* DoRatioPlot(TH1D *, std::vector<TH1D*> , string , string );

void resolutionCalculator(string train = "698032", string detector = "TPC", bool runByRun = false) {
    LoadStyle();
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

        TH2F *hR2RebinSPAB = (TH2F*) hR2SPAB->Clone("hR2RebinSPAB");
        TH2F *hR2RebinSPAC = (TH2F*) hR2SPAC->Clone("hR2RebinSPAC");
        TH2F *hR2RebinSPBC = (TH2F*) hR2SPBC->Clone("hR2RebinSPBC");

        hR2RebinSPAB->RebinX(5);
        hR2RebinSPAC->RebinX(5);
        hR2RebinSPBC->RebinX(5);

        TH1D *hSumR2SP = (TH1D*) ComputeResolution(hR2SPAB, hR2SPAC, hR2SPBC, "");
        hSumR2SP->SetMarkerStyle(20);
        hSumR2SP->SetMarkerColor(kBlack);
        hSumR2SP->SetLineColor(kBlack);

        TH1D *hSumR2RebinSP = (TH1D*) ComputeResolution(hR2RebinSPAB, hR2RebinSPAC, hR2RebinSPBC, "");
        hSumR2RebinSP->SetMarkerStyle(20);
        hSumR2RebinSP->SetMarkerColor(kBlack);
        hSumR2RebinSP->SetLineColor(kBlack);

        // Compute the dimuon weight
        THashList *hListSEPM = (THashList*) fIn->Get("analysis-same-event-pairing/output");
        TList *listSEPM = (TList*) hListSEPM->FindObject("PairsMuonSEPM_muonQualityCutsMUONStandalone");
        THnSparseF *histSEPM = (THnSparseF*) listSEPM->FindObject("Mass_Pt_centrFT0C_V2");

        //double minMassBin = histSEPM->GetAxis(0)->FindBin(2.1);
        //double maxMassBin = histSEPM->GetAxis(0)->FindBin(4.9);
        //histSEPM->GetAxis(0)->SetRange(minMassBin, maxMassBin);
        TH1D *hWeight = (TH1D*) histSEPM->Projection(3, "hWeight");

        // Extract the resolution in the required centrality bins
        double minCentrBins[] = {0, 20, 60};
        double maxCentrBins[] = {20, 60, 90};
        double centrBins[] = {0, 20, 60, 90};
        TH1D *hMeanR2SP = new TH1D("hMeanR2SP", ";Centrality FT0C(%);#it{R}_{2}^{SP}", 3, centrBins);
        hMeanR2SP->SetMarkerStyle(20);
        hMeanR2SP->SetMarkerColor(kRed);
        hMeanR2SP->SetLineColor(kRed);

        TH1D *hWmeanR2SP = new TH1D("hWmeanR2SP", ";Centrality FT0C(%);#it{R}_{2}^{SP}", 3, centrBins);
        hWmeanR2SP->SetMarkerStyle(25);
        hWmeanR2SP->SetMarkerColor(kRed);
        hWmeanR2SP->SetLineColor(kRed);

        TH1D *hDimuWmeanR2SP = new TH1D("hDimuWmeanR2SP", ";Centrality FT0C(%);#it{R}_{2}^{SP}", 3, centrBins);
        hDimuWmeanR2SP->SetMarkerStyle(24);
        hDimuWmeanR2SP->SetMarkerColor(kRed);
        hDimuWmeanR2SP->SetLineColor(kRed);


        for (int iCentrBin = 0;iCentrBin < 3;iCentrBin++) {
            int binBegin = hSumR2RebinSP->FindBin(minCentrBins[iCentrBin]);
            int binEnd = hSumR2RebinSP->FindBin(maxCentrBins[iCentrBin]);

            double meanR2SP = 0, wMeanR2SP = 0, wSum =0, dimuWmeanR2SP = 0, dimuWsum = 0;
            for (int iBin = binBegin;iBin < binEnd;iBin++) {
                meanR2SP += hSumR2RebinSP->GetBinContent(iBin);

                wMeanR2SP += hSumR2RebinSP->GetBinContent(iBin) / (hSumR2RebinSP->GetBinError(iBin) * hSumR2RebinSP->GetBinError(iBin));
                wSum += 1 / (hSumR2RebinSP->GetBinError(iBin) * hSumR2RebinSP->GetBinError(iBin));

                dimuWmeanR2SP += hSumR2RebinSP->GetBinContent(iBin)*hWeight->GetBinContent(iBin);
                dimuWsum += hWeight->GetBinContent(iBin);

                //std::cout << hSumR2RebinSP->GetBinContent(iBin) << " " << hWeight->GetBinContent(iBin) << std::endl;
            }
            meanR2SP = meanR2SP/(binEnd-binBegin);
            hMeanR2SP->SetBinContent(iCentrBin+1, meanR2SP);
            hMeanR2SP->SetBinError(iCentrBin+1, 0);

            wMeanR2SP = wMeanR2SP/wSum;
            hWmeanR2SP->SetBinContent(iCentrBin+1, wMeanR2SP);
            hWmeanR2SP->SetBinError(iCentrBin+1, 1./TMath::Sqrt(wSum));

            dimuWmeanR2SP = dimuWmeanR2SP/dimuWsum;
            hDimuWmeanR2SP->SetBinContent(iCentrBin+1, dimuWmeanR2SP);
            hDimuWmeanR2SP->SetBinError(iCentrBin+1, 1./TMath::Sqrt(dimuWsum));
            std::cout << "mean = " << meanR2SP << " ; w. mean = " << wMeanR2SP << " ; dimu. w. mean = " << dimuWmeanR2SP << std::endl; 
        }

        TCanvas *canvasR2SP = new TCanvas("canvasR2SP", "", 800, 600);
        hSumR2SP->GetYaxis()->SetRangeUser(0, 0.7);
        hSumR2SP->Draw("EP");

        TCanvas *canvasR2RebinSP = new TCanvas("canvasR2RebinSP", "", 800, 600);
        hSumR2RebinSP->GetYaxis()->SetRangeUser(0, 0.7);
        hSumR2RebinSP->Draw("EP");
        hMeanR2SP->Draw("EP SAME");
        hWmeanR2SP->Draw("EP SAME");
        hDimuWmeanR2SP->Draw("EP SAME");

        TLegend *legendR2RebinSP = new TLegend(0.55, 0.68, 0.85, 0.88, " ", "brNDC");
        SetLegend(legendR2RebinSP);
        legendR2RebinSP -> SetTextSize(0.035);
        legendR2RebinSP -> AddEntry(hSumR2RebinSP, "Resolution", "PL");
        legendR2RebinSP -> AddEntry(hMeanR2SP, "Normal average", "PL");
        legendR2RebinSP -> AddEntry(hWmeanR2SP, "Weighted average", "PL");
        legendR2RebinSP -> AddEntry(hDimuWmeanR2SP, "Weighted average with #mu^{+}#mu^{-}", "PL");
        legendR2RebinSP -> Draw();
        
        TFile *fOut = new TFile(Form("output/resolution_%s_train_%s.root", detector.c_str(), train.c_str()), "RECREATE");
        hSumR2SP->Write("hR2SP_original");
        hSumR2RebinSP->Write("hR2SP_original_rebin");
        hMeanR2SP->Write("hMeanR2SP");
        hWmeanR2SP->Write("hWmeanR2SP");
        hDimuWmeanR2SP->Write("hDimuWmeanR2SP");
        fOut->Close();

        canvasR2RebinSP->SaveAs(Form("figures/resolution_%s_train_%s.pdf", detector.c_str(), train.c_str()));
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetHist(TH1D *hist, Color_t color, int size, bool scale = false) {
    hist->SetMarkerStyle(size);
    hist->SetMarkerColor(color);
    hist->SetMarkerSize(0.5);
    hist->SetLineColor(color);
    if (scale) {hist->Scale(1./hist->Integral());}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
    histSum->GetYaxis()->SetRangeUser(0.02, 0.7);
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetLegend(TLegend *legend) {
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.045);
}