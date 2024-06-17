void SetProfile(TProfile *, int , int , int , double );

void spectator_plane() {
    string pathToFiles = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/data/pass3/LHC23zzh_small/train_225020_spectator_plane";
    TFile *fIn = new TFile(Form("%s/AnalysisResults.root", pathToFiles.c_str()));

    TList *list1 = (TList*) fIn -> Get("d-q-event-qvector/outputQA");
    TList *list2 = (TList*) list1 -> FindObject("Event_AfterCuts_centralFW");

    THnSparseF* histSparseZNA = (THnSparseF*) list2 -> FindObject("Q1ZNAX_Q1ZNAY_CentFT0C");
    THnSparseF* histSparseZNC = (THnSparseF*) list2 -> FindObject("Q1ZNCX_Q1ZNCY_CentFT0C");

    TH2D *histQ1ZNAXCentFT0C = (TH2D*) histSparseZNA -> Projection(0, 2, "");
    TH2D *histQ1ZNAYCentFT0C = (TH2D*) histSparseZNA -> Projection(1, 2, "");
    TH2D *histQ1ZNCXCentFT0C = (TH2D*) histSparseZNC -> Projection(0, 2, "");
    TH2D *histQ1ZNCYCentFT0C = (TH2D*) histSparseZNC -> Projection(1, 2, "");

    TH2D *histQ1ZNACXXCentFT0C = (TH2D*)list2 -> FindObject("Q1ZNACXX_CentFT0C");
    TH2D *histQ1ZNACYYCentFT0C = (TH2D*)list2 -> FindObject("Q1ZNACYY_CentFT0C");
    TH2D *histQ1ZNACYXCentFT0C = (TH2D*)list2 -> FindObject("Q1ZNACYX_CentFT0C");
    TH2D *histQ1ZNACXYCentFT0C = (TH2D*)list2 -> FindObject("Q1ZNACXY_CentFT0C");

    TH2D *histIntercalibZNACentFT0C = (TH2D*)list2 -> FindObject("IntercalibZNA_CentFT0C");
    TH2D *histIntercalibZNCCentFT0C = (TH2D*)list2 -> FindObject("IntercalibZNC_CentFT0C");

    histQ1ZNAXCentFT0C -> GetYaxis() -> SetRangeUser(-2, 2);
    histQ1ZNAYCentFT0C -> GetYaxis() -> SetRangeUser(-2, 2);
    histQ1ZNCXCentFT0C -> GetYaxis() -> SetRangeUser(-2, 2);
    histQ1ZNCYCentFT0C -> GetYaxis() -> SetRangeUser(-2, 2);

    histQ1ZNACXXCentFT0C -> GetYaxis() -> SetRangeUser(-2, 2);
    histQ1ZNACYYCentFT0C -> GetYaxis() -> SetRangeUser(-2, 2);
    histQ1ZNACYXCentFT0C -> GetYaxis() -> SetRangeUser(-2, 2);
    histQ1ZNACXYCentFT0C -> GetYaxis() -> SetRangeUser(-2, 2);

    TProfile *profQ1ZNAXCentFT0C = (TProfile*) histQ1ZNAXCentFT0C -> ProfileX("profQ1ZNAXCentFT0C");
    TProfile *profQ1ZNAYCentFT0C = (TProfile*) histQ1ZNAYCentFT0C -> ProfileX("profQ1ZNAYCentFT0C");
    TProfile *profQ1ZNCXCentFT0C = (TProfile*) histQ1ZNCXCentFT0C -> ProfileX("profQ1ZNCXCentFT0C");
    TProfile *profQ1ZNCYCentFT0C = (TProfile*) histQ1ZNCYCentFT0C -> ProfileX("profQ1ZNCYCentFT0C");

    TProfile *profQ1ZNACXXCentFT0C = (TProfile*) histQ1ZNACXXCentFT0C -> ProfileX("profQ1ZNACXXCentFT0C");
    TProfile *profQ1ZNACYYCentFT0C = (TProfile*) histQ1ZNACYYCentFT0C -> ProfileX("profQ1ZNACYYCentFT0C");
    TProfile *profQ1ZNACYXCentFT0C = (TProfile*) histQ1ZNACYXCentFT0C -> ProfileX("profQ1ZNACYXCentFT0C");
    TProfile *profQ1ZNACXYCentFT0C = (TProfile*) histQ1ZNACXYCentFT0C -> ProfileX("profQ1ZNACXYCentFT0C");

    TProfile *profIntercalibZNACentFT0C = (TProfile*) histIntercalibZNACentFT0C -> ProfileX("profIntercalibZNACentFT0C");
    TProfile *profIntercalibZNCCentFT0C = (TProfile*) histIntercalibZNCCentFT0C -> ProfileX("profIntercalibZNCCentFT0C");

    SetProfile(profQ1ZNAXCentFT0C, 1, 1, 20, 0.8);
    SetProfile(profQ1ZNAYCentFT0C, 1, 1, 20, 0.8);
    SetProfile(profQ1ZNCXCentFT0C, 1, 1, 20, 0.8);
    SetProfile(profQ1ZNCYCentFT0C, 1, 1, 20, 0.8);

    SetProfile(profQ1ZNACXXCentFT0C, 1, 1, 20, 0.5);
    SetProfile(profQ1ZNACYYCentFT0C, 1, 1, 20, 0.5);
    SetProfile(profQ1ZNACYXCentFT0C, 1, 1, 20, 0.5);
    SetProfile(profQ1ZNACXYCentFT0C, 1, 1, 20, 0.5);

    SetProfile(profIntercalibZNACentFT0C, 1, 1, 20, 0.5);
    SetProfile(profIntercalibZNCCentFT0C, 1, 1, 20, 0.5);

    TLine *lineUnity = new TLine(0, 0, 90, 0);
    lineUnity -> SetLineColor(kGray+2);
    lineUnity -> SetLineStyle(kDashed);

    TCanvas *canvasQ1ZNCentFT0C = new TCanvas("canvasQ1ZNCentFT0C", "", 1200, 1200);
    canvasQ1ZNCentFT0C -> Divide(2, 2);

    canvasQ1ZNCentFT0C -> cd(1);
    histQ1ZNAXCentFT0C -> Draw("COLZ");
    lineUnity -> Draw("SAME");
    profQ1ZNAXCentFT0C -> Draw("EP SAME");

    canvasQ1ZNCentFT0C -> cd(2);
    histQ1ZNAYCentFT0C -> Draw("COLZ");
    lineUnity -> Draw("SAME");
    profQ1ZNAYCentFT0C -> Draw("EP SAME");

    canvasQ1ZNCentFT0C -> cd(3);
    histQ1ZNCXCentFT0C -> Draw("COLZ");
    lineUnity -> Draw("SAME");
    profQ1ZNCXCentFT0C -> Draw("EP SAME");

    canvasQ1ZNCentFT0C -> cd(4);
    histQ1ZNCYCentFT0C -> Draw("COLZ");
    lineUnity -> Draw("SAME");
    profQ1ZNCYCentFT0C -> Draw("EP SAME");


    TCanvas *canvasQ1ZNACCentFT0C = new TCanvas("canvasQ1ZNACCentFT0C", "", 1200, 1200);
    canvasQ1ZNACCentFT0C -> Divide(2, 2);

    canvasQ1ZNACCentFT0C -> cd(1);
    histQ1ZNACXXCentFT0C -> Draw("COLZ");
    lineUnity -> Draw("SAME");
    profQ1ZNACXXCentFT0C -> Draw("EP SAME");

    canvasQ1ZNACCentFT0C -> cd(2);
    histQ1ZNACYYCentFT0C -> Draw("COLZ");
    lineUnity -> Draw("SAME");
    profQ1ZNACYYCentFT0C -> Draw("EP SAME");

    canvasQ1ZNACCentFT0C -> cd(3);
    histQ1ZNACYXCentFT0C -> Draw("COLZ");
    lineUnity -> Draw("SAME");
    profQ1ZNACYXCentFT0C -> Draw("EP SAME");

    canvasQ1ZNACCentFT0C -> cd(4);
    histQ1ZNACXYCentFT0C -> Draw("COLZ");
    lineUnity -> Draw("SAME");
    profQ1ZNACXYCentFT0C -> Draw("EP SAME");


    double minCentrRange[] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
    double maxCentrRange[] = {10, 20, 30, 40, 50, 60, 70, 80, 90};

    TCanvas *canvasQ1ZNA = new TCanvas("canvasQ1ZNA", "", 1800, 1800);
    canvasQ1ZNA -> Divide(3, 3);

    for (int iCentr = 0;iCentr < 9;iCentr++) {
        double minCentrBin = histSparseZNA -> GetAxis(2) -> FindBin(minCentrRange[iCentr]);
        double maxCentrBin = histSparseZNA -> GetAxis(2) -> FindBin(maxCentrRange[iCentr] - 0.01);
        histSparseZNA -> GetAxis(2) -> SetRange(minCentrBin, maxCentrBin);
        TH2D *histProj = (TH2D*) histSparseZNA -> Projection(1, 0, "");
        histProj -> SetName(Form("Centr. FT0C %1.0f - %1.0f %%", minCentrRange[iCentr], maxCentrRange[iCentr]));
        histProj -> SetTitle(Form("Centr. FT0C %1.0f - %1.0f %%", minCentrRange[iCentr], maxCentrRange[iCentr]));
        canvasQ1ZNA -> cd(iCentr+1);
        histProj -> GetXaxis() -> SetRangeUser(-3., 3.);
        histProj -> GetYaxis() -> SetRangeUser(-3., 3.);
        histProj -> Draw("COLZ");
    }

    TCanvas *canvasQ1ZNC = new TCanvas("canvasQ1ZNC", "", 1800, 1800);
    canvasQ1ZNC -> Divide(3, 3);

    for (int iCentr = 0;iCentr < 9;iCentr++) {
        double minCentrBin = histSparseZNC -> GetAxis(2) -> FindBin(minCentrRange[iCentr]);
        double maxCentrBin = histSparseZNC -> GetAxis(2) -> FindBin(maxCentrRange[iCentr] - 0.01);
        histSparseZNC -> GetAxis(2) -> SetRange(minCentrBin, maxCentrBin);
        TH2D *histProj = (TH2D*) histSparseZNC -> Projection(1, 0, "");
        histProj -> SetName(Form("Centr. FT0C %1.0f - %1.0f %%", minCentrRange[iCentr], maxCentrRange[iCentr]));
        histProj -> SetTitle(Form("Centr. FT0C %1.0f - %1.0f %%", minCentrRange[iCentr], maxCentrRange[iCentr]));
        canvasQ1ZNC -> cd(iCentr+1);
        histProj -> GetXaxis() -> SetRangeUser(-3., 3.);
        histProj -> GetYaxis() -> SetRangeUser(-3., 3.);
        histProj -> Draw("COLZ");
    }

    TCanvas *canvasIntercalibZNACentFT0C = new TCanvas("canvasIntercalibZNACentFT0C", "", 800, 600);
    histIntercalibZNACentFT0C -> Draw("COLZ");
    lineUnity -> Draw("SAME");
    profIntercalibZNACentFT0C -> Draw("EP SAME");

    TCanvas *canvasIntercalibZNCCentFT0C = new TCanvas("canvasIntercalibZNCCentFT0C", "", 800, 600);
    histIntercalibZNCCentFT0C -> Draw("COLZ");
    lineUnity -> Draw("SAME");
    profIntercalibZNCCentFT0C -> Draw("EP SAME");

    canvasQ1ZNCentFT0C -> SaveAs("plots/Q1ZNCentFT0C.pdf");
    canvasQ1ZNACCentFT0C -> SaveAs("plots/Q1ZNACCentFT0C.pdf");
    canvasQ1ZNA -> SaveAs("plots/Q1ZNA.pdf");
    canvasQ1ZNC -> SaveAs("plots/Q1ZNC.pdf");
    canvasIntercalibZNACentFT0C -> SaveAs("plots/IntercalibZNACentFT0C.pdf");
    canvasIntercalibZNCCentFT0C -> SaveAs("plots/IntercalibZNCCentFT0C.pdf");
}
////////////////////////////////////////////////////////////////////////////////
void SetProfile(TProfile *prof, int color, int lineWidth, int markerStyle, double markerSize) {
    prof -> SetLineColor(color);
    prof -> SetLineWidth(lineWidth);
    prof -> SetMarkerColor(color);
    prof -> SetMarkerStyle(markerStyle);
    prof -> SetMarkerSize(markerSize);
}