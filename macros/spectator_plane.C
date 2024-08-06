void SetProfile(TProfile *, int , int , int , double );
void CouputeCalibrations(TProfile *, TProfile *, TH1D *);

void spectator_plane() {
    //string pathToFiles = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/data/pass3/LHC23zzh_small/train_225020_spectator_plane/544122";
    //TFile *fIn = new TFile(Form("%s/AnalysisResults.root", pathToFiles.c_str()));
    string pathToFiles = "/Users/lucamicheletti/Downloads";
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

    TH2D *histEnergyCommonZNACentFT0C = (TH2D*)list2 -> FindObject("EnergyCommonZNA");
    TH2D *histEnergyCommonZNCCentFT0C = (TH2D*)list2 -> FindObject("EnergyCommonZNC");

    TH2D *histEnergyZNA1CentFT0C = (TH2D*)list2 -> FindObject("EnergyZNA1");
    TH2D *histEnergyZNA2CentFT0C = (TH2D*)list2 -> FindObject("EnergyZNA2");
    TH2D *histEnergyZNA3CentFT0C = (TH2D*)list2 -> FindObject("EnergyZNA3");
    TH2D *histEnergyZNA4CentFT0C = (TH2D*)list2 -> FindObject("EnergyZNA4");

    TH2D *histEnergyZNC1CentFT0C = (TH2D*)list2 -> FindObject("EnergyZNC1");
    TH2D *histEnergyZNC2CentFT0C = (TH2D*)list2 -> FindObject("EnergyZNC2");
    TH2D *histEnergyZNC3CentFT0C = (TH2D*)list2 -> FindObject("EnergyZNC3");
    TH2D *histEnergyZNC4CentFT0C = (TH2D*)list2 -> FindObject("EnergyZNC4");

    histQ1ZNAXCentFT0C -> GetYaxis() -> SetRangeUser(-2, 2);
    histQ1ZNAYCentFT0C -> GetYaxis() -> SetRangeUser(-2, 2);
    histQ1ZNCXCentFT0C -> GetYaxis() -> SetRangeUser(-2, 2);
    histQ1ZNCYCentFT0C -> GetYaxis() -> SetRangeUser(-2, 2);

    histQ1ZNACXXCentFT0C -> GetYaxis() -> SetRangeUser(-2, 2);
    histQ1ZNACYYCentFT0C -> GetYaxis() -> SetRangeUser(-2, 2);
    histQ1ZNACYXCentFT0C -> GetYaxis() -> SetRangeUser(-2, 2);
    histQ1ZNACXYCentFT0C -> GetYaxis() -> SetRangeUser(-2, 2);

    histEnergyCommonZNACentFT0C -> GetYaxis() -> SetRangeUser(0, 250);
    histEnergyCommonZNCCentFT0C -> GetYaxis() -> SetRangeUser(0, 250);

    histEnergyZNA1CentFT0C -> GetYaxis() -> SetRangeUser(0, 250);
    histEnergyZNA2CentFT0C -> GetYaxis() -> SetRangeUser(0, 250);
    histEnergyZNA3CentFT0C -> GetYaxis() -> SetRangeUser(0, 250);
    histEnergyZNA4CentFT0C -> GetYaxis() -> SetRangeUser(0, 250);

    histEnergyZNC1CentFT0C -> GetYaxis() -> SetRangeUser(0, 250);
    histEnergyZNC2CentFT0C -> GetYaxis() -> SetRangeUser(0, 250);
    histEnergyZNC3CentFT0C -> GetYaxis() -> SetRangeUser(0, 250);
    histEnergyZNC4CentFT0C -> GetYaxis() -> SetRangeUser(0, 250);

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

    TProfile *profEnergyCommonZNACentFT0C = (TProfile*) histEnergyCommonZNACentFT0C -> ProfileX("profEnergyCommonZNACentFT0C");
    TProfile *profEnergyCommonZNCCentFT0C = (TProfile*) histEnergyCommonZNCCentFT0C -> ProfileX("profEnergyCommonZNCCentFT0C");

    TProfile *profEnergyZNA1CentFT0C = (TProfile*) histEnergyZNA1CentFT0C -> ProfileX("profEnergyZNA1CentFT0C");
    TProfile *profEnergyZNA2CentFT0C = (TProfile*) histEnergyZNA2CentFT0C -> ProfileX("profEnergyZNA2CentFT0C");
    TProfile *profEnergyZNA3CentFT0C = (TProfile*) histEnergyZNA3CentFT0C -> ProfileX("profEnergyZNA3CentFT0C");
    TProfile *profEnergyZNA4CentFT0C = (TProfile*) histEnergyZNA4CentFT0C -> ProfileX("profEnergyZNA4CentFT0C");

    TProfile *profEnergyZNC1CentFT0C = (TProfile*) histEnergyZNC1CentFT0C -> ProfileX("profEnergyZNC1CentFT0C");
    TProfile *profEnergyZNC2CentFT0C = (TProfile*) histEnergyZNC2CentFT0C -> ProfileX("profEnergyZNC2CentFT0C");
    TProfile *profEnergyZNC3CentFT0C = (TProfile*) histEnergyZNC3CentFT0C -> ProfileX("profEnergyZNC3CentFT0C");
    TProfile *profEnergyZNC4CentFT0C = (TProfile*) histEnergyZNC4CentFT0C -> ProfileX("profEnergyZNC4CentFT0C");

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

    SetProfile(profEnergyCommonZNACentFT0C, 1, 1, 20, 0.5);
    SetProfile(profEnergyCommonZNCCentFT0C, 1, 1, 20, 0.5);

    SetProfile(profEnergyZNA1CentFT0C, 1, 1, 20, 0.5);
    SetProfile(profEnergyZNA2CentFT0C, 1, 1, 20, 0.5);
    SetProfile(profEnergyZNA3CentFT0C, 1, 1, 20, 0.5);
    SetProfile(profEnergyZNA4CentFT0C, 1, 1, 20, 0.5);

    SetProfile(profEnergyZNC1CentFT0C, 1, 1, 20, 0.5);
    SetProfile(profEnergyZNC2CentFT0C, 1, 1, 20, 0.5);
    SetProfile(profEnergyZNC3CentFT0C, 1, 1, 20, 0.5);
    SetProfile(profEnergyZNC4CentFT0C, 1, 1, 20, 0.5);

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

    TCanvas *canvasEnergyZDCCentFT0C = new TCanvas("canvasEnergyCommonZNACentFT0C", "", 3000, 1200);
    canvasEnergyZDCCentFT0C -> Divide(5, 2);

    canvasEnergyZDCCentFT0C -> cd(1);
    histEnergyCommonZNACentFT0C -> Draw("COLZ");
    profEnergyCommonZNACentFT0C -> Draw("EP SAME");
    canvasEnergyZDCCentFT0C -> cd(2);
    histEnergyZNA1CentFT0C -> Draw("COLZ");
    profEnergyZNA1CentFT0C -> Draw("EP SAME");
    canvasEnergyZDCCentFT0C -> cd(3);
    histEnergyZNA2CentFT0C -> Draw("COLZ");
    profEnergyZNA2CentFT0C -> Draw("EP SAME");
    canvasEnergyZDCCentFT0C -> cd(4);
    histEnergyZNA3CentFT0C -> Draw("COLZ");
    profEnergyZNA3CentFT0C -> Draw("EP SAME");
    canvasEnergyZDCCentFT0C -> cd(5);
    histEnergyZNA4CentFT0C -> Draw("COLZ");
    profEnergyZNA4CentFT0C -> Draw("EP SAME");

    canvasEnergyZDCCentFT0C -> cd(6);
    histEnergyCommonZNCCentFT0C -> Draw("COLZ");
    profEnergyCommonZNCCentFT0C -> Draw("EP SAME");
    canvasEnergyZDCCentFT0C -> cd(7);
    histEnergyZNC1CentFT0C -> Draw("COLZ");
    profEnergyZNC1CentFT0C -> Draw("EP SAME");
    canvasEnergyZDCCentFT0C -> cd(8);
    histEnergyZNC2CentFT0C -> Draw("COLZ");
    profEnergyZNC2CentFT0C -> Draw("EP SAME");
    canvasEnergyZDCCentFT0C -> cd(9);
    histEnergyZNC3CentFT0C -> Draw("COLZ");
    profEnergyZNC3CentFT0C -> Draw("EP SAME");
    canvasEnergyZDCCentFT0C -> cd(10);
    histEnergyZNC4CentFT0C -> Draw("COLZ");
    profEnergyZNC4CentFT0C -> Draw("EP SAME");


    // Compute ZDC energy calibrations
    TH1D *histEnergyCalibZNA1CentFT0C = new TH1D("histEnergyCalibZNA1CentFT0C", "", 90, 0, 90);
    CouputeCalibrations(profEnergyCommonZNACentFT0C, profEnergyZNA1CentFT0C, histEnergyCalibZNA1CentFT0C);

    TH1D *histEnergyCalibZNA2CentFT0C = new TH1D("histEnergyCalibZNA2CentFT0C", "", 90, 0, 90);
    CouputeCalibrations(profEnergyCommonZNACentFT0C, profEnergyZNA2CentFT0C, histEnergyCalibZNA2CentFT0C);

    TH1D *histEnergyCalibZNA3CentFT0C = new TH1D("histEnergyCalibZNA3CentFT0C", "", 90, 0, 90);
    CouputeCalibrations(profEnergyCommonZNACentFT0C, profEnergyZNA3CentFT0C, histEnergyCalibZNA3CentFT0C);

    TH1D *histEnergyCalibZNA4CentFT0C = new TH1D("histEnergyCalibZNA4CentFT0C", "", 90, 0, 90);
    CouputeCalibrations(profEnergyCommonZNACentFT0C, profEnergyZNA4CentFT0C, histEnergyCalibZNA4CentFT0C);

    TH1D *histEnergyCalibZNC1CentFT0C = new TH1D("histEnergyCalibZNC1CentFT0C", "", 90, 0, 90);
    CouputeCalibrations(profEnergyCommonZNCCentFT0C, profEnergyZNC1CentFT0C, histEnergyCalibZNC1CentFT0C);

    TH1D *histEnergyCalibZNC2CentFT0C = new TH1D("histEnergyCalibZNC2CentFT0C", "", 90, 0, 90);
    CouputeCalibrations(profEnergyCommonZNCCentFT0C, profEnergyZNC2CentFT0C, histEnergyCalibZNC2CentFT0C);

    TH1D *histEnergyCalibZNC3CentFT0C = new TH1D("histEnergyCalibZNC3CentFT0C", "", 90, 0, 90);
    CouputeCalibrations(profEnergyCommonZNCCentFT0C, profEnergyZNC3CentFT0C, histEnergyCalibZNC3CentFT0C);

    TH1D *histEnergyCalibZNC4CentFT0C = new TH1D("histEnergyCalibZNC4CentFT0C", "", 90, 0, 90);
    CouputeCalibrations(profEnergyCommonZNCCentFT0C, profEnergyZNC4CentFT0C, histEnergyCalibZNC4CentFT0C);

    /*
    TProfile *profEnergyCalibZNA1CentFT0C = (TProfile*) profEnergyCommonZNACentFT0C -> Clone("profEnergyCalibZNA1CentFT0C");
    profEnergyCalibZNA1CentFT0C -> Sumw2();
    profEnergyCalibZNA1CentFT0C -> Scale(1. / 4.);
    profEnergyCalibZNA1CentFT0C -> Divide(profEnergyZNA1CentFT0C);
    //profEnergyCalibZNA1CentFT0C -> GetYaxis() -> SetRangeUser(0, 2);

    TProfile *profEnergyCalibZNA2CentFT0C = (TProfile*) profEnergyCommonZNACentFT0C -> Clone("profEnergyCalibZNA2CentFT0C");
    profEnergyCalibZNA2CentFT0C -> Scale(1. / 4.);
    profEnergyCalibZNA2CentFT0C -> Divide(profEnergyZNA2CentFT0C);

    TProfile *profEnergyCalibZNA3CentFT0C = (TProfile*) profEnergyCommonZNACentFT0C -> Clone("profEnergyCalibZNA3CentFT0C");
    profEnergyCalibZNA3CentFT0C -> Scale(1. / 4.);
    profEnergyCalibZNA3CentFT0C -> Divide(profEnergyZNA3CentFT0C);

    TProfile *profEnergyCalibZNA4CentFT0C = (TProfile*) profEnergyCommonZNACentFT0C -> Clone("profEnergyCalibZNA4CentFT0C");
    profEnergyCalibZNA4CentFT0C -> Scale(1. / 4.);
    profEnergyCalibZNA4CentFT0C -> Divide(profEnergyZNA4CentFT0C);


    TProfile *profEnergyCalibZNC1CentFT0C = (TProfile*) profEnergyCommonZNCCentFT0C -> Clone("profEnergyCalibZNC1CentFT0C");
    profEnergyCalibZNC1CentFT0C -> Scale(1. / 4.);
    profEnergyCalibZNC1CentFT0C -> Divide(profEnergyZNC1CentFT0C);

    TProfile *profEnergyCalibZNC2CentFT0C = (TProfile*) profEnergyCommonZNCCentFT0C -> Clone("profEnergyCalibZNC2CentFT0C");
    profEnergyCalibZNC2CentFT0C -> Scale(1. / 4.);
    profEnergyCalibZNC2CentFT0C -> Divide(profEnergyZNC2CentFT0C);

    TProfile *profEnergyCalibZNC3CentFT0C = (TProfile*) profEnergyCommonZNCCentFT0C -> Clone("profEnergyCalibZNC3CentFT0C");
    profEnergyCalibZNC3CentFT0C -> Scale(1. / 4.);
    profEnergyCalibZNC3CentFT0C -> Divide(profEnergyZNC3CentFT0C);

    TProfile *profEnergyCalibZNC4CentFT0C = (TProfile*) profEnergyCommonZNCCentFT0C -> Clone("profEnergyCalibZNC4CentFT0C");
    profEnergyCalibZNC4CentFT0C -> Scale(1. / 4.);
    profEnergyCalibZNC4CentFT0C -> Divide(profEnergyZNC4CentFT0C);
    */


    TCanvas *canvasEnergyCalibZDCCentFT0C = new TCanvas("canvasEnergyCalibZDCCentFT0C", "", 2400, 1200);
    canvasEnergyCalibZDCCentFT0C -> Divide(4, 2);

    canvasEnergyCalibZDCCentFT0C -> cd(1);
    histEnergyCalibZNA1CentFT0C -> Draw("EP SAME");
    canvasEnergyCalibZDCCentFT0C -> cd(2);
    histEnergyCalibZNA2CentFT0C -> Draw("EP SAME");
    canvasEnergyCalibZDCCentFT0C -> cd(3);
    histEnergyCalibZNA3CentFT0C -> Draw("EP SAME");
    canvasEnergyCalibZDCCentFT0C -> cd(4);
    histEnergyCalibZNA4CentFT0C -> Draw("EP SAME");

    canvasEnergyCalibZDCCentFT0C -> cd(5);
    histEnergyCalibZNC1CentFT0C -> Draw("EP SAME");
    canvasEnergyCalibZDCCentFT0C -> cd(6);
    histEnergyCalibZNC2CentFT0C -> Draw("EP SAME");
    canvasEnergyCalibZDCCentFT0C -> cd(7);
    histEnergyCalibZNC3CentFT0C -> Draw("EP SAME");
    canvasEnergyCalibZDCCentFT0C -> cd(8);
    histEnergyCalibZNC4CentFT0C -> Draw("EP SAME");



    return;
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
////////////////////////////////////////////////////////////////////////////////
void CouputeCalibrations(TProfile *profEnergyCommon, TProfile *profEnergy, TH1D *histCalib) {

    histCalib -> GetYaxis() -> SetRangeUser(0, 2);
    histCalib -> SetLineColor(kBlack);
    histCalib -> SetLineWidth(1);
    histCalib -> SetMarkerColor(kBlack);
    histCalib -> SetMarkerStyle(20);
    histCalib -> SetMarkerSize(0.8);

    for (int iCentr = 0;iCentr < 90;iCentr++) {
        double energyCommon = profEnergyCommon -> GetBinContent(iCentr+1) / 4.;
        double errEnergyCommon = profEnergyCommon -> GetBinError(iCentr+1) / 4.;
        double energy = profEnergy -> GetBinContent(iCentr+1);
        double errEnergy = profEnergy -> GetBinError(iCentr+1);
        histCalib -> SetBinContent(iCentr+1, energyCommon / energy);
        histCalib -> SetBinError(iCentr+1, energyCommon / energy * TMath::Sqrt(TMath::Power(errEnergyCommon / energyCommon, 2.) + TMath::Power(errEnergy / energy, 2.)));
    }
}