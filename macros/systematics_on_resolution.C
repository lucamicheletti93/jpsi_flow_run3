void LoadStyle();
void SetLegend(TLegend *);

inline void EvalXaxis(int nVarBins, double minVarX[], double maxVarX[], double varCentr[], double varWidthStat[], double varWidthSyst[], double systWidth) {
    for (int iVar = 0;iVar < nVarBins;iVar++) {
        varCentr[iVar] = (maxVarX[iVar] + minVarX[iVar]) / 2.;
        varWidthStat[iVar] = (maxVarX[iVar] - minVarX[iVar]) / 2.;
        varWidthSyst[iVar] = systWidth;
    }
}

inline void SetGraph(TGraphErrors *graph, Color_t mkrCol, int mkrSty, Color_t lnCol, int lnWidth, int fillSty) {
    graph -> SetMarkerColor(mkrCol);
    graph -> SetMarkerStyle(mkrSty);
    graph -> SetLineColor(lnCol);
    graph -> SetLineWidth(lnWidth);
    graph -> SetFillStyle(fillSty);
}

void systematics_on_resolution() {
    LoadStyle();

    TLatex *latexTitle = new TLatex();
    latexTitle -> SetTextSize(0.050);
    latexTitle -> SetNDC();
    latexTitle -> SetTextFont(42);

    TLine *lineZeroV2Pt = new TLine(-0.5, 0, 15.5, 0);
    lineZeroV2Pt -> SetLineStyle(kDashed);
    lineZeroV2Pt -> SetLineWidth(2);
    lineZeroV2Pt -> SetLineColor(kGray+1);

    TLine *lineZeroV2Centr = new TLine(-0.5, 0, 90.5, 0);
    lineZeroV2Centr -> SetLineStyle(kDashed);
    lineZeroV2Centr -> SetLineWidth(2);
    lineZeroV2Centr -> SetLineColor(kGray+1);

    //////////////////////////////////////////////////////
    const int nCentrBinsInt = 3;
    double centrMinInt[] = {10, 30, 50};
    double centrMaxInt[] = {30, 50, 80};
    double centrCentrInt[nCentrBinsInt], centrWidthStatInt[nCentrBinsInt], centrWidthSystint[nCentrBinsInt];
    EvalXaxis(nCentrBinsInt, centrMinInt, centrMaxInt, centrCentrInt, centrWidthStatInt, centrWidthSystint, 0.15);
    
    double resoIntInt[] = {1. / 1.16148, 1. / 1.15516, 1. / 1.43847};
    double resoIntPt05[] = {1. / 1.16148, 1. / 1.15515, 1. / 1.43842};
    double resoIntPt515[] = {1. / 1.15923, 1. / 1.15679, 1. / 1.46334};

    double centrCentrIntPt05[nCentrBinsInt], centrCentrIntPt515[nCentrBinsInt];
    for(int iCentr = 0;iCentr < nCentrBinsInt;iCentr++) {
        centrCentrIntPt05[iCentr] = centrCentrInt[iCentr] - 1.5;
        centrCentrIntPt515[iCentr] = centrCentrInt[iCentr] + 1.5;

        double invResoIntInt = 1. / resoIntInt[iCentr];
        double invResoIntIntPt515 = 1. / resoIntPt515[iCentr];
        std::cout << Form("%1.0f-%1.0f | max diff = %f", centrMinInt[iCentr], centrMaxInt[iCentr], TMath::Abs(invResoIntInt - invResoIntIntPt515) / invResoIntInt) << std::endl;
    }

    TGraphErrors *graResoIntInt = new TGraphErrors(nCentrBinsInt, centrCentrInt, resoIntInt, centrWidthStatInt, 0);
    SetGraph(graResoIntInt, kRed+1, 20, kBlack, 2, 0);

    TGraphErrors *graResoIntPt05 = new TGraphErrors(nCentrBinsInt, centrCentrIntPt05, resoIntPt05, 0, 0);
    SetGraph(graResoIntPt05, kAzure+4, 20, kAzure+4, 2, 0);

    TGraphErrors *graResoIntPt515 = new TGraphErrors(nCentrBinsInt, centrCentrIntPt515, resoIntPt515, 0, 0);
    SetGraph(graResoIntPt515, kGreen+2, 20, kGreen+2, 2, 0);

    const int nCentrBinsDiff = 8;
    double centrMinDiff[] = {0, 10, 20, 30, 40, 50, 60, 70};
    double centrMaxDiff[] = {10, 20, 30, 40, 50, 60, 70, 80};
    double centrCentrDiffInt[nCentrBinsDiff], centrWidthStatDiffInt[nCentrBinsDiff], centrWidthSystDiffInt[nCentrBinsDiff];
    EvalXaxis(nCentrBinsDiff, centrMinDiff, centrMaxDiff, centrCentrDiffInt, centrWidthStatDiffInt, centrWidthSystDiffInt, 0.15);

    double resoDiffInt[] = {1. / 1.56041, 1. / 1.17495, 1. / 1.1283, 1. / 1.14024, 1. / 1.1949, 1. / 1.31888, 1. / 1.5944, 1. / 2.25324};
    double resoDiffIntPt05[] = {1. / 1.56042, 1. / 1.17495, 1. / 1.1283, 1. / 1.14024, 1. / 1.19499, 1. / 1.31887, 1. / 1.59439, 1. / 2.25323};
    double resoDiffIntPt515[] = {1. / 1.5468, 1. / 1.17394, 1. / 1.12826, 1. / 1.14049, 1. / 1.19562, 1. / 1.3211, 1. / 1.59965, 1. / 2.25422};

    std::cout << "---------------------------------" << std::endl;
    double centrCentrDiffIntPt05[nCentrBinsDiff], centrCentrDiffIntPt515[nCentrBinsDiff];
    for(int iCentr = 0;iCentr < nCentrBinsDiff;iCentr++) {
        centrCentrDiffIntPt05[iCentr] = centrCentrDiffInt[iCentr] - 1.;
        centrCentrDiffIntPt515[iCentr] = centrCentrDiffInt[iCentr] + 1.;

        double invResoDiffInt = 1. / resoDiffInt[iCentr];
        double invResoDiffIntPt515 = 1. / resoDiffIntPt515[iCentr];
        std::cout << Form("%1.0f-%1.0f | max diff = %f", centrMinDiff[iCentr], centrMaxDiff[iCentr], TMath::Abs(invResoDiffInt - invResoDiffIntPt515) / invResoDiffInt) << std::endl;
    }



    TGraphErrors *graResoDiffInt = new TGraphErrors(nCentrBinsDiff, centrCentrDiffInt, resoDiffInt, centrWidthStatDiffInt, 0);
    SetGraph(graResoDiffInt, kRed+1, 33, kBlack, 1, 0);

    TGraphErrors *graResoDiffIntPt05 = new TGraphErrors(nCentrBinsDiff, centrCentrDiffIntPt05, resoDiffIntPt05, 0, 0);
    SetGraph(graResoDiffIntPt05, kAzure+4, 33, kAzure+4, 2, 0);

    TGraphErrors *graResoDiffIntPt515 = new TGraphErrors(nCentrBinsDiff, centrCentrDiffIntPt515, resoDiffIntPt515, 0, 0);
    SetGraph(graResoDiffIntPt515, kGreen+2, 33, kGreen+2, 2, 0);

    TH2D *histGridReso = new TH2D("histGridReso", "", 100, -0.5, 80, 100, 0, 1);
    histGridReso -> GetXaxis() -> SetTitle("Centrality (%)");
    histGridReso -> GetYaxis() -> SetTitle("R_{EP}");

    TCanvas *canvasReso = new TCanvas("canvasReso", "", 800, 600);

    TLine *lineCentr10 = new TLine(10, 0.3, 10, 1);
    lineCentr10 -> SetLineColor(kGray+1);
    lineCentr10 -> SetLineStyle(kDashed);
    
    TLine *lineCentr30 = new TLine(30, 0.3, 30, 1);
    lineCentr30 -> SetLineColor(kGray+1);
    lineCentr30 -> SetLineStyle(kDashed);
    
    TLine *lineCentr50 = new TLine(50, 0.3, 50, 1);
    lineCentr50 -> SetLineColor(kGray+1);
    lineCentr50 -> SetLineStyle(kDashed);
    
    histGridReso -> Draw();
    graResoIntInt -> Draw("EP SAME");
    graResoIntPt05 -> Draw("EP SAME");
    graResoIntPt515 -> Draw("EP SAME");
    graResoDiffInt -> Draw("EP SAME");
    graResoDiffIntPt05 -> Draw("EP SAME");
    graResoDiffIntPt515 -> Draw("EP SAME");

    lineCentr10 -> Draw("SAME");
    lineCentr30 -> Draw("SAME");
    lineCentr50 -> Draw("SAME");

    TLegend *legendReso = new TLegend(0.25, 0.20, 0.45, 0.40, " ", "brNDC");
    SetLegend(legendReso);
    legendReso -> SetTextSize(0.05);
    legendReso -> AddEntry(graResoIntInt, "Integrated", "PL");
    legendReso -> AddEntry(graResoIntPt05, "0 < #it{p}_{T} < 5 GeV/#it{c}", "P");
    legendReso -> AddEntry(graResoIntPt515, "5 < #it{p}_{T} < 15 GeV/#it{c}", "P");
    legendReso -> Draw();

    canvasReso -> SaveAs("systematic_on_resolution.pdf");
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