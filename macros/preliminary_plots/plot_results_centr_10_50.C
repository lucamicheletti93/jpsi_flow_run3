void LoadStyle();
void SetLegend(TLegend *);

void plot_results_centr_10_50() {
    LoadStyle();

    const int nPtBins = 10;
    double ptMinRun3[] = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12};
    double ptMaxRun3[] = {1, 2, 3, 4, 5, 6, 8, 10, 12, 15};
    double ptCentrRun3[nPtBins];
    double ptWidthStatRun3[nPtBins];
    double ptWidthSystRun3[nPtBins];

    for (int iPt = 0;iPt < 10;iPt++) {
        ptCentrRun3[iPt] = (ptMaxRun3[iPt] + ptMinRun3[iPt]) / 2.;
        ptWidthStatRun3[iPt] = (ptMaxRun3[iPt] - ptMinRun3[iPt]) / 2.;
        ptWidthSystRun3[iPt] = 0.30;
    }

    double v2JpsiVsPtCentr1050Run3[] = {0.00272,0.03622, 0.08255, 0.08029, 0.10396, 0.13124, 0.06305, 0.08404, 0.03606, 0.09862};
    double statV2JpsiVsPtCentr1050Run3[] = {0.00619, 0.00436, 0.00474, 0.00630, 0.00871, 0.01167, 0.01176, 0.01972, 0.03491, 0.04140};
    double systV2JpsiVsPtCentr1050Run3[] = {0.00574, 0.00304, 0.00231, 0.00457, 0.00405, 0.00505, 0.00751, 0.01160, 0.00780, 0.02168};

    TGraphErrors *graStatV2JpsiVsPtCentr1050Run3 = new TGraphErrors(nPtBins, ptCentrRun3, v2JpsiVsPtCentr1050Run3, ptWidthStatRun3, statV2JpsiVsPtCentr1050Run3);
    graStatV2JpsiVsPtCentr1050Run3 -> SetLineColor(kRed+1);
    graStatV2JpsiVsPtCentr1050Run3 -> SetLineWidth(2);
    graStatV2JpsiVsPtCentr1050Run3 -> SetMarkerStyle(20);
    graStatV2JpsiVsPtCentr1050Run3 -> SetMarkerColor(kRed+1);

    TGraphErrors *graSystV2JpsiVsPtCentr1050Run3 = new TGraphErrors(nPtBins, ptCentrRun3, v2JpsiVsPtCentr1050Run3, ptWidthSystRun3, systV2JpsiVsPtCentr1050Run3);
    graSystV2JpsiVsPtCentr1050Run3 -> SetLineColor(kRed+1);
    graSystV2JpsiVsPtCentr1050Run3 -> SetLineWidth(2);
    graSystV2JpsiVsPtCentr1050Run3 -> SetMarkerStyle(20);
    graSystV2JpsiVsPtCentr1050Run3 -> SetMarkerColor(kRed+1);
    graSystV2JpsiVsPtCentr1050Run3 -> SetFillStyle(0);


    // J/psi v2 vs pT in 10-50%
    TCanvas *canvasV2JpsiVsPtCentr1050 = new TCanvas("canvasV2JpsiVsPtCentr1050", "", 800, 600);
    canvasV2JpsiVsPtCentr1050 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiVsPtCentr1050 = new TH2D("histGridV2JpsiVsPtCentr1050", "", 100, -0.5, 15.5, 100, -0.05, 0.4);
    histGridV2JpsiVsPtCentr1050 -> GetXaxis() -> SetTitle("#it{p}_{T} GeV/#it{c}");
    //histGridV2JpsiVsPtCentr1050 -> GetYaxis() -> SetTitle("v_{2} {EP, |#Delta#eta > 2.5|}");
    histGridV2JpsiVsPtCentr1050 -> GetYaxis() -> SetTitle("v_{2}");
    histGridV2JpsiVsPtCentr1050 -> Draw();

    TLine *lineZero = new TLine(-0.5, 0, 15.5, 0);
    lineZero -> SetLineStyle(kDashed);
    lineZero -> SetLineWidth(2);
    lineZero -> SetLineColor(kGray+1);
    lineZero -> Draw();

    graStatV2JpsiVsPtCentr1050Run3 -> Draw("EP SAME");
    graSystV2JpsiVsPtCentr1050Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsPtCentr1050 = new TLegend(0.65, 0.60, 0.85, 0.80, " ", "brNDC");
    SetLegend(legendV2JpsiVsPtCentr1050);
    legendV2JpsiVsPtCentr1050 -> SetTextSize(0.05);
    legendV2JpsiVsPtCentr1050 -> AddEntry(graStatV2JpsiVsPtCentr1050Run3, "Stat. Uncert.", "PL");
    legendV2JpsiVsPtCentr1050 -> AddEntry(graSystV2JpsiVsPtCentr1050Run3, "Syst. Uncert.", "F");
    legendV2JpsiVsPtCentr1050 -> Draw();

    TLatex *latexTitle = new TLatex();
    latexTitle -> SetTextSize(0.050);
    latexTitle -> SetNDC();
    latexTitle -> SetTextFont(42);
    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.18, 0.77, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4, 10#minus50\%");



    // J/psi v2 vs pT vs Run 2
    double ptCentrRun2[] = {0.64, 1.49, 2.47, 3.46, 4.45, 5.45, 6.819, 8.835, 10.84, 14.25};
    double ptWidthStatRun2[] = {0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
    double ptWidthSystRun2[] = {0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30};

    double v2JpsiVsPtCentr1030Run2[] = {0.011, 0.043, 0.074, 0.088, 0.085, 0.103, 0.083, 0.100, 0.049, 0.022};
    double statV2JpsiVsPtCentr1030Run2[] = {0.0085, 0.0069, 0.0069, 0.0077, 0.009, 0.011, 0.011, 0.018, 0.028, 0.032};
    double systV2JpsiVsPtCentr1030Run2[] = {0.0038391, 0.0036633, 0.004898, 0.0035068, 0.0037855, 0.0029726, 0.0036802, 0.0075789, 0.0093488, 0.0091828};

    TGraphErrors *graStatV2JpsiVsPtCentr1030Run2 = new TGraphErrors(nPtBins, ptCentrRun2, v2JpsiVsPtCentr1030Run2, ptWidthStatRun2, statV2JpsiVsPtCentr1030Run2);
    graStatV2JpsiVsPtCentr1030Run2 -> SetLineColor(kOrange+1);
    graStatV2JpsiVsPtCentr1030Run2 -> SetLineWidth(2);
    graStatV2JpsiVsPtCentr1030Run2 -> SetMarkerStyle(20);
    graStatV2JpsiVsPtCentr1030Run2 -> SetMarkerColor(kOrange+1);

    TGraphErrors *graSystV2JpsiVsPtCentr1030Run2 = new TGraphErrors(nPtBins, ptCentrRun2, v2JpsiVsPtCentr1030Run2, ptWidthSystRun2, systV2JpsiVsPtCentr1030Run2);
    graSystV2JpsiVsPtCentr1030Run2 -> SetLineColor(kOrange+1);
    graSystV2JpsiVsPtCentr1030Run2 -> SetLineWidth(2);
    graSystV2JpsiVsPtCentr1030Run2 -> SetMarkerStyle(20);
    graSystV2JpsiVsPtCentr1030Run2 -> SetMarkerColor(kOrange+1);
    graSystV2JpsiVsPtCentr1030Run2 -> SetFillStyle(0);


    double v2JpsiVsPtCentr3050Run2[] = {0.0008, 0.029, 0.067, 0.099, 0.098, 0.101, 0.098, 0.092, 0.055, 0.026};
    double statV2JpsiVsPtCentr3050Run2[] = {0.011, 0.0091, 0.0089, 0.0095, 0.011, 0.013, 0.023, 0.022, 0.037, 0.039};
    double systV2JpsiVsPtCentr3050Run2[] = {0.0032389, 0.0032904, 0.003359, 0.0043267, 0.0065719, 0.0066256, 0.0065651, 0.0067724, 0.0075293, 0.0093145};

    TGraphErrors *graStatV2JpsiVsPtCentr3050Run2 = new TGraphErrors(nPtBins, ptCentrRun2, v2JpsiVsPtCentr3050Run2, ptWidthStatRun2, statV2JpsiVsPtCentr3050Run2);
    graStatV2JpsiVsPtCentr3050Run2 -> SetLineColor(kAzure+1);
    graStatV2JpsiVsPtCentr3050Run2 -> SetLineWidth(2);
    graStatV2JpsiVsPtCentr3050Run2 -> SetMarkerStyle(20);
    graStatV2JpsiVsPtCentr3050Run2 -> SetMarkerColor(kAzure+1);

    TGraphErrors *graSystV2JpsiVsPtCentr3050Run2 = new TGraphErrors(nPtBins, ptCentrRun2, v2JpsiVsPtCentr3050Run2, ptWidthSystRun2, systV2JpsiVsPtCentr3050Run2);
    graSystV2JpsiVsPtCentr3050Run2 -> SetLineColor(kAzure+1);
    graSystV2JpsiVsPtCentr3050Run2 -> SetLineWidth(2);
    graSystV2JpsiVsPtCentr3050Run2 -> SetMarkerStyle(20);
    graSystV2JpsiVsPtCentr3050Run2 -> SetMarkerColor(kAzure+1);
    graSystV2JpsiVsPtCentr3050Run2 -> SetFillStyle(0);

    TCanvas *canvasV2JpsiVsPtRun2VsRun3 = new TCanvas("canvasV2JpsiVsPtRun2VsRun3", "", 800, 600);
    canvasV2JpsiVsPtRun2VsRun3 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiVsPtRun2VsRun3 = new TH2D("histGridV2JpsiVsPtRun2VsRun3", "", 100, -0.5, 15.5, 100, -0.05, 0.4);
    histGridV2JpsiVsPtRun2VsRun3 -> GetXaxis() -> SetTitle("#it{p}_{T} GeV/#it{c}");
    //histGridV2JpsiVsPtRun2VsRun3 -> GetYaxis() -> SetTitle("v_{2} {EP, |#Delta#eta > 2.5|}");
    histGridV2JpsiVsPtRun2VsRun3 -> GetYaxis() -> SetTitle("v_{2}");
    histGridV2JpsiVsPtRun2VsRun3 -> Draw();

    lineZero -> Draw();

    graStatV2JpsiVsPtCentr1030Run2 -> Draw("EP SAME");
    graSystV2JpsiVsPtCentr1030Run2 -> Draw("E2 SAME");
    graStatV2JpsiVsPtCentr3050Run2 -> Draw("EP SAME");
    graSystV2JpsiVsPtCentr3050Run2 -> Draw("E2 SAME");
    graStatV2JpsiVsPtCentr1050Run3 -> Draw("EP SAME");
    graSystV2JpsiVsPtCentr1050Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsPtRun2 = new TLegend(0.20, 0.60, 0.45, 0.80, " ", "brNDC");
    SetLegend(legendV2JpsiVsPtRun2);
    legendV2JpsiVsPtRun2 -> SetTextSize(0.05);
    legendV2JpsiVsPtRun2 -> AddEntry(graSystV2JpsiVsPtCentr1030Run2, "10#minus30\%", "PF");
    legendV2JpsiVsPtRun2 -> AddEntry(graSystV2JpsiVsPtCentr3050Run2, "30#minus50\%", "PF");
    legendV2JpsiVsPtRun2 -> Draw();

    TLegend *legendV2JpsiVsPtRun2VsRun3 = new TLegend(0.50, 0.66, 0.75, 0.80, " ", "brNDC");
    SetLegend(legendV2JpsiVsPtRun2VsRun3);
    legendV2JpsiVsPtRun2VsRun3 -> SetTextSize(0.05);
    legendV2JpsiVsPtRun2VsRun3 -> AddEntry(graSystV2JpsiVsPtCentr1050Run3, "10#minus50\%", "PF");
    legendV2JpsiVsPtRun2VsRun3 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb, J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4");
    latexTitle -> DrawLatex(0.20, 0.77, "#sqrt{#it{s}_{NN}} = 5.02 TeV");
    latexTitle -> DrawLatex(0.50, 0.77, "#sqrt{#it{s}_{NN}} = 5.36 TeV");


    



    canvasV2JpsiVsPtCentr1050 -> SaveAs("v2JpsiVsPtCentr1050.pdf");
    canvasV2JpsiVsPtRun2VsRun3 -> SaveAs("v2JpsiVsPtRun2VsRun3.pdf");
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