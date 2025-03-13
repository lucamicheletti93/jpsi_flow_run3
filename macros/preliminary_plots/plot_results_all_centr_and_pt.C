void LoadStyle();
void SetLegend(TLegend *);

inline void EvalXaxis(int nVarBins, double minVarX[], double maxVarX[], double varCentr[], double varWidthStat[], double varWidthSyst[]) {
    for (int iVar = 0;iVar < nVarBins;iVar++) {
        varCentr[iVar] = (maxVarX[iVar] + minVarX[iVar]) / 2.;
        varWidthStat[iVar] = (maxVarX[iVar] - minVarX[iVar]) / 2.;
        varWidthSyst[iVar] = 0.15;
    }
}

void plot_results_all_centr_and_pt() {
    LoadStyle();

    TLatex *latexTitle = new TLatex();
    latexTitle -> SetTextSize(0.050);
    latexTitle -> SetNDC();
    latexTitle -> SetTextFont(42);

    TLine *lineZero = new TLine(-0.5, 0, 15.5, 0);
    lineZero -> SetLineStyle(kDashed);
    lineZero -> SetLineWidth(2);
    lineZero -> SetLineColor(kGray+1);

    // 10-30%
    const int nPtBins1030Run3 = 13;
    double ptMin1030Run3[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0};
    double ptMax1030Run3[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0};
    double ptCentr1030Run3[nPtBins1030Run3], ptWidthStat1030Run3[nPtBins1030Run3], ptWidthSyst1030Run3[nPtBins1030Run3];
    EvalXaxis(nPtBins1030Run3, ptMin1030Run3, ptMax1030Run3, ptCentr1030Run3, ptWidthStat1030Run3, ptWidthSyst1030Run3);

    double v2JpsiVsPt1030Run3[] = {-0.00756602, 0.0130263, 0.0235715, 0.047108, 0.0616868, 0.0853649, 0.074347, 0.0923566, 0.116014, 0.068199, 0.0728007, 0.0361536, 0.08162}; 
    double statV2JpsiVsPt1030Run3[] = {0.0133131, 0.00996566, 0.00768733, 0.00713583, 0.00743519, 0.00841299, 0.00748719, 0.0106387, 0.0148189, 0.0138707, 0.0249091, 0.0402871, 0.0512455};
    double systV2JpsiVsPt1030Run3[] = {0.00894376, 0.00836294, 0.00285496, 0.00460982, 0.00338182, 0.00408964, 0.00336151, 0.00498123, 0.00504415, 0.0141917, 0.0167561, 0.00926357, 0.032};

    TGraphErrors *graStatV2JpsiVsPt1030Run3 = new TGraphErrors(nPtBins1030Run3, ptCentr1030Run3, v2JpsiVsPt1030Run3, ptWidthStat1030Run3, statV2JpsiVsPt1030Run3);
    graStatV2JpsiVsPt1030Run3 -> SetLineColor(kRed+1);
    graStatV2JpsiVsPt1030Run3 -> SetLineWidth(2);
    graStatV2JpsiVsPt1030Run3 -> SetMarkerStyle(20);
    graStatV2JpsiVsPt1030Run3 -> SetMarkerColor(kRed+1);

    TGraphErrors *graSystV2JpsiVsPt1030Run3 = new TGraphErrors(nPtBins1030Run3, ptCentr1030Run3, v2JpsiVsPt1030Run3, ptWidthSyst1030Run3, systV2JpsiVsPt1030Run3);
    graSystV2JpsiVsPt1030Run3 -> SetLineColor(kRed+1);
    graSystV2JpsiVsPt1030Run3 -> SetLineWidth(2);
    graSystV2JpsiVsPt1030Run3 -> SetMarkerStyle(20);
    graSystV2JpsiVsPt1030Run3 -> SetMarkerColor(kRed+1);
    graSystV2JpsiVsPt1030Run3 -> SetFillStyle(0);

    // 30-50%
    const int nPtBins3050Run3 = 13;
    double ptMin3050Run3[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0};
    double ptMax3050Run3[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0};
    double ptCentr3050Run3[nPtBins3050Run3], ptWidthStat3050Run3[nPtBins3050Run3], ptWidthSyst3050Run3[nPtBins3050Run3];
    EvalXaxis(nPtBins3050Run3, ptMin3050Run3, ptMax3050Run3, ptCentr3050Run3, ptWidthStat3050Run3, ptWidthSyst3050Run3);

    double v2JpsiVsPt3050Run3[] = {-0.00690106, 0.0187907, 0.0342579, 0.0392987, 0.0806316, 0.0655535, 0.0935674, 0.0960825, 0.0863652, 0.06843, 0.0704275, 0.0748724, 0.0816124};
    double statV2JpsiVsPt3050Run3[] = {0.015088, 0.011354, 0.00969041, 0.00951278, 0.00989837, 0.0106802, 0.00926315, 0.0122852, 0.0152306, 0.0158748, 0.0266773, 0.0487097, 0.0584373};
    double systV2JpsiVsPt3050Run3[] = {0.0104899, 0.00335386, 0.00721273, 0.00270693, 0.00310409, 0.00783956, 0.00191871, 0.00281299, 0.0054744, 0.00802774, 0.0107482, 0.018323, 0.0316189};

    TGraphErrors *graStatV2JpsiVsPt3050Run3 = new TGraphErrors(nPtBins3050Run3, ptCentr3050Run3, v2JpsiVsPt3050Run3, ptWidthStat3050Run3, statV2JpsiVsPt3050Run3);
    graStatV2JpsiVsPt3050Run3 -> SetLineColor(kRed+1);
    graStatV2JpsiVsPt3050Run3 -> SetLineWidth(2);
    graStatV2JpsiVsPt3050Run3 -> SetMarkerStyle(20);
    graStatV2JpsiVsPt3050Run3 -> SetMarkerColor(kRed+1);

    TGraphErrors *graSystV2JpsiVsPt3050Run3 = new TGraphErrors(nPtBins3050Run3, ptCentr3050Run3, v2JpsiVsPt3050Run3, ptWidthSyst3050Run3, systV2JpsiVsPt3050Run3);
    graSystV2JpsiVsPt3050Run3 -> SetLineColor(kRed+1);
    graSystV2JpsiVsPt3050Run3 -> SetLineWidth(2);
    graSystV2JpsiVsPt3050Run3 -> SetMarkerStyle(20);
    graSystV2JpsiVsPt3050Run3 -> SetMarkerColor(kRed+1);
    graSystV2JpsiVsPt3050Run3 -> SetFillStyle(0);

    // 50-80%
    const int nPtBins5080Run3 = 8;
    double ptMin5080Run3[] = {0.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
    double ptMax5080Run3[] = {2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0};
    double ptCentr5080Run3[nPtBins5080Run3], ptWidthStat5080Run3[nPtBins5080Run3], ptWidthSyst5080Run3[nPtBins5080Run3];
    EvalXaxis(nPtBins5080Run3, ptMin5080Run3, ptMax5080Run3, ptCentr5080Run3, ptWidthStat5080Run3, ptWidthSyst5080Run3);

    double v2JpsiVsPt5080Run3[] = {0.0117977, 0.0578502, 0.0509676, 0.073072, 0.10159, 0.0593964, 0.079027, 0.0210381};
    double statV2JpsiVsPt5080Run3[] = {0.00877944, 0.012978, 0.0161448, 0.0199122, 0.0262047, 0.0260302, 0.0455554, 0.0587834};
    double systV2JpsiVsPt5080Run3[] = {0.00685833, 0.00766383, 0.00216662, 0.00835467, 0.00831863, 0.012291, 0.022382, 0.0299183};

    TGraphErrors *graStatV2JpsiVsPt5080Run3 = new TGraphErrors(nPtBins5080Run3, ptCentr5080Run3, v2JpsiVsPt5080Run3, ptWidthStat5080Run3, statV2JpsiVsPt5080Run3);
    graStatV2JpsiVsPt5080Run3 -> SetLineColor(kRed+1);
    graStatV2JpsiVsPt5080Run3 -> SetLineWidth(2);
    graStatV2JpsiVsPt5080Run3 -> SetMarkerStyle(20);
    graStatV2JpsiVsPt5080Run3 -> SetMarkerColor(kRed+1);

    TGraphErrors *graSystV2JpsiVsPt5080Run3 = new TGraphErrors(nPtBins5080Run3, ptCentr5080Run3, v2JpsiVsPt5080Run3, ptWidthSyst5080Run3, systV2JpsiVsPt5080Run3);
    graSystV2JpsiVsPt5080Run3 -> SetLineColor(kRed+1);
    graSystV2JpsiVsPt5080Run3 -> SetLineWidth(2);
    graSystV2JpsiVsPt5080Run3 -> SetMarkerStyle(20);
    graSystV2JpsiVsPt5080Run3 -> SetMarkerColor(kRed+1);
    graSystV2JpsiVsPt5080Run3 -> SetFillStyle(0);


    // J/psi v2 vs pT
    // 10-30%
    TCanvas *canvasV2JpsiVsPtCentr1030 = new TCanvas("canvasV2JpsiVsPtCentr1030", "", 800, 600);
    canvasV2JpsiVsPtCentr1030 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiVsPtCentr1030 = new TH2D("histGridV2JpsiVsPtCentr1030", "", 100, -0.5, 15.5, 100, -0.05, 0.3);
    histGridV2JpsiVsPtCentr1030 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    //histGridV2JpsiVsPtCentr1030 -> GetYaxis() -> SetTitle("v_{2} {EP, |#Delta#eta > 2.5|}");
    histGridV2JpsiVsPtCentr1030 -> GetYaxis() -> SetTitle("v_{2}");
    histGridV2JpsiVsPtCentr1030 -> Draw();
    lineZero -> Draw();
    graStatV2JpsiVsPt1030Run3 -> Draw("EP SAME");
    graSystV2JpsiVsPt1030Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsPtCentr1030 = new TLegend(0.65, 0.60, 0.85, 0.80, " ", "brNDC");
    SetLegend(legendV2JpsiVsPtCentr1030);
    legendV2JpsiVsPtCentr1030 -> SetTextSize(0.05);
    legendV2JpsiVsPtCentr1030 -> AddEntry(graStatV2JpsiVsPt1030Run3, "Stat. Uncert.", "PL");
    legendV2JpsiVsPtCentr1030 -> AddEntry(graSystV2JpsiVsPt1030Run3, "Syst. Uncert.", "F");
    legendV2JpsiVsPtCentr1030 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.18, 0.77, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4, 10#minus30\%");

    // 30-50%
    TCanvas *canvasV2JpsiVsPt3050 = new TCanvas("canvasV2JpsiVsPt3050", "", 800, 600);
    canvasV2JpsiVsPt3050 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiVsPt3050 = new TH2D("histGridV2JpsiVsPt3050", "", 100, -0.5, 15.5, 100, -0.05, 0.3);
    histGridV2JpsiVsPt3050 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    //histGridV2JpsiVsPt3050 -> GetYaxis() -> SetTitle("v_{2} {EP, |#Delta#eta > 2.5|}");
    histGridV2JpsiVsPt3050 -> GetYaxis() -> SetTitle("v_{2}");
    histGridV2JpsiVsPt3050 -> Draw();
    lineZero -> Draw();
    graStatV2JpsiVsPt3050Run3 -> Draw("EP SAME");
    graSystV2JpsiVsPt3050Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsPt3050 = new TLegend(0.65, 0.60, 0.85, 0.80, " ", "brNDC");
    SetLegend(legendV2JpsiVsPt3050);
    legendV2JpsiVsPt3050 -> SetTextSize(0.05);
    legendV2JpsiVsPt3050 -> AddEntry(graStatV2JpsiVsPt3050Run3, "Stat. Uncert.", "PL");
    legendV2JpsiVsPt3050 -> AddEntry(graSystV2JpsiVsPt3050Run3, "Syst. Uncert.", "F");
    legendV2JpsiVsPt3050 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.18, 0.77, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4, 30#minus50\%");

    // 50-80%
    TCanvas *canvasV2JpsiVsPtCentr5080 = new TCanvas("canvasV2JpsiVsPtCentr5080", "", 800, 600);
    canvasV2JpsiVsPtCentr5080 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiVsPtCentr5080 = new TH2D("histGridV2JpsiVsPtCentr5080", "", 100, -0.5, 15.5, 100, -0.05, 0.3);
    histGridV2JpsiVsPtCentr5080 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    //histGridV2JpsiVsPtCentr5080 -> GetYaxis() -> SetTitle("v_{2} {EP, |#Delta#eta > 2.5|}");
    histGridV2JpsiVsPtCentr5080 -> GetYaxis() -> SetTitle("v_{2}");
    histGridV2JpsiVsPtCentr5080 -> Draw();
    lineZero -> Draw();
    graStatV2JpsiVsPt5080Run3 -> Draw("EP SAME");
    graSystV2JpsiVsPt5080Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsPtCentr5080 = new TLegend(0.65, 0.60, 0.85, 0.80, " ", "brNDC");
    SetLegend(legendV2JpsiVsPtCentr5080);
    legendV2JpsiVsPtCentr5080 -> SetTextSize(0.05);
    legendV2JpsiVsPtCentr5080 -> AddEntry(graStatV2JpsiVsPt5080Run3, "Stat. Uncert.", "PL");
    legendV2JpsiVsPtCentr5080 -> AddEntry(graSystV2JpsiVsPt5080Run3, "Syst. Uncert.", "F");
    legendV2JpsiVsPtCentr5080 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.18, 0.77, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4, 50#minus80\%");


    /////////////////////////////
    // J/psi v2 vs pT vs Run 2 //
    /////////////////////////////
    const int nPtBinsRun2 = 10;
    double ptCentrRun2[] = {0.64, 1.49, 2.47, 3.46, 4.45, 5.45, 6.819, 8.835, 10.84, 14.25};
    double ptWidthStatRun2[] = {0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
    double ptWidthSystRun2[] = {0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30};

    double v2JpsiVsPtCentr1030Run2[] = {0.011, 0.043, 0.074, 0.088, 0.085, 0.103, 0.083, 0.100, 0.049, 0.022};
    double statV2JpsiVsPtCentr1030Run2[] = {0.0085, 0.0069, 0.0069, 0.0077, 0.009, 0.011, 0.011, 0.018, 0.028, 0.032};
    double systV2JpsiVsPtCentr1030Run2[] = {0.0038391, 0.0036633, 0.004898, 0.0035068, 0.0037855, 0.0029726, 0.0036802, 0.0075789, 0.0093488, 0.0091828};

    TGraphErrors *graStatV2JpsiVsPtCentr1030Run2 = new TGraphErrors(nPtBinsRun2, ptCentrRun2, v2JpsiVsPtCentr1030Run2, ptWidthStatRun2, statV2JpsiVsPtCentr1030Run2);
    graStatV2JpsiVsPtCentr1030Run2 -> SetLineColor(kOrange+1);
    graStatV2JpsiVsPtCentr1030Run2 -> SetLineWidth(2);
    graStatV2JpsiVsPtCentr1030Run2 -> SetMarkerStyle(20);
    graStatV2JpsiVsPtCentr1030Run2 -> SetMarkerColor(kOrange+1);

    TGraphErrors *graSystV2JpsiVsPtCentr1030Run2 = new TGraphErrors(nPtBinsRun2, ptCentrRun2, v2JpsiVsPtCentr1030Run2, ptWidthSystRun2, systV2JpsiVsPtCentr1030Run2);
    graSystV2JpsiVsPtCentr1030Run2 -> SetLineColor(kOrange+1);
    graSystV2JpsiVsPtCentr1030Run2 -> SetLineWidth(2);
    graSystV2JpsiVsPtCentr1030Run2 -> SetMarkerStyle(20);
    graSystV2JpsiVsPtCentr1030Run2 -> SetMarkerColor(kOrange+1);
    graSystV2JpsiVsPtCentr1030Run2 -> SetFillStyle(0);


    double v2JpsiVsPt3050Run2[] = {0.0008, 0.029, 0.067, 0.099, 0.098, 0.101, 0.098, 0.092, 0.055, 0.026};
    double statV2JpsiVsPt3050Run2[] = {0.011, 0.0091, 0.0089, 0.0095, 0.011, 0.013, 0.023, 0.022, 0.037, 0.039};
    double systV2JpsiVsPt3050Run2[] = {0.0032389, 0.0032904, 0.003359, 0.0043267, 0.0065719, 0.0066256, 0.0065651, 0.0067724, 0.0075293, 0.0093145};

    TGraphErrors *graStatV2JpsiVsPt3050Run2 = new TGraphErrors(nPtBinsRun2, ptCentrRun2, v2JpsiVsPt3050Run2, ptWidthStatRun2, statV2JpsiVsPt3050Run2);
    graStatV2JpsiVsPt3050Run2 -> SetLineColor(kOrange+1);
    graStatV2JpsiVsPt3050Run2 -> SetLineWidth(2);
    graStatV2JpsiVsPt3050Run2 -> SetMarkerStyle(20);
    graStatV2JpsiVsPt3050Run2 -> SetMarkerColor(kOrange+1);

    TGraphErrors *graSystV2JpsiVsPt3050Run2 = new TGraphErrors(nPtBinsRun2, ptCentrRun2, v2JpsiVsPt3050Run2, ptWidthSystRun2, systV2JpsiVsPt3050Run2);
    graSystV2JpsiVsPt3050Run2 -> SetLineColor(kOrange+1);
    graSystV2JpsiVsPt3050Run2 -> SetLineWidth(2);
    graSystV2JpsiVsPt3050Run2 -> SetMarkerStyle(20);
    graSystV2JpsiVsPt3050Run2 -> SetMarkerColor(kOrange+1);
    graSystV2JpsiVsPt3050Run2 -> SetFillStyle(0);

    TCanvas *canvasV2JpsiVsPtCentr1030Run2VsRun3 = new TCanvas("canvasV2JpsiVsPtCentr1030Run2VsRun3", "", 800, 600);
    canvasV2JpsiVsPtCentr1030Run2VsRun3 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiVsPtCentr1030Run2VsRun3 = new TH2D("histGridV2JpsiVsPtCentr1030Run2VsRun3", "", 100, -0.5, 15.5, 100, -0.05, 0.3);
    histGridV2JpsiVsPtCentr1030Run2VsRun3 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    //histGridV2JpsiVsPtCentr1030Run2VsRun3 -> GetYaxis() -> SetTitle("v_{2} {EP, |#Delta#eta > 2.5|}");
    histGridV2JpsiVsPtCentr1030Run2VsRun3 -> GetYaxis() -> SetTitle("v_{2}");
    histGridV2JpsiVsPtCentr1030Run2VsRun3 -> Draw();
    lineZero -> Draw();
    graStatV2JpsiVsPtCentr1030Run2 -> Draw("EP SAME");
    graSystV2JpsiVsPtCentr1030Run2 -> Draw("E2 SAME");
    graStatV2JpsiVsPt1030Run3 -> Draw("EP SAME");
    graSystV2JpsiVsPt1030Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsPtCentr1030Run2VsRun3 = new TLegend(0.20, 0.60, 0.45, 0.80, " ", "brNDC");
    SetLegend(legendV2JpsiVsPtCentr1030Run2VsRun3);
    legendV2JpsiVsPtCentr1030Run2VsRun3 -> SetTextSize(0.05);
    legendV2JpsiVsPtCentr1030Run2VsRun3 -> AddEntry(graSystV2JpsiVsPtCentr1030Run2, "#sqrt{#it{s}_{NN}} = 5.02 TeV", "PF");
    legendV2JpsiVsPtCentr1030Run2VsRun3 -> AddEntry(graSystV2JpsiVsPt1030Run3, "#sqrt{#it{s}_{NN}} = 5.36 TeV", "PF"); 
    legendV2JpsiVsPtCentr1030Run2VsRun3 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary");
    latexTitle -> DrawLatex(0.18, 0.80, "Pb#minusPb, J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4, 10#minus30\%");


    TCanvas *canvasV2JpsiVsPt3050Run2VsRun3 = new TCanvas("canvasV2JpsiVsPt3050Run2VsRun3", "", 800, 600);
    canvasV2JpsiVsPt3050Run2VsRun3 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiVsPt3050Run2VsRun3 = new TH2D("histGridV2JpsiVsPt3050Run2VsRun3", "", 100, -0.5, 15.5, 100, -0.05, 0.3);
    histGridV2JpsiVsPt3050Run2VsRun3 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    //histGridV2JpsiVsPt3050Run2VsRun3 -> GetYaxis() -> SetTitle("v_{2} {EP, |#Delta#eta > 2.5|}");
    histGridV2JpsiVsPt3050Run2VsRun3 -> GetYaxis() -> SetTitle("v_{2}");
    histGridV2JpsiVsPt3050Run2VsRun3 -> Draw();
    lineZero -> Draw();
    graStatV2JpsiVsPt3050Run2 -> Draw("EP SAME");
    graSystV2JpsiVsPt3050Run2 -> Draw("E2 SAME");
    graStatV2JpsiVsPt3050Run3 -> Draw("EP SAME");
    graSystV2JpsiVsPt3050Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsPt3050Run2VsRun3 = new TLegend(0.20, 0.60, 0.45, 0.80, " ", "brNDC");
    SetLegend(legendV2JpsiVsPt3050Run2VsRun3);
    legendV2JpsiVsPt3050Run2VsRun3 -> SetTextSize(0.05);
    legendV2JpsiVsPt3050Run2VsRun3 -> AddEntry(graSystV2JpsiVsPt3050Run2, "#sqrt{#it{s}_{NN}} = 5.02 TeV", "PF");
    legendV2JpsiVsPt3050Run2VsRun3 -> AddEntry(graSystV2JpsiVsPt3050Run3, "#sqrt{#it{s}_{NN}} = 5.36 TeV", "PF"); 
    legendV2JpsiVsPt3050Run2VsRun3 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary");
    latexTitle -> DrawLatex(0.18, 0.80, "Pb#minusPb, J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4, 30#minus50\%");

    ///////////////////////////////
    // J/psi v2 vs pT collection //
    ///////////////////////////////
    TGraphErrors *graCollStatV2JpsiVsPt1030Run3 = new TGraphErrors(nPtBins1030Run3, ptCentr1030Run3, v2JpsiVsPt1030Run3, ptWidthStat1030Run3, statV2JpsiVsPt1030Run3);
    graCollStatV2JpsiVsPt1030Run3 -> SetLineColor(kRed+1);
    graCollStatV2JpsiVsPt1030Run3 -> SetLineWidth(2);
    graCollStatV2JpsiVsPt1030Run3 -> SetMarkerStyle(20);
    graCollStatV2JpsiVsPt1030Run3 -> SetMarkerColor(kRed+1);

    TGraphErrors *graCollSystV2JpsiVsPt1030Run3 = new TGraphErrors(nPtBins1030Run3, ptCentr1030Run3, v2JpsiVsPt1030Run3, ptWidthSyst1030Run3, systV2JpsiVsPt1030Run3);
    graCollSystV2JpsiVsPt1030Run3 -> SetLineColor(kRed+1);
    graCollSystV2JpsiVsPt1030Run3 -> SetLineWidth(2);
    graCollSystV2JpsiVsPt1030Run3 -> SetMarkerStyle(20);
    graCollSystV2JpsiVsPt1030Run3 -> SetMarkerColor(kRed+1);
    graCollSystV2JpsiVsPt1030Run3 -> SetFillStyle(0);

    TGraphErrors *graCollStatV2JpsiVsPt3050Run3 = new TGraphErrors(nPtBins3050Run3, ptCentr3050Run3, v2JpsiVsPt3050Run3, ptWidthStat3050Run3, statV2JpsiVsPt3050Run3);
    graCollStatV2JpsiVsPt3050Run3 -> SetLineColor(kGreen+2);
    graCollStatV2JpsiVsPt3050Run3 -> SetLineWidth(2);
    graCollStatV2JpsiVsPt3050Run3 -> SetMarkerStyle(20);
    graCollStatV2JpsiVsPt3050Run3 -> SetMarkerColor(kGreen+2);

    TGraphErrors *graCollSystV2JpsiVsPt3050Run3 = new TGraphErrors(nPtBins3050Run3, ptCentr3050Run3, v2JpsiVsPt3050Run3, ptWidthSyst3050Run3, systV2JpsiVsPt3050Run3);
    graCollSystV2JpsiVsPt3050Run3 -> SetLineColor(kGreen+2);
    graCollSystV2JpsiVsPt3050Run3 -> SetLineWidth(2);
    graCollSystV2JpsiVsPt3050Run3 -> SetMarkerStyle(20);
    graCollSystV2JpsiVsPt3050Run3 -> SetMarkerColor(kGreen+2);
    graCollSystV2JpsiVsPt3050Run3 -> SetFillStyle(0);

    TGraphErrors *graCollStatV2JpsiVsPt5080Run3 = new TGraphErrors(nPtBins5080Run3, ptCentr5080Run3, v2JpsiVsPt5080Run3, ptWidthStat5080Run3, statV2JpsiVsPt5080Run3);
    graCollStatV2JpsiVsPt5080Run3 -> SetLineColor(kAzure+4);
    graCollStatV2JpsiVsPt5080Run3 -> SetLineWidth(2);
    graCollStatV2JpsiVsPt5080Run3 -> SetMarkerStyle(20);
    graCollStatV2JpsiVsPt5080Run3 -> SetMarkerColor(kAzure+4);

    TGraphErrors *graCollSystV2JpsiVsPt5080Run3 = new TGraphErrors(nPtBins5080Run3, ptCentr5080Run3, v2JpsiVsPt5080Run3, ptWidthSyst5080Run3, systV2JpsiVsPt5080Run3);
    graCollSystV2JpsiVsPt5080Run3 -> SetLineColor(kAzure+4);
    graCollSystV2JpsiVsPt5080Run3 -> SetLineWidth(2);
    graCollSystV2JpsiVsPt5080Run3 -> SetMarkerStyle(20);
    graCollSystV2JpsiVsPt5080Run3 -> SetMarkerColor(kAzure+4);
    graCollSystV2JpsiVsPt5080Run3 -> SetFillStyle(0);

    TCanvas *canvasV2JpsiVsPtColl = new TCanvas("canvasV2JpsiVsPtColl", "", 800, 600);
    canvasV2JpsiVsPtColl -> SetTicks(1, 1);

    TH2D *histGridV2JpsiVsPtColl = new TH2D("histGridV2JpsiVsPtColl", "", 100, -0.5, 15.5, 100, -0.05, 0.3);
    histGridV2JpsiVsPtColl -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    //histGridV2JpsiVsPtColl -> GetYaxis() -> SetTitle("v_{2} {EP, |#Delta#eta > 2.5|}");
    histGridV2JpsiVsPtColl -> GetYaxis() -> SetTitle("v_{2}");
    histGridV2JpsiVsPtColl -> Draw();
    lineZero -> Draw();
    graCollStatV2JpsiVsPt1030Run3 -> Draw("EP SAME");
    graCollSystV2JpsiVsPt1030Run3 -> Draw("E2 SAME");
    graCollStatV2JpsiVsPt3050Run3 -> Draw("EP SAME");
    graCollSystV2JpsiVsPt3050Run3 -> Draw("E2 SAME");
    graCollStatV2JpsiVsPt5080Run3 -> Draw("EP SAME");
    graCollSystV2JpsiVsPt5080Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsPtColl = new TLegend(0.20, 0.60, 0.45, 0.80, " ", "brNDC");
    SetLegend(legendV2JpsiVsPtColl);
    legendV2JpsiVsPtColl -> SetTextSize(0.05);
    legendV2JpsiVsPtColl -> AddEntry(graCollSystV2JpsiVsPt1030Run3, "10#minus30\%", "PF");
    legendV2JpsiVsPtColl -> AddEntry(graCollSystV2JpsiVsPt3050Run3, "30#minus50\%", "PF");
    legendV2JpsiVsPtColl -> AddEntry(graCollSystV2JpsiVsPt5080Run3, "50#minus80\%", "PF");
    legendV2JpsiVsPtColl -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.18, 0.80, "Pb#minusPb, J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4");

    canvasV2JpsiVsPtCentr1030 -> SaveAs("v2JpsiVsPtCentr1030.pdf");
    canvasV2JpsiVsPtCentr1030Run2VsRun3 -> SaveAs("v2JpsiVsPtCentr1030Run2VsRun3.pdf");
    canvasV2JpsiVsPt3050 -> SaveAs("v2JpsiVsPt3050.pdf");
    canvasV2JpsiVsPt3050Run2VsRun3 -> SaveAs("v2JpsiVsPt3050Run2VsRun3.pdf");
    canvasV2JpsiVsPtCentr5080 -> SaveAs("v2JpsiVsPtCentr5080.pdf");
    canvasV2JpsiVsPtColl -> SaveAs("v2JpsiVsPtCollection.pdf");
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