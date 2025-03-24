void LoadStyle();
void SetLegend(TLegend *);

inline void EvalXaxis(int nVarBins, double minVarX[], double maxVarX[], double varCentr[], double varWidthStat[], double varWidthSyst[], double systWidth) {
    for (int iVar = 0;iVar < nVarBins;iVar++) {
        varCentr[iVar] = (maxVarX[iVar] + minVarX[iVar]) / 2.;
        varWidthStat[iVar] = (maxVarX[iVar] - minVarX[iVar]) / 2.;
        varWidthSyst[iVar] = systWidth;
    }
}

inline void SetGraph(auto *graph, Color_t mkrCol, int mkrSty, Color_t lnCol, int lnWidth, int fillSty) {
    graph -> SetMarkerColor(mkrCol);
    graph -> SetMarkerStyle(mkrSty);
    graph -> SetLineColor(lnCol);
    graph -> SetLineWidth(lnWidth);
    graph -> SetFillStyle(fillSty);
}

void plot_results_all_centr_and_pt() {
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

    // 10-30%
    const int nPtBins1030Run3 = 13;
    double ptMin1030Run3[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0};
    double ptMax1030Run3[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0};
    double ptCentr1030Run3[nPtBins1030Run3], ptWidthStat1030Run3[nPtBins1030Run3], ptWidthSyst1030Run3[nPtBins1030Run3];
    EvalXaxis(nPtBins1030Run3, ptMin1030Run3, ptMax1030Run3, ptCentr1030Run3, ptWidthStat1030Run3, ptWidthSyst1030Run3, 0.15);

    double v2JpsiVsPt1030Run3[] = {-0.00758901, 0.0132946, 0.0235542, 0.0470634, 0.0616915, 0.0854803, 0.0744066, 0.0937042, 0.115946, 0.0671242, 0.0761203, 0.0370125, 0.0803783}; 
    double statV2JpsiVsPt1030Run3[] = {0.0133299, 0.00989667, 0.0076973, 0.00713118, 0.00744479, 0.00842379, 0.007468, 0.0109966, 0.0148175, 0.0136262, 0.0274283, 0.0445005, 0.0524056};
    double systV2JpsiVsPt1030Run3[] = {0.00896739, 0.00841674, 0.00287965, 0.00464068, 0.00344948, 0.00411742, 0.00333794, 0.00928915, 0.00494801, 0.0152623, 0.0156752, 0.00966239, 0.031102};

    // WARNING: add the systematic on the resolution
    const double systResoVsPt1030Run3 = 0.01;
    for (int iPt = 0;iPt < nPtBins1030Run3;iPt++) {
        double systReso = v2JpsiVsPt1030Run3[iPt] * systResoVsPt1030Run3;
        systV2JpsiVsPt1030Run3[iPt] = TMath::Sqrt(systV2JpsiVsPt1030Run3[iPt]*systV2JpsiVsPt1030Run3[iPt] + systReso*systReso);
    }

    TGraphErrors *graStatV2JpsiVsPt1030Run3 = new TGraphErrors(nPtBins1030Run3, ptCentr1030Run3, v2JpsiVsPt1030Run3, ptWidthStat1030Run3, statV2JpsiVsPt1030Run3);
    SetGraph(graStatV2JpsiVsPt1030Run3, kRed+1, 20, kRed+1, 2, 0);

    TGraphErrors *graSystV2JpsiVsPt1030Run3 = new TGraphErrors(nPtBins1030Run3, ptCentr1030Run3, v2JpsiVsPt1030Run3, ptWidthSyst1030Run3, systV2JpsiVsPt1030Run3);
    SetGraph(graSystV2JpsiVsPt1030Run3, kRed+1, 20, kRed+1, 2, 0);

    // 30-50%
    const int nPtBins3050Run3 = 13;
    double ptMin3050Run3[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0};
    double ptMax3050Run3[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0};
    double ptCentr3050Run3[nPtBins3050Run3], ptWidthStat3050Run3[nPtBins3050Run3], ptWidthSyst3050Run3[nPtBins3050Run3];
    EvalXaxis(nPtBins3050Run3, ptMin3050Run3, ptMax3050Run3, ptCentr3050Run3, ptWidthStat3050Run3, ptWidthSyst3050Run3, 0.15);

    /*double v2JpsiVsPt3050Run3[] = {-0.00690106, 0.0187907, 0.0342579, 0.0392987, 0.0806316, 0.0655535, 0.0935674, 0.0960825, 0.0863652, 0.06843, 0.0704275, 0.0748724, 0.0816124};
    double statV2JpsiVsPt3050Run3[] = {0.015088, 0.011354, 0.00969041, 0.00951278, 0.00989837, 0.0106802, 0.00926315, 0.0122852, 0.0152306, 0.0158748, 0.0266773, 0.0487097, 0.0584373};
    double systV2JpsiVsPt3050Run3[] = {0.0104899, 0.00335386, 0.00721273, 0.00270693, 0.00310409, 0.00783956, 0.00191871, 0.00281299, 0.0054744, 0.00802774, 0.0107482, 0.018323, 0.0316189};*/

    double v2JpsiVsPt3050Run3[] = {-0.00638724, 0.0187648, 0.0342108, 0.0392447, 0.0805207, 0.0654634, 0.0934388, 0.0959505, 0.0862465, 0.068336, 0.0703307, 0.074402, 0.0815052};
    double statV2JpsiVsPt3050Run3[] = {0.0150673, 0.0113384, 0.0096771, 0.00949971, 0.00988476, 0.0106655, 0.00925043, 0.0122683, 0.0152096, 0.015853, 0.0266407, 0.0486423, 0.058357};
    double systV2JpsiVsPt3050Run3[] = {0.0104801, 0.00334925, 0.00720282, 0.0027032, 0.00309982, 0.00782878, 0.00191604, 0.00280911, 0.00546688, 0.00801671, 0.0107334, 0.0185521, 0.0314803};

    // WARNING: add the systematic on the resolution
    const double systResoVsPt3050Run3 = 0.01;
    for (int iPt = 0;iPt < nPtBins3050Run3;iPt++) {
        double systReso = v2JpsiVsPt3050Run3[iPt] * systResoVsPt3050Run3;
        systV2JpsiVsPt3050Run3[iPt] = TMath::Sqrt(systV2JpsiVsPt3050Run3[iPt]*systV2JpsiVsPt3050Run3[iPt] + systReso*systReso);
    }

    TGraphErrors *graStatV2JpsiVsPt3050Run3 = new TGraphErrors(nPtBins3050Run3, ptCentr3050Run3, v2JpsiVsPt3050Run3, ptWidthStat3050Run3, statV2JpsiVsPt3050Run3);
    SetGraph(graStatV2JpsiVsPt3050Run3, kRed+1, 20, kRed+1, 2, 0);

    TGraphErrors *graSystV2JpsiVsPt3050Run3 = new TGraphErrors(nPtBins3050Run3, ptCentr3050Run3, v2JpsiVsPt3050Run3, ptWidthSyst3050Run3, systV2JpsiVsPt3050Run3);
    SetGraph(graSystV2JpsiVsPt3050Run3, kRed+1, 20, kRed+1, 2, 0);

    // 50-80%
    const int nPtBins5080Run3 = 8;
    double ptMin5080Run3[] = {0.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
    double ptMax5080Run3[] = {2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0};
    double ptCentr5080Run3[nPtBins5080Run3], ptWidthStat5080Run3[nPtBins5080Run3], ptWidthSyst5080Run3[nPtBins5080Run3];
    EvalXaxis(nPtBins5080Run3, ptMin5080Run3, ptMax5080Run3, ptCentr5080Run3, ptWidthStat5080Run3, ptWidthSyst5080Run3, 0.15);

    /*double v2JpsiVsPt5080Run3[] = {0.0126049, 0.0625654, 0.0552993, 0.0781961, 0.110691, 0.0636064, 0.0832959, 0.0249238};
    double statV2JpsiVsPt5080Run3[] = {0.00948747, 0.014073, 0.0174842, 0.0209127, 0.02838, 0.0281896, 0.0490464, 0.0634836};
    double systV2JpsiVsPt5080Run3[] = {0.00739106, 0.00833212, 0.00254726, 0.00659095, 0.00758594, 0.0112154, 0.0215367, 0.0327847};*/

    double v2JpsiVsPt5080Run3[] = {0.0103796, 0.05169, 0.0456227, 0.0646736, 0.0916401, 0.052272, 0.0691233, 0.0185192};
    double statV2JpsiVsPt5080Run3[] = {0.00787057, 0.0115652, 0.0143992, 0.0172359, 0.0233886, 0.0232925, 0.0406867, 0.0531321};
    double systV2JpsiVsPt5080Run3[] = {0.00614408, 0.00681589, 0.00196597, 0.00543031, 0.00625659, 0.00943061, 0.018147, 0.026985};

    // WARNING: add the systematic on the resolution
    const double systResoVsPt5080Run3 = 0.017;
    for (int iPt = 0;iPt < nPtBins5080Run3;iPt++) {
        double systReso = v2JpsiVsPt5080Run3[iPt] * systResoVsPt5080Run3;
        systV2JpsiVsPt5080Run3[iPt] = TMath::Sqrt(systV2JpsiVsPt5080Run3[iPt]*systV2JpsiVsPt5080Run3[iPt] + systReso*systReso);
    }

    TGraphErrors *graStatV2JpsiVsPt5080Run3 = new TGraphErrors(nPtBins5080Run3, ptCentr5080Run3, v2JpsiVsPt5080Run3, ptWidthStat5080Run3, statV2JpsiVsPt5080Run3);
    SetGraph(graStatV2JpsiVsPt5080Run3, kRed+1, 20, kRed+1, 2, 0);

    TGraphErrors *graSystV2JpsiVsPt5080Run3 = new TGraphErrors(nPtBins5080Run3, ptCentr5080Run3, v2JpsiVsPt5080Run3, ptWidthSyst5080Run3, systV2JpsiVsPt5080Run3);
    SetGraph(graSystV2JpsiVsPt5080Run3, kRed+1, 20, kRed+1, 2, 0);

    // J/psi v2 vs pT
    // 10-30%
    TCanvas *canvasV2JpsiVsPtCentr1030 = new TCanvas("canvasV2JpsiVsPtCentr1030", "", 800, 600);
    canvasV2JpsiVsPtCentr1030 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiVsPtCentr1030 = new TH2D("histGridV2JpsiVsPtCentr1030", "", 100, -0.5, 15.5, 100, -0.05, 0.3);
    histGridV2JpsiVsPtCentr1030 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    //histGridV2JpsiVsPtCentr1030 -> GetYaxis() -> SetTitle("v_{2} {EP, |#Delta#eta > 2.5|}");
    histGridV2JpsiVsPtCentr1030 -> GetYaxis() -> SetTitle("v_{2}");
    histGridV2JpsiVsPtCentr1030 -> Draw();
    lineZeroV2Pt -> Draw();
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
    lineZeroV2Pt -> Draw();
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
    lineZeroV2Pt -> Draw();
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
    SetGraph(graStatV2JpsiVsPtCentr1030Run2, kGray+2, 20, kGray+2, 2, 0);

    TGraphErrors *graSystV2JpsiVsPtCentr1030Run2 = new TGraphErrors(nPtBinsRun2, ptCentrRun2, v2JpsiVsPtCentr1030Run2, ptWidthSystRun2, systV2JpsiVsPtCentr1030Run2);
    SetGraph(graSystV2JpsiVsPtCentr1030Run2, kGray+2, 20, kGray+2, 2, 0);


    double v2JpsiVsPt3050Run2[] = {0.0008, 0.029, 0.067, 0.099, 0.098, 0.101, 0.098, 0.092, 0.055, 0.026};
    double statV2JpsiVsPt3050Run2[] = {0.011, 0.0091, 0.0089, 0.0095, 0.011, 0.013, 0.023, 0.022, 0.037, 0.039};
    double systV2JpsiVsPt3050Run2[] = {0.0032389, 0.0032904, 0.003359, 0.0043267, 0.0065719, 0.0066256, 0.0065651, 0.0067724, 0.0075293, 0.0093145};

    TGraphErrors *graStatV2JpsiVsPt3050Run2 = new TGraphErrors(nPtBinsRun2, ptCentrRun2, v2JpsiVsPt3050Run2, ptWidthStatRun2, statV2JpsiVsPt3050Run2);
    SetGraph(graStatV2JpsiVsPt3050Run2, kGray+2, 20, kGray+2, 2, 0);

    TGraphErrors *graSystV2JpsiVsPt3050Run2 = new TGraphErrors(nPtBinsRun2, ptCentrRun2, v2JpsiVsPt3050Run2, ptWidthSystRun2, systV2JpsiVsPt3050Run2);
    SetGraph(graSystV2JpsiVsPt3050Run2, kGray+2, 20, kGray+2, 2, 0);

    TCanvas *canvasV2JpsiVsPtCentr1030Run2VsRun3 = new TCanvas("canvasV2JpsiVsPtCentr1030Run2VsRun3", "", 800, 600);
    canvasV2JpsiVsPtCentr1030Run2VsRun3 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiVsPtCentr1030Run2VsRun3 = new TH2D("histGridV2JpsiVsPtCentr1030Run2VsRun3", "", 100, -0.5, 15.5, 100, -0.05, 0.3);
    histGridV2JpsiVsPtCentr1030Run2VsRun3 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    //histGridV2JpsiVsPtCentr1030Run2VsRun3 -> GetYaxis() -> SetTitle("v_{2} {EP, |#Delta#eta > 2.5|}");
    histGridV2JpsiVsPtCentr1030Run2VsRun3 -> GetYaxis() -> SetTitle("v_{2}");
    histGridV2JpsiVsPtCentr1030Run2VsRun3 -> Draw();
    lineZeroV2Pt -> Draw();
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
    lineZeroV2Pt -> Draw();
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
    SetGraph(graCollStatV2JpsiVsPt1030Run3, kRed+1, 20, kRed+1, 2, 0);

    TGraphErrors *graCollSystV2JpsiVsPt1030Run3 = new TGraphErrors(nPtBins1030Run3, ptCentr1030Run3, v2JpsiVsPt1030Run3, ptWidthSyst1030Run3, systV2JpsiVsPt1030Run3);
    SetGraph(graCollSystV2JpsiVsPt1030Run3, kRed+1, 20, kRed+1, 2, 0);

    TGraphErrors *graCollStatV2JpsiVsPt3050Run3 = new TGraphErrors(nPtBins3050Run3, ptCentr3050Run3, v2JpsiVsPt3050Run3, ptWidthStat3050Run3, statV2JpsiVsPt3050Run3);
    SetGraph(graCollStatV2JpsiVsPt3050Run3, kGreen+2, 20, kGreen+2, 2, 0);

    TGraphErrors *graCollSystV2JpsiVsPt3050Run3 = new TGraphErrors(nPtBins3050Run3, ptCentr3050Run3, v2JpsiVsPt3050Run3, ptWidthSyst3050Run3, systV2JpsiVsPt3050Run3);
    SetGraph(graCollSystV2JpsiVsPt3050Run3, kGreen+2, 20, kGreen+2, 2, 0);

    TGraphErrors *graCollStatV2JpsiVsPt5080Run3 = new TGraphErrors(nPtBins5080Run3, ptCentr5080Run3, v2JpsiVsPt5080Run3, ptWidthStat5080Run3, statV2JpsiVsPt5080Run3);
    SetGraph(graCollStatV2JpsiVsPt5080Run3, kAzure+4, 20, kAzure+4, 2, 0);

    TGraphErrors *graCollSystV2JpsiVsPt5080Run3 = new TGraphErrors(nPtBins5080Run3, ptCentr5080Run3, v2JpsiVsPt5080Run3, ptWidthSyst5080Run3, systV2JpsiVsPt5080Run3);
    SetGraph(graCollSystV2JpsiVsPt5080Run3, kAzure+4, 20, kAzure+4, 2, 0);

    TCanvas *canvasV2JpsiVsPtColl = new TCanvas("canvasV2JpsiVsPtColl", "", 800, 600);
    canvasV2JpsiVsPtColl -> SetTicks(1, 1);

    TH2D *histGridV2JpsiVsPtColl = new TH2D("histGridV2JpsiVsPtColl", "", 100, -0.5, 15.5, 100, -0.05, 0.3);
    histGridV2JpsiVsPtColl -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    //histGridV2JpsiVsPtColl -> GetYaxis() -> SetTitle("v_{2} {EP, |#Delta#eta > 2.5|}");
    histGridV2JpsiVsPtColl -> GetYaxis() -> SetTitle("v_{2}");
    histGridV2JpsiVsPtColl -> Draw();
    lineZeroV2Pt -> Draw();
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



    // Collection plot with the same binning of 50-80%
    //double v2JpsiVsPt1030Run3LargeBins[] = {0.0251153, 0.071785, 0.0744066, 0.0937042, 0.115946, 0.0671242, 0.0761203, 0.0508805};
    //double statV2JpsiVsPt1030Run3LargeBins[] = {0.00410048, 0.00551796, 0.007468, 0.0109966, 0.0148175, 0.0136262, 0.0274283, 0.0333535};
    //double systV2JpsiVsPt1030Run3LargeBins[] = {0.00245636, 0.00224522, 0.00333794, 0.00928915, 0.00494801, 0.0152623, 0.0156752, 0.0113642};

    double v2JpsiVsPt1030Run3LargeBins[] = {0.0252683, 0.0721892, 0.0744066, 0.0937042, 0.115946, 0.0671242, 0.0761203, 0.0508805};
    double statV2JpsiVsPt1030Run3LargeBins[] = {0.0041272, 0.00555273, 0.007468, 0.0109966, 0.0148175, 0.0136262, 0.0274283, 0.0333535};
    double systV2JpsiVsPt1030Run3LargeBins[] = {0.00251499, 0.00225653, 0.00333794, 0.00928915, 0.00494801, 0.0152623, 0.0156752, 0.0113642};

    for (int iPt = 0;iPt < nPtBins1030Run3;iPt++) {
        double systReso = v2JpsiVsPt1030Run3LargeBins[iPt] * systResoVsPt1030Run3;
        systV2JpsiVsPt1030Run3LargeBins[iPt] = TMath::Sqrt(systV2JpsiVsPt1030Run3LargeBins[iPt]*systV2JpsiVsPt1030Run3LargeBins[iPt] + systReso*systReso);
    }

    //double v2JpsiVsPt3050Run3LargeBins[] = {0.0266931, 0.0742065, 0.0935674, 0.0960825, 0.0863652, 0.06843, 0.0704275, 0.0798585};
    //double statV2JpsiVsPt3050Run3LargeBins[] = {0.00541945, 0.00723374, 0.00926315, 0.0122852, 0.0152306, 0.0158748, 0.0266773, 0.0369357};
    //double systV2JpsiVsPt3050Run3LargeBins[] = {0.00222263, 0.00379158, 0.00191871, 0.00281299, 0.0054744, 0.00802774, 0.0107482, 0.0204162};

    double v2JpsiVsPt3050Run3LargeBins[] = {0.0266564, 0.0741045, 0.0934388, 0.0959505, 0.0862465, 0.068336, 0.0703307, 0.0797487};
    double statV2JpsiVsPt3050Run3LargeBins[] = {0.005412, 0.0072238, 0.00925043, 0.0122683, 0.0152096, 0.015853, 0.0266407, 0.0368849};
    double systV2JpsiVsPt3050Run3LargeBins[] = {0.00221957, 0.00378637, 0.00191604, 0.00280911, 0.00546688, 0.00801671, 0.0107334, 0.0203881};

    for (int iPt = 0;iPt < nPtBins3050Run3;iPt++) {
        double systReso = v2JpsiVsPt3050Run3LargeBins[iPt] * systResoVsPt3050Run3;
        systV2JpsiVsPt3050Run3LargeBins[iPt] = TMath::Sqrt(systV2JpsiVsPt3050Run3LargeBins[iPt]*systV2JpsiVsPt3050Run3LargeBins[iPt] + systReso*systReso);
    }

    TGraphErrors *graCollStatV2JpsiVsPt1030Run3LargeBins = new TGraphErrors(nPtBins5080Run3, ptCentr5080Run3, v2JpsiVsPt1030Run3LargeBins, ptWidthStat5080Run3, statV2JpsiVsPt1030Run3LargeBins);
    SetGraph(graCollStatV2JpsiVsPt1030Run3LargeBins, kRed+1, 20, kRed+1, 2, 0);

    TGraphErrors *graCollSystV2JpsiVsPt1030Run3LargeBins = new TGraphErrors(nPtBins5080Run3, ptCentr5080Run3, v2JpsiVsPt1030Run3LargeBins, ptWidthSyst5080Run3, systV2JpsiVsPt1030Run3LargeBins);
    SetGraph(graCollSystV2JpsiVsPt1030Run3LargeBins, kRed+1, 20, kRed+1, 2, 0);

    TGraphErrors *graCollStatV2JpsiVsPt3050Run3LargeBins = new TGraphErrors(nPtBins5080Run3, ptCentr5080Run3, v2JpsiVsPt3050Run3LargeBins, ptWidthStat5080Run3, statV2JpsiVsPt3050Run3LargeBins);
    SetGraph(graCollStatV2JpsiVsPt3050Run3LargeBins, kGreen+2, 20, kGreen+2, 2, 0);

    TGraphErrors *graCollSystV2JpsiVsPt3050Run3LargeBins = new TGraphErrors(nPtBins5080Run3, ptCentr5080Run3, v2JpsiVsPt3050Run3LargeBins, ptWidthSyst5080Run3, systV2JpsiVsPt3050Run3LargeBins);
    SetGraph(graCollSystV2JpsiVsPt3050Run3LargeBins, kGreen+2, 20, kGreen+2, 2, 0);

    TCanvas *canvasV2JpsiVsPtCollLargeBins = new TCanvas("canvasV2JpsiVsPtCollLargeBins", "", 800, 600);
    canvasV2JpsiVsPtCollLargeBins -> SetTicks(1, 1);
    histGridV2JpsiVsPtColl -> Draw();
    lineZeroV2Pt -> Draw();
    graCollStatV2JpsiVsPt1030Run3LargeBins -> Draw("EP SAME");
    graCollSystV2JpsiVsPt1030Run3LargeBins -> Draw("E2 SAME");
    graCollStatV2JpsiVsPt3050Run3LargeBins -> Draw("EP SAME");
    graCollSystV2JpsiVsPt3050Run3LargeBins -> Draw("E2 SAME");
    graCollStatV2JpsiVsPt5080Run3 -> Draw("EP SAME");
    graCollSystV2JpsiVsPt5080Run3 -> Draw("E2 SAME");
    legendV2JpsiVsPtColl -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.18, 0.80, "Pb#minusPb, J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4");

    ////////////////////////////////////////
    ////////////////////////////////////////
    // pT < 5 GeV/c
    const int nCentrBinsPt05Run3 = 8;
    double centrMinPt05Run3[] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0};
    double centrMaxPt05Run3[] = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0};
    double centrPt05Run3[nCentrBinsPt05Run3], centrWidthStatPt05Run3[nCentrBinsPt05Run3], centrWidthSystPt05Run3[nCentrBinsPt05Run3];
    EvalXaxis(nCentrBinsPt05Run3, centrMinPt05Run3, centrMaxPt05Run3, centrPt05Run3, centrWidthStatPt05Run3, centrWidthSystPt05Run3, 1.5);

    //double v2JpsiVsCentrPt05Run3[] = {0.0159439, 0.0454924, 0.0509332, 0.049155, 0.0623088, 0.0379151, 0.0210997, 0.0199327}; 
    //double statV2JpsiVsCentrPt05Run3[] = {0.0044295, 0.00389648, 0.00427198, 0.00487373, 0.00582226, 0.00743649, 0.0108959, 0.0197447};
    //double systV2JpsiVsCentrPt05Run3[] = {0.00652571, 0.00186915, 0.000983468, 0.000880583, 0.00215478, 0.0045007, 0.00313358, 0.00399369};

    double v2JpsiVsCentrPt05Run3[] = {0.0171446, 0.0459331, 0.0511079, 0.0491207, 0.0621999, 0.0375835, 0.0208883, 0.0199414}; 
    double statV2JpsiVsCentrPt05Run3[] = {0.00449126, 0.00394339, 0.00430864, 0.00487536, 0.00580195, 0.00736735, 0.0105835, 0.019498};  
    double systV2JpsiVsCentrPt05Run3[] = {0.00679839, 0.00168885, 0.000889025, 0.000873776, 0.00222375, 0.00447077, 0.0030507, 0.00387123};

    const double systResoVsCentrPt05Run3[nCentrBinsPt05Run3] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
    for (int iCentr = 0;iCentr < nCentrBinsPt05Run3;iCentr++) {
        double systReso = v2JpsiVsCentrPt05Run3[iCentr] * systResoVsCentrPt05Run3[iCentr];
        systV2JpsiVsCentrPt05Run3[iCentr] = TMath::Sqrt(systV2JpsiVsCentrPt05Run3[iCentr]*systV2JpsiVsCentrPt05Run3[iCentr] + systReso*systReso);
    }

    TGraphErrors *graStatV2JpsiVsCentrPt05Run3 = new TGraphErrors(nCentrBinsPt05Run3, centrPt05Run3, v2JpsiVsCentrPt05Run3, centrWidthStatPt05Run3, statV2JpsiVsCentrPt05Run3);
    SetGraph(graStatV2JpsiVsCentrPt05Run3, kRed+1, 20, kRed+1, 2, 0);

    TGraphErrors *graSystV2JpsiVsCentrPt05Run3 = new TGraphErrors(nCentrBinsPt05Run3, centrPt05Run3, v2JpsiVsCentrPt05Run3, centrWidthSystPt05Run3, systV2JpsiVsCentrPt05Run3);
    SetGraph(graSystV2JpsiVsCentrPt05Run3, kRed+1, 20, kRed+1, 2, 0);

    const int nCentrBinsPt05Run2 = 7;
    double centrCentrPt05Run2[] = {5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 75.0};
    double centrWidthStatPt05Run2[] = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 15.0};
    double centrWidthSystPt05Run2[] = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};

    double v2JpsiVsCentrPt05Run2[] = {0.036, 0.05, 0.056, 0.055, 0.031, 0.065, 0.045};
    double statV2JpsiVsCentrPt05Run2[] = {0.0093, 0.0081, 0.0084, 0.0099, 0.011, 0.014, 0.019};
    double systV2JpsiVsCentrPt05Run2[] = {0.0031078, 0.0044159, 0.0057615, 0.0053348, 0.0065608, 0.0065192, 0.0082189};

    TGraphErrors *graStatV2JpsiVsCentrPt05Run2 = new TGraphErrors(nCentrBinsPt05Run2, centrCentrPt05Run2, v2JpsiVsCentrPt05Run2, centrWidthStatPt05Run2, statV2JpsiVsCentrPt05Run2);
    SetGraph(graStatV2JpsiVsCentrPt05Run2, kGray+2, 20, kGray+2, 2, 0);

    TGraphErrors *graSystV2JpsiVsCentrPt05Run2 = new TGraphErrors(nCentrBinsPt05Run2, centrCentrPt05Run2, v2JpsiVsCentrPt05Run2, centrWidthSystPt05Run2, systV2JpsiVsCentrPt05Run2);
    SetGraph(graSystV2JpsiVsCentrPt05Run2, kGray+2, 20, kGray+2, 2, 0);

    TCanvas *canvasV2JpsiVsCentrPt05Run2VsRun3 = new TCanvas("canvasV2JpsiVsCentrPt05Run2VsRun3", "", 800, 600);
    canvasV2JpsiVsCentrPt05Run2VsRun3 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiVsCentrPt05Run2VsRun3 = new TH2D("histGridV2JpsiVsCentrPt05Run2VsRun3", "", 100, -0.5, 90.5, 100, -0.05, 0.3);
    histGridV2JpsiVsCentrPt05Run2VsRun3 -> GetXaxis() -> SetTitle("Centrality (%)");
    //histGridV2JpsiVsCentrPt05Run2VsRun3 -> GetYaxis() -> SetTitle("v_{2} {EP, |#Delta#eta > 2.5|}");
    histGridV2JpsiVsCentrPt05Run2VsRun3 -> GetYaxis() -> SetTitle("v_{2}");
    histGridV2JpsiVsCentrPt05Run2VsRun3 -> Draw();
    lineZeroV2Centr -> Draw();
    graStatV2JpsiVsCentrPt05Run2 -> Draw("EP SAME");
    graSystV2JpsiVsCentrPt05Run2 -> Draw("E2 SAME");
    graStatV2JpsiVsCentrPt05Run3 -> Draw("EP SAME");
    graSystV2JpsiVsCentrPt05Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsCentrPt05Run2VsRun3 = new TLegend(0.20, 0.60, 0.45, 0.80, " ", "brNDC");
    SetLegend(legendV2JpsiVsCentrPt05Run2VsRun3);
    legendV2JpsiVsCentrPt05Run2VsRun3 -> SetTextSize(0.05);
    legendV2JpsiVsCentrPt05Run2VsRun3 -> AddEntry(graSystV2JpsiVsCentrPt05Run2, "#sqrt{#it{s}_{NN}} = 5.02 TeV", "PF");
    legendV2JpsiVsCentrPt05Run2VsRun3 -> AddEntry(graSystV2JpsiVsCentrPt05Run3, "#sqrt{#it{s}_{NN}} = 5.36 TeV", "PF"); 
    legendV2JpsiVsCentrPt05Run2VsRun3 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary");
    latexTitle -> DrawLatex(0.18, 0.80, "Pb#minusPb, J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4, 0 < #it{p}_{T} < 5 GeV/#it{c}");


    // 5 < pT < 15 GeV/c
    const int nCentrBinsPt515Run3 = 8;
    double centrMinPt515Run3[] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0};
    double centrMaxPt515Run3[] = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0};
    double centrPt515Run3[nCentrBinsPt515Run3], centrWidthStatPt515Run3[nCentrBinsPt515Run3], centrWidthSystPt515Run3[nCentrBinsPt515Run3];
    EvalXaxis(nCentrBinsPt515Run3, centrMinPt515Run3, centrMaxPt515Run3, centrPt515Run3, centrWidthStatPt515Run3, centrWidthSystPt515Run3, 1.5);

    //double v2JpsiVsCentrPt515Run3[] = {0.0219178, 0.103274, 0.072014, 0.0797112, 0.0761756, 0.0956604, 0.0234934, 0.0455937}; 
    //double statV2JpsiVsCentrPt515Run3[] = {0.0178678, 0.0133656, 0.0127442, 0.0128951, 0.0149092, 0.0187528, 0.0286451, 0.0512898};
    //double systV2JpsiVsCentrPt515Run3[] = {0.01357, 0.0105526, 0.00794372, 0.00636919, 0.00863161, 0.00227993, 0.0210715, 0.009457};

    double v2JpsiVsCentrPt515Run3[] = {0.0233286, 0.103797, 0.0721063, 0.0795816, 0.0758477, 0.0948522, 0.0231561, 0.0448654}; 
    double statV2JpsiVsCentrPt515Run3[] = {0.0190183, 0.0134334, 0.0128311, 0.0128743, 0.0148452, 0.0185945, 0.0282124, 0.049162}; 
    double systV2JpsiVsCentrPt515Run3[] = {0.0144432, 0.0106041, 0.00787289, 0.00635375, 0.00859144, 0.00226504, 0.0207722, 0.00931585}; 

    const double systResoVsCentrPt515Run3[nCentrBinsPt515Run3] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
    for (int iCentr = 0;iCentr < nCentrBinsPt515Run3;iCentr++) {
        double systReso = v2JpsiVsCentrPt515Run3[iCentr] * systResoVsCentrPt515Run3[iCentr];
        systV2JpsiVsCentrPt515Run3[iCentr] = TMath::Sqrt(systV2JpsiVsCentrPt515Run3[iCentr]*systV2JpsiVsCentrPt515Run3[iCentr] + systReso*systReso);
    }

    TGraphErrors *graStatV2JpsiVsCentrPt515Run3 = new TGraphErrors(nCentrBinsPt515Run3, centrPt515Run3, v2JpsiVsCentrPt515Run3, centrWidthStatPt515Run3, statV2JpsiVsCentrPt515Run3);
    SetGraph(graStatV2JpsiVsCentrPt515Run3, kRed+1, 20, kRed+1, 2, 0);

    TGraphErrors *graSystV2JpsiVsCentrPt515Run3 = new TGraphErrors(nCentrBinsPt515Run3, centrPt515Run3, v2JpsiVsCentrPt515Run3, centrWidthSystPt515Run3, systV2JpsiVsCentrPt515Run3);
    SetGraph(graSystV2JpsiVsCentrPt515Run3, kRed+1, 20, kRed+1, 2, 0);

    const int nCentrBinsPt515Run2 = 7;
    double centrCentrPt515Run2[] = {5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 75.0};
    double centrWidthStatPt515Run2[] = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 15.0};
    double centrWidthSystPt515Run2[] = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};

    double v2JpsiVsCentrPt515Run2[] = {0.04, 0.064, 0.105, 0.096, 0.092, 0.106, 0.1};
    double statV2JpsiVsCentrPt515Run2[] = {0.02, 0.015, 0.015, 0.016, 0.019, 0.025, 0.036};
    double systV2JpsiVsCentrPt515Run2[] = {0.0044385, 0.0044843, 0.0065419, 0.0064786, 0.0083954, 0.0085108, 0.011057};

    TGraphErrors *graStatV2JpsiVsCentrPt515Run2 = new TGraphErrors(nCentrBinsPt515Run2, centrCentrPt515Run2, v2JpsiVsCentrPt515Run2, centrWidthStatPt515Run2, statV2JpsiVsCentrPt515Run2);
    SetGraph(graStatV2JpsiVsCentrPt515Run2, kGray+2, 20, kGray+2, 2, 0);

    TGraphErrors *graSystV2JpsiVsCentrPt515Run2 = new TGraphErrors(nCentrBinsPt515Run2, centrCentrPt515Run2, v2JpsiVsCentrPt515Run2, centrWidthSystPt515Run2, systV2JpsiVsCentrPt515Run2);
    SetGraph(graSystV2JpsiVsCentrPt515Run2, kGray+2, 20, kGray+2, 2, 0);

    TCanvas *canvasV2JpsiVsCentrPt515Run2VsRun3 = new TCanvas("canvasV2JpsiVsCentrPt515Run2VsRun3", "", 800, 600);
    canvasV2JpsiVsCentrPt515Run2VsRun3 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiVsCentrPt515Run2VsRun3 = new TH2D("histGridV2JpsiVsCentrPt515Run2VsRun3", "", 100, -0.5, 90.5, 100, -0.05, 0.3);
    histGridV2JpsiVsCentrPt515Run2VsRun3 -> GetXaxis() -> SetTitle("Centrality (%)");
    //histGridV2JpsiVsCentrPt515Run2VsRun3 -> GetYaxis() -> SetTitle("v_{2} {EP, |#Delta#eta > 2.5|}");
    histGridV2JpsiVsCentrPt515Run2VsRun3 -> GetYaxis() -> SetTitle("v_{2}");
    histGridV2JpsiVsCentrPt515Run2VsRun3 -> Draw();
    lineZeroV2Centr -> Draw();
    graStatV2JpsiVsCentrPt515Run2 -> Draw("EP SAME");
    graSystV2JpsiVsCentrPt515Run2 -> Draw("E2 SAME");
    graStatV2JpsiVsCentrPt515Run3 -> Draw("EP SAME");
    graSystV2JpsiVsCentrPt515Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsCentrPt515Run2VsRun3 = new TLegend(0.20, 0.60, 0.45, 0.80, " ", "brNDC");
    SetLegend(legendV2JpsiVsCentrPt515Run2VsRun3);
    legendV2JpsiVsCentrPt515Run2VsRun3 -> SetTextSize(0.05);
    legendV2JpsiVsCentrPt515Run2VsRun3 -> AddEntry(graSystV2JpsiVsCentrPt515Run2, "#sqrt{#it{s}_{NN}} = 5.02 TeV, 5 < #it{p}_{T} < 20 GeV/#it{c}", "PF");
    legendV2JpsiVsCentrPt515Run2VsRun3 -> AddEntry(graSystV2JpsiVsCentrPt515Run3, "#sqrt{#it{s}_{NN}} = 5.36 TeV, 5 < #it{p}_{T} < 15 GeV/#it{c}", "PF"); 
    legendV2JpsiVsCentrPt515Run2VsRun3 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary");
    latexTitle -> DrawLatex(0.18, 0.80, "Pb#minusPb, J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4");


    canvasV2JpsiVsPtCentr1030 -> SaveAs("QM2025_preliminary/v2JpsiVsPtCentr1030.pdf");
    canvasV2JpsiVsPtCentr1030Run2VsRun3 -> SaveAs("QM2025_preliminary/v2JpsiVsPtCentr1030Run2VsRun3.pdf");
    canvasV2JpsiVsPt3050 -> SaveAs("QM2025_preliminary/v2JpsiVsPt3050.pdf");
    canvasV2JpsiVsPt3050Run2VsRun3 -> SaveAs("QM2025_preliminary/v2JpsiVsPt3050Run2VsRun3.pdf");
    canvasV2JpsiVsPtCentr5080 -> SaveAs("QM2025_preliminary/v2JpsiVsPtCentr5080.pdf");
    canvasV2JpsiVsPtColl -> SaveAs("QM2025_preliminary/v2JpsiVsPtCollection.pdf");
    canvasV2JpsiVsPtCollLargeBins -> SaveAs("QM2025_preliminary/v2JpsiVsPtCollectionLargeBins.pdf");
    canvasV2JpsiVsCentrPt05Run2VsRun3 -> SaveAs("QM2025_preliminary/v2JpsiVsCentr05Run2VsRun3.pdf");
    canvasV2JpsiVsCentrPt515Run2VsRun3 -> SaveAs("QM2025_preliminary/v2JpsiVsCentr515Run2VsRun3.pdf");
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