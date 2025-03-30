void LoadStyle();
void SetLegend(TLegend *);

inline void EvalXaxis(int nVarBins, double minVarX[], double maxVarX[], double varCentr[], double varWidthStat[], double varWidthSyst[], double systWidth) {
    for (int iVar = 0;iVar < nVarBins;iVar++) {
        varCentr[iVar] = (maxVarX[iVar] + minVarX[iVar]) / 2.;
        varWidthStat[iVar] = (maxVarX[iVar] - minVarX[iVar]) / 2.;
        varWidthSyst[iVar] = systWidth;
    }
}

inline void EvalXaxis2(int nVarBins, double minVarX[], double maxVarX[], double varCentr[], double varWidthStat[], double varWidthSyst[], double systWidth) {
    for (int iVar = 0;iVar < nVarBins;iVar++) {
        varCentr[iVar] = (maxVarX[iVar] + minVarX[iVar]) / 2.;
        varWidthStat[iVar] = (maxVarX[iVar] - minVarX[iVar]) / 2.;
        varWidthSyst[iVar] = systWidth;
    }
}

inline void SetGraph(auto *graph, Color_t mkrCol, int mkrSty, Color_t lnCol, int lnWidth, int fillSty, double alpha = 1) {
    graph -> SetMarkerColorAlpha(mkrCol, alpha);
    graph -> SetMarkerStyle(mkrSty);
    graph -> SetLineColorAlpha(lnCol, alpha);
    graph -> SetLineWidth(lnWidth);
    graph -> SetFillStyle(fillSty);
}

inline void GetDataFromTxt(string fInName, 
                           std::vector<double>& vecCentrX, std::vector<double>& vecStatX, std::vector<double>& vecSystX, 
                           std::vector<double>& vecCentrY, std::vector<double>& vecStatY, std::vector<double>& vecSystY, double systWidth) {

    std::ifstream fIn(fInName.c_str());
    if (!fIn) { std::cerr << "Error while opening fIn!" << std::endl;}

    double minX, maxX, centrY, statY, systY;
    while (fIn >> minX >> maxX >> centrY >> statY >> systY) {
        std::cout << minX << " " << maxX << " " << centrY << " " << statY << " " << systY << std::endl;
        vecCentrX.push_back((maxX + minX) / 2.);
        vecStatX.push_back((maxX - minX) / 2.);
        vecSystX.push_back(systWidth);
        vecCentrY.push_back(centrY);
        vecStatY.push_back(statY);
        vecSystY.push_back(systY);
    }

    fIn.close();
}

void plot_results_all_centr_and_pt() {
    // EP, List of input files
    string pathNarrowBinsEpPt1030 = "../../results/QM2025/EP/jpsi_v2_vs_pt_10_30_narrow_binning.txt";
    string pathNarrowBinsEpPt3050 = "../../results/QM2025/EP/jpsi_v2_vs_pt_30_50_narrow_binning.txt";
    string pathNarrowBinsEpPt5080 = "../../results/QM2025/EP/jpsi_v2_vs_pt_50_80_large_binning.txt";
    string pathLargeBinsEpPt1030 = "../../results/QM2025/EP/jpsi_v2_vs_pt_10_30_large_binning.txt";
    string pathLargeBinsEpPt3050 = "../../results/QM2025/EP/jpsi_v2_vs_pt_30_50_large_binning.txt";
    string pathLargeBinsEpPt5080 = "../../results/QM2025/EP/jpsi_v2_vs_pt_50_80_large_binning.txt";
    string pathEpCentr05 = "../../results/QM2025/EP/jpsi_v2_vs_centr_0_5.txt";
    string pathEpCentr515 = "../../results/QM2025/EP/jpsi_v2_vs_centr_5_15.txt";

    // SP, List of input files
    string pathNarrowBinsSpPt010 = "../../results/QM2025/SP/jpsi_v2_vs_pt_0_10_narrow_binning.txt";
    string pathNarrowBinsSpPt1030 = "../../results/QM2025/SP/jpsi_v2_vs_pt_10_30_narrow_binning.txt";
    string pathNarrowBinsSpPt3050 = "../../results/QM2025/SP/jpsi_v2_vs_pt_30_50_narrow_binning.txt";
    string pathSpCentr05 = "../../results/QM2025/SP/jpsi_v2_vs_centr_0_5.txt";
    string pathSpCentr520 = "../../results/QM2025/SP/jpsi_v2_vs_centr_5_20.txt";

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

    //////////////////////////////////////////
    // J/psi v2 vs pT / centrality vs Run 2 //
    //////////////////////////////////////////
    const int nPtBinsRun2 = 10;
    double ptCentrRun2[] = {0.64, 1.49, 2.47, 3.46, 4.45, 5.45, 6.819, 8.835, 10.84, 14.25};
    double ptWidthStatRun2[] = {0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
    double ptWidthSystRun2[] = {0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30};

    double v2JpsiVsPtCentr010Run2[] = {0.03,0.034, 0.022, 0.053, 0.043, 0.05, 0.045, 0.0006, 0.068, 0.002};
    double statV2JpsiVsPtCentr010Run2[] = {0.013, 0.01, 0.011, 0.013, 0.016, 0.019, 0.019, 0.028, 0.046, 0.059};
    double systV2JpsiVsPtCentr010Run2[] = {0.0041243, 0.0031563, 0.0040117, 0.0030156, 0.0017944, 0.0025338, 0.0033808, 0.0051196, 0.0059211, 0.0051197};

    TGraphErrors *graStatV2JpsiVsPt010Run2 = new TGraphErrors(nPtBinsRun2, ptCentrRun2, v2JpsiVsPtCentr010Run2, ptWidthStatRun2, statV2JpsiVsPtCentr010Run2);
    SetGraph(graStatV2JpsiVsPt010Run2, kGray+2, 20, kGray+2, 2, 0, 0.7);
    TGraphErrors *graSystV2JpsiVsPt010Run2 = new TGraphErrors(nPtBinsRun2, ptCentrRun2, v2JpsiVsPtCentr010Run2, ptWidthSystRun2, systV2JpsiVsPtCentr010Run2);
    SetGraph(graSystV2JpsiVsPt010Run2, kGray+2, 20, kGray+2, 2, 0, 0.7);

    double v2JpsiVsPtCentr1030Run2[] = {0.011, 0.043, 0.074, 0.088, 0.085, 0.103, 0.083, 0.100, 0.049, 0.022};
    double statV2JpsiVsPtCentr1030Run2[] = {0.0085, 0.0069, 0.0069, 0.0077, 0.009, 0.011, 0.011, 0.018, 0.028, 0.032};
    double systV2JpsiVsPtCentr1030Run2[] = {0.0038391, 0.0036633, 0.004898, 0.0035068, 0.0037855, 0.0029726, 0.0036802, 0.0075789, 0.0093488, 0.0091828};

    TGraphErrors *graStatV2JpsiVsPt1030Run2 = new TGraphErrors(nPtBinsRun2, ptCentrRun2, v2JpsiVsPtCentr1030Run2, ptWidthStatRun2, statV2JpsiVsPtCentr1030Run2);
    SetGraph(graStatV2JpsiVsPt1030Run2, kGray+2, 20, kGray+2, 2, 0, 0.7);
    TGraphErrors *graSystV2JpsiVsPt1030Run2 = new TGraphErrors(nPtBinsRun2, ptCentrRun2, v2JpsiVsPtCentr1030Run2, ptWidthSystRun2, systV2JpsiVsPtCentr1030Run2);
    SetGraph(graSystV2JpsiVsPt1030Run2, kGray+2, 20, kGray+2, 2, 0, 0.7);

    double v2JpsiVsPt3050Run2[] = {0.0008, 0.029, 0.067, 0.099, 0.098, 0.101, 0.098, 0.092, 0.055, 0.026};
    double statV2JpsiVsPt3050Run2[] = {0.011, 0.0091, 0.0089, 0.0095, 0.011, 0.013, 0.023, 0.022, 0.037, 0.039};
    double systV2JpsiVsPt3050Run2[] = {0.0032389, 0.0032904, 0.003359, 0.0043267, 0.0065719, 0.0066256, 0.0065651, 0.0067724, 0.0075293, 0.0093145};

    TGraphErrors *graStatV2JpsiVsPt3050Run2 = new TGraphErrors(nPtBinsRun2, ptCentrRun2, v2JpsiVsPt3050Run2, ptWidthStatRun2, statV2JpsiVsPt3050Run2);
    SetGraph(graStatV2JpsiVsPt3050Run2, kGray+2, 20, kGray+2, 2, 0, 0.7);
    TGraphErrors *graSystV2JpsiVsPt3050Run2 = new TGraphErrors(nPtBinsRun2, ptCentrRun2, v2JpsiVsPt3050Run2, ptWidthSystRun2, systV2JpsiVsPt3050Run2);
    SetGraph(graSystV2JpsiVsPt3050Run2, kGray+2, 20, kGray+2, 2, 0, 0.7);

    const int nCentrBins05Run2 = 7;
    double centrCentr05Run2[] = {5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 75.0};
    double centrWidthStatPt05Run2[] = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 15.0};
    double centrWidthSystPt05Run2[] = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};

    double v2JpsiVsCentr05Run2[] = {0.036, 0.05, 0.056, 0.055, 0.031, 0.065, 0.045};
    double statV2JpsiVsCentr05Run2[] = {0.0093, 0.0081, 0.0084, 0.0099, 0.011, 0.014, 0.019};
    double systV2JpsiVsCentr05Run2[] = {0.0031078, 0.0044159, 0.0057615, 0.0053348, 0.0065608, 0.0065192, 0.0082189};

    TGraphErrors *graStatV2JpsiVsCentr05Run2 = new TGraphErrors(nCentrBins05Run2, centrCentr05Run2, v2JpsiVsCentr05Run2, centrWidthStatPt05Run2, statV2JpsiVsCentr05Run2);
    SetGraph(graStatV2JpsiVsCentr05Run2, kGray+2, 20, kGray+2, 2, 0, 0.7);
    TGraphErrors *graSystV2JpsiVsCentr05Run2 = new TGraphErrors(nCentrBins05Run2, centrCentr05Run2, v2JpsiVsCentr05Run2, centrWidthSystPt05Run2, systV2JpsiVsCentr05Run2);
    SetGraph(graSystV2JpsiVsCentr05Run2, kGray+2, 20, kGray+2, 2, 0, 0.7);

    const int nCentrBins515Run2 = 7;
    double centrCentr515Run2[] = {5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 75.0};
    double centrWidthStatPt515Run2[] = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 15.0};
    double centrWidthSystPt515Run2[] = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};

    double v2JpsiVsCentr515Run2[] = {0.04, 0.064, 0.105, 0.096, 0.092, 0.106, 0.1};
    double statV2JpsiVsCentr515Run2[] = {0.02, 0.015, 0.015, 0.016, 0.019, 0.025, 0.036};
    double systV2JpsiVsCentr515Run2[] = {0.0044385, 0.0044843, 0.0065419, 0.0064786, 0.0083954, 0.0085108, 0.011057};

    TGraphErrors *graStatV2JpsiVsCentr515Run2 = new TGraphErrors(nCentrBins515Run2, centrCentr515Run2, v2JpsiVsCentr515Run2, centrWidthStatPt515Run2, statV2JpsiVsCentr515Run2);
    SetGraph(graStatV2JpsiVsCentr515Run2, kGray+2, 20, kGray+2, 2, 0, 0.7);
    TGraphErrors *graSystV2JpsiVsCentr515Run2 = new TGraphErrors(nCentrBins515Run2, centrCentr515Run2, v2JpsiVsCentr515Run2, centrWidthSystPt515Run2, systV2JpsiVsCentr515Run2);
    SetGraph(graSystV2JpsiVsCentr515Run2, kGray+2, 20, kGray+2, 2, 0, 0.7);

    /////////////////////////////////////
    // J/psi v2 vs pT / centr vs Run 2 //
    /////////////////////////////////////
    // 0-10%
    // SP, narrow binning
    std::vector<double> ptCentr010SpRun3, ptWidthStat010SpRun3, ptWidthSyst010SpRun3, v2JpsiVsPt010SpRun3, statV2JpsiVsPt010SpRun3, systV2JpsiVsPt010SpRun3;
    GetDataFromTxt(pathNarrowBinsSpPt010.c_str(), ptCentr010SpRun3, ptWidthStat010SpRun3, ptWidthSyst010SpRun3, v2JpsiVsPt010SpRun3, statV2JpsiVsPt010SpRun3, systV2JpsiVsPt010SpRun3, 0.15);
    TGraphErrors *graStatV2JpsiSpVsPt010Run3 = new TGraphErrors(ptCentr010SpRun3.size(), &(ptCentr010SpRun3[0]), &(v2JpsiVsPt010SpRun3[0]), &(ptWidthStat010SpRun3[0]), &(statV2JpsiVsPt010SpRun3[0]));
    SetGraph(graStatV2JpsiSpVsPt010Run3, kRed+1, 20, kRed+1, 2, 0);
    TGraphErrors *graSystV2JpsiSpVsPt010Run3 = new TGraphErrors(ptCentr010SpRun3.size(), &(ptCentr010SpRun3[0]), &(v2JpsiVsPt010SpRun3[0]), &(ptWidthSyst010SpRun3[0]), &(systV2JpsiVsPt010SpRun3[0]));
    SetGraph(graSystV2JpsiSpVsPt010Run3, kRed+1, 20, kRed+1, 2, 0);

    TGraphErrors *graStatV2JpsiSpVsPt010Run3CollPlot = new TGraphErrors(ptCentr010SpRun3.size(), &(ptCentr010SpRun3[0]), &(v2JpsiVsPt010SpRun3[0]), &(ptWidthStat010SpRun3[0]), &(statV2JpsiVsPt010SpRun3[0]));
    SetGraph(graStatV2JpsiSpVsPt010Run3CollPlot, kOrange+7, 20, kOrange+7, 2, 0);
    TGraphErrors *graSystV2JpsiSpVsPt010Run3CollPlot = new TGraphErrors(ptCentr010SpRun3.size(), &(ptCentr010SpRun3[0]), &(v2JpsiVsPt010SpRun3[0]), &(ptWidthSyst010SpRun3[0]), &(systV2JpsiVsPt010SpRun3[0]));
    SetGraph(graSystV2JpsiSpVsPt010Run3CollPlot, kOrange+7, 20, kOrange+7, 2, 0);

    // 10-30%
    // EP, narrow / large binning
    std::vector<double> ptCentr1030EpRun3, ptWidthStat1030EpRun3, ptWidthSyst1030EpRun3, v2JpsiVsPt1030EpRun3, statV2JpsiVsPt1030EpRun3, systV2JpsiVsPt1030EpRun3;
    GetDataFromTxt(pathNarrowBinsEpPt1030.c_str(), ptCentr1030EpRun3, ptWidthStat1030EpRun3, ptWidthSyst1030EpRun3, v2JpsiVsPt1030EpRun3, statV2JpsiVsPt1030EpRun3, systV2JpsiVsPt1030EpRun3, 0.15);
    TGraphErrors *graStatV2JpsiEpVsPt1030Run3 = new TGraphErrors(ptCentr1030EpRun3.size(), &(ptCentr1030EpRun3[0]), &(v2JpsiVsPt1030EpRun3[0]), &(ptWidthStat1030EpRun3[0]), &(statV2JpsiVsPt1030EpRun3[0]));
    SetGraph(graStatV2JpsiEpVsPt1030Run3, kRed+1, 20, kRed+1, 2, 0);
    TGraphErrors *graSystV2JpsiEpVsPt1030Run3 = new TGraphErrors(ptCentr1030EpRun3.size(), &(ptCentr1030EpRun3[0]), &(v2JpsiVsPt1030EpRun3[0]), &(ptWidthSyst1030EpRun3[0]), &(systV2JpsiVsPt1030EpRun3[0]));
    SetGraph(graSystV2JpsiEpVsPt1030Run3, kRed+1, 20, kRed+1, 2, 0);

    std::vector<double> ptCentr1030EpRun3LargeBins, ptWidthStat1030EpRun3LargeBins, ptWidthSyst1030EpRun3LargeBins, v2JpsiVsPt1030EpRun3LargeBins, statV2JpsiVsPt1030EpRun3LargeBins, systV2JpsiVsPt1030EpRun3LargeBins;
    GetDataFromTxt(pathLargeBinsEpPt1030.c_str(), ptCentr1030EpRun3LargeBins, ptWidthStat1030EpRun3LargeBins, ptWidthSyst1030EpRun3LargeBins, v2JpsiVsPt1030EpRun3LargeBins, statV2JpsiVsPt1030EpRun3LargeBins, systV2JpsiVsPt1030EpRun3LargeBins, 0.15);
    TGraphErrors *graStatV2JpsiEpVsPt1030Run3LargeBins = new TGraphErrors(ptCentr1030EpRun3LargeBins.size(), &(ptCentr1030EpRun3LargeBins[0]), &(v2JpsiVsPt1030EpRun3LargeBins[0]), &(ptWidthStat1030EpRun3LargeBins[0]), &(statV2JpsiVsPt1030EpRun3LargeBins[0]));
    SetGraph(graStatV2JpsiEpVsPt1030Run3LargeBins, kRed+1, 20, kRed+1, 2, 0);
    TGraphErrors *graSystV2JpsiEpVsPt1030Run3LargeBins = new TGraphErrors(ptCentr1030EpRun3LargeBins.size(), &(ptCentr1030EpRun3LargeBins[0]), &(v2JpsiVsPt1030EpRun3LargeBins[0]), &(ptWidthSyst1030EpRun3LargeBins[0]), &(systV2JpsiVsPt1030EpRun3LargeBins[0]));
    SetGraph(graSystV2JpsiEpVsPt1030Run3LargeBins, kRed+1, 20, kRed+1, 2, 0);

    TGraphErrors *graStatV2JpsiEpVsPt1030Run3LargeBinsCollPlot = new TGraphErrors(ptCentr1030EpRun3LargeBins.size(), &(ptCentr1030EpRun3LargeBins[0]), &(v2JpsiVsPt1030EpRun3LargeBins[0]), &(ptWidthStat1030EpRun3LargeBins[0]), &(statV2JpsiVsPt1030EpRun3LargeBins[0]));
    SetGraph(graStatV2JpsiEpVsPt1030Run3LargeBinsCollPlot, kRed+1, 20, kRed+1, 2, 0);
    TGraphErrors *graSystV2JpsiEpVsPt1030Run3LargeBinsCollPlot = new TGraphErrors(ptCentr1030EpRun3LargeBins.size(), &(ptCentr1030EpRun3LargeBins[0]), &(v2JpsiVsPt1030EpRun3LargeBins[0]), &(ptWidthSyst1030EpRun3LargeBins[0]), &(systV2JpsiVsPt1030EpRun3LargeBins[0]));
    SetGraph(graSystV2JpsiEpVsPt1030Run3LargeBinsCollPlot, kRed+1, 20, kRed+1, 2, 0);

    // SP, narrow binning
    std::vector<double> ptCentr1030SpRun3, ptWidthStat1030SpRun3, ptWidthSyst1030SpRun3, v2JpsiVsPt1030SpRun3, statV2JpsiVsPt1030SpRun3, systV2JpsiVsPt1030SpRun3;
    GetDataFromTxt(pathNarrowBinsSpPt1030.c_str(), ptCentr1030SpRun3, ptWidthStat1030SpRun3, ptWidthSyst1030SpRun3, v2JpsiVsPt1030SpRun3, statV2JpsiVsPt1030SpRun3, systV2JpsiVsPt1030SpRun3, 0.15);
    TGraphErrors *graStatV2JpsiSpVsPt1030Run3 = new TGraphErrors(ptCentr1030SpRun3.size(), &(ptCentr1030SpRun3[0]), &(v2JpsiVsPt1030SpRun3[0]), &(ptWidthStat1030SpRun3[0]), &(statV2JpsiVsPt1030SpRun3[0]));
    SetGraph(graStatV2JpsiSpVsPt1030Run3, kRed+1, 20, kRed+1, 2, 0);
    TGraphErrors *graSystV2JpsiSpVsPt1030Run3 = new TGraphErrors(ptCentr1030SpRun3.size(), &(ptCentr1030SpRun3[0]), &(v2JpsiVsPt1030SpRun3[0]), &(ptWidthSyst1030SpRun3[0]), &(systV2JpsiVsPt1030SpRun3[0]));
    SetGraph(graSystV2JpsiSpVsPt1030Run3, kRed+1, 20, kRed+1, 2, 0);

    TGraphErrors *graStatV2JpsiSpVsPt1030Run3CollPlot = new TGraphErrors(ptCentr1030SpRun3.size(), &(ptCentr1030SpRun3[0]), &(v2JpsiVsPt1030SpRun3[0]), &(ptWidthStat1030SpRun3[0]), &(statV2JpsiVsPt1030SpRun3[0]));
    SetGraph(graStatV2JpsiSpVsPt1030Run3CollPlot, kRed+1, 20, kRed+1, 2, 0);
    TGraphErrors *graSystV2JpsiSpVsPt1030Run3CollPlot = new TGraphErrors(ptCentr1030SpRun3.size(), &(ptCentr1030SpRun3[0]), &(v2JpsiVsPt1030SpRun3[0]), &(ptWidthSyst1030SpRun3[0]), &(systV2JpsiVsPt1030SpRun3[0]));
    SetGraph(graSystV2JpsiSpVsPt1030Run3CollPlot, kRed+1, 20, kRed+1, 2, 0);

    // 30-50% 
    // narrow / large binning
    std::vector<double> ptCentr3050EpRun3, ptWidthStat3050EpRun3, ptWidthSyst3050EpRun3, v2JpsiVsPt3050EpRun3, statV2JpsiVsPt3050EpRun3, systV2JpsiVsPt3050EpRun3;
    GetDataFromTxt(pathNarrowBinsEpPt3050.c_str(), ptCentr3050EpRun3, ptWidthStat3050EpRun3, ptWidthSyst3050EpRun3, v2JpsiVsPt3050EpRun3, statV2JpsiVsPt3050EpRun3, systV2JpsiVsPt3050EpRun3, 0.15);
    TGraphErrors *graStatV2JpsiEpVsPt3050Run3 = new TGraphErrors(ptCentr3050EpRun3.size(), &(ptCentr3050EpRun3[0]), &(v2JpsiVsPt3050EpRun3[0]), &(ptWidthStat3050EpRun3[0]), &(statV2JpsiVsPt3050EpRun3[0]));
    SetGraph(graStatV2JpsiEpVsPt3050Run3, kRed+1, 20, kRed+1, 2, 0);
    TGraphErrors *graSystV2JpsiEpVsPt3050Run3 = new TGraphErrors(ptCentr3050EpRun3.size(), &(ptCentr3050EpRun3[0]), &(v2JpsiVsPt3050EpRun3[0]), &(ptWidthSyst3050EpRun3[0]), &(systV2JpsiVsPt3050EpRun3[0]));
    SetGraph(graSystV2JpsiEpVsPt3050Run3, kRed+1, 20, kRed+1, 2, 0);

    std::vector<double> ptCentr3050EpRun3LargeBins, ptWidthStat3050EpRun3LargeBins, ptWidthSyst3050EpRun3LargeBins, v2JpsiVsPt3050EpRun3LargeBins, statV2JpsiVsPt3050EpRun3LargeBins, systV2JpsiVsPt3050EpRun3LargeBins;
    GetDataFromTxt(pathLargeBinsEpPt3050.c_str(), ptCentr3050EpRun3LargeBins, ptWidthStat3050EpRun3LargeBins, ptWidthSyst3050EpRun3LargeBins, v2JpsiVsPt3050EpRun3LargeBins, statV2JpsiVsPt3050EpRun3LargeBins, systV2JpsiVsPt3050EpRun3LargeBins, 0.15);
    TGraphErrors *graStatV2JpsiEpVsPt3050Run3LargeBins = new TGraphErrors(ptCentr3050EpRun3LargeBins.size(), &(ptCentr3050EpRun3LargeBins[0]), &(v2JpsiVsPt3050EpRun3LargeBins[0]), &(ptWidthStat3050EpRun3LargeBins[0]), &(statV2JpsiVsPt3050EpRun3LargeBins[0]));
    SetGraph(graStatV2JpsiEpVsPt3050Run3LargeBins, kRed+1, 20, kRed+1, 2, 0);
    TGraphErrors *graSystV2JpsiEpVsPt3050Run3LargeBins = new TGraphErrors(ptCentr3050EpRun3LargeBins.size(), &(ptCentr3050EpRun3LargeBins[0]), &(v2JpsiVsPt3050EpRun3LargeBins[0]), &(ptWidthSyst3050EpRun3LargeBins[0]), &(systV2JpsiVsPt3050EpRun3LargeBins[0]));
    SetGraph(graSystV2JpsiEpVsPt3050Run3LargeBins, kRed+1, 20, kRed+1, 2, 0);

    TGraphErrors *graStatV2JpsiEpVsPt3050Run3LargeBinsCollPlot = new TGraphErrors(ptCentr3050EpRun3LargeBins.size(), &(ptCentr3050EpRun3LargeBins[0]), &(v2JpsiVsPt3050EpRun3LargeBins[0]), &(ptWidthStat3050EpRun3LargeBins[0]), &(statV2JpsiVsPt3050EpRun3LargeBins[0]));
    SetGraph(graStatV2JpsiEpVsPt3050Run3LargeBinsCollPlot, kAzure+4, 20, kAzure+4, 2, 0);
    TGraphErrors *graSystV2JpsiEpVsPt3050Run3LargeBinsCollPlot = new TGraphErrors(ptCentr3050EpRun3LargeBins.size(), &(ptCentr3050EpRun3LargeBins[0]), &(v2JpsiVsPt3050EpRun3LargeBins[0]), &(ptWidthSyst3050EpRun3LargeBins[0]), &(systV2JpsiVsPt3050EpRun3LargeBins[0]));
    SetGraph(graSystV2JpsiEpVsPt3050Run3LargeBinsCollPlot, kAzure+4, 20, kAzure+4, 2, 0);

    // SP, narrow binning
    std::vector<double> ptCentr3050SpRun3, ptWidthStat3050SpRun3, ptWidthSyst3050SpRun3, v2JpsiVsPt3050SpRun3, statV2JpsiVsPt3050SpRun3, systV2JpsiVsPt3050SpRun3;
    GetDataFromTxt(pathNarrowBinsSpPt3050.c_str(), ptCentr3050SpRun3, ptWidthStat3050SpRun3, ptWidthSyst3050SpRun3, v2JpsiVsPt3050SpRun3, statV2JpsiVsPt3050SpRun3, systV2JpsiVsPt3050SpRun3, 0.15);
    TGraphErrors *graStatV2JpsiSpVsPt3050Run3 = new TGraphErrors(ptCentr3050SpRun3.size(), &(ptCentr3050SpRun3[0]), &(v2JpsiVsPt3050SpRun3[0]), &(ptWidthStat3050SpRun3[0]), &(statV2JpsiVsPt3050SpRun3[0]));
    SetGraph(graStatV2JpsiSpVsPt3050Run3, kRed+1, 20, kRed+1, 2, 0);
    TGraphErrors *graSystV2JpsiSpVsPt3050Run3 = new TGraphErrors(ptCentr3050SpRun3.size(), &(ptCentr3050SpRun3[0]), &(v2JpsiVsPt3050SpRun3[0]), &(ptWidthSyst3050SpRun3[0]), &(systV2JpsiVsPt3050SpRun3[0]));
    SetGraph(graSystV2JpsiSpVsPt3050Run3, kRed+1, 20, kRed+1, 2, 0);

    TGraphErrors *graStatV2JpsiSpVsPt3050Run3CollPlot = new TGraphErrors(ptCentr3050SpRun3.size(), &(ptCentr3050SpRun3[0]), &(v2JpsiVsPt3050SpRun3[0]), &(ptWidthStat3050SpRun3[0]), &(statV2JpsiVsPt3050SpRun3[0]));
    SetGraph(graStatV2JpsiSpVsPt3050Run3CollPlot, kAzure+4, 20, kAzure+4, 2, 0);
    TGraphErrors *graSystV2JpsiSpVsPt3050Run3CollPlot = new TGraphErrors(ptCentr3050SpRun3.size(), &(ptCentr3050SpRun3[0]), &(v2JpsiVsPt3050SpRun3[0]), &(ptWidthSyst3050SpRun3[0]), &(systV2JpsiVsPt3050SpRun3[0]));
    SetGraph(graSystV2JpsiSpVsPt3050Run3CollPlot, kAzure+4, 20, kAzure+4, 2, 0);

    // 50-80% 
    // narrow / large binning
    std::vector<double> ptCentr5080EpRun3, ptWidthStat5080EpRun3, ptWidthSyst5080EpRun3, v2JpsiVsPt5080EpRun3, statV2JpsiVsPt5080EpRun3, systV2JpsiVsPt5080EpRun3;
    GetDataFromTxt(pathNarrowBinsEpPt5080.c_str(), ptCentr5080EpRun3, ptWidthStat5080EpRun3, ptWidthSyst5080EpRun3, v2JpsiVsPt5080EpRun3, statV2JpsiVsPt5080EpRun3, systV2JpsiVsPt5080EpRun3, 0.15);
    TGraphErrors *graStatV2JpsiEpVsPt5080Run3 = new TGraphErrors(ptCentr5080EpRun3.size(), &(ptCentr5080EpRun3[0]), &(v2JpsiVsPt5080EpRun3[0]), &(ptWidthStat5080EpRun3[0]), &(statV2JpsiVsPt5080EpRun3[0]));
    SetGraph(graStatV2JpsiEpVsPt5080Run3, kRed+1, 20, kRed+1, 2, 0);
    TGraphErrors *graSystV2JpsiEpVsPt5080Run3 = new TGraphErrors(ptCentr5080EpRun3.size(), &(ptCentr5080EpRun3[0]), &(v2JpsiVsPt5080EpRun3[0]), &(ptWidthSyst5080EpRun3[0]), &(systV2JpsiVsPt5080EpRun3[0]));
    SetGraph(graSystV2JpsiEpVsPt5080Run3, kRed+1, 20, kRed+1, 2, 0);

    TGraphErrors *graStatV2JpsiEpVsPt5080Run3CollPlot = new TGraphErrors(ptCentr5080EpRun3.size(), &(ptCentr5080EpRun3[0]), &(v2JpsiVsPt5080EpRun3[0]), &(ptWidthStat5080EpRun3[0]), &(statV2JpsiVsPt5080EpRun3[0]));
    SetGraph(graStatV2JpsiEpVsPt5080Run3CollPlot, kGreen+2, 20, kGreen+2, 2, 0);
    TGraphErrors *graSystV2JpsiEpVsPt5080Run3CollPlot = new TGraphErrors(ptCentr5080EpRun3.size(), &(ptCentr5080EpRun3[0]), &(v2JpsiVsPt5080EpRun3[0]), &(ptWidthSyst5080EpRun3[0]), &(systV2JpsiVsPt5080EpRun3[0]));
    SetGraph(graSystV2JpsiEpVsPt5080Run3CollPlot, kGreen+2, 20, kGreen+2, 2, 0);

    // 0 < pT < 5 GeV/c
    // EP
    std::vector<double> centrCenCentr05EpRun3, ptWidthStCentr05EpRun3, ptWidthSyCentr05EpRun3, v2JpsiVsCentr05EpRun3, statV2JpsiVsCentr05EpRun3, systV2JpsiVsCentr05EpRun3;
    GetDataFromTxt(pathEpCentr05.c_str(), centrCenCentr05EpRun3, ptWidthStCentr05EpRun3, ptWidthSyCentr05EpRun3, v2JpsiVsCentr05EpRun3, statV2JpsiVsCentr05EpRun3, systV2JpsiVsCentr05EpRun3, 1.5);
    TGraphErrors *graStatV2JpsiEpVsCentr05Run3 = new TGraphErrors(centrCenCentr05EpRun3.size(), &(centrCenCentr05EpRun3[0]), &(v2JpsiVsCentr05EpRun3[0]), &(ptWidthStCentr05EpRun3[0]), &(statV2JpsiVsCentr05EpRun3[0]));
    SetGraph(graStatV2JpsiEpVsCentr05Run3, kRed+1, 20, kRed+1, 2, 0);
    TGraphErrors *graSystV2JpsiEpVsCentr05Run3 = new TGraphErrors(centrCenCentr05EpRun3.size(), &(centrCenCentr05EpRun3[0]), &(v2JpsiVsCentr05EpRun3[0]), &(ptWidthSyCentr05EpRun3[0]), &(systV2JpsiVsCentr05EpRun3[0]));
    SetGraph(graSystV2JpsiEpVsCentr05Run3, kRed+1, 20, kRed+1, 2, 0);

    // SP
    std::vector<double> centrCenCentr05SpRun3, ptWidthStCentr05SpRun3, ptWidthSyCentr05SpRun3, v2JpsiVsCentr05SpRun3, statV2JpsiVsCentr05SpRun3, systV2JpsiVsCentr05SpRun3;
    GetDataFromTxt(pathSpCentr05.c_str(), centrCenCentr05SpRun3, ptWidthStCentr05SpRun3, ptWidthSyCentr05SpRun3, v2JpsiVsCentr05SpRun3, statV2JpsiVsCentr05SpRun3, systV2JpsiVsCentr05SpRun3, 1.5);
    TGraphErrors *graStatV2JpsiSpVsCentr05Run3 = new TGraphErrors(centrCenCentr05SpRun3.size(), &(centrCenCentr05SpRun3[0]), &(v2JpsiVsCentr05SpRun3[0]), &(ptWidthStCentr05SpRun3[0]), &(statV2JpsiVsCentr05SpRun3[0]));
    SetGraph(graStatV2JpsiSpVsCentr05Run3, kRed+1, 20, kRed+1, 2, 0);
    TGraphErrors *graSystV2JpsiSpVsCentr05Run3 = new TGraphErrors(centrCenCentr05SpRun3.size(), &(centrCenCentr05SpRun3[0]), &(v2JpsiVsCentr05SpRun3[0]), &(ptWidthSyCentr05SpRun3[0]), &(systV2JpsiVsCentr05SpRun3[0]));
    SetGraph(graSystV2JpsiSpVsCentr05Run3, kRed+1, 20, kRed+1, 2, 0);

    // 5 < pT < 15 GeV/c, EP
    std::vector<double> centrCenCentr515EpRun3, ptWidthStCentr515EpRun3, ptWidthSyCentr515EpRun3, v2JpsiVsCentr515EpRun3, statV2JpsiVsCentr515EpRun3, systV2JpsiVsCentr515EpRun3;
    GetDataFromTxt(pathEpCentr515.c_str(), centrCenCentr515EpRun3, ptWidthStCentr515EpRun3, ptWidthSyCentr515EpRun3, v2JpsiVsCentr515EpRun3, statV2JpsiVsCentr515EpRun3, systV2JpsiVsCentr515EpRun3, 1.5);
    TGraphErrors *graStatV2JpsiEpVsCentr515Run3 = new TGraphErrors(centrCenCentr515EpRun3.size(), &(centrCenCentr515EpRun3[0]), &(v2JpsiVsCentr515EpRun3[0]), &(ptWidthStCentr515EpRun3[0]), &(statV2JpsiVsCentr515EpRun3[0]));
    SetGraph(graStatV2JpsiEpVsCentr515Run3, kRed+1, 20, kRed+1, 2, 0);
    TGraphErrors *graSystV2JpsiEpVsCentr515Run3 = new TGraphErrors(centrCenCentr515EpRun3.size(), &(centrCenCentr515EpRun3[0]), &(v2JpsiVsCentr515EpRun3[0]), &(ptWidthSyCentr515EpRun3[0]), &(systV2JpsiVsCentr515EpRun3[0]));
    SetGraph(graSystV2JpsiEpVsCentr515Run3, kRed+1, 20, kRed+1, 2, 0);

    // 5 < pT < 20 GeV/c, SP
    std::vector<double> centrCenCentr520SpRun3, ptWidthStCentr520SpRun3, ptWidthSyCentr520SpRun3, v2JpsiVsCentr520SpRun3, statV2JpsiVsCentr520SpRun3, systV2JpsiVsCentr520SpRun3;
    GetDataFromTxt(pathSpCentr520.c_str(), centrCenCentr520SpRun3, ptWidthStCentr520SpRun3, ptWidthSyCentr520SpRun3, v2JpsiVsCentr520SpRun3, statV2JpsiVsCentr520SpRun3, systV2JpsiVsCentr520SpRun3, 1.5);
    TGraphErrors *graStatV2JpsiSpVsCentr520Run3 = new TGraphErrors(centrCenCentr520SpRun3.size(), &(centrCenCentr520SpRun3[0]), &(v2JpsiVsCentr520SpRun3[0]), &(ptWidthStCentr520SpRun3[0]), &(statV2JpsiVsCentr520SpRun3[0]));
    SetGraph(graStatV2JpsiSpVsCentr520Run3, kRed+1, 20, kRed+1, 2, 0);
    TGraphErrors *graSystV2JpsiSpVsCentr520Run3 = new TGraphErrors(centrCenCentr520SpRun3.size(), &(centrCenCentr520SpRun3[0]), &(v2JpsiVsCentr520SpRun3[0]), &(ptWidthSyCentr520SpRun3[0]), &(systV2JpsiVsCentr520SpRun3[0]));
    SetGraph(graSystV2JpsiSpVsCentr520Run3, kRed+1, 20, kRed+1, 2, 0);

    //******************************************************************************************//
    // J/psi v2 vs pT
    //******************************************************************************************//
    // 0-10%, SP
    TCanvas *canvasV2JpsiSpVsPt010 = new TCanvas("canvasV2JpsiSpVsPt010", "", 800, 600);
    canvasV2JpsiSpVsPt010 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiSpVsPt010 = new TH2D("histGridV2JpsiSpVsPt010", "", 100, -0.5, 15.5, 100, -0.07, 0.3);
    histGridV2JpsiSpVsPt010 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histGridV2JpsiSpVsPt010 -> GetYaxis() -> SetTitle("#it{v}_{2} {SP, |#Delta#eta| > 1.7}");
    histGridV2JpsiSpVsPt010 -> Draw();
    lineZeroV2Pt -> Draw();
    graStatV2JpsiSpVsPt010Run3 -> Draw("EP SAME");
    graSystV2JpsiSpVsPt010Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsPt010 = new TLegend(0.65, 0.62, 0.85, 0.82, " ", "brNDC");
    SetLegend(legendV2JpsiVsPt010);
    legendV2JpsiVsPt010 -> AddEntry(graStatV2JpsiSpVsPt010Run3, "Stat. Uncert.", "EP");
    legendV2JpsiVsPt010 -> AddEntry(graSystV2JpsiSpVsPt010Run3, "Syst. Uncert.", "F");
    legendV2JpsiVsPt010 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary");
    latexTitle -> DrawLatex(0.18, 0.80, "Pb#minusPb,#kern[0.3]{#sqrt{#it{s}_{NN}}} = 5.36 TeV, 0#minus10\%");
    latexTitle -> DrawLatex(0.18, 0.73, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4");

    // 10-30% 
    // EP
    TCanvas *canvasV2JpsiEpVsPt1030 = new TCanvas("canvasV2JpsiEpVsPt1030", "", 800, 600);
    canvasV2JpsiEpVsPt1030 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiEpVsPt1030 = new TH2D("histGridV2JpsiEpVsPt1030", "", 100, -0.5, 15.5, 100, -0.07, 0.3);
    histGridV2JpsiEpVsPt1030 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histGridV2JpsiEpVsPt1030 -> GetYaxis() -> SetTitle("#it{v}_{2} {EP, |#Delta#eta| > 1.7}");
    histGridV2JpsiEpVsPt1030 -> Draw();
    lineZeroV2Pt -> Draw();
    graStatV2JpsiEpVsPt1030Run3 -> Draw("EP SAME");
    graSystV2JpsiEpVsPt1030Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsPt1030 = new TLegend(0.65, 0.62, 0.85, 0.82, " ", "brNDC");
    SetLegend(legendV2JpsiVsPt1030);
    legendV2JpsiVsPt1030 -> AddEntry(graStatV2JpsiEpVsPt1030Run3, "Stat. Uncert.", "EP");
    legendV2JpsiVsPt1030 -> AddEntry(graSystV2JpsiEpVsPt1030Run3, "Syst. Uncert.", "F");
    legendV2JpsiVsPt1030 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary");
    latexTitle -> DrawLatex(0.18, 0.80, "Pb#minusPb,#kern[0.3]{#sqrt{#it{s}_{NN}}} = 5.36 TeV, 10#minus30\%");
    latexTitle -> DrawLatex(0.18, 0.73, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4");

    // SP
    TCanvas *canvasV2JpsiSpVsPt1030 = new TCanvas("canvasV2JpsiSpVsPt1030", "", 800, 600);
    canvasV2JpsiSpVsPt1030 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiSpVsPt1030 = new TH2D("histGridV2JpsiSpVsPt1030", "", 100, -0.5, 15.5, 100, -0.07, 0.3);
    histGridV2JpsiSpVsPt1030 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histGridV2JpsiSpVsPt1030 -> GetYaxis() -> SetTitle("#it{v}_{2} {SP, |#Delta#eta| > 1.7}");
    histGridV2JpsiSpVsPt1030 -> Draw();
    lineZeroV2Pt -> Draw();
    graStatV2JpsiSpVsPt1030Run3 -> Draw("EP SAME");
    graSystV2JpsiSpVsPt1030Run3 -> Draw("E2 SAME");
    legendV2JpsiVsPt1030 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary");
    latexTitle -> DrawLatex(0.18, 0.80, "Pb#minusPb,#kern[0.3]{#sqrt{#it{s}_{NN}}} = 5.36 TeV, 10#minus30\%");
    latexTitle -> DrawLatex(0.18, 0.73, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4");

    // 30-50%
    // EP
    TCanvas *canvasV2JpsiEpVsPt3050 = new TCanvas("canvasV2JpsiEpVsPt3050", "", 800, 600);
    canvasV2JpsiEpVsPt3050 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiEpVsPt3050 = new TH2D("histGridV2JpsiEpVsPt3050", "", 100, -0.5, 15.5, 100, -0.07, 0.3);
    histGridV2JpsiEpVsPt3050 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histGridV2JpsiEpVsPt3050 -> GetYaxis() -> SetTitle("#it{v}_{2} {EP, |#Delta#eta| > 1.7}");
    histGridV2JpsiEpVsPt3050 -> Draw();
    lineZeroV2Pt -> Draw();
    graStatV2JpsiEpVsPt3050Run3 -> Draw("EP SAME");
    graSystV2JpsiEpVsPt3050Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsPt3050 = new TLegend(0.65, 0.62, 0.85, 0.82, " ", "brNDC");
    SetLegend(legendV2JpsiVsPt3050);
    legendV2JpsiVsPt3050 -> AddEntry(graStatV2JpsiEpVsPt3050Run3, "Stat. Uncert.", "EP");
    legendV2JpsiVsPt3050 -> AddEntry(graSystV2JpsiEpVsPt3050Run3, "Syst. Uncert.", "F");
    legendV2JpsiVsPt3050 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary");
    latexTitle -> DrawLatex(0.18, 0.80, "Pb#minusPb,#kern[0.3]{#sqrt{#it{s}_{NN}}} = 5.36 TeV, 30#minus50\%");
    latexTitle -> DrawLatex(0.18, 0.73, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4");

    // SP
    TCanvas *canvasV2JpsiSpVsPt3050 = new TCanvas("canvasV2JpsiSpVsPt3050", "", 800, 600);
    canvasV2JpsiSpVsPt3050 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiSpVsPt3050 = new TH2D("histGridV2JpsiSpVsPt3050", "", 100, -0.5, 15.5, 100, -0.07, 0.3);
    histGridV2JpsiSpVsPt3050 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histGridV2JpsiSpVsPt3050 -> GetYaxis() -> SetTitle("#it{v}_{2} {SP, |#Delta#eta| > 1.7}");
    histGridV2JpsiSpVsPt3050 -> Draw();
    lineZeroV2Pt -> Draw();
    graStatV2JpsiSpVsPt3050Run3 -> Draw("EP SAME");
    graSystV2JpsiSpVsPt3050Run3 -> Draw("E2 SAME");
    legendV2JpsiVsPt3050 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary");
    latexTitle -> DrawLatex(0.18, 0.80, "Pb#minusPb,#kern[0.3]{#sqrt{#it{s}_{NN}}} = 5.36 TeV, 30#minus50\%");
    latexTitle -> DrawLatex(0.18, 0.73, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4");

    // 50-80%
    TCanvas *canvasV2JpsiEpVsPt5080 = new TCanvas("canvasV2JpsiEpVsPt5080", "", 800, 600);
    canvasV2JpsiEpVsPt5080 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiVsPt5080 = new TH2D("histGridV2JpsiVsPt5080", "", 100, -0.5, 15.5, 100, -0.07, 0.3);
    histGridV2JpsiVsPt5080 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histGridV2JpsiVsPt5080 -> GetYaxis() -> SetTitle("#it{v}_{2} {EP, |#Delta#eta| > 1.7}");
    histGridV2JpsiVsPt5080 -> Draw();
    lineZeroV2Pt -> Draw();
    graStatV2JpsiEpVsPt5080Run3 -> Draw("EP SAME");
    graSystV2JpsiEpVsPt5080Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsPt5080 = new TLegend(0.65, 0.62, 0.85, 0.82, " ", "brNDC");
    SetLegend(legendV2JpsiVsPt5080);
    legendV2JpsiVsPt5080 -> AddEntry(graStatV2JpsiEpVsPt5080Run3, "Stat. Uncert.", "EP");
    legendV2JpsiVsPt5080 -> AddEntry(graSystV2JpsiEpVsPt5080Run3, "Syst. Uncert.", "F");
    legendV2JpsiVsPt5080 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary");
    latexTitle -> DrawLatex(0.18, 0.80, "Pb#minusPb,#kern[0.3]{#sqrt{#it{s}_{NN}}} = 5.36 TeV, 50#minus80\%");
    latexTitle -> DrawLatex(0.18, 0.73, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4");

    //******************************************************************************************//
    // J/psi v2 vs centrality
    //******************************************************************************************//
    // 0 < pT < 5 GeV/c
    // EP
    TCanvas *canvasV2JpsiEpVsCentr05 = new TCanvas("canvasV2JpsiEpVsCentr05", "", 800, 600);
    canvasV2JpsiEpVsCentr05 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiEpVsCentr05 = new TH2D("histGridV2JpsiEpVsCentr05", "", 100, -0.5, 90.5, 100, -0.05, 0.20);
    histGridV2JpsiEpVsCentr05 -> GetXaxis() -> SetTitle("Centrality (%)");
    histGridV2JpsiEpVsCentr05 -> GetYaxis() -> SetTitle("#it{v}_{2} {EP, |#Delta#eta| > 1.7}");
    histGridV2JpsiEpVsCentr05 -> Draw();
    lineZeroV2Centr -> Draw();
    graStatV2JpsiEpVsCentr05Run3 -> Draw("EP SAME");
    graSystV2JpsiEpVsCentr05Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsCentr05 = new TLegend(0.20, 0.60, 0.45, 0.80, " ", "brNDC");
    SetLegend(legendV2JpsiVsCentr05);
    legendV2JpsiVsCentr05 -> AddEntry(graStatV2JpsiEpVsCentr05Run3, "Stat. Uncert.", "EP");
    legendV2JpsiVsCentr05 -> AddEntry(graSystV2JpsiEpVsCentr05Run3, "Syst. Uncert.", "F");
    legendV2JpsiVsCentr05 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb");
    latexTitle -> DrawLatex(0.18, 0.80, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4, #it{p}_{T} < 5 GeV/#it{c}");

    // SP
    TCanvas *canvasV2JpsiSpVsCentr05 = new TCanvas("canvasV2JpsiSpVsCentr05", "", 800, 600);
    canvasV2JpsiSpVsCentr05 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiSpVsCentr05 = new TH2D("histGridV2JpsiSpVsCentr05", "", 100, -0.5, 90.5, 100, -0.05, 0.20);
    histGridV2JpsiSpVsCentr05 -> GetXaxis() -> SetTitle("Centrality (%)");
    histGridV2JpsiSpVsCentr05 -> GetYaxis() -> SetTitle("#it{v}_{2} {SP, |#Delta#eta| > 1.7}");
    histGridV2JpsiSpVsCentr05 -> Draw();
    lineZeroV2Centr -> Draw();
    graStatV2JpsiSpVsCentr05Run3 -> Draw("EP SAME");
    graSystV2JpsiSpVsCentr05Run3 -> Draw("E2 SAME");
    legendV2JpsiVsCentr05 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb");
    latexTitle -> DrawLatex(0.18, 0.80, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4, #it{p}_{T} < 5 GeV/#it{c}");

    // 5 < pT < 15 GeV/c, EP
    TCanvas *canvasV2JpsiEpVsCentr515 = new TCanvas("canvasV2JpsiEpVsCentr515", "", 800, 600);
    canvasV2JpsiEpVsCentr515 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiEpVsCentr515 = new TH2D("histGridV2JpsiEpVsCentr515", "", 100, -0.5, 90.5, 100, -0.05, 0.30);
    histGridV2JpsiEpVsCentr515 -> GetXaxis() -> SetTitle("Centrality (%)");
    histGridV2JpsiEpVsCentr515 -> GetYaxis() -> SetTitle("#it{v}_{2} {EP, |#Delta#eta| > 1.7}");
    histGridV2JpsiEpVsCentr515 -> Draw();
    lineZeroV2Centr -> Draw();
    graStatV2JpsiEpVsCentr515Run3 -> Draw("EP SAME");
    graSystV2JpsiEpVsCentr515Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsCentr515 = new TLegend(0.20, 0.60, 0.45, 0.80, " ", "brNDC");
    SetLegend(legendV2JpsiVsCentr515);
    legendV2JpsiVsCentr515 -> AddEntry(graStatV2JpsiEpVsCentr515Run3, "Stat. Uncert.", "EP");
    legendV2JpsiVsCentr515 -> AddEntry(graSystV2JpsiEpVsCentr515Run3, "Syst. Uncert.", "F");
    legendV2JpsiVsCentr515 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb");
    latexTitle -> DrawLatex(0.18, 0.80, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4, 5 < #it{p}_{T} < 15 GeV/#it{c}");

    // 5 < pT < 20 GeV/c, SP
    TCanvas *canvasV2JpsiSpVsCentr520 = new TCanvas("canvasV2JpsiSpVsCentr520", "", 800, 600);
    canvasV2JpsiSpVsCentr520 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiSpVsCentr520 = new TH2D("histGridV2JpsiSpVsCentr520", "", 100, -0.5, 90.5, 100, -0.05, 0.30);
    histGridV2JpsiSpVsCentr520 -> GetXaxis() -> SetTitle("Centrality (%)");
    histGridV2JpsiSpVsCentr520 -> GetYaxis() -> SetTitle("#it{v}_{2} {SP, |#Delta#eta| > 1.7}");
    histGridV2JpsiSpVsCentr520 -> Draw();
    lineZeroV2Centr -> Draw();
    graStatV2JpsiSpVsCentr520Run3 -> Draw("EP SAME");
    graSystV2JpsiSpVsCentr520Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsCentr520 = new TLegend(0.20, 0.60, 0.45, 0.80, " ", "brNDC");
    SetLegend(legendV2JpsiVsCentr520);
    legendV2JpsiVsCentr520 -> AddEntry(graStatV2JpsiSpVsCentr520Run3, "Stat. Uncert.", "EP");
    legendV2JpsiVsCentr520 -> AddEntry(graSystV2JpsiSpVsCentr520Run3, "Syst. Uncert.", "F");
    legendV2JpsiVsCentr520 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb");
    latexTitle -> DrawLatex(0.18, 0.80, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4, 5 < #it{p}_{T} < 20 GeV/#it{c}");

    //******************************************************************************************//
    // J/psi v2 vs pT vs Run 2
    //******************************************************************************************//
    // 0-10%, SP 
    TCanvas *canvasV2JpsiSpVsPt010Run2VsRun3 = new TCanvas("canvasV2JpsiSpVsPt010Run2VsRun3", "", 800, 600);
    canvasV2JpsiSpVsPt010Run2VsRun3 -> SetTicks(1, 1);

    histGridV2JpsiSpVsPt010 -> Draw();
    lineZeroV2Pt -> Draw();
    graStatV2JpsiVsPt010Run2 -> Draw("EP SAME");
    graSystV2JpsiVsPt010Run2 -> Draw("E2 SAME");
    graStatV2JpsiSpVsPt010Run3 -> Draw("EP SAME");
    graSystV2JpsiSpVsPt010Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsPt010Run2VsSpRun3 = new TLegend(0.20, 0.70, 0.85, 0.85, " ", "brNDC");
    SetLegend(legendV2JpsiVsPt010Run2VsSpRun3);
    legendV2JpsiVsPt010Run2VsSpRun3 -> SetNColumns(2);
    legendV2JpsiVsPt010Run2VsSpRun3 -> AddEntry(graSystV2JpsiVsPt010Run2, "#sqrt{#it{s}_{NN}} = 5.02 TeV", "PF");
    legendV2JpsiVsPt010Run2VsSpRun3 -> AddEntry(graSystV2JpsiSpVsPt010Run3, "#sqrt{#it{s}_{NN}} = 5.36 TeV", "PF"); 
    legendV2JpsiVsPt010Run2VsSpRun3 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb");
    latexTitle -> DrawLatex(0.18, 0.80, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4, 0#minus10\%");

    // 10-30% 
    // EP
    TCanvas *canvasV2JpsiEpVsPt1030Run2VsRun3 = new TCanvas("canvasV2JpsiEpVsPt1030Run2VsRun3", "", 800, 600);
    canvasV2JpsiEpVsPt1030Run2VsRun3 -> SetTicks(1, 1);

    histGridV2JpsiEpVsPt1030 -> Draw();
    lineZeroV2Pt -> Draw();
    graStatV2JpsiVsPt1030Run2 -> Draw("EP SAME");
    graSystV2JpsiVsPt1030Run2 -> Draw("E2 SAME");
    graStatV2JpsiEpVsPt1030Run3 -> Draw("EP SAME");
    graSystV2JpsiEpVsPt1030Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsPt1030Run2VsEpRun3 = new TLegend(0.20, 0.70, 0.85, 0.85, " ", "brNDC");
    SetLegend(legendV2JpsiVsPt1030Run2VsEpRun3);
    legendV2JpsiVsPt1030Run2VsEpRun3 -> SetNColumns(2);
    legendV2JpsiVsPt1030Run2VsEpRun3 -> AddEntry(graSystV2JpsiVsPt1030Run2, "#sqrt{#it{s}_{NN}} = 5.02 TeV", "PF");
    legendV2JpsiVsPt1030Run2VsEpRun3 -> AddEntry(graSystV2JpsiEpVsPt1030Run3, "#sqrt{#it{s}_{NN}} = 5.36 TeV", "PF"); 
    legendV2JpsiVsPt1030Run2VsEpRun3 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb");
    latexTitle -> DrawLatex(0.18, 0.80, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4, 10#minus30\%");

    // SP
    TCanvas *canvasV2JpsiSpVsPt1030Run2VsRun3 = new TCanvas("canvasV2JpsiSpVsPt1030Run2VsRun3", "", 800, 600);
    canvasV2JpsiSpVsPt1030Run2VsRun3 -> SetTicks(1, 1);

    histGridV2JpsiSpVsPt1030 -> Draw();
    lineZeroV2Pt -> Draw();
    graStatV2JpsiVsPt1030Run2 -> Draw("EP SAME");
    graSystV2JpsiVsPt1030Run2 -> Draw("E2 SAME");
    graStatV2JpsiSpVsPt1030Run3 -> Draw("EP SAME");
    graSystV2JpsiSpVsPt1030Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsPt1030Run2VsSpRun3 = new TLegend(0.20, 0.70, 0.85, 0.85, " ", "brNDC");
    SetLegend(legendV2JpsiVsPt1030Run2VsSpRun3);
    legendV2JpsiVsPt1030Run2VsSpRun3 -> SetNColumns(2);
    legendV2JpsiVsPt1030Run2VsSpRun3 -> AddEntry(graSystV2JpsiVsPt1030Run2, "#sqrt{#it{s}_{NN}} = 5.02 TeV", "PF");
    legendV2JpsiVsPt1030Run2VsSpRun3 -> AddEntry(graSystV2JpsiSpVsPt1030Run3, "#sqrt{#it{s}_{NN}} = 5.36 TeV", "PF"); 
    legendV2JpsiVsPt1030Run2VsSpRun3 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb");
    latexTitle -> DrawLatex(0.18, 0.80, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4, 10#minus30\%");

    // 30-50% 
    // EP
    TCanvas *canvasV2JpsiEpVsPt3050Run2VsRun3 = new TCanvas("canvasV2JpsiEpVsPt3050Run2VsRun3", "", 800, 600);
    canvasV2JpsiEpVsPt3050Run2VsRun3 -> SetTicks(1, 1);

    histGridV2JpsiEpVsPt3050 -> Draw();
    lineZeroV2Pt -> Draw();
    graStatV2JpsiVsPt3050Run2 -> Draw("EP SAME");
    graSystV2JpsiVsPt3050Run2 -> Draw("E2 SAME");
    graStatV2JpsiEpVsPt3050Run3 -> Draw("EP SAME");
    graSystV2JpsiEpVsPt3050Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsPt3050Run2VsEpRun3 = new TLegend(0.20, 0.70, 0.85, 0.85, " ", "brNDC");
    SetLegend(legendV2JpsiVsPt3050Run2VsEpRun3);
    legendV2JpsiVsPt3050Run2VsEpRun3 -> SetNColumns(2);
    legendV2JpsiVsPt3050Run2VsEpRun3 -> AddEntry(graSystV2JpsiVsPt3050Run2, "#sqrt{#it{s}_{NN}} = 5.02 TeV", "PF");
    legendV2JpsiVsPt3050Run2VsEpRun3 -> AddEntry(graSystV2JpsiEpVsPt3050Run3, "#sqrt{#it{s}_{NN}} = 5.36 TeV", "PF"); 
    legendV2JpsiVsPt3050Run2VsEpRun3 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb");
    latexTitle -> DrawLatex(0.18, 0.80, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4, 30#minus50\%");

    // SP
    TCanvas *canvasV2JpsiSpVsPt3050Run2VsRun3 = new TCanvas("canvasV2JpsiSpVsPt3050Run2VsRun3", "", 800, 600);
    canvasV2JpsiSpVsPt3050Run2VsRun3 -> SetTicks(1, 1);

    histGridV2JpsiSpVsPt3050 -> Draw();
    lineZeroV2Pt -> Draw();
    graStatV2JpsiVsPt3050Run2 -> Draw("EP SAME");
    graSystV2JpsiVsPt3050Run2 -> Draw("E2 SAME");
    graStatV2JpsiSpVsPt3050Run3 -> Draw("EP SAME");
    graSystV2JpsiSpVsPt3050Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsPt3050Run2VsSpRun3 = new TLegend(0.20, 0.70, 0.85, 0.85, " ", "brNDC");
    SetLegend(legendV2JpsiVsPt3050Run2VsSpRun3);
    legendV2JpsiVsPt3050Run2VsSpRun3 -> SetNColumns(2);
    legendV2JpsiVsPt3050Run2VsSpRun3 -> AddEntry(graSystV2JpsiVsPt3050Run2, "#sqrt{#it{s}_{NN}} = 5.02 TeV", "PF");
    legendV2JpsiVsPt3050Run2VsSpRun3 -> AddEntry(graSystV2JpsiSpVsPt3050Run3, "#sqrt{#it{s}_{NN}} = 5.36 TeV", "PF"); 
    legendV2JpsiVsPt3050Run2VsSpRun3 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb");
    latexTitle -> DrawLatex(0.18, 0.80, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4, 30#minus50\%");
    
    //******************************************************************************************//
    // J/psi v2 vs pT collection plot
    //******************************************************************************************//
    // EP
    TCanvas *canvasV2JpsiEpVsPtLargeBinsCollPlot = new TCanvas("canvasV2JpsiEpVsPtLargeBinsCollPlot", "", 800, 600);
    canvasV2JpsiEpVsPtLargeBinsCollPlot -> SetTicks(1, 1);

    TH2D *histGridV2JpsiEpVsPtCollPlot = new TH2D("histGridV2JpsiEpVsPtCollPlot", "", 100, -0.5, 15.5, 100, -0.07, 0.3);
    histGridV2JpsiEpVsPtCollPlot -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histGridV2JpsiEpVsPtCollPlot -> GetYaxis() -> SetTitle("#it{v}_{2} {EP, |#Delta#eta| > 1.7}");
    histGridV2JpsiEpVsPtCollPlot -> Draw();

    lineZeroV2Pt -> Draw();
    graStatV2JpsiEpVsPt1030Run3LargeBinsCollPlot -> Draw("EP SAME");
    graSystV2JpsiEpVsPt1030Run3LargeBinsCollPlot -> Draw("E2 SAME");
    graStatV2JpsiEpVsPt3050Run3LargeBinsCollPlot -> Draw("EP SAME");
    graSystV2JpsiEpVsPt3050Run3LargeBinsCollPlot -> Draw("E2 SAME");
    graStatV2JpsiEpVsPt5080Run3CollPlot -> Draw("EP SAME");
    graSystV2JpsiEpVsPt5080Run3CollPlot -> Draw("E2 SAME");

    TLegend *legendV2JpsiEpVsPtColl = new TLegend(0.20, 0.70, 0.85, 0.85, " ", "brNDC");
    SetLegend(legendV2JpsiEpVsPtColl);
    legendV2JpsiEpVsPtColl -> SetNColumns(3);
    legendV2JpsiEpVsPtColl -> AddEntry(graSystV2JpsiEpVsPt1030Run3LargeBinsCollPlot, "10#minus30\%", "PF");
    legendV2JpsiEpVsPtColl -> AddEntry(graSystV2JpsiEpVsPt3050Run3LargeBinsCollPlot, "30#minus50\%", "PF");
    legendV2JpsiEpVsPtColl -> AddEntry(graSystV2JpsiEpVsPt5080Run3CollPlot, "50#minus80\%", "PF");
    legendV2JpsiEpVsPtColl -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb,#kern[0.3]{#sqrt{#it{s}_{NN}}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.18, 0.80, "J/#psi#rightarrow#mu^{#plus}#mu^{#minus}, 2.5 <#kern[0.7]{#it{y}} < 4");

    // SP
    TCanvas *canvasV2JpsiSpVsPtCollPlot = new TCanvas("canvasV2JpsiSpVsPtCollPlot", "", 800, 600);
    canvasV2JpsiSpVsPtCollPlot -> SetTicks(1, 1);

    TH2D *histGridV2JpsiSpVsPtCollPlot = new TH2D("histGridV2JpsiSpVsPtCollPlot", "", 100, -0.5, 15.5, 100, -0.07, 0.3);
    histGridV2JpsiSpVsPtCollPlot -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histGridV2JpsiSpVsPtCollPlot -> GetYaxis() -> SetTitle("#it{v}_{2} {SP, |#Delta#eta| > 1.7}");
    histGridV2JpsiSpVsPtCollPlot -> Draw();

    lineZeroV2Pt -> Draw();
    graStatV2JpsiSpVsPt010Run3CollPlot -> Draw("EP SAME");
    graSystV2JpsiSpVsPt010Run3CollPlot -> Draw("E2 SAME");
    graStatV2JpsiSpVsPt1030Run3CollPlot -> Draw("EP SAME");
    graSystV2JpsiSpVsPt1030Run3CollPlot -> Draw("E2 SAME");
    graStatV2JpsiSpVsPt3050Run3CollPlot -> Draw("EP SAME");
    graSystV2JpsiSpVsPt3050Run3CollPlot -> Draw("E2 SAME");

    TLegend *legendV2JpsiSpVsPtColl = new TLegend(0.20, 0.70, 0.85, 0.85, " ", "brNDC");
    SetLegend(legendV2JpsiSpVsPtColl);
    legendV2JpsiSpVsPtColl -> SetNColumns(3);
    legendV2JpsiSpVsPtColl -> AddEntry(graSystV2JpsiSpVsPt010Run3CollPlot, "0#minus10\%", "PF");
    legendV2JpsiSpVsPtColl -> AddEntry(graSystV2JpsiSpVsPt1030Run3CollPlot, "10#minus30\%", "PF");
    legendV2JpsiSpVsPtColl -> AddEntry(graSystV2JpsiSpVsPt3050Run3CollPlot, "30#minus50\%", "PF");
    legendV2JpsiSpVsPtColl -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb,#kern[0.3]{#sqrt{#it{s}_{NN}}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.18, 0.80, "J/#psi#rightarrow#mu^{#plus}#mu^{#minus}, 2.5 <#kern[0.7]{#it{y}} < 4");

    //******************************************************************************************//
    // J/psi v2 vs centrality vs Run 2
    //******************************************************************************************//
    // 0 < pT < 5 GeV/c
    // EP
    TCanvas *canvasV2JpsiEpVsCentr05Run2VsRun3 = new TCanvas("canvasV2JpsiEpVsCentr05Run2VsRun3", "", 800, 600);
    canvasV2JpsiEpVsCentr05Run2VsRun3 -> SetTicks(1, 1);

    histGridV2JpsiEpVsCentr05 -> Draw();
    lineZeroV2Centr -> Draw();
    graStatV2JpsiVsCentr05Run2 -> Draw("EP SAME");
    graSystV2JpsiVsCentr05Run2 -> Draw("E2 SAME");
    graStatV2JpsiEpVsCentr05Run3 -> Draw("EP SAME");
    graSystV2JpsiEpVsCentr05Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsCentr05Run2VsRun3 = new TLegend(0.20, 0.70, 0.85, 0.85, " ", "brNDC");
    SetLegend(legendV2JpsiVsCentr05Run2VsRun3);
    legendV2JpsiVsCentr05Run2VsRun3 -> SetNColumns(2);
    legendV2JpsiVsCentr05Run2VsRun3 -> AddEntry(graSystV2JpsiVsCentr05Run2, "#sqrt{#it{s}_{NN}} = 5.02 TeV", "PF");
    legendV2JpsiVsCentr05Run2VsRun3 -> AddEntry(graSystV2JpsiEpVsCentr05Run3, "#sqrt{#it{s}_{NN}} = 5.36 TeV", "PF"); 
    legendV2JpsiVsCentr05Run2VsRun3 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb");
    latexTitle -> DrawLatex(0.18, 0.80, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4, #it{p}_{T} < 5 GeV/#it{c}");

    // SP
    TCanvas *canvasV2JpsiSpVsCentr05Run2VsRun3 = new TCanvas("canvasV2JpsiSpVsCentr05Run2VsRun3", "", 800, 600);
    canvasV2JpsiSpVsCentr05Run2VsRun3 -> SetTicks(1, 1);

    histGridV2JpsiSpVsCentr05 -> Draw();
    lineZeroV2Centr -> Draw();
    graStatV2JpsiVsCentr05Run2 -> Draw("EP SAME");
    graSystV2JpsiVsCentr05Run2 -> Draw("E2 SAME");
    graStatV2JpsiSpVsCentr05Run3 -> Draw("EP SAME");
    graSystV2JpsiSpVsCentr05Run3 -> Draw("E2 SAME");
    legendV2JpsiVsCentr05Run2VsRun3 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb");
    latexTitle -> DrawLatex(0.18, 0.80, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4, #it{p}_{T} < 5 GeV/#it{c}");


    // 5 < pT < 15 GeV/c, EP
    TCanvas *canvasV2JpsiEpVsCentr515Run2VsRun3 = new TCanvas("canvasV2JpsiEpVsCentr515Run2VsRun3", "", 800, 600);
    canvasV2JpsiEpVsCentr515Run2VsRun3 -> SetTicks(1, 1);

    histGridV2JpsiEpVsCentr515 -> Draw();
    lineZeroV2Centr -> Draw();
    graStatV2JpsiVsCentr515Run2 -> Draw("EP SAME");
    graSystV2JpsiVsCentr515Run2 -> Draw("E2 SAME");
    graStatV2JpsiEpVsCentr515Run3 -> Draw("EP SAME");
    graSystV2JpsiEpVsCentr515Run3 -> Draw("E2 SAME");

    TLegend *legendV2JpsiVsCentr515Run2VsRun3 = new TLegend(0.20, 0.70, 0.85, 0.85, " ", "brNDC");
    SetLegend(legendV2JpsiVsCentr515Run2VsRun3);
    legendV2JpsiVsCentr515Run2VsRun3 -> SetNColumns(2);
    legendV2JpsiVsCentr515Run2VsRun3 -> AddEntry(graSystV2JpsiVsCentr515Run2, "#sqrt{#it{s}_{NN}} = 5.02 TeV", "PF");
    legendV2JpsiVsCentr515Run2VsRun3 -> AddEntry(graSystV2JpsiEpVsCentr515Run3, "#sqrt{#it{s}_{NN}} = 5.36 TeV", "PF"); 
    legendV2JpsiVsCentr515Run2VsRun3 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb");
    latexTitle -> DrawLatex(0.18, 0.80, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4, 5 < #it{p}_{T} < 15 GeV/#it{c}");

    // 5 < pT < 20 GeV/c, SP
    TCanvas *canvasV2JpsiSpVsCentr520Run2VsRun3 = new TCanvas("canvasV2JpsiSpVsCentr520Run2VsRun3", "", 800, 600);
    canvasV2JpsiSpVsCentr520Run2VsRun3 -> SetTicks(1, 1);

    histGridV2JpsiSpVsCentr520 -> Draw();
    lineZeroV2Centr -> Draw();
    graStatV2JpsiVsCentr515Run2 -> Draw("EP SAME");
    graSystV2JpsiVsCentr515Run2 -> Draw("E2 SAME");
    graStatV2JpsiSpVsCentr520Run3 -> Draw("EP SAME");
    graSystV2JpsiSpVsCentr520Run3 -> Draw("E2 SAME");
    legendV2JpsiVsCentr515Run2VsRun3 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb");
    latexTitle -> DrawLatex(0.18, 0.80, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4, 5 < #it{p}_{T} < 20 GeV/#it{c}");

    TCanvas *canvasV2JpsiVsCentr515Run2VsRun3 = new TCanvas("canvasV2JpsiVsCentr515Run2VsRun3", "", 800, 600);
    canvasV2JpsiVsCentr515Run2VsRun3 -> SetTicks(1, 1);

    TH2D *histGridV2JpsiVsCentr515Run2VsRun3 = new TH2D("histGridV2JpsiVsCentr515Run2VsRun3", "", 100, -0.5, 90.5, 100, -0.07, 0.3);
    histGridV2JpsiVsCentr515Run2VsRun3 -> GetXaxis() -> SetTitle("Centrality (%)");
    histGridV2JpsiVsCentr515Run2VsRun3 -> GetYaxis() -> SetTitle("#it{v}_{2} {EP, |#Delta#eta| > 1.7}");
    histGridV2JpsiVsCentr515Run2VsRun3 -> Draw();
    lineZeroV2Centr -> Draw();
    graStatV2JpsiVsCentr515Run2 -> Draw("EP SAME");
    graSystV2JpsiVsCentr515Run2 -> Draw("E2 SAME");
    graStatV2JpsiEpVsCentr515Run3 -> Draw("EP SAME");
    graSystV2JpsiEpVsCentr515Run3 -> Draw("E2 SAME");

    legendV2JpsiVsCentr515Run2VsRun3 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb");
    latexTitle -> DrawLatex(0.18, 0.80, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 <#kern[0.7]{#it{y}} < 4, 5 < #it{p}_{T} < 15 GeV/#it{c}");

    // EP
    canvasV2JpsiEpVsPt1030 -> SaveAs("QM2025_preliminary/v2JpsiEpVsPt1030.pdf");
    canvasV2JpsiEpVsPt3050 -> SaveAs("QM2025_preliminary/v2JpsiEpVsPt3050.pdf");
    canvasV2JpsiEpVsPt5080 -> SaveAs("QM2025_preliminary/v2JpsiEpVsPt5080.pdf");
    canvasV2JpsiEpVsPt1030Run2VsRun3 -> SaveAs("QM2025_preliminary/v2JpsiEpVsPt1030Run2VsRun3.pdf");
    canvasV2JpsiEpVsPt3050Run2VsRun3 -> SaveAs("QM2025_preliminary/v2JpsiEpVsPt3050Run2VsRun3.pdf");
    canvasV2JpsiEpVsPtLargeBinsCollPlot -> SaveAs("QM2025_preliminary/v2JpsiEpVsPtLargeBinsCollPlot.pdf");

    canvasV2JpsiEpVsCentr05 -> SaveAs("QM2025_preliminary/v2JpsiEpVsCentr05.pdf");
    canvasV2JpsiEpVsCentr515 -> SaveAs("QM2025_preliminary/v2JpsiEpVsCentr515.pdf");
    canvasV2JpsiEpVsCentr05Run2VsRun3 -> SaveAs("QM2025_preliminary/v2JpsiEpVsCentr05Run2VsRun3.pdf");
    canvasV2JpsiEpVsCentr515Run2VsRun3 -> SaveAs("QM2025_preliminary/v2JpsiEpVsCentr515Run2VsRun3.pdf");

    // SP
    canvasV2JpsiSpVsPt010 -> SaveAs("QM2025_preliminary/v2JpsiSpVsPt010.pdf");
    canvasV2JpsiSpVsPt1030 -> SaveAs("QM2025_preliminary/v2JpsiSpVsPt1030.pdf");
    canvasV2JpsiSpVsPt3050 -> SaveAs("QM2025_preliminary/v2JpsiSpVsPt3050.pdf");
    canvasV2JpsiSpVsPt010Run2VsRun3 -> SaveAs("QM2025_preliminary/v2JpsiSpVsPt010Run2VsRun3.pdf");
    canvasV2JpsiSpVsPt1030Run2VsRun3 -> SaveAs("QM2025_preliminary/v2JpsiSpVsPt1030Run2VsRun3.pdf");
    canvasV2JpsiSpVsPt3050Run2VsRun3 -> SaveAs("QM2025_preliminary/v2JpsiSpVsPt3050Run2VsRun3.pdf");
    canvasV2JpsiSpVsPtCollPlot -> SaveAs("QM2025_preliminary/v2JpsiSpVsPtCollPlot.pdf");

    canvasV2JpsiSpVsCentr05 -> SaveAs("QM2025_preliminary/v2JpsiSpVsCentr05.pdf");
    canvasV2JpsiSpVsCentr520 -> SaveAs("QM2025_preliminary/v2JpsiSpVsCentr520.pdf");
    canvasV2JpsiSpVsCentr05Run2VsRun3 -> SaveAs("QM2025_preliminary/v2JpsiSpVsCentr05Run2VsRun3.pdf");
    canvasV2JpsiSpVsCentr520Run2VsRun3 -> SaveAs("QM2025_preliminary/v2JpsiSpVsCentr520Run2VsRun3.pdf");
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