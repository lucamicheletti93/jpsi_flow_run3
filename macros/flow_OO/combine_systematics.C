void combine_systematics(string fixVar = "centrality", double minFixVar = 0, double maxFixVar = 20, string binning = "narrow", string strEtaGap = "etaGap17") {
    string dirInPathStd = Form("/Users/lucamicheletti/GITHUB/jpsi_flow_run3/macros/flow_OO/systematics_std_fit/%s/%s_%1.0f_%1.0f", strEtaGap.c_str(), fixVar.c_str(), minFixVar, maxFixVar);
    string dirInPathMix = Form("/Users/lucamicheletti/GITHUB/jpsi_flow_run3/macros/flow_OO/systematics_mix_fit/%s/%s_%1.0f_%1.0f", strEtaGap.c_str(), fixVar.c_str(), minFixVar, maxFixVar);
    string dirOutPath = Form("combined_systematics/%s/%s_%1.0f_%1.0f", strEtaGap.c_str(), fixVar.c_str(), minFixVar, maxFixVar);

    const int nTrials = 48;
    string varAxisTitle, varName, varFixName;
    int nVarBins = -999;
    vector <double> vecMinVarBins, vecMaxVarBins, vecVarBins;

    //WARNING! 1% systematic applied to all centrality for pT integrated results
    std::map<std::pair<double, double>, double> mapResoRelErr = {{{0, 20}, 0.01}, {{20, 60}, 0.01}, {{10, 30}, 0.01}, {{30, 50}, 0.01}, {{50, 80}, 0.017}, {{0, 5}, 0.01}, {{5, 15}, 0.01}};

    if (fixVar == "centrality") {
        // Pt dependence
        varAxisTitle = "#it{p}_{T} (GeV/#it{c})";
        varName = "Pt";
        varFixName = "centrality";

        if (binning == "narrow") {
            nVarBins = 6;
            double minVarBins[] = {0.0, 1.0, 2.0, 3.0, 4.0, 6.0};
            double maxVarBins[] = {1.0, 2.0, 3.0, 4.0, 6.0, 8.0};
            double varBins[] = {0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0};

            for (int iPt = 0;iPt < nVarBins+1;iPt++) {
                if (iPt < nVarBins) {
                    vecMinVarBins.push_back(minVarBins[iPt]);
                    vecMaxVarBins.push_back(maxVarBins[iPt]);
                }
                vecVarBins.push_back(varBins[iPt]);
            }
        } else {
            nVarBins = 4;
            double minVarBins[] = {0.0, 2.0, 4.0, 6.0};
            double maxVarBins[] = {2.0, 4.0, 6.0, 8.0};
            double varBins[] = {0.0, 2.0, 4.0, 6.0, 8.0};

            for (int iPt = 0;iPt < nVarBins+1;iPt++) {
                if (iPt < nVarBins) {
                    vecMinVarBins.push_back(minVarBins[iPt]);
                    vecMaxVarBins.push_back(maxVarBins[iPt]);
                }
                vecVarBins.push_back(varBins[iPt]);
            }
        }
    } else {
        // Centrality dependence
        varAxisTitle = "Centrality (%)";
        varName = "Centr";
        varFixName = "pt";

        nVarBins = 8;
        double minVarBins[] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0};
        double maxVarBins[] = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0};
        double varBins[] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0};

        for (int iCentr = 0;iCentr < nVarBins+1;iCentr++) {
            if (iCentr < nVarBins) {
                vecMinVarBins.push_back(minVarBins[iCentr]);
                vecMaxVarBins.push_back(maxVarBins[iCentr]);
            }
            vecVarBins.push_back(varBins[iCentr]);
        }
    }

    TCanvas *canvasChi2Ndf = new TCanvas("canvasChi2Ndf", "", 3000, 1800);
    canvasChi2Ndf -> Divide(3, 2);

    TH1D *histStatJpsiV2 = new TH1D("histStatJpsiV2", "", nVarBins, &(vecVarBins[0]));
    histStatJpsiV2 -> GetXaxis() -> SetTitle(varAxisTitle.c_str());
    histStatJpsiV2 -> GetYaxis() -> SetRangeUser(-0.05, 0.4);
    histStatJpsiV2 -> GetYaxis() -> SetTitle("#it{v}_{2}");
    histStatJpsiV2 -> SetMarkerStyle(20);
    histStatJpsiV2 -> SetMarkerColor(kRed+1);
    histStatJpsiV2 -> SetLineColor(kRed+1);

    TH1D *histSystJpsiV2 = new TH1D("histSystJpsiV2", "", nVarBins, &(vecVarBins[0]));
    histSystJpsiV2 -> SetMarkerStyle(20);
    histSystJpsiV2 -> SetMarkerColor(kRed+1);
    histSystJpsiV2 -> SetLineColor(kRed+1);
    histSystJpsiV2 -> SetFillStyle(0);

    TH1D *histStatJpsiMixV2 = new TH1D("histStatJpsiMixV2", "", nVarBins, &(vecVarBins[0]));
    histStatJpsiMixV2 -> GetXaxis() -> SetTitle(varAxisTitle.c_str());
    histStatJpsiMixV2 -> GetYaxis() -> SetRangeUser(-0.05, 0.4);
    histStatJpsiMixV2 -> GetYaxis() -> SetTitle("#it{v}_{2}");
    histStatJpsiMixV2 -> SetMarkerStyle(20);
    histStatJpsiMixV2 -> SetMarkerColor(kBlue+1);
    histStatJpsiMixV2 -> SetLineColor(kBlue+1);

    TH1D *histSystJpsiMixV2 = new TH1D("histSystJpsiMixV2", "", nVarBins, &(vecVarBins[0]));
    histSystJpsiMixV2 -> SetMarkerStyle(20);
    histSystJpsiMixV2 -> SetMarkerColor(kBlue+1);
    histSystJpsiMixV2 -> SetLineColor(kBlue+1);
    histSystJpsiMixV2 -> SetFillStyle(0);

    for (int iVar = 0;iVar < nVarBins;iVar++) {
        TFile *fInStdFit = new TFile(Form("%s/fitResults_%s_%2.1f_%2.1f.root", dirInPathStd.c_str(), varName.c_str(), vecMinVarBins[iVar], vecMaxVarBins[iVar]));
        TFile *fInMixFit = new TFile(Form("%s/fitResults_%s_%2.1f_%2.1f.root", dirInPathMix.c_str(), varName.c_str(), vecMinVarBins[iVar], vecMaxVarBins[iVar]));

        TCanvas *canvasStdFit = (TCanvas*) fInStdFit -> Get("canvasSyst");
        TH1D *histStdFit = (TH1D*) canvasStdFit -> GetPrimitive("histSyst");
        TH1D *histChi2NdfStdFit = (TH1D*) fInStdFit -> Get("histChi2Ndf");

        TCanvas *canvasMixFit = (TCanvas*) fInMixFit -> Get("canvasSyst");
        TH1D *histMixFit = (TH1D*) canvasMixFit -> GetPrimitive("histSyst");
        TH1D *histChi2NdfMixFit = (TH1D*) fInMixFit -> Get("histChi2Ndf");

        vector<double> jpsiV2s;
        vector<double> errJpsiV2s;

        TH1D *histCombinedSysts = new TH1D("histCombinedSysts", "", nTrials, 0, nTrials);
        histCombinedSysts -> SetMarkerStyle(20);
        histCombinedSysts -> SetMarkerColor(kBlack);
        histCombinedSysts -> SetLineColor(kBlack);

        vector<double> jpsiMixV2s;
        vector<double> errJpsiMixV2s;

        TH1D *histMixSysts = new TH1D("histMixSysts", "", nTrials/2, 0, nTrials/2);
        histMixSysts -> SetMarkerStyle(20);
        histMixSysts -> SetMarkerColor(kBlack);
        histMixSysts -> SetLineColor(kBlack);

        TH1D *histChi2NdfAll = new TH1D("histChi2NdfAll", Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c} ; #chi^{2} / NDF", vecMinVarBins[iVar], vecMaxVarBins[iVar]), 100, 0, 5);
        int nTrueTrials = 0;
        int nTrueTrialsMix = 0;
        for (int iTrial = 0;iTrial < nTrials;iTrial++) {
            if (iTrial < nTrials/2) {
                double chi2Ndf = histChi2NdfStdFit -> GetBinContent(iTrial+1);
                histChi2NdfAll -> Fill(chi2Ndf);
                double v2Jpsi = histStdFit -> GetBinContent(iTrial+1);
                if (chi2Ndf > 4 || chi2Ndf == 0 || v2Jpsi < -0.15 || v2Jpsi > 0.3) {
                    continue;
                }
                nTrueTrials++;
                jpsiV2s.push_back(histStdFit -> GetBinContent(iTrial+1));
                errJpsiV2s.push_back(histStdFit -> GetBinError(iTrial+1));

                histCombinedSysts -> SetBinContent(iTrial+1, histStdFit -> GetBinContent(iTrial+1));
                histCombinedSysts -> SetBinError(iTrial+1, histStdFit -> GetBinError(iTrial+1));
                histCombinedSysts -> GetXaxis() -> SetBinLabel(iTrial+1, histStdFit -> GetXaxis() -> GetBinLabel(iTrial+1));
            } else {
                double chi2Ndf = histChi2NdfMixFit -> GetBinContent(iTrial+1-(nTrials/2));
                histChi2NdfAll -> Fill(chi2Ndf);
                double v2Jpsi = histMixFit -> GetBinContent(iTrial+1-(nTrials/2));
                if (chi2Ndf > 4 || chi2Ndf == 0 || v2Jpsi < -0.15 || v2Jpsi > 0.3) {
                    continue;
                }
                nTrueTrials++;
                jpsiV2s.push_back(histMixFit -> GetBinContent(iTrial+1-(nTrials/2)));
                errJpsiV2s.push_back(histMixFit -> GetBinError(iTrial+1-(nTrials/2)));
                string trialName = histMixFit -> GetXaxis() -> GetBinLabel(iTrial+1-(nTrials/2));

                histCombinedSysts -> SetBinContent(iTrial+1, histMixFit -> GetBinContent(iTrial+1-(nTrials/2)));
                histCombinedSysts -> SetBinError(iTrial+1, histMixFit -> GetBinError(iTrial+1-(nTrials/2)));
                histCombinedSysts -> GetXaxis() -> SetBinLabel(iTrial+1, Form("#color[864]{%s}", trialName.c_str()));

                nTrueTrialsMix++;
                jpsiMixV2s.push_back(histMixFit -> GetBinContent(iTrial+1-(nTrials/2)));
                errJpsiMixV2s.push_back(histMixFit -> GetBinError(iTrial+1-(nTrials/2)));

                histMixSysts -> SetBinContent(iTrial+1-(nTrials/2), histMixFit -> GetBinContent(iTrial+1-(nTrials/2)));
                histMixSysts -> SetBinError(iTrial+1-(nTrials/2), histMixFit -> GetBinError(iTrial+1-(nTrials/2)));
                histMixSysts -> GetXaxis() -> SetBinLabel(iTrial+1-(nTrials/2), Form("%s", trialName.c_str()));
            }
        }

        // combined std and mix fits
        // compute the mean v2, the mean v2 error and the systematic
        double sumJpsiV2s = 0; 
        double statErrJpsiV2s = 0;

        for (int iTrial = 0;iTrial < nTrueTrials;iTrial++) {
            sumJpsiV2s = sumJpsiV2s + jpsiV2s[iTrial];
            statErrJpsiV2s = statErrJpsiV2s + errJpsiV2s[iTrial];
        }
        double meanJpsiV2 = sumJpsiV2s / nTrueTrials;
        double meanStatErrJpsiV2 = statErrJpsiV2s / nTrueTrials;

        double sumSyst = 0;
        for (int iTrial = 0;iTrial < nTrueTrials;iTrial++) {
        double dev = (jpsiV2s[iTrial] - meanJpsiV2);
            sumSyst = sumSyst + dev * dev;
        }
        double meanSystErrJpsiV2 = TMath::Sqrt(sumSyst / nTrueTrials);

        histStatJpsiV2 -> SetBinContent(iVar+1, meanJpsiV2);
        histStatJpsiV2 -> SetBinError(iVar+1, meanStatErrJpsiV2);
        histSystJpsiV2 -> SetBinContent(iVar+1, meanJpsiV2);
        histSystJpsiV2 -> SetBinError(iVar+1, meanSystErrJpsiV2);

        TLine *lineMean = new TLine(0, meanJpsiV2, nTrials, meanJpsiV2);
        lineMean -> SetLineStyle(1);
        lineMean -> SetLineColor(kBlue);
        lineMean -> SetLineWidth(2);

        TLine *lineSystUp = new TLine(0, meanJpsiV2 + meanSystErrJpsiV2, nTrials, meanJpsiV2 + meanSystErrJpsiV2);
        lineSystUp-> SetLineStyle(2);
        lineSystUp-> SetLineColor(kBlue);
        lineSystUp-> SetLineWidth(2);

        TLine *lineSystLw = new TLine(0, meanJpsiV2 - meanSystErrJpsiV2, nTrials, meanJpsiV2 - meanSystErrJpsiV2);
        lineSystLw -> SetLineStyle(2);
        lineSystLw -> SetLineColor(kBlue);
        lineSystLw -> SetLineWidth(2);

        TLine *lineStatUp = new TLine(0, meanJpsiV2 + meanStatErrJpsiV2, nTrials, meanJpsiV2 + meanStatErrJpsiV2);
        lineStatUp-> SetLineStyle(2);
        lineStatUp-> SetLineColor(kGray+1);
        lineStatUp-> SetLineWidth(2);

        TLine *lineStatLw = new TLine(0, meanJpsiV2 - meanStatErrJpsiV2, nTrials, meanJpsiV2 - meanStatErrJpsiV2);
        lineStatLw -> SetLineStyle(2);
        lineStatLw -> SetLineColor(kGray+1);
        lineStatLw -> SetLineWidth(2);

        TCanvas *canvasCombinedSysts = new TCanvas("canvasCombinedSysts", "", 1400, 900);
        canvasCombinedSysts -> SetTopMargin(0.05);
        canvasCombinedSysts -> SetBottomMargin(0.5);
        histCombinedSysts -> SetStats(0);
        histCombinedSysts -> GetXaxis() -> LabelsOption("v");
        histCombinedSysts -> GetYaxis() -> SetRangeUser(meanJpsiV2 -  0.12, meanJpsiV2 + 0.12);
        /*if (meanJpsiV2 > 0.01) {
            histCombinedSysts -> GetYaxis() -> SetRangeUser(meanJpsiV2 -  (2*meanJpsiV2), meanJpsiV2 + (2*meanJpsiV2));
        } else if (meanJpsiV2 < 0){ 
            histCombinedSysts -> GetYaxis() -> SetRangeUser(-0.15, 0.1);
        }else {
            histCombinedSysts -> GetYaxis() -> SetRangeUser(-0.05, 0.05);
        }*/
        
        if (vecMinVarBins[iVar] >= 10) {
            //histCombinedSysts -> GetYaxis() -> SetRangeUser(-0.05, 0.2);
            histCombinedSysts -> GetYaxis() -> SetRangeUser(meanJpsiV2 -  (2*meanJpsiV2), meanJpsiV2 + (2*meanJpsiV2));
        }

        histCombinedSysts -> Draw();
        lineMean -> Draw();
        lineSystUp-> Draw();
        lineSystLw -> Draw();
        lineStatUp-> Draw();
        lineStatLw -> Draw();

        TPaveText *display1 = new TPaveText(0.20, 0.88, 0.80, 0.90, "blNDC");
        display1 -> SetTextFont(42);
        display1 -> SetTextSize(0.036);
        display1 -> SetTextColor(kBlack);
        display1 -> SetBorderSize(0);
        display1 -> SetFillColor(0);

        TText *text1 = display1 -> AddText(Form(" v_{2}^{J/#psi} [%1.0f < #it{p}_{T} < %1.0f GeV/#it{c}] = %5.4f #pm %5.4f #pm %5.4f", vecMinVarBins[iVar], vecMaxVarBins[iVar], meanJpsiV2, meanStatErrJpsiV2, meanSystErrJpsiV2));
        display1 -> Draw("same");


        // combined mix fits only
        // compute the mean v2, the mean v2 error and the systematic
        double sumMixJpsiV2s = 0; 
        double statErrMixJpsiV2s = 0;

        for (int iTrial = 0;iTrial < nTrueTrialsMix;iTrial++) {
            sumMixJpsiV2s = sumMixJpsiV2s + jpsiMixV2s[iTrial];
            statErrMixJpsiV2s = statErrMixJpsiV2s + errJpsiMixV2s[iTrial];
        }
        double meanJpsiMixV2 = sumMixJpsiV2s / nTrueTrialsMix;
        double meanStatErrJpsiMixV2 = statErrMixJpsiV2s / nTrueTrialsMix;

        double sumMixSyst = 0;
        for (int iTrial = 0;iTrial < nTrueTrialsMix;iTrial++) {
        double dev = (jpsiMixV2s[iTrial] - meanJpsiV2);
            sumMixSyst = sumMixSyst + dev * dev;
        }
        double meanSystErrJpsiMixV2 = TMath::Sqrt(sumMixSyst / nTrueTrialsMix);

        histStatJpsiMixV2 -> SetBinContent(iVar+1, meanJpsiMixV2);
        histStatJpsiMixV2 -> SetBinError(iVar+1, meanStatErrJpsiMixV2);
        histSystJpsiMixV2 -> SetBinContent(iVar+1, meanJpsiMixV2);
        histSystJpsiMixV2 -> SetBinError(iVar+1, meanSystErrJpsiMixV2);

        TLine *lineMeanMix = new TLine(0, meanJpsiMixV2, nTrials/2, meanJpsiMixV2);
        lineMeanMix -> SetLineStyle(1);
        lineMeanMix -> SetLineColor(kBlue);
        lineMeanMix -> SetLineWidth(2);

        TLine *lineSystUpMix = new TLine(0, meanJpsiMixV2 + meanSystErrJpsiMixV2, nTrials/2, meanJpsiMixV2 + meanSystErrJpsiMixV2);
        lineSystUpMix-> SetLineStyle(2);
        lineSystUpMix-> SetLineColor(kBlue);
        lineSystUpMix-> SetLineWidth(2);

        TLine *lineSystLwMix = new TLine(0, meanJpsiMixV2 - meanSystErrJpsiMixV2, nTrials/2, meanJpsiMixV2 - meanSystErrJpsiMixV2);
        lineSystLwMix -> SetLineStyle(2);
        lineSystLwMix -> SetLineColor(kBlue);
        lineSystLwMix -> SetLineWidth(2);

        TLine *lineStatUpMix = new TLine(0, meanJpsiMixV2 + meanStatErrJpsiMixV2, nTrials/2, meanJpsiMixV2 + meanStatErrJpsiMixV2);
        lineStatUpMix-> SetLineStyle(2);
        lineStatUpMix-> SetLineColor(kGray+1);
        lineStatUpMix-> SetLineWidth(2);

        TLine *lineStatLwMix = new TLine(0, meanJpsiMixV2 - meanStatErrJpsiMixV2, nTrials/2, meanJpsiMixV2 - meanStatErrJpsiMixV2);
        lineStatLwMix -> SetLineStyle(2);
        lineStatLwMix -> SetLineColor(kGray+1);
        lineStatLwMix -> SetLineWidth(2);

        TCanvas *canvasMixSysts = new TCanvas("canvasMixSysts", "", 1400, 900);
        canvasMixSysts -> SetTopMargin(0.05);
        canvasMixSysts -> SetBottomMargin(0.5);
        histMixSysts -> SetStats(0);
        histMixSysts -> GetXaxis() -> LabelsOption("v");
        histMixSysts -> GetYaxis() -> SetRangeUser(meanJpsiV2 -  0.12, meanJpsiV2 + 0.12);
        histMixSysts -> Draw();
        canvasMixSysts -> Update();

        lineMeanMix -> Draw();
        lineSystUpMix-> Draw();
        lineSystLwMix -> Draw();
        lineStatUpMix -> Draw();
        lineStatLwMix -> Draw();

        TPaveText *displayMix1 = new TPaveText(0.20, 0.88, 0.80, 0.90, "blNDC");
        displayMix1 -> SetTextFont(42);
        displayMix1 -> SetTextSize(0.036);
        displayMix1 -> SetTextColor(kBlack);
        displayMix1 -> SetBorderSize(0);
        displayMix1 -> SetFillColor(0);

        TText *textMix1 = displayMix1 -> AddText(Form(" v_{2}^{J/#psi} [%1.0f < #it{p}_{T} < %1.0f GeV/#it{c}] = %5.4f #pm %5.4f #pm %5.4f", vecMinVarBins[iVar], vecMaxVarBins[iVar], meanJpsiMixV2, meanStatErrJpsiMixV2, meanSystErrJpsiMixV2));
        displayMix1 -> Draw("same");

        if (varName == "Centr") {
            canvasCombinedSysts -> SaveAs(Form("%s/v2_sys_cent_%1.0f_%1.0f_pt_%1.0f_%1.0f.pdf", dirOutPath.c_str(), vecMinVarBins[iVar], vecMaxVarBins[iVar], minFixVar, maxFixVar));
        }
        if (varName == "Pt") {
            canvasCombinedSysts -> SaveAs(Form("%s/v2_sys_pt_%1.0f_%1.0f_cent_%1.0f_%1.0f.pdf", dirOutPath.c_str(), vecMinVarBins[iVar], vecMaxVarBins[iVar], minFixVar, maxFixVar));
            canvasMixSysts -> SaveAs(Form("%s/v2_sys_pt_%1.0f_%1.0f_cent_%1.0f_%1.0f_mixing_only.pdf", dirOutPath.c_str(), vecMinVarBins[iVar], vecMaxVarBins[iVar], minFixVar, maxFixVar));
        }

        canvasChi2Ndf -> cd(iVar+1);
        histChi2NdfAll -> Draw("HIST");

        delete canvasCombinedSysts;
    }

    canvasChi2Ndf -> SaveAs(Form("%s/Chi2Summary.pdf", dirOutPath.c_str()));

    TCanvas *canvasJpsiV2 = new TCanvas("canvasJpsiV2", "", 800, 600);
    histStatJpsiV2 -> Draw("EP");
    histSystJpsiV2 -> Draw("E2P SAME");
    histStatJpsiMixV2 -> Draw("EP SAME");
    histSystJpsiMixV2 -> Draw("E2P SAME");

    cout << "----------------------------------------------------------------------------------" << endl;
    std::pair<int, int> key = {minFixVar, maxFixVar};
    double resoRelErr;
    auto tmpResoRelErr = mapResoRelErr.find(key);
    if (tmpResoRelErr != mapResoRelErr.end()) {
        resoRelErr = tmpResoRelErr -> second;
    }

    cout << Form("******** Adding in quadrature the syst. on resolution (%f) *******", resoRelErr) << endl;
    cout << "x_min x_max val stat syst " << endl;
    for (int iVar = 0;iVar < nVarBins;iVar++) {
        double systSigExtr = histSystJpsiV2 -> GetBinError(iVar+1);
        double systReso = histStatJpsiV2 -> GetBinContent(iVar+1) * resoRelErr;
        double systAll = TMath::Sqrt(systSigExtr*systSigExtr + systReso*systReso);
        //Printf("%3.2f %3.2f %6.5f %6.5f %6.5f ", vecMinVarBins[iVar], vecMaxVarBins[iVar], histStatJpsiV2 -> GetBinContent(iVar+1), histStatJpsiV2 -> GetBinError(iVar+1), histSystJpsiV2 -> GetBinError(iVar+1));
        Printf("%3.2f %3.2f %6.5f %6.5f %6.5f ", vecMinVarBins[iVar], vecMaxVarBins[iVar], histStatJpsiV2 -> GetBinContent(iVar+1), histStatJpsiV2 -> GetBinError(iVar+1), systAll);
    }
    cout << "----------------------------------------------------------------------------------" << endl;
    cout << "* * * Combining standard and mixing fits * * *" << endl;
    for (int iVar = 0;iVar < nVarBins;iVar++) {
        std::cout << histStatJpsiV2 -> GetBinContent(iVar+1) << ", ";
    }
    std::cout << std::endl;
    for (int iVar = 0;iVar < nVarBins;iVar++) {
        std::cout << histStatJpsiV2 -> GetBinError(iVar+1) << ", ";
    }
    std::cout << std::endl;
    for (int iVar = 0;iVar < nVarBins;iVar++) {
        double systSigExtr = histSystJpsiV2 -> GetBinError(iVar+1);
        double systReso = histStatJpsiV2 -> GetBinContent(iVar+1) * resoRelErr;
        double systAll = TMath::Sqrt(systSigExtr*systSigExtr + systReso*systReso);
        std::cout << systAll << ", ";
    }
    std::cout << std::endl;
    cout << "----------------------------------------------------------------------------------" << endl;

    cout << "----------------------------------------------------------------------------------" << endl;
    cout << "* * * Mixing fits only * * *" << endl;
    for (int iVar = 0;iVar < nVarBins;iVar++) {
        std::cout << histStatJpsiMixV2 -> GetBinContent(iVar+1) << ", ";
    }
    std::cout << std::endl;
    for (int iVar = 0;iVar < nVarBins;iVar++) {
        std::cout << histStatJpsiMixV2 -> GetBinError(iVar+1) << ", ";
    }
    std::cout << std::endl;
    for (int iVar = 0;iVar < nVarBins;iVar++) {
        double systSigExtr = histSystJpsiMixV2 -> GetBinError(iVar+1);
        double systReso = histStatJpsiMixV2 -> GetBinContent(iVar+1) * resoRelErr;
        double systAll = TMath::Sqrt(systSigExtr*systSigExtr + systReso*systReso);
        std::cout << systAll << ", ";
    }
    std::cout << std::endl;
    cout << "----------------------------------------------------------------------------------" << endl;
}