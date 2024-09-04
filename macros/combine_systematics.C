void combine_systematics() {
    const int nTrials = 72;

    const int nPtBins = 10;
    double minPtBins[] = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12};
    double maxPtBins[] = {1, 2, 3, 4, 5, 6, 8, 10, 12, 15};
    double ptBins[] = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15};
    //double minPtBins[] = {2};
    //double maxPtBins[] = {3};

    TH1D *histStatJpsiV2 = new TH1D("histStatJpsiV2", "", 10, ptBins);
    histStatJpsiV2 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histStatJpsiV2 -> GetYaxis() -> SetRangeUser(-0.05, 0.4);
    histStatJpsiV2 -> GetYaxis() -> SetTitle("#it{v}_{2}");
    histStatJpsiV2 -> SetMarkerStyle(20);
    histStatJpsiV2 -> SetMarkerColor(kRed+1);
    histStatJpsiV2 -> SetLineColor(kRed+1);

    TH1D *histSystJpsiV2 = new TH1D("histSystJpsiV2", "", 10, ptBins);
    histSystJpsiV2 -> SetMarkerStyle(20);
    histSystJpsiV2 -> SetMarkerColor(kRed+1);
    histSystJpsiV2 -> SetLineColor(kRed+1);
    histSystJpsiV2 -> SetFillStyle(0);

    for (int iPt = 0;iPt < nPtBins;iPt++) {
        TFile *fInStdFit = new TFile(Form("/Users/lucamicheletti/GITHUB/jpsi_flow_run3/macros/systematics_std_fit/fitResults_Pt_%2.1f_%2.1f.root", minPtBins[iPt], maxPtBins[iPt]));
        TFile *fInMixFit = new TFile(Form("/Users/lucamicheletti/GITHUB/jpsi_flow_run3/macros/systematics_mix_fit/fitResults_Pt_%2.1f_%2.1f.root", minPtBins[iPt], maxPtBins[iPt]));

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

        int nTrueTrials = 0;
        for (int iTrial = 0;iTrial < nTrials;iTrial++) {
            if (iTrial < nTrials/2) {
                double chi2Ndf = histChi2NdfStdFit -> GetBinContent(iTrial+1);
                if (chi2Ndf > 4 || chi2Ndf == 0) {
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
                if (chi2Ndf > 4 || chi2Ndf == 0) {
                    continue;
                }
                nTrueTrials++;
                jpsiV2s.push_back(histMixFit -> GetBinContent(iTrial+1-(nTrials/2)));
                errJpsiV2s.push_back(histMixFit -> GetBinError(iTrial+1-(nTrials/2)));
                string trialName = histMixFit -> GetXaxis() -> GetBinLabel(iTrial+1-(nTrials/2));

                histCombinedSysts -> SetBinContent(iTrial+1, histMixFit -> GetBinContent(iTrial+1-(nTrials/2)));
                histCombinedSysts -> SetBinError(iTrial+1, histMixFit -> GetBinError(iTrial+1-(nTrials/2)));
                histCombinedSysts -> GetXaxis() -> SetBinLabel(iTrial+1, Form("#color[864]{%s}", trialName.c_str()));
            }
        }

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

        histStatJpsiV2 -> SetBinContent(iPt+1, meanJpsiV2);
        histStatJpsiV2 -> SetBinError(iPt+1, meanStatErrJpsiV2);
        histSystJpsiV2 -> SetBinContent(iPt+1, meanJpsiV2);
        histSystJpsiV2 -> SetBinError(iPt+1, meanSystErrJpsiV2);

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
        if (meanJpsiV2 > 0.01) {
            histCombinedSysts -> GetYaxis() -> SetRangeUser(meanJpsiV2 -  (0.75*meanJpsiV2), meanJpsiV2 + (0.75*meanJpsiV2));
        } else {
            histCombinedSysts -> GetYaxis() -> SetRangeUser(-0.05, 0.05);
        } 
        if (minPtBins[iPt] >= 10) {
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

        TText *text1 = display1 -> AddText(Form(" v_{2}^{J/#psi} [%1.0f < #it{p}_{T} < %1.0f GeV/#it{c}] = %5.4f #pm %5.4f #pm %5.4f", minPtBins[iPt], maxPtBins[iPt], meanJpsiV2, meanStatErrJpsiV2, meanSystErrJpsiV2));
        display1 -> Draw("same");

        canvasCombinedSysts -> SaveAs(Form("combined_systematics/v2_sys_pt_%1.0f_%1.0f_cent_10_50.pdf", minPtBins[iPt], maxPtBins[iPt]));

        delete canvasCombinedSysts;
    }

    TCanvas *canvasJpsiV2 = new TCanvas("canvasJpsiV2", "", 800, 600);
    histStatJpsiV2 -> Draw("EP");
    histSystJpsiV2 -> Draw("E2P SAME");

    cout << "-------------------------" << endl;
    cout << "x_min x_max val stat syst " << endl;
    for (int iPt = 0;iPt < nPtBins;iPt++) {
        Printf("%3.2f %3.2f %6.5f %6.5f %6.5f ", minPtBins[iPt], maxPtBins[iPt], histStatJpsiV2 -> GetBinContent(iPt+1), histStatJpsiV2 -> GetBinError(iPt+1), histSystJpsiV2 -> GetBinError(iPt+1));
    }

}