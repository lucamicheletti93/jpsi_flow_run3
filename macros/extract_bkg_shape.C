double FuncBkg1(double *, double *);
double FuncBkg2(double *, double *);
double FuncBkg3(double *, double *);
double FuncBkg4(double *, double *);
double FuncPol2(double *, double *);
double FuncGaus(double *, double *);
double FuncVWG(double *, double *);
double FuncWeight1(double *, double *);
double FuncWeight2(double *, double *);

void extract_bkg_shape() {
    // Fixed variable -> Centrality
    double minFixVarBin = 10;
    double maxFixVarBin = 30;

    // Variable -> pT
    const int nVarBins = 13;
    double minVarBins[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0};
    double maxVarBins[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0};
    double varBins[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0};

    TFile *fInRatio = new TFile("/Users/lucamicheletti/Downloads/Histograms_Fullpass4PbPbQualitymatchedMchMid_Ratio_Mass_10_80.root", "READ");
    TFile *fInCut00 = new TFile("/Users/lucamicheletti/Downloads/Histograms_Fullpass4PbPbQualitymatchedMchMid_MchMid_CentBins_withMixing_10_80.root", "READ"); // MCH-MID
    TFile *fInCut07 = new TFile("/Users/lucamicheletti/Downloads/Histograms_Fullpass4PbPbQualitymatchedMchMid_withpt07_CentBins_withMixing_10_80.root", "READ"); // muon pT > 0.7 GeV/c


    TCanvas *canvasEventMixingBkg = new TCanvas("canvasEventMixingBkg", "", 1200, 1200);
    canvasEventMixingBkg -> Divide(5, 3);

    TCanvas *canvasRatio = new TCanvas("canvasRatio", "", 1200, 1200);
    canvasRatio -> Divide(5, 3);

    for (int iVar = 0;iVar < nVarBins;iVar++) {
        TH1D *histRatioMass = (TH1D*) fInRatio -> Get(Form("histMassMEPM_%2.1f_%2.1f__%1.0f_%1.0f", minVarBins[iVar], maxVarBins[iVar], minFixVarBin, maxFixVarBin));
        histRatioMass -> GetXaxis() -> SetRangeUser(2, 5);
        histRatioMass -> SetLineColor(kBlack);
        histRatioMass -> SetLineColor(kBlack);

        TH1D *histMassProj = (TH1D*) fInCut00 -> Get(Form("histMassMEPM_%2.1f_%2.1f__%1.0f_%1.0f", minVarBins[iVar], maxVarBins[iVar], minFixVarBin, maxFixVarBin));
        histMassProj -> GetXaxis() -> SetRangeUser(2, 5);
        histMassProj -> SetLineColor(kBlue);
        histMassProj -> SetLineColor(kBlue);

        TH1D *histMassProjPtCutRun3 = (TH1D*) fInCut07 -> Get(Form("histMassMEPM_%2.1f_%2.1f__%1.0f_%1.0f", minVarBins[iVar], maxVarBins[iVar], minFixVarBin, maxFixVarBin));
        histMassProjPtCutRun3 -> GetXaxis() -> SetRangeUser(2, 5);
        histMassProjPtCutRun3 -> SetLineColor(kRed);
        histMassProjPtCutRun3 -> SetLineColor(kRed);

        TF1 *funcWeight2 = new TF1("funcWeight2", FuncWeight2, 1.9, 4.5, 5);
        funcWeight2 -> SetParameter(0, 1);
        funcWeight2 -> SetParameter(1, 1);
        funcWeight2 -> SetParameter(2, 1);
        funcWeight2 -> SetParameter(3, 1);
        funcWeight2 -> SetParameter(4, 1);

        histRatioMass -> Fit(funcWeight2, "R0");
        funcWeight2 -> SetLineColor(kViolet);

        TF1 *funcVWG = new TF1("FuncVWG", FuncVWG, 1.9, 4.5, 4);
        funcVWG -> SetParameter(0, 1.6432);
        funcVWG -> SetParameter(1, 4.5482e-01);
        funcVWG -> SetParameter(2, 3.6494e-01);
        funcVWG -> SetParameter(3, 1e6);

        histMassProj -> Fit(funcVWG, "R0");
        funcVWG -> SetLineColor(kRed);

        TF1 *funcBkg4 = new TF1("funcBkg4", FuncBkg4, 1.9, 4.5, 9);
        funcBkg4 -> FixParameter(0, funcWeight2 -> GetParameter(0));
        funcBkg4 -> FixParameter(1, funcWeight2 -> GetParameter(1));
        funcBkg4 -> FixParameter(2, funcWeight2 -> GetParameter(2));
        funcBkg4 -> FixParameter(3, funcWeight2 -> GetParameter(3));
        funcBkg4 -> FixParameter(4, funcWeight2 -> GetParameter(4));
        funcBkg4 -> SetParameter(5, 1e6);
        funcBkg4 -> SetParameter(6, funcVWG -> GetParameter(0));
        funcBkg4 -> SetParameter(7, funcVWG -> GetParameter(1));
        funcBkg4 -> SetParameter(8, funcVWG -> GetParameter(2));

        histMassProjPtCutRun3 -> Fit(funcBkg4, "R0");
        funcBkg4 -> SetLineColor(kViolet);


        canvasRatio -> cd(iVar+1);
        gPad -> SetLogy(true);
        histRatioMass -> Draw("EP");
        funcWeight2 -> Draw("SAME");

        canvasEventMixingBkg -> cd(iVar+1);
        gPad -> SetLogy(true);
        histMassProj -> Draw("EP");
        histMassProjPtCutRun3 -> Draw("EP SAME");
        funcVWG -> Draw("SAME");
        funcBkg4 -> Draw("SAME");
    }
}
///////////////////////////////////////////////////
double FuncBkg1(double *x, double *par) {
    double xx = x[0];
    double alpha = par[0];

    if (xx > alpha - 0.2 && xx < alpha + 0.2) {
        return par[1] + par[2]*xx + par[3]*xx*xx;
    } else {
        if (xx > alpha + 0.2) {
            return par[4] + par[5]*xx + par[6]*xx*xx;
        }
        if (xx < alpha - 0.2) {
            return par[7] + par[8]*xx + par[9]*xx*xx;
        }
    }
    return 0;
}
///////////////////////////////////////////////////
double FuncBkg2(double *x, double *par) {
    double xx = x[0];
    double mean = par[0];
    double sigma = par[1];

    if (xx > (mean - sigma) && xx < (mean + sigma)) {
        return par[2] + par[3]*xx + par[4]*xx*xx + par[5] + par[6]*xx + par[7]*xx*xx;
    }
    if (xx >= (mean + sigma)) {
        return par[2] + par[3]*xx + par[4]*xx*xx;
    }
    if (xx < (mean - sigma)) {
        return par[5] + par[6]*xx + par[7]*xx*xx;
    }
    return 0;
}
///////////////////////////////////////////////////
double FuncBkg3(double *x, double *par) {
    double xx = x[0];

    double mean1 = par[1];
    double sigma1 = par[2];
    double arg1 = - ((xx - mean1) * (xx - mean1)) / (2 * sigma1 * sigma1);

    double mean2 = par[5];
    double sigma2 = par[6];
    double arg2 = - ((xx - mean2) * (xx - mean2)) / (2 * sigma2 * sigma2);

    return (par[4] * TMath::Exp(arg2)) / (par[3] + par[0] * TMath::Exp(arg1));
}
///////////////////////////////////////////////////
double FuncBkg4(double *x, double *par) {
    double xx = x[0];

    double mean1 = par[1];
    double sigma1 = par[2] + par[3] * ((xx - mean1) / mean1);
    double arg1 = - ((xx - mean1) * (xx - mean1)) / (2 * sigma1 * sigma1);

    double mean2 = par[6];
    double sigma2 = par[7] + par[8] * ((xx - mean2) / mean2);
    double arg2 = - (xx - mean2) * (xx - mean2) / (2. * sigma2 * sigma2);

    return (par[5] * TMath::Exp(arg2)) / (par[4] + par[0] * TMath::Exp(arg1));
}
///////////////////////////////////////////////////
double FuncPol2(double *x, double *par) {
    double xx = x[0];
    return par[0] + par[1]*xx + par[2]*xx*xx;
}
///////////////////////////////////////////////////
double FuncGaus(double *x, double *par) {
    double xx = x[0];
    double mean = par[0];
    double sigma = par[1];

    double arg = - ((xx - mean) * (xx - mean)) / (2 * sigma * sigma);
    return par[2] * TMath::Exp(arg);
}
///////////////////////////////////////////////////
double FuncVWG(double *x, double *par) {
    double xx = x[0];
    double mean = par[0];
    double sigma = par[1] + par[2] * ((xx - mean) / mean);
    double arg = - (xx - mean) * (xx - mean) / (2. * sigma * sigma);

    return par[3] * TMath::Exp(arg);
}
///////////////////////////////////////////////////
double FuncWeight1(double *x, double *par) {
    double xx = x[0];
    double mean = par[1];
    double sigma = par[2];

    double arg = - ((xx - mean) * (xx - mean)) / (2 * sigma * sigma);
    return par[3] + par[0] * TMath::Exp(arg);
}
///////////////////////////////////////////////////
double FuncWeight2(double *x, double *par) {
    double xx = x[0];
    double mean = par[1];
    double sigma = par[2] + par[3] * ((xx - mean) / mean);
    double arg = - (xx - mean) * (xx - mean) / (2. * sigma * sigma);

    return par[4] + par[0] * TMath::Exp(arg);
}
