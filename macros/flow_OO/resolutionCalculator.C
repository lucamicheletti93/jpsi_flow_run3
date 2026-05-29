void resolutionCalculator() {
    TFile *fIn = TFile::Open("/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/LHC25ae/AnalysisResults_687739.root");

    THashList *hList = (THashList*) fIn->Get("analysis-event-selection/output");
    TList *subList = (TList*) hList->FindObject("Event_AfterCuts");
    TH2F *hs_R2SPAB = (TH2F*) subList->FindObject("R2SP_TPCFT0A_CentFT0C");
    TH2F *hs_R2SPAC = (TH2F*) subList->FindObject("R2SP_TPCFT0C_CentFT0C");
    TH2F *hs_R2SPBC = (TH2F*) subList->FindObject("R2SP_FT0AFT0C_CentFT0C");

    TAxis *axis = (TAxis*) hs_R2SPAB->GetXaxis();
    int nCentrBins = axis->GetNbins();
    double *centrBins = new double[nCentrBins + 1];
    axis->GetLowEdge(centrBins);
    centrBins[nCentrBins] = axis->GetBinUpEdge(nCentrBins);

    TH1D *hist_r2sp = new TH1D("R2SP_Cent", "R_{2}^{SP}", nCentrBins, centrBins);
    hist_r2sp->GetXaxis()->SetTitle("Centrality FT0C(%)");
    hist_r2sp->GetYaxis()->SetTitle("R_{2}{SP}");

    for (int iBin = 0; iBin < nCentrBins; iBin++) {
        std::cout << centrBins[iBin] << " " << centrBins[iBin+1] << std::endl;
        hs_R2SPAB->GetXaxis()->SetRangeUser(centrBins[iBin], centrBins[iBin+1]);
        hs_R2SPAC->GetXaxis()->SetRangeUser(centrBins[iBin], centrBins[iBin+1]);
        hs_R2SPBC->GetXaxis()->SetRangeUser(centrBins[iBin], centrBins[iBin+1]);

        double R2SPAB = hs_R2SPAB->GetMean(2);
        double R2SPAC = hs_R2SPAC->GetMean(2);
        double R2SPBC = hs_R2SPBC->GetMean(2);
        double R2SPABe = hs_R2SPAB->GetMean(12);
        double R2SPACe = hs_R2SPAC->GetMean(12);
        double R2SPBCe = hs_R2SPBC->GetMean(12);

        double R22SP = R2SPBC != 0 ? R2SPAB * R2SPAC / R2SPBC : 0.0;
        //double R22SP = R2SPBC != 0 ? R2SPBC * R2SPAC / R2SPAB : 0.0;
        double R2SP = R22SP > 0 ? TMath::Sqrt(R22SP) : 0.0;
        double R2SPe =
            R2SPAB * R2SPAC * R2SPBC == 0
                ? 0.0
                : TMath::Sqrt(
                      1. / 4 * (R2SPAC / (R2SPAB * R2SPBC)) * pow(R2SPABe, 2.) +
                      1. / 4 * (R2SPAB / (R2SPAC * R2SPBC)) * pow(R2SPACe, 2.) +
                      1. / 4 * R2SPAC * R2SPAB / pow(R2SPBC, 3.) *
                          pow(R2SPBCe, 2.));
        R2SPe = isnan(R2SPe) || isinf(R2SPe) ? 0. : R2SPe;

        std::cout << R2SP << " +/- " << R2SPe << std::endl;

        hist_r2sp->SetBinContent(iBin+1, R2SP);
        hist_r2sp->SetBinError(iBin+1, R2SPe);
    }
    hist_r2sp->Draw();
}