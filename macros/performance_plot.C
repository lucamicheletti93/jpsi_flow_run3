void LoadStyle();
void SetLegend(TLegend *);

void performance_plot() {
    LoadStyle();
    gStyle -> SetOptStat(0);

    TFile *fIn = new TFile("/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/data/pass2/LHC23_golden/dimuon_train_201824_no_cuts/Histograms_matchedMchMid_centr_10_50.root");
    fIn -> ls();
    TH1D *histMassSEPM = (TH1D*) fIn -> Get("histMassSEPM_2_3__10_50");
    TH1D *histV2SEPM = (TH1D*) fIn -> Get("histV2SEPM_2_3__10_50");
    TH1D *histPtSEPM = (TH1D*) fIn -> Get("histPtSEPM_2_3__10_50");

    TH1D *histMassMEPM = (TH1D*) fIn -> Get("histMassMEPM_2_3__10_50");
    TH1D *histV2MEPM = (TH1D*) fIn -> Get("histV2MEPM_2_3__10_50");
    TH1D *histPtMEPM = (TH1D*) fIn -> Get("histPtMEPM_2_3__10_50");

    histMassSEPM -> SetTitle("");
    histV2SEPM -> SetTitle("");
    histPtSEPM -> SetTitle("");

    histMassSEPM -> GetXaxis() -> SetRangeUser(2., 4.9);
    histMassSEPM -> GetYaxis() -> SetLabelSize(0.04);
    histMassSEPM -> GetYaxis() -> SetTitle("d#it{N} / d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");
    histMassSEPM -> GetYaxis() -> SetTitleSize(0.05);

    histV2SEPM -> GetXaxis() -> SetRangeUser(2., 4.9);
    histV2SEPM -> GetXaxis() -> SetLabelSize(0.06);
    histV2SEPM -> GetXaxis() -> SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    histV2SEPM -> GetXaxis() -> SetTitleSize(0.07);
    histV2SEPM -> GetYaxis() -> SetLabelSize(0.06);
    histV2SEPM -> GetYaxis() -> CenterTitle(true);
    histV2SEPM -> GetYaxis() -> SetTitle("#it{v}_{2} (EP)");
    histV2SEPM -> GetYaxis() -> SetTitleOffset(0.75);
    histV2SEPM -> GetYaxis() -> SetTitleSize(0.07);

    histPtSEPM -> GetXaxis() -> SetRangeUser(2., 4.9);

    histV2MEPM -> SetLineColor(kAzure+2);
    histMassSEPM -> SetLineColor(kAzure+2);

    TCanvas *canvasMassV2Pt = new TCanvas("canvasMassV2Pt", "", 800, 1800);
    
    canvasMassV2Pt -> cd();
    gStyle -> SetOptStat(0);
    TPad *pad1 = new TPad("pad1", "pad1", 0.05, 0.4, 0.95, 0.95);
    pad1 -> SetBottomMargin(0);
    pad1 -> Draw();
    pad1 -> cd();
    gPad -> SetLogy(true);
    histMassSEPM -> Draw("EP");
    histMassMEPM -> Draw("H SAME");

    TLatex *latexTitle = new TLatex();
    latexTitle -> SetTextSize(0.06);
    latexTitle -> SetNDC();
    latexTitle -> SetTextFont(42);
    latexTitle -> DrawLatex(0.22, 0.87, "ALICE performance, Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.50, 0.77, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < y < 4");
    latexTitle -> DrawLatex(0.50, 0.67, "10#minus50\%, 2 < #it{p}_{T} < 3 GeV/#it{c}");

    canvasMassV2Pt -> cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.05, 0.05, 0.95, 0.4);
    pad2 -> SetTopMargin(0);
    pad2 -> Draw();
    pad2 -> cd();
    histV2SEPM -> Draw("EP");
    histV2MEPM -> Draw("H SAME");

    canvasMassV2Pt -> Update();
    canvasMassV2Pt -> SaveAs("performance_plot_test.pdf");

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