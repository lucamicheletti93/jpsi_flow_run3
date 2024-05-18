void LoadStyle();
void SetLegend(TLegend *);

void performance_plot() {
    LoadStyle();
    gStyle -> SetOptStat(0);

    const double resolution = 0.65200005;

    TFile *fIn = new TFile("/Users/lucamicheletti/GITHUB/jpsi_flow_run3/macros/Histograms_10_50_pT_2-3.root");

    TH1D *histMassSEPM = (TH1D*) fIn -> Get("histMassSEPM_2_3__10_50");
    TH1D *histV2SEPM = (TH1D*) fIn -> Get("histV2SEPM_2_3__10_50");
    histV2SEPM -> Scale(1. / resolution);

    TH1D *histMassMEPM = (TH1D*) fIn -> Get("histMassMEPM_2_3__10_50");
    TH1D *histV2MEPM = (TH1D*) fIn -> Get("histV2MEPM_2_3__10_50");
    histV2MEPM -> Scale(1. / resolution);

    std::cout << histV2MEPM -> GetBinContent(histV2MEPM -> FindBin(3.)) << std::endl;

    histMassSEPM -> SetTitle("");
    histV2SEPM -> SetTitle("");

    histMassSEPM -> GetXaxis() -> SetRangeUser(2.5, 4.5);
    histMassSEPM -> GetYaxis() -> SetRangeUser(5e3, 1e6);
    histMassSEPM -> GetYaxis() -> SetLabelSize(0.05);
    histMassSEPM -> GetYaxis() -> SetTitle("d#it{N} / d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");
    histMassSEPM -> GetYaxis() -> SetTitleSize(0.06);

    histV2SEPM -> GetXaxis() -> SetRangeUser(2.5, 4.5);
    histV2SEPM -> GetXaxis() -> SetLabelSize(0.08);
    histV2SEPM -> GetXaxis() -> SetTitle("");
    histV2SEPM -> GetXaxis() -> SetTitleSize(0.07);
    histV2SEPM -> GetYaxis() -> SetRangeUser(0.016, 0.074);
    histV2SEPM -> GetYaxis() -> SetLabelSize(0.07);
    histV2SEPM -> GetYaxis() -> CenterTitle(true);
    histV2SEPM -> GetYaxis() -> SetTitle("#it{v}_{2} {EP, |#Delta#eta| > 1.1}");
    histV2SEPM -> GetYaxis() -> SetTitleOffset(0.75);
    histV2SEPM -> GetYaxis() -> SetTitleSize(0.08);

    histMassSEPM -> SetLineColor(kBlack);
    histMassMEPM -> SetLineColor(kBlack);
    histMassMEPM -> SetMarkerStyle(24);
    histMassMEPM -> SetMarkerSize(0.80);
    histV2SEPM -> SetLineColor(kBlack);
    histV2MEPM -> SetLineColor(kBlack);
    histV2MEPM -> SetMarkerStyle(24);
    histV2MEPM -> SetMarkerSize(0.80);

    TFile *fFuncIn = new TFile("Fit_Function_Cent10-50_Pt2_3Bin.root");
    TF1 *funcMassSigBkg = (TF1*) fFuncIn -> Get("Mass_Signal_plus_bkg");
    funcMassSigBkg -> SetLineColor(kRed+1);
    TF1 *funcMassBkg = (TF1*) fFuncIn -> Get("Mass_bkg");
    funcMassBkg -> SetLineColor(kBlue+1);
    TF1 *funcV2SigBkg = (TF1*) fFuncIn -> Get("v2_Signal_plus_bkg");
    funcV2SigBkg -> SetLineColor(kRed+1);
    TF1 *funcV2Bkg = (TF1*) fFuncIn -> Get("v2_bkg");
    funcV2Bkg -> SetLineColor(kBlue+1);

    TCanvas *canvasMassV2Pt = new TCanvas("canvasMassV2Pt", "", 800, 1800);
    
    canvasMassV2Pt -> cd();
    gStyle -> SetOptStat(0);
    TPad *pad1 = new TPad("pad1", "pad1", 0.05, 0.4, 0.95, 0.95);
    pad1 -> SetBottomMargin(0);
    pad1 -> Draw();
    pad1 -> cd();
    gPad -> SetLogy(true);
    histMassSEPM -> Draw("EP");
    histMassMEPM -> Draw("EP SAME");
    funcMassSigBkg -> Draw("SAME");
    //funcMassBkg -> Draw("SAME");

    TLegend *legendHist = new TLegend(0.58, 0.35, 0.75, 0.60, " ", "brNDC");
    SetLegend(legendHist);
    legendHist -> SetTextSize(0.055);
    legendHist -> SetHeader("Opposite-sign pairs");
    legendHist -> AddEntry(histMassSEPM,"Same-event", "PL");
    legendHist -> AddEntry(histMassMEPM,"Mixed-event", "PL");
    legendHist -> Draw();

    TLegend *legendFunc = new TLegend(0.20, 0.15, 0.37, 0.40, " ", "brNDC");
    SetLegend(legendFunc);
    legendFunc -> SetTextSize(0.055);
    legendFunc -> AddEntry(funcV2SigBkg,"Signal", "L");
    legendFunc -> AddEntry(funcV2Bkg,"Background", "L");
    legendFunc -> Draw();

    TLatex *latexTitle = new TLatex();
    latexTitle -> SetTextSize(0.06);
    latexTitle -> SetNDC();
    latexTitle -> SetTextFont(42);
    latexTitle -> DrawLatex(0.22, 0.87, "ALICE performance, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.50, 0.77, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < y < 4");
    latexTitle -> DrawLatex(0.50, 0.67, "10#minus50\%, 2 < #it{p}_{T} < 3 GeV/#it{c}");

    canvasMassV2Pt -> cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.05, 0.05, 0.95, 0.4);
    pad2 -> SetTopMargin(0);
    pad2 -> Draw();
    pad2 -> cd();
    histV2SEPM -> Draw("EP");
    histV2MEPM -> Draw("EP SAME");
    funcV2Bkg -> Draw("SAME");
    funcV2SigBkg -> Draw("SAME");

    canvasMassV2Pt -> cd();
    TLatex *latexAxis = new TLatex();
    latexAxis -> SetTextSize(0.04);
    latexAxis -> SetNDC();
    latexAxis -> SetTextFont(42);
    latexAxis -> DrawLatex(0.70, 0.030, "#it{m}_{#mu#mu} (GeV/#it{c}^{2})");

    canvasMassV2Pt -> Update();
    canvasMassV2Pt -> SaveAs("flow_performance.pdf");
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