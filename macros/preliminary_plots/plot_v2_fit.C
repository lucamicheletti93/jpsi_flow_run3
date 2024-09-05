void LoadStyle();
void SetLegend(TLegend *);

void plot_v2_fit(double minPtBin = 2, double maxPtBin = 3) {
    LoadStyle();
    gStyle -> SetOptStat(0);

    TFile *fIn = new TFile(Form("/Users/lucamicheletti/GITHUB/jpsi_flow_run3/macros/systematics_mix_fit/fitResults_Pt_%2.1f_%2.1f.root", minPtBin, maxPtBin));
    TCanvas *canvasFit = (TCanvas*) fIn -> Get("v2_fit_2.3_4.7_CB2_VWG_data_Pol2_v2bkg");

    TPad *padMassFit = (TPad*) canvasFit -> GetPrimitive("pad1");
    TH1D *histMassSEPM = (TH1D*) padMassFit -> GetPrimitive(Form("histMassSEPM_%1.0f_%1.0f__10_50", minPtBin, maxPtBin));
    TH1D *histMassMEPM = (TH1D*) padMassFit -> GetPrimitive(Form("histMassMEPM_%1.0f_%1.0f__10_50", minPtBin, maxPtBin));

    TF1 *funcMassSigBkg = (TF1*) padMassFit -> GetPrimitive("funcMassSigBkg");
    TF1 *funcMassSig = (TF1*) padMassFit -> GetPrimitive("funcMassSig");
    TF1 *funcMassBkg = (TF1*) padMassFit -> GetPrimitive("funcMassBkg");

    TPad *padFlowFit = (TPad*) canvasFit -> GetPrimitive("pad2");
    TF1 *funcFlowBkg = (TF1*) padFlowFit -> GetPrimitive("bck2");

    // The histogram has no names -> necessary to retrieve with primitives
    TList* primitives = padFlowFit -> GetListOfPrimitives();
    TIter next(primitives);
    TObject* obj = nullptr;
    TH1D *histFlowSEPM;
    TH1D *histFlowMEPM;

    int counter = 0;
    while ((obj = next())) {
        if (obj->InheritsFrom(TH1::Class())) {
            TH1* hist = (TH1*) obj;
            if (strcmp(hist -> GetName(), "") == 0) {
                if (counter == 0) {
                    histFlowSEPM = (TH1D*) hist;
                }
                if (counter == 1) {
                    histFlowMEPM = (TH1D*) hist;
                }
                counter++;
            }
        }
    }

    // Retrieve the fit result
    TCanvas *canvasSyst = (TCanvas*) fIn -> Get("canvasSyst");
    TH1D *histSyst = (TH1D*) canvasSyst -> GetPrimitive("histSyst");
    double jpsiV2 = histSyst -> GetBinContent(histSyst -> GetXaxis() -> FindBin("CB2 + VWG + data tails, 2.3 - 4.7, Pol2[v2 bkg]"));
    double errJpsiV2 = histSyst -> GetBinError(histSyst -> GetXaxis() -> FindBin("CB2 + VWG + data tails, 2.3 - 4.7, Pol2[v2 bkg]"));

    ///////////////////////////
    // Plot the fit
    ///////////////////////////
    TCanvas *canvasMassFlow = new TCanvas("canvasMassFlow", "", 800, 1800);
    
    canvasMassFlow -> cd();
    gStyle -> SetOptStat(0);
    
    // Pad 1 [J/psi mass fit]
    TPad *pad1 = new TPad("pad1", "pad1", 0.05, 0.5, 0.95, 0.95);
    pad1 -> SetBottomMargin(0);
    pad1 -> SetTopMargin(0.07);
    pad1 -> Draw();
    pad1 -> cd();
    pad1 -> SetTicks(1, 1);

    histMassSEPM -> SetTitle("");
    double minMassBin, maxMassBin;
    if (minPtBin < 4) {
        minMassBin = 0.01 * histMassSEPM -> GetBinContent(histMassSEPM -> FindBin(2.5));
        maxMassBin = 1.10 * histMassSEPM -> GetBinContent(histMassSEPM -> FindBin(2.5));
    } else {
        minMassBin = 0.01 * histMassSEPM -> GetBinContent(histMassSEPM -> FindBin(2.5));
        maxMassBin = 1.50 * histMassSEPM -> GetBinContent(histMassSEPM -> FindBin(2.5));
    }
    histMassSEPM -> GetXaxis() -> SetRangeUser(2.48, 4.52);
    histMassSEPM -> GetYaxis() -> SetRangeUser(minMassBin, maxMassBin);
    histMassSEPM -> GetYaxis() -> SetLabelSize(0.07);
    histMassSEPM -> GetYaxis() -> SetTitleOffset(0.90);
    histMassSEPM -> GetYaxis() -> SetTitle("d#it{N} / d#it{m}_{#mu#mu} (GeV/#it{c}^{2})^{-1}");
    histMassSEPM -> GetYaxis() -> SetTitleSize(0.08);
    histMassSEPM -> GetYaxis() -> SetMaxDigits(3);
    histMassSEPM -> SetMarkerSize(1);
    histMassSEPM -> SetLineColor(kBlack);

    histMassSEPM -> Draw("EP");

    funcMassBkg -> SetLineColor(kGray+1);
    funcMassBkg -> Draw("SAME");

    funcMassSigBkg -> SetLineColor(kBlue);
    funcMassSigBkg -> Draw("SAME");

    TLegend *legendMassFit = new TLegend(0.65, 0.33, 0.85, 0.68, " ", "brNDC");
    SetLegend(legendMassFit);
    legendMassFit -> SetTextSize(0.07);
    legendMassFit -> AddEntry(histMassSEPM, "Data", "PL");
    legendMassFit -> AddEntry(funcMassSigBkg, "Total fit", "L");
    legendMassFit -> AddEntry(funcMassBkg, "Background", "L");
    legendMassFit -> Draw();

    TLatex *latexTitle = new TLatex();
    latexTitle -> SetTextSize(0.07);
    latexTitle -> SetNDC();
    latexTitle -> SetTextFont(42);
    latexTitle -> DrawLatex(0.27, 0.85, "ALICE Preliminary");
    latexTitle -> DrawLatex(0.27, 0.75, "Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV, 10#minus50\%");
    latexTitle -> DrawLatex(0.27, 0.65, Form("J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4, %1.0f < #it{p}_{T} < %1.0f GeV/#it{c}", minPtBin, maxPtBin));

    canvasMassFlow -> cd();

    // Pad 1 [J/psi flow fit]
    TPad *pad2 = new TPad("pad2", "pad2", 0.05, 0.05, 0.95, 0.5);
    pad2 -> SetTopMargin(0);
    pad2 -> Draw();
    pad2 -> cd();
    pad2 -> SetTicks(1, 1);

    histFlowMEPM -> SetTitle("");
    double minFlowBin = 0.30 * histFlowSEPM -> GetBinContent(histFlowSEPM -> FindBin(2.5));
    double maxFlowBin = 1.30 * histFlowSEPM -> GetBinContent(histFlowSEPM -> FindBin(2.5));
    histFlowMEPM -> GetXaxis() -> SetRangeUser(2.48, 4.52);
    histFlowMEPM -> GetXaxis() -> SetLabelSize(0.08);
    histFlowMEPM -> GetXaxis() -> SetTitle("");
    histFlowMEPM -> GetXaxis() -> SetTitleSize(0.07);
    histFlowMEPM -> GetYaxis() -> SetRangeUser(minFlowBin, maxFlowBin);
    histFlowMEPM -> GetYaxis() -> SetLabelSize(0.07);
    histFlowMEPM -> GetYaxis() -> CenterTitle(true);
    histFlowMEPM -> GetYaxis() -> SetTitle("#it{v}_{2} {EP, |#Delta#eta| > 2.5}");
    histFlowMEPM -> GetYaxis() -> SetTitleOffset(0.90);
    histFlowMEPM -> GetYaxis() -> SetTitleSize(0.08);
    histFlowMEPM -> SetLineColor(kBlack);
    histFlowMEPM -> SetMarkerStyle(24);
    histFlowMEPM -> SetMarkerSize(1);
    histFlowMEPM -> Draw("EP");

    funcFlowBkg -> SetLineColor(kGray+1);
    funcFlowBkg -> Draw("SAME");

    histFlowSEPM -> SetMarkerSize(1);
    histFlowSEPM -> SetLineColor(kBlack);
    histFlowSEPM -> Draw("EP SAME");

    TLegend *legendFlowFit = new TLegend(0.65, 0.85, 0.88, 0.97, " ", "brNDC");
    SetLegend(legendFlowFit);
    legendFlowFit -> SetTextSize(0.07);
    legendFlowFit -> AddEntry(histFlowMEPM, "Mixed event", "PL");
    legendFlowFit -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, Form("v_{2}^{Sig}(J/#psi) = %4.3f#kern[0.5]{#pm} %4.3f", jpsiV2, errJpsiV2));
    //latexTitle -> DrawLatex(0.28, 0.77, "#chi^{2}/NDF = 1.1");

    canvasMassFlow -> cd();
    TLatex *latexAxis = new TLatex();
    latexAxis -> SetTextSize(0.04);
    latexAxis -> SetNDC();
    latexAxis -> SetTextFont(42);
    latexAxis -> DrawLatex(0.70, 0.030, "#it{m}_{#mu#mu} (GeV/#it{c}^{2})");

    canvasMassFlow -> Update();
    canvasMassFlow -> SaveAs(Form("fitFlowVsMassPt_%1.0f_%1.0f_MixedEvent.pdf", minPtBin, maxPtBin));
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