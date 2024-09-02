void LoadStyle();
void SetLegend(TLegend *);

void plot_v2() {
    LoadStyle();
    gROOT -> SetStyle("Plain");
    gStyle -> SetOptStat(0);

    const int nMassBins = 50;
    const double minMassRange = 2;
    const double maxMassRange = 5;

    TLatex *latexTitle = new TLatex(); 
    latexTitle -> SetTextSize(0.08); 
    latexTitle -> SetNDC(); 
    latexTitle -> SetTextFont(42);

    TLatex *latexAxisMass = new TLatex(); 
    latexAxisMass -> SetTextSize(0.03); 
    latexAxisMass -> SetNDC(); 
    latexAxisMass -> SetTextFont(42);

    TLine *lineZero = new TLine(minMassRange, 0, maxMassRange, 0);
    lineZero -> SetLineStyle(kDashed);
    lineZero -> SetLineColor(kGray+2);
    lineZero -> SetLineWidth(2);

    TFile *fIn = new TFile("v2_results.root", "READ");

    TH1D *histProjMass = (TH1D*) fIn -> Get("histProjMass");
    histProjMass -> SetMarkerStyle(20);
    histProjMass -> SetMarkerSize(0.8);
    histProjMass -> SetMarkerColor(kBlack);
    histProjMass -> SetLineColor(kBlack);
    histProjMass -> GetYaxis() -> SetLabelSize(0.07);

    TH1D *histProjPt = (TH1D*) fIn -> Get("histProjPt");
    histProjPt -> SetMarkerStyle(24);
    histProjPt -> SetMarkerSize(0.8);
    histProjPt -> SetMarkerColor(kBlack);
    histProjPt -> SetLineColor(kBlack);
    histProjPt -> GetXaxis() -> SetLabelSize(0.1);
    histProjPt -> GetYaxis() -> SetLabelSize(0.07);
    histProjPt -> GetYaxis() -> CenterTitle(true);
    histProjPt -> GetYaxis() -> SetTitleOffset(0.4);
    histProjPt -> GetYaxis() -> SetTitleSize(0.1);
    histProjPt -> GetYaxis() -> SetTitle("<#it{p}_{T}>");

    TH1F *histV2SP = (TH1F*) fIn -> Get("histV2SP");
    histV2SP -> SetMarkerStyle(24);
    histV2SP -> SetMarkerSize(0.8);
    histV2SP -> SetMarkerColor(kBlack);
    histV2SP -> SetLineColor(kBlack);
    histV2SP -> GetXaxis() -> SetLabelSize(0.07);
    histV2SP -> GetYaxis() -> SetRangeUser(-0.95, 0.95);
    histV2SP -> GetYaxis() -> SetLabelSize(0.07);
    histV2SP -> GetYaxis() -> CenterTitle(true);
    histV2SP -> GetYaxis() -> SetTitleOffset(0.4);
    histV2SP -> GetYaxis() -> SetTitleSize(0.1);
    histV2SP -> GetYaxis() -> SetTitle("<#it{v}_{2} SP>");

    TH1F *histV2EP = (TH1F*) fIn -> Get("histV2EP");
    histV2EP -> SetMarkerStyle(24);
    histV2EP -> SetMarkerSize(0.8);
    histV2EP -> SetMarkerColor(kBlack);
    histV2EP -> SetLineColor(kBlack);
    histV2EP -> GetXaxis() -> SetLabelSize(0.07);
    histV2EP -> GetYaxis() -> SetRangeUser(-0.18, 0.18);
    histV2EP -> GetYaxis() -> SetLabelSize(0.07);
    histV2EP -> GetYaxis() -> CenterTitle(true);
    histV2EP -> GetYaxis() -> SetTitleOffset(0.4);
    histV2EP -> GetYaxis() -> SetTitleSize(0.1);
    histV2EP -> GetYaxis() -> SetTitle("<#it{v}_{2} EP>");

    // Scalar-Product results
    TCanvas *canvasMassVsV2SP = new TCanvas("canvasMassVsV2SP", "", 800, 700);
    gStyle -> SetOptTitle(0);
    
    TPad *padMassSP = new TPad("padMassSP", "padMassSP", 0.1, 0.6, 0.9, 0.9, 0, 0, 0);
    padMassSP -> SetBottomMargin(0);
    padMassSP -> Draw();
    
    TPad *padV2SP = new TPad("padV2SP", "padV2SP", 0.1, 0.35, 0.9, 0.6, 0, 0, 0);
    padV2SP -> SetTopMargin(0);
    padV2SP -> SetBottomMargin(0);
    padV2SP -> Draw();

    TPad *padPtSP = new TPad("padPtSP", "padPtSP", 0.1, 0.1, 0.9, 0.35, 0, 0, 0);
    padPtSP -> SetTopMargin(0);
    padPtSP -> Draw();

    
    padMassSP -> cd(); 
    histProjMass -> Draw("EP");
    latexTitle -> DrawLatex(0.55, 0.80, "ALICE, pp #sqrt{s_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.60, 0.70, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < y < 4");
    latexTitle -> DrawLatex(0.65, 0.60, "0#minus50\%, Scalar Product");

    padV2SP -> cd(); 
    histV2SP -> Draw("EP"); 
    lineZero -> Draw("SAME");

    padPtSP -> cd(); 
    histProjPt -> Draw("EP");

    canvasMassVsV2SP -> cd();
    latexAxisMass -> DrawLatex(0.70, 0.07, "#it{M}_{#mu#mu} (GeV/c^{2})");

    // Event-Plane results
    TCanvas *canvasMassVsV2EP = new TCanvas("canvasMassVsV2EP", "", 800, 700);
    gStyle -> SetOptTitle(0);
    
    TPad *padMassEP = new TPad("padMassEP", "padMassEP", 0.1, 0.6, 0.9, 0.9, 0, 0, 0);
    padMassEP -> SetBottomMargin(0);
    padMassEP -> Draw();
    
    TPad *padV2EP = new TPad("padV2EP", "padV2EP", 0.1, 0.35, 0.9, 0.6, 0, 0, 0);
    padV2EP -> SetTopMargin(0);
    padV2EP -> SetBottomMargin(0);
    padV2EP -> Draw();

    TPad *padPtEP = new TPad("padPtEP", "padPtEP", 0.1, 0.1, 0.9, 0.35, 0, 0, 0);
    padPtEP -> SetTopMargin(0);
    padPtEP -> Draw();

    
    padMassEP -> cd(); 
    histProjMass -> Draw("EP");
    latexTitle -> DrawLatex(0.55, 0.80, "ALICE, pp #sqrt{s_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.60, 0.70, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < y < 4");
    latexTitle -> DrawLatex(0.65, 0.60, "0#minus50\%, Event Plane");

    padV2EP -> cd(); 
    histV2EP -> Draw("EP"); 
    lineZero -> Draw("SAME");

    padPtEP -> cd(); 
    histProjPt -> Draw("EP");

    canvasMassVsV2EP -> cd();
    latexAxisMass -> DrawLatex(0.70, 0.07, "#it{M}_{#mu#mu} (GeV/c^{2})");

    canvasMassVsV2SP -> SaveAs("V2SP.pdf");
    canvasMassVsV2EP -> SaveAs("V2EP.pdf");

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