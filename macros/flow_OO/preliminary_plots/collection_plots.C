void LoadStyle();
void SetLegend(TLegend *);
void SetHist(TH1D *, Color_t , double , int , double , bool );
void SetGraph(TGraphErrors *, Color_t , double , int , double , bool );
void SetAsymmGraph(TGraphAsymmErrors *, Color_t , double , int , double , bool );
TGraphErrors* MakeSystGraph(TGraphErrors* graStat, TH1D* histSyst);

void collection_plots() {
    LoadStyle();

    // OO collisions: forward rapidity
    const int nPtBins = 4;
    double minPtBins[] = {0, 2, 4, 6};
    double maxPtBins[] = {2, 4, 6, 8};
    double ptJpsiFwdWidth[] = {0, 0, 0, 0};
    double ptJpsiFwdSyst[] = {0.15, 0.15, 0.15, 0.15};

    double ptJpsiOOFwdCentr020[] = {1.28, 2.87, 4.82, 6.84};
    double v2JpsiOOFwdCentr020Vals[] = {0.0143719, 0.0321096, 0.072157, 0.0819288};
    double v2JpsiOOFwdCentr020Stats[] = {0.0179254, 0.0187689, 0.0300562, 0.0488713};
    double v2JpsiOOFwdCentr020Systs[] = {0.00310963, 0.0101671, 0.0069517, 0.0301406};

    TGraphErrors *graStatV2JpsiOOFwdCentr020 = new TGraphErrors(nPtBins, ptJpsiOOFwdCentr020, v2JpsiOOFwdCentr020Vals, ptJpsiFwdWidth, v2JpsiOOFwdCentr020Stats);
    SetGraph(graStatV2JpsiOOFwdCentr020, kRed+1, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiOOFwdCentr020 = new TGraphErrors(nPtBins, ptJpsiOOFwdCentr020, v2JpsiOOFwdCentr020Vals, ptJpsiFwdSyst, v2JpsiOOFwdCentr020Systs);
    SetGraph(graSystV2JpsiOOFwdCentr020, kRed+1, 1.5, 20, 1, true);

    // ***************************************************************************************** //
    // LF results
    // ***************************************************************************************** //
    // Pion
    TFile *fInPiKp = new TFile("v2_piKp_OO_nch_0_5.root", "READ");
    TGraphErrors *graStatV2PiMidCentr010 = (TGraphErrors*) fInPiKp->Get("gv2_pi_OO_nch_0_5");
    SetGraph(graStatV2PiMidCentr010, kBlack, 1.5, 20, 1, false);

    TH1D *hSystV2PiMidCentr010 = (TH1D*) fInPiKp->Get("hv2_OO_pi_nch_0_5_syst_err");
    TGraphErrors *graSystV2PiMidCentr010 = MakeSystGraph(graStatV2PiMidCentr010, hSystV2PiMidCentr010);
    SetGraph(graSystV2PiMidCentr010, kBlack, 1.5, 20, 1, true);

    // Kaon
    TGraphErrors *graStatV2KaMidCentr010 = (TGraphErrors*) fInPiKp->Get("gv2_pi_OO_nch_0_5");
    SetGraph(graStatV2KaMidCentr010, kBlack, 1.5, 20, 1, false);

    TH1D *hSystV2KaMidCentr010 = (TH1D*) fInPiKp->Get("hv2_OO_pi_nch_0_5_syst_err");
    TGraphErrors *graSystV2KaMidCentr010 = MakeSystGraph(graStatV2PiMidCentr010, hSystV2KaMidCentr010);
    SetGraph(graSystV2KaMidCentr010, kBlack, 1.5, 20, 1, true);

    // Proton
    TGraphErrors *graStatV2PrMidCentr010 = (TGraphErrors*) fInPiKp->Get("gv2_pi_OO_nch_0_5");
    SetGraph(graStatV2PrMidCentr010, kBlack, 1.5, 20, 1, false);

    TH1D *hSystV2PrMidCentr010 = (TH1D*) fInPiKp->Get("hv2_OO_pi_nch_0_5_syst_err");
    TGraphErrors *graSystV2PrMidCentr010 = MakeSystGraph(graStatV2PiMidCentr010, hSystV2PrMidCentr010);
    SetGraph(graSystV2PrMidCentr010, kBlack, 1.5, 20, 1, true);

    // Lambda
    TFile *fInLambda = new TFile("v2Lambda.root", "READ");
    TGraphErrors *graStatV2LambdaMidCentr010 = (TGraphErrors*) fInLambda->Get("gist_lambda_010");
    SetGraph(graStatV2LambdaMidCentr010, kGreen+2, 1.5, 20, 1, false);

    TGraphErrors *graSyst1V2LambdaMidCentr010 = (TGraphErrors*) fInLambda->Get("syst_lambda_010");
    SetGraph(graSyst1V2LambdaMidCentr010, kGreen+2, 1.5, 20, 1, true);

    // KzeroS
    TFile *fInKzero = new TFile("v2K0.root", "READ");
    TGraphErrors *graStatV2KzeroMidCentr010 = (TGraphErrors*) fInKzero->Get("gist_k0_010");
    SetGraph(graStatV2KzeroMidCentr010, kOrange+7, 1.5, 20, 1, false);

    TGraphErrors *graSyst1V2KzeroMidCentr010 = (TGraphErrors*) fInKzero->Get("syst_k0_010");
    SetGraph(graSyst1V2KzeroMidCentr010, kOrange+7, 1.5, 20, 1, true);
    // ***************************************************************************************** //
    // D-meson results
    // ***************************************************************************************** //
    TFile *fInDzero = new TFile("v2Dzero_020.root", "READ");
    TGraphAsymmErrors *graStatV2DzeroMidCentr020 = (TGraphAsymmErrors*) fInDzero->Get("gvn_prompt_stat");
    SetAsymmGraph(graStatV2DzeroMidCentr020, kAzure+4, 1.5, 20, 1, false);

    TGraphAsymmErrors *graSyst1V2DzeroMidCentr020 = (TGraphAsymmErrors*) fInDzero->Get("tot_syst");
    SetAsymmGraph(graSyst1V2DzeroMidCentr020, kAzure+4, 1.5, 20, 1, true);

    for (int iPoint = 0;iPoint < 11;iPoint++) {
        graSyst1V2DzeroMidCentr020->SetPointEXhigh(iPoint, 0.15);
        graSyst1V2DzeroMidCentr020->SetPointEXlow(iPoint, 0.15);
    }

    // ======================================================================= //
    // Plots
    // ======================================================================= //
    TLine *lineUnity = new TLine(0, 0, 8, 0);
    lineUnity->SetLineColor(kGray+2);
    lineUnity->SetLineWidth(2);
    lineUnity->SetLineStyle(kDashed);

    TLatex *latexTitle = new TLatex();
    latexTitle->SetTextSize(0.045);
    latexTitle->SetNDC();
    latexTitle->SetTextFont(42);

    TLatex latexTable;

    TCanvas *canvasV2JpsiAllSpecies = new TCanvas("canvasV2JpsiAllSpecies", "", 800, 600);
    TH2D *histGridV2JpsiAllSpecies = new TH2D("histGridV2JpsiAllSpecies", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}", 100, 0, 8, 100, -0.03, 0.37);
    histGridV2JpsiAllSpecies->Draw();

    TLegend *legendV2JpsiAllSpecies1 = new TLegend(0.20,0.72,0.81,0.80);
    SetLegend(legendV2JpsiAllSpecies1);
    legendV2JpsiAllSpecies1->SetNColumns(2);
    legendV2JpsiAllSpecies1->AddEntry(graStatV2KzeroMidCentr010,"K_{S}^{0}","P");
    legendV2JpsiAllSpecies1->AddEntry(graStatV2LambdaMidCentr010,"#Lambda","P");
    legendV2JpsiAllSpecies1->Draw();

    TLegend *legendV2JpsiAllSpecies2 = new TLegend(0.20,0.60,0.82,0.68);
    SetLegend(legendV2JpsiAllSpecies2);
    legendV2JpsiAllSpecies2->SetNColumns(2);
    legendV2JpsiAllSpecies2->AddEntry(graStatV2DzeroMidCentr020,"Prompt D^{0}, |#Delta#it{#eta}| > 1.3","P");
    legendV2JpsiAllSpecies2->AddEntry(graStatV2JpsiOOFwdCentr020,"J/#psi, |#Delta#it{#eta}| > 1.7","P");
    legendV2JpsiAllSpecies2->Draw();

    latexTitle->DrawLatex(0.20, 0.88, "ALICE Preliminary, OO,  #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle->DrawLatex(0.20, 0.80, "2PC, 0#minus10%, 1.2 < |#Delta#it{#eta}| < 1.8");
    latexTitle->DrawLatex(0.20, 0.68, "SP, 0#minus20%");

    lineUnity->Draw("SAME");
    graStatV2KzeroMidCentr010->Draw("P SAME");
    graSyst1V2KzeroMidCentr010->Draw("E2P SAME");
    graStatV2LambdaMidCentr010->Draw("P SAME");
    graSyst1V2LambdaMidCentr010->Draw("E2P SAME");
    graStatV2DzeroMidCentr020->Draw("P SAME");
    graSyst1V2DzeroMidCentr020->Draw("E2P SAME");
    graStatV2JpsiOOFwdCentr020->Draw("P SAME");
    graSystV2JpsiOOFwdCentr020->Draw("E2P SAME");
    //graStatV2PiMidCentr010->Draw("P SAME");
    //graSystV2PiMidCentr010->Draw("E2P SAME");

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LoadStyle() {
    int font = 42;
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(10);
    gStyle->SetCanvasColor(10);
    gStyle->SetTitleFillColor(10);
    gStyle->SetTitleBorderSize(1);
    gStyle->SetStatColor(10);
    gStyle->SetStatBorderSize(1);
    gStyle->SetLegendBorderSize(1);
    gStyle->SetDrawBorder(0);
    gStyle->SetTextFont(font);
    gStyle->SetStatFontSize(0.05);
    gStyle->SetStatX(0.97);
    gStyle->SetStatY(0.98);
    gStyle->SetStatH(0.03);
    gStyle->SetStatW(0.3);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTickLength(0.02,"y");
    gStyle->SetEndErrorSize(3);
    gStyle->SetLabelSize(0.05,"xyz");
    gStyle->SetLabelFont(font,"xyz");
    gStyle->SetLabelOffset(0.01,"xyz");
    gStyle->SetTitleFont(font,"xyz");
    gStyle->SetTitleOffset(0.9,"x");
    gStyle->SetTitleOffset(1.02,"y");
    gStyle->SetTitleSize(0.05,"xyz");
    gStyle->SetMarkerSize(1.3);
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(0);
    gStyle->SetCanvasPreferGL(kTRUE);
    gStyle->SetHatchesSpacing(0.5);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(0.0);
    gStyle->SetTitleSize(0.05,"X");
    gStyle->SetTitleSize(0.045,"Y");
    gStyle->SetLabelSize(0.045,"X");
    gStyle->SetLabelSize(0.045,"Y");
    gStyle->SetTitleOffset(1.2,"X");
    gStyle->SetTitleOffset(1.35,"Y");
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetLegend(TLegend *legend) {
    legend->SetBorderSize(0);
    legend->SetFillColor(10);
    legend->SetFillStyle(1);
    legend->SetLineStyle(0);
    legend->SetLineColor(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.045);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetHist(TH1D *hist, Color_t color, double size, int style, double alpha = 1, bool fillStyle = false) {
    hist->SetMarkerStyle(style);
    hist->SetMarkerColorAlpha(color, alpha);
    hist->SetMarkerSize(size);
    hist->SetLineColorAlpha(color, alpha);
    hist->SetLineWidth(2);
    if (fillStyle) {hist->SetFillStyle(0);}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetGraph(TGraphErrors *gra, Color_t color, double size, int style, double alpha = 1, bool fillStyle = false) {
    gra->SetMarkerStyle(style);
    gra->SetMarkerColorAlpha(color, alpha);
    gra->SetMarkerSize(size);
    gra->SetLineColorAlpha(color, alpha);
    gra->SetLineWidth(2);
    if (fillStyle) {gra->SetFillStyle(0);}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetAsymmGraph(TGraphAsymmErrors *gra, Color_t color, double size, int style, double alpha = 1, bool fillStyle = false) {
    gra->SetMarkerStyle(style);
    gra->SetMarkerColorAlpha(color, alpha);
    gra->SetMarkerSize(size);
    gra->SetLineColorAlpha(color, alpha);
    gra->SetLineWidth(2);
    if (fillStyle) {gra->SetFillStyle(0);}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors* MakeSystGraph(TGraphErrors* graStat, TH1D* histSyst) {
    if (!graStat || !histSyst) return nullptr;

    int nBins = graStat->GetN();
    if (histSyst->GetNbinsX() < nBins) {
        return nullptr;
    }

    TGraphErrors *graSyst = new TGraphErrors(nBins);
    for (int i = 0;i < nBins;i++) {
        double x, y;
        graStat->GetPoint(i,x,y);
        double ex = graStat->GetErrorX(i);
        double ey = histSyst->GetBinContent(i+1);
        graSyst->SetPoint(i,x,y);
        graSyst->SetPointError(i,ex,ey);
    }

    return graSyst;
}