void LoadStyle();
void SetLegend(TLegend *);

void plot_results(TString model = "THU") {
    LoadStyle();
    gStyle -> SetLineStyleString(9, "80 20");

    double xMin, xMax, yCentral, yMin, yMax;

    std::ifstream fInFwd(Form("/Users/lucamicheletti/GITHUB/dq_run3_analyses/charmonia_production_pO_OO_NeNe/theory/%s/forward/data_v2_jpsi_010_y254.dat",model.Data()));

    std::vector<double> xCentrals, exMins, exMaxs, exZeros, yFwdCentrals, eyFwdMins, eyFwdMaxs, eyFwdZeros;
    while (fInFwd >> xMin >> xMax >> yCentral >> yMin >> yMax) {
        double xCentral = (xMax + xMin) / 2.;

        xCentrals.push_back(xCentral);
        exMins.push_back(xCentral - xMin);
        exMaxs.push_back(xMax - xCentral);
        exZeros.push_back(0);
        yFwdCentrals.push_back(yCentral / 100.);
        eyFwdMins.push_back((yCentral - yMin) / 100.);
        eyFwdMaxs.push_back((yMax - yCentral) / 100.);
        eyFwdZeros.push_back(0);
    }

    TFile *fInDzero = new TFile("v2Dzero.root", "READ");
    TGraphAsymmErrors *graStatV2DzeroMidCentr020 = (TGraphAsymmErrors*) fInDzero -> Get("gvn_prompt_stat");
    graStatV2DzeroMidCentr020 -> SetMarkerColor(kAzure+4);
    graStatV2DzeroMidCentr020 -> SetLineColor(kAzure+4);
    graStatV2DzeroMidCentr020 -> SetMarkerStyle(20);
    graStatV2DzeroMidCentr020 -> SetMarkerSize(1.5);

    TGraphAsymmErrors *graSyst1V2DzeroMidCentr020 = (TGraphAsymmErrors*) fInDzero -> Get("tot_syst");
    graSyst1V2DzeroMidCentr020 -> SetMarkerColor(kAzure+4);
    graSyst1V2DzeroMidCentr020 -> SetMarkerStyle(20);
    graSyst1V2DzeroMidCentr020 -> SetFillStyle(0);
    graSyst1V2DzeroMidCentr020 -> SetLineColor(kAzure+4);
    graSyst1V2DzeroMidCentr020 -> SetMarkerSize(1.5);

    for (int iPoint = 0;iPoint < 11;iPoint++) {
        graSyst1V2DzeroMidCentr020 -> SetPointEXhigh(iPoint, 0.15);
        graSyst1V2DzeroMidCentr020 -> SetPointEXlow(iPoint, 0.15);
    }

    TGraphAsymmErrors *graSyst2V2DzeroMidCentr020 = (TGraphAsymmErrors*) fInDzero -> Get("tot_syst_wnonflow");
    graSyst2V2DzeroMidCentr020 -> SetMarkerColor(kAzure+1);
    graSyst2V2DzeroMidCentr020 -> SetMarkerStyle(20);
    graSyst2V2DzeroMidCentr020 -> SetFillStyle(0);
    graSyst2V2DzeroMidCentr020 -> SetLineColor(kAzure+1);
    graSyst2V2DzeroMidCentr020 -> SetMarkerSize(1.5);

    TGraphAsymmErrors *graTheorFwdCentrVal = new TGraphAsymmErrors(xCentrals.size(), &(xCentrals[0]), &(yFwdCentrals[0]), &(exZeros[0]), &(exZeros[0]), &(exZeros[0]), &(eyFwdZeros[0]));
    gStyle -> SetLineStyleString(9,"80 20");
    graTheorFwdCentrVal -> SetLineColor(kOrange+7);
    graTheorFwdCentrVal -> SetLineWidth(2);
    graTheorFwdCentrVal -> SetLineStyle(9);

    double ptJpsiFwdCentrBins[] = {0.5, 1.5, 2.5, 3.5, 5, 7};
    double ptJpsiFwdWidthBins[] = {0.5, 0.5, 0.5, 0.5, 1, 1};
    double ptJpsiFwdSystBins[] = {0.15, 0.15, 0.15, 0.15, 0.15, 0.15};
    
    double v2JpsiFwdCentr020Vals[] = {-0.02049, 0.01836, 0.00113, 0.02118, 0.03563, 0.06178};
    double v2JpsiFwdCentr020Stats[] = {0.02700, 0.01943, 0.02102, 0.02533, 0.02749, 0.04121};
    double v2JpsiFwdCentr020Systs[] = {0.00243, 0.00103, 0.00316, 0.00343, 0.00300, 0.00543};

    TGraphErrors *graStatV2JpsiFwdCentr020 = new TGraphErrors(6, ptJpsiFwdCentrBins, v2JpsiFwdCentr020Vals, ptJpsiFwdWidthBins, v2JpsiFwdCentr020Stats);
    graStatV2JpsiFwdCentr020 -> SetMarkerColor(kRed+1);
    graStatV2JpsiFwdCentr020 -> SetLineColor(kRed+1);
    graStatV2JpsiFwdCentr020 -> SetLineWidth(2);
    graStatV2JpsiFwdCentr020 -> SetMarkerStyle(20);
    graStatV2JpsiFwdCentr020 -> SetMarkerSize(1.5);
    
    TGraphErrors *graSystV2JpsiFwdCentr020 = new TGraphErrors(6, ptJpsiFwdCentrBins, v2JpsiFwdCentr020Vals, ptJpsiFwdSystBins, v2JpsiFwdCentr020Systs);
    graSystV2JpsiFwdCentr020 -> SetMarkerColor(kRed+1);
    graSystV2JpsiFwdCentr020 -> SetMarkerStyle(20);
    graSystV2JpsiFwdCentr020 -> SetFillStyle(0);
    graSystV2JpsiFwdCentr020 -> SetLineColor(kRed+1);
    graSystV2JpsiFwdCentr020 -> SetLineWidth(2);
    graSystV2JpsiFwdCentr020 -> SetMarkerSize(1.5);

    double v2JpsiFwdCentr010Vals[] = {-0.00640, 0.04541, 0.03017, 0.04114, 0.05747, 0.04378};
    double v2JpsiFwdCentr010Stats[] = {0.03459, 0.02572, 0.02809, 0.03341, 0.03487, 0.05330};
    double v2JpsiFwdCentr010Systs[] = {0.00258, 0.00194, 0.00590, 0.00411, 0.00545, 0.01101};

    TGraphErrors *graStatV2JpsiFwdCentr010 = new TGraphErrors(6, ptJpsiFwdCentrBins, v2JpsiFwdCentr010Vals, ptJpsiFwdWidthBins, v2JpsiFwdCentr010Stats);
    graStatV2JpsiFwdCentr010 -> SetMarkerColor(kOrange+7);
    graStatV2JpsiFwdCentr010 -> SetLineColor(kOrange+7);
    graStatV2JpsiFwdCentr010 -> SetLineWidth(2);
    graStatV2JpsiFwdCentr010 -> SetMarkerStyle(20);
    graStatV2JpsiFwdCentr010 -> SetMarkerSize(1.5);
    
    TGraphErrors *graSystV2JpsiFwdCentr010 = new TGraphErrors(6, ptJpsiFwdCentrBins, v2JpsiFwdCentr010Vals, ptJpsiFwdSystBins, v2JpsiFwdCentr010Systs);    
    graSystV2JpsiFwdCentr010 -> SetMarkerColor(kOrange+7);
    graSystV2JpsiFwdCentr010 -> SetMarkerStyle(20);
    graSystV2JpsiFwdCentr010 -> SetFillStyle(0);
    graSystV2JpsiFwdCentr010 -> SetLineColor(kOrange+7);
    graSystV2JpsiFwdCentr010 -> SetLineWidth(2);
    graSystV2JpsiFwdCentr010 -> SetMarkerSize(1.5);

    double ptJpsiMidCentrBins[] = {1.230, 2.870, 5.300};
    double ptJpsiMidWidthBins[] = {0. , 0., 0.};
    double ptJpsiMidSystBins[] = {0.15, 0.15, 0.15};
    
    double v2JpsiMidCentr020Vals[] = {-0.180000, -0.018714, 0.059150};
    double v2JpsiMidCentr020Stats[] = {0.102775, 0.102947, 0.058341};
    double v2JpsiMidCentr020Systs[] = {0.050000, 0.090000, 0.060000};

    TGraphErrors *graStatV2JpsiMidCentr020 = new TGraphErrors(3, ptJpsiFwdCentrBins, v2JpsiMidCentr020Vals, ptJpsiFwdWidthBins, v2JpsiMidCentr020Stats);
    graStatV2JpsiMidCentr020 -> SetMarkerColor(kGreen+2);
    graStatV2JpsiMidCentr020 -> SetLineColor(kGreen+2);
    graStatV2JpsiMidCentr020 -> SetLineWidth(2);
    graStatV2JpsiMidCentr020 -> SetMarkerStyle(20);
    graStatV2JpsiMidCentr020 -> SetMarkerSize(1.5);
    
    TGraphErrors *graSystV2JpsiMidCentr020 = new TGraphErrors(3, ptJpsiFwdCentrBins, v2JpsiMidCentr020Vals, ptJpsiFwdSystBins, v2JpsiMidCentr020Systs);    
    graSystV2JpsiMidCentr020 -> SetMarkerColor(kGreen+2);
    graSystV2JpsiMidCentr020 -> SetMarkerStyle(20);
    graSystV2JpsiMidCentr020 -> SetFillStyle(0);
    graSystV2JpsiMidCentr020 -> SetLineColor(kGreen+2);
    graSystV2JpsiMidCentr020 -> SetLineWidth(2);
    graSystV2JpsiMidCentr020 -> SetMarkerSize(1.5);


    TLine *lineUnity = new TLine(0, 0, 8, 0);
    lineUnity -> SetLineColor(kGray+2);
    lineUnity -> SetLineWidth(2);
    lineUnity -> SetLineStyle(kDashed);

    TLatex *latexTitle = new TLatex();
    latexTitle -> SetTextSize(0.050);
    latexTitle -> SetNDC();
    latexTitle -> SetTextFont(42);

    // Jpsi Fwd vs theory
    TCanvas *canvasJpsiFwdVsTheor = new TCanvas("canvasJpsiFwdVsTheor", "", 800, 600);
    TH2D *histGridJpsiFwdVsTheor  = new TH2D("histGridJpsiFwdVsTheor", ";#it{p}_{T} (GeV/#it{c});#it{v}_{2}", 100, 0, 8, 100, -0.07, 0.25);
    histGridJpsiFwdVsTheor -> Draw();
    lineUnity -> Draw("SAME");
    graTheorFwdCentrVal -> Draw("L SAME");
    graSystV2JpsiFwdCentr020 -> Draw("E2P");
    graStatV2JpsiFwdCentr020 -> Draw("P");
    graSystV2JpsiFwdCentr010 -> Draw("E2P");
    graStatV2JpsiFwdCentr010 -> Draw("P");


    TLegend *legend1 = new TLegend(0.45,0.60,0.65,0.80);
    SetLegend(legend1);
    legend1 -> AddEntry(graTheorFwdCentrVal,"THU model, 0#minus10%","L");
    legend1 -> AddEntry(graSystV2JpsiFwdCentr010,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4, 0#minus10%","P");
    legend1 -> AddEntry(graSystV2JpsiFwdCentr020,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4, 0#minus20%","P");
    legend1 -> Draw();

    latexTitle -> DrawLatex(0.20, 0.88, "ALICE Work in Progress, OO, #sqrt{#it{s}_{NN}} = 5.36 TeV");

    // Jpsi Fwd vs Jpsi Mid
    TCanvas *canvasJpsiFwdVsJpsiMid = new TCanvas("canvasJpsiFwdVsJpsiMid", "", 800, 600);
    TH2D *histGridJpsiFwdVsJpsiMid  = new TH2D("histGridJpsiFwdVsJpsiMid", ";#it{p}_{T} (GeV/#it{c});#it{v}_{2}", 100, 0, 8, 100, -0.3, 0.30);
    histGridJpsiFwdVsJpsiMid -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiMidCentr020 -> Draw("E2P");
    graStatV2JpsiMidCentr020 -> Draw("P");
    graSystV2JpsiFwdCentr020 -> Draw("E2P");
    graStatV2JpsiFwdCentr020 -> Draw("P");


    TLegend *legend2 = new TLegend(0.45,0.20,0.65,0.40);
    SetLegend(legend2);
    legend2 -> AddEntry(graSystV2JpsiFwdCentr020,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4, 0#minus20%","P");
    legend2 -> AddEntry(graSystV2JpsiMidCentr020,"J/#psi #rightarrow e^{+}e^{-}, |#it{y}| < 0.9, 0#minus20%","P");
    legend2 -> Draw();

    latexTitle -> DrawLatex(0.20, 0.88, "ALICE Work in Progress, OO, #sqrt{#it{s}_{NN}} = 5.36 TeV");


    // Jpsi Fwd vs Dzero Mid
    TCanvas *canvasJpsiFwdVsDzeroMid = new TCanvas("canvasJpsiFwdVsDzeroMid", "", 800, 600);
    TH2D *histGridJpsiFwdVsDzeroMid  = new TH2D("histGridJpsiFwdVsDzeroMid", ";#it{p}_{T} (GeV/#it{c});#it{v}_{2}", 100, 0, 8, 100, -0.07, 0.25);
    histGridJpsiFwdVsDzeroMid -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiFwdCentr020 -> Draw("E2P");
    graStatV2JpsiFwdCentr020 -> Draw("P");
    graSyst1V2DzeroMidCentr020 -> Draw("E2P");
    //graSyst2V2DzeroMidCentr020 -> Draw("E2P");
    graStatV2DzeroMidCentr020 -> Draw("EP");

    TLegend *legend3 = new TLegend(0.18,0.70,0.48,0.85);
    SetLegend(legend3);
    legend3 -> SetTextSize(0.045);
    legend3 -> AddEntry(graSystV2JpsiFwdCentr020,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4, 0#minus20%","P");
    legend3 -> AddEntry(graStatV2DzeroMidCentr020,"D^{0} #rightarrow K^{-}#pi^{+} and charge conj., |#it{y}| < 0.8, 0#minus20%","P");
    legend3 -> Draw();

    latexTitle -> DrawLatex(0.20, 0.88, "ALICE Work in Progress, OO, #sqrt{#it{s}_{NN}} = 5.36 TeV");


    canvasJpsiFwdVsTheor -> SaveAs("figures/jpsiFwd_vs_models.pdf");
    canvasJpsiFwdVsJpsiMid -> SaveAs("figures/jpsiFwd_vs_jpsiMid.pdf");
    canvasJpsiFwdVsDzeroMid -> SaveAs("figures/jpsiFwd_vs_dzeroMid.pdf");

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetLegend(TLegend *legend){
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.045);
}