void LoadStyle();
void SetLegend(TLegend *);
void SetGraph(TGraphErrors *, Color_t , double , int , double , bool );
void SetAsymmGraph(TGraphAsymmErrors *, Color_t , double , int , double , bool );
TGraphAsymmErrors* DoGraphFromTheory(string);

void compare_results() {
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

    double ptJpsiOOFwdCentr2060[] = {1.24, 2.87, 4.83, 6.76};
    double v2JpsiOOFwdCentr2060Vals[] = {-0.00880563, 0.00844508, 0.0413613, 0.119999};
    double v2JpsiOOFwdCentr2060Stats[] = {0.0213792, 0.0220216, 0.0348384, 0.0635753};
    double v2JpsiOOFwdCentr2060Systs[] = {0.0069356, 0.00718644, 0.00866651, 0.0273761};

    TGraphErrors *graStatV2JpsiOOFwdCentr2060 = new TGraphErrors(nPtBins, ptJpsiOOFwdCentr2060, v2JpsiOOFwdCentr2060Vals, ptJpsiFwdWidth, v2JpsiOOFwdCentr2060Stats);
    SetGraph(graStatV2JpsiOOFwdCentr2060, kRed+1, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiOOFwdCentr2060 = new TGraphErrors(nPtBins, ptJpsiOOFwdCentr2060, v2JpsiOOFwdCentr2060Vals, ptJpsiFwdSyst, v2JpsiOOFwdCentr2060Systs);
    SetGraph(graSystV2JpsiOOFwdCentr2060, kRed+1, 1.5, 20, 1, true);


    double v2JpsiOOFwdCentr020ValsEm[] = {0.01037, 0.03796, 0.07001, 0.11965};
    double v2JpsiOOFwdCentr020StatsEm[] = {0.01667, 0.01740, 0.02914, 0.04731};
    double v2JpsiOOFwdCentr020SystsEm[] = {0.00101, 0.00378, 0.00744, 0.00405};

    TGraphErrors *graStatV2JpsiOOFwdCentr020Em = new TGraphErrors(nPtBins, ptJpsiOOFwdCentr020, v2JpsiOOFwdCentr020ValsEm, ptJpsiFwdWidth, v2JpsiOOFwdCentr020StatsEm);
    SetGraph(graStatV2JpsiOOFwdCentr020Em, kBlue+1, 1.5, 20, 1, false);

    TGraphErrors *graSystV2JpsiOOFwdCentr020Em = new TGraphErrors(nPtBins, ptJpsiOOFwdCentr020, v2JpsiOOFwdCentr020ValsEm, ptJpsiFwdSyst, v2JpsiOOFwdCentr020SystsEm);
    SetGraph(graSystV2JpsiOOFwdCentr020Em, kBlue+1, 1.5, 20, 1, true);


    double v2JpsiOOFwdCentr2060ValsEm[] = {-0.01060, 0.02185, 0.05566, 0.13816};
    double v2JpsiOOFwdCentr2060StatsEm[] = {0.02062, 0.02166, 0.03534, 0.06060};
    double v2JpsiOOFwdCentr2060SystsEm[] = {0.00085, 0.00481, 0.00325, 0.00536};

    TGraphErrors *graStatV2JpsiOOFwdCentr2060Em = new TGraphErrors(nPtBins, ptJpsiOOFwdCentr2060, v2JpsiOOFwdCentr2060ValsEm, ptJpsiFwdWidth, v2JpsiOOFwdCentr2060StatsEm);
    SetGraph(graStatV2JpsiOOFwdCentr2060Em, kBlue+1, 1.5, 20, 1, false);

    TGraphErrors *graSystV2JpsiOOFwdCentr2060Em = new TGraphErrors(nPtBins, ptJpsiOOFwdCentr2060, v2JpsiOOFwdCentr2060ValsEm, ptJpsiFwdSyst, v2JpsiOOFwdCentr2060SystsEm);
    SetGraph(graSystV2JpsiOOFwdCentr2060Em, kBlue+1, 1.5, 20, 1, true);

    TH2D *histGridV2JpsiOO = new TH2D("histGridV2JpsiOO", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}^{SP} {|#Delta#it{#eta}| > 1.7}", 100, 0, 8, 100, -0.05, 0.20);

    TLine *lineUnity = new TLine(0, 0, 8, 0);
    lineUnity->SetLineColor(kGray+2);
    lineUnity->SetLineWidth(2);
    lineUnity->SetLineStyle(kDashed);

    TF1 *funcPol0Centr020 = new TF1("funcPol0Centr020", "[0]", 0, 8);
    graSystV2JpsiOOFwdCentr020->Fit(funcPol0Centr020, "R");

    TF1 *funcPol1Centr020 = new TF1("funcPol1Centr020", "[0] + [1]*x", 0, 8);
    graSystV2JpsiOOFwdCentr020->Fit(funcPol1Centr020, "R");

    TCanvas *canvasV2JpsiOOFwdCentr020 = new TCanvas("canvasV2JpsiOOFwdCentr020", "", 800, 600);
    histGridV2JpsiOO->Draw();
    lineUnity->Draw("SAME");
    graSystV2JpsiOOFwdCentr020->Draw("E2P SAME");
    graStatV2JpsiOOFwdCentr020->Draw("P SAME");
    graSystV2JpsiOOFwdCentr020Em->Draw("E2P SAME");
    graStatV2JpsiOOFwdCentr020Em->Draw("P SAME");

    TLegend *legendV2JpsiOOFwdCentr020 = new TLegend(0.20,0.65,0.40,0.78);
    SetLegend(legendV2JpsiOOFwdCentr020);
    legendV2JpsiOOFwdCentr020->AddEntry(graSystV2JpsiOOFwdCentr020,"Preliminary results","P");
    legendV2JpsiOOFwdCentr020->AddEntry(graStatV2JpsiOOFwdCentr020Em,"New event mixing","P");
    legendV2JpsiOOFwdCentr020->Draw();

    canvasV2JpsiOOFwdCentr020-> SaveAs("comparison_0_20.pdf");


    TCanvas *canvasV2JpsiOOFwdCentr2060 = new TCanvas("canvasV2JpsiOOFwdCentr2060", "", 800, 600);
    histGridV2JpsiOO->Draw();
    lineUnity->Draw("SAME");
    graSystV2JpsiOOFwdCentr2060->Draw("E2P SAME");
    graStatV2JpsiOOFwdCentr2060->Draw("P SAME");
    graSystV2JpsiOOFwdCentr2060Em->Draw("E2P SAME");
    graStatV2JpsiOOFwdCentr2060Em->Draw("P SAME");

    TLegend *legendV2JpsiOOFwdCentr2060 = new TLegend(0.20,0.65,0.40,0.78);
    SetLegend(legendV2JpsiOOFwdCentr2060);
    legendV2JpsiOOFwdCentr2060->AddEntry(graSystV2JpsiOOFwdCentr2060,"Preliminary results","P");
    legendV2JpsiOOFwdCentr2060->AddEntry(graStatV2JpsiOOFwdCentr2060Em,"New event mixing","P");
    legendV2JpsiOOFwdCentr2060->Draw();

    canvasV2JpsiOOFwdCentr2060-> SaveAs("comparison_20_60.pdf");
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphAsymmErrors* DoGraphFromTheory(string fInName) {
    double xMin, xMax, yCentral, yMin, yMax;

    std::ifstream fInFwd(fInName.c_str());

    std::vector<double> xCentrals, exMins, exMaxs, exZeros, yCentrals, eyMins, eyMaxs, eyZeros;
    while (fInFwd >> xMin >> xMax >> yCentral >> yMin >> yMax) {
        double xCentral = (xMax + xMin) / 2.;
        //std::cout << yCentral << std::endl;

        xCentrals.push_back(xCentral);
        exMins.push_back(xCentral - xMin);
        exMaxs.push_back(xMax - xCentral);
        exZeros.push_back(0);
        yCentrals.push_back(yCentral / 100.);
        eyMins.push_back((yCentral - yMin) / 100.);
        eyMaxs.push_back((yMax - yCentral) / 100.);
        eyZeros.push_back(0);
    }

    TGraphAsymmErrors *graAsymm = new TGraphAsymmErrors(xCentrals.size(), &(xCentrals[0]), &(yCentrals[0]), &(exZeros[0]), &(exZeros[0]), &(eyZeros[0]), &(eyZeros[0]));
    return graAsymm;
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