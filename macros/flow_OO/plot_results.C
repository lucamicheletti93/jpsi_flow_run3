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

    const int nPtBins = 6;
    double minPtBins[] = {0, 1, 2, 3, 4, 6};
    double maxPtBins[] = {1, 2, 3, 4, 6, 8};
    double ptJpsiFwdCentrBins[] = {0.5, 1.5, 2.5, 3.5, 5, 7};
    double ptJpsiFwdWidthBins[] = {0.5, 0.5, 0.5, 0.5, 1, 1};
    double ptJpsiFwdSystBins[] = {0.15, 0.15, 0.15, 0.15, 0.15, 0.15};
    
    // ***************************************************************************************** //
    // Centrality 0-20%
    // ***************************************************************************************** //
    // Results shown at the pag of 17/04/26
    //double v2JpsiFwdCentr020Vals[] = {-0.02049, 0.01836, 0.00113, 0.02118, 0.03563, 0.06178};
    //double v2JpsiFwdCentr020Stats[] = {0.02700, 0.01943, 0.02102, 0.02533, 0.02749, 0.04121};
    //double v2JpsiFwdCentr020Systs[] = {0.00243, 0.00103, 0.00316, 0.00343, 0.00300, 0.00543};
    // TPC-POS only Q-vectors
    //double v2JpsiFwdCentr020Vals[] = {-0.0290694, 0.0310567, 0.0101512, 0.0546731, 0.0819529, 0.0789928};
    //double v2JpsiFwdCentr020Stats[] = {0.0326985, 0.023158, 0.025574, 0.031402, 0.0316347, 0.0515139};
    //double v2JpsiFwdCentr020Systs[] = {0.00706876, 0.00573476, 0.0138293, 0.00413571, 0.00587116, 0.042338};
    // TPC-ALL Q-vectors
    double v2JpsiFwdCentr020Vals[] = {-0.0405239, 0.0301976, 0.00319262, 0.0744494, 0.0698103, 0.0844539};
    double v2JpsiFwdCentr020Stats[] = {0.0312219, 0.0221446, 0.0243449, 0.0300693, 0.0303064, 0.049283};
    double v2JpsiFwdCentr020Systs[] = {0.00287178, 0.00929753, 0.0180505, 0.00464562, 0.0194577, 0.0269159};

    Printf("J/psi v2 vs pT in centrality 0-20");
    for (int iPt = 0;iPt < nPtBins;iPt++) {
        Printf("%1.0f - %1.0f & %4.3f $?pm$ %4.3f $?pm$ %4.3f ??", minPtBins[iPt], maxPtBins[iPt], v2JpsiFwdCentr020Vals[iPt], v2JpsiFwdCentr020Stats[iPt], v2JpsiFwdCentr020Systs[iPt]);
    }

    TGraphErrors *graStatV2JpsiFwdCentr020 = new TGraphErrors(nPtBins, ptJpsiFwdCentrBins, v2JpsiFwdCentr020Vals, ptJpsiFwdWidthBins, v2JpsiFwdCentr020Stats);
    graStatV2JpsiFwdCentr020 -> SetMarkerColor(kRed+1);
    graStatV2JpsiFwdCentr020 -> SetLineColor(kRed+1);
    graStatV2JpsiFwdCentr020 -> SetLineWidth(2);
    graStatV2JpsiFwdCentr020 -> SetMarkerStyle(20);
    graStatV2JpsiFwdCentr020 -> SetMarkerSize(1.5);
    
    TGraphErrors *graSystV2JpsiFwdCentr020 = new TGraphErrors(nPtBins, ptJpsiFwdCentrBins, v2JpsiFwdCentr020Vals, ptJpsiFwdSystBins, v2JpsiFwdCentr020Systs);
    graSystV2JpsiFwdCentr020 -> SetMarkerColor(kRed+1);
    graSystV2JpsiFwdCentr020 -> SetMarkerStyle(20);
    graSystV2JpsiFwdCentr020 -> SetFillStyle(0);
    graSystV2JpsiFwdCentr020 -> SetLineColor(kRed+1);
    graSystV2JpsiFwdCentr020 -> SetLineWidth(2);
    graSystV2JpsiFwdCentr020 -> SetMarkerSize(1.5);

    const int nLargePtBins = 4;
    double minLargePtBins[] = {0, 2, 4, 6};
    double maxLargePtBins[] = {2, 4, 6, 8};
    double ptJpsiFwdCentrLargeBins[] = {1, 3, 5, 7};
    double ptJpsiFwdWidthLargeBins[] = {1, 1, 1, 1};
    double ptJpsiFwdSystLargeBins[] = {0.15, 0.15, 0.15, 0.15};
    
    double v2JpsiFwdCentr020LargeBinsVals[] = {0.0144323, 0.0329074, 0.0698103, 0.0844539};
    double v2JpsiFwdCentr020LargeBinsStats[] = {0.0180728, 0.0188962, 0.0303064, 0.049283};
    double v2JpsiFwdCentr020LargeBinsSysts[] = {0.00312387, 0.00928278, 0.0194577, 0.0269159};

    TGraphErrors *graStatV2JpsiFwdCentr020LargeBins = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr020LargeBinsVals, ptJpsiFwdWidthLargeBins, v2JpsiFwdCentr020LargeBinsStats);
    graStatV2JpsiFwdCentr020LargeBins -> SetMarkerColor(kOrange+7);
    graStatV2JpsiFwdCentr020LargeBins -> SetLineColor(kOrange+7);
    graStatV2JpsiFwdCentr020LargeBins -> SetLineWidth(2);
    graStatV2JpsiFwdCentr020LargeBins -> SetMarkerStyle(20);
    graStatV2JpsiFwdCentr020LargeBins -> SetMarkerSize(1.5);
    
    TGraphErrors *graSystV2JpsiFwdCentr020LargeBins = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr020LargeBinsVals, ptJpsiFwdSystLargeBins, v2JpsiFwdCentr020LargeBinsSysts);
    graSystV2JpsiFwdCentr020LargeBins -> SetMarkerColor(kOrange+7);
    graSystV2JpsiFwdCentr020LargeBins -> SetMarkerStyle(20);
    graSystV2JpsiFwdCentr020LargeBins -> SetFillStyle(0);
    graSystV2JpsiFwdCentr020LargeBins -> SetLineColor(kOrange+7);
    graSystV2JpsiFwdCentr020LargeBins -> SetLineWidth(2);
    graSystV2JpsiFwdCentr020LargeBins -> SetMarkerSize(1.5);

    // ***************************************************************************************** //
    // Centrality 0-20%
    // ***************************************************************************************** //
    // TPC-POS only Q-vectors
    //double v2JpsiFwdCentr2060Vals[] = {-0.0162254, -0.0178981, 0.00161603, 0.0205455, 0.0736224, 0.104662};
    //double v2JpsiFwdCentr2060Stats[] = {0.0388019, 0.0276514, 0.0301216, 0.0367978, 0.0371583, 0.066764};
    //double v2JpsiFwdCentr2060Systs[] = {0.0106743, 0.0289646, 0.00971417, 0.0148252, 0.00970835, 0.0282096};
    // TPC-ALL only Q-vectors
    double v2JpsiFwdCentr2060Vals[] = {-0.017018, -0.0210286, 0.000239199, 0.0185226, 0.0418244, 0.12257};
    double v2JpsiFwdCentr2060Stats[] = {0.0370604, 0.0263352, 0.0286848, 0.0347268, 0.0350661, 0.0639947};
    double v2JpsiFwdCentr2060Systs[] = {0.0145221, 0.0218097, 0.0039893, 0.00969075, 0.00867721, 0.0234387};

    Printf("J/psi v2 vs pT in centrality 20-60");
    for (int iPt = 0;iPt < nPtBins;iPt++) {
        Printf("%1.0f - %1.0f & %4.3f $?pm$ %4.3f $?pm$ %4.3f ??", minPtBins[iPt], maxPtBins[iPt], v2JpsiFwdCentr2060Vals[iPt], v2JpsiFwdCentr2060Stats[iPt], v2JpsiFwdCentr2060Systs[iPt]);
    }

    TGraphErrors *graStatV2JpsiFwdCentr2060 = new TGraphErrors(6, ptJpsiFwdCentrBins, v2JpsiFwdCentr2060Vals, ptJpsiFwdWidthBins, v2JpsiFwdCentr2060Stats);
    graStatV2JpsiFwdCentr2060 -> SetMarkerColor(kBlue+1);
    graStatV2JpsiFwdCentr2060 -> SetLineColor(kBlue+1);
    graStatV2JpsiFwdCentr2060 -> SetLineWidth(2);
    graStatV2JpsiFwdCentr2060 -> SetMarkerStyle(20);
    graStatV2JpsiFwdCentr2060 -> SetMarkerSize(1.5);
    
    TGraphErrors *graSystV2JpsiFwdCentr2060 = new TGraphErrors(6, ptJpsiFwdCentrBins, v2JpsiFwdCentr2060Vals, ptJpsiFwdSystBins, v2JpsiFwdCentr2060Systs);
    graSystV2JpsiFwdCentr2060 -> SetMarkerColor(kBlue+1);
    graSystV2JpsiFwdCentr2060 -> SetMarkerStyle(20);
    graSystV2JpsiFwdCentr2060 -> SetFillStyle(0);
    graSystV2JpsiFwdCentr2060 -> SetLineColor(kBlue+1);
    graSystV2JpsiFwdCentr2060 -> SetLineWidth(2);
    graSystV2JpsiFwdCentr2060 -> SetMarkerSize(1.5);


    double v2JpsiFwdCentr2060LargeBinsVals[] = {-0.00894195, 0.00854762, 0.0418244, 0.12257};
    double v2JpsiFwdCentr2060LargeBinsStats[] = {0.0214885, 0.0221865, 0.0350661, 0.0639947};
    double v2JpsiFwdCentr2060LargeBinsSysts[] = {0.00690191, 0.00723872, 0.00867721, 0.0234387};

    TGraphErrors *graStatV2JpsiFwdCentr2060LargeBins = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr2060LargeBinsVals, ptJpsiFwdWidthLargeBins, v2JpsiFwdCentr2060LargeBinsStats);
    graStatV2JpsiFwdCentr2060LargeBins -> SetMarkerColor(kAzure+2);
    graStatV2JpsiFwdCentr2060LargeBins -> SetLineColor(kAzure+2);
    graStatV2JpsiFwdCentr2060LargeBins -> SetLineWidth(2);
    graStatV2JpsiFwdCentr2060LargeBins -> SetMarkerStyle(20);
    graStatV2JpsiFwdCentr2060LargeBins -> SetMarkerSize(1.5);
    
    TGraphErrors *graSystV2JpsiFwdCentr2060LargeBins = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr2060LargeBinsVals, ptJpsiFwdSystLargeBins, v2JpsiFwdCentr2060LargeBinsSysts);
    graSystV2JpsiFwdCentr2060LargeBins -> SetMarkerColor(kAzure+2);
    graSystV2JpsiFwdCentr2060LargeBins -> SetMarkerStyle(20);
    graSystV2JpsiFwdCentr2060LargeBins -> SetFillStyle(0);
    graSystV2JpsiFwdCentr2060LargeBins -> SetLineColor(kAzure+2);
    graSystV2JpsiFwdCentr2060LargeBins -> SetLineWidth(2);
    graSystV2JpsiFwdCentr2060LargeBins -> SetMarkerSize(1.5);



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

    TGraphErrors *graStatV2JpsiMidCentr020 = new TGraphErrors(3, ptJpsiMidCentrBins, v2JpsiMidCentr020Vals, ptJpsiMidWidthBins, v2JpsiMidCentr020Stats);
    graStatV2JpsiMidCentr020 -> SetMarkerColor(kGreen+2);
    graStatV2JpsiMidCentr020 -> SetLineColor(kGreen+2);
    graStatV2JpsiMidCentr020 -> SetLineWidth(2);
    graStatV2JpsiMidCentr020 -> SetMarkerStyle(20);
    graStatV2JpsiMidCentr020 -> SetMarkerSize(1.5);
    
    TGraphErrors *graSystV2JpsiMidCentr020 = new TGraphErrors(3, ptJpsiMidCentrBins, v2JpsiMidCentr020Vals, ptJpsiMidSystBins, v2JpsiMidCentr020Systs);    
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

    // Jpsi Fwd 0-20%
    TCanvas *canvasJpsiFwdCentr020 = new TCanvas("canvasJpsiFwdCentr020", "", 800, 600);
    TH2D *histGridJpsiFwdCentr020  = new TH2D("histGridJpsiFwdCentr020", ";#it{p}_{T} (GeV/#it{c});#it{v}_{2}", 100, 0, 8, 100, -0.10, 0.20);
    histGridJpsiFwdCentr020 -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiFwdCentr020 -> Draw("E2P");
    graStatV2JpsiFwdCentr020 -> Draw("P SAME");
    graSystV2JpsiFwdCentr020LargeBins -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr020LargeBins -> Draw("P SAME");
    latexTitle -> DrawLatex(0.20, 0.88, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, 0#minus20%");

    // Jpsi Fwd 20-60%
    TCanvas *canvasJpsiFwdCentr2060 = new TCanvas("canvasJpsiFwdCentr2060", "", 800, 600);
    TH2D *histGridJpsiFwdCentr2060  = new TH2D("histGridJpsiFwdCentr2060", ";#it{p}_{T} (GeV/#it{c});#it{v}_{2}", 100, 0, 8, 100, -0.10, 0.20);
    histGridJpsiFwdCentr2060 -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiFwdCentr2060 -> Draw("E2P");
    graStatV2JpsiFwdCentr2060 -> Draw("P SAME");
    graSystV2JpsiFwdCentr2060LargeBins -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr2060LargeBins -> Draw("P SAME");
    latexTitle -> DrawLatex(0.20, 0.88, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, 20#minus60%");

    return;

    // Jpsi Fwd 0-20 vs 20-60%
    TCanvas *canvasJpsiFwd = new TCanvas("canvasJpsiFwd", "", 800, 600);
    TH2D *histGridJpsiFwd  = new TH2D("histGridJpsiFwd", ";#it{p}_{T} (GeV/#it{c});#it{v}_{2}", 100, 0, 8, 100, -0.10, 0.20);
    histGridJpsiFwd -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiFwdCentr020 -> Draw("E2P");
    graStatV2JpsiFwdCentr020 -> Draw("P SAME");
    graSystV2JpsiFwdCentr2060 -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr2060 -> Draw("P SAME");
    latexTitle -> DrawLatex(0.20, 0.88, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, 20#minus60%");

    // Jpsi Fwd vs theory
    TCanvas *canvasJpsiFwdVsTheor = new TCanvas("canvasJpsiFwdVsTheor", "", 800, 600);
    TH2D *histGridJpsiFwdVsTheor  = new TH2D("histGridJpsiFwdVsTheor", ";#it{p}_{T} (GeV/#it{c});#it{v}_{2}", 100, 0, 8, 100, -0.30, 0.30);
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
    TH2D *histGridJpsiFwdVsJpsiMid  = new TH2D("histGridJpsiFwdVsJpsiMid", ";#it{p}_{T} (GeV/#it{c});#it{v}_{2}", 100, 0, 8, 100, -0.30, 0.30);
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
    TH2D *histGridJpsiFwdVsDzeroMid  = new TH2D("histGridJpsiFwdVsDzeroMid", ";#it{p}_{T} (GeV/#it{c});#it{v}_{2}", 100, 0, 8, 100, -0.30, 0.30);
    histGridJpsiFwdVsDzeroMid -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiFwdCentr020 -> Draw("E2P");
    graStatV2JpsiFwdCentr020 -> Draw("P");
    graSyst1V2DzeroMidCentr020 -> Draw("E2P");
    //graSyst2V2DzeroMidCentr020 -> Draw("E2P");
    graStatV2DzeroMidCentr020 -> Draw("EP");

    TLegend *legend3 = new TLegend(0.18,0.20,0.48,0.35);
    SetLegend(legend3);
    legend3 -> SetTextSize(0.045);
    legend3 -> AddEntry(graSystV2JpsiFwdCentr020,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4, 0#minus20%","P");
    legend3 -> AddEntry(graStatV2DzeroMidCentr020,"D^{0} #rightarrow K^{-}#pi^{+} and charge conj., |#it{y}| < 0.8, 0#minus20%","P");
    legend3 -> Draw();

    latexTitle -> DrawLatex(0.20, 0.88, "ALICE Work in Progress, OO, #sqrt{#it{s}_{NN}} = 5.36 TeV");

    canvasJpsiFwd -> SaveAs("figures/jpsiFwd.pdf");
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