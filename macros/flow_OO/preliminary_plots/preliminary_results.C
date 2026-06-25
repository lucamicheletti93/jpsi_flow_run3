void LoadStyle();
void SetLegend(TLegend *);
void SetGraph(TGraphErrors *, Color_t , double , int , double , bool );
void SetAsymmGraph(TGraphAsymmErrors *, Color_t , double , int , double , bool );
TGraphAsymmErrors* DoGraphFromTheory(string);

void preliminary_results() {
    LoadStyle();
    // ***************************************************************************************** //
    // Theory predictions from THU
    // ***************************************************************************************** //
    string fInNameCentr020 = "/Users/lucamicheletti/GITHUB/dq_run3_analyses/charmonia_production_pO_OO_NeNe/theory/THU/forward/data_v2_jpsi_020_y254.dat";
    TGraphAsymmErrors *graTheorFwdCentr020CentrVal = DoGraphFromTheory(fInNameCentr020);
    gStyle->SetLineStyleString(9,"80 20");
    graTheorFwdCentr020CentrVal->SetLineColor(kOrange+7);
    graTheorFwdCentr020CentrVal->SetLineWidth(2);
    graTheorFwdCentr020CentrVal->SetLineStyle(9);

    string fInNameCentr2060 = "/Users/lucamicheletti/GITHUB/dq_run3_analyses/charmonia_production_pO_OO_NeNe/theory/THU/forward/data_v2_jpsi_2060_y254.dat";
    TGraphAsymmErrors *graTheorFwdCentr2060CentrVal = DoGraphFromTheory(fInNameCentr2060);
    gStyle->SetLineStyleString(9,"80 20");
    graTheorFwdCentr2060CentrVal->SetLineColor(kOrange+7);
    graTheorFwdCentr2060CentrVal->SetLineWidth(2);
    graTheorFwdCentr2060CentrVal->SetLineStyle(9);

    // ***************************************************************************************** //
    // LF results
    // ***************************************************************************************** //
    TFile *fInLambda = new TFile("v2Lambda.root", "READ");
    TGraphErrors *graStatV2LambdaMidCentr010 = (TGraphErrors*) fInLambda->Get("gist_lambda_010");
    SetGraph(graStatV2LambdaMidCentr010, kGreen+2, 1.5, 20, 1, false);

    TGraphErrors *graSyst1V2LambdaMidCentr010 = (TGraphErrors*) fInLambda->Get("syst_lambda_010");
    SetGraph(graSyst1V2LambdaMidCentr010, kGreen+2, 1.5, 20, 1, true);

    TFile *fInKzero = new TFile("v2K0.root", "READ");
    TGraphErrors *graStatV2KzeroMidCentr010 = (TGraphErrors*) fInKzero->Get("gist_k0_010");
    SetGraph(graStatV2KzeroMidCentr010, kOrange+7, 1.5, 20, 1, false);

    TGraphErrors *graSyst1V2KzeroMidCentr010 = (TGraphErrors*) fInKzero->Get("syst_k0_010");
    SetGraph(graSyst1V2KzeroMidCentr010, kOrange+7, 1.5, 20, 1, true);
    // ***************************************************************************************** //
    // D-meson results
    // ***************************************************************************************** //
    TFile *fInDzero = new TFile("v2Dzero.root", "READ");
    TGraphAsymmErrors *graStatV2DzeroMidCentr020 = (TGraphAsymmErrors*) fInDzero->Get("gvn_prompt_stat");
    SetAsymmGraph(graStatV2DzeroMidCentr020, kAzure+4, 1.5, 20, 1, false);

    TGraphAsymmErrors *graSyst1V2DzeroMidCentr020 = (TGraphAsymmErrors*) fInDzero->Get("tot_syst");
    SetAsymmGraph(graSyst1V2DzeroMidCentr020, kAzure+4, 1.5, 20, 1, true);

    for (int iPoint = 0;iPoint < 11;iPoint++) {
        graSyst1V2DzeroMidCentr020->SetPointEXhigh(iPoint, 0.15);
        graSyst1V2DzeroMidCentr020->SetPointEXlow(iPoint, 0.15);
    }

    TGraphAsymmErrors *graSyst2V2DzeroMidCentr020 = (TGraphAsymmErrors*) fInDzero->Get("tot_syst_wnonflow");
    SetAsymmGraph(graSyst2V2DzeroMidCentr020, kAzure+1, 1.5, 20, 1, true);

    // ***************************************************************************************** //
    // J/psi results
    // ***************************************************************************************** //
    // Pb-Pb collisions
    const int nPtBinsPbPb = 13;
    double minPtBinsPbPb[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0};
    double maxPtBinsPbPb[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0};
    double ptJpsiFwdCentrBinsPbPb[] = {0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.5, 4.5, 5.5, 7.0, 9.0, 11.0, 13.5};
    double ptJpsiFwdWidthBinsPbPb[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double ptJpsiFwdSystBinsPbPb[] = {0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125};

    double v2JpsiPbPbFwdCentr010Vals[] = {0.0228447, 0.012234, 0.00634267, 0.0157593, 0.0355709, 0.0336122, 0.0614967, 0.0387399, 0.0310459, 0.0278144, 0.0358532, 0.000178693, -0.0188813};
    double v2JpsiPbPbFwdCentr010Stats[] = {0.0161063, 0.0142315, 0.00820513, 0.0123487, 0.00994421, 0.0108773, 0.0093194, 0.0124067, 0.0159081, 0.0160622, 0.0257618, 0.040254, 0.0496009};
    double v2JpsiPbPbFwdCentr010Systs[] = {0.00794097, 0.00727943, 0.00279877, 0.0022176, 0.00722703, 0.00257224, 0.00135059, 0.0050472, 0.00415074, 0.00587882, 0.00474813, 0.00980131, 0.00566724};

    TGraphErrors *graStatV2JpsiPbPbFwdCentr010 = new TGraphErrors(nPtBinsPbPb, ptJpsiFwdCentrBinsPbPb, v2JpsiPbPbFwdCentr010Vals, ptJpsiFwdWidthBinsPbPb, v2JpsiPbPbFwdCentr010Stats);
    SetGraph(graStatV2JpsiPbPbFwdCentr010, kGray+1, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiPbPbFwdCentr010 = new TGraphErrors(nPtBinsPbPb, ptJpsiFwdCentrBinsPbPb, v2JpsiPbPbFwdCentr010Vals, ptJpsiFwdSystBinsPbPb, v2JpsiPbPbFwdCentr010Systs);
    SetGraph(graSystV2JpsiPbPbFwdCentr010, kGray+1, 1.5, 20, 1, true);

    double v2JpsiPbPbFwdCentr3050Vals[] = {-0.00751166, 0.0324682, 0.0483128, 0.0623836, 0.0790007, 0.0750308, 0.0821814, 0.0844028, 0.110533, 0.0894459, 0.0864354, 0.0662619, 0.105772};
    double v2JpsiPbPbFwdCentr3050Stats[] = {0.0138548, 0.0119431, 0.00818823, 0.010507, 0.00958244, 0.00985444, 0.00793826, 0.00942195, 0.0112991, 0.011538, 0.019468, 0.0290814, 0.0375577};
    double v2JpsiPbPbFwdCentr3050Systs[] = {0.00456438, 0.00532624, 0.00657618, 0.00767329, 0.00376474, 0.00244701, 0.00190689, 0.00266505, 0.00643058, 0.00989524, 0.00819967, 0.0119924, 0.012523};

    TGraphErrors *graStatV2JpsiPbPbFwdCentr3050 = new TGraphErrors(nPtBinsPbPb, ptJpsiFwdCentrBinsPbPb, v2JpsiPbPbFwdCentr3050Vals, ptJpsiFwdWidthBinsPbPb, v2JpsiPbPbFwdCentr3050Stats);
    SetGraph(graStatV2JpsiPbPbFwdCentr3050, kGray+1, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiPbPbFwdCentr3050 = new TGraphErrors(nPtBinsPbPb, ptJpsiFwdCentrBinsPbPb, v2JpsiPbPbFwdCentr3050Vals, ptJpsiFwdSystBinsPbPb, v2JpsiPbPbFwdCentr3050Systs);
    SetGraph(graSystV2JpsiPbPbFwdCentr3050, kGray+1, 1.5, 20, 1, true);

    // p-Pb & Pb-p collisions
    const int nPtBinspPb = 5;
    double minPtBinspPb[] = {0.0, 2.0, 3.0, 4.0, 6.0};
    double maxPtBinspPb[] = {2.0, 3.0, 4.0, 6.0, 8.0};
    double ptJpsiFwdCentrBinspPb[] = {1.24, 2.48, 3.48, 4.87, 6.83};
    double ptJpsiFwdWidthBinspPb[] = {0, 0, 0, 0, 0};
    double ptJpsiFwdSystBinspPb[] = {0.125, 0.125, 0.125, 0.125, 0.125};

    double v2JpsipPbFwdVals[] = {-0.014, 0.004, 0.038, 0.092, 0.033};
    double v2JpsipPbFwdStats[] = {0.023, 0.026, 0.029, 0.026, 0.039};
    double v2JpsipPbFwdSysts1Minus[] = {-0.015, -0.016, -0.016, -0.015, -0.026};
    double v2JpsipPbFwdSysts1Plus[] = {0.015, 0.016, 0.016, 0.013, 0.019};
    double v2JpsipPbFwdSysts2[] = {-0.001, 0.000, 0.002, 0.005, 0.002};
    double v2JpsipPbFwdSystsMin[nPtBinspPb], v2JpsipPbFwdSystsMax[nPtBinspPb];

    for (int iPt = 0;iPt < nPtBinspPb;iPt++) {
        v2JpsipPbFwdSystsMin[iPt] = TMath::Sqrt(v2JpsipPbFwdSysts1Minus[iPt]*v2JpsipPbFwdSysts1Minus[iPt] + v2JpsipPbFwdSysts2[iPt]*v2JpsipPbFwdSysts2[iPt]);
        v2JpsipPbFwdSystsMax[iPt] = TMath::Sqrt(v2JpsipPbFwdSysts1Plus[iPt]*v2JpsipPbFwdSysts1Plus[iPt] + v2JpsipPbFwdSysts2[iPt]*v2JpsipPbFwdSysts2[iPt]);
    }

    TGraphErrors *graStatV2JpsipPbFwd = new TGraphErrors(nPtBinspPb, ptJpsiFwdCentrBinspPb, v2JpsipPbFwdVals, ptJpsiFwdWidthBinspPb, v2JpsipPbFwdStats);
    SetGraph(graStatV2JpsipPbFwd, kGreen+2, 1.5, 20, 1, false);

    TGraphAsymmErrors *graSystV2JpsipPbFwd = new TGraphAsymmErrors(nPtBinspPb, ptJpsiFwdCentrBinspPb, v2JpsipPbFwdVals, ptJpsiFwdSystBinspPb, ptJpsiFwdSystBinspPb, v2JpsipPbFwdSystsMin, v2JpsipPbFwdSystsMax);
    SetAsymmGraph(graSystV2JpsipPbFwd, kGreen+2, 1.5, 20, 1, true);

    const int nPtBinsPbp = 5;
    double minPtBinsPbp[] = {0.0, 2.0, 3.0, 4.0, 6.0};
    double maxPtBinsPbp[] = {2.0, 3.0, 4.0, 6.0, 8.0};
    double ptJpsiFwdCentrBinsPbp[] = {1.23, 2.48, 3.46, 4.84, 6.81};
    double ptJpsiFwdWidthBinsPbp[] = {0, 0, 0, 0, 0};
    double ptJpsiFwdSystBinsPbp[] = {0.125, 0.125, 0.125, 0.125, 0.125};

    double v2JpsiPbpFwdVals[] = {0.026, 0.013, 0.058, 0.062, 0.083};
    double v2JpsiPbpFwdStats[] = {0.017, 0.020, 0.024, 0.022, 0.037};
    double v2JpsiPbpFwdSysts1Minus[] = {-0.017, -0.017, -0.017, -0.015, -0.022};
    double v2JpsiPbpFwdSysts1Plus[] = {0.017, 0.017, 0.017, 0.015, 0.022};
    double v2JpsiPbpFwdSysts2[] = {0.001, 0.001, 0.003, 0.003, 0.004};
    double v2JpsiPbpFwdSystsMin[nPtBinsPbp], v2JpsiPbpFwdSystsMax[nPtBinsPbp];

    for (int iPt = 0;iPt < nPtBinsPbp;iPt++) {
        v2JpsiPbpFwdSystsMin[iPt] = TMath::Sqrt(v2JpsiPbpFwdSysts1Minus[iPt]*v2JpsiPbpFwdSysts1Minus[iPt] + v2JpsiPbpFwdSysts2[iPt]*v2JpsiPbpFwdSysts2[iPt]);
        v2JpsiPbpFwdSystsMax[iPt] = TMath::Sqrt(v2JpsiPbpFwdSysts1Plus[iPt]*v2JpsiPbpFwdSysts1Plus[iPt] + v2JpsiPbpFwdSysts2[iPt]*v2JpsiPbpFwdSysts2[iPt]);
    }

    TGraphErrors *graStatV2JpsiPbpFwd = new TGraphErrors(nPtBinsPbp, ptJpsiFwdCentrBinsPbp, v2JpsiPbpFwdVals, ptJpsiFwdWidthBinsPbp, v2JpsiPbpFwdStats);
    SetGraph(graStatV2JpsiPbpFwd, kAzure+2, 1.5, 20, 1, false);

    TGraphAsymmErrors *graSystV2JpsiPbpFwd = new TGraphAsymmErrors(nPtBinsPbp, ptJpsiFwdCentrBinsPbp, v2JpsiPbpFwdVals, ptJpsiFwdSystBinsPbp, ptJpsiFwdSystBinsPbp, v2JpsiPbpFwdSystsMin, v2JpsiPbpFwdSystsMax);
    SetAsymmGraph(graSystV2JpsiPbpFwd, kAzure+2, 1.5, 20, 1, true);

    // OO collisions: midrapidity
    double ptJpsiMidWidth[] = {0. , 0., 0.};
    double ptJpsiMidSyst[] = {0.15, 0.15, 0.15};
    
    double ptJpsiOOMidCentr020[] = {1.230, 2.870, 5.300};
    double v2JpsiMidCentr020Vals[] = {-0.180000, -0.018714, 0.059150};
    double v2JpsiMidCentr020Stats[] = {0.102775, 0.102947, 0.058341};
    double v2JpsiMidCentr020Systs[] = {0.050000, 0.090000, 0.060000};

    TGraphErrors *graStatV2JpsiOOMidCentr020 = new TGraphErrors(3, ptJpsiOOMidCentr020, v2JpsiMidCentr020Vals, ptJpsiMidWidth, v2JpsiMidCentr020Stats);
    SetGraph(graStatV2JpsiOOMidCentr020, kAzure+2, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiOOMidCentr020 = new TGraphErrors(3, ptJpsiOOMidCentr020, v2JpsiMidCentr020Vals, ptJpsiMidSyst, v2JpsiMidCentr020Systs);
    SetGraph(graSystV2JpsiOOMidCentr020, kAzure+2, 1.5, 20, 1, true);

    double ptJpsiOOMidCentr2060[] = {1.190, 2.920 , 5.270};
    double v2JpsiMidCentr2060Vals[] = {-0.032969, 0.070266, 0.145450};
    double v2JpsiMidCentr2060Stats[] = {0.097961, 0.079300, 0.067461};
    double v2JpsiMidCentr2060Systs[] = {0.070000, 0.080000, 0.060000};

    TGraphErrors *graStatV2JpsiOOMidCentr2060 = new TGraphErrors(3, ptJpsiOOMidCentr2060, v2JpsiMidCentr2060Vals, ptJpsiMidWidth, v2JpsiMidCentr2060Stats);
    SetGraph(graStatV2JpsiOOMidCentr2060, kAzure+2, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiOOMidCentr2060 = new TGraphErrors(3, ptJpsiOOMidCentr2060, v2JpsiMidCentr2060Vals, ptJpsiMidSyst, v2JpsiMidCentr2060Systs);
    SetGraph(graSystV2JpsiOOMidCentr2060, kAzure+2, 1.5, 20, 1, true);

    // OO collisions: forward rapidity
    const int nPtBins = 4;
    double minPtBins[] = {0, 2, 4, 6};
    double maxPtBins[] = {2, 4, 6, 8};
    double ptJpsiFwdWidth[] = {0, 0, 0, 0};
    double ptJpsiFwdSyst[] = {0.15, 0.15, 0.15, 0.15};

    double ptJpsiOOFwdCentr020[] = {1.28, 2.87, 4.82, 6.84};
    double v2JpsiOOFwdCentr020Vals[] = {0.0144323, 0.0329074, 0.0698103, 0.0844539};
    double v2JpsiOOFwdCentr020Stats[] = {0.0180728, 0.0188962, 0.0303064, 0.049283};
    double v2JpsiOOFwdCentr020Systs[] = {0.00312387, 0.00928278, 0.0194577, 0.0269159};

    TGraphErrors *graStatV2JpsiOOFwdCentr020 = new TGraphErrors(nPtBins, ptJpsiOOFwdCentr020, v2JpsiOOFwdCentr020Vals, ptJpsiFwdWidth, v2JpsiOOFwdCentr020Stats);
    SetGraph(graStatV2JpsiOOFwdCentr020, kRed+1, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiOOFwdCentr020 = new TGraphErrors(nPtBins, ptJpsiOOFwdCentr020, v2JpsiOOFwdCentr020Vals, ptJpsiFwdSyst, v2JpsiOOFwdCentr020Systs);
    SetGraph(graSystV2JpsiOOFwdCentr020, kRed+1, 1.5, 20, 1, true);

    double ptJpsiOOFwdCentr2060[] = {1.24, 2.87, 4.83, 6.76};
    double v2JpsiOOFwdCentr2060Vals[] = {-0.00894195, 0.00854762, 0.0418244, 0.12257};
    double v2JpsiOOFwdCentr2060Stats[] = {0.0214885, 0.0221865, 0.0350661, 0.0639947};
    double v2JpsiOOFwdCentr2060Systs[] = {0.00690191, 0.00723872, 0.00867721, 0.0234387};

    TGraphErrors *graStatV2JpsiOOFwdCentr2060 = new TGraphErrors(nPtBins, ptJpsiOOFwdCentr2060, v2JpsiOOFwdCentr2060Vals, ptJpsiFwdWidth, v2JpsiOOFwdCentr2060Stats);
    SetGraph(graStatV2JpsiOOFwdCentr2060, kRed+1, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiOOFwdCentr2060 = new TGraphErrors(nPtBins, ptJpsiOOFwdCentr2060, v2JpsiOOFwdCentr2060Vals, ptJpsiFwdSyst, v2JpsiOOFwdCentr2060Systs);
    SetGraph(graSystV2JpsiOOFwdCentr2060, kRed+1, 1.5, 20, 1, true);

    TGraphErrors *graStatV2JpsiOOFwdCentr2060_clone = (TGraphErrors*) graStatV2JpsiOOFwdCentr2060 -> Clone("graStatV2JpsiOOFwdCentr2060_clone");
    SetGraph(graStatV2JpsiOOFwdCentr2060_clone, kAzure+2, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiOOFwdCentr2060_clone = (TGraphErrors*) graSystV2JpsiOOFwdCentr2060 -> Clone("graSystV2JpsiOOFwdCentr2060_clone");
    SetGraph(graSystV2JpsiOOFwdCentr2060_clone, kAzure+2, 1.5, 20, 1, true);

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

    TH2D *histGridV2JpsiOO = new TH2D("histGridV2JpsiOO", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}^{SP} {|#Delta#it{#eta}| > 1.7}", 100, 0, 8, 100, -0.05, 0.20);
    TH2D *histGridV2JpsiAllSystems = new TH2D("histGridV2JpsiAllSystems", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}^{J/#psi}", 100, 0, 8, 100, -0.05, 0.30);
    TH2D *histGridV2JpsiAllSpecies = new TH2D("histGridV2JpsiAllSpecies", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}", 100, 0, 8, 100, -0.03, 0.37);
    TH2D *histGridV2JpsiOOvsPbPb = new TH2D("histGridV2JpsiOOvsPbPb", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}^{SP} {|#Delta#it{#eta}| > 1.7}", 100, 0, 8, 100, -0.05, 0.20);
    TH2D *histGridV2JpsiOOVsTheor = new TH2D("histGridV2JpsiOOVsTheor", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}^{J/#psi}", 100, 0, 8, 100, -0.05, 0.20);
    TH2D *histGridV2JpsiFwdVsMid = new TH2D("histGridV2JpsiFwdVsMid", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}^{J/#psi}", 100, 0, 8, 100, -0.30, 0.30);

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // J/psi: Oxygen-Oxygen
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    TCanvas *canvasV2JpsiOO = new TCanvas("canvasV2JpsiOO", "", 800, 600);
    histGridV2JpsiOO->Draw();
    lineUnity->Draw("SAME");
    graSystV2JpsiOOFwdCentr020->Draw("E2P SAME");
    graStatV2JpsiOOFwdCentr020->Draw("P SAME");
    graSystV2JpsiOOFwdCentr2060_clone->Draw("E2P SAME");
    graStatV2JpsiOOFwdCentr2060_clone->Draw("P SAME");

    latexTitle->DrawLatex(0.20, 0.88, "ALICE Preliminary");
    latexTitle->DrawLatex(0.20, 0.82, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");

    TLegend *legendV2JpsiOO = new TLegend(0.20,0.65,0.40,0.78);
    SetLegend(legendV2JpsiOO);
    legendV2JpsiOO->AddEntry(graSystV2JpsiOOFwdCentr020,"0#minus20%","P");
    legendV2JpsiOO->AddEntry(graSystV2JpsiOOFwdCentr2060_clone,"20#minus60%","P");
    legendV2JpsiOO->Draw();


    TCanvas *canvasV2JpsiOOCentral = new TCanvas("canvasV2JpsiOOCentral", "", 800, 600);
    histGridV2JpsiOO->Draw();
    lineUnity->Draw("SAME");
    graSystV2JpsiOOFwdCentr020->Draw("E2P SAME");
    graStatV2JpsiOOFwdCentr020->Draw("P SAME");

    latexTitle->DrawLatex(0.20, 0.88, "ALICE Preliminary, OO, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle->DrawLatex(0.20, 0.82, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4, 0#minus20%");


    TCanvas *canvasV2JpsiOOSemientral = new TCanvas("canvasV2JpsiOOSemientral", "", 800, 600);
    histGridV2JpsiOO->Draw();
    lineUnity->Draw("SAME");
    graSystV2JpsiOOFwdCentr2060->Draw("E2P SAME");
    graStatV2JpsiOOFwdCentr2060->Draw("P SAME");

    latexTitle->DrawLatex(0.20, 0.88, "ALICE Preliminary, OO, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle->DrawLatex(0.20, 0.82, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4, 20#minus60%");

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // J/psi: all collision systems comparison
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    TCanvas *canvasV2JpsiAllSystems = new TCanvas("canvasV2JpsiAllSystems", "", 800, 600);
    histGridV2JpsiAllSystems->Draw();

    TLegend *legendV2Jpsi1 = new TLegend(0.20,0.72,0.80,0.82);
    SetLegend(legendV2Jpsi1);
    legendV2Jpsi1->SetNColumns(2);
    legendV2Jpsi1->AddEntry(graSystV2JpsipPbFwd,"2.03 < #it{y} < 3.53","P");
    legendV2Jpsi1->AddEntry(graSystV2JpsiPbpFwd,"#minus4.46 < #it{y} < #minus2.96","P");
    legendV2Jpsi1->Draw();

    TLegend *legendV2Jpsi2 = new TLegend(0.20,0.57,0.82,0.67);
    SetLegend(legendV2Jpsi2);
    legendV2Jpsi2->SetNColumns(2);
    legendV2Jpsi2->AddEntry(graSystV2JpsiOOFwdCentr2060,"OO, 20#minus60%","P");
    legendV2Jpsi2->AddEntry(graSystV2JpsiPbPbFwdCentr3050,"Pb#minusPb, 30#minus50%","P");
    legendV2Jpsi2->Draw();

    latexTitle->DrawLatex(0.20, 0.88, "ALICE Preliminary, J/#psi #rightarrow #mu^{+}#mu^{-}");
    latexTitle->DrawLatex(0.20, 0.82, "p#minusPb, #sqrt{#it{s}_{NN}} = 8.16 TeV, PLB 780 (2018)");
    latexTitle->DrawLatex(0.20, 0.67, "#sqrt{#it{s}_{NN}} = 5.36 TeV, 2.5 < #it{y} < 4");

    lineUnity->Draw("SAME");
    graSystV2JpsipPbFwd->Draw("E2P");
    graStatV2JpsipPbFwd->Draw("P SAME");
    graSystV2JpsiPbpFwd->Draw("E2P SAME");
    graStatV2JpsiPbpFwd->Draw("P SAME");
    graSystV2JpsiPbPbFwdCentr3050->Draw("E2P SAME");
    graStatV2JpsiPbPbFwdCentr3050->Draw("P SAME");
    graStatV2JpsiOOFwdCentr2060->Draw("P SAME");
    graSystV2JpsiOOFwdCentr2060->Draw("E2P SAME");

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // J/psi: OO vs Pb-Pb
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    TCanvas *canvasV2JpsiOOvsPbPbCentral = new TCanvas("canvasV2JpsiOOvsPbPbCentral", "", 800, 600);
    histGridV2JpsiOOvsPbPb->Draw();
    lineUnity->Draw("SAME");
    graSystV2JpsiPbPbFwdCentr010->Draw("E2P SAME");
    graStatV2JpsiPbPbFwdCentr010->Draw("P SAME");
    graSystV2JpsiOOFwdCentr020->Draw("E2P SAME");
    graStatV2JpsiOOFwdCentr020->Draw("P SAME");
    
    TLegend *legendV2JpsiOOvsPbPbCentral = new TLegend(0.20,0.65,0.40,0.78);
    SetLegend(legendV2JpsiOOvsPbPbCentral);
    legendV2JpsiOOvsPbPbCentral->AddEntry(graSystV2JpsiPbPbFwdCentr010,"Pb#minusPb, 0#minus10%","P");
    legendV2JpsiOOvsPbPbCentral->AddEntry(graSystV2JpsiOOFwdCentr020,"OO, 0#minus20%","P");
    legendV2JpsiOOvsPbPbCentral->Draw();

    latexTitle->DrawLatex(0.20, 0.88, "ALICE Preliminary");
    latexTitle->DrawLatex(0.20, 0.82, "#sqrt{#it{s}_{NN}} = 5.36 TeV, J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");


    TCanvas *canvasV2JpsiOOvsPbPbSemicentral = new TCanvas("canvasV2JpsiOOvsPbPbSemicentral", "", 800, 600);
    histGridV2JpsiOOvsPbPb->Draw();
    lineUnity->Draw("SAME");
    graSystV2JpsiPbPbFwdCentr3050->Draw("E2P SAME");
    graStatV2JpsiPbPbFwdCentr3050->Draw("P SAME");
    graSystV2JpsiOOFwdCentr2060->Draw("E2P SAME");
    graStatV2JpsiOOFwdCentr2060->Draw("P SAME");
    
    TLegend *legendV2JpsiOOvsPbPbSemicentral = new TLegend(0.20,0.65,0.40,0.78);
    SetLegend(legendV2JpsiOOvsPbPbSemicentral);
    legendV2JpsiOOvsPbPbSemicentral->AddEntry(graSystV2JpsiPbPbFwdCentr3050,"Pb#minusPb, 30#minus50%","P");
    legendV2JpsiOOvsPbPbSemicentral->AddEntry(graSystV2JpsiOOFwdCentr2060,"OO, 20#minus60%","P");
    legendV2JpsiOOvsPbPbSemicentral->Draw();

    latexTitle->DrawLatex(0.20, 0.88, "ALICE Preliminary");
    latexTitle->DrawLatex(0.20, 0.82, "#sqrt{#it{s}_{NN}} = 5.36 TeV, J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // J/psi: model comparison
    TCanvas *canvasV2JpsiOOVsTheorCentral = new TCanvas("canvasV2JpsiOOVsTheorCentral", "", 800, 600);
    histGridV2JpsiOOVsTheor->Draw();
    lineUnity->Draw("SAME");
    graTheorFwdCentr020CentrVal->Draw("L SAME");
    graSystV2JpsiOOFwdCentr020->Draw("E2P SAME");
    graStatV2JpsiOOFwdCentr020->Draw("P SAME");

    TLegend *legendV2JpsiOOVsTheorCentral = new TLegend(0.20,0.60,0.40,0.73);
    SetLegend(legendV2JpsiOOVsTheorCentral);
    legendV2JpsiOOVsTheorCentral->AddEntry(graTheorFwdCentr020CentrVal,"Tsinghua Transport + MUSIC Hydro.","L");
    legendV2JpsiOOVsTheorCentral->AddEntry(graStatV2JpsiOOFwdCentr020,"Data, SP, |#Delta#it{#eta}| > 1.7","P");
    legendV2JpsiOOVsTheorCentral->Draw();

    latexTitle->DrawLatex(0.20, 0.88, "ALICE Preliminary");
    latexTitle->DrawLatex(0.20, 0.82, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle->DrawLatex(0.20, 0.76, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4, 0#minus20%");


    TCanvas *canvasV2JpsiOOVsTheorSemicentral = new TCanvas("canvasV2JpsiOOVsTheorSemicentral", "", 800, 600);
    histGridV2JpsiOOVsTheor->Draw();
    lineUnity->Draw("SAME");
    graTheorFwdCentr2060CentrVal->Draw("L SAME");
    graSystV2JpsiOOFwdCentr2060->Draw("E2P SAME");
    graStatV2JpsiOOFwdCentr2060->Draw("P SAME");

    TLegend *legendV2JpsiOOVsTheorSemientral = new TLegend(0.20,0.60,0.40,0.73);
    SetLegend(legendV2JpsiOOVsTheorSemientral);
    legendV2JpsiOOVsTheorSemientral->AddEntry(graTheorFwdCentr020CentrVal,"Tsinghua Transport + MUSIC Hydro.","L");
    legendV2JpsiOOVsTheorSemientral->AddEntry(graStatV2JpsiOOFwdCentr020,"Data, SP, |#Delta#it{#eta}| > 1.7","P");
    legendV2JpsiOOVsTheorSemientral->Draw();

    latexTitle->DrawLatex(0.20, 0.88, "ALICE Preliminary");
    latexTitle->DrawLatex(0.20, 0.82, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle->DrawLatex(0.20, 0.76, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4, 20#minus60%");

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // J/psi: mid vs forward rapidity
    TCanvas *canvasV2JpsiFwdVsMidCentral = new TCanvas("canvasV2JpsiFwdVsMidCentral", "", 800, 600);
    histGridV2JpsiFwdVsMid->Draw();
    lineUnity->Draw("SAME");
    graSystV2JpsiOOMidCentr020->Draw("E2P SAME");
    graStatV2JpsiOOMidCentr020->Draw("P SAME");
    graSystV2JpsiOOFwdCentr020->Draw("E2P SAME");
    graStatV2JpsiOOFwdCentr020->Draw("P SAME");

    TLegend *legendV2JpsiFwdVsMidCentral = new TLegend(0.50,0.27,0.80,0.40);
    SetLegend(legendV2JpsiFwdVsMidCentral);
    legendV2JpsiFwdVsMidCentral->AddEntry(graSystV2JpsiOOMidCentr020,"J/#psi #rightarrow e^{+}e^{-}, |#it{y}| < 0.9","P");
    legendV2JpsiFwdVsMidCentral->AddEntry(graSystV2JpsiOOFwdCentr020,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4","P");
    legendV2JpsiFwdVsMidCentral->Draw();

    latexTitle->DrawLatex(0.20, 0.88, "ALICE Preliminary");
    latexTitle->DrawLatex(0.20, 0.82, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, 0#minus20%");


    TCanvas *canvasV2JpsiFwdVsMidSemicentral = new TCanvas("canvasV2JpsiFwdVsMidSemicentral", "", 800, 600);
    histGridV2JpsiFwdVsMid->Draw();
    lineUnity->Draw("SAME");
    graSystV2JpsiOOMidCentr2060->Draw("E2P SAME");
    graStatV2JpsiOOMidCentr2060->Draw("P SAME");
    graSystV2JpsiOOFwdCentr2060->Draw("E2P SAME");
    graStatV2JpsiOOFwdCentr2060->Draw("P SAME");

    TLegend *legendV2JpsiFwdVsMidSemicentral = new TLegend(0.50,0.27,0.80,0.40);
    SetLegend(legendV2JpsiFwdVsMidSemicentral);
    legendV2JpsiFwdVsMidSemicentral->AddEntry(graSystV2JpsiOOMidCentr2060,"J/#psi #rightarrow e^{+}e^{-}, |#it{y}| < 0.9","P");
    legendV2JpsiFwdVsMidSemicentral->AddEntry(graSystV2JpsiOOFwdCentr2060,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4","P");
    legendV2JpsiFwdVsMidSemicentral->Draw();

    latexTitle->DrawLatex(0.20, 0.88, "ALICE Preliminary");
    latexTitle->DrawLatex(0.20, 0.82, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, 20#minus60%");

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // J/psi: all collision systems comparison
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    TCanvas *canvasV2JpsiAllSpecies = new TCanvas("canvasV2JpsiAllSpecies", "", 800, 600);
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



    ////////////////////////////////////////////////////////////////////////////////////
    canvasV2JpsiOO->SaveAs("ICHEP2026/approved_plots/jpsi_v2_OO.pdf");
    canvasV2JpsiOOCentral->SaveAs("ICHEP2026/approved_plots/jpsi_v2_OO_central.pdf");
    canvasV2JpsiOOSemientral->SaveAs("ICHEP2026/approved_plots/jpsi_v2_OO_semicentral.pdf");
    canvasV2JpsiAllSystems->SaveAs("ICHEP2026/approved_plots/jpsi_v2_all_systems.pdf");
    canvasV2JpsiOOvsPbPbCentral->SaveAs("ICHEP2026/approved_plots/jpsi_v2_OO_vs_PbPb_central.pdf");
    canvasV2JpsiOOvsPbPbSemicentral->SaveAs("ICHEP2026/approved_plots/jpsi_v2_OO_vs_PbPb_semicentral.pdf");
    canvasV2JpsiOOVsTheorCentral->SaveAs("ICHEP2026/approved_plots/jpsi_v2_OO_vs_theory_central.pdf");
    canvasV2JpsiOOVsTheorSemicentral->SaveAs("ICHEP2026/approved_plots/jpsi_v2_OO_vs_theory_semicentral.pdf");
    canvasV2JpsiFwdVsMidCentral->SaveAs("ICHEP2026/approved_plots/jpsi_v2_fwd_vs_mid_central.pdf");
    canvasV2JpsiFwdVsMidSemicentral->SaveAs("ICHEP2026/approved_plots/jpsi_v2_fwd_vs_mid_semicentral.pdf");
    canvasV2JpsiAllSpecies->SaveAs("ICHEP2026/approved_plots/jpsi_v2_all_species.pdf");

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