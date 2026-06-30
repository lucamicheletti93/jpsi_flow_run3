void LoadStyle();
void SetLegend(TLegend *);
void SetGraph(TGraphErrors *, Color_t , double , int , double , bool );
TGraphAsymmErrors* DoGraphFromTheory(string);

void plot_results(TString model = "THU") {
    LoadStyle();
    gStyle -> SetLineStyleString(9, "80 20");

    // ***************************************************************************************** //
    // Theory predictions
    // ***************************************************************************************** //
    /*string fInNameCentr010 = Form("/Users/lucamicheletti/GITHUB/dq_run3_analyses/charmonia_production_pO_OO_NeNe/theory/%s/forward/data_v2_jpsi_010_y254.dat",model.Data());
    TGraphAsymmErrors *graTheorFwdCentr010CentrVal = DoGraphFromTheory(fInNameCentr010);
    gStyle -> SetLineStyleString(9,"80 20");
    graTheorFwdCentr010CentrVal -> SetLineColor(kOrange+7);
    graTheorFwdCentr010CentrVal -> SetLineWidth(2);
    graTheorFwdCentr010CentrVal -> SetLineStyle(9);*/

    string fInNameCentr020 = Form("/Users/lucamicheletti/GITHUB/dq_run3_analyses/charmonia_production_pO_OO_NeNe/theory/%s/forward/data_v2_jpsi_020_y254.dat",model.Data());
    TGraphAsymmErrors *graTheorFwdCentr020CentrVal = DoGraphFromTheory(fInNameCentr020);
    gStyle -> SetLineStyleString(9,"80 20");
    graTheorFwdCentr020CentrVal -> SetLineColor(kOrange+7);
    graTheorFwdCentr020CentrVal -> SetLineWidth(2);
    graTheorFwdCentr020CentrVal -> SetLineStyle(9);

    /*string fInNameCentr1050 = Form("/Users/lucamicheletti/GITHUB/dq_run3_analyses/charmonia_production_pO_OO_NeNe/theory/%s/forward/data_v2_jpsi_1050_y254.dat",model.Data());
    TGraphAsymmErrors *graTheorFwdCentr1050CentrVal = DoGraphFromTheory(fInNameCentr1050);
    gStyle -> SetLineStyleString(9,"80 20");
    graTheorFwdCentr1050CentrVal -> SetLineColor(kOrange+7);
    graTheorFwdCentr1050CentrVal -> SetLineWidth(2);
    graTheorFwdCentr1050CentrVal -> SetLineStyle(9);*/
    std::cout << "++++++++" << std::endl;

    string fInNameCentr2060 = Form("/Users/lucamicheletti/GITHUB/dq_run3_analyses/charmonia_production_pO_OO_NeNe/theory/%s/forward/data_v2_jpsi_2060_y254.dat",model.Data());
    TGraphAsymmErrors *graTheorFwdCentr2060CentrVal = DoGraphFromTheory(fInNameCentr2060);
    gStyle -> SetLineStyleString(9,"80 20");
    graTheorFwdCentr2060CentrVal -> SetLineColor(kOrange+7);
    graTheorFwdCentr2060CentrVal -> SetLineWidth(2);
    graTheorFwdCentr2060CentrVal -> SetLineStyle(9);

    // ***************************************************************************************** //
    // D-meson results
    // ***************************************************************************************** //
    TFile *fInDzero = new TFile("preliminary_plots/v2Dzero.root", "READ");
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
    //double v2JpsiFwdCentr020Vals[] = {-0.0405239, 0.0301976, 0.00319262, 0.0744494, 0.0698103, 0.0844539};
    //double v2JpsiFwdCentr020Stats[] = {0.0312219, 0.0221446, 0.0243449, 0.0300693, 0.0303064, 0.049283};
    //double v2JpsiFwdCentr020Systs[] = {0.00287178, 0.00929753, 0.0180505, 0.00464562, 0.0194577, 0.0269159};
    // TPC-ALL Q-vectors - ICHEP preliminaries
    double v2JpsiFwdCentr020Vals[] = {-0.0402795, 0.0301301, 0.00247346, 0.0744542, 0.072157, 0.0819288};
    double v2JpsiFwdCentr020Stats[] = {0.0310093, 0.0220362, 0.0241573, 0.0298396, 0.0300562, 0.0488713};
    double v2JpsiFwdCentr020Systs[] = {0.00285989, 0.00925626, 0.0183232, 0.00216452, 0.0069517, 0.0301406};

    Printf("J/psi v2 vs pT in centrality 0-20");
    for (int iPt = 0;iPt < nPtBins;iPt++) {
        Printf("%1.0f - %1.0f & %4.3f $?pm$ %4.3f $?pm$ %4.3f ??", minPtBins[iPt], maxPtBins[iPt], v2JpsiFwdCentr020Vals[iPt], v2JpsiFwdCentr020Stats[iPt], v2JpsiFwdCentr020Systs[iPt]);
    }

    TGraphErrors *graStatV2JpsiFwdCentr020 = new TGraphErrors(nPtBins, ptJpsiFwdCentrBins, v2JpsiFwdCentr020Vals, ptJpsiFwdWidthBins, v2JpsiFwdCentr020Stats);
    SetGraph(graStatV2JpsiFwdCentr020, kRed+1, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiFwdCentr020 = new TGraphErrors(nPtBins, ptJpsiFwdCentrBins, v2JpsiFwdCentr020Vals, ptJpsiFwdSystBins, v2JpsiFwdCentr020Systs);
    SetGraph(graSystV2JpsiFwdCentr020, kRed+1, 1.5, 20, 1, true);

    const int nLargePtBins = 4;
    double minLargePtBins[] = {0, 2, 4, 6};
    double maxLargePtBins[] = {2, 4, 6, 8};
    double ptJpsiFwdCentrLargeBins[] = {1, 3, 5, 7};
    double ptJpsiFwdWidthLargeBins[] = {1, 1, 1, 1};
    double ptJpsiFwdSystLargeBins[] = {0.15, 0.15, 0.15, 0.15};
    
    // TPC-ALL Q-vectors
    //double v2JpsiFwdCentr020LargeBinsVals[] = {0.0144323, 0.0329074, 0.0698103, 0.0844539};
    //double v2JpsiFwdCentr020LargeBinsStats[] = {0.0180728, 0.0188962, 0.0303064, 0.049283};
    //double v2JpsiFwdCentr020LargeBinsSysts[] = {0.00312387, 0.00928278, 0.0194577, 0.0269159};
    // TPC-ALL Q-vectors - ICHEP preliminaries
    double v2JpsiFwdCentr020LargeBinsVals[] = {0.0143719, 0.0321096, 0.072157, 0.0819288};
    double v2JpsiFwdCentr020LargeBinsStats[] = {0.0179254, 0.0187689, 0.0300562, 0.0488713};
    double v2JpsiFwdCentr020LargeBinsSysts[] = {0.00310963, 0.0101671, 0.0069517, 0.0301406};

    TGraphErrors *graStatV2JpsiFwdCentr020LargeBins = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr020LargeBinsVals, ptJpsiFwdWidthLargeBins, v2JpsiFwdCentr020LargeBinsStats);
    SetGraph(graStatV2JpsiFwdCentr020LargeBins, kRed+1, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiFwdCentr020LargeBins = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr020LargeBinsVals, ptJpsiFwdSystLargeBins, v2JpsiFwdCentr020LargeBinsSysts);
    SetGraph(graSystV2JpsiFwdCentr020LargeBins, kRed+1, 1.5, 20, 1, true);

    TGraphErrors *graStatV2JpsiFwdCentr020LargeBins_clone = (TGraphErrors*) graStatV2JpsiFwdCentr020LargeBins -> Clone("graStatV2JpsiFwdCentr020LargeBins_clone");
    SetGraph(graStatV2JpsiFwdCentr020LargeBins_clone, kOrange+7, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiFwdCentr020LargeBins_clone = (TGraphErrors*) graSystV2JpsiFwdCentr020LargeBins -> Clone("graSystV2JpsiFwdCentr020LargeBins_clone");
    SetGraph(graSystV2JpsiFwdCentr020LargeBins_clone, kOrange+7, 1.5, 20, 1, true);

    // Time association test
    double v2JpsiFwdCentr020LargeBinsTimeAssocVals[] = {0.012291, 0.0297956, 0.0698693, 0.0804606};
    double v2JpsiFwdCentr020LargeBinsTimeAssocStats[] = {0.0179296, 0.0187449, 0.0302308, 0.0486299}; 
    double v2JpsiFwdCentr020LargeBinsTimeAssocSysts[] = {0.00320326, 0.0107928, 0.0164601, 0.0237762};

    TGraphErrors *graStatV2JpsiFwdCentr020LargeBinsTimeAssoc = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr020LargeBinsTimeAssocVals, ptJpsiFwdWidthLargeBins, v2JpsiFwdCentr020LargeBinsTimeAssocStats);
    SetGraph(graStatV2JpsiFwdCentr020LargeBinsTimeAssoc, kBlack, 1.5, 24, 1, false);
    
    TGraphErrors *graSystV2JpsiFwdCentr020LargeBinsTimeAssoc = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr020LargeBinsTimeAssocVals, ptJpsiFwdSystLargeBins, v2JpsiFwdCentr020LargeBinsTimeAssocSysts);
    SetGraph(graSystV2JpsiFwdCentr020LargeBinsTimeAssoc, kBlack, 1.5, 24, 1, true);


    // Eta gap 2.5
    double v2JpsiFwdCentr020LargeEtaGapVals[] = {-0.0204004, 0.0205158, 0.0133768, 0.0813665, 0.0794812, 0.112315};
    double v2JpsiFwdCentr020LargeEtaGapStats[] = {0.039146, 0.0280243, 0.0303973, 0.0388386, 0.0377464, 0.060836};
    double v2JpsiFwdCentr020LargeEtaGapSysts[] = {0.0170701, 0.0106744, 0.0243818, 0.00436408, 0.00791315, 0.0621036};

    TGraphErrors *graStatV2JpsiFwdCentr020LargeEtaGap = new TGraphErrors(nPtBins, ptJpsiFwdCentrBins, v2JpsiFwdCentr020LargeEtaGapVals, ptJpsiFwdWidthBins, v2JpsiFwdCentr020LargeEtaGapStats);
    SetGraph(graStatV2JpsiFwdCentr020LargeEtaGap, kMagenta, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiFwdCentr020LargeEtaGap = new TGraphErrors(nPtBins, ptJpsiFwdCentrBins, v2JpsiFwdCentr020LargeEtaGapVals, ptJpsiFwdSystBins, v2JpsiFwdCentr020LargeEtaGapSysts);
    SetGraph(graSystV2JpsiFwdCentr020LargeEtaGap, kMagenta, 1.5, 20, 1, true);

    double v2JpsiFwdCentr020LargeBinsLargeEtaGapVals[] = {0.0135507, 0.0458208, 0.0794812, 0.112495};
    double v2JpsiFwdCentr020LargeBinsLargeEtaGapStats[] = {0.0227494, 0.0239035, 0.0377464, 0.0607986};
    double v2JpsiFwdCentr020LargeBinsLargeEtaGapSysts[] = {0.00571276, 0.010273, 0.00791315, 0.0619603};

    TGraphErrors *graStatV2JpsiFwdCentr020LargeBinsLargeEtaGap = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr020LargeBinsLargeEtaGapVals, ptJpsiFwdWidthLargeBins, v2JpsiFwdCentr020LargeBinsLargeEtaGapStats);
    SetGraph(graStatV2JpsiFwdCentr020LargeBinsLargeEtaGap, kGreen+2, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiFwdCentr020LargeBinsLargeEtaGap = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr020LargeBinsLargeEtaGapVals, ptJpsiFwdSystLargeBins, v2JpsiFwdCentr020LargeBinsLargeEtaGapSysts);
    SetGraph(graSystV2JpsiFwdCentr020LargeBinsLargeEtaGap, kGreen+2, 1.5, 20, 1, true);

    // Comparison 1.7 vs 2.5 eta gap
    double v2JpsiFwdCentr020EtaGap17[] = {0.0172906, 0.0235403, 0.0772568, 0.055712};
    double v2JpsiFwdCentr020EtaGap17Stats[] = {0.0166245, 0.0173056, 0.0274043, 0.0454146}; 
    double v2JpsiFwdCentr020EtaGap17Systs[] = {0.00296767, 0.009143, 0.00586852, 0.0317969};

    TGraphErrors *graStatV2JpsiFwdCentr020EtaGap17 = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr020EtaGap17, ptJpsiFwdWidthLargeBins, v2JpsiFwdCentr020EtaGap17Stats);
    SetGraph(graStatV2JpsiFwdCentr020EtaGap17, kRed+1, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiFwdCentr020EtaGap17 = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr020EtaGap17, ptJpsiFwdSystLargeBins, v2JpsiFwdCentr020EtaGap17Systs);
    SetGraph(graSystV2JpsiFwdCentr020EtaGap17, kRed+1, 1.5, 20, 1, true);

    double v2JpsiFwdCentr020EtaGap25[] = {0.00802074, 0.036169, 0.0748037, 0.0497307,};
    double v2JpsiFwdCentr020EtaGap25Stats[] = {0.0211309, 0.0220904, 0.0344394, 0.0563183};
    double v2JpsiFwdCentr020EtaGap25Systs[] = {0.00315599, 0.00286456, 0.00950001, 0.0416407};

    TGraphErrors *graStatV2JpsiFwdCentr020EtaGap25 = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr020EtaGap25, ptJpsiFwdWidthLargeBins, v2JpsiFwdCentr020EtaGap25Stats);
    SetGraph(graStatV2JpsiFwdCentr020EtaGap25, kBlack, 1.5, 24, 1, false);
    
    TGraphErrors *graSystV2JpsiFwdCentr020EtaGap25 = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr020EtaGap25, ptJpsiFwdSystLargeBins, v2JpsiFwdCentr020EtaGap25Systs);
    SetGraph(graSystV2JpsiFwdCentr020EtaGap25, kBlack, 1.5, 24, 1, true);


    // ***************************************************************************************** //
    // Centrality 20-60%
    // ***************************************************************************************** //
    // TPC-POS only Q-vectors
    //double v2JpsiFwdCentr2060Vals[] = {-0.0162254, -0.0178981, 0.00161603, 0.0205455, 0.0736224, 0.104662};
    //double v2JpsiFwdCentr2060Stats[] = {0.0388019, 0.0276514, 0.0301216, 0.0367978, 0.0371583, 0.066764};
    //double v2JpsiFwdCentr2060Systs[] = {0.0106743, 0.0289646, 0.00971417, 0.0148252, 0.00970835, 0.0282096};
    // TPC-ALL only Q-vectors
    //double v2JpsiFwdCentr2060Vals[] = {-0.017018, -0.0210286, 0.000239199, 0.0185226, 0.0418244, 0.12257};
    //double v2JpsiFwdCentr2060Stats[] = {0.0370604, 0.0263352, 0.0286848, 0.0347268, 0.0350661, 0.0639947};
    //double v2JpsiFwdCentr2060Systs[] = {0.0145221, 0.0218097, 0.0039893, 0.00969075, 0.00867721, 0.0234387};
    // TPC-ALL only Q-vectors - ICHEP preliminaries
    double v2JpsiFwdCentr2060Vals[] = {-0.0165498, -0.0203789, 0.000216791, 0.0184516, 0.0413613, 0.119999};
    double v2JpsiFwdCentr2060Stats[] = {0.0368338, 0.0261904, 0.0284674, 0.0344189, 0.0348384, 0.0635753};
    double v2JpsiFwdCentr2060Systs[] = {0.0140023, 0.0215822, 0.00395956, 0.00944929, 0.00866651, 0.0273761};

    Printf("J/psi v2 vs pT in centrality 20-60");
    for (int iPt = 0;iPt < nPtBins;iPt++) {
        Printf("%1.0f - %1.0f & %4.3f $?pm$ %4.3f $?pm$ %4.3f ??", minPtBins[iPt], maxPtBins[iPt], v2JpsiFwdCentr2060Vals[iPt], v2JpsiFwdCentr2060Stats[iPt], v2JpsiFwdCentr2060Systs[iPt]);
    }

    TGraphErrors *graStatV2JpsiFwdCentr2060 = new TGraphErrors(6, ptJpsiFwdCentrBins, v2JpsiFwdCentr2060Vals, ptJpsiFwdWidthBins, v2JpsiFwdCentr2060Stats);
    SetGraph(graStatV2JpsiFwdCentr2060, kRed+1, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiFwdCentr2060 = new TGraphErrors(6, ptJpsiFwdCentrBins, v2JpsiFwdCentr2060Vals, ptJpsiFwdSystBins, v2JpsiFwdCentr2060Systs);
    SetGraph(graSystV2JpsiFwdCentr2060, kRed+1, 1.5, 20, 1, true);

    TGraphErrors *graStatV2JpsiFwdCentr2060_clone = (TGraphErrors*) graStatV2JpsiFwdCentr2060 -> Clone("graStatV2JpsiFwdCentr2060_clone");
    SetGraph(graStatV2JpsiFwdCentr2060_clone, kBlue+1, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiFwdCentr2060_clone = (TGraphErrors*) graSystV2JpsiFwdCentr2060 -> Clone("graStatV2JpsiFwdCentr2060_clone");
    SetGraph(graSystV2JpsiFwdCentr2060_clone, kBlue+1, 1.5, 20, 1, true);

    // TPC-ALL only Q-vectors
    //double v2JpsiFwdCentr2060LargeBinsVals[] = {-0.00894195, 0.00854762, 0.0418244, 0.12257};
    //double v2JpsiFwdCentr2060LargeBinsStats[] = {0.0214885, 0.0221865, 0.0350661, 0.0639947};
    //double v2JpsiFwdCentr2060LargeBinsSysts[] = {0.00690191, 0.00723872, 0.00867721, 0.0234387};
    // TPC-ALL only Q-vectors - ICHEP preliminaries
    double v2JpsiFwdCentr2060LargeBinsVals[] = {-0.00880563, 0.00844508, 0.0413613, 0.119999};
    double v2JpsiFwdCentr2060LargeBinsStats[] = {0.0213792, 0.0220216, 0.0348384, 0.0635753};
    double v2JpsiFwdCentr2060LargeBinsSysts[] = {0.0069356, 0.00718644, 0.00866651, 0.0273761};

    TGraphErrors *graStatV2JpsiFwdCentr2060LargeBins = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr2060LargeBinsVals, ptJpsiFwdWidthLargeBins, v2JpsiFwdCentr2060LargeBinsStats);
    SetGraph(graStatV2JpsiFwdCentr2060LargeBins, kRed+1, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiFwdCentr2060LargeBins = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr2060LargeBinsVals, ptJpsiFwdSystLargeBins, v2JpsiFwdCentr2060LargeBinsSysts);
    SetGraph(graSystV2JpsiFwdCentr2060LargeBins, kRed+1, 1.5, 20, 1, true);

    TGraphErrors *graStatV2JpsiFwdCentr2060LargeBins_clone = (TGraphErrors*) graStatV2JpsiFwdCentr2060LargeBins -> Clone("graStatV2JpsiFwdCentr2060LargeBins_clone");
    SetGraph(graStatV2JpsiFwdCentr2060LargeBins_clone, kAzure+2, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiFwdCentr2060LargeBins_clone = (TGraphErrors*) graSystV2JpsiFwdCentr2060LargeBins -> Clone("graSystV2JpsiFwdCentr2060LargeBins_clone");
    SetGraph(graSystV2JpsiFwdCentr2060LargeBins_clone, kAzure+2, 1.5, 20, 1, true);

    // Time association test
    double v2JpsiFwdCentr2060LargeBinsTimeAssocVals[] = {-0.00785662, 0.00842857, 0.0366479, 0.124399};
    double v2JpsiFwdCentr2060LargeBinsTimeAssocStats[] = {0.0213507, 0.0220006, 0.0348177, 0.063143}; 
    double v2JpsiFwdCentr2060LargeBinsTimeAssocSysts[] = {0.00932627, 0.00554435, 0.0210392, 0.0219922};

    TGraphErrors *graStatV2JpsiFwdCentr2060LargeBinsTimeAssoc = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr2060LargeBinsTimeAssocVals, ptJpsiFwdWidthLargeBins, v2JpsiFwdCentr2060LargeBinsTimeAssocStats);
    SetGraph(graStatV2JpsiFwdCentr2060LargeBinsTimeAssoc, kBlack, 1.5, 24, 1, false);
    
    TGraphErrors *graSystV2JpsiFwdCentr2060LargeBinsTimeAssoc = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr2060LargeBinsTimeAssocVals, ptJpsiFwdSystLargeBins, v2JpsiFwdCentr2060LargeBinsTimeAssocSysts);
    SetGraph(graSystV2JpsiFwdCentr2060LargeBinsTimeAssoc, kBlack, 1.5, 24, 1, true);

    // Eta gap 2.5
    double v2JpsiFwdCentr2060LargeEtaGapVals[] = {0.0343215, 0.0386001, 0.0270601, -0.0420709, 0.060201, 0.21096};
    double v2JpsiFwdCentr2060LargeEtaGapStats[] = {0.0487787, 0.0350896, 0.0381595, 0.0459997, 0.0463647, 0.0810375};
    double v2JpsiFwdCentr2060LargeEtaGapSysts[] = {0.0318106, 0.0145527, 0.0283817, 0.00834783, 0.0302841, 0.0427581};

    TGraphErrors *graStatV2JpsiFwdCentr2060LargeEtaGap = new TGraphErrors(nPtBins, ptJpsiFwdCentrBins, v2JpsiFwdCentr2060LargeEtaGapVals, ptJpsiFwdWidthBins, v2JpsiFwdCentr2060LargeEtaGapStats);
    SetGraph(graStatV2JpsiFwdCentr2060LargeEtaGap, kBlue+2, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiFwdCentr2060LargeEtaGap = new TGraphErrors(nPtBins, ptJpsiFwdCentrBins, v2JpsiFwdCentr2060LargeEtaGapVals, ptJpsiFwdSystBins, v2JpsiFwdCentr2060LargeEtaGapSysts);
    SetGraph(graSystV2JpsiFwdCentr2060LargeEtaGap, kBlue+2, 1.5, 20, 1, true);


    double v2JpsiFwdCentr2060LargeBinsLargeEtaGapVals[] = {0.0270983, 0.00583639, 0.0613767, 0.210388};
    double v2JpsiFwdCentr2060LargeBinsLargeEtaGapStats[] = {0.0283058, 0.0294719, 0.0463647, 0.0810375};
    double v2JpsiFwdCentr2060LargeBinsLargeEtaGapSysts[] = {0.0286739, 0.0125549, 0.030383, 0.0444629};

    TGraphErrors *graStatV2JpsiFwdCentr2060LargeBinsLargeEtaGap = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr2060LargeBinsLargeEtaGapVals, ptJpsiFwdWidthLargeBins, v2JpsiFwdCentr2060LargeBinsLargeEtaGapStats);
    SetGraph(graStatV2JpsiFwdCentr2060LargeBinsLargeEtaGap, kGreen+2, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiFwdCentr2060LargeBinsLargeEtaGap = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr2060LargeBinsLargeEtaGapVals, ptJpsiFwdSystLargeBins, v2JpsiFwdCentr2060LargeBinsLargeEtaGapSysts);
    SetGraph(graSystV2JpsiFwdCentr2060LargeBinsLargeEtaGap, kGreen+2, 1.5, 20, 1, true);

    // Comparison 1.7 vs 2.5 eta gap
    double v2JpsiFwdCentr2060EtaGap17[] = {-0.0154992, 0.0146533, 0.0430052, 0.100778};
    double v2JpsiFwdCentr2060EtaGap17Stats[] = {0.0200534, 0.0205939, 0.032431, 0.058807}; 
    double v2JpsiFwdCentr2060EtaGap17Systs[] = {0.00678459, 0.00773385, 0.0256593, 0.0242779};

    TGraphErrors *graStatV2JpsiFwdCentr2060EtaGap17 = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr2060EtaGap17, ptJpsiFwdWidthLargeBins, v2JpsiFwdCentr2060EtaGap17Stats);
    SetGraph(graStatV2JpsiFwdCentr2060EtaGap17, kRed+1, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiFwdCentr2060EtaGap17 = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr2060EtaGap17, ptJpsiFwdSystLargeBins, v2JpsiFwdCentr2060EtaGap17Systs);
    SetGraph(graSystV2JpsiFwdCentr2060EtaGap17, kRed+1, 1.5, 20, 1, true);

    double v2JpsiFwdCentr2060EtaGap25[] = {-0.00104163, 0.017433, 0.0912705, 0.16855};
    double v2JpsiFwdCentr2060EtaGap25Stats[] = {0.0266955, 0.0275955, 0.0432319, 0.0753419};
    double v2JpsiFwdCentr2060EtaGap25Systs[] = {0.0287657, 0.0119787, 0.0299155, 0.0458988};

    TGraphErrors *graStatV2JpsiFwdCentr2060EtaGap25 = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr2060EtaGap25, ptJpsiFwdWidthLargeBins, v2JpsiFwdCentr2060EtaGap25Stats);
    SetGraph(graStatV2JpsiFwdCentr2060EtaGap25, kBlack, 1.5, 24, 1, false);
    
    TGraphErrors *graSystV2JpsiFwdCentr2060EtaGap25 = new TGraphErrors(nLargePtBins, ptJpsiFwdCentrLargeBins, v2JpsiFwdCentr2060EtaGap25, ptJpsiFwdSystLargeBins, v2JpsiFwdCentr2060EtaGap25Systs);
    SetGraph(graSystV2JpsiFwdCentr2060EtaGap25, kBlack, 1.5, 24, 1, true);


    double ptJpsiOOFwdCentr020[] = {1.28, 2.87, 4.82, 6.84};
    double ptJpsiOOFwdCentr2060[] = {1.24, 2.87, 4.83, 6.76};
    Printf("J/psi v2 vs pT in centrality 0-20 & 20-60");
    for (int iPt = 0;iPt < nPtBins;iPt++) {
        double minPt = minPtBins[iPt];
        double maxPt = maxPtBins[iPt];
        double meanPt020 = (minPt + maxPt)/2.;
        double val020 = v2JpsiFwdCentr020Vals[iPt];
        double stat020 = v2JpsiFwdCentr020Stats[iPt];
        double syst020 = v2JpsiFwdCentr020Systs[iPt];
        double meanPt2060 = (minPt + maxPt)/2.;
        double val2060 = v2JpsiFwdCentr2060Vals[iPt];
        double stat2060 = v2JpsiFwdCentr2060Stats[iPt];
        double syst2060 = v2JpsiFwdCentr2060Systs[iPt];
        Printf("%1.0f - %1.0f & %1.2f & %4.3f $?pm$ %4.3f $?pm$ %4.3f & %1.2f & %4.3f $?pm$ %4.3f $?pm$ %4.3f ??", minPt, maxPt, meanPt020, val020, stat020, syst020, meanPt2060, val2060, stat2060, syst2060);
    }
    for (int iPt = 0;iPt < nLargePtBins;iPt++) {
        double minPt = minLargePtBins[iPt];
        double maxPt = maxLargePtBins[iPt];
        double meanPt020 = (minPt + maxPt)/2.;
        double val020 = v2JpsiFwdCentr020LargeBinsVals[iPt];
        double stat020 = v2JpsiFwdCentr020LargeBinsStats[iPt];
        double syst020 = v2JpsiFwdCentr020LargeBinsSysts[iPt];
        double meanPt2060 = (minPt + maxPt)/2.;
        double val2060 = v2JpsiFwdCentr2060LargeBinsVals[iPt];
        double stat2060 = v2JpsiFwdCentr2060LargeBinsStats[iPt];
        double syst2060 = v2JpsiFwdCentr2060LargeBinsSysts[iPt];
        Printf("%1.0f - %1.0f & %1.2f & %4.3f $?pm$ %4.3f $?pm$ %4.3f & %1.2f & %4.3f $?pm$ %4.3f $?pm$ %4.3f ??", minPt, maxPt, meanPt020, val020, stat020, syst020, meanPt2060, val2060, stat2060, syst2060);
    }

    // ***************************************************************************************** //
    // Preparing other results...
    // ***************************************************************************************** //

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
    SetGraph(graStatV2JpsiPbPbFwdCentr010, kGray+2, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiPbPbFwdCentr010 = new TGraphErrors(nPtBinsPbPb, ptJpsiFwdCentrBinsPbPb, v2JpsiPbPbFwdCentr010Vals, ptJpsiFwdSystBinsPbPb, v2JpsiPbPbFwdCentr010Systs);
    SetGraph(graSystV2JpsiPbPbFwdCentr010, kGray+2, 1.5, 20, 1, true);


    double v2JpsiPbPbFwdCentr3050Vals[] = {-0.00751166, 0.0324682, 0.0483128, 0.0623836, 0.0790007, 0.0750308, 0.0821814, 0.0844028, 0.110533, 0.0894459, 0.0864354, 0.0662619, 0.105772};
    double v2JpsiPbPbFwdCentr3050Stats[] = {0.0138548, 0.0119431, 0.00818823, 0.010507, 0.00958244, 0.00985444, 0.00793826, 0.00942195, 0.0112991, 0.011538, 0.019468, 0.0290814, 0.0375577};
    double v2JpsiPbPbFwdCentr3050Systs[] = {0.00456438, 0.00532624, 0.00657618, 0.00767329, 0.00376474, 0.00244701, 0.00190689, 0.00266505, 0.00643058, 0.00989524, 0.00819967, 0.0119924, 0.012523};

    TGraphErrors *graStatV2JpsiPbPbFwdCentr3050 = new TGraphErrors(nPtBinsPbPb, ptJpsiFwdCentrBinsPbPb, v2JpsiPbPbFwdCentr3050Vals, ptJpsiFwdWidthBinsPbPb, v2JpsiPbPbFwdCentr3050Stats);
    SetGraph(graStatV2JpsiPbPbFwdCentr3050, kGray+2, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiPbPbFwdCentr3050 = new TGraphErrors(nPtBinsPbPb, ptJpsiFwdCentrBinsPbPb, v2JpsiPbPbFwdCentr3050Vals, ptJpsiFwdSystBinsPbPb, v2JpsiPbPbFwdCentr3050Systs);
    SetGraph(graSystV2JpsiPbPbFwdCentr3050, kGray+2, 1.5, 20, 1, true);

    //=======================================================================//
    double ptJpsiMidWidthBins[] = {0. , 0., 0.};
    double ptJpsiMidSystBins[] = {0.15, 0.15, 0.15};
    
    double ptJpsiMidCentr020[] = {1.230, 2.870, 5.300};
    double v2JpsiMidCentr020Vals[] = {-0.180000, -0.018714, 0.059150};
    double v2JpsiMidCentr020Stats[] = {0.102775, 0.102947, 0.058341};
    double v2JpsiMidCentr020Systs[] = {0.050000, 0.090000, 0.060000};

    TGraphErrors *graStatV2JpsiMidCentr020 = new TGraphErrors(3, ptJpsiMidCentr020, v2JpsiMidCentr020Vals, ptJpsiMidWidthBins, v2JpsiMidCentr020Stats);
    SetGraph(graStatV2JpsiMidCentr020, kGreen+2, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiMidCentr020 = new TGraphErrors(3, ptJpsiMidCentr020, v2JpsiMidCentr020Vals, ptJpsiMidSystBins, v2JpsiMidCentr020Systs);
    SetGraph(graSystV2JpsiMidCentr020, kGreen+2, 1.5, 20, 1, true);

    double ptJpsiMidCentr2060[] = {1.190, 2.920 , 5.270};
    double v2JpsiMidCentr2060Vals[] = {-0.032969, 0.070266, 0.145450};
    double v2JpsiMidCentr2060Stats[] = {0.097961, 0.079300, 0.067461};
    double v2JpsiMidCentr2060Systs[] = {0.070000, 0.080000, 0.060000};

    TGraphErrors *graStatV2JpsiMidCentr2060 = new TGraphErrors(3, ptJpsiMidCentr2060, v2JpsiMidCentr2060Vals, ptJpsiMidWidthBins, v2JpsiMidCentr2060Stats);
    SetGraph(graStatV2JpsiMidCentr2060, kGreen+2, 1.5, 20, 1, false);
    
    TGraphErrors *graSystV2JpsiMidCentr2060 = new TGraphErrors(3, ptJpsiMidCentr2060, v2JpsiMidCentr2060Vals, ptJpsiMidSystBins, v2JpsiMidCentr2060Systs);
    SetGraph(graSystV2JpsiMidCentr2060, kGreen+2, 1.5, 20, 1, true);

    
    // ======================================================================= //
    // Plots
    // ======================================================================= //

    TLine *lineUnity = new TLine(0, 0, 8, 0);
    lineUnity -> SetLineColor(kGray+2);
    lineUnity -> SetLineWidth(2);
    lineUnity -> SetLineStyle(kDashed);

    TLatex *latexTitle = new TLatex();
    latexTitle -> SetTextSize(0.050);
    latexTitle -> SetNDC();
    latexTitle -> SetTextFont(42);

    // ************************************************************************** //
    // Jpsi Fwd 0-20%
    // Default
    TCanvas *canvasJpsiFwdCentr020 = new TCanvas("canvasJpsiFwdCentr020", "", 800, 600);
    TH2D *histGridJpsiFwdCentr020  = new TH2D("histGridJpsiFwdCentr020", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}^{SP} {|#Delta#it{#eta}| > 1.7}", 100, 0, 8, 100, -0.20, 0.30);
    histGridJpsiFwdCentr020 -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiFwdCentr020LargeBins_clone -> Draw("E2P");
    graStatV2JpsiFwdCentr020LargeBins_clone -> Draw("P SAME");
    graSystV2JpsiFwdCentr020 -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr020 -> Draw("P SAME");
    latexTitle -> DrawLatex(0.20, 0.88, "ALICE work in progress");
    latexTitle -> DrawLatex(0.20, 0.82, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, 0#minus20%");
    latexTitle -> DrawLatex(0.20, 0.76, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");

    // Default vs time assoc
    TCanvas *canvasJpsiFwdCentr020vsTimeAssoc = new TCanvas("canvasJpsiFwdCentr020vsTimeAssoc", "", 800, 600);
    histGridJpsiFwdCentr020 -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiFwdCentr020LargeBins -> Draw("E2P");
    graStatV2JpsiFwdCentr020LargeBins -> Draw("P SAME");
    graSystV2JpsiFwdCentr020LargeBinsTimeAssoc -> Draw("P SAME");
    graStatV2JpsiFwdCentr020LargeBinsTimeAssoc -> Draw("P SAME");
    latexTitle -> DrawLatex(0.20, 0.88, "ALICE work in progress");
    latexTitle -> DrawLatex(0.20, 0.82, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, 0#minus20%");
    latexTitle -> DrawLatex(0.20, 0.76, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");

    TLegend *legendJpsiFwdCentr020vsTimeAssoc = new TLegend(0.60,0.25,0.85,0.40);
    SetLegend(legendJpsiFwdCentr020vsTimeAssoc);
    legendJpsiFwdCentr020vsTimeAssoc -> AddEntry(graSystV2JpsiFwdCentr020LargeBins,"Standard assoc.","FP");
    legendJpsiFwdCentr020vsTimeAssoc -> AddEntry(graSystV2JpsiFwdCentr020LargeBinsTimeAssoc,"Time assoc.","FP");
    legendJpsiFwdCentr020vsTimeAssoc -> Draw();

    // OO vs Pb-Pb
    TCanvas *canvasJpsiFwdCentr020vsPbPb = new TCanvas("canvasJpsiFwdCentr020vsPbPb", "", 800, 600);
    histGridJpsiFwdCentr020 -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiPbPbFwdCentr010 -> Draw("E2P");
    graStatV2JpsiPbPbFwdCentr010 -> Draw("P SAME");
    graSystV2JpsiFwdCentr020 -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr020 -> Draw("P SAME");
    latexTitle -> DrawLatex(0.20, 0.88, "ALICE work in progress");
    latexTitle -> DrawLatex(0.20, 0.82, "#sqrt{#it{s}_{NN}} = 5.36 TeV, J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");

    TLegend *legendJpsiFwdCentr020vsPbPb = new TLegend(0.60,0.25,0.85,0.40);
    SetLegend(legendJpsiFwdCentr020vsPbPb);
    legendJpsiFwdCentr020vsPbPb -> AddEntry(graSystV2JpsiFwdCentr020,"OO, 0#minus20%","FP");
    legendJpsiFwdCentr020vsPbPb -> AddEntry(graSystV2JpsiPbPbFwdCentr010,"Pb-Pb, 0#minus10%","FP");
    legendJpsiFwdCentr020vsPbPb -> Draw();

    // OO vs Pb-Pb, large bins
    TCanvas *canvasJpsiFwdCentr020vsPbPbLargeBins = new TCanvas("canvasJpsiFwdCentr020vsPbPbLargeBins", "", 800, 600);
    histGridJpsiFwdCentr020 -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiPbPbFwdCentr010 -> Draw("E2P");
    graStatV2JpsiPbPbFwdCentr010 -> Draw("P SAME");
    graSystV2JpsiFwdCentr020LargeBins -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr020LargeBins -> Draw("P SAME");
    latexTitle -> DrawLatex(0.20, 0.88, "ALICE work in progress");
    latexTitle -> DrawLatex(0.20, 0.82, "#sqrt{#it{s}_{NN}} = 5.36 TeV, J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");
    legendJpsiFwdCentr020vsPbPb -> Draw();

    TCanvas *canvasJpsiFwdVsJpsiMidCentr020 = new TCanvas("canvasJpsiFwdVsJpsiMidCentr020", "", 800, 600);
    TH2D *histGridJpsiFwdVsJpsiMid  = new TH2D("histGridJpsiFwdVsJpsiMid", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}^{SP}", 100, 0, 8, 100, -0.30, 0.30);
    histGridJpsiFwdVsJpsiMid -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiMidCentr020 -> Draw("E2P");
    graStatV2JpsiMidCentr020 -> Draw("P");
    graSystV2JpsiFwdCentr020 -> Draw("E2P");
    graStatV2JpsiFwdCentr020 -> Draw("P");

    TLegend *legendJpsiFwdVsJpsiMidCentr020 = new TLegend(0.55,0.20,0.75,0.35);
    SetLegend(legendJpsiFwdVsJpsiMidCentr020);
    legendJpsiFwdVsJpsiMidCentr020 -> AddEntry(graSystV2JpsiFwdCentr020,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4","P");
    legendJpsiFwdVsJpsiMidCentr020 -> AddEntry(graSystV2JpsiMidCentr020,"J/#psi #rightarrow e^{+}e^{-}, |#it{y}| < 0.9","P");
    legendJpsiFwdVsJpsiMidCentr020 -> Draw();

    latexTitle -> DrawLatex(0.20, 0.88, "ALICE work in progress");
    latexTitle -> DrawLatex(0.20, 0.82, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, 0#minus20%");

    TCanvas *canvasJpsiFwdVsJpsiMidCentr020LargeBins = new TCanvas("canvasJpsiFwdVsJpsiMidCentr020LargeBins", "", 800, 600);
    histGridJpsiFwdVsJpsiMid -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiMidCentr020 -> Draw("E2P");
    graStatV2JpsiMidCentr020 -> Draw("P");
    graSystV2JpsiFwdCentr020LargeBins -> Draw("E2P");
    graStatV2JpsiFwdCentr020LargeBins -> Draw("P");
    legendJpsiFwdVsJpsiMidCentr020 -> Draw();

    latexTitle -> DrawLatex(0.20, 0.88, "ALICE work in progress");
    latexTitle -> DrawLatex(0.20, 0.82, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, 0#minus20%");

    // ************************************************************************** //
    // Jpsi Fwd 20-60%
    // Default
    TCanvas *canvasJpsiFwdCentr2060 = new TCanvas("canvasJpsiFwdCentr2060", "", 800, 600);
    TH2D *histGridJpsiFwdCentr2060  = new TH2D("histGridJpsiFwdCentr2060", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}^{SP} {|#Delta#it{#eta}| > 1.7}", 100, 0, 8, 100, -0.20, 0.30);
    histGridJpsiFwdCentr2060 -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiFwdCentr2060LargeBins_clone -> Draw("E2P");
    graStatV2JpsiFwdCentr2060LargeBins_clone -> Draw("P SAME");
    graSystV2JpsiFwdCentr2060_clone -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr2060_clone -> Draw("P SAME");
    latexTitle -> DrawLatex(0.20, 0.88, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, 20#minus60%");

    // Default vs time assoc
    TCanvas *canvasJpsiFwdCentr2060vsTimeAssoc = new TCanvas("canvasJpsiFwdCentr2060vsTimeAssoc", "", 800, 600);
    histGridJpsiFwdCentr2060 -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiFwdCentr2060LargeBins -> Draw("E2P");
    graStatV2JpsiFwdCentr2060LargeBins -> Draw("P SAME");
    graSystV2JpsiFwdCentr2060LargeBinsTimeAssoc -> Draw("P SAME");
    graStatV2JpsiFwdCentr2060LargeBinsTimeAssoc -> Draw("P SAME");
    latexTitle -> DrawLatex(0.20, 0.88, "ALICE work in progress");
    latexTitle -> DrawLatex(0.20, 0.82, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, 20#minus60%");
    latexTitle -> DrawLatex(0.20, 0.76, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");

    TLegend *legendJpsiFwdCentr2060vsTimeAssoc = new TLegend(0.60,0.25,0.85,0.40);
    SetLegend(legendJpsiFwdCentr2060vsTimeAssoc);
    legendJpsiFwdCentr2060vsTimeAssoc -> AddEntry(graSystV2JpsiFwdCentr2060LargeBins,"Standard assoc.","FP");
    legendJpsiFwdCentr2060vsTimeAssoc -> AddEntry(graSystV2JpsiFwdCentr2060LargeBinsTimeAssoc,"Time assoc.","FP");
    legendJpsiFwdCentr2060vsTimeAssoc -> Draw();

    TCanvas *canvasJpsiFwdCentr2060vsPbPb = new TCanvas("canvasJpsiFwdCentr2060vsPbPb", "", 800, 600);
    histGridJpsiFwdCentr2060 -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiPbPbFwdCentr3050 -> Draw("E2P");
    graStatV2JpsiPbPbFwdCentr3050 -> Draw("P SAME");
    graSystV2JpsiFwdCentr2060 -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr2060 -> Draw("P SAME");
    latexTitle -> DrawLatex(0.20, 0.88, "ALICE work in progress");
    latexTitle -> DrawLatex(0.20, 0.82, "#sqrt{#it{s}_{NN}} = 5.36 TeV, J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");

    TLegend *legendJpsiFwdCentr2060vsPbPb = new TLegend(0.60,0.25,0.85,0.40);
    SetLegend(legendJpsiFwdCentr2060vsPbPb);
    legendJpsiFwdCentr2060vsPbPb -> AddEntry(graSystV2JpsiFwdCentr2060,"OO, 20#minus60%","FP");
    legendJpsiFwdCentr2060vsPbPb -> AddEntry(graSystV2JpsiPbPbFwdCentr010,"Pb-Pb, 30#minus50%","FP");
    legendJpsiFwdCentr2060vsPbPb -> Draw();

    TCanvas *canvasJpsiFwdCentr2060vsPbPbLargeBins = new TCanvas("canvasJpsiFwdCentr2060vsPbPbLargeBins", "", 800, 600);
    histGridJpsiFwdCentr2060 -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiPbPbFwdCentr3050 -> Draw("E2P");
    graStatV2JpsiPbPbFwdCentr3050 -> Draw("P SAME");
    graSystV2JpsiFwdCentr2060LargeBins -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr2060LargeBins -> Draw("P SAME");
    latexTitle -> DrawLatex(0.20, 0.88, "ALICE work in progress");
    latexTitle -> DrawLatex(0.20, 0.82, "#sqrt{#it{s}_{NN}} = 5.36 TeV, J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");
    legendJpsiFwdCentr2060vsPbPb->Draw();


    TCanvas *canvasJpsiFwdVsJpsiMidCentr2060 = new TCanvas("canvasJpsiFwdVsJpsiMidCentr2060", "", 800, 600);
    histGridJpsiFwdVsJpsiMid -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiMidCentr2060 -> Draw("E2P");
    graStatV2JpsiMidCentr2060 -> Draw("P");
    graSystV2JpsiFwdCentr2060 -> Draw("E2P");
    graStatV2JpsiFwdCentr2060 -> Draw("P");

    TLegend *legendJpsiFwdVsJpsiMidCentr2060 = new TLegend(0.55,0.20,0.75,0.35);
    SetLegend(legendJpsiFwdVsJpsiMidCentr2060);
    legendJpsiFwdVsJpsiMidCentr2060 -> AddEntry(graSystV2JpsiFwdCentr2060,"J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4","P");
    legendJpsiFwdVsJpsiMidCentr2060 -> AddEntry(graSystV2JpsiMidCentr2060,"J/#psi #rightarrow e^{+}e^{-}, |#it{y}| < 0.9","P");
    legendJpsiFwdVsJpsiMidCentr2060 -> Draw();

    latexTitle -> DrawLatex(0.20, 0.88, "ALICE work in progress");
    latexTitle -> DrawLatex(0.20, 0.82, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, 20#minus60%");

    TCanvas *canvasJpsiFwdVsJpsiMidCentr2060LargeBins = new TCanvas("canvasJpsiFwdVsJpsiMidCentr2060LargeBins", "", 800, 600);
    histGridJpsiFwdVsJpsiMid -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiMidCentr2060 -> Draw("E2P");
    graStatV2JpsiMidCentr2060 -> Draw("P");
    graSystV2JpsiFwdCentr2060LargeBins -> Draw("E2P");
    graStatV2JpsiFwdCentr2060LargeBins -> Draw("P");
    legendJpsiFwdVsJpsiMidCentr2060 -> Draw();

    latexTitle -> DrawLatex(0.20, 0.88, "ALICE work in progress");
    latexTitle -> DrawLatex(0.20, 0.82, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, 20#minus60%");

    // Jpsi Fwd 0-20 vs 20-60%
    TCanvas *canvasJpsiFwd = new TCanvas("canvasJpsiFwd", "", 800, 600);
    TH2D *histGridJpsiFwd  = new TH2D("histGridJpsiFwd", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}^{SP} {|#Delta#it{#eta}| > 1.7}", 100, 0, 8, 100, -0.20, 0.30);
    histGridJpsiFwd -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiFwdCentr020 -> Draw("E2P");
    graStatV2JpsiFwdCentr020 -> Draw("P SAME");
    graSystV2JpsiFwdCentr2060_clone -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr2060_clone -> Draw("P SAME");
    latexTitle -> DrawLatex(0.20, 0.88, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");

    TLegend *legendJpsiFwd = new TLegend(0.20,0.65,0.40,0.80);
    SetLegend(legendJpsiFwd);
    legendJpsiFwd -> AddEntry(graSystV2JpsiFwdCentr020,"0#minus20%","P");
    legendJpsiFwd -> AddEntry(graSystV2JpsiFwdCentr2060_clone,"20#minus60%","P");
    legendJpsiFwd -> Draw();

    TCanvas *canvasJpsiFwdLargeBins = new TCanvas("canvasJpsiFwdLargeBins", "", 800, 600);
    histGridJpsiFwd -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiFwdCentr020LargeBins -> Draw("E2P");
    graStatV2JpsiFwdCentr020LargeBins -> Draw("P SAME");
    graSystV2JpsiFwdCentr2060LargeBins_clone -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr2060LargeBins_clone -> Draw("P SAME");
    latexTitle -> DrawLatex(0.20, 0.88, "OO, #sqrt{#it{s}_{NN}} = 5.36 TeV, J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4");

    TLegend *legendJpsiFwdLargeBins = new TLegend(0.20,0.65,0.40,0.80);
    SetLegend(legendJpsiFwdLargeBins);
    legendJpsiFwdLargeBins -> AddEntry(graSystV2JpsiFwdCentr020LargeBins,"0#minus20%","P");
    legendJpsiFwdLargeBins -> AddEntry(graSystV2JpsiFwdCentr2060LargeBins_clone,"20#minus60%","P");
    legendJpsiFwdLargeBins -> Draw();


    // Jpsi Fwd vs theory
    TCanvas *canvasJpsiFwdVsTheorCentr020 = new TCanvas("canvasJpsiFwdVsTheorCentr020", "", 800, 600);
    TH2D *histGridJpsiFwdVsTheor  = new TH2D("histGridJpsiFwdVsTheor", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}^{SP}", 100, 0, 8, 100, -0.30, 0.30);
    histGridJpsiFwdVsTheor -> Draw();
    lineUnity -> Draw("SAME");
    graTheorFwdCentr020CentrVal -> Draw("L SAME");
    graSystV2JpsiFwdCentr020 -> Draw("E2P");
    graStatV2JpsiFwdCentr020 -> Draw("P");

    TLegend *legendJpsiFwdVsTheorCentr020 = new TLegend(0.20,0.25,0.45,0.40);
    SetLegend(legendJpsiFwdVsTheorCentr020);
    legendJpsiFwdVsTheorCentr020 -> AddEntry(graSystV2JpsiFwdCentr020LargeBins,"Data, SP, |#Delta#it{#eta}| > 1.7","FP");
    legendJpsiFwdVsTheorCentr020 -> AddEntry(graTheorFwdCentr020CentrVal,"Tsinghua Transport + MUSIC Hydro.","L");
    legendJpsiFwdVsTheorCentr020 -> Draw();

    latexTitle -> DrawLatex(0.20, 0.88, "ALICE Preliminary, OO, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.20, 0.82, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4, 0#minus20%");

    TCanvas *canvasJpsiFwdVsTheorCentr020LargeBins = new TCanvas("canvasJpsiFwdVsTheorCentr020LargeBins", "", 800, 600);
    histGridJpsiFwdVsTheor -> Draw();
    lineUnity -> Draw("SAME");
    graTheorFwdCentr020CentrVal -> Draw("L SAME");
    graSystV2JpsiFwdCentr020LargeBins -> Draw("E2P");
    graStatV2JpsiFwdCentr020LargeBins -> Draw("P");
    legendJpsiFwdVsTheorCentr020 -> Draw();

    latexTitle -> DrawLatex(0.20, 0.88, "ALICE Preliminary, OO, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.20, 0.82, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4, 0#minus20%");


    TCanvas *canvasJpsiFwdVsTheorCentr2060 = new TCanvas("canvasJpsiFwdVsTheorCentr2060", "", 800, 600);
    histGridJpsiFwdVsTheor -> Draw();
    lineUnity -> Draw("SAME");
    graTheorFwdCentr2060CentrVal -> Draw("L SAME");
    graSystV2JpsiFwdCentr2060 -> Draw("E2P");
    graStatV2JpsiFwdCentr2060 -> Draw("P");

    TLegend *legendJpsiFwdVsTheorCentr2060 = new TLegend(0.20,0.25,0.45,0.40);
    SetLegend(legendJpsiFwdVsTheorCentr2060);
    legendJpsiFwdVsTheorCentr2060 -> AddEntry(graSystV2JpsiFwdCentr2060LargeBins,"Data, SP, |#Delta#it{#eta}| > 1.7","FP");
    legendJpsiFwdVsTheorCentr2060 -> AddEntry(graTheorFwdCentr2060CentrVal,"Tsinghua Transport + MUSIC Hydro.","L");
    legendJpsiFwdVsTheorCentr2060 -> Draw();

    latexTitle -> DrawLatex(0.20, 0.88, "ALICE Preliminary, OO, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.20, 0.82, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4, 20#minus60%");

    TCanvas *canvasJpsiFwdVsTheorCentr2060LargeBins = new TCanvas("canvasJpsiFwdVsTheorCentr2060LargeBins", "", 800, 600);
    histGridJpsiFwdVsTheor -> Draw();
    lineUnity -> Draw("SAME");
    graTheorFwdCentr2060CentrVal -> Draw("L SAME");
    graSystV2JpsiFwdCentr2060LargeBins -> Draw("E2P");
    graStatV2JpsiFwdCentr2060LargeBins -> Draw("P");
    legendJpsiFwdVsTheorCentr2060 -> Draw();

    latexTitle -> DrawLatex(0.20, 0.88, "ALICE Preliminary, OO, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.20, 0.82, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4, 20#minus60%");
    

    TCanvas *canvasJpsiFwdCentr020LargeEtaGap = new TCanvas("canvasJpsiFwdCentr020LargeEtaGap", "", 800, 600);
    TH2D *histGridJpsiFwdCentr020LargeEtaGap  = new TH2D("histGridJpsiFwdCentr020LargeEtaGap", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}^{SP}", 100, 0, 8, 100, -0.20, 0.30);
    histGridJpsiFwdCentr020LargeEtaGap -> Draw();
    lineUnity -> Draw("SAME");
    //graSystV2JpsiFwdCentr020 -> Draw("E2P");
    //graStatV2JpsiFwdCentr020 -> Draw("P SAME");
    graSystV2JpsiFwdCentr020LargeBins -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr020LargeBins -> Draw("P SAME");
    //graSystV2JpsiFwdCentr020LargeEtaGap -> Draw("E2P SAME");
    //graStatV2JpsiFwdCentr020LargeEtaGap -> Draw("P SAME");
    graSystV2JpsiFwdCentr020LargeBinsLargeEtaGap -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr020LargeBinsLargeEtaGap -> Draw("P SAME");

    TLegend *legendJpsiFwdCentr020LargeEtaGap = new TLegend(0.20,0.25,0.45,0.40);
    SetLegend(legendJpsiFwdCentr020LargeEtaGap);
    legendJpsiFwdCentr020LargeEtaGap -> AddEntry(graSystV2JpsiFwdCentr020LargeBins,"|#Delta#it{#eta}| > 1.7","FP");
    legendJpsiFwdCentr020LargeEtaGap -> AddEntry(graSystV2JpsiFwdCentr020LargeBinsLargeEtaGap,"|#Delta#it{#eta}| > 2.5","FP");
    legendJpsiFwdCentr020LargeEtaGap -> Draw();

    latexTitle -> DrawLatex(0.20, 0.88, "ALICE Preliminary, OO, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.20, 0.82, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4, 0#minus20%");


    TCanvas *canvasJpsiFwdCentr2060LargeEtaGap = new TCanvas("canvasJpsiFwdCentr2060LargeEtaGap", "", 800, 600);
    TH2D *histGridJpsiFwdCentr2060LargeEtaGap  = new TH2D("histGridJpsiFwdCentr2060LargeEtaGap", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}^{SP}", 100, 0, 8, 100, -0.20, 0.30);
    histGridJpsiFwdCentr2060LargeEtaGap -> Draw();
    lineUnity -> Draw("SAME");
    //graSystV2JpsiFwdCentr2060 -> Draw("E2P");
    //graStatV2JpsiFwdCentr2060 -> Draw("P SAME");
    graSystV2JpsiFwdCentr2060LargeBins -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr2060LargeBins -> Draw("P SAME");
    //graSystV2JpsiFwdCentr2060LargeEtaGap -> Draw("E2P SAME");
    //graStatV2JpsiFwdCentr2060LargeEtaGap -> Draw("P SAME");
    graSystV2JpsiFwdCentr2060LargeBinsLargeEtaGap -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr2060LargeBinsLargeEtaGap -> Draw("P SAME");

    TLegend *legendJpsiFwdCentr2060LargeEtaGap = new TLegend(0.20,0.25,0.45,0.40);
    SetLegend(legendJpsiFwdCentr2060LargeEtaGap);
    legendJpsiFwdCentr2060LargeEtaGap -> AddEntry(graSystV2JpsiFwdCentr2060LargeBins,"|#Delta#it{#eta}| > 1.7","FP");
    legendJpsiFwdCentr2060LargeEtaGap -> AddEntry(graSystV2JpsiFwdCentr2060LargeBinsLargeEtaGap,"|#Delta#it{#eta}| > 2.5","FP");
    legendJpsiFwdCentr2060LargeEtaGap -> Draw();

    latexTitle -> DrawLatex(0.20, 0.88, "ALICE Preliminary, OO, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.20, 0.82, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4, 20#minus60%");




    TCanvas *canvasJpsiFwdCentr020EtaGap17vsEtaGap25 = new TCanvas("canvasJpsiFwdCentr020EtaGap17vsEtaGap25", "", 800, 600);
    TH2D *histGridJpsiFwdCentr020EtaGap17vsEtaGap25  = new TH2D("histGridJpsiFwdCentr020EtaGap17vsEtaGap25", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}^{SP}", 100, 0, 8, 100, -0.20, 0.30);
    histGridJpsiFwdCentr020EtaGap17vsEtaGap25 -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiFwdCentr020EtaGap17 -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr020EtaGap17 -> Draw("P SAME");
    graSystV2JpsiFwdCentr020EtaGap25 -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr020EtaGap25 -> Draw("P SAME");

    TLegend *legendJpsiFwdCentr020EtaGap17vsEtaGap25 = new TLegend(0.20,0.25,0.45,0.40);
    SetLegend(legendJpsiFwdCentr020EtaGap17vsEtaGap25);
    legendJpsiFwdCentr020EtaGap17vsEtaGap25 -> AddEntry(graSystV2JpsiFwdCentr020EtaGap17,"|#Delta#it{#eta}| > 1.7","FP");
    legendJpsiFwdCentr020EtaGap17vsEtaGap25 -> AddEntry(graSystV2JpsiFwdCentr020EtaGap25,"|#Delta#it{#eta}| > 2.5","FP");
    legendJpsiFwdCentr020EtaGap17vsEtaGap25 -> Draw();

    TCanvas *canvasJpsiFwdCentr2060EtaGap17vsEtaGap25 = new TCanvas("canvasJpsiFwdCentr2060EtaGap17vsEtaGap25", "", 800, 600);
    TH2D *histGridJpsiFwdCentr2060EtaGap17vsEtaGap25  = new TH2D("histGridJpsiFwdCentr2060EtaGap17vsEtaGap25", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}^{SP}", 100, 0, 8, 100, -0.20, 0.30);
    histGridJpsiFwdCentr2060EtaGap17vsEtaGap25 -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiFwdCentr2060EtaGap17 -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr2060EtaGap17 -> Draw("P SAME");
    graSystV2JpsiFwdCentr2060EtaGap25 -> Draw("E2P SAME");
    graStatV2JpsiFwdCentr2060EtaGap25 -> Draw("P SAME");

    TLegend *legendJpsiFwdCentr2060EtaGap17vsEtaGap25 = new TLegend(0.20,0.25,0.45,0.40);
    SetLegend(legendJpsiFwdCentr2060EtaGap17vsEtaGap25);
    legendJpsiFwdCentr2060EtaGap17vsEtaGap25 -> AddEntry(graSystV2JpsiFwdCentr2060EtaGap17,"|#Delta#it{#eta}| > 1.7","FP");
    legendJpsiFwdCentr2060EtaGap17vsEtaGap25 -> AddEntry(graSystV2JpsiFwdCentr2060EtaGap25,"|#Delta#it{#eta}| > 2.5","FP");
    legendJpsiFwdCentr2060EtaGap17vsEtaGap25 -> Draw();

    latexTitle -> DrawLatex(0.20, 0.88, "ALICE Preliminary, OO, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.20, 0.82, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < #it{y} < 4, 20#minus60%");

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

    // Jpsi Fwd vs Dzero Mid
    TCanvas *canvasJpsiFwdVsDzeroMid = new TCanvas("canvasJpsiFwdVsDzeroMid", "", 800, 600);
    TH2D *histGridJpsiFwdVsDzeroMid  = new TH2D("histGridJpsiFwdVsDzeroMid", ";#it{p}_{T} (GeV/#it{c});#it{#nu}_{2}^{SP}", 100, 0, 8, 100, -0.30, 0.30);
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

    TCanvas *canvasJpsiFwdVsDzeroMidLargeBins = new TCanvas("canvasJpsiFwdVsDzeroMidLargeBins", "", 800, 600);
    histGridJpsiFwdVsDzeroMid -> Draw();
    lineUnity -> Draw("SAME");
    graSystV2JpsiFwdCentr020LargeBins -> Draw("E2P");
    graStatV2JpsiFwdCentr020LargeBins -> Draw("P");
    graSyst1V2DzeroMidCentr020 -> Draw("E2P");
    //graSyst2V2DzeroMidCentr020 -> Draw("E2P");
    graStatV2DzeroMidCentr020 -> Draw("EP");
    legend3 -> Draw();

    latexTitle -> DrawLatex(0.20, 0.88, "ALICE Work in Progress, OO, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    

    canvasJpsiFwdCentr020 -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_centr_0_20.pdf");
    canvasJpsiFwdCentr020vsTimeAssoc -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_centr_0_20_vs_time_assoc.pdf");
    canvasJpsiFwdCentr020vsPbPb -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_centr_0_20_vs_PbPb.pdf");
    canvasJpsiFwdCentr020vsPbPbLargeBins -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_centr_0_20_vs_PbPb_large_bins.pdf");
    canvasJpsiFwdCentr2060 -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_centr_20_60.pdf");
    canvasJpsiFwdCentr2060vsTimeAssoc -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_centr_20_60_vs_time_assoc.pdf");
    canvasJpsiFwdCentr2060vsPbPb -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_centr_20_60_vs_PbPb.pdf");
    canvasJpsiFwdCentr2060vsPbPbLargeBins -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_centr_20_60_vs_PbPb_large_bins.pdf");
    canvasJpsiFwd -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd.pdf");
    canvasJpsiFwdLargeBins -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_large_bins.pdf");
    canvasJpsiFwdVsJpsiMidCentr020 -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_vs_jpsiMid_centr_0_20.pdf");
    canvasJpsiFwdVsJpsiMidCentr020LargeBins -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_vs_jpsiMid_centr_0_20_large_bins.pdf");
    canvasJpsiFwdVsJpsiMidCentr2060 -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_vs_jpsiMid_centr_20_60.pdf");
    canvasJpsiFwdVsJpsiMidCentr2060LargeBins -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_vs_jpsiMid_centr_20_60_large_bins.pdf");
    canvasJpsiFwdVsTheorCentr020 -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_vs_models_centr_0_20.pdf");
    canvasJpsiFwdVsTheorCentr020LargeBins -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_vs_models_centr_0_20_large_bins.pdf");
    canvasJpsiFwdVsTheorCentr2060 -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_vs_models_centr_20_60.pdf");
    canvasJpsiFwdVsTheorCentr2060LargeBins -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_vs_models_centr_20_60_large_bins.pdf");
    canvasJpsiFwdCentr020LargeEtaGap -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_centr_0_20_large_delta_eta.pdf");
    canvasJpsiFwdCentr2060LargeEtaGap -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_centr_20_60_large_delta_eta.pdf");
    canvasJpsiFwdVsDzeroMid -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_vs_dzeroMid.pdf");
    canvasJpsiFwdVsDzeroMidLargeBins -> SaveAs("preliminary_plots/ICHEP2026/jpsiFwd_vs_dzeroMid_large_bins.pdf");

    

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphAsymmErrors* DoGraphFromTheory(string fInName) {
    double xMin, xMax, yCentral, yMin, yMax;

    std::ifstream fInFwd(fInName.c_str());

    std::vector<double> xCentrals, exMins, exMaxs, exZeros, yCentrals, eyMins, eyMaxs, eyZeros;
    while (fInFwd >> xMin >> xMax >> yCentral >> yMin >> yMax) {
        double xCentral = (xMax + xMin) / 2.;
        std::cout << yCentral << std::endl;

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
void SetLegend(TLegend *legend) {
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.045);
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