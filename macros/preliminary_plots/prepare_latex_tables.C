inline void GetDataFromTxt(string fInName, 
    std::vector<double>& vecCentrX, std::vector<double>& vecStatX, std::vector<double>& vecSystX, 
    std::vector<double>& vecCentrY, std::vector<double>& vecStatY, std::vector<double>& vecSystY, double systWidth) {

std::ifstream fIn(fInName.c_str());
if (!fIn) { std::cerr << "Error while opening fIn!" << std::endl;}

double minX, maxX, centrY, statY, systY;
while (fIn >> minX >> maxX >> centrY >> statY >> systY) {
vecCentrX.push_back((maxX + minX) / 2.);
vecStatX.push_back((maxX - minX) / 2.);
vecSystX.push_back(systWidth);
vecCentrY.push_back(centrY);
vecStatY.push_back(statY);
vecSystY.push_back(systY);
}

fIn.close();
}

void prepare_latex_tables() {
    // EP, List of input files
    string pathNarrowBinsEpPt1030 = "../../results/QM2025/EP/jpsi_v2_vs_pt_10_30_narrow_binning.txt";
    string pathNarrowBinsEpPt3050 = "../../results/QM2025/EP/jpsi_v2_vs_pt_30_50_narrow_binning.txt";
    string pathNarrowBinsEpPt5080 = "../../results/QM2025/EP/jpsi_v2_vs_pt_50_80_large_binning.txt";
    string pathLargeBinsEpPt1030 = "../../results/QM2025/EP/jpsi_v2_vs_pt_10_30_large_binning.txt";
    string pathLargeBinsEpPt3050 = "../../results/QM2025/EP/jpsi_v2_vs_pt_30_50_large_binning.txt";
    string pathLargeBinsEpPt5080 = "../../results/QM2025/EP/jpsi_v2_vs_pt_50_80_large_binning.txt";
    string pathEpCentr05 = "../../results/QM2025/EP/jpsi_v2_vs_centr_0_5.txt";
    string pathEpCentr515 = "../../results/QM2025/EP/jpsi_v2_vs_centr_5_15.txt";

    // v2 vs pT
    // Narrow Binning
    std::vector<double> ptCentr1030EpRun3, ptWidthStat1030EpRun3, ptWidthSyst1030EpRun3, v2JpsiVsPt1030EpRun3, statV2JpsiVsPt1030EpRun3, systV2JpsiVsPt1030EpRun3;
    GetDataFromTxt(pathNarrowBinsEpPt1030.c_str(), ptCentr1030EpRun3, ptWidthStat1030EpRun3, ptWidthSyst1030EpRun3, v2JpsiVsPt1030EpRun3, statV2JpsiVsPt1030EpRun3, systV2JpsiVsPt1030EpRun3, 0.15);

    std::vector<double> ptCentr3050EpRun3, ptWidthStat3050EpRun3, ptWidthSyst3050EpRun3, v2JpsiVsPt3050EpRun3, statV2JpsiVsPt3050EpRun3, systV2JpsiVsPt3050EpRun3;
    GetDataFromTxt(pathNarrowBinsEpPt3050.c_str(), ptCentr3050EpRun3, ptWidthStat3050EpRun3, ptWidthSyst3050EpRun3, v2JpsiVsPt3050EpRun3, statV2JpsiVsPt3050EpRun3, systV2JpsiVsPt3050EpRun3, 0.15);

    // Large Binning
    std::vector<double> ptCentr1030EpRun3LargeBins, ptWidthStat1030EpRun3LargeBins, ptWidthSyst1030EpRun3LargeBins, v2JpsiVsPt1030EpRun3LargeBins, statV2JpsiVsPt1030EpRun3LargeBins, systV2JpsiVsPt1030EpRun3LargeBins;
    GetDataFromTxt(pathLargeBinsEpPt1030.c_str(), ptCentr1030EpRun3LargeBins, ptWidthStat1030EpRun3LargeBins, ptWidthSyst1030EpRun3LargeBins, v2JpsiVsPt1030EpRun3LargeBins, statV2JpsiVsPt1030EpRun3LargeBins, systV2JpsiVsPt1030EpRun3LargeBins, 0.15);

    std::vector<double> ptCentr3050EpRun3LargeBins, ptWidthStat3050EpRun3LargeBins, ptWidthSyst3050EpRun3LargeBins, v2JpsiVsPt3050EpRun3LargeBins, statV2JpsiVsPt3050EpRun3LargeBins, systV2JpsiVsPt3050EpRun3LargeBins;
    GetDataFromTxt(pathLargeBinsEpPt3050.c_str(), ptCentr3050EpRun3LargeBins, ptWidthStat3050EpRun3LargeBins, ptWidthSyst3050EpRun3LargeBins, v2JpsiVsPt3050EpRun3LargeBins, statV2JpsiVsPt3050EpRun3LargeBins, systV2JpsiVsPt3050EpRun3LargeBins, 0.15);

    std::vector<double> ptCentr5080EpRun3LargeBins, ptWidthStat5080EpRun3LargeBins, ptWidthSyst5080EpRun3LargeBins, v2JpsiVsPt5080EpRun3LargeBins, statV2JpsiVsPt5080EpRun3LargeBins, systV2JpsiVsPt5080EpRun3LargeBins;
    GetDataFromTxt(pathLargeBinsEpPt5080.c_str(), ptCentr5080EpRun3LargeBins, ptWidthStat5080EpRun3LargeBins, ptWidthSyst5080EpRun3LargeBins, v2JpsiVsPt5080EpRun3LargeBins, statV2JpsiVsPt5080EpRun3LargeBins, systV2JpsiVsPt5080EpRun3LargeBins, 0.15);

    std::cout << "*** Narrow pT binning ***" << std::endl;
    for (int iVar = 0;iVar < ptCentr1030EpRun3.size();iVar++) {
        double minX = ptCentr1030EpRun3[iVar] - ptWidthStat1030EpRun3[iVar];
        double maxX = ptCentr1030EpRun3[iVar] + ptWidthStat1030EpRun3[iVar];
        Printf(
            "%3.2f - %3.2f & %4.3f $?pm$ %4.3f $?pm$ %4.3f & %4.3f $?pm$ %4.3f $?pm$ %4.3f ??", minX, maxX, 
            v2JpsiVsPt1030EpRun3[iVar], statV2JpsiVsPt1030EpRun3[iVar], systV2JpsiVsPt1030EpRun3[iVar],
            v2JpsiVsPt3050EpRun3[iVar], statV2JpsiVsPt3050EpRun3[iVar], systV2JpsiVsPt3050EpRun3[iVar]
        );
    }
    std::cout << "*** Large pT binning ***" << std::endl;
    for (int iVar = 0;iVar < ptCentr1030EpRun3LargeBins.size();iVar++) {
        double minX = ptCentr1030EpRun3LargeBins[iVar] - ptWidthStat1030EpRun3LargeBins[iVar];
        double maxX = ptCentr1030EpRun3LargeBins[iVar] + ptWidthStat1030EpRun3LargeBins[iVar];
        Printf(
            "%3.2f - %3.2f & %4.3f $?pm$ %4.3f $?pm$ %4.3f & %4.3f $?pm$ %4.3f $?pm$ %4.3f & %4.3f $?pm$ %4.3f $?pm$ %4.3f ??", minX, maxX, 
            v2JpsiVsPt1030EpRun3LargeBins[iVar], statV2JpsiVsPt1030EpRun3LargeBins[iVar], systV2JpsiVsPt1030EpRun3LargeBins[iVar],
            v2JpsiVsPt3050EpRun3LargeBins[iVar], statV2JpsiVsPt3050EpRun3LargeBins[iVar], systV2JpsiVsPt3050EpRun3LargeBins[iVar],
            v2JpsiVsPt5080EpRun3LargeBins[iVar], statV2JpsiVsPt5080EpRun3LargeBins[iVar], systV2JpsiVsPt5080EpRun3LargeBins[iVar]
        );
    }

    // v2 vs centrality
    std::vector<double> centrCenCentr05EpRun3, ptWidthStCentr05EpRun3, ptWidthSyCentr05EpRun3, v2JpsiVsCentr05EpRun3, statV2JpsiVsCentr05EpRun3, systV2JpsiVsCentr05EpRun3;
    GetDataFromTxt(pathEpCentr05.c_str(), centrCenCentr05EpRun3, ptWidthStCentr05EpRun3, ptWidthSyCentr05EpRun3, v2JpsiVsCentr05EpRun3, statV2JpsiVsCentr05EpRun3, systV2JpsiVsCentr05EpRun3, 1.5);

    std::vector<double> centrCenCentr515EpRun3, ptWidthStCentr515EpRun3, ptWidthSyCentr515EpRun3, v2JpsiVsCentr515EpRun3, statV2JpsiVsCentr515EpRun3, systV2JpsiVsCentr515EpRun3;
    GetDataFromTxt(pathEpCentr515.c_str(), centrCenCentr515EpRun3, ptWidthStCentr515EpRun3, ptWidthSyCentr515EpRun3, v2JpsiVsCentr515EpRun3, statV2JpsiVsCentr515EpRun3, systV2JpsiVsCentr515EpRun3, 1.5);

    std::cout << "*** Narrow centrality binning ***" << std::endl;
    for (int iVar = 0;iVar < centrCenCentr05EpRun3.size();iVar++) {
        double minX = centrCenCentr05EpRun3[iVar] - ptWidthStCentr05EpRun3[iVar];
        double maxX = centrCenCentr05EpRun3[iVar] + ptWidthStCentr05EpRun3[iVar];
        Printf(
            "%3.2f - %3.2f & %4.3f $?pm$ %4.3f $?pm$ %4.3f & %4.3f $?pm$ %4.3f $?pm$ %4.3f ??", minX, maxX, 
            v2JpsiVsCentr05EpRun3[iVar], statV2JpsiVsCentr05EpRun3[iVar], systV2JpsiVsCentr05EpRun3[iVar],
            v2JpsiVsCentr515EpRun3[iVar], statV2JpsiVsCentr515EpRun3[iVar], systV2JpsiVsCentr515EpRun3[iVar]
        );
    }
}