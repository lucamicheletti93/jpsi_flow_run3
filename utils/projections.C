#include "TProfile.h"
#include <fstream>
#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLatex.h"
#include "TStyle.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TGraphErrors.h"

void LoadStyle();
void SetLegend(TLegend *);

TH1D* ProjectTHnSparse(THnSparseD *, double , double , double , double , double , double );
TH1D* ProjectTH2(TH2D *, double , double );

void projections() {
    //string pathToFiles = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/data/pass3/LHC23_golden/train_231992"; // 20% golden Pb-Pb
    string pathToFiles = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/data/pass3/LHC23_full/train_235036"; // full Pb-Pb
    string histName = "Mass_Pt_centrFT0C_V2";
    string cuts[] = {"muonLowPt210SigmaPDCA"};

    ifstream fRunList(Form("%s/run_list.txt", pathToFiles.c_str()));
    int runNumber;
    vector<int> runList;
    while(!fRunList.eof()) {
        fRunList >> runNumber;
        runList.push_back(runNumber);
    }

    //const int nRuns = 28;
    /*int runList[] = {544013, 544028, 544032, 544091, 544095, 544098, 544116, 544121, 544122, 544123, 
                     544124, 544184, 544185, 544389, 544390, 544391, 544392, 544451, 544454, 544474, 
                     544475, 544476, 544477, 544490, 544491, 544492, 544508, 544510};*/

    const int nPtBins1 = 12;
    double minPtBins1[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20};
    double maxPtBins1[] = {1, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20, 30};
    const int nRapBins1 = 1;
    double minRapBins1[] = {2.5};
    double maxRapBins1[] = {4};
    const int nCentrBins1 = 1;
    double minCentrBins1[] = {0};
    double maxCentrBins1[] = {90};

    const int nPtBins2 = 1;
    double minPtBins2[] = {0};
    double maxPtBins2[] = {30};
    const int nRapBins2 = 1;
    double minRapBins2[] = {2.5};
    double maxRapBins2[] = {4};
    const int nCentrBins2 = 9;
    double minCentrBins2[] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
    double maxCentrBins2[] = {10, 20, 30, 40, 50, 60, 70, 80, 90};

    for (auto const& run : runList) {
        std::cout << run << std::endl;
        TFile *fIn = new TFile(Form("%s/%i/AnalysisResults.root", pathToFiles.c_str(), run));
        TList *listSE = (TList*) fIn -> Get("analysis-same-event-pairing/output");
        TList *listME = (TList*) fIn -> Get("analysis-event-mixing/output");
        for (auto& cut : cuts) {
            TFile *fOut = new TFile(Form("%s/%i/Histograms_%s.root", pathToFiles.c_str(), run, cut.c_str()), "RECREATE");
            if (fIn -> IsZombie()) continue;
            TList *listSEPM = (TList*) listSE -> FindObject(Form("PairsMuonSEPM_%s", cut.c_str()));
            TList *listSEPP = (TList*) listSE -> FindObject(Form("PairsMuonSEPP_%s", cut.c_str()));
            TList *listSEMM = (TList*) listSE -> FindObject(Form("PairsMuonSEMM_%s", cut.c_str()));
            TList *listMEPM = (TList*) listME -> FindObject(Form("PairsMuonMEPM_%s", cut.c_str()));
            TList *listMEPP = (TList*) listME -> FindObject(Form("PairsMuonMEPP_%s", cut.c_str()));
            TList *listMEMM = (TList*) listME -> FindObject(Form("PairsMuonMEMM_%s", cut.c_str()));

            THnSparseD *histSparseSEPM = (THnSparseD*) listSEPM -> FindObject(histName.c_str());
            THnSparseD *histSparseSEPP = (THnSparseD*) listSEPP -> FindObject(histName.c_str());
            THnSparseD *histSparseSEMM = (THnSparseD*) listSEMM -> FindObject(histName.c_str());
            THnSparseD *histSparseMEPM = (THnSparseD*) listMEPM -> FindObject(histName.c_str());
            THnSparseD *histSparseMEPP = (THnSparseD*) listMEPP -> FindObject(histName.c_str());
            THnSparseD *histSparseMEMM = (THnSparseD*) listMEMM -> FindObject(histName.c_str());

            TH1D *histSEPMProjInt = (TH1D*) ProjectTHnSparse(histSparseSEPM, 0, 30, 2.5, 4, 0, 90);
            TH1D *histSEPPProjInt = (TH1D*) ProjectTHnSparse(histSparseSEPP, 0, 30, 2.5, 4, 0, 90);
            TH1D *histSEMMProjInt = (TH1D*) ProjectTHnSparse(histSparseSEMM, 0, 30, 2.5, 4, 0, 90);
            TH1D *histMEPMProjInt = (TH1D*) ProjectTHnSparse(histSparseMEPM, 0, 30, 2.5, 4, 0, 90);
            TH1D *histMEPPProjInt = (TH1D*) ProjectTHnSparse(histSparseMEPP, 0, 30, 2.5, 4, 0, 90);
            TH1D *histMEMMProjInt = (TH1D*) ProjectTHnSparse(histSparseMEMM, 0, 30, 2.5, 4, 0, 90);

            histSEPMProjInt -> Write("Mass_Int_SEPM");
            histSEPPProjInt -> Write("Mass_Int_SEPP");
            histSEMMProjInt -> Write("Mass_Int_SEMM");
            histMEPMProjInt -> Write("Mass_Int_MEPM");
            histMEPPProjInt -> Write("Mass_Int_MEPP");
            histMEMMProjInt -> Write("Mass_Int_MEMM");

            for (int iPt = 0;iPt < nPtBins1;iPt++) {
                TH1D *histSEPMProjPt = (TH1D*) ProjectTHnSparse(histSparseSEPM, minPtBins1[iPt], maxPtBins1[iPt], minRapBins1[0], maxRapBins1[0], minCentrBins1[0], maxCentrBins1[0]);
                TH1D *histSEPPProjPt = (TH1D*) ProjectTHnSparse(histSparseSEPP, minPtBins1[iPt], maxPtBins1[iPt], minRapBins1[0], maxRapBins1[0], minCentrBins1[0], maxCentrBins1[0]);
                TH1D *histSEMMProjPt = (TH1D*) ProjectTHnSparse(histSparseSEMM, minPtBins1[iPt], maxPtBins1[iPt], minRapBins1[0], maxRapBins1[0], minCentrBins1[0], maxCentrBins1[0]);
                TH1D *histMEPMProjPt = (TH1D*) ProjectTHnSparse(histSparseMEPM, minPtBins1[iPt], maxPtBins1[iPt], minRapBins1[0], maxRapBins1[0], minCentrBins1[0], maxCentrBins1[0]);
                TH1D *histMEPPProjPt = (TH1D*) ProjectTHnSparse(histSparseMEPP, minPtBins1[iPt], maxPtBins1[iPt], minRapBins1[0], maxRapBins1[0], minCentrBins1[0], maxCentrBins1[0]);
                TH1D *histMEMMProjPt = (TH1D*) ProjectTHnSparse(histSparseMEMM, minPtBins1[iPt], maxPtBins1[iPt], minRapBins1[0], maxRapBins1[0], minCentrBins1[0], maxCentrBins1[0]);

                histSEPMProjPt -> Write(Form("Mass_Pt_%2.1f_%2.1f_SEPM", minPtBins1[iPt], maxPtBins1[iPt]));
                histSEPPProjPt -> Write(Form("Mass_Pt_%2.1f_%2.1f_SEPP", minPtBins1[iPt], maxPtBins1[iPt]));
                histSEMMProjPt -> Write(Form("Mass_Pt_%2.1f_%2.1f_SEMM", minPtBins1[iPt], maxPtBins1[iPt]));
                histMEPMProjPt -> Write(Form("Mass_Pt_%2.1f_%2.1f_MEPM", minPtBins1[iPt], maxPtBins1[iPt]));
                histMEPPProjPt -> Write(Form("Mass_Pt_%2.1f_%2.1f_MEPP", minPtBins1[iPt], maxPtBins1[iPt]));
                histMEMMProjPt -> Write(Form("Mass_Pt_%2.1f_%2.1f_MEMM", minPtBins1[iPt], maxPtBins1[iPt]));
            }

            for (int iCentr = 0;iCentr < nCentrBins2;iCentr++) {
                TH1D *histSEPMProjCentr = (TH1D*) ProjectTHnSparse(histSparseSEPM, minPtBins2[0], maxPtBins2[0], minRapBins2[0], maxRapBins2[0], minCentrBins2[iCentr], maxCentrBins2[iCentr]);
                TH1D *histSEPPProjCentr = (TH1D*) ProjectTHnSparse(histSparseSEPP, minPtBins2[0], maxPtBins2[0], minRapBins2[0], maxRapBins2[0], minCentrBins2[iCentr], maxCentrBins2[iCentr]);
                TH1D *histSEMMProjCentr = (TH1D*) ProjectTHnSparse(histSparseSEMM, minPtBins2[0], maxPtBins2[0], minRapBins2[0], maxRapBins2[0], minCentrBins2[iCentr], maxCentrBins2[iCentr]);
                TH1D *histMEPMProjCentr = (TH1D*) ProjectTHnSparse(histSparseMEPM, minPtBins2[0], maxPtBins2[0], minRapBins2[0], maxRapBins2[0], minCentrBins2[iCentr], maxCentrBins2[iCentr]);
                TH1D *histMEPPProjCentr = (TH1D*) ProjectTHnSparse(histSparseMEPP, minPtBins2[0], maxPtBins2[0], minRapBins2[0], maxRapBins2[0], minCentrBins2[iCentr], maxCentrBins2[iCentr]);
                TH1D *histMEMMProjCentr = (TH1D*) ProjectTHnSparse(histSparseMEMM, minPtBins2[0], maxPtBins2[0], minRapBins2[0], maxRapBins2[0], minCentrBins2[iCentr], maxCentrBins2[iCentr]);

                histSEPMProjCentr -> Write(Form("Mass_CentrFT0C_%2.1f_%2.1f_SEPM", minCentrBins2[iCentr], maxCentrBins2[iCentr]));
                histSEPPProjCentr -> Write(Form("Mass_CentrFT0C_%2.1f_%2.1f_SEPP", minCentrBins2[iCentr], maxCentrBins2[iCentr]));
                histSEMMProjCentr -> Write(Form("Mass_CentrFT0C_%2.1f_%2.1f_SEMM", minCentrBins2[iCentr], maxCentrBins2[iCentr]));
                histMEPMProjCentr -> Write(Form("Mass_CentrFT0C_%2.1f_%2.1f_MEPM", minCentrBins2[iCentr], maxCentrBins2[iCentr]));
                histMEPPProjCentr -> Write(Form("Mass_CentrFT0C_%2.1f_%2.1f_MEPP", minCentrBins2[iCentr], maxCentrBins2[iCentr]));
                histMEMMProjCentr -> Write(Form("Mass_CentrFT0C_%2.1f_%2.1f_MEMM", minCentrBins2[iCentr], maxCentrBins2[iCentr]));
            }

            fOut -> Close();
        }
    }
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
////////////////////////////////////////////////////////////////////////////////
TH1D* ProjectTHnSparse(THnSparseD *histSparse, double minVar1, double maxVar1, double minVar2, double maxVar2, double minVar3, double maxVar3) {
    double minVar1Bin = histSparse -> GetAxis(1) -> FindBin(minVar1);
    double maxVar1Bin = histSparse -> GetAxis(1) -> FindBin(maxVar1 - 0.01);
    double minVar2Bin = histSparse -> GetAxis(2) -> FindBin(minVar2);
    double maxVar2Bin = histSparse -> GetAxis(2) -> FindBin(maxVar2 - 0.01);
    double minVar3Bin = histSparse -> GetAxis(3) -> FindBin(minVar3);
    double maxVar3Bin = histSparse -> GetAxis(3) -> FindBin(maxVar3 - 0.01);

    histSparse -> GetAxis(1) -> SetRange(minVar1Bin, maxVar1Bin);
    histSparse -> GetAxis(2) -> SetRange(minVar2Bin, maxVar2Bin);
    histSparse -> GetAxis(3) -> SetRange(minVar3Bin, maxVar3Bin);

    TH1D *histProj = (TH1D*) histSparse -> Projection(0, Form("histProj_%1.0f_%1.0f__%1.0f_%1.0f__%1.0f_%1.0f", minVar1, maxVar1, minVar2, maxVar2, minVar3, maxVar3));
    return histProj;
}
////////////////////////////////////////////////////////////////////////////////
TH1D* ProjectTH2(TH2D *hist2D, double minPtRange, double maxPtRange) {
    double minPtBin = hist2D -> GetYaxis() -> FindBin(minPtRange);
    double maxPtBin = hist2D -> GetYaxis() -> FindBin(maxPtRange - 0.01);

    Printf("minPtBin = %1.0f, maxPtBin = %1.0f", minPtBin, maxPtBin);
    hist2D -> GetYaxis() -> SetRange(minPtBin, maxPtBin);

    TH1D *histProj = (TH1D*) hist2D -> ProjectionX(Form("histProj_%1.0f_%1.0f", minPtRange, maxPtRange));
    return histProj;
}