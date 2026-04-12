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

TH1D* projectHistogram(THnSparseF *, double , double , double , double );
TProfile* projectProfile(THnSparseF *, double , double , double , double , int);

void plot_dimuon_dist()
{
   bool integrated = false;
   //string pathToFiles = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/data/pass2/LHC23_golden/dimuon_train_198686_with_cuts";
    string pathToFiles = "/Users/dhananjaya/Desktop/ALICE_PbPb_Analysis/jpsi_flow_run3/macros";
    //string pathToFiles = "/afs/cern.ch/user/d/dhthakur/pass4_PbPbQuality";
    //string muonCut = "muonLowPt10SigmaPDCA";
    //string muonCut = "muonLowPt210SigmaPDCA";
    string muonCut = "matchedMchMid";
    string suffix = "";
    //string suffix = "_new_EP_calculation";
    string methodName = "SP";
    int method = 0;

    if (methodName == "EP") {
        method = 5;
    } else {
        method = 4;
    }

    //TFile *fIn  = new TFile(Form("%s%s/AnalysisResults.root", pathToFiles.c_str(), suffix.c_str()), "READ");
    //TFile *fIn  = new TFile("AnalysisResults_Fullpass3_singlerun.root", "READ");
    TFile *fIn2  = new TFile("AnalysisResults_OO_pass2_SE_eta05.root", "READ");
    fIn2->ls();
    
    // Same event pairing
    TList *list1 = (TList*) fIn2 -> Get("analysis-same-event-pairing/output");
    TList *listSEPM = (TList*) list1 -> FindObject(Form("PairsMuonSEPM_%s", muonCut.c_str()));
    THnSparseF* histSEPM = (THnSparseF*) listSEPM->FindObject("Mass_Pt_centrFT0C_V2"); // 0: mass; 1: pT ; 2:FT0C_Centrality; 3:u2q2_A; 4:cos2(phi-psi)

    TList *listSEPP = (TList*) list1 -> FindObject(Form("PairsMuonSEPP_%s", muonCut.c_str()));
    THnSparseF* histSEPP = (THnSparseF*) listSEPP -> FindObject("Mass_Pt_centrFT0C_V2");

    TList *listSEMM = (TList*) list1 -> FindObject(Form("PairsMuonSEMM_%s", muonCut.c_str()));
    THnSparseF* histSEMM = (THnSparseF*) listSEMM -> FindObject("Mass_Pt_centrFT0C_V2");
    

    //TFile *fIn3  = new TFile("AnalysisResults_ME123456789.root", "READ");
    TFile *fIn3  = new TFile("AnalysisResults_OO_pass2_ME_eta05.root", "READ");
    // Mixed event pairing
    TList *list2 = (TList*) fIn3 -> Get("analysis-event-mixing/output");
    TList *listMEPM = (TList*) list2 -> FindObject(Form("PairsMuonMEPM_%s", muonCut.c_str()));
    THnSparseF* histMEPM = (THnSparseF*) listMEPM->FindObject("Mass_Pt_centrFT0C_V2");

    TList *listMEPP = (TList*) list2 -> FindObject(Form("PairsMuonMEPP_%s", muonCut.c_str()));
    THnSparseF* histMEPP = (THnSparseF*) listMEPP -> FindObject("Mass_Pt_centrFT0C_V2");

    TList *listMEMM = (TList*) list2 -> FindObject(Form("PairsMuonMEMM_%s", muonCut.c_str()));
    THnSparseF* histMEMM = (THnSparseF*) listMEMM -> FindObject("Mass_Pt_centrFT0C_V2");
    
    // v2 and pT distributions
    if (!integrated) {
      //const int nPtBins = 12;
        //double minCentrBins[] = {10,0,5,10,20,30,40,50,60,70};
        //double maxCentrBins[] = {50,5,10,20,30,40,50,60,70,80};

      //int nCent = 5;
      //double minCentrBins[] = {10,30,50,10,10};
      //double maxCentrBins[] = {30,50,80,50,80};

      int nCent = 1;
      double minCentrBins[] = {0};
      double maxCentrBins[] = {20};
      
      //int nCent = 26;
      //double minCentrBins[] = {0,10,20,30,40,50,60,70,0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85};
      //double maxCentrBins[] = {10,20,30,40,50,60,70,80,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90};

      //const int nPtBins = 8;
      //double minPtBins[] = {0, 0, 1, 2, 3, 4, 6,8};
      //double maxPtBins[] = {10, 1, 2, 3, 4, 6, 8,15};

      //const int nPtBins = 13;
      //double minPtBins[] = {0,1,2,3,4,5,6,8.,10.,12.,15.,0.,10.};
      //double maxPtBins[] = {1,2,3,4,5,6,8,10.,12.,15.,20.,2.,15.};

      const int nPtBins = 8;
      double minPtBins[] = {0.0,6.0,0.0,0.0, 2.0, 4.0, 6.0, 8.0};
      double maxPtBins[] = {6.0,15.0,15.0,2.0, 4.0, 6.0, 8.0, 15.0};


      //const int nPtBins = 20;
      //double minPtBins[] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,0,1,2,3,4,5,6,8.,10.,12.,15.};
      //double maxPtBins[] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,1,2,3,4,5,6,8,10.,12.,15.,20.};

	//const int nPtBins = 17;
	//double minPtBins[] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,8.0,10.0,12.0,15.0};
	//double maxPtBins[] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,8.0,10.0,12.0,15.0,20.0};

      //const int nPtBins = 9;
      //double minPtBins[] = {0,1,2,3,4,5,6,8,10};
      //double maxPtBins[] = {1,2,3,4,5,6,8,10,15};
          
	TFile *fOut = new TFile(Form("%s/Histograms_O_O_%s_CentBins_MchMid_%1.f_%1.f_SP.root", pathToFiles.c_str(), muonCut.c_str(), minCentrBins[0], maxCentrBins[0]), "RECREATE");
	for (int icent = 0;icent < nCent;icent++)
	  {	   
	    //TFile *fOut = new TFile(Form("%s/Histograms_%s_centr_%1.f_%1.f.root", pathToFiles.c_str(), muonCut.c_str(), minCentrBins[icent], maxCentrBins[icent]), "RECREATE");
        for (int iPt = 0;iPt < nPtBins;iPt++) {
            TH1D *histMassSEPM = (TH1D*) projectHistogram(histSEPM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]);
            histMassSEPM -> SetName(Form("histMassSEPM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
	    cout<<" Get entries ==========:::::: "<<histMassSEPM->GetEntries()<<endl;
	    TH1D *histMassSEPP = (TH1D*) projectHistogram(histSEPP, minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]);
            histMassSEPP -> SetName(Form("histMassSEPP_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
	    TH1D *histMassSEMM = (TH1D*) projectHistogram(histSEMM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]);
            histMassSEMM -> SetName(Form("histMassSEMM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
	    
	    TH1D *histMassMEPM = (TH1D*) projectHistogram(histMEPM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]);
            histMassMEPM -> SetName(Form("histMassMEPM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
	    TH1D *histMassMEPP = (TH1D*) projectHistogram(histMEPP, minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]);
            histMassMEPP -> SetName(Form("histMassMEPP_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
	    TH1D *histMassMEMM = (TH1D*) projectHistogram(histMEMM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]);
            histMassMEMM -> SetName(Form("histMassMEMM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));

	    TProfile *histV2SEPM = (TProfile*) projectProfile(histSEPM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent], method);
            histV2SEPM -> SetName(Form("histV2SEPM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
            TProfile *histV2MEPM = (TProfile*) projectProfile(histMEPM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent], method);
            histV2MEPM -> SetName(Form("histV2MEPM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));

	    TProfile *histV2SEPP = (TProfile*) projectProfile(histSEPP, minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent], method);
            histV2SEPP -> SetName(Form("histV2SEPP_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
            TProfile *histV2MEPP = (TProfile*) projectProfile(histMEPP, minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent], method);
            histV2MEPP -> SetName(Form("histV2MEPP_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));

	    TProfile *histV2SEMM = (TProfile*) projectProfile(histSEMM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent], method);
            histV2SEMM -> SetName(Form("histV2SEMM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
            TProfile *histV2MEMM = (TProfile*) projectProfile(histMEMM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent], method);
            histV2MEMM -> SetName(Form("histV2MEMM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
	    
            TProfile *histPtSEPM = (TProfile*) projectProfile(histSEPM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent], 1);
            histPtSEPM -> SetName(Form("histPtSEPM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
            TProfile *histPtMEPM = (TProfile*) projectProfile(histMEPM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent], 1);
            histPtMEPM -> SetName(Form("histPtMEPM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));


	    TProfile *histU2Q2SEPM = (TProfile*) projectProfile(histSEPM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent], 4);
            histU2Q2SEPM -> SetName(Form("histU2Q2SEPM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
            TProfile *histU2Q2MEPM = (TProfile*) projectProfile(histMEPM, minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent], 4);
            histU2Q2MEPM -> SetName(Form("histU2Q2MEPM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));

            histMassSEPM -> SetTitle(Form("SE - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
            histMassMEPM -> SetTitle(Form("ME - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
            histV2SEPM -> SetTitle(Form("V2 SE - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
            histV2MEPM -> SetTitle(Form("V2 ME - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
	    histV2SEPP -> SetTitle(Form("V2 SEPP - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
            histV2MEPP -> SetTitle(Form("V2 MEPP - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
	    histV2SEMM -> SetTitle(Form("V2 SEMM - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
            histV2MEMM -> SetTitle(Form("V2 MEMM - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
            histPtSEPM -> SetTitle(Form("<pT> SE - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
            histPtMEPM -> SetTitle(Form("<pT> ME - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
	    histU2Q2SEPM -> SetTitle(Form("U2Q2 SE - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));
            histU2Q2MEPM -> SetTitle(Form("U2Q2 ME - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], minCentrBins[icent], maxCentrBins[icent]));

	    fOut -> cd();
            histMassSEPM -> Write();
	    histMassSEPP -> Write();
	    histMassSEMM -> Write();   
            histMassMEPM -> Write();
	    histMassMEPP -> Write();
	    histMassMEMM -> Write();

	    histV2SEPM -> Write();
            histV2MEPM -> Write();
	    histV2SEPP -> Write();
            histV2MEPP -> Write();
	    histV2SEMM -> Write();
            histV2MEMM -> Write();
            histPtSEPM -> Write();
            histPtMEPM -> Write();
	    histU2Q2SEPM -> Write();
            histU2Q2MEPM -> Write();
	    
	    
            histMassSEPM -> SetLineColor(kBlack);
            histMassSEPM -> SetMarkerStyle(20);
            histMassSEPM -> SetMarkerSize(0.8);
            histMassSEPM -> SetMarkerColor(kBlack);

            histV2SEPM -> SetLineColor(kBlack);
            histV2SEPM -> SetMarkerStyle(20);
            histV2SEPM -> SetMarkerSize(0.8);
            histV2SEPM -> SetMarkerColor(kBlack);

	    histU2Q2SEPM-> SetLineColor(kBlack);
	    histU2Q2SEPM-> SetMarkerStyle(20);
	    histU2Q2SEPM-> SetMarkerSize(0.8);
	    histU2Q2SEPM-> SetMarkerColor(kBlack);

            histPtSEPM -> SetLineColor(kBlack);
            histPtSEPM -> SetMarkerStyle(20);
            histPtSEPM -> SetMarkerSize(0.8);
            histPtSEPM -> SetMarkerColor(kBlack);

            histMassMEPM -> SetLineColor(kAzure+4);
            histV2MEPM -> SetLineColor(kAzure+4);
            histPtMEPM -> SetLineColor(kAzure+4);
	    histU2Q2MEPM -> SetLineColor(kAzure+4);

            Printf("Compute normalization...\n");
            double minMassBin1 = histMassSEPM -> GetXaxis() -> FindBin(2.00);
            double maxMassBin1 = histMassSEPM -> GetXaxis() -> FindBin(2.50);
            double minMassBin2 = histMassSEPM -> GetXaxis() -> FindBin(3.72);
            double maxMassBin2 = histMassSEPM -> GetXaxis() -> FindBin(5.00);

            double normMassSEPM = histMassSEPM -> Integral(minMassBin1, maxMassBin1) + histMassSEPM->Integral(minMassBin2, maxMassBin2);
            double normMassMEPM = histMassMEPM -> Integral(minMassBin1, maxMassBin1) + histMassMEPM->Integral(minMassBin2, maxMassBin2);

            Printf("Mass normalization factor = %f\n", normMassSEPM / normMassMEPM);


            TCanvas *canvasMassV2Pt = new TCanvas("canvasMassV2Pt", "", 600, 1800);
            canvasMassV2Pt -> Divide(1, 4);

            canvasMassV2Pt -> cd(1);
            histMassSEPM -> GetXaxis() -> SetRangeUser(2, 5);
            histMassMEPM -> GetXaxis() -> SetRangeUser(2, 5);
            histMassSEPM -> Draw("EP");
            //histMassMEPM -> Scale(normMassSEPM / normMassMEPM);
            histMassMEPM -> Draw("EP SAME");

            TH1D *histBkgSubtrSEMP = (TH1D*) histMassSEPM -> Clone(Form("histMassBkgSubtrSEPM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
            histBkgSubtrSEMP -> Add(histMassMEPM, -1.);
            histBkgSubtrSEMP -> SetLineColor(kRed+1);
            histBkgSubtrSEMP -> SetMarkerStyle(24);
            histBkgSubtrSEMP -> SetMarkerColor(kRed+1);
            histBkgSubtrSEMP -> Draw("EP SAME");

            canvasMassV2Pt -> cd(2);
            histV2SEPM -> GetXaxis() -> SetRangeUser(2, 5);
            histV2SEPM -> Draw("EP");
            histV2MEPM -> Scale(histV2SEPM -> GetBinContent(histV2SEPM -> GetXaxis() -> FindBin(2)) / histV2MEPM -> GetBinContent(histV2MEPM -> GetXaxis() -> FindBin(2)));
            histV2MEPM -> Draw("H SAME");

            canvasMassV2Pt -> cd(3);
            histPtSEPM -> GetXaxis() -> SetRangeUser(2, 5);
            histPtSEPM -> Draw("EP");
            histPtMEPM -> Draw("H SAME");

	    canvasMassV2Pt -> cd(4);
	    histU2Q2SEPM -> GetXaxis() -> SetRangeUser(2, 5);
            histU2Q2SEPM -> Draw("EP");

            canvasMassV2Pt -> SaveAs(Form("LHC23_golden/flow/%s_distrib_centr_%s_%0.1f_%0.1f__%1.0f_%1.0f%s.pdf", muonCut.c_str(), methodName.c_str(), minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0], suffix.c_str()));
        
       
        }  //Pt Loop
      
	} //Centrality loop
	fOut -> Close();
    }

    if (integrated) {
        TFile *fOut = new TFile(Form("%s/Histograms_%s.root", pathToFiles.c_str(), muonCut.c_str()), "RECREATE");
        // pT scan
        //const int nPtBins = 9;
        //double minPtBins[] = {0, 1, 2, 3, 5, 7, 10, 15, 20};
        //double maxPtBins[] = {1, 2, 3, 5, 7, 10, 15, 20, 30};
        const int nPtBins = 4;
        double minPtBins[] = {0,0, 2, 4, 6};
        double maxPtBins[] = {100,2, 4, 6, 12};

        TCanvas *canvasMassPtScan = new TCanvas("canvasMassPtScan", "", 1800, 1800);
        canvasMassPtScan -> Divide(3, 3);

        for (int iPt = 0;iPt < nPtBins;iPt++) {
            TH1D *histMassSEPM = (TH1D*) projectHistogram(histSEPM, minPtBins[iPt], maxPtBins[iPt], 0., 90.);
            TH1D *histMassMEPM = (TH1D*) projectHistogram(histMEPM, minPtBins[iPt], maxPtBins[iPt], 0., 90.);

            histMassSEPM -> SetName(Form("histMassSEPM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], 0., 90.));
            histMassMEPM -> SetName(Form("histMassMEPM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], 0., 90.));
            histMassSEPM -> SetTitle(Form("SE - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], 0., 90.));
            histMassMEPM -> SetTitle(Form("ME - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", minPtBins[iPt], maxPtBins[iPt], 0., 90.));

            histMassSEPM -> SetLineColor(kBlack);
            histMassSEPM -> SetMarkerStyle(20);
            histMassSEPM -> SetMarkerSize(0.8);
            histMassSEPM -> SetMarkerColor(kBlack);
            histMassMEPM -> SetLineColor(kAzure+4);

            Printf("Compute normalization...\n");
            double minMassBin1 = histMassSEPM -> GetXaxis() -> FindBin(1.00);
            double maxMassBin1 = histMassSEPM -> GetXaxis() -> FindBin(2.50);
            double minMassBin2 = histMassSEPM -> GetXaxis() -> FindBin(3.72);
            double maxMassBin2 = histMassSEPM -> GetXaxis() -> FindBin(5.00);

            double normMassSEPM = histMassSEPM -> Integral(minMassBin1, maxMassBin1) + histMassSEPM->Integral(minMassBin2, maxMassBin2);
            double normMassMEPM = histMassMEPM -> Integral(minMassBin1, maxMassBin1) + histMassMEPM->Integral(minMassBin2, maxMassBin2);

            Printf("Normalization factor = %f\n", normMassSEPM / normMassMEPM);

            canvasMassPtScan -> cd(iPt+1);
            gPad -> SetLogy(1);
            histMassSEPM -> GetXaxis() -> SetRangeUser(2, 5);
            histMassSEPM -> GetYaxis() -> SetRangeUser(1, 1e7);
            histMassMEPM -> GetXaxis() -> SetRangeUser(2, 5);
            histMassSEPM -> Draw("EP");
            histMassMEPM -> Scale(normMassSEPM / normMassMEPM);
            histMassMEPM -> Draw("EP SAME");

            TH1D *histBkgSubtrSEMP = (TH1D*) histMassSEPM -> Clone(Form("histMassBkgSubtrSEPM_%0.1f_%0.1f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], 0., 90.));
            histBkgSubtrSEMP -> Add(histMassMEPM, -1.);
            histBkgSubtrSEMP -> SetLineColor(kRed+1);
            histBkgSubtrSEMP -> SetMarkerStyle(24);
            histBkgSubtrSEMP -> SetMarkerSize(0.8);
            histBkgSubtrSEMP -> SetMarkerColor(kRed+1);
            histBkgSubtrSEMP -> Draw("EP SAME");

            fOut -> cd();
            histMassSEPM -> Write();
            histMassMEPM -> Write();
            histBkgSubtrSEMP -> Write();
        }

        // Centrality scan
        const int nCentrBins = 9;
        double minCentrBins[] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
        double maxCentrBins[] = {10, 20, 30, 40, 50, 60, 70, 80, 90};

        TCanvas *canvasCentrScan = new TCanvas("canvasCentrScan", "", 1800, 1800);
        canvasCentrScan -> Divide(3, 3);

        TCanvas *canvasV2CentrScan = new TCanvas("canvasV2CentrScan", "", 1800, 1800);
        canvasV2CentrScan -> Divide(3, 3);


        for (int iCentr = 0;iCentr < nCentrBins;iCentr++) {
            TH1D *histMassSEPM = (TH1D*) projectHistogram(histSEPM, 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]);
            histMassSEPM -> SetName(Form("histMassSEPM_%0.1f_%0.1f__%1.0f_%1.0f", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));
            TH1D *histMassMEPM = (TH1D*) projectHistogram(histMEPM, 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]);
            histMassMEPM -> SetName(Form("histMassMEPM_%0.1f_%0.1f__%1.0f_%1.0f", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));
            TProfile *histV2SEPM = (TProfile*) projectProfile(histSEPM, 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr], method);
            histV2SEPM -> SetName(Form("histV2SEPM_%0.1f_%0.1f__%1.0f_%1.0f", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));
            TProfile *histV2MEPM = (TProfile*) projectProfile(histMEPM, 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr], method);
            histV2MEPM -> SetName(Form("histV2MEPM_%0.1f_%0.1f__%1.0f_%1.0f", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));

            histMassSEPM -> SetTitle(Form("SE - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));
            histMassMEPM -> SetTitle(Form("ME - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));
            histV2SEPM -> SetTitle(Form("V2 SE - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));
            histV2MEPM -> SetTitle(Form("V2 ME - %0.1f < #it{p}_{T} < %0.1f GeV/#it{c}, %1.0f #minus %1.0f %%", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));

            histMassSEPM -> SetLineColor(kBlack);
            histMassSEPM -> SetMarkerStyle(20);
            histMassSEPM -> SetMarkerSize(0.8);
            histMassSEPM -> SetMarkerColor(kBlack);
            histMassMEPM -> SetLineColor(kAzure+4);

            histV2SEPM -> SetLineColor(kBlack);
            histV2SEPM -> SetMarkerStyle(20);
            histV2SEPM -> SetMarkerSize(0.8);
            histV2SEPM -> SetMarkerColor(kBlack);
            histV2MEPM -> SetLineColor(kAzure+4);

            Printf("Compute normalization...\n");
            double minMassBin1 = histMassSEPM -> GetXaxis() -> FindBin(1.00);
            double maxMassBin1 = histMassSEPM -> GetXaxis() -> FindBin(2.50);
            double minMassBin2 = histMassSEPM -> GetXaxis() -> FindBin(3.72);
            double maxMassBin2 = histMassSEPM -> GetXaxis() -> FindBin(5.00);

            double minV2Bin1 = histV2SEPM -> GetXaxis() -> FindBin(1.00);
            double maxV2Bin1 = histV2SEPM -> GetXaxis() -> FindBin(2.50);
            double minV2Bin2 = histV2SEPM -> GetXaxis() -> FindBin(3.72);
            double maxV2Bin2 = histV2SEPM -> GetXaxis() -> FindBin(5.00);

            double normMassSEPM = histMassSEPM -> Integral(minMassBin1, maxMassBin1) + histMassSEPM->Integral(minMassBin2, maxMassBin2);
            double normMassMEPM = histMassMEPM -> Integral(minMassBin1, maxMassBin1) + histMassMEPM->Integral(minMassBin2, maxMassBin2);
            double normV2SEPM = histV2SEPM -> Integral(minV2Bin1, maxV2Bin1) + histV2SEPM->Integral(minV2Bin2, maxV2Bin2);
            double normV2MEPM = histV2MEPM -> Integral(minV2Bin1, maxV2Bin1) + histV2MEPM->Integral(minV2Bin2, maxV2Bin2);

            Printf("Mass normalization factor = %f\n", normMassSEPM / normMassMEPM);
            Printf("V2 normalization factor = %f (%f - %f)\n", normV2SEPM / normV2MEPM, normV2SEPM, normV2MEPM);

            canvasCentrScan -> cd(iCentr+1);
            gPad -> SetLogy(1);
            histMassSEPM -> GetXaxis() -> SetRangeUser(2, 5);
            histMassSEPM -> GetYaxis() -> SetRangeUser(1, 1e7);
            histMassMEPM -> GetXaxis() -> SetRangeUser(2, 5);
            histMassSEPM -> Draw("EP");
            histMassMEPM -> Scale(normMassSEPM / normMassMEPM);
            histMassMEPM -> Draw("EP SAME");

            TH1D *histBkgSubtrSEMP = (TH1D*) histMassSEPM -> Clone(Form("histMassBkgSubtrSEPM_%0.1f_%0.1f__%1.0f_%1.0f", 0., 30., minCentrBins[iCentr], maxCentrBins[iCentr]));
            histBkgSubtrSEMP -> Add(histMassMEPM, -1.);
            histBkgSubtrSEMP -> SetLineColor(kRed+1);
            histBkgSubtrSEMP -> SetMarkerStyle(24);
            histBkgSubtrSEMP -> SetMarkerSize(0.8);
            histBkgSubtrSEMP -> SetMarkerColor(kRed+1);
            histBkgSubtrSEMP -> Draw("EP SAME");

            canvasV2CentrScan -> cd(iCentr+1);
            histV2SEPM -> GetXaxis() -> SetRangeUser(1, 5);
            histV2SEPM -> Draw("EP");
            //histV2MEPM -> Scale(histV2SEPM -> GetBinContent(histV2SEPM -> GetXaxis() -> FindBin(2)) / histV2MEPM -> GetBinContent(histV2MEPM -> GetXaxis() -> FindBin(2)));
            histV2MEPM -> Scale(normV2SEPM / normV2MEPM);
            histV2MEPM -> Draw("H SAME");

            fOut -> cd();
            histMassSEPM -> Write();
            histMassMEPM -> Write();
            histBkgSubtrSEMP -> Write();
        }

        canvasMassPtScan -> SaveAs(Form("LHC23_golden/flow/plot_mass_pt_scan_%s.pdf", muonCut.c_str()));
        canvasCentrScan -> SaveAs(Form("LHC23_golden/flow/plot_mass_centrality_scan_%s.pdf", muonCut.c_str()));
        canvasV2CentrScan -> SaveAs(Form("LHC23_golden/flow/plot_v2_centrality_scan_%s.pdf", muonCut.c_str()));

        fOut -> ls();
        fOut -> Close();
    }
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

TH1D* projectHistogram(THnSparseF *histSparse, double minPtRange, double maxPtRange, double minCentrRange, double maxCentrRange) {
    double minPtBin = histSparse -> GetAxis(1) -> FindBin(minPtRange);
    double maxPtBin = histSparse -> GetAxis(1) -> FindBin(maxPtRange - 0.01);
    double minCentrBin = histSparse -> GetAxis(3) -> FindBin(minCentrRange);
    double maxCentrBin = histSparse -> GetAxis(3) -> FindBin(maxCentrRange - 0.01);

    Printf("minPtBin = %0.1f, maxPtBin = %0.1f, minCentrBin = %0.1f, maxCentrBin = %0.1f", minPtBin, maxPtBin, minCentrBin, maxCentrBin);
    histSparse -> GetAxis(1) -> SetRange(minPtBin, maxPtBin);
    histSparse -> GetAxis(3) -> SetRange(minCentrBin, maxCentrBin);

    TH1D *histProj = (TH1D*) histSparse -> Projection(0, Form("histProj_%0.1f_%0.1f__%1.0f_%1.0f", minPtRange, maxPtRange, minCentrRange, maxCentrRange));
    return histProj;
}

TProfile* projectProfile(THnSparseF *histSparse, double minPtRange, double maxPtRange, double minCentrRange, double maxCentrRange, int var) {
    double minPtBin = histSparse -> GetAxis(1) -> FindBin(minPtRange);
    double maxPtBin = histSparse -> GetAxis(1) -> FindBin(maxPtRange - 0.01);
    double minCentrBin = histSparse -> GetAxis(3) -> FindBin(minCentrRange);
    double maxCentrBin = histSparse -> GetAxis(3) -> FindBin(maxCentrRange - 0.01);

    Printf("minPtBin = %0.1f, maxPtBin = %0.1f, minCentrBin = %0.1f, maxCentrBin = %0.1f", minPtBin, maxPtBin, minCentrBin, maxCentrBin);
    histSparse -> GetAxis(1) -> SetRange(minPtBin, maxPtBin);
    histSparse -> GetAxis(3) -> SetRange(minCentrBin, maxCentrBin);

    TH2D *hist2DProj = (TH2D*) histSparse -> Projection(var, 0, "Projection");
    hist2DProj -> Rebin2D(2);
    TProfile *histProj = (TProfile*) hist2DProj -> ProfileX(Form("histProj_%i_%0.1f_%0.1f__%1.0f_%1.0f", var, minPtRange, maxPtRange, minCentrRange, maxCentrRange));
    return histProj;
}
