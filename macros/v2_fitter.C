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
#include "TLatex.h"
#include "TStyle.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TPaveText.h"
#include <vector>

Double_t FitFunctionBackgroundPol2(Double_t *x, Double_t *par);
Double_t FitFunctionBackgroundPol3(Double_t *x, Double_t *par);
Double_t alphaCB2VWG(Double_t *x, Double_t *par);
Double_t FitFunctionFlowS2CB2VWGPOL2(Double_t *x, Double_t *par);
void SetFitRejectRange(Double_t a, Double_t b);
Double_t BackgroundVWG(Double_t *x, Double_t *par);
Double_t CrystalBallExtended(Double_t *x,Double_t *par);
Double_t FitFunctionSignalCrystalBallExtended(Double_t *x,Double_t *par);
Double_t FitFunctionBackgroundVWG(Double_t *x, Double_t *par);
Double_t FitFunctionFlowS2CB2POL2EXPPOL2(Double_t *x, Double_t *par);
Double_t alphaCB2POL2EXP(Double_t*x, Double_t* par);
Double_t alphaCB2POL4EXP(Double_t*x, Double_t* par);
Double_t FitFunctionBackgroundPol2Exp(Double_t *x, Double_t *par);
Double_t FitFunctionBackgroundPol4Exp(Double_t *x, Double_t *par);
Double_t FitFunctionNA60New(Double_t *x,Double_t *par);
Double_t FitFunctionFlowS2NA60NEWPOL2EXPPOL2(Double_t *x, Double_t *par);
Double_t FitFunctionFlowS2NA60NEWVWGPOL2(Double_t *x, Double_t *par);
Double_t alphaNA60NEWVWG(Double_t*x, Double_t* par);
Double_t alphaNA60NEWPOL2EXP(Double_t*x, Double_t* par);
Double_t alphaNA60NEWPOL4EXP(Double_t*x, Double_t* par);
Double_t FitFunctionBackgroundPol2Cheb(Double_t *x, Double_t *par);
Double_t FitFunctionFlowS2CB2VWGPol3Cheb(Double_t *x, Double_t *par);
Double_t FitFunctionFlowS2NA60NEWPOL4EXPPOL2(Double_t *x, Double_t *par);
Double_t FitFunctionFlowS2CB2POL4EXPPOL2(Double_t *x, Double_t *par);
Double_t fitFunctionCB2VWG(Double_t *x, Double_t *par);
Double_t fitFunctionCB2Pol4Exp(Double_t *x, Double_t *par);
Double_t fitFunctionNA60NEWVWG(Double_t *x, Double_t *par);
Double_t fitFunctionNA60NEWPol4Exp(Double_t *x, Double_t *par);

TF1* v2Fit(int iPt,TProfile *v2prof, int type,Double_t parLimit[3], Double_t minFitRange, Double_t maxFitRange);
double* DoFlowFit(TH1D *, TH1D *, TProfile *, double , double , bool , string , TFile *);

Double_t xmin= 2.3;
Double_t xmax= 4.8;

const Double_t v2xmin= 2.3;
const Double_t v2xmax= 4.8;

vector<double> fitPars;
TFile *fOut;
TFile *fv2fit;
Double_t Jpsiv2 ;
Double_t Jpsiv2_err ;
Double_t v2_ch_ndf;

vector<double> v2par;
vector<double> v2par_er;
std::vector<TH1*> histSystematics;

std::string Fit_subStr1 = "CB2_VWG";
std::string Fit_subStr2 = "CB2_Pol4Exp";
std::string Fit_subStr3 = "NA60_Pol4Exp";
std::string Fit_subStr4 = "NA60_VWG";


const int nMassFitTrials = 1;
string massSigTrials[] = {"CB2", "CB2", "CB2", "CB2", "NA60", "NA60"};
string massBkgTrials[] = {"VWG_data", "VWG_MC", "Pol4Exp_data", "Pol4Exp_MC", "VWG_MC", "Pol4Exp_MC"};

string func_sig[] = {"CB2", "CB2", "CB2", "CB2", "NA60", "NA60"};
string func_bkg[] = {"VWG", "VWG", "Pol4Exp", "Pol4Exp", "VWG", "Pol4Exp"};

const int nSigFuncs = 1;
string sigFuncs[] = {"CB2", "NA60"};
const int nBkgFuncs = 1;
string bkgFuncs[] = {"VWG", "Pol4Exp"};
const int nTailSets = 2;
string tailSets[] = {"data", "MC"};

double minCentrBins[] = {10};
double maxCentrBins[] = {50};

const int nPtBins = 1;
//double minPtBins[] = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12};
//double maxPtBins[] = {1, 2, 3, 4, 5, 6, 8, 10, 12, 15};
double minPtBins[] = {2};
double maxPtBins[] = {3};

const int nFitRanges = 3;
double minFitRanges[] = {2.3, 2.4, 2.5};
double maxFitRanges[] = {4.7, 4.6, 4.5};

const int nbck_type = 1;
int bck_type[2] = {0,1};
string func_v2bck[] = {"Pol2", "Cheb"};

int nSysbins = nMassFitTrials * nFitRanges * nbck_type;

double resolution = 1.28004;

double norm_sig[] = {2.12944e+04,3.61983e+04,2.49511e+04,1.24006e+04,6.77359e+03,3.69119e+03,3.02467e+03,9.84673e+02,2.90510e+02,100,10};
double norm_bkg[] = {8.12173e+06,1.13264e+07,1.48518e+06,2.21428e+05,4.73500e+04,1.51403e+04,8.22036e+03,2.02200e+03,5.60000e+02,100,10};

TProfile *prof_v2EP[20];
TProfile *prof_MEv2EP[20];
TH1D *hT_mass[20];
TH1D *hT_mass2[20];
TCanvas *m_plot[20];
TF1 *fitv2[20][20][20];

double par[19];
string sig_plus_bkg[] = {"bkg","aa","bb","cc","sig_Jpsi","mean_Jpsi","width_Jpsi","a_Jpsi","b_Jpsi","c_Jpsi","d_Jpsi"};
string cb2sig_plus_pol4expbkg[] = {"bkg","aa","bb","cc","dd","ee","ff","sig_Jpsi","mean_Jpsi","width_Jpsi","a_Jpsi","b_Jpsi","c_Jpsi","d_Jpsi"};
string signa60_plus_bkg[] = {"bkg","aa","bb","cc","sig_Jpsi","mean_Jpsi","width_Jpsi","b_Jpsi","c_Jpsi","d_Jpsi","f_Jpsi","g_Jpsi","h_Jpsi","a_Jpsi","e_Jpsi"};
string signa60_plus_pol4expbkg[] = {"bkg","aa","bb","cc","dd","ee","ff","sig_Jpsi","mean_Jpsi","width_Jpsi","b_Jpsi","c_Jpsi","d_Jpsi","f_Jpsi","g_Jpsi","h_Jpsi","a_Jpsi","e_Jpsi"};

std::vector<TH1*> Uncertainities(std::vector<TH1*> histlist, vector<double> parameter, vector<double> parameter_er, double value[nPtBins][3]);

void v2_fitter() {
    for (int iPt = 0;iPt < nPtBins;iPt++) {
        double jpsiV2;
        double errJpsiV2;
        vector<double> jpsiV2s;
        vector<double> errJpsiV2s;

        double minPtBin = minPtBins[iPt];
        double maxPtBin = maxPtBins[iPt];

        TFile *fIn = new TFile("/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/data/pass3/LHC23_full/Histograms_Fullpass3matchedMchMid_centr_Add_Pt_bin_10_50.root");
        TFile *fOut = new TFile(Form("systematics/fitResults_Pt_%2.1f_%2.1f.root", minPtBin, maxPtBin), "RECREATE");

        for (int iSigFunc = 0;iSigFunc < nSigFuncs;iSigFunc++) {
            for (int iBkgFunc = 0;iBkgFunc < nSigFuncs;iBkgFunc++) {
                for (int iTailSet = 0;iTailSet < nTailSets;iTailSet++) {
                    string sigFunc = sigFuncs[iSigFunc];
                    string bkgFunc = bkgFuncs[iBkgFunc];
                    string tailSet = tailSets[iTailSet];
                    TFile *fInMassFit = new TFile(Form("/Users/lucamicheletti/GITHUB/dq_fitter/analysis/LHC23_pass3_full/centrality_10_50/pt_dependence/Pt_%2.1f_%2.1f/multi_trial_%s_%s_%s_tails.root", minPtBin, maxPtBin, sigFunc.c_str(), bkgFunc.c_str(), tailSet.c_str()));

                    for (int iFitRange = 0;iFitRange < nFitRanges;iFitRange++) {
                        double minFitRange = minFitRanges[iFitRange];
                        double maxFitRange = maxFitRanges[iFitRange];

                        TH1D *histMassFitPars = (TH1D*) fInMassFit -> Get(Form("fit_results_CB2_VWG__%2.1f_%2.1f_histMassSEPM_%1.0f_%1.0f__10_50", minFitRange, maxFitRange, minPtBin, maxPtBin));
                        TH1D *histMass = (TH1D*) fIn -> Get(Form("histMassSEPM_%1.0f_%1.0f__10_50", minPtBin, maxPtBin));
                        TProfile *profV2 = (TProfile*) fIn -> Get(Form("histV2SEPM_%1.0f_%1.0f__10_50", minPtBin, maxPtBin));

                        // Plot options
                        histMass -> GetXaxis() -> SetRangeUser(minFitRange, maxFitRange);
                        histMass -> GetXaxis() -> SetTitle("#it{M}_{#mu#mu} (GeV/#it{c^{2}})");
                        histMass -> SetMarkerStyle(20);
                        histMass -> SetMarkerSize(0.8);
                        histMass -> SetMarkerColor(kBlack);
                        histMass -> SetLineColor(kBlack);

                        profV2 -> SetTitle("");
                        profV2 -> GetXaxis() -> SetRangeUser(minFitRange, maxFitRange);
                        profV2 -> GetXaxis() -> SetTitle("#it{M}_{#mu#mu} (GeV/#it{c^{2}})");
                        profV2 -> SetMarkerStyle(20);
                        profV2 -> SetMarkerSize(0.8);
                        profV2 -> SetMarkerColor(kBlack);
                        profV2 -> SetLineColor(kBlack);

                        double *results = DoFlowFit(histMassFitPars, histMass, profV2, minFitRange, maxFitRange, kFALSE, Form("%s_%s_%s", sigFunc.c_str(), bkgFunc.c_str(), tailSet.c_str()), fOut);

                        jpsiV2s.push_back(results[0]);
                        errJpsiV2s.push_back(results[1]);
                        std::cout << results[0] << " +/- " << results[1] << std::endl;
                    }
                }
            }
        }

        //double value[nPtBins][3];
        //std::vector<TH1*> myHist = Uncertainities(histSystematics, jpsiV2s, errJpsiV2s, value);


        fOut -> Close();
    }

    /*
    // Get signal shape from previous fit
    vector<double> fitPars;
    for (int i = 0; i < 11; i++) {
        fitPars.push_back(histMassFitPars -> GetBinContent(histMassFitPars -> GetXaxis() -> FindBin(Form("%s",sig_plus_bkg[i].data()))));
    }

    TF1 *funcMassCB2VWG = new TF1("funcMassCB2VWG", fitFunctionCB2VWG, minFitRange, maxFitRange, 11);
    for (int i = 1; i < 4; i++) {
        funcMassCB2VWG->FixParameter(i, fitPars[i]);
    }
    for (int i = 5; i < 11; i++) {
        funcMassCB2VWG->FixParameter(i, fitPars[i]);
    }
    funcMassCB2VWG->SetParameter(0, 1e6); // Bkg normalization
    funcMassCB2VWG->SetParameter(4, 300000); // Jpsi normalization

    histMass ->Fit(funcMassCB2VWG, "RL0"); // Get the normalization from the fit


    TF1 *funcAlphaCB2VWG = new TF1("funcAlphaCB2VWG", alphaCB2VWG, minFitRange, maxFitRange, 11);
    for (int iPar = 0;iPar < 11;iPar++) {
        funcAlphaCB2VWG -> FixParameter(iPar, funcMassCB2VWG -> GetParameter(iPar));
        //std::cout << "Par " << iPar << " -> " << funcAlphaCB2VWG -> GetParameter(iPar) << std::endl;
    }

    TF1 *funcMassVWG = new TF1("funcMassVWG", FitFunctionBackgroundVWG, minFitRange, maxFitRange, 4);
    funcMassVWG -> SetLineColor(kGray+1);
    funcMassVWG -> SetLineStyle(kDashed);
    for (int iPar = 0;iPar < 4;iPar++) {
        funcMassVWG -> FixParameter(iPar, funcMassCB2VWG -> GetParameter(iPar));
    }

    TF1 *funcMassCB2 = new TF1("funcMassCB2", FitFunctionSignalCrystalBallExtended, minFitRange, maxFitRange, 7);
    funcMassCB2 -> SetLineColor(kAzure+2);
    funcMassCB2 -> SetLineStyle(kSolid);
    for (int iPar = 0;iPar < 7;iPar++) {
        funcMassCB2 -> FixParameter(iPar, funcMassCB2VWG -> GetParameter(iPar+4));
    }

    TF1 *funcFlowCB2VWG  = new TF1("funcFlowCB2VWG ", FitFunctionFlowS2CB2VWGPOL2, minFitRange, maxFitRange, 17);
    funcFlowCB2VWG -> SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
    funcFlowCB2VWG -> SetParName(11, "kPsiP");
    funcFlowCB2VWG -> SetParName(12, "v_{2} JPsi");
    funcFlowCB2VWG -> SetParName(13, "v_{2} BG0");
    funcFlowCB2VWG -> SetParName(14, "v_{2} BG1");
    funcFlowCB2VWG -> SetParName(15, "v_{2} BG2");
    funcFlowCB2VWG -> SetParName(16, "type");
    funcFlowCB2VWG -> SetLineColor(kBlue);

    for (int iPar = 0;iPar < 11;iPar++) {
        funcFlowCB2VWG -> FixParameter(iPar, funcMassCB2VWG -> GetParameter(iPar));
    }
    funcFlowCB2VWG -> FixParameter(11, 0);
    funcFlowCB2VWG -> SetParameter(12, 0.072);
    funcFlowCB2VWG -> SetParameter(13, 0.09);
    funcFlowCB2VWG -> SetParameter(14, -0.0252894);
    funcFlowCB2VWG -> SetParameter(15, 0.0022179);
    funcFlowCB2VWG -> FixParameter(16, 0);
    profV2 -> Fit (funcFlowCB2VWG, "RI"); // Get the normalization from the fit


    TCanvas *canvasAlpha = new TCanvas("canvasAlpha", "", 800, 600);
    TH2D *histGridAlpha = new TH2D("histGridAlpha", "", 100, minFitRange, maxFitRange, 100, 0, 0.3);
    histGridAlpha -> Draw();
    funcAlphaCB2VWG -> Draw("SAME");
    funcFlowCB2VWG -> Draw("SAME");

    TCanvas *canvasFit = new TCanvas("canvasFit", "", 600, 1200);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 1, 1);  // Pad superiore
    TPad *pad2 = new TPad("pad2", "pad1", 0, 0, 1, 0.5);  // Pad inferiore

    pad1 -> SetBottomMargin(0);
    pad2 -> SetTopMargin(0);
    pad2 -> SetBottomMargin(0.1);

    canvasFit -> cd();
    pad1 -> Draw();
    pad2 -> Draw();

    pad1 -> cd();
    gPad -> SetLogy(kTRUE);
    histMass -> Draw();
    funcMassVWG -> Draw("SAME");
    funcMassCB2 -> Draw("SAME");
    funcMassCB2VWG -> Draw("SAME");

    pad2 -> cd();
    gStyle -> SetOptStat(0);
    gStyle -> SetOptFit(1111);
    profV2 -> Draw();
    funcFlowCB2VWG -> Draw("SAME");
    */

}
////////////////////////////////////////////////////////////
double* DoFlowFit(TH1D *histMassFitPars, TH1D *histMass, TProfile *profV2, double minFitRange, double maxFitRange, bool debug, string fitFuncs, TFile *fOut) {
    // Get signal shape from previous fit
    vector<double> fitPars;


    for (int i = 0; i < 11; i++) {
        fitPars.push_back(histMassFitPars -> GetBinContent(histMassFitPars -> GetXaxis() -> FindBin(Form("%s",sig_plus_bkg[i].data()))));
    }

    // STARTING THE FIT TO THE INVARIANT MASS
    TF1 *funcMassCB2VWG = new TF1("funcMassCB2VWG", fitFunctionCB2VWG, minFitRange, maxFitRange, 11);
    for (int i = 1; i < 4; i++) {
        funcMassCB2VWG->FixParameter(i, fitPars[i]);
    }
    for (int i = 5; i < 11; i++) {
        funcMassCB2VWG->FixParameter(i, fitPars[i]);
    }
    funcMassCB2VWG->SetParameter(0, 1e6); // Bkg normalization
    funcMassCB2VWG->SetParameter(4, 300000); // Jpsi normalization

    histMass ->Fit(funcMassCB2VWG, "RL0"); // Get the normalization from the fit

    TF1 *funcAlphaCB2VWG = new TF1("funcAlphaCB2VWG", alphaCB2VWG, minFitRange, maxFitRange, 11);
    for (int iPar = 0;iPar < 11;iPar++) {
        funcAlphaCB2VWG -> FixParameter(iPar, funcMassCB2VWG -> GetParameter(iPar));
        //std::cout << "Par " << iPar << " -> " << funcAlphaCB2VWG -> GetParameter(iPar) << std::endl;
    }

    TF1 *funcMassVWG = new TF1("funcMassVWG", FitFunctionBackgroundVWG, minFitRange, maxFitRange, 4);
    funcMassVWG -> SetLineColor(kGray+1);
    funcMassVWG -> SetLineStyle(kDashed);
    for (int iPar = 0;iPar < 4;iPar++) {
        funcMassVWG -> FixParameter(iPar, funcMassCB2VWG -> GetParameter(iPar));
    }

    TF1 *funcMassCB2 = new TF1("funcMassCB2", FitFunctionSignalCrystalBallExtended, minFitRange, maxFitRange, 7);
    funcMassCB2 -> SetLineColor(kAzure+2);
    funcMassCB2 -> SetLineStyle(kSolid);
    for (int iPar = 0;iPar < 7;iPar++) {
        funcMassCB2 -> FixParameter(iPar, funcMassCB2VWG -> GetParameter(iPar+4));
    }
    
    TF1 *funcFlowCB2VWG  = new TF1("funcFlowCB2VWG ", FitFunctionFlowS2CB2VWGPOL2, minFitRange, maxFitRange, 17);
    funcFlowCB2VWG -> SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
    funcFlowCB2VWG -> SetParName(11, "kPsiP");
    funcFlowCB2VWG -> SetParName(12, "v_{2} JPsi");
    funcFlowCB2VWG -> SetParName(13, "v_{2} BG0");
    funcFlowCB2VWG -> SetParName(14, "v_{2} BG1");
    funcFlowCB2VWG -> SetParName(15, "v_{2} BG2");
    funcFlowCB2VWG -> SetParName(16, "type");
    funcFlowCB2VWG -> SetLineColor(kBlue);

    if (debug) {
        TCanvas *canvasAlpha = new TCanvas("canvasAlpha", "", 800, 600);
        TH2D *histGridAlpha = new TH2D("histGridAlpha", "", 100, minFitRange, maxFitRange, 100, 0, 0.3);
        histGridAlpha -> Draw();
        funcAlphaCB2VWG -> Draw("SAME");
        funcFlowCB2VWG -> Draw("SAME");
    }

    for (int iPar = 0;iPar < 11;iPar++) {
        funcFlowCB2VWG -> FixParameter(iPar, funcMassCB2VWG -> GetParameter(iPar));
    }
    funcFlowCB2VWG -> FixParameter(11, 0);
    funcFlowCB2VWG -> SetParameter(12, 0.072);
    funcFlowCB2VWG -> SetParameter(13, 0.09);
    funcFlowCB2VWG -> SetParameter(14, -0.0252894);
    funcFlowCB2VWG -> SetParameter(15, 0.0022179);
    funcFlowCB2VWG -> FixParameter(16, 0);
    profV2 -> Fit (funcFlowCB2VWG, "RI"); // Get the normalization from the fit

    double* results = new double[2];
    results[0] = funcFlowCB2VWG -> GetParameter(12);
    results[1] = funcFlowCB2VWG -> GetParError(12);

    TCanvas *canvasFit = new TCanvas("canvasFit", "", 600, 1200);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 1, 1);  // Pad superiore
    TPad *pad2 = new TPad("pad2", "pad1", 0, 0, 1, 0.5);  // Pad inferiore

    pad1 -> SetBottomMargin(0);
    pad2 -> SetTopMargin(0);
    pad2 -> SetBottomMargin(0.1);

    canvasFit -> cd();
    pad1 -> Draw();
    pad2 -> Draw();

    pad1 -> cd();
    gPad -> SetLogy(kTRUE);
    histMass -> Draw();
    funcMassVWG -> Draw("SAME");
    funcMassCB2 -> Draw("SAME");
    funcMassCB2VWG -> Draw("SAME");

    pad2 -> cd();
    gStyle -> SetOptStat(0);
    gStyle -> SetOptFit(1111);
    profV2 -> Draw();
    funcFlowCB2VWG -> Draw("SAME");

    fOut -> cd();
    canvasFit -> Write(Form("v2_fit_%2.1f_%2.1f_%s", minFitRange, maxFitRange, fitFuncs.c_str()));

    return results;
}
















////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
std::vector<TH1*> Uncertainities(std::vector<TH1*> histlist, vector<double> parameter, vector<double> parameter_er, double value[nPtBins][3]) {
    std::vector<TH1*> histograms;
    int ncombination = parameter.size() / nPtBins;
    int nbin = 0;
    for (int i = 0; i < parameter.size(); i++) {
        if ((i) % ncombination != 0) {
            continue;
        }
        if ((i) % ncombination == 0) {
            cout<<" HERE IS THE DIVISIONS "<<i<<endl;
        }

        double sum_v2 = 0; 
        double sum_v2er = 0;
        int ibin =0;
        for (int j = i; j < i+ncombination; j++) {
            histlist[nbin]->SetBinContent(ibin+1,parameter[j]);
	        cout<<" v2 parameter inside  "<<i<<" "<< j<<" "<<parameter[j]<<" +/- "<<parameter_er[j]<<endl;
	        //cout<<" v2 parameter inside2  "<<histlist[nbin]->GetBinContent(ibin+1)<<endl;
	        sum_v2 = sum_v2 + parameter[j];
	        sum_v2er = sum_v2er + parameter_er[j];
	        ibin++;
	    }

        histograms.push_back(histlist[nbin]);
        double mean_v2 = sum_v2/ncombination;
        double mean_v2er = sum_v2er/ncombination;
      
        Double_t sum2 =0.0;
        for (int j = i; j < i+ncombination; j++) {
            Double_t dev = (parameter[i] - mean_v2);
	        sum2 = sum2 + dev*dev;	     
	    }

        double_t sys = TMath::Sqrt(sum2/ncombination);
        value[nbin][0] = mean_v2;
        value[nbin][1] = mean_v2er;
        value[nbin][2] = sys;
      
        Double_t per_stat = (mean_v2er/mean_v2)*100;
        Double_t per_sys = (sys/mean_v2)*100;

        cout<<" Average ==== "<<mean_v2<< " +/- " <<mean_v2er<<" ("<<per_stat<<"\%)"<<" +/- "<<sys<<" ("<<per_sys<<"\%)"<<endl;
        //cout<<" Size "<<histlist.size()<<" Bin Content ==== "<<j<<" "<<histlist[0]->GetBinContent(j)<<endl;
        nbin++;
    }
    return histograms;
}


void SetFitRejectRange(Double_t a, Double_t b)
{
  /// Set a range the fit function(s) can ignore
  if ( a <= TMath::Limits<Double_t>::Max() && b <= TMath::Limits<Double_t>::Max() ) TF1::RejectPoint();
}
//------------------------------------------------------------------------------
Double_t FitFunctionBackgroundPol2(Double_t *x, Double_t *par)
{
  int type = par[3];
  // pol2 3 params

  //if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  //if (x[0] > 3.5 && x[0] < 3.7) TF1::RejectPoint();
  if(type==0) return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];  //return pol2 function

  //Return chebyshev Pol2 function
  else return FitFunctionBackgroundPol2Cheb(x,par); 
}

//------------------------------------------------------------------------------
Double_t FitFunctionBackgroundPol3(Double_t *x, Double_t *par)
{
  // pol2 3 params

  //if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  //if (x[0] > 3.5 && x[0] < 3.7) TF1::RejectPoint();
  return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0];
}

//____________________________________________________________________________
Double_t FitFunctionBackgroundPol4Exp(Double_t *x, Double_t *par)
{
  // pol4 x exp : 6 params
 //if(x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
 //return par[0]*(par[1]+par[2]*x[0]+par[3]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]+par[5]*x[0]*x[0]*x[0]*x[0])*TMath::Exp(par[6]/x[0]);
 //return (par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]*x[0])*TMath::Exp(par[5]/x[0]);
 //return par[0]*(par[1]+par[2]*x[0]+par[3]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]+par[5]*x[0]*x[0]*x[0]*x[0])*TMath::Exp(par[6]*x[0]);
  return par[0]*(par[1]+par[2]*x[0]+par[3]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]+par[5]*x[0]*x[0]*x[0]*x[0])*TMath::Exp(-par[6]*x[0]);
  
}

//____________________________________________________________________________
Double_t FitFunctionBackgroundPol2Cheb(Double_t *x, Double_t *par)
{
  if(x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  Double_t xmin = 2.2;
  Double_t xmax = 4.7;
  double xx = (2.0 * x[0] - xmin -xmax)/(xmax-xmin);
  const int order = 2;
  Double_t fT[order+1] = {0};
  if (order == 1) return par[0];
  if (order == 2) return par[0] + xx*par[1];
  // build the polynomials
  fT[0] = 1;
  fT[1] = xx;
  for (int i = 1; i< order; ++i) {
    fT[i+1] =  2 *xx * fT[i] - fT[i-1];
  }
  double sum = par[0]*fT[0];
  for (int i = 1; i<= order; ++i) {
    sum += par[i] * fT[i];
  }
  return sum;
}
//____________________________________________________________________________

Double_t alphaCB2VWG(Double_t*x, Double_t* par)
{
  return FitFunctionSignalCrystalBallExtended(x, &par[4])/(FitFunctionSignalCrystalBallExtended(x, &par[4]) + FitFunctionBackgroundVWG(x,par));
}

//____________________________________________________________________________
Double_t FitFunctionSignalCrystalBallExtended(Double_t *x,Double_t *par)
{
  // Extended Crystal Ball : 7 parameters
  //
  // par[0] = Normalization
  // par[1] = mean
  // par[2] = sigma
  // par[3] = alpha
  // par[4] = n
  // par[5] = alpha'
  // par[6] = n'

  Double_t t = (x[0]-par[1])/par[2];
  if (par[3] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);

  if (t >= -absAlpha && t < absAlpha2) // gaussian core
  {
    return par[0]*(exp(-0.5*t*t));
  }

  if (t < -absAlpha) //left tail
  {
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[4]/absAlpha - absAlpha;
    return par[0]*(a/TMath::Power(b - t, par[4]));
  }

  if (t >= absAlpha2) //right tail
  {

    Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
    Double_t d = par[6]/absAlpha2 - absAlpha2;
    return par[0]*(c/TMath::Power(d + t, par[6]));
  }

  return 0. ;
}

//-------------------------------------------------------------------------------
Double_t FitFunctionBackgroundVWG(Double_t *x, Double_t *par)
{
  // gaussian variable width : 4 params

  //if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  //if (x[0] > 3.5 && x[0] < 3.7) TF1::RejectPoint();
  Double_t sigma = par[2]+par[3]*((x[0]-par[1])/par[1]);
  return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));
}


//-------------------------------------------------------------------------------
Double_t FitFunctionFlowS2CB2VWGPOL2(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = 1.0341;

  Double_t par2[11] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[11], //kPsi'
    par[5]+(3.68609-3.096916),
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10]
  };

  return alphaCB2VWG(x,par)*par[12] + (1. - alphaCB2VWG(x,par))*FitFunctionBackgroundPol2(x,&par[13]);
}



//-------------------------------------------------------------------------------
Double_t FitFunctionFlowS2CB2VWGPol3Cheb(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = 1.0341;

  Double_t par2[11] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[11], //kPsi'
    par[5]+(3.68609-3.096916),
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10]
  };


  return alphaCB2VWG(x,par)*par[12] + (1. - alphaCB2VWG(x,par))*FitFunctionBackgroundPol2Cheb(x,&par[13]);
}


//-------------------------------------------------------------------------------
Double_t FitFunctionFlowS2CB2POL2EXPPOL2(Double_t *x, Double_t *par)
{
  Double_t SPsiPFactor = 1.0341;
  return alphaCB2POL2EXP(x,par)*par[12]  + (1. - alphaCB2POL2EXP(x,par))*FitFunctionBackgroundPol2(x,&par[13]);

}

//-------------------------------------------------------------------------------
Double_t FitFunctionFlowS2CB2POL4EXPPOL2(Double_t *x, Double_t *par)
{
  return alphaCB2POL4EXP(x,par)*par[15]  + (1. - alphaCB2POL4EXP(x,par))*FitFunctionBackgroundPol2(x,&par[16]);

}

//-------------------------------------------------------------------------------
Double_t alphaCB2POL2EXP(Double_t*x, Double_t* par)
{
  return FitFunctionSignalCrystalBallExtended(x, &par[4])/(FitFunctionSignalCrystalBallExtended(x, &par[4]) + FitFunctionBackgroundPol2Exp(x,par));
}

//-------------------------------------------------------------------------------
Double_t alphaCB2POL4EXP(Double_t*x, Double_t* par)
{
  return FitFunctionSignalCrystalBallExtended(x, &par[7])/(FitFunctionSignalCrystalBallExtended(x, &par[7]) + FitFunctionBackgroundPol4Exp(x,par));
}

//____________________________________________________________________________
Double_t FitFunctionBackgroundPol2Exp(Double_t *x, Double_t *par)
{
  // pol2 x exp : 4 params
  //if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  //if (x[0] > 3.5 && x[0] < 3.7) TF1::RejectPoint();
//  return par[0]*(par[1]+par[2]*x[0]+par[3]*x[0]*x[0])*TMath::Exp(par[4]/x[0]);
  return (par[0]+par[1]*x[0]+par[2]*x[0]*x[0])*TMath::Exp(par[3]*x[0]);
}



//____________________________________________________________________________
Double_t FitFunctionNA60New(Double_t *x,Double_t *par)
{
  // New Formulation of NA60 : 11 parameters
  //
  // par[0] = Normalization
  // par[1] = mean
  // par[2] = sigma
  // par[3] = p1Left
  // par[4] = p2Left
  // par[5] = p3Left
  // par[6] = p1Right
  // par[7] = p2Right
  // par[8] = p3Right
  // par[9] = alphaLeft
  // par[10] = alphaRight


  const Double_t t = (x[0]-par[1])/par[2];

  Double_t sigmaRatio(0.);
  if( t < par[9] ) sigmaRatio = ( 1.0 + TMath::Power( par[3]*(par[9]-t), par[4]-par[5]*TMath::Sqrt(par[9] - t) ) );
  else if( t >= par[9] && t < par[10] ) sigmaRatio = 1;
  else if( t >= par[10] ) sigmaRatio = ( 1.0 + TMath::Power( par[6]*(t-par[10]), par[7]-par[8]*TMath::Sqrt(t - par[10]) ) );

  return par[0]*TMath::Exp( -(1/2.)*TMath::Power(t/sigmaRatio,2.));

}

//Mass Fits
//------------------------------------------------------------------------------
Double_t fitFunctionCB2VWG(Double_t *x, Double_t *par)
{
  return FitFunctionBackgroundVWG(x, par) + FitFunctionSignalCrystalBallExtended(x, &par[4]);
}

//------------------------------------------------------------------------------
Double_t fitFunctionNA60NEWVWG(Double_t *x, Double_t *par)
{
  return FitFunctionBackgroundVWG(x, par) + FitFunctionNA60New(x, &par[4]);
}

//------------------------------------------------------------------------------
Double_t fitFunctionCB2Pol4Exp(Double_t *x, Double_t *par)
{
  return FitFunctionBackgroundPol4Exp(x, par) + FitFunctionSignalCrystalBallExtended(x, &par[7]);
}

//------------------------------------------------------------------------------
Double_t fitFunctionNA60NEWPol4Exp(Double_t *x, Double_t *par)
{
  return FitFunctionBackgroundPol4Exp(x, par) + FitFunctionNA60New(x, &par[7]);
}

//------------------------------------------------------------------------------
Double_t alphaNA60NEWVWG(Double_t*x, Double_t* par)
{
  return FitFunctionNA60New(x, &par[4])/(FitFunctionNA60New(x, &par[4]) + FitFunctionBackgroundVWG(x,par));
}

//------------------------------------------------------------------------------
Double_t alphaNA60NEWPOL2EXP(Double_t*x, Double_t* par)
{
  return FitFunctionNA60New(x, &par[4])/(FitFunctionNA60New(x, &par[4]) + FitFunctionBackgroundPol2Exp(x,par));
}

//------------------------------------------------------------------------------
Double_t alphaNA60NEWPOL4EXP(Double_t*x, Double_t* par)
{
  return FitFunctionNA60New(x, &par[7])/(FitFunctionNA60New(x, &par[7]) + FitFunctionBackgroundPol4Exp(x,par));
}

//____________________________________________________________________________
Double_t FitFunctionFlowS2NA60NEWPOL2EXPPOL2(Double_t *x, Double_t *par)
{
  Double_t SPsiPFactor = 1.0341;
  Double_t par2[15] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[15],
    par[5]+(3.68609-3.096916),
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10],
    par[11],
    par[12],
    par[13],
    par[14],
  };

  return alphaNA60NEWPOL2EXP(x,par)*par[16] + (1. - alphaNA60NEWPOL2EXP(x,par))*FitFunctionBackgroundPol2(x,&par[17]);

}

//____________________________________________________________________________
Double_t FitFunctionFlowS2NA60NEWPOL4EXPPOL2(Double_t *x, Double_t *par)
{
  return alphaNA60NEWPOL4EXP(x,par)*par[19] + (1. - alphaNA60NEWPOL4EXP(x,par))*FitFunctionBackgroundPol2(x,&par[20]);
}


//____________________________________________________________________________
Double_t FitFunctionFlowS2NA60NEWVWGPOL2(Double_t *x, Double_t *par)
{
  Double_t SPsiPFactor = 1.0341;

  Double_t par2[15] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[15],
    par[5]+(3.68609-3.096916),
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10],
    par[11],
    par[12],
    par[13],
    par[14],
  };

  return alphaNA60NEWVWG(x,par)*par[16] + (1. - alphaNA60NEWVWG(x,par))*FitFunctionBackgroundPol2(x,&par[17]);

}

