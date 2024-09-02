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

void LoadStyle();
void SetLegend(TLegend *);

Double_t FitFunctionBackgroundPol2(Double_t *x, Double_t *par);
Double_t FitFunctionBackgroundPol3(Double_t *x, Double_t *par);
Double_t FitFunctionSignalCrystalBallExtended(Double_t *x,Double_t *par);
Double_t alphaCB2VWG(Double_t *x, Double_t *par);
Double_t FitFunctionBackgroundVWG(Double_t *x, Double_t *par);
Double_t FitFunctionFlowS2CB2VWGPOL2(Double_t *x, Double_t *par);
Double_t FitFunctionMeanPtS2CB2VWGPOL3(Double_t *x, Double_t *par);
void SetFitRejectRange(Double_t a, Double_t b);

Double_t func2CB2VWG(Double_t *x, Double_t *par);
Double_t BackgroundVWG(Double_t *x, Double_t *par);
Double_t CrystalBallExtended(Double_t *x,Double_t *par);
void ProcessMinvFit(TH1D *fHisto, TFitResultPtr& fitResult, TF1* fitTotal, TF1* bckInit, const char* fitOption, Int_t iParKPsip, Int_t iLastParBkg);
void SetFitRejectRange(Double_t a, Double_t b);
void v2Fit(TProfile *v2prof, Double_t par[12],Double_t parLimit[3], TFitResultPtr r);
void MeanPtFit(TProfile *Ptprof, Double_t par[12],Double_t parLimit[3], TFitResultPtr r);

      Double_t xmin= 2.3;
      Double_t xmax= 4.8;

const Double_t v2xmin= 2.3;
const Double_t v2xmax= 4.8;

const Double_t Ptxmin= 2.3;
const Double_t Ptxmax= 4.8;

TFile *fOut;
Double_t Jpsiv2 ;
Double_t Jpsiv2_err ;
Double_t v2_ch_ndf;

void Mass_Fit_pTBin()
{
  LoadStyle();
  gStyle -> SetOptStat(0);
  float pT_bin[11] = {0.000000, 1.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000,12.000000,20.000000};
  int cent[2] = {1, 2}; // 10 - 20 %
 
  TCanvas *c0 = new TCanvas("c0", "c0",75,70,1309,696);
  TCanvas* c1 = new TCanvas("c1", "c1",75,70,1309,696);
  c0->Divide(8,2,0.0);
  c1->Divide(8,2,0.0);

  TH2D *th2 = new TH2D("th2","", 360, 1.5, 5, 6000, -0.02, 0.08);
  th2->SetStats(0);
  th2->GetXaxis()->SetLabelFont(43);
  th2->GetXaxis()->SetLabelSize(12); // labels will be 14 pixels
  th2->GetYaxis()->SetLabelFont(43);
  th2->GetYaxis()->SetLabelSize(12); // labels will be 14 pixels
  th2->GetYaxis()->SetTitle("Mass");
  th2->GetXaxis()->SetTitle("v_{2}"); 
  //th2->Draw("0");
  
  TH1D *hT_mass[10];
  TH1D *hT_mass_MEPM[10];
  TProfile *prof_v2EP[10];
  TProfile *prof_MEv2EP[10];
  TProfile *prof_v2SP[10];
  TProfile *prof_MEv2SP[10];
  TProfile *prof_meanPtSE[10];
  TProfile *prof_meanPtME[10];			  

  double par[12];
  //const int nPtBins = 5;
  double minCentrBins[] = {10};
  double maxCentrBins[] = {50};
  //double minPtBins[] = {0, 1, 2, 4, 5};
  //double maxPtBins[] = {1, 2, 4, 5, 10};

  const int nPtBins = 8;
  double minPtBins[] = {0.0,1.0,2,3,4,5,6,8};
  double maxPtBins[] = {1.0,2.0,3,4,5,6,8,10};

  //double minPtBins[] = {0};
  //double maxPtBins[] = {5};

  fOut = new TFile("Histogram_Cent10-50.root","RECREATE");
  TFile *f = new TFile("Histograms_Fullpass3matchedMchMid_3cent.root");
  if (f->IsOpen()) cout << "File opened successfully" << endl;

  for (int iPt = 0;iPt < nPtBins;iPt++) {
    
  hT_mass[iPt] = (TH1D*) f->Get(Form("histMassSEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
  hT_mass_MEPM[iPt] = (TH1D*) f->Get(Form("histMassMEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
  prof_v2EP[iPt] = (TProfile*) f->Get(Form("histV2SEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
  prof_MEv2EP[iPt] = (TProfile*) f->Get(Form("histV2MEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
  prof_v2SP[iPt] = (TProfile*) f->Get(Form("histU2Q2SEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
  prof_MEv2SP[iPt] = (TProfile*) f->Get(Form("histU2Q2MEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
  prof_meanPtSE[iPt] = (TProfile*) f->Get(Form("histPtSEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
  prof_meanPtME[iPt] = (TProfile*) f->Get(Form("histPtSEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));

  //hT_mass[iPt]->Rebin(2);
  //hT_mass_MEPM[iPt]->Rebin(2);

  fOut->cd();
  //hT_mass[iPt]->Write();
  

  hT_mass[iPt]->SetMarkerColor(kBlack);
  hT_mass[iPt]->SetMarkerStyle(20);
  hT_mass[iPt]->SetMarkerSize(0.6); 
  hT_mass_MEPM[iPt]->SetMarkerColor(kBlue);
  hT_mass_MEPM[iPt]->SetMarkerStyle(20);
  hT_mass_MEPM[iPt]->SetMarkerSize(0.6);

  prof_v2EP[iPt]->SetMarkerColor(kBlack);
  prof_v2EP[iPt]->SetMarkerStyle(20);
  prof_v2EP[iPt]->SetMarkerSize(0.6);
  prof_MEv2EP[iPt]->SetMarkerColor(kBlue);
  prof_MEv2EP[iPt]->SetMarkerStyle(20);
  prof_MEv2EP[iPt]->SetMarkerSize(0.6);
 
   
  //fOut->WriteObject(hT_mass[iPt],"histMassSEPM_2_3__10_50");
  //fOut->WriteObject(hT_mass_MEPM[iPt],"histMassMEPM_2_3__10_50");
  //fOut->WriteObject(prof_v2EP[iPt],"histV2SEPM_2_3__10_50");
  //fOut->WriteObject(prof_MEv2EP[iPt],"histV2MEPM_2_3__10_50");

  
  //fOut->WriteObject(hT_mass[iPt],Form("histMassSEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
  //fOut->WriteObject(prof_v2EP[iPt],Form("histV2SEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
  //fOut->WriteObject(hT_mass_MEPM[iPt],Form("histMassMEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
  //fOut->WriteObject(prof_MEv2EP[iPt],Form("histV2MEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
 
  //if(iPt==0) {xmin = 2.0 ; xmax = 5.0;}
  //else {xmin = 2.5 ; xmax = 4.8;}
  //hT_mass[i]->Draw();

  /*
  TF1* bck = new TF1("bck",BackgroundVWG,2.2,4.8,4);
  bck->SetLineColor(kBlue);
  bck->SetLineStyle(2);
  bck->SetParameter(0, 5000000.);
  bck->SetParameter(1, 1.9);
  bck->SetParameter(2, 0.5);
  bck->SetParLimits(2, 0., 100.);
  bck->SetParameter(3, 0.3);
  bck->SetParLimits(3, 0., 100.);
  SetFitRejectRange(2.9,3.2);
  hT_mass[iPt]->Fit(bck,"REM","",2.0,5.0);
  */
  
  hT_mass[iPt]->GetXaxis()->SetRangeUser(xmin,xmax);
  TF1* fitTotal = new TF1("fitTotal",func2CB2VWG,xmin,xmax,12);
  fitTotal->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kPsi","mPsi","sPsi","alPsi","nlPsi","auPsi","nuPsi");
  fitTotal->SetParName(11, "kPsi'");

  fitTotal->SetParameter(0, 5000000.);
  fitTotal->SetParameter(1, 1.9);
  fitTotal->SetParameter(2, 0.5);
  fitTotal->SetParLimits(2, 0., 100.);
  fitTotal->SetParameter(3, 0.3);
  fitTotal->SetParLimits(3, 0., 100.);
  fitTotal->SetParameter(4, 100.);
  fitTotal->SetParameter(5, 3.1);
  fitTotal->SetParLimits(5, 3.085, 3.2);
  fitTotal->SetParameter(6, 0.072);
  fitTotal->SetParLimits(6, 0.06, 0.15);
  
  //r = FitJpsi2CB2VWG(*hminv,0.93,5.59,2.32,3.39);

  
    fitTotal->FixParameter(7,1.0169);
    //fitTotal->SetParLimits(7,0.1,10.0);
 
    fitTotal->FixParameter(8,3.7064);
    //fitTotal->SetParLimits(8,0.0,10.0);
  
    fitTotal->FixParameter(9, 2.4183);
    //fitTotal->SetParLimits(9,0.1,10.0);
  
    fitTotal->FixParameter(10,7.9655);
    //fitTotal->SetParLimits(10,0.0,10.0);
  
    fitTotal->FixParameter(11, 1.);

  const char* fitOption = "QSER"; //+";
  
  //TFitResultPtr fitResult = hfit->Fit(fitTotal,fitOption,"");

  c0->cd(iPt+1);
  hT_mass[iPt]->Draw("E");
  TFitResultPtr r = hT_mass[iPt]->Fit(fitTotal,"RESI");
  TMatrixDSym cov = r->GetCovarianceMatrix();
  fitTotal->GetParameters(&par[0]);
  fOut->WriteObject(fitTotal,"Mass_Signal_plus_bkg");
  
  TF1 *fitFctCB2 = new TF1("fitFctCB2", CrystalBallExtended,xmin,xmax, 7);
  fitFctCB2->SetLineColor(kRed);
  fitFctCB2->SetLineStyle(2);
  fitFctCB2->SetParNames("kPsi","mPsi","sPsi","alPsi","nlPsi","auPsi","nuPsi");
  fitFctCB2->FixParameter(0, fitTotal->GetParameter(4)); 
  fitFctCB2->FixParameter(1, fitTotal->GetParameter(5)); 
  fitFctCB2->FixParameter(2, fitTotal->GetParameter(6));
  fitFctCB2->FixParameter(3, fitTotal->GetParameter(7));   
  fitFctCB2->FixParameter(4, fitTotal->GetParameter(8));   
  fitFctCB2->FixParameter(5, fitTotal->GetParameter(9));   
  fitFctCB2->FixParameter(6, fitTotal->GetParameter(10));
  //fitFctCB2->Draw("same");

  TF1* bck2 = new TF1("bck2",BackgroundVWG,xmin,xmax,4);
  for ( Int_t i = 0; i < 4; ++i )
  {
    bck2->SetParameter(i,fitTotal->GetParameter(i));
  }
  bck2->SetLineColor(kGreen);
  bck2->SetLineStyle(2);
  //bck2->Draw("same");
  //bck->Draw("same");
  
  fOut->WriteObject(bck2,"Mass_bkg_noscale");

  double SE =hT_mass[iPt]->Integral(hT_mass[iPt]->GetXaxis()->FindBin(1.0),hT_mass[iPt]->GetXaxis()->FindBin(2.5)) + hT_mass[iPt]->Integral(hT_mass[iPt]->GetXaxis()->FindBin(3.72),hT_mass[iPt]->GetXaxis()->FindBin(5.0));

  double ME = hT_mass_MEPM[iPt]->Integral(hT_mass_MEPM[iPt]->GetXaxis()->FindBin(1.0),hT_mass_MEPM[iPt]->GetXaxis()->FindBin(2.5)) + hT_mass_MEPM[iPt]->Integral(hT_mass_MEPM[iPt]->GetXaxis()->FindBin(3.72),hT_mass_MEPM[iPt]->GetXaxis()->FindBin(5.0));

 cout<<" Normalization factor "<<ME/SE <<endl;


  hT_mass_MEPM[iPt]->Scale(SE/ME);
  //hT_mass_MEPM[iPt]->Draw("phsame");

  TLegend *legend1_1 = new TLegend(0.8125205,0.5806881,0.8948445,0.7828001,NULL,"brNDC");
  legend1_1->SetTextFont(42);
  legend1_1->SetTextSize(0.06);
  legend1_1->SetLineColor(0);
  legend1_1->SetFillColor(0);
  //legend1_1->AddEntry(prof_v2EP[iPt],"Same Event","Ep");
  //legend1_1->AddEntry(hT_mass_MEPM[iPt],"Mixed Event","Ep");
  if(iPt == 0) legend1_1->AddEntry(fitTotal,"Fit","l");
  legend1_1->Draw();
  
  //==============================================
  ///////////////////////////////////////////////////////////////
  std::cout << "FitResult = " << static_cast<int>(r) << std::endl;
  std::cout << "CovMatrixStatus = " << r->CovMatrixStatus() << std::endl;

  
  if ( ( static_cast<int>(r) && static_cast<int>(r)!=4000 ) ||  static_cast<int>(r->CovMatrixStatus())!=3 ) ProcessMinvFit(hT_mass[iPt],r,fitTotal,bck2,"RESLI",12,4);
  //Set("FitResult",static_cast<int>(r),0);
  //Set("CovMatrixStatus",static_cast<int>(r->CovMatrixStatus()),0);
  //printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n");
  //printf(" Fit Status : %d <-> Cov. Mat. : %d ",static_cast<int>(r),r->CovMatrixStatus());
  //printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n\n");


   //====================J/Psi Signal Estimation==========================================================
 
  // Set("NofJPsi",njpsi,nerr);
  
  double m =fitTotal->GetParameter(5);
  double m_err =fitTotal->GetParError(5);
  double s =fitTotal->GetParameter(6);
  double s_err =fitTotal->GetParError(6);
  cout<<" Mass "<<m<<endl;
  double njpsi3s = fitTotal->Integral(m-3*s,m+3*s)/hT_mass[iPt]->GetBinWidth(1);
  double nerr3s = fitTotal->IntegralError(m-3*s,m+3*s,fitTotal->GetParameters(),cov.GetMatrixArray())/hT_mass[iPt]->GetBinWidth(1);
  cout<< "No of J/Psi == " <<  njpsi3s <<"+/-" << nerr3s<<endl;
  cout<< "Sigma == " << s <<"+/-" <<s_err<<endl;
  Double_t chi=fitTotal->GetChisquare();
  Double_t ndf=fitTotal->GetNDF();
  Double_t ch_ndf=chi/ndf;
  
  //===============================================================================
  //Computation of bin significance and signal over background
  
  double nbck3s = bck2->Integral(m-3.*s,m+3.*s)/hT_mass[iPt]->GetBinWidth(1);
  //double nbck3sErr = exp2->IntegralError(m-3*s,m+3*s,r2->GetParams(),r2->GetCovarianceMatrix().GetMatrixArray() )/hT->GetBinWidth(1);
  double nbck3sErr = bck2->IntegralError(m-3*s,m+3*s,bck2->GetParameters(),cov.GetMatrixArray())/hT_mass[iPt]->GetBinWidth(1);
  double sOverB3s = njpsi3s / nbck3s;
  double sOverB3sErr = sOverB3s*TMath::Sqrt(TMath::Power(nerr3s/njpsi3s,2.) + TMath::Power(nbck3sErr/nbck3s,2.));

  cout<<" S/B ==="<< sOverB3s << "+/-" << sOverB3sErr <<endl;
  //Set("SignalOverBkg3s",sOverB3s,sOverB3sErr);
  
  double sig = njpsi3s/TMath::Sqrt(njpsi3s + nbck3s);
  double sigErr = TMath::Sqrt( TMath::Power((1. - (1./2.)*njpsi3s/(njpsi3s + nbck3s) )*nerr3s/TMath::Sqrt(njpsi3s + nbck3s),2.) +
                               TMath::Power(njpsi3s*nbck3sErr/(2.*TMath::Power(njpsi3s + nbck3s,3./2.)),2.) );


  cout<<" s/sqrt(s+b) " <<sig << " +/- "<<sigErr<<endl;  
  
  //==========================================================================================================
  

  cout<<"Double_t Fit_Par[12] = {";
  for(int i=0;i<12;i++) {
    cout<<fitTotal->GetParameter(i);
    if(i!=11) cout<<", ";
  }
  cout<<" }; "<<endl;

  TLatex *latexTitle = new TLatex();
  if(iPt == 0)
    {
      latexTitle -> SetTextSize(0.06);
      latexTitle -> SetNDC();
      latexTitle -> SetTextFont(42);
      //latexTitle -> DrawLatex(0.22, 0.87, "ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
      //latexTitle -> DrawLatex(0.22, 0.77, "J/#psi #rightarrow #mu^{+}#mu^{-}, 2.5 < y < 4");
      latexTitle -> DrawLatex(0.22, 0.87, Form("%1.0f#minus%1.0f %%",minCentrBins[0],maxCentrBins[0]));
      latexTitle -> DrawLatex(0.22, 0.77, Form("%1.0f < #it{p}_{T} < %1.0f GeV/#it{c}",minPtBins[iPt], maxPtBins[iPt]));
      //latexTitle -> DrawLatex(0.22, 0.67, Form("N_{J/#psi}(3#sigma) = %3.0f +/- %3.0f", njpsi3s,nerr3s));
      //latexTitle -> DrawLatex(0.25, 0.57, Form("#sigma_{J/#psi}(3#sigma) = %f +/- %f GeV",s,s_err));
      //latexTitle -> DrawLatex(0.15, 0.57, Form("S/B = %f +/- %f GeV",sOverB3s,sOverB3sErr));
      //latexTitle -> DrawLatex(0.5, 0.57, Form("S/(S+B) = %f +/- %f GeV",sig,sigErr));
      //latexTitle -> DrawLatex(0.15, 0.27, Form("#chi^{2}/ndf =%3.0f",ch_ndf));
    }
  else if(iPt == 4)
    {
      latexTitle -> SetTextSize(0.06);
      latexTitle -> SetNDC();
      latexTitle -> SetTextFont(42);
      latexTitle -> DrawLatex(0.06, 0.47, Form("%1.0f < #it{p}_{T} < %1.0f GeV/#it{c}",minPtBins[iPt], maxPtBins[iPt]));
      latexTitle -> DrawLatex(0.06, 0.37, Form("#chi^{2}/ndf =%3.0f",ch_ndf));
      //latexTitle -> DrawLatex(0.06, 0.27, Form("N_{J/#psi}(3#sigma) = %3.0f +/- %3.0f", njpsi3s, nerr3s));
      //latexTitle -> DrawLatex(0.06, 0.17, Form("#sigma_{J/#psi} = %f +/- %f GeV",s,s_err));
      //latexTitle -> DrawLatex(0.06, 0.17, Form("S/B = %f +/- %f GeV",sOverB3s,sOverB3sErr));
    }
  
  else
    {
      latexTitle -> SetTextSize(0.06);
      latexTitle -> SetNDC();
      latexTitle -> SetTextFont(42);
      latexTitle -> DrawLatex(0.22, 0.87, Form("%1.0f < #it{p}_{T} < %1.0f GeV/#it{c}",minPtBins[iPt], maxPtBins[iPt]));
      latexTitle -> DrawLatex(0.22, 0.77, Form("#chi^{2}/ndf =%3.0f",ch_ndf));
      //latexTitle -> DrawLatex(0.22, 0.67, Form("N_{J/#psi}(3#sigma) = %3.0f +/- %3.0f", njpsi3s, nerr3s));
      //latexTitle -> DrawLatex(0.22, 0.57, Form("#sigma_{J/#psi} = %f +/- %f GeV",s,s_err));
      //latexTitle -> DrawLatex(0.22, 0.57, Form("S/B = %f +/- %f GeV",sOverB3s,sOverB3sErr));
    }

  
  //v2 Fitting Limits ======
  Double_t parLimit0[3];
  parLimit0[0] = 0.005; 
  parLimit0[1] = 0.004;
  parLimit0[2] = 0.008;
  
  Double_t parLimit1[3];
  parLimit1[0] = 0.045; 
  parLimit1[1] = 0.042;
  parLimit1[2] = 0.048;

  Double_t parLimit2[3];
  parLimit2[0] = 0.084; 
  parLimit2[1] = 0.082;
  parLimit2[2] = 0.09;

  Double_t parLimit3[3];
  parLimit3[0] = 0.094; 
  parLimit3[1] = 0.082;
  parLimit3[2] = 0.099;

  Double_t parLimit4[3];
  parLimit4[0] = 0.13; 
  parLimit4[1] = 0.11;
  parLimit4[2] = 0.16;

   Double_t parLimit5[3];
  parLimit5[0] = 0.09; 
  parLimit5[1] = 0.05;
  parLimit5[2] = 0.13;

  Double_t parLimit6[3];
  parLimit6[0] = 0.06; 
  parLimit6[1] = 0.04;
  parLimit6[2] = 0.08;
  
  
  c0->cd(iPt+9);
  //th2->Draw("0");
  prof_v2EP[iPt]->Rebin(2);
  prof_v2EP[iPt]->Draw("p");
  prof_MEv2EP[iPt]->Draw("psame");
  prof_MEv2EP[iPt]->Rebin(2);
  
 //if(iPt ==1)
 //  {
 //   prof_v2EP[iPt]->SetMarkerColor(kRed);
 //    prof_v2EP[iPt]->Draw("p");
 //    prof_v2SP[iPt]->Draw("psame");
 //  }
  
  if(iPt ==0) v2Fit(prof_v2EP[iPt], par,parLimit0,r);
  else if(iPt ==1) v2Fit(prof_v2EP[iPt], par,parLimit1,r);
  else if(iPt ==2) v2Fit(prof_v2EP[iPt], par,parLimit2,r);
  else if(iPt ==3) v2Fit(prof_v2EP[iPt], par,parLimit3,r);
  else if(iPt ==4) v2Fit(prof_v2EP[iPt], par,parLimit4,r);
  else if(iPt ==5) v2Fit(prof_v2EP[iPt], par,parLimit5,r);
  else v2Fit(prof_v2EP[iPt], par,parLimit6,r);

 
  prof_MEv2EP[iPt]->Rebin(1);
  //prof_MEv2EP[iPt]->Scale(1.53);
  prof_MEv2EP[iPt]->Scale(2.7);
  // prof_MEv2EP[iPt]->Draw("psame");

  
  if(iPt ==0)
    {
  TLegend *legend2_1 = new TLegend(0.6125205,0.5806881,0.8948445,0.7828001,NULL,"brNDC");
  legend2_1->SetTextFont(42);
  legend2_1->SetTextSize(0.06);
  legend2_1->SetLineColor(0);
  legend2_1->SetFillColor(0);
  legend2_1->AddEntry(prof_v2EP[iPt],"Same event dimuon","p");
  legend2_1->AddEntry(prof_MEv2EP[iPt],"Mixed event dimuon","p");
  legend2_1->Draw();
    }
  
  cout<<" Pt Bins "<<minPtBins[iPt]<< " "<<maxPtBins[iPt]<<endl;
  cout<<" J/Psi v2  "<< Jpsiv2<<endl;
  cout<<" J/Psi v2 Error "<< Jpsiv2_err<<endl;

  TLatex *latexTitle2 = new TLatex();
   if(iPt == 4)
     {
   latexTitle2 -> SetTextSize(0.06);
   latexTitle2 -> SetNDC();
   latexTitle2 -> SetTextFont(42);
   latexTitle2 -> DrawLatex(0.22, 0.27, Form("#chi^{2}/ndf =%3.0f",v2_ch_ndf));
   latexTitle2 -> DrawLatex(0.22, 0.17, Form("v_{2} = %f +/- %f ",Jpsiv2,Jpsiv2_err));
     }

   else {
   latexTitle2 -> SetTextSize(0.06);
   latexTitle2 -> SetNDC();
   latexTitle2 -> SetTextFont(42);
   latexTitle2 -> DrawLatex(0.22, 0.27, Form("#chi^{2}/ndf =%3.0f",v2_ch_ndf));
   latexTitle2 -> DrawLatex(0.22, 0.17, Form("v_{2} = %f +/- %f ",Jpsiv2,Jpsiv2_err));
     }
	
	
  //c0->cd(iPt+11);
  //prof_meanPtSE[iPt]->Draw("p"); 
  //MeanPtFit(prof_meanPtSE[iPt], par,parLimit2, r);

   
    }


  
}



void v2Fit(TProfile *v2prof, Double_t par[12], Double_t parLimit[3], TFitResultPtr r)
{
  //v2prof->Scale(1/0.652);
  v2prof->Scale(1.280);
  v2prof->GetXaxis()->SetRangeUser(xmin,xmax);
  v2prof->Draw();
  //=============================== v2 Fitting =================
  TF1* bck = new TF1("bck",FitFunctionBackgroundPol2,xmin,xmax,3);
  bck->SetLineColor(kBlue);
  bck->SetLineStyle(9);
  bck->SetParameters(0.0191306,-0.00836616,0.00102991);
  SetFitRejectRange(3.0,4.0);
  v2prof->Fit(bck,"SRIE","",xmin,xmax);
  //bck->Draw("same");
  //SetFitRejectRange();

  TF1* fitv2 = new TF1("fitv2",FitFunctionFlowS2CB2VWGPOL2,v2xmin,v2xmax,16);
  fitv2->SetLineColor(kRed);
  fitv2->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
  fitv2->SetParName(11,"kPsiP");
  fitv2->SetParName(12,"v_{2} JPsi");
  fitv2->SetParName(13,"v_{2} BG0");
  fitv2->SetParName(14,"v_{2} BG1");
  fitv2->SetParName(15,"v_{2} BG2");

  //Double_t Fit_Par[12] = {3.97736e+07, 0.000772307, 0.44601, 0.000106707, 14438.4, 3.07397, 0.0770682, 9.75573, 6.75144, 9.81753, 7.71034, -113.444 };
  //Double_t Fit_Par[12] = {124166, 2.32721, 0.818141, 0.341954, 43926.9, 3.08071, 0.0847431, 9.75573, 6.75144, 9.81753, 7.71034, -113.444};
  //if ( ( static_cast<int>(r) && static_cast<int>(r)==4000 ) ||  static_cast<int>(r->CovMatrixStatus())==3 )
  {
  for( Int_t i = 0; i < 12; ++i ) fitv2->FixParameter(i,par[i]);
  }

  /*
  else
    {
      //if i == 0 ========
  Double_t Fit_Par[12] = {3.97736e+07, 0.000772307, 0.44601, 0.000106707, 14438.4, 3.07397, 0.0770682, 9.75573, 6.75144, 9.81753, 7.71034, -113.444 };
  for( Int_t i = 0; i < 12; ++i ) fitv2->FixParameter(i,Fit_Par[i]);
    }
  */
  
  fitv2->SetParameter(12, parLimit[0]);
  fitv2->SetParLimits(12, parLimit[1],parLimit[2]);
  
  
 for ( Int_t i = 0; i < 3; ++i )
  {
    fitv2->SetParameter(i + 13, bck->GetParameter(i));
  }

 TFitResultPtr r2 = v2prof->Fit(fitv2,"RESI");

 //fOut->WriteObject(fitv2,"v2_Signal_plus_bkg");
 
  if((fitv2->GetParError(12) - fitv2->GetParameter(12)) > 0 )
   {
     for( Int_t i = 0; i < 12; ++i ) fitv2->FixParameter(i,par[i]);
     //fitv2->SetParameter(12, fitv2->GetParameter(12)*0.9);
     //fitv2->SetParLimits(12,fitv2->GetParameter(12)*0.88,fitv2->GetParameter(12)*0.92);
     TFitResultPtr r3 = v2prof->Fit(fitv2,"RESI");
   }
 
  TF1* bck2 = new TF1("bck2",FitFunctionBackgroundPol2,v2xmin,v2xmax,3);
  for ( Int_t i = 0; i < 3; ++i )
  {
    bck2->SetParameter(i, fitv2->GetParameter(i + 13));
  }
  bck2->SetLineColor(kBlue);
  bck2->SetLineStyle(2);
  bck2->Draw("same");
  //fOut->WriteObject(bck2,"v2_bkg");

  Jpsiv2 =  fitv2->GetParameter(12);
  Jpsiv2_err =  fitv2->GetParError(12);
  v2_ch_ndf = fitv2->GetChisquare()/fitv2->GetNDF();
  //const char* fitOption = "SER";
}



void MeanPtFit(TProfile *Ptprof, Double_t par[12],Double_t parLimit[3], TFitResultPtr r)
{
  Ptprof->Draw("p");
  TF1* bck = new TF1("bck",FitFunctionBackgroundPol3,Ptxmin,Ptxmax,4);
  bck->SetLineColor(kBlue);
  bck->SetLineStyle(9);
  bck->SetParameters(6.71209,-0.879906,0.340363,-0.0333421);
  //bck->SetParLimits(0, 0.,5.0);
  //SetFitRejectRange(2.6,4.0);
  Ptprof->Fit(bck,"SRL","",Ptxmin,Ptxmax);
  bck->Draw();
  //SetFitRejectRange();

  TF1* fitMeanpt = new TF1("fitMeanpt",FitFunctionMeanPtS2CB2VWGPOL3,Ptxmin,Ptxmax,17);
  fitMeanpt->SetLineColor(kRed);
  fitMeanpt->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
  fitMeanpt->SetParName(11,"kPsiP");
  fitMeanpt->SetParName(12,"<pt>JPsi");
  fitMeanpt->SetParName(13,"<pt>BG0");
  fitMeanpt->SetParName(14,"<pt>BG1");
  fitMeanpt->SetParName(15,"<pt>BG2");
  fitMeanpt->SetParName(16,"<pt>BG3");

  Double_t Fit_Par[12] = {3.97736e+07, 0.000772307, 0.44601, 0.000106707, 14438.4, 3.07397, 0.0770682, 9.75573, 6.75144, 9.81753, 7.71034, -113.444 };
  if ( ( static_cast<int>(r) && static_cast<int>(r)==4000 ) ||  static_cast<int>(r->CovMatrixStatus())==3 )
    {
      for( Int_t i = 0; i < 12; ++i ) fitMeanpt->FixParameter(i,par[i]);
    }

    else
    {
      //if i == 0 ========
      Double_t Fit_Par[12] = {3.97736e+07, 0.000772307, 0.44601, 0.000106707, 14438.4, 3.07397, 0.0770682, 9.75573, 6.75144, 9.81753, 7.71034, -113.444 };
      for( Int_t i = 0; i < 12; ++i ) fitMeanpt->FixParameter(i,Fit_Par[i]);
    }

  fitMeanpt->SetParameter(12, 0.7);
  fitMeanpt->SetParLimits(12, 0.01,10.0);

 for ( Int_t i = 0; i < 4; ++i )
  {
    fitMeanpt->SetParameter(i + 13, bck->GetParameter(i));
  }


 TFitResultPtr r2 = Ptprof->Fit(fitMeanpt,"SREI");

  
  TF1* bck2 = new TF1("bck2",FitFunctionBackgroundPol2,Ptxmin,Ptxmax,4);
  for ( Int_t i = 0; i < 4; ++i )
  {
    bck2->SetParameter(i, bck->GetParameter(i));
  }
  bck2->SetLineColor(kBlue);
  bck2->SetLineStyle(2);
  bck2->Draw("same");
}



Double_t BackgroundVWG(Double_t *x, Double_t *par)
  {
    // gaussian variable width
    Double_t sigma = par[2]+par[3]*((x[0]-par[1])/par[1]);
    return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));
    
  }
  

 Double_t CrystalBallExtended(Double_t *x,Double_t *par)
  {
    //par[0] = Normalization
    //par[1] = mean
    //par[2] = sigma
    //par[3] = alpha
    //par[4] = n
    //par[5] = alpha'
    //par[6] = n'
    
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
    
    return 0.;
  }
  
  Double_t func2CB2VWG(Double_t *x, Double_t *par)
  {
    /// 2 extended crystal balls + variable width gaussian
    /// width of the second CB related to the first (free) one.
    
    Double_t par2[7] = {
      par[11],
      3.68609+(par[5]-3.096916)/3.096916*3.68609,
      par[6]/3.096916*3.68609,
      par[7],
      par[8],
      par[9],
      par[10]
    };
    //return BackgroundVWG(x, par) + CrystalBallExtended(x, &par[4]) + CrystalBallExtended(x, par2); //Jpsi and psi Fit
    return BackgroundVWG(x, par) + CrystalBallExtended(x, &par[4]); //Only J/psi
  }






void ProcessMinvFit(TH1D *fHisto, TFitResultPtr& fitResult, TF1* fitTotal, TF1* bckInit, const char* fitOption, Int_t iParKPsip, Int_t iLastParBkg)
{
  // If a Minv fit fails this algorithm changes some initial parameters to get the fit converged

  Int_t bin(0);

  if ( (static_cast<int>(fitResult) && (static_cast<int>(fitResult)!=4000||static_cast<int>(fitResult)!=0)  ) || static_cast<int>(fitResult->CovMatrixStatus())!=3 /*|| static_cast<int>(fitResult->CovMatrixStatus())!=2*/)
  {
    if ( (0.5*fitTotal->GetParameter(iParKPsip) <= fitTotal->GetParError(iParKPsip))) { //kPsi'
      std::cout << "-----------------------------------------------" << std::endl;
      std::cout << "------- Setting Psi'norm= Psi' norm*0.8 -------" << std::endl;
      std::cout << "-----------------------------------------------" << std::endl;
      bin = fHisto->FindBin(3.68);
      // fitTotal->SetParLimits(iParKPsip, 0.,fHisto->GetBinContent(bin)*1.5); // we further restrict the range of psi' norm
      fitTotal->SetParameter(iParKPsip, fHisto->GetBinContent(bin)*0.8);
    }

    if ( 0.5*fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ) {
      std::cout << "-----------------------------------------------" << std::endl;
      std::cout << "-------       Setting p0=MAX/2          -------" << std::endl;
      std::cout << "-----------------------------------------------" << std::endl;
    }

    std::cout << "================================\\" << std::endl;
    std::cout << "======== Refitting again =======\\" << std::endl;
    std::cout << "================================\\" << std::endl;
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

    //Check if there are poles somewhere/
    //if(iLastParBkg == 6 )CheckRoots(fitResult,fitTotal,3,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6),fitOption);
    //if(iLastParBkg == 5 )CheckRoots(fitResult,fitTotal,2,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),0.,fitOption);
   }

  if ( (static_cast<int>(fitResult) && (static_cast<int>(fitResult)!=4000||static_cast<int>(fitResult)!=0)) || static_cast<int>(fitResult->CovMatrixStatus())!=3 /*|| static_cast<int>(fitResult->CovMatrixStatus())!=2*/)
  {
    if ( (0.5*fitTotal->GetParameter(iParKPsip) <= fitTotal->GetParError(iParKPsip))  ){ //kPsi'
      std::cout << "------------------------------------------------" << std::endl;
      std::cout << "------- Setting Psi'norm= Psi' norm*0.5) -------" << std::endl;
      std::cout << "------------------------------------------------" << std::endl;
      bin = fHisto->FindBin(3.68);
      // fitTotal->SetParLimits(iParKPsip, 0.,fHisto->GetBinContent(bin)*0.9); // we further restrict the range of psi' norm
      fitTotal->SetParameter(iParKPsip, fHisto->GetBinContent(bin)*0.5);
    }
    if ( 0.5*fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ) {
      std::cout << "------------------------------------------------" << std::endl;
      std::cout << "-------         Setting p0=MAX/2)        -------" << std::endl;
      std::cout << "------------------------------------------------" << std::endl;
      fitTotal->SetParameter(0, fHisto->GetMaximum()*2.); // kVWG
    }

    std::cout << "================================\\" << std::endl;
    std::cout << "======== Refitting again =======\\" << std::endl;
    std::cout << "================================\\" << std::endl;
    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

     //Check if there are poles somewhere/
    // if(iLastParBkg == 6 )CheckRoots(fitResult,fitTotal,3,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6),fitOption) ;
    //if(iLastParBkg == 5 )CheckRoots(fitResult,fitTotal,2,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),0.,fitOption);
  }

  if ( (static_cast<int>(fitResult) && (static_cast<int>(fitResult)!=4000||static_cast<int>(fitResult)!=0)) || static_cast<int>(fitResult->CovMatrixStatus())!=3 /*|| static_cast<int>(fitResult->CovMatrixStatus())!=2*/) {
    std::cout << "============================================================================================\\" << std::endl;
    std::cout << "======== Refitting bkg again (setting range rejected 2.5-3.7, and fit range 1.7-4.5) =======\\" << std::endl;
    std::cout << "============================================================================================\\" << std::endl;
    SetFitRejectRange(2.5,3.7);
    TFitResultPtr fitResultInit = fHisto->Fit(bckInit,fitOption,"",1.7,4.5);
    //SetFitRejectRange();

    for ( Int_t i = 0; i < iLastParBkg+1 ; ++i ) fitTotal->SetParameter(i, bckInit->GetParameter(i)); //set initial background parameters

    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

     //Check if there are poles somewhere/
    //if(iLastParBkg == 6 )CheckRoots(fitResult,fitTotal,3,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6),fitOption) ;
    // if(iLastParBkg == 5 )CheckRoots(fitResult,fitTotal,2,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),0.,fitOption);
  }

  if ( (static_cast<int>(fitResult) && (static_cast<int>(fitResult)!=4000||static_cast<int>(fitResult)!=0)) || static_cast<int>(fitResult->CovMatrixStatus())!=3 /*|| static_cast<int>(fitResult->CovMatrixStatus())!=2*/) {
    std::cout << "============================================================================================\\" << std::endl;
    std::cout << "======== Refitting bkg again (setting range rejected 2.5-3.5, and fit range 1.5-5)  ========\\" << std::endl;
    std::cout << "============================================================================================\\" << std::endl;

    SetFitRejectRange(2.5,3.5);
    TFitResultPtr fitResultInit = fHisto->Fit(bckInit,fitOption,"",1.5,5.);
    //SetFitRejectRange();

    for ( Int_t i = 0; i < iLastParBkg+1 ; ++i ) fitTotal->SetParameter(i, bckInit->GetParameter(i)); //set initial background parameters

    fitResult = fHisto->Fit(fitTotal,fitOption,"");
    std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

    //if(iLastParBkg == 6 )CheckRoots(fitResult,fitTotal,3,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6),fitOption) ;
    // if(iLastParBkg == 5 )CheckRoots(fitResult,fitTotal,2,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),0.,fitOption);
  }

  if ( (static_cast<int>(fitResult) && (static_cast<int>(fitResult)!=4000||static_cast<int>(fitResult)!=0)) || static_cast<int>(fitResult->CovMatrixStatus())!=3.)
  {

    for ( Int_t i = 0; i < iLastParBkg+1 ; ++i ) fitTotal->SetParameter(i, bckInit->GetParameter(i));

    if ( (0.5*fitTotal->GetParameter(iParKPsip) <= fitTotal->GetParError(iParKPsip))  ) { //kPsi'
      std::cout << "------------------------------------------------" << std::endl;
      std::cout << "------- Setting Psi'norm= Psi' norm*0.3) -------" << std::endl;
      std::cout << "------------------------------------------------" << std::endl;
      bin = fHisto->FindBin(3.68);
      // fitTotal->SetParLimits(iParKPsip, 0.,fHisto->GetBinContent(bin)*0.7); // we further restrict the range of psi' norm
      fitTotal->SetParameter(iParKPsip, fHisto->GetBinContent(bin)*0.3);
    }

    if ( 0.5*fitTotal->GetParameter(0) <= fitTotal->GetParError(0) ){
      std::cout << "------------------------------------------------" << std::endl;
      std::cout << "-------         Setting p0=MAX*0.        -------" << std::endl;
      std::cout << "------------------------------------------------" << std::endl;
      // fitTotal->SetParLimits(0,bckInit->GetParameter(0)*0.1,bckInit->GetParameter(0)*1.5);
      fitTotal->SetParameter(0, fHisto->GetMaximum()*0.6); // kVWG
    }

    std::cout << "================================\\" << std::endl;
    std::cout << "======== Refitting again =======\\" << std::endl;
    std::cout << "================================\\" << std::endl;    fitResult = fHisto->Fit(fitTotal,fitOption,"");

    std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
    std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;

    //if(iLastParBkg == 6 )CheckRoots(fitResult,fitTotal,3,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),fitTotal->GetParameter(6),fitOption) ;
    //if(iLastParBkg == 5 )CheckRoots(fitResult,fitTotal,2,fitTotal->GetParameter(3),fitTotal->GetParameter(4),fitTotal->GetParameter(5),0.,fitOption);
  }

  // if ( !static_cast<int>(fitResult) && static_cast<int>(fitResult->CovMatrixStatus())!=3.){
  //   // change fit option if the only problem is the error calculation
  //   std::cout << "//------- Erro estimation problem, changing fit option" << std::endl;
  //   std::cout << "//======== Refitting again =======\\" << std::endl;

  //   fitResult = fHisto->Fit(fitTotal,Form("%sM",fitOption),"");
  //   std::cout << "FitResult = " << static_cast<int>(fitResult) << std::endl;
  //   std::cout << "CovMatrixStatus = " << fitResult->CovMatrixStatus() << std::endl;
  // }

  if ( (static_cast<int>(fitResult) && (static_cast<int>(fitResult)!=4000||static_cast<int>(fitResult)!=0)) || static_cast<int>(fitResult->CovMatrixStatus())!=3.){
    std::cout << "===========================================================\\" << std::endl;
    std::cout << "======== Cannot fit properly, try something else... =======\\" << std::endl;
    std::cout << "===========================================================\\" << std::endl;
  }
}


void SetFitRejectRange(Double_t a, Double_t b)
{
  /// Set a range the fit function(s) can ignore
  if ( a <= TMath::Limits<Double_t>::Max() && b <= TMath::Limits<Double_t>::Max() ) TF1::RejectPoint();

  /*fRejectFitPoints = kFALSE;

  fFitRejectRangeLow = a;
  fFitRejectRangeHigh = b;
  if ( a <= TMath::Limits<Double_t>::Max() && b <= TMath::Limits<Double_t>::Max() )
  {
    fRejectFitPoints = kTRUE;
  }
  */
}
//------------------------------------------------------------------------------
Double_t FitFunctionBackgroundPol2(Double_t *x, Double_t *par)
{
  // pol2 3 params

  if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  //if (x[0] > 3.5 && x[0] < 3.7) TF1::RejectPoint();
  return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
}

//------------------------------------------------------------------------------
Double_t FitFunctionBackgroundPol3(Double_t *x, Double_t *par)
{
  // pol2 3 params

  if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  //if (x[0] > 3.5 && x[0] < 3.7) TF1::RejectPoint();
  return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0];
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

//____________________________________________________________________________
Double_t alphaCB2VWG(Double_t*x, Double_t* par)
{
  return FitFunctionSignalCrystalBallExtended(x, &par[4])/(FitFunctionSignalCrystalBallExtended(x, &par[4]) + FitFunctionBackgroundVWG(x,par));
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

Double_t FitFunctionMeanPtS2CB2VWGPOL3(Double_t *x, Double_t *par)
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


  return alphaCB2VWG(x,par)*par[12] + (1. - alphaCB2VWG(x,par))*FitFunctionBackgroundPol3(x,&par[13]);
}

//------------------------------------------------------------------------------
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
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(2.4,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(font,"xyz");  
  gStyle->SetTitleOffset(0.9,"x");  
  gStyle->SetTitleOffset(0.9,"y");  
  gStyle->SetTitleSize(0.06,"x");  
  gStyle->SetTitleSize(0.06,"yz");  
  gStyle->SetMarkerSize(2.3); 
  gStyle->SetPalette(1,0); 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(10);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetEndErrorSize(0);
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
