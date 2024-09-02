#include "TProfile.h"
#include <fstream>
#include <iostream>
#include <string>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TH2D.h"
#include "TVirtualFitter.h"
#include "TPaveText.h"
#include <TLatex.h>
#include <TStyle.h>
#include"TMatrixDSym.h"
#include"TFitResult.h"
#include"TROOT.h"
Bool_t reject;
void SetStyle(Bool_t graypalette=kTRUE);
Double_t FitFunctionBackgroundPol2(Double_t *x, Double_t *par);
Double_t FitFunctionSignalCrystalBallExtended(Double_t *x,Double_t *par);
Double_t alphaCB2VWG(Double_t *x, Double_t *par);
Double_t FitFunctionBackgroundVWG(Double_t *x, Double_t *par);
Double_t FitFunctionMeanPtS2CB2VWGPOL2(Double_t *x, Double_t *par);
void SetFitRejectRange(Double_t a, Double_t b);

Double_t Fit_Par[12] = {3.97736e+07, 0.000772307, 0.44601, 0.000106707, 14438.4, 3.07397, 0.0770682, 9.75573, 6.75144, 9.81753, 7.71034, -113.444 };
  
void v2Fit()
{
  TCanvas *c = new TCanvas("c","c",20,20,600,600);
  //c->SetLogy();
  //SetStyle();
  gStyle->SetOptStat(111);
  //gStyle->SetOptFit(0);
  //gStyle->SetOptTitle(0);
  //gStyle->SetCanvasColor(10);
  //gStyle->SetFrameFillColor(10);


   TPaveText *display1 = new TPaveText(4.2867031,3915.316,5.0608323,6109.354,"BLARC");
  display1->SetTextFont(62);
  display1->SetTextSize(0.03);
  display1->SetTextColor(kBlack);
  display1->SetBorderSize(0);
  display1->SetFillColor(0);

  TPaveText *display2 = new TPaveText(4.3318969,963.6719,4.8850661,1293.534,"BLARC");
  display2->SetTextFont(62);
  display2->SetTextSize(0.03);
  display2->SetTextColor(kBlack);
  display2->SetBorderSize(0);
  display2->SetFillColor(0);


  TPaveText *display3 = new TPaveText(4.03318969,1300.6719,4.850661,1893.534,"BLARC");
  display3->SetTextFont(62);
  display3->SetTextSize(0.03);
  display3->SetTextColor(kBlack);
  display3->SetBorderSize(0);
  display3->SetFillColor(0);

  
  TH2F *f0=new TH2F("f","M_{#mu#mu}(GeV/c^{2})",120,2.0,5.0,1.e+07,0,1.e+07);
  f0->GetXaxis()->SetTitle("M_{#mu#mu}(GeV/c^{2})");
  f0->GetYaxis()->SetTitle("dN/dM_{#mu#mu} (c^{2}/GeV)");
  f0->GetXaxis()->CenterTitle();
  f0->GetYaxis()->CenterTitle();

  Double_t fitRangeLow = 2.2;
  Double_t fitRangeHigh = 4.7;

  TFile *file = new TFile("v2_results.root");
  if (file->IsOpen()) cout << "File opened successfully" << endl;

  TFile *file2 = new TFile("v2_Bkg.root");
  if (file->IsOpen()) cout << "File opened successfully" << endl;

  TH1D *hT = (TH1D*) file->Get("histV2SP");
  hT->SetMarkerStyle(20);
  hT->SetMarkerColor(kBlack);
  hT->SetLineColor(kBlack);
  hT->SetMarkerSize(0.7);
  hT->Draw();
  
  TProfile *hT_bkg = (TProfile*) file2->Get("v2Bkg_vs_mass;1");
  hT_bkg->SetMarkerStyle(20);
  hT_bkg->SetMarkerColor(kBlue);
  hT_bkg->SetLineColor(kBlue);
  hT_bkg->SetMarkerSize(0.7);
  hT_bkg->Scale(0.12);
  hT_bkg->Rebin(1);
  //hT_bkg->Draw("psame");
  ///////////////Fitting //////////////////////////////////

  TF1* bck = new TF1("bck",FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3);
  bck->SetLineColor(kBlue);
  bck->SetLineStyle(9);
  bck->SetParameters(20.,1.,1.);
  //bck->SetParLimits(0, 0.,50.0);
  //SetFitRejectRange(3.0,4.0);
  hT->Fit(bck,"SRL","",fitRangeLow,fitRangeHigh);
  //bck->Draw("same");
  //SetFitRejectRange();


  TF1* fitMeanpt = new TF1("fitMeanpt",FitFunctionMeanPtS2CB2VWGPOL2,fitRangeLow,fitRangeHigh,16);
  fitMeanpt->SetLineColor(kRed);
  fitMeanpt->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
  fitMeanpt->SetParName(11,"kPsiP");
  fitMeanpt->SetParName(12,"v_{2} JPsi");
  fitMeanpt->SetParName(13,"v_{2} BG0");
  fitMeanpt->SetParName(14,"v_{2} BG1");
  fitMeanpt->SetParName(15,"v_{2} BG2");


  for( Int_t i = 0; i < 12; ++i ) fitMeanpt->FixParameter(i,Fit_Par[i]);

  /*
  fitMeanpt->FixParameter(0,1.02600e+04);
  fitMeanpt->FixParameter(1,2.36336e+00);
  fitMeanpt->FixParameter(2,5.42282e-01);
  fitMeanpt->FixParameter(3,2.66333e-01);
  fitMeanpt->FixParameter(4,1.30535e+03);
  fitMeanpt->FixParameter(5,3.08194e+00);
  fitMeanpt->FixParameter(6,6.25909e-02);
  fitMeanpt->FixParameter(7,8.64747e-01);
  fitMeanpt->FixParameter(8,4.06186e-01);
  fitMeanpt->FixParameter(9,1.79065e+00);
  fitMeanpt->FixParameter(10,4.94294e-02);
  fitMeanpt->FixParameter(11,-7.66891e+00);
  */

  fitMeanpt->SetParameter(12, 0.04);
  fitMeanpt->SetParLimits(12, 0.002,0.08);

 for ( Int_t i = 0; i < 3; ++i )
  {
    fitMeanpt->SetParameter(i + 13, bck->GetParameter(i));
  }


  TFitResultPtr r = hT->Fit(fitMeanpt,"REMI");
  

  TF1* bck2 = new TF1("bck2",FitFunctionBackgroundPol2,fitRangeLow,fitRangeHigh,3);
  for ( Int_t i = 0; i < 3; ++i )
  {
    bck2->SetParameter(i, fitMeanpt->GetParameter(i + 13));
  }
  bck2->SetLineColor(kBlue);
  bck2->SetLineStyle(2);
  bck2->Draw("same");


  //const char* fitOption = "SER";
  //TFitResultPtr fitResult = hT->Fit(fitMeanpt,fitOption,"");
  //TFitResultPtr r = hT->Fit(fitMeanpt,"REMI");
  //fitMeanpt->Draw();
  //hT->Draw();
 
  ///////////////////////////////////////////////////////////////

  TLegend *legend2_1 = new TLegend(0.6125205,0.5806881,0.8948445,0.7828001,NULL,"brNDC");
  legend2_1->SetTextFont(42);
  legend2_1->SetTextSize(0.04);
  legend2_1->SetLineColor(0);
  legend2_1->SetFillColor(0);
  legend2_1->AddEntry(hT,"Same event dimuon","p");
  legend2_1->AddEntry(hT_bkg,"Mixed event dimuon","p");
  legend2_1->Draw();



  c->Update();
  c->SaveAs("Spectra6.root");   
 
}

//-----------------------------------------------------------------------
void SetFitRejectRange(Double_t a, Double_t b)
{
  /// Set a range the fit function(s) can ignore

  Bool_t fRejectFitPoints = kFALSE;

  Double_t fFitRejectRangeLow = a;
  Double_t fFitRejectRangeHigh = b;
  if ( a <= TMath::Limits<Double_t>::Max() && b <= TMath::Limits<Double_t>::Max() )
  {
    fRejectFitPoints = kTRUE;
  }
}
//------------------------------------------------------------------------------
Double_t FitFunctionBackgroundPol2(Double_t *x, Double_t *par)
{
  // pol2 3 params

  if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  if (x[0] > 3.5 && x[0] < 3.7) TF1::RejectPoint();
  return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
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

  if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  if (x[0] > 3.5 && x[0] < 3.7) TF1::RejectPoint();
  Double_t sigma = par[2]+par[3]*((x[0]-par[1])/par[1]);
  return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));
}

//-------------------------------------------------------------------------------

Double_t FitFunctionMeanPtS2CB2VWGPOL2(Double_t *x, Double_t *par)
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

//------------------------------------------------------------------------------

void SetStyle(Bool_t graypalette) {
  //cout << "Setting style!" << endl;
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(1);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kBlack);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.019,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetTitleOffset(1.31,"y");
  gStyle->SetTitleOffset(1.05,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(28);
  gStyle->SetTextFont(42);
  gStyle->SetTickLength(0.03,"X");
  gStyle->SetTickLength(0.03,"Y"); 
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
}
