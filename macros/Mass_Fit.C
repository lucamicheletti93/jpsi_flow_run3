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

Double_t func2CB2VWG(Double_t *x, Double_t *par);
Double_t BackgroundVWG(Double_t *x, Double_t *par);
Double_t CrystalBallExtended(Double_t *x,Double_t *par);
void ProcessMinvFit(TH1D *fHisto, TFitResultPtr& fitResult, TF1* fitTotal, TF1* bckInit, const char* fitOption, Int_t iParKPsip, Int_t iLastParBkg);
void SetFitRejectRange(Double_t a, Double_t b);
void Mass_Fit()
{

   TPaveText *display1 = new TPaveText(4.2867031,3915.316,5.0608323,6109.354,"BLARC");
  display1->SetTextFont(42);
  display1->SetTextSize(0.041);
  display1->SetTextColor(kBlack);
  display1->SetBorderSize(0);
  display1->SetFillColor(0);

  TPaveText *display2 = new TPaveText(4.3318969,963.6719,4.8850661,1293.534,"BLARC");
  display2->SetTextFont(42);
  display2->SetTextSize(0.041);
  display2->SetTextColor(kBlack);
  display2->SetBorderSize(0);
  display2->SetFillColor(0);


  TPaveText *display3 = new TPaveText(4.03318969,1300.6719,4.850661,1893.534,"BLARC");
  display3->SetTextFont(62);
  display3->SetTextSize(0.03);
  display3->SetTextColor(kBlack);
  display3->SetBorderSize(0);
  display3->SetFillColor(0);


  /*
  TFile *file = new TFile("v2_results.root");
  if (file->IsOpen()) cout << "File opened successfully" << endl;
  
  TH1D *h_Mass_SEPM = (TH1D*) file->Get("histProjMass");
  h_Mass_SEPM->SetMarkerStyle(20);
  h_Mass_SEPM->SetMarkerColor(kBlack);
  h_Mass_SEPM->SetMarkerSize(0.7);
  h_Mass_SEPM->Draw();
  */

  TFile *file = new TFile("AnalysisResults_Mass.root");
  TList *l = (TList*) file->Get("analysis-same-event-pairing/output");
  TList *l2 = (TList*) file->Get("analysis-event-mixing/output");
  TList *l_SEPM = (TList*)l->FindObject("PairsMuonSEPM_muonLowPt10SigmaPDCA");
  TList *l_MEPM = (TList*)l2->FindObject("PairsMuonMEPM_muonLowPt10SigmaPDCA");
  TH1D *h_Mass_SEPM = (TH1D*)l_SEPM->FindObject("Mass");
  TH1D *h_Mass_MEPM = (TH1D*)l_MEPM->FindObject("Mass");

  TFile *fOut = new TFile("mass_results.root", "RECREATE");
  h_Mass_SEPM ->GetXaxis()->SetRangeUser(1.5, 6);
  h_Mass_SEPM -> Write();
  fOut -> Close();

  h_Mass_SEPM->SetMarkerStyle(20);
  h_Mass_SEPM->SetMarkerColor(kBlack);
  h_Mass_SEPM->SetLineColor(kBlack);
  h_Mass_SEPM->SetMarkerSize(0.4);
 
  h_Mass_MEPM->SetMarkerStyle(20);
  h_Mass_MEPM->SetMarkerColor(kBlue);
  h_Mass_MEPM->SetLineColor(kBlue);
  h_Mass_MEPM->SetMarkerSize(0.4);
   
  const Double_t xmin=2.2;
  const Double_t xmax=4.7;
  
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
  fitTotal->SetParLimits(5, 3.07, 3.2);
  fitTotal->SetParameter(6, 0.08);
  fitTotal->SetParLimits(6, 0.05, 0.15);
  
//  r = FitJpsi2CB2VWG(*hminv,0.93,5.59,2.32,3.39);

  
    fitTotal->SetParameter(7,0.9);
    fitTotal->SetParLimits(7,0.1,10.0);
 
    fitTotal->SetParameter(8,5.0);
    fitTotal->SetParLimits(8,0.0,10.0);
  
    fitTotal->SetParameter(9, 2.0);
    fitTotal->SetParLimits(9,0.1,10.0);
  
    fitTotal->SetParameter(10,3.0);
    fitTotal->SetParLimits(10,0.0,10.0);
  
    fitTotal->SetParameter(11, 10.);

  const char* fitOption = "QSER"; //+";
  
  //TFitResultPtr fitResult = hfit->Fit(fitTotal,fitOption,"");

  double SE = h_Mass_SEPM->Integral(h_Mass_SEPM->GetXaxis()->FindBin(1.0),h_Mass_SEPM->GetXaxis()->FindBin(2.5)) + h_Mass_SEPM->Integral(h_Mass_SEPM->GetXaxis()->FindBin(3.72),h_Mass_SEPM->GetXaxis()->FindBin(5.0));

  double ME = h_Mass_MEPM->Integral(h_Mass_MEPM->GetXaxis()->FindBin(1.0),h_Mass_MEPM->GetXaxis()->FindBin(2.5)) + h_Mass_MEPM->Integral(h_Mass_MEPM->GetXaxis()->FindBin(3.72),h_Mass_MEPM->GetXaxis()->FindBin(5.0));

  cout<<" Normalization factor "<<ME/SE <<endl;

  h_Mass_SEPM->GetXaxis()->SetRangeUser(2, 5);
  h_Mass_SEPM->Fit(fitTotal,"REM");
  h_Mass_SEPM->Draw("E");
  //h_Mass_MEPM->Rebin(2);
  //h_Mass_MEPM->Scale(0.165564);
  //h_Mass_MEPM->Draw("Esame");
  
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
  fitFctCB2->Draw("same");


  TF1* bck2 = new TF1("bck2",BackgroundVWG,xmin,xmax,4);
  for ( Int_t i = 0; i < 4; ++i )
  {
    bck2->SetParameter(i,fitTotal->GetParameter(i));
  }
  bck2->SetLineColor(kGreen);
  bck2->SetLineStyle(2);
  //bck2->Draw("same");

  
  TFitResultPtr r = h_Mass_SEPM->Fit(fitTotal,"RESLI");
  TMatrixDSym cov = r->GetCovarianceMatrix();



  
  //==============================================
  ///////////////////////////////////////////////////////////////
  std::cout << "FitResult = " << static_cast<int>(r) << std::endl;
  std::cout << "CovMatrixStatus = " << r->CovMatrixStatus() << std::endl;

  
  if ( ( static_cast<int>(r) && static_cast<int>(r)!=4000 ) ||  static_cast<int>(r->CovMatrixStatus())!=3 ) ProcessMinvFit(h_Mass_SEPM,r,fitTotal,bck2,"RESLI",12,4);
  //Set("FitResult",static_cast<int>(r),0);
  //Set("CovMatrixStatus",static_cast<int>(r->CovMatrixStatus()),0);
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n");
  printf(" Fit Status : %d <-> Cov. Mat. : %d ",static_cast<int>(r),r->CovMatrixStatus());
  printf("\n -_-_-_-_-_-_-_-_-_-_-_-_-_-_\n\n");
  

  cout<<"Double_t Fit_Par[12] = {";
  for(int i=0;i<12;i++) {
    cout<<fitTotal->GetParameter(i);
    if(i!=11) cout<<", ";
  }
  cout<<" }; "<<endl;


  cout<<"mass Jpsi" << fitTotal->GetParameter(5) << " +/-" << fitTotal->GetParError(5) << endl;
  cout<<"sigma Jpsi" << fitTotal->GetParameter(6) << " +/-" << fitTotal->GetParError(6) << endl;;
  cout << " mass of psi' "<< (3.68609+(fitTotal->GetParameter(5)-3.096916)/3.096916*3.68609) << endl;
  cout << "sigma psi' " << (fitTotal->GetParameter(6)/3.096916*3.68609) << endl;

  //====================J/Psi Signal Estimation==========================================================
 
  // Set("NofJPsi",njpsi,nerr);
  
  double m =fitTotal->GetParameter(5);
  double m_err =fitTotal->GetParError(5);
  double s =fitTotal->GetParameter(6);
  double s_err =fitTotal->GetParError(6);
  double njpsi3s = fitTotal->Integral(m-3*s,m+3*s)/h_Mass_SEPM->GetBinWidth(1);
  double nerr3s = fitTotal->IntegralError(m-3*s,m+3*s,fitTotal->GetParameters(),cov.GetMatrixArray())/h_Mass_SEPM->GetBinWidth(1);
  cout<< "No of J/Psi == " <<  njpsi3s <<"+/-" << nerr3s<<endl;
  Double_t chi=fitTotal->GetChisquare();
  Double_t ndf=fitTotal->GetNDF();
  Double_t ch_ndf=chi/ndf;
  
  //===============================================================================
  //Computation of bin significance and signal over background
  
  double nbck3s = bck2->Integral(m-3.*s,m+3.*s)/h_Mass_SEPM->GetBinWidth(1);
  //double nbck3sErr = exp2->IntegralError(m-3*s,m+3*s,r2->GetParams(),r2->GetCovarianceMatrix().GetMatrixArray() )/hT->GetBinWidth(1);
  double nbck3sErr = bck2->IntegralError(m-3*s,m+3*s,bck2->GetParameters(),cov.GetMatrixArray())/h_Mass_SEPM->GetBinWidth(1);
  double sOverB3s = njpsi3s / nbck3s;
  double sOverB3sErr = sOverB3s*TMath::Sqrt(TMath::Power(nerr3s/njpsi3s,2.) + TMath::Power(nbck3sErr/nbck3s,2.));

  cout<<" S/B ==="<< sOverB3s << "+/-" << sOverB3sErr <<endl;
  //Set("SignalOverBkg3s",sOverB3s,sOverB3sErr);
  
  double sig = njpsi3s/TMath::Sqrt(njpsi3s + nbck3s);
  double sigErr = TMath::Sqrt( TMath::Power((1. - (1./2.)*njpsi3s/(njpsi3s + nbck3s) )*nerr3s/TMath::Sqrt(njpsi3s + nbck3s),2.) +
                               TMath::Power(njpsi3s*nbck3sErr/(2.*TMath::Power(njpsi3s + nbck3s,3./2.)),2.) );


  cout<<" s/sqrt(s+b) " <<sig << " +/- "<<sigErr<<endl;  
  
  //==========================================================================================================

  TText *text0 = display1->AddText(Form("#chi^{2}/ndf =%f",ch_ndf));
  TText *text1 = display1->AddText(Form("M_{J/#psi}(3#sigma) =%f +/- %f GeV",m,m_err));
  TText *text2 = display1->AddText(Form("#sigma_{J/#psi}(3#sigma) = %f +/- %f GeV",s,s_err));
  TText *text3 = display2->AddText(Form("N_{J/#psi}(3#sigma) = %f +/- %f", njpsi3s, nerr3s));
  TText *text4 = display2->AddText(Form("#frac{S}{B}_{J/#psi}(3#sigma) = %f +/- %f",sOverB3s,sOverB3sErr));
  //TText *text5 = display2->AddText(Form("#frac{S}{#sqrt{S+B}}_{J/#psi}(3#sigma) = %f +/- %f",sig,sigErr));

  TLegend *legend2_1 = new TLegend(0.6125205,0.5806881,0.8948445,0.7828001,NULL,"brNDC");
  legend2_1->SetTextFont(42);
  legend2_1->SetTextSize(0.04);
  legend2_1->SetLineColor(0);
  legend2_1->SetFillColor(0);
  legend2_1->AddEntry(h_Mass_SEPM,"Same event dimuon","p");
  legend2_1->AddEntry(h_Mass_MEPM,"Mixed event dimuon","p");
  legend2_1->Draw();


  display1->Draw("same");
  display2->Draw("same");
  
  
  
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
    return BackgroundVWG(x, par) + CrystalBallExtended(x, &par[4]) + CrystalBallExtended(x, par2);
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
