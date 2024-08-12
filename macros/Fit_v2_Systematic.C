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
Double_t FitFunctionBackgroundPol2Exp(Double_t *x, Double_t *par);
Double_t FitFunctionNA60New(Double_t *x,Double_t *par);
Double_t FitFunctionFlowS2NA60NEWPOL2EXPPOL2(Double_t *x, Double_t *par);
Double_t FitFunctionFlowS2NA60NEWVWGPOL2(Double_t *x, Double_t *par);
Double_t alphaNA60NEWVWG(Double_t*x, Double_t* par);
Double_t alphaNA60NEWPOL2EXP(Double_t*x, Double_t* par);
Double_t FitFunctionBackgroundPol3Cheb(Double_t *x, Double_t *par);
Double_t FitFunctionFlowS2CB2VWGPol3Cheb(Double_t *x, Double_t *par);

TF1* v2Fit(TProfile *v2prof, Double_t par[12],Double_t parLimit[3]);

Double_t xmin= 2.5;
Double_t xmax= 4.8;

const Double_t v2xmin= 2.2;
const Double_t v2xmax= 4.8;

TFile *fOut;
TFile *fv2fit;
Double_t Jpsiv2 ;
Double_t Jpsiv2_err ;
Double_t v2_ch_ndf;


const int nFits_fun = 5;
string mass_sig[] = {"CB2","CB2","CB2","CB2","NA60","NA60"};
string mass_bkg[] = {"VWG_data","Pol4Exp_data","VWG_MC","Pol4Exp_MC","Pol4Exp_MC","VWG_MC"};

string func_sig[] = {"CB2","CB2","CB2","CB2","NA60","NA60"};
string func_bkg[] = {"VWG","Pol4Exp","VWG","Pol4Exp","Pol4Exp","VWG_MC"};

const int nRanges = 3;
double min_range[] = {2.3,2.4,2.5};
double max_range[] = {4.7,4.6,4.5};

double resolution = 1.3093;

void Fit_v2_Systematic()
{
  TCanvas *c1 = new TCanvas("c1", "c1",75,70,1309,696);
  TProfile *prof_v2EP[20];
  TProfile *prof_MEv2EP[20];
  TH1D *hT_mass[20];
  TF1 *fitv2[20][20][20];
  
  double par[16];

  string sig_plus_bkg[] = {"bkg","aa","bb","cc","sig_Jpsi","mean_Jpsi","width_Jpsi","a_Jpsi","b_Jpsi","c_Jpsi","d_Jpsi"};
                           //("kVWG","mVWG","sVWG1","sVWG2","kPsi","mPsi","sPsi","alPsi","nlPsi","auPsi","nuPsi");

 string signa60_plus_bkg[] = {"bkg","aa","bb","cc","sig_Jpsi","mean_Jpsi","width_Jpsi","a_Jpsi","b_Jpsi","c_Jpsi","d_Jpsi","e_Jpsi","f_Jpsi","g_Jpsi","h_Jpsi"};
                         //("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi", "p1LJPsi","p2LJPsi","p3LJPsi","p1RJPsi","p2RJPsi","p3RJPsi","aLJPsi","aRJPsi");
  
  
  double minCentrBins[] = {10};
  double maxCentrBins[] = {50};

  const int nPtBins = 5;
  double minPtBins[] = {0,1.0,2.0,4.0,5.0};
  double maxPtBins[] = {1.0,2.0,4.0,5.0,10.0};

  //const int nPtBins = 1;
  //double minPtBins[] = {1.0};
  //double maxPtBins[] = {2.0};

  Double_t mean_v2[nPtBins];
  Double_t mean_v2er[nPtBins];

  TH1F *hSys[nPtBins]; 
  TCanvas *cPt[nPtBins];
  //cPt[nPtBins] = new TCanvas("cPt", "cPt",75,70,1309,696);

 
  
  //v2 Fitting Limits ======
  Double_t parLimit0[3];
  parLimit0[0] = 0.0040; 
  parLimit0[1] = 0.0037;
  parLimit0[2] = 0.0045;
  
  Double_t parLimit1[3];
  parLimit1[0] = 0.035; 
  parLimit1[1] = 0.032;
  parLimit1[2] = 0.039;

  Double_t parLimit2[3];
  parLimit2[0] = 0.0828; 
  parLimit2[1] = 0.078;
  parLimit2[2] = 0.085;

  Double_t parLimit3[3];
  parLimit3[0] = 0.115; 
  parLimit3[1] = 0.11;
  parLimit3[2] = 0.122;

  Double_t parLimit4[3];
  parLimit4[0] = 0.159; 
  parLimit4[1] = 0.152;
  parLimit4[2] = 0.169;
  

  int nSysbins = nFits_fun*nRanges;
  

  TFile *f_hist = new TFile("Histograms_Fullpass3matchedMchMid_centr_10_50.root");
  if (f_hist->IsOpen()) cout << "File opened successfully" << endl;
  
 for(int iPt = 0;iPt<nPtBins;iPt++) //Loop over mass and v2 histograms in different pT ranges
    {
  hSys[iPt] = new TH1F("hSys", " ",nSysbins,0,nSysbins);
  int ncombination = 0;
  Double_t sum_v2 =0.0;
  Double_t sum_v2er =0.0;
 
  for(int kfun = 0; kfun < nFits_fun; kfun++)  //Loop over Signal and backgrounds shape
    {
      TFile *f_fit = new TFile(Form("mass_fit_results/Pt_%1.0f_%1.0f/multi_trial_%s_%s_tails.root",minPtBins[iPt], maxPtBins[iPt],mass_sig[kfun].data(),mass_bkg[kfun].data()));
      if (f_fit->IsOpen()) cout << "File opened successfully" << endl;
      //f_fit->ls();

      for(int iR = 0; iR < nRanges ; iR++) //Loop over Signal Range
	{
 
	   prof_v2EP[iPt] = (TProfile*) f_hist->Get(Form("histV2SEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
	   prof_MEv2EP[iPt] = (TProfile*) f_hist->Get(Form("histV2MEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
	   hT_mass[iPt] = (TH1D*)f_fit->Get(Form("fit_results_%s_%s__%0.1f_%0.1f_histMassSEPM_%1.0f_%1.0f__%1.0f_%1.0f",func_sig[kfun].data(),func_bkg[kfun].data(), min_range[iR],max_range[iR],minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));

           //NA60 
	   if(kfun > 3) {for(int i = 0; i < 15 ; i++) par[i] = hT_mass[iPt]->GetBinContent(hT_mass[iPt]->GetXaxis()->FindBin(Form("%s",signa60_plus_bkg[i].data())));
	   par[15] = 1.0;}

	   //CB2
	   else{for(int i = 0; i < 11 ; i++) par[i] = hT_mass[iPt]->GetBinContent(hT_mass[iPt]->GetXaxis()->FindBin(Form("%s",sig_plus_bkg[i].data())));
           par[11] = 1.0;}
	   //hT_mass[iPt]->Draw();

	   prof_v2EP[iPt]->SetName(Form("fitv2_Sig_Bkg:%s_%s__FitRange:%0.1f_%0.1f__Pt:%1.0f_%1.0f__CENT:%1.0f_%1.0f",mass_sig[kfun].data(),mass_bkg[kfun].data(), min_range[iR],max_range[iR],minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
	   prof_v2EP[iPt]->SetTitle(Form("fitv2_Sig_Bkg:%s_%s__FitRange:%0.1f_%0.1f",mass_sig[kfun].data(),mass_bkg[kfun].data(), min_range[iR],max_range[iR]));
	  
	  if(iPt ==0) fitv2[iPt][kfun][iR] = v2Fit(prof_v2EP[iPt], par,parLimit0);
	  if(iPt ==1) fitv2[iPt][kfun][iR] = v2Fit(prof_v2EP[iPt], par,parLimit1);
	  if(iPt ==2) fitv2[iPt][kfun][iR] = v2Fit(prof_v2EP[iPt], par,parLimit2);
	  if(iPt ==3) fitv2[iPt][kfun][iR] = v2Fit(prof_v2EP[iPt], par,parLimit3);
	  if(iPt ==4) fitv2[iPt][kfun][iR] = v2Fit(prof_v2EP[iPt], par,parLimit4);
	   
	   fitv2[iPt][kfun][iR]->SetName(Form("fitv2_Sig_Bkg:%s_%s__FitRange:%0.1f_%0.1f__Pt:%1.0f_%1.0f__CENT:%1.0f_%1.0f",mass_sig[kfun].data(),mass_bkg[kfun].data(), min_range[iR],max_range[iR],minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
	   fitv2[iPt][kfun][iR]->SetTitle(Form("fitv2_Sig_Bkg:%s_%s__FitRange:%0.1f_%0.1f",mass_sig[kfun].data(),mass_bkg[kfun].data(), min_range[iR],max_range[iR]));
           
	   hSys[iPt]->GetXaxis()->SetBinLabel(ncombination+1,Form("fitv2_Sig_Bkg:%s_%s__FitRange:%0.1f_%0.1f",mass_sig[kfun].data(),mass_bkg[kfun].data(), min_range[iR],max_range[iR]));
           hSys[iPt]->SetBinContent(ncombination+1,fitv2[iPt][kfun][iR]->GetParameter(12));

	   sum_v2 = sum_v2 + fitv2[iPt][kfun][iR]->GetParameter(12);
           sum_v2er = sum_v2er + fitv2[iPt][kfun][iR]->GetParError(12);

	   ncombination++;
	   cout<<" ncombination "<<ncombination<<endl;
	   
	  }

	}

   mean_v2[iPt] = sum_v2/ncombination;
   mean_v2er[iPt] = sum_v2er/ncombination;

  cout<<" Mean + Statistical "<<mean_v2[iPt]<<" +/- "<<mean_v2er[iPt]<<endl;
    
  hSys[iPt]->SetTitle(Form("SystV2SEPM_Pt:%1.0f_%1.0f__CENT:%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
  hSys[iPt]->SetName(Form("SystV2SEPM_Pt:%1.0f_%1.0f__CENT:%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
  hSys[iPt]->SetMarkerColor(kRed);
  hSys[iPt]->SetMarkerStyle(20);
  hSys[iPt]->SetMarkerSize(0.9);

  //hSys[iPt]->Draw();
  
    }


 for(int iPt = 0;iPt<nPtBins;iPt++) //Loop over v2_signal and background
   {
     Double_t sum2 =0.0;
     int NN = 0;
     for(int kfun = 0; kfun < nFits_fun; kfun++)  //Loop over Signal shape
       {
	 for(int iR = 0; iR < nRanges ; iR++) //Loop over Signal Range
	   {

	     Double_t dev = (fitv2[iPt][kfun][iR]->GetParameter(12) - mean_v2[iPt]);
	     //cout<<fitv2[iPt][kfun][iR]->GetParameter(12) <<"  "<< mean_v2[iPt]<<"  "<<fitv2[iPt][kfun][iR]->GetParameter(12) - mean_v2[iPt]<<endl;
	     sum2 = sum2 + dev*dev;
	     NN++;
	   }
       }

     //cout<<" Sum "<<sum2<<" NN "<<NN<<endl;
     Double_t sys = TMath::Sqrt(sum2/NN);

     Double_t per_stat = (mean_v2er[iPt]/mean_v2[iPt])*100;
     Double_t per_sys = (sys/mean_v2[iPt])*100;

     cPt[iPt] = new TCanvas("c1", "c1",75,70,1309,696);
     cPt[iPt]->SetName(Form("SystV2SEPM_Pt:%1.0f_%1.0f__CENT:%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
     cPt[iPt]->cd();
     hSys[iPt]->Draw("p");

     TLine *line = new TLine(0,mean_v2[iPt],NN,mean_v2[iPt]);
     line->SetLineStyle(1);
     line->SetLineColor(kBlue);
     line->SetLineWidth(2);
     line->Draw();

     TLine *line1 = new TLine(0,sys+mean_v2[iPt],NN,sys+mean_v2[iPt]);
     line1->SetLineStyle(2);
     line1->SetLineColor(kBlue);
     line1->SetLineWidth(2);
     line1->Draw();

     TLine *line2 = new TLine(0,mean_v2[iPt]-sys,NN,mean_v2[iPt]-sys);
     line2->SetLineStyle(2);
     line2->SetLineColor(kBlue);
     line2->SetLineWidth(2);
     line2->Draw();

     TPaveText *display1 = new TPaveText(0.1963037,0.7368947,0.8036963,0.795,"blNDC");
     display1->SetTextFont(42);
     display1->SetTextSize(0.046);
     display1->SetTextColor(kBlack);
     display1->SetBorderSize(0);
     display1->SetFillColor(0);

    
     TText *text1 = display1->AddText(Form(" v_{2} {J/#psi} =%.3f +/- %.3f (%.2f %%) +/- %.3f (%.2f %% ) ",mean_v2[iPt],mean_v2er[iPt],per_stat,sys,per_sys));
     display1->Draw("same");
   
     cout<<" Average ==== "<<mean_v2[iPt]<< " +/- " <<mean_v2er[iPt]<<" ("<<per_stat<<"\%)"<<" +/- "<<sys<<" ("<<per_sys<<"\%)"<<endl;

   }

 
 fOut = new TFile("Systematic_Cent10-50.root","RECREATE");
 cPt[0]->Write();
 cPt[1]->Write();
 cPt[2]->Write();
 cPt[3]->Write();
 cPt[4]->Write();
 fOut->Close();
 
 
}


TF1* v2Fit(TProfile *v2prof, Double_t par[12],Double_t parLimit[3])
{
  v2prof->Rebin(2);
  TCanvas *c3 = new TCanvas("c3", "c3",75,70,1309,696);
  //v2prof->Scale(1.405);
  v2prof->Scale(resolution);
  v2prof->GetXaxis()->SetRangeUser(xmin,xmax);
  v2prof->Draw();

  //=============================== v2 Fitting =================
  TF1* bck = new TF1("bck",FitFunctionBackgroundPol2,xmin,xmax,3);
  bck->SetLineColor(kBlue);
  bck->SetLineStyle(9);
  //bck->SetParameters(0.0191306,-0.00836616,0.00102991);
  bck->SetParameters(0.106962,-0.0252894,0.0022179);
  //SetFitRejectRange(2.9,3.3);
  v2prof->Fit(bck,"SRLI","",v2xmin,v2xmax);
  //bck->Draw("same");
  //SetFitRejectRange();

 
  std::string Fit_mainStr = v2prof->GetName();
  cout<<" Fit Function ========================= "<<Fit_mainStr.data()<<endl;
  std::string Fit_subStr1 = "CB2_VWG";
  std::string Fit_subStr2 = "CB2_Pol4Exp";
  std::string Fit_subStr3 = "NA60_Pol4Exp";
  std::string Fit_subStr4 = "NA60_VWG";
  size_t found_cb2_vwg = Fit_mainStr.find(Fit_subStr1);
  size_t found_cb2_pol4exp = Fit_mainStr.find(Fit_subStr2);
  size_t found_na60_pol4exp = Fit_mainStr.find(Fit_subStr3);
  size_t found_na60_vwg = Fit_mainStr.find(Fit_subStr4);

 TF1* fitv2;
 if (found_cb2_pol4exp != std::string::npos) 
    {
  cout<<"CB2_Pol4Exp EXIST"<<endl;
  fitv2 = new TF1("fitv2",FitFunctionFlowS2CB2POL2EXPPOL2,v2xmin,v2xmax,16);
  fitv2->SetParNames("pol0","pol1","pol2","exp","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
    }

 if (found_na60_vwg != std::string::npos) 
   {
     cout<<"NA60_VWG EXIST"<<endl;
     fitv2 = new TF1("fitv2",FitFunctionFlowS2NA60NEWVWGPOL2,v2xmin,v2xmax,20);
     fitv2->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","p1LJPsi","p2LJPsi","p3LJPsi","p1RJPsi");
     fitv2->SetParName(11,"p2RJPsi");
     fitv2->SetParName(12,"p3RJPsi");
     fitv2->SetParName(13,"aLJPsi");
     fitv2->SetParName(14,"aRJPsi");
     fitv2->SetParName(15,"kPsiP");
     fitv2->SetParName(16,"v_{2} JPsi");
     fitv2->SetParName(17,"v_{2} >BG0");
     fitv2->SetParName(18,"v_{2} BG1");
     fitv2->SetParName(19,"v_{2} BG2");
   }

 if (found_na60_pol4exp != std::string::npos) 
   {
     cout<<"NA60_EXPPol EXIST"<<endl;
     fitv2 = new TF1("fitv2",FitFunctionFlowS2NA60NEWPOL2EXPPOL2,v2xmin,v2xmax,20);
     fitv2->SetParNames("pol0","pol1","pol2","exp","kJPsi","mJPsi","sJPsi","p1LJPsi","p2LJPsi","p3LJPsi","p1RJPsi");
     fitv2->SetParName(11,"p2RJPsi");
     fitv2->SetParName(12,"p3RJPsi");
     fitv2->SetParName(13,"aLJPsi");
     fitv2->SetParName(14,"aRJPsi");
     fitv2->SetParName(15,"kPsiP");
     fitv2->SetParName(16,"v_{2} JPsi");
     fitv2->SetParName(17,"v_{2} BG0");
     fitv2->SetParName(18,"v_{2} BG1");
     fitv2->SetParName(19,"v_{2} BG2");
   }

 else
   {
  fitv2 = new TF1("fitv2",FitFunctionFlowS2CB2VWGPOL2,v2xmin,v2xmax,16);
  fitv2->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
   }

  fitv2->SetLineColor(kRed);
  fitv2->SetParName(11,"kPsiP");
  fitv2->SetParName(12,"v_{2} JPsi");
  fitv2->SetParName(13,"v_{2} BG0");
  fitv2->SetParName(14,"v_{2} BG1");
  fitv2->SetParName(15,"v_{2} BG2");

  for( Int_t i = 0; i < 12; ++i ) fitv2->FixParameter(i,par[i]);
  fitv2->SetParameter(12, parLimit[0]);
  fitv2->SetParLimits(12, parLimit[1],parLimit[2]);
  
  for ( Int_t i = 0; i < 3; ++i ) { fitv2->SetParameter(i + 13, bck->GetParameter(i));}

  TFitResultPtr r2 = v2prof->Fit(fitv2,"RESLI");

  TF1* bck2 = new TF1("bck2",FitFunctionBackgroundPol2,v2xmin,v2xmax,3);
  for ( Int_t i = 0; i < 3; ++i ) { bck2->SetParameter(i, fitv2->GetParameter(i + 13));}
  bck2->SetLineColor(kBlue);
  bck2->SetLineStyle(2);
  bck2->Draw("same");

 bck->SetName("Fit_mainStr.data()");
 fv2fit = new TFile("FitResults_v2-50.root","RECREATE");
 v2prof->Write();
 bck2->Write();
 fv2fit->Close();

  Jpsiv2 =  fitv2->GetParameter(12);
  Jpsiv2_err =  fitv2->GetParError(12);
  v2_ch_ndf = fitv2->GetChisquare()/fitv2->GetNDF();

  return fitv2;
}



void SetFitRejectRange(Double_t a, Double_t b)
{
  /// Set a range the fit function(s) can ignore
  if ( a <= TMath::Limits<Double_t>::Max() && b <= TMath::Limits<Double_t>::Max() ) TF1::RejectPoint();
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
Double_t FitFunctionBackgroundPol3Cheb(Double_t *x, Double_t *par)
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

  if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
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


  return alphaCB2VWG(x,par)*par[12] + (1. - alphaCB2VWG(x,par))*FitFunctionBackgroundPol3Cheb(x,&par[13]);
}


//-------------------------------------------------------------------------------
Double_t FitFunctionFlowS2CB2POL2EXPPOL2(Double_t *x, Double_t *par)
{
  Double_t SPsiPFactor = 1.0341;
  return alphaCB2POL2EXP(x,par)*par[12]  + (1. - alphaCB2POL2EXP(x,par))*FitFunctionBackgroundPol2(x,&par[13]);

}

//-------------------------------------------------------------------------------
Double_t alphaCB2POL2EXP(Double_t*x, Double_t* par)
{
  return FitFunctionSignalCrystalBallExtended(x, &par[4])/(FitFunctionSignalCrystalBallExtended(x, &par[4]) + FitFunctionBackgroundPol2Exp(x,par));
}

//____________________________________________________________________________
Double_t FitFunctionBackgroundPol2Exp(Double_t *x, Double_t *par)
{
  // pol2 x exp : 4 params
  if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
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


