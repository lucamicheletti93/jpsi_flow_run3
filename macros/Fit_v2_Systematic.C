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

TF1* v2Fit(int iPt,TProfile *v2prof, int type,Double_t parLimit[3]);

Double_t xmin= 2.5;
Double_t xmax= 4.8;

const Double_t v2xmin= 2.5;
const Double_t v2xmax= 4.8;

vector<double> Fitpar;
TFile *fOut;
TFile *fv2fit;
Double_t Jpsiv2 ;
Double_t Jpsiv2_err ;
Double_t v2_ch_ndf;

vector<double> v2par;
vector<double> v2par_er;
std::vector<TH1*> histsys;

std::string Fit_subStr1 = "CB2_VWG";
std::string Fit_subStr2 = "CB2_Pol4Exp";
std::string Fit_subStr3 = "NA60_Pol4Exp";
std::string Fit_subStr4 = "NA60_VWG";


const int nFits_fun = 6;
string mass_sig[] = {"CB2","CB2","CB2","CB2","NA60","NA60"};
string mass_bkg[] = {"VWG_data","VWG_MC","Pol4Exp_data","Pol4Exp_MC","VWG_MC","Pol4Exp_MC"};

string func_sig[] = {"CB2","CB2","CB2","CB2","NA60","NA60"};
string func_bkg[] = {"VWG","VWG","Pol4Exp","Pol4Exp","VWG","Pol4Exp"};


/*
  const int nFits_fun = 1;
  string mass_sig[] = {"NA60"};
  string mass_bkg[] = {"Pol4Exp_data"};

  string func_sig[] = {"NA60"};
  string func_bkg[] = {"Pol4Exp"};
*/

//CB2+VWG, 2.3-4.7, datatail, fixed-mixing
////CB2+VWG, 2.3-4.7, same combination
double minCentrBins[] = {10};
double maxCentrBins[] = {50};


const int nPtBins = 10;
double minPtBins[] = {0.0,1.0,2,3,4,5,6,8,10,12};
double maxPtBins[] = {1.0,2.0,3,4,5,6,8,10,12,15};

//const int nPtBins = 1;
//double minPtBins[] = {12};
//double maxPtBins[] = {15};

const int nRanges = 3;
double min_range[] = {2.3,2.4,2.5};
double max_range[] = {4.7,4.6,4.5};

const int nbck_type =2;
int bck_type[2] = {0,1};
string func_v2bck[] = {"Pol2", "Cheb"};

int nSysbins = nFits_fun*nRanges*nbck_type;

double resolution = 1.3093;

double norm_sig[] = {2.12944e+04,3.61983e+04,2.49511e+04,1.24006e+04,6.77359e+03,3.69119e+03,3.02467e+03,9.84673e+02,2.90510e+02,100,10};
double norm_bkg[] = {8.12173e+06,1.13264e+07,1.48518e+06,2.21428e+05,4.73500e+04,1.51403e+04,8.22036e+03,2.02200e+03,5.60000e+02,100,10};

TProfile *prof_v2EP[20];
TProfile *prof_MEv2EP[20];
TH1D *hT_mass[20];
TH1D *hT_mass2[20];
TCanvas *m_plot[20];
TF1 *fitv2[20][20][20];

std::vector<TH1*> Uncertainities(std::vector<TH1*> histlist, vector<double> parameter, vector<double> parameter_er, double value[nPtBins][3]);

void Fit_v2_Systematic()
{
  TCanvas *c1 = new TCanvas("c1", "c1",75,70,1309,696);
  double par[19];

  string sig_plus_bkg[] = {"bkg","aa","bb","cc","sig_Jpsi","mean_Jpsi","width_Jpsi","a_Jpsi","b_Jpsi","c_Jpsi","d_Jpsi"};
                           //("kVWG","mVWG","sVWG1","sVWG2","kPsi","mPsi","sPsi","alPsi","nlPsi","auPsi","nuPsi");

  string cb2sig_plus_pol4expbkg[] = {"bkg","aa","bb","cc","dd","ee","ff","sig_Jpsi","mean_Jpsi","width_Jpsi","a_Jpsi","b_Jpsi","c_Jpsi","d_Jpsi"};
                           //("kVWG","mVWG","sVWG1","sVWG2","kPsi","mPsi","sPsi","alPsi","nlPsi","auPsi","nuPsi");

  string signa60_plus_bkg[] = {"bkg","aa","bb","cc","sig_Jpsi","mean_Jpsi","width_Jpsi","b_Jpsi","c_Jpsi","d_Jpsi","f_Jpsi","g_Jpsi","h_Jpsi","a_Jpsi","e_Jpsi"};
                         //("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi", "p1LJPsi","p2LJPsi","p3LJPsi","p1RJPsi","p2RJPsi","p3RJPsi","aLJPsi","aRJPsi");

  string signa60_plus_pol4expbkg[] = {"bkg","aa","bb","cc","dd","ee","ff","sig_Jpsi","mean_Jpsi","width_Jpsi","b_Jpsi","c_Jpsi","d_Jpsi","f_Jpsi","g_Jpsi","h_Jpsi","a_Jpsi","e_Jpsi"};
                         //("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi", "p1LJPsi","p2LJPsi","p3LJPsi","p1RJPsi","p2RJPsi","p3RJPsi","aLJPsi","aRJPsi");
  

  Double_t mean_v2[nPtBins];
  Double_t mean_v2er[nPtBins];

  TH1F *hSys[nPtBins]; 
  TCanvas *cPt[nPtBins];
  //cPt[nPtBins] = new TCanvas("cPt", "cPt",75,70,1309,696);


  //v2 Fitting Limits ======
  Double_t parLimit0[3];
  parLimit0[0] = 0.0040; 
  parLimit0[1] = 0.0034;
  parLimit0[2] = 0.0041;
  
  Double_t parLimit1[3];
  parLimit1[0] = 0.036405; 
  parLimit1[1] = 0.034;
  parLimit1[2] = 0.038;

  Double_t parLimit2[3];
  parLimit2[0] = 0.077; 
  parLimit2[1] = 0.072;
  parLimit2[2] = 0.087;

  Double_t parLimit3[3];
  parLimit3[0] = 0.091; 
  parLimit3[1] = 0.09014;
  parLimit3[2] = 0.096;

  Double_t parLimit4[3];
  parLimit4[0] = 0.096; 
  parLimit4[1] = 0.094002;
  parLimit4[2] = 0.098;

   Double_t parLimit5[3];
  parLimit5[0] = 0.13; 
  parLimit5[1] = 0.12;
  parLimit5[2] = 0.130103;

  Double_t parLimit6[3];
  parLimit6[0] = 0.079; 
  parLimit6[1] = 0.072623;
  parLimit6[2] = 0.082;

  Double_t parLimit7[3];
  parLimit7[0] = 0.09; 
  parLimit7[1] = 0.088;
  parLimit7[2] = 0.092;

  Double_t parLimit8[3];
  parLimit8[0] = 0.066; 
  parLimit8[1] = 0.058;
  parLimit8[2] = 0.072;

  Double_t parLimit9[3];
  parLimit9[0] = 0.036; 
  parLimit9[1] = 0.028;
  parLimit9[2] = 0.072;
  
  
  TFile *f_hist = new TFile("Histograms_Fullpass3matchedMchMid_centr_Add_Pt_bin_10_50.root");
  if (f_hist->IsOpen()) cout << "File opened successfully" << endl;
 
 for(int iPt = 0;iPt<nPtBins;iPt++) //Loop over mass and v2 histograms in different pT ranges
    {
  hSys[iPt] = new TH1F("hSys", " ",nSysbins,0,nSysbins);

  int ncombination = 0;
  for(int type = 0; type < 2; type++)  //Loop over background function type
    {
      
  for(int kfun = 0; kfun < nFits_fun; kfun++)  //Loop over Signal and backgrounds shape
    {
      TFile *f_fit = new TFile(Form("mass_fit_results/Pt_%0.1f_%0.1f/multi_trial_%s_%s_tails.root",minPtBins[iPt], maxPtBins[iPt],mass_sig[kfun].data(),mass_bkg[kfun].data()));
      if (f_fit->IsOpen()) cout << "File opened successfully" << endl;
      f_fit->ls();

      for(int iR = 0; iR < nRanges ; iR++) //Loop over Signal Range
	{
	   cout<<" Strig name "<<func_sig[iPt].data()<<endl;
	   hT_mass2[iPt] = (TProfile*) f_hist->Get(Form("histMassSEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
	   prof_v2EP[iPt] = (TProfile*) f_hist->Get(Form("histV2SEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
	   prof_MEv2EP[iPt] = (TProfile*) f_hist->Get(Form("histV2MEPM_%1.0f_%1.0f__%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
	   hT_mass[iPt] = (TH1D*)f_fit->Get(Form("fit_results_%s_%s__%0.1f_%0.1f_histMassSEPM_%1.0f_%1.0f__%1.0f_%1.0f",func_sig[kfun].data(),func_bkg[kfun].data(), min_range[iR],max_range[iR],minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));

	   m_plot[iPt] = (TCanvas*)f_fit->Get(Form("fit_plot_%s_%s__%0.1f_%0.1f_histMassSEPM_%1.0f_%1.0f__%1.0f_%1.0f",func_sig[kfun].data(),func_bkg[kfun].data(), min_range[iR],max_range[iR],minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));

	   //hT_mass[iPt]->Draw();
	  
	   std::string Fit_massStr = hT_mass[iPt]->GetTitle();
	   cout<<" Mass Fit Histogram ========================= "<<Fit_massStr.data()<<endl;
	   size_t found_cb2_vwg = Fit_massStr.find(Fit_subStr1);
	   size_t found_cb2_pol4exp = Fit_massStr.find(Fit_subStr2);
	   size_t found_na60_pol4exp = Fit_massStr.find(Fit_subStr3);
	   size_t found_na60_vwg = Fit_massStr.find(Fit_subStr4);

	   
	   if(found_cb2_vwg != std::string::npos)
	     {
	       for (int i = 0; i < 11; i++) 
	       Fitpar.push_back(hT_mass[iPt]->GetBinContent(hT_mass[iPt]->GetXaxis()->FindBin(Form("%s",sig_plus_bkg[i].data()))));
	       Fitpar.push_back(0.0);
	     }

	   if(found_cb2_pol4exp != std::string::npos)
	     {
	       for (int i = 0; i < 14; i++)
	       Fitpar.push_back(hT_mass[iPt]->GetBinContent(hT_mass[iPt]->GetXaxis()->FindBin(Form("%s",cb2sig_plus_pol4expbkg[i].data()))));
	       Fitpar.push_back(0.0);
	       //parLimit5[3] = {0.12,0.1,0.14};
	     }

	   if(found_na60_vwg != std::string::npos)
	     {
	       for (int i = 0; i < 15; i++)
	       Fitpar.push_back(hT_mass[iPt]->GetBinContent(hT_mass[iPt]->GetXaxis()->FindBin(Form("%s",signa60_plus_bkg[i].data()))));
	       Fitpar.push_back(0.0);
	     }

	   if(found_na60_pol4exp != std::string::npos)
	     {
	       for (int i = 0; i < 18; i++)
	       Fitpar.push_back(hT_mass[iPt]->GetBinContent(hT_mass[iPt]->GetXaxis()->FindBin(Form("%s",signa60_plus_pol4expbkg[i].data()))));
	       Fitpar.push_back(0.0);
	     }
	   
	   //hT_mass[iPt]->Draw();

	   prof_v2EP[iPt]->SetName(Form("fitv2_Sig_Bkg:%s_%s__FitRange:%0.1f_%0.1f__Pt:%1.0f_%1.0f__CENT:%1.0f_%1.0f",mass_sig[kfun].data(),mass_bkg[kfun].data(), min_range[iR],max_range[iR],minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
	   prof_v2EP[iPt]->SetTitle(Form("fitv2_Sig_Bkg:%s_%s__Pt:%1.0f_%1.0f",mass_sig[kfun].data(),mass_bkg[kfun].data(),minPtBins[iPt], maxPtBins[iPt]));
	  
	   if(iPt ==0) fitv2[iPt][kfun][iR] = v2Fit(iPt,prof_v2EP[iPt], type,parLimit0);
	  else if(iPt ==1) fitv2[iPt][kfun][iR] = v2Fit(iPt, prof_v2EP[iPt], type,parLimit1);
	  else if(iPt ==2) fitv2[iPt][kfun][iR] = v2Fit(iPt,prof_v2EP[iPt], type,parLimit2);
	  else if(iPt ==3) fitv2[iPt][kfun][iR] = v2Fit(iPt,prof_v2EP[iPt], type,parLimit3);
	  else if(iPt ==4) fitv2[iPt][kfun][iR] = v2Fit(iPt,prof_v2EP[iPt], type,parLimit4);
	  else if(iPt ==5) fitv2[iPt][kfun][iR] = v2Fit(iPt,prof_v2EP[iPt], type,parLimit5);
	  else if(iPt ==6) fitv2[iPt][kfun][iR] = v2Fit(iPt,prof_v2EP[iPt], type,parLimit6);
	  else if(iPt ==7) fitv2[iPt][kfun][iR] = v2Fit(iPt,prof_v2EP[iPt], type,parLimit7);
	  else if(iPt ==8) fitv2[iPt][kfun][iR] = v2Fit(iPt,prof_v2EP[iPt], type,parLimit8);
	  else if(iPt ==9) fitv2[iPt][kfun][iR] = v2Fit(iPt,prof_v2EP[iPt], type,parLimit9);

	   prof_v2EP[iPt]->Draw("p");
           //prof_MEv2EP[iPt]->Draw("psame");
	   
	   fitv2[iPt][kfun][iR]->SetName(Form("fitv2_Sig_Bkg:%s_%s__FitRange:%0.1f_%0.1f__Pt:%1.0f_%1.0f__CENT:%1.0f_%1.0f",mass_sig[kfun].data(),mass_bkg[kfun].data(), min_range[iR],max_range[iR],minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
	   fitv2[iPt][kfun][iR]->SetTitle(Form("fitv2_Sig_Bkg:%s_%s__FitRange:%0.1f_%0.1f",mass_sig[kfun].data(),mass_bkg[kfun].data(), min_range[iR],max_range[iR]));
           
	  hSys[iPt]->GetXaxis()->SetBinLabel(ncombination+1,Form("fitv2_Sig_Bkg:%s_%s__FitRange:%0.1f_%0.1f_%s",mass_sig[kfun].data(),mass_bkg[kfun].data(), min_range[iR],max_range[iR],func_v2bck[type].data()));
	  
	   cout<<" v2 ========================== /////// "<<Jpsiv2<<" +/- "<<Jpsiv2_err<<endl;
	   v2par.push_back(Jpsiv2);
	   v2par_er.push_back(Jpsiv2_err);

	   ncombination++;
	   cout<<" ncombination "<<ncombination<<endl;
	   
	   
	  }

	}

    }  //background function type
    
  hSys[iPt]->SetTitle(Form("SystV2SEPM_Pt:%1.0f_%1.0f__CENT:%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
  hSys[iPt]->SetName(Form("SystV2SEPM_Pt:%1.0f_%1.0f__CENT:%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
  hSys[iPt]->SetMarkerColor(kRed);
  hSys[iPt]->SetMarkerStyle(20);
  hSys[iPt]->SetMarkerSize(0.9);

  histsys.push_back(hSys[iPt]);
  //hSys[iPt]->Draw();
    }


 //=======Systematic uncertainities ===============
 double value[nPtBins][3];
 std::vector<TH1*> myHist = Uncertainities(histsys,v2par, v2par_er, value);

 for(int iPt = 0;iPt<nPtBins;iPt++) //Loop over v2_signal and background
   {
 double fmean_v2 = value[iPt][0];
 double fmean_v2er = value[iPt][1];
 double fsys = value[iPt][2];
 cout<<iPt<<" iPt "<< fmean_v2<<endl;
 Double_t fper_stat = (fmean_v2er/fmean_v2)*100;
 Double_t fper_sys = (fsys/fmean_v2)*100;
 //for(int j=0; j <18; j++) cout<<" v2 parameter outside "<<myHist[iPt]->GetBinContent(j)<<endl;
 cout<<" Average Again ==== "<<fmean_v2<< " $\pm$ " <<fmean_v2er<<" ("<<fper_stat<<"\%)"<<" $\pm$ "<<fsys<<" ("<<fper_sys<<"\%)"<<endl;
      
 cPt[iPt] = new TCanvas("c1", "c1",75,70,1309,696);
 cPt[iPt]->SetName(Form("SystV2SEPM_Pt:%1.0f_%1.0f__CENT:%1.0f_%1.0f", minPtBins[iPt], maxPtBins[iPt], minCentrBins[0], maxCentrBins[0]));
 cPt[iPt]->cd();

 myHist[iPt]->SetMarkerStyle(20);
 myHist[iPt]->SetMarkerColor(kRed);
 myHist[iPt]->SetMarkerSize(1.2);
 myHist[iPt]->Draw("p");

 TLine *line = new TLine(0,fmean_v2,nSysbins,fmean_v2);
 line->SetLineStyle(1);
 line->SetLineColor(kBlue);
 line->SetLineWidth(2);
 line->Draw();

 TLine *line1 = new TLine(0,fsys+fmean_v2,nSysbins,fsys+fmean_v2);
 line1->SetLineStyle(2);
 line1->SetLineColor(kBlue);
 line1->SetLineWidth(2);
 line1->Draw();

 TLine *line2 = new TLine(0,fmean_v2-fsys,nSysbins,fmean_v2-fsys);
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

 TText *text1 = display1->AddText(Form(" v_{2} {J/#psi} =%.3f +/- %.3f (%.2f %%) +/- %.3f (%.2f %% ) ",fmean_v2,fmean_v2er,fper_stat,fsys,fper_sys));
 display1->Draw("same");

   }

 
 //=======Plotting ==================
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(1);
 fOut = new TFile("Systematic_Cent10-50.root","RECREATE");
 for(int i =0; i<nPtBins; i++)
 cPt[i]->Write();
 fOut->Close();

 fv2fit = new TFile("FitResults_v2-50.root","RECREATE");
 for(int i =0; i<nPtBins; i++)
 prof_v2EP[i]->Write();
 fv2fit->Close();

 
TCanvas *c2 = new TCanvas("c2", "c2",75,70,1309,696);
 c2->Divide(5,2);
 for(int i =0; i<nPtBins; i++)
   {
     c2->cd(i+1);
     prof_v2EP[i]->SetMarkerColor(kBlack);
     prof_v2EP[i]->SetMarkerStyle(20);
     prof_v2EP[i]->SetMarkerSize(0.7);
     prof_v2EP[i]->GetYaxis()->SetTitle("v_{2}^{#mu#mu}");
     prof_v2EP[i]->SetTitle(Form("Pt: %1.0f - %1.0f GeV/c ",minPtBins[i], maxPtBins[i]));

     prof_MEv2EP[i]->SetMarkerColor(kRed);
     prof_MEv2EP[i]->SetMarkerStyle(20);
     prof_MEv2EP[i]->SetMarkerSize(0.7);
     //prof_MEv2EP[i]->Rebin(2);
     prof_MEv2EP[i]->Scale(2.8);
     prof_v2EP[i]->Draw("e");
     //prof_MEv2EP[i]->Draw("psame");
     //m_plot[i]->Draw();
}
 

 
}


TF1* v2Fit(int iPt,TProfile *v2prof, int type, Double_t parLimit[3])
{
  cout << "Output of para begin and end: "<<endl;
  for (int i = 0; i < Fitpar.size(); i++)
   cout << Fitpar[i] << " ";
  cout <<endl;
  v2prof->Rebin(2);
  TCanvas *c3 = new TCanvas("c3", "c3",75,70,1309,696);
  v2prof->Scale(resolution);
  v2prof->GetXaxis()->SetRangeUser(xmin,xmax);
  v2prof->Draw();

  //=============================== v2 Fitting =================
  TF1* bck = new TF1("bck",FitFunctionBackgroundPol2,xmin,xmax,3);
  bck->FixParameter(3,type);
  bck->SetLineColor(kBlue);
  bck->SetLineStyle(9);
  //bck->SetParameters(0.0191306,-0.00836616,0.00102991);
  bck->SetParameters(0.106962,-0.0252894,0.0022179);
  //SetFitRejectRange(2.9,3.3);
  v2prof->Fit(bck,"SERI","",v2xmin,v2xmax);
  //bck->Draw("same");
  //SetFitRejectRange();

 
  std::string Fit_v2Str = v2prof->GetName();
  cout<<" Fit Function ========================= "<<Fit_v2Str.data()<<endl;
  size_t found_cb2_vwg = Fit_v2Str.find(Fit_subStr1);
  size_t found_cb2_pol4exp = Fit_v2Str.find(Fit_subStr2);
  size_t found_na60_pol4exp = Fit_v2Str.find(Fit_subStr3);
  size_t found_na60_vwg = Fit_v2Str.find(Fit_subStr4);

 TF1* fitv2;

 if(found_cb2_vwg != std::string::npos) 
   {

  //Fit mass spectra and get normalization parameters ..........
  TF1 *fitFctCB2VWG = new TF1("fitFctCB2VWG", fitFunctionCB2VWG,v2xmin,v2xmax,11);
  fitFctCB2VWG->SetParameter(0,norm_bkg[iPt]);
  for (int i = 1; i < 4; i++) fitFctCB2VWG->FixParameter(i,Fitpar[i]);
  fitFctCB2VWG->SetParameter(4,norm_sig[iPt]);
  for (int i = 5; i < 11; i++) fitFctCB2VWG->FixParameter(i,Fitpar[i]);
  hT_mass2[iPt]->Fit(fitFctCB2VWG,"REMSI");
  //for (int i = 0; i < Fitpar.size(); i++) cout<<" PARAMETERS COMPARISION \\\\\\\\\\  "<<fitFctCB2VWG->GetParameter(i)<<"   "<<Fitpar[i]<<endl;
  
  fitv2 = new TF1("fitv2",FitFunctionFlowS2CB2VWGPOL2,v2xmin,v2xmax,17);
  fitv2->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
  fitv2->SetLineColor(kRed);
  fitv2->SetParName(11,"kPsiP");
  fitv2->SetParName(12,"v_{2} JPsi");
  fitv2->SetParName(13,"v_{2} BG0");
  fitv2->SetParName(14,"v_{2} BG1");
  fitv2->SetParName(15,"v_{2} BG2");
  fitv2->FixParameter(16,type);
  
  for (int i = 0; i < Fitpar.size(); i++) fitv2->FixParameter(i,Fitpar[i]);
  //for (int i = 0; i < Fitpar.size(); i++) cout<<" ================== \\\\\\\\\\  "<<fitv2->GetParameter(i)<<endl;
  //for( Int_t i = 0; i < 12; ++i ) fitv2->FixParameter(i,par[i]);
  fitv2->SetParameter(12, parLimit[0]);
  fitv2->SetParLimits(12, parLimit[1],parLimit[2]);
  fitv2->FixParameter(0,fitFctCB2VWG->GetParameter(0));
  fitv2->FixParameter(4,fitFctCB2VWG->GetParameter(4));
  
  for ( Int_t i = 0; i < 3; ++i ) { fitv2->SetParameter(i + 13, bck->GetParameter(i));}
 
  TFitResultPtr r2 = v2prof->Fit(fitv2,"REI");

  TF1* bck2 = new TF1("bck2",FitFunctionBackgroundPol2,v2xmin,v2xmax,4);
  bck2->FixParameter(3,type);
  for ( Int_t i = 0; i < 3; ++i ) { bck2->FixParameter(i, fitv2->GetParameter(i + 13));}
  bck2->SetLineColor(kBlue);
  bck2->SetLineStyle(2);
  v2prof->Fit(bck2,"RE+");
  bck2->Draw("same");

  Jpsiv2 =  fitv2->GetParameter(12);
  Jpsiv2_err =  fitv2->GetParError(12);
  
   }
 
 if (found_cb2_pol4exp != std::string::npos) 
    {
      cout<<"CB2_Pol4Exp EXIST"<<endl;
      //Fit mass spectra and get normalization parameters ..........
      TF1 *fitFctCB2Pol4Exp = new TF1("fitFctCB2Pol4Exp",fitFunctionCB2Pol4Exp,v2xmin,v2xmax,14);
      fitFctCB2Pol4Exp->SetParameter(0,norm_bkg[iPt]);
      for (int i = 1; i < 7; i++) fitFctCB2Pol4Exp->FixParameter(i,Fitpar[i]);
      fitFctCB2Pol4Exp->SetParameter(7,norm_sig[iPt]);
      for (int i = 8; i < 14; i++) fitFctCB2Pol4Exp->FixParameter(i,Fitpar[i]);
      hT_mass2[iPt]->Fit(fitFctCB2Pol4Exp,"REMSI");
      //for (int i = 0; i < Fitpar.size(); i++) cout<<" PARAMETERS COMPARISION \\\\\\\\\\  "<<fitFctCB2Pol4Exp->GetParameter(i)<<"   "<<Fitpar[i]<<endl;

      fitv2 = new TF1("fitv2",FitFunctionFlowS2CB2POL4EXPPOL2,v2xmin,v2xmax,20);
      fitv2->SetParNames("pol0","pol1","pol2","pol3","exp1","exp2","exp3","kJPsi","mJPsi","sJPsi","alJPsi");
      fitv2->SetParName(11,"nlJPsi");
      fitv2->SetParName(12,"auJPsi");
      fitv2->SetParName(13,"nuJPsi");
      fitv2->SetParName(14,"kPsiP");
      fitv2->SetParName(15,"v_{2} JPsi");
      fitv2->SetParName(16,"v_{2} BG0");
      fitv2->SetParName(17,"v_{2} BG1");
      fitv2->SetParName(18,"v_{2} BG2");
      fitv2->FixParameter(19,type);
      
      for (int i = 1; i < Fitpar.size(); i++) fitv2->FixParameter(i,Fitpar[i]);
      for (int i = 0; i < Fitpar.size(); i++) cout<<" outFitparCB2 "<<Fitpar[i]<<", "<<endl;
      //fitv2->FixParameter(0,0.06);
      //for( Int_t i = 0; i < 15; ++i ) fitv2->FixParameter(i,par[i]);
      fitv2->SetParameter(15,parLimit[0]);
      fitv2->SetParLimits(15,parLimit[1],parLimit[2]);
      fitv2->FixParameter(0,fitFctCB2Pol4Exp->GetParameter(0));
      fitv2->FixParameter(7,fitFctCB2Pol4Exp->GetParameter(7));
  
     for ( Int_t i = 0; i < 3; ++i ) { fitv2->SetParameter(i + 16, bck->GetParameter(i));}

     TFitResultPtr r2 = v2prof->Fit(fitv2,"REI");

     TF1* bck2 = new TF1("bck2",FitFunctionBackgroundPol2,v2xmin,v2xmax,4);
     for ( Int_t i = 0; i < 3; ++i ) { bck2->FixParameter(i, fitv2->GetParameter(i + 16));}
     bck2->SetLineColor(kBlue);
     bck2->SetLineStyle(2);
     v2prof->Fit(bck2,"RE+");
     //bck2->Draw("same");

     Jpsiv2 =  fitv2->GetParameter(15);
     Jpsiv2_err =  fitv2->GetParError(15);

     cout<<" EXP jpsi v2 "<<Jpsiv2<<" "<<Jpsiv2_err<<endl;
    }

 if (found_na60_vwg != std::string::npos) 
   {
     cout<<"NA60_VWG EXIST"<<endl;
     //Fit mass spectra and get normalization parameters NA60..........
     TF1 *fitFctna60vwg = new TF1("fitFctna60vwg",fitFunctionNA60NEWVWG,v2xmin,v2xmax,15);
     fitFctna60vwg->SetParameter(0,norm_bkg[iPt]);
     for (int i = 1; i < 4; i++) fitFctna60vwg->FixParameter(i,Fitpar[i]);
     fitFctna60vwg->SetParameter(4,norm_sig[iPt]);
     for (int i = 5; i < 15; i++) fitFctna60vwg->FixParameter(i,Fitpar[i]);
     hT_mass2[iPt]->Fit(fitFctna60vwg,"REMSI");
     //for (int i = 0; i < Fitpar.size(); i++) cout<<" PARAMETERS COMPARISION NA60\\\\\\\\\\  "<<fitFctna60vwg->GetParameter(i)<<"   "<<Fitpar[i]<<endl;

       
     fitv2 = new TF1("fitv2",FitFunctionFlowS2NA60NEWVWGPOL2,v2xmin,v2xmax,21);
     fitv2->SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","p1LJPsi","p2LJPsi","p3LJPsi","p1RJPsi");
     fitv2->SetParName(11,"p2RJPsi");
     fitv2->SetParName(12,"p3RJPsi");
     fitv2->SetParName(13,"aLJPsi");
     fitv2->SetParName(14,"aRJPsi");
     fitv2->SetParName(15,"kPsiP");
     fitv2->SetParName(16,"v_{2} JPsi");
     fitv2->SetParName(17,"v_{2} BG0");
     fitv2->SetParName(18,"v_{2} BG1");
     fitv2->SetParName(19,"v_{2} BG2");
     fitv2->FixParameter(20,type);

     for (int i = 0; i < Fitpar.size(); i++) fitv2->FixParameter(i,Fitpar[i]);
     //for( Int_t i = 0; i < 16; ++i ) fitv2->FixParameter(i,par[i]);
     fitv2->SetParameter(16,parLimit[0]);
     fitv2->SetParLimits(16,parLimit[1],parLimit[2]);
     fitv2->FixParameter(0,fitFctna60vwg->GetParameter(0));
     fitv2->FixParameter(4,fitFctna60vwg->GetParameter(4));
     
     for ( Int_t i = 0; i < 3; ++i ) { fitv2->SetParameter(i + 17, bck->GetParameter(i));}

     TFitResultPtr r2 = v2prof->Fit(fitv2,"REI");

     TF1* bck2 = new TF1("bck2",FitFunctionBackgroundPol2,v2xmin,v2xmax,4);
     for ( Int_t i = 0; i < 3; ++i ) { bck2->FixParameter(i, fitv2->GetParameter(i + 17));}
     bck2->SetLineColor(kBlue);
     bck2->SetLineStyle(2);
     //v2prof->Fit(bck2,"RE+");
     //bck2->Draw("same");

     Jpsiv2 =  fitv2->GetParameter(16);
     Jpsiv2_err =  fitv2->GetParError(16);
   }

 if (found_na60_pol4exp != std::string::npos) 
   {
     cout<<"NA60_EXPPol EXIST"<<endl;
     //Fit mass spectra and get normalization parameters NA60..........
     TF1 *fitFctna60pol4exp = new TF1("fitFctna60pol4exp",fitFunctionNA60NEWPol4Exp,v2xmin,v2xmax,18);
     fitFctna60pol4exp->SetParameter(0,norm_bkg[iPt]);
     for (int i = 1; i < 7; i++) fitFctna60pol4exp->FixParameter(i,Fitpar[i]);
     fitFctna60pol4exp->SetParameter(7,norm_sig[iPt]);
     for (int i = 8; i < 18; i++) fitFctna60pol4exp->FixParameter(i,Fitpar[i]);
     hT_mass2[iPt]->Fit(fitFctna60pol4exp,"REMSI");
     //for (int i = 0; i < Fitpar.size(); i++) cout<<" PARAMETERS COMPARISION NA60PLO4EXP\\\\\\\\\\  "<<fitFctna60pol4exp->GetParameter(i)<<"   "<<Fitpar[i]<<endl;

     
     fitv2 = new TF1("fitv2",FitFunctionFlowS2NA60NEWPOL4EXPPOL2,v2xmin,v2xmax,24);
     fitv2->SetParNames("pol0","pol1","pol2","pol3","exp1","exp2","exp3","kJPsi","mJPsi","sJPsi","p1LJPsi");
     fitv2->SetParName(11,"p2LJPsi");
     fitv2->SetParName(12,"p3LJPsi");
     fitv2->SetParName(13,"p1RJPsi");   
     fitv2->SetParName(14,"p2RJPsi");
     fitv2->SetParName(15,"p3RJPsi");
     fitv2->SetParName(16,"aLJPsi");
     fitv2->SetParName(17,"aRJPsi");
     fitv2->SetParName(18,"kPsiP");
     fitv2->SetParName(19,"v_{2} JPsi");
     fitv2->SetParName(20,"v_{2} BG0");
     fitv2->SetParName(21,"v_{2} BG1");
     fitv2->SetParName(22,"v_{2} BG2");
     fitv2->FixParameter(23,type);

     //for( Int_t i = 0; i < 19; ++i ) fitv2->FixParameter(i,par[i]);
     for (int i = 0; i < Fitpar.size(); i++) fitv2->FixParameter(i,Fitpar[i]);
     //fitv2->FixParameter(0,0.06);
     fitv2->SetParameter(19, parLimit[0]);
     fitv2->SetParLimits(19, parLimit[1],parLimit[2]);
     fitv2->FixParameter(0,fitFctna60pol4exp->GetParameter(0));
     fitv2->FixParameter(7,fitFctna60pol4exp->GetParameter(7));
  
     for ( Int_t i = 0; i < 3; ++i ) { fitv2->SetParameter(i + 20, bck->GetParameter(i));}

     TFitResultPtr r2 = v2prof->Fit(fitv2,"REI");

     TF1* bck2 = new TF1("bck2",FitFunctionBackgroundPol2,v2xmin,v2xmax,4);
     for ( Int_t i = 0; i < 3; ++i ) { bck2->FixParameter(i, fitv2->GetParameter(i + 20));}
     bck2->SetLineColor(kBlue);
     bck2->SetLineStyle(2);
     //v2prof->Fit(bck2,"RE+");
     //bck2->Draw("same");

     Jpsiv2 =  fitv2->GetParameter(19);
     Jpsiv2_err =  fitv2->GetParError(19);
   }

 
  bck->SetName("Fit_mainStr.data()");
 
  v2_ch_ndf = fitv2->GetChisquare()/fitv2->GetNDF();

  Fitpar.clear();
  return fitv2;
}


std::vector<TH1*> Uncertainities(std::vector<TH1*> histlist, vector<double> parameter, vector<double> parameter_er, double value[nPtBins][3])
{
  std::vector<TH1*> histograms;
  int ncombination = parameter.size()/nPtBins;

  int nbin =0;
  for (int i = 0; i < parameter.size(); i++)
    {
      if ((i) % ncombination != 0) continue;
      if ((i) % ncombination == 0) cout<<" HERE IS THE DIVISIONS "<<i<<endl;
      double sum_v2 =0; double sum_v2er =0;
      int ibin =0;
      for (int j = i; j < i+ncombination; j++)
	{
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
      for (int j = i; j < i+ncombination; j++)
	{
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

  if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  //if (x[0] > 3.5 && x[0] < 3.7) TF1::RejectPoint();
  if(type==0) return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];  //return pol2 function

  //Return chebyshev Pol2 function
  else return FitFunctionBackgroundPol2Cheb(x,par); 
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


