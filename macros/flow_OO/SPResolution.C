#include "TAxis.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include <TCanvas.h>
#include <TChain.h>
#include <TClass.h>
#include <TCollection.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphMultiErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <THashList.h>
#include <THnSparse.h>
#include <TKey.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TMatrixD.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TParameter.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "FlowAnalysis_Helper.h"
//#include "Framework/Logger.h"

using namespace std;
vector<double> GetMeanError(int size, double *sample);
void LoadStyle();
void SetLegend(TLegend *);

void SPResolution(
    string Label, string FileName = "input_AnalysisResults.txt") 
{

   //gSystem->CompileMacro("FlowAnalysis_Helper.cxx");
   gSystem->Load("FlowAnalysis_Helper_cxx.so");    

  FlowAnalysis_Helper helper;
  fstream InputFiles;
  vector<vector<double>> stats;
  THStack *hs_stack = new THStack("Run_by_Run_ResolutionFactors", "");
  InputFiles.open(FileName, ios::in);
  TList *ls = new TList();
  TList *ls_run_sp = new TList();
  TList *ls_run_ep = new TList();
  int Nrow = 0;
  int Ncol = 0;
  int NBins_cent = 0;
  int NBins_cent_profile = 0;
  double *Bin_cent = nullptr;
  double *Bin_cent_profile = nullptr;
  TFile *f = new TFile("Sp_resolution_dummy.root","RECREATE");
  if (InputFiles.is_open()) {
    string File;
    cout << "Start reading input AnalysisResults.root ..." << endl;
    //TH1D *hist_r2sp[10];
    while (getline(InputFiles, File)) {
      cout << "Reading input from: " << File << endl;
      TH2F *hs_R2SPAB, *hs_R2SPAC, *hs_R2SPBC;
      TProfile *tp_R2SPAB, *tp_R2SPAC, *tp_R2SPBC;
      TProfile *tp_R2SPAB_Im, *tp_R2SPAC_Im, *tp_R2SPBC_Im;
      TProfile *tp_R2EPAB, *tp_R2EPAC, *tp_R2EPBC;
      TProfile *tp_R2EPAB_Im, *tp_R2EPAC_Im, *tp_R2EPBC_Im;
      vector<string> File_string = helper.tokenize(File);
      string run_number = File_string[File_string.size() - 2];
      helper.LoadReso(File.c_str(), hs_R2SPAB, hs_R2SPAC, hs_R2SPBC);

      // From histograms
      Bin_cent = helper.CreateBinsFromAxis(hs_R2SPAB->GetXaxis());
      NBins_cent = hs_R2SPAB->GetXaxis()->GetNbins();
      Ncol = NBins_cent;
      TH1D *hist_r2sp =
          new TH1D(Form("R2SP_Cent_%s", run_number.c_str()),
                   Form("R_{2}{SP} for run: %s", run_number.c_str()),
                   NBins_cent, Bin_cent);
      hist_r2sp->GetXaxis()->SetTitle("Centrality FT0C(%)");
      hist_r2sp->GetYaxis()->SetTitle("R_{2}{SP}");

      vector<double> row;
      for (int i = 0; i < NBins_cent; i++) {
        // Resolution factor for SP
        TH2F *hs_R2SPAB_cp = dynamic_cast<TH2F *>(hs_R2SPAB->Clone(
            Form("R2SPAB_Cent_Copy_%g_%g", Bin_cent[i], Bin_cent[i + 1])));
        TH2F *hs_R2SPAC_cp = dynamic_cast<TH2F *>(hs_R2SPAC->Clone(
            Form("R2SPAC_Cent_Copy_%g_%g", Bin_cent[i], Bin_cent[i + 1])));
        TH2F *hs_R2SPBC_cp = dynamic_cast<TH2F *>(hs_R2SPBC->Clone(
            Form("R2SPBC_Cent_Copy_%g_%g", Bin_cent[i], Bin_cent[i + 1])));
        hs_R2SPAB_cp->GetXaxis()->SetRangeUser(Bin_cent[i], Bin_cent[i + 1]);
        hs_R2SPAC_cp->GetXaxis()->SetRangeUser(Bin_cent[i], Bin_cent[i + 1]);
        hs_R2SPBC_cp->GetXaxis()->SetRangeUser(Bin_cent[i], Bin_cent[i + 1]);

        double R2SPAB = hs_R2SPAB_cp->GetMean(2);
        double R2SPAC = hs_R2SPAC_cp->GetMean(2);
        double R2SPBC = hs_R2SPBC_cp->GetMean(2);
        double R2SPABe = hs_R2SPAB_cp->GetMean(12);
        double R2SPACe = hs_R2SPAC_cp->GetMean(12);
        double R2SPBCe = hs_R2SPBC_cp->GetMean(12);

        double R22SP = R2SPBC != 0 ? R2SPAB * R2SPAC / R2SPBC : 0.0;
        double R2SP = R22SP > 0 ? TMath::Sqrt(R22SP) : 0.0;
        double R2SPe =
            R2SPAB * R2SPAC * R2SPBC == 0
                ? 0.0
                : TMath::Sqrt(
                      1. / 4 * (R2SPAC / (R2SPAB * R2SPBC)) * pow(R2SPABe, 2.) +
                      1. / 4 * (R2SPAB / (R2SPAC * R2SPBC)) * pow(R2SPACe, 2.) +
                      1. / 4 * R2SPAC * R2SPAB / pow(R2SPBC, 3.) *
                          pow(R2SPBCe, 2.));
        R2SPe = isnan(R2SPe) || isinf(R2SPe) ? 0. : R2SPe;

        hist_r2sp->SetBinContent(i + 1, R2SP);
        hist_r2sp->SetBinError(i + 1, R2SPe);

        //hist_r2sp->Rebin(2);
        //if(ifile == 0) hist_r2sp[ifile]->Draw();
        // else hist_r2sp[ifile]->Draw("same");

        //if(ifile > 0) hist_r2sp[0]->Add(hist_r2sp[ifile]);
        row.emplace_back(R2SP);
        delete hs_R2SPAB_cp;
        delete hs_R2SPAC_cp;
        delete hs_R2SPBC_cp;
       
      }

      ls->Add(hist_r2sp);
      hs_stack->Add(hist_r2sp);
      stats.emplace_back(row);
      Nrow++;
    }

}
 vector<double *> stats_trans;
  for (int i = 0; i < Ncol; i++) {
    stats_trans.emplace_back(new double[Nrow]);
  }
  for (int i = 0; i < Nrow; i++) {
    for (int j = 0; j < Ncol; j++) {
      stats_trans[j][i] = stats[i][j];
    }
  }

  TH1D *hist_final =
      new TH1D(Form("R2SP_All_%s", Label.c_str()),
               Form("R2SP_All_%s", Label.c_str()), NBins_cent, Bin_cent);
  hist_final->GetXaxis()->SetTitle("Centrality FT0C (%)");
  hist_final->GetYaxis()->SetTitle("R_{2}{SP}");
  for (int i = 0; i < Ncol; i++) {
    vector<double> final_stats = GetMeanError(Nrow, stats_trans[i]);
    hist_final->SetBinContent(i + 1, final_stats[0]);
    hist_final->SetBinError(i + 1, final_stats[1]);
  }
  TList *ls_final = new TList();
  ls_final->Add(hist_final);
  ls_final->Add(hs_stack);

  TFile froot(Form("ResolutionsAll_%s.root", Label.c_str()), "RECREATE");
  ls->Write("RunByRun_Histogram", TObject::kSingleKey);
  ls_run_sp->Write("RunByRun_HistogramFromProfile_SP", TObject::kSingleKey);
  ls_run_ep->Write("RunByRun_HistogramFromProfile_EP", TObject::kSingleKey);
  ls_final->Write("FinalStatistics", TObject::kSingleKey);
  f->cd();
  hist_final->Write();
  f->Close();
  froot.Close();
  InputFiles.close();

  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<double> GetMeanError(int size, double *sample) {
  vector<double> results;
  double mean = 0.;
  for (int i = 0; i < size; i++) {
    mean += sample[i] / size;
  }
  results.emplace_back(mean);

  double sum2 = 0;
  for (int i = 0; i < size; i++) {
    sum2 += pow(sample[i] - mean, 2.) / size;
  }
  double rms = pow(sum2, 0.5);
  results.emplace_back(rms);
  return results;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void plot_resolution() {
  LoadStyle();
  gStyle -> SetOptStat(false);
  TFile *fIn = new TFile("ResolutionsAll_LHC25ae_pass2.root", "READ");
  TList *list = (TList*) fIn -> Get("RunByRun_Histogram");
  TH1D *histResolution = (TH1D*) list -> FindObject("R2SP_Cent_LHC25ae");
  histResolution -> GetYaxis() -> SetRangeUser(0, 0.7);
  histResolution -> GetYaxis() -> SetTitle("#it{R}_{2}^{SP}");
  histResolution -> SetTitle("");
  histResolution -> SetMarkerStyle(20);
  histResolution -> SetMarkerColor(kBlack);
  histResolution -> SetLineColor(kBlack);

  TLine *line020 = new TLine(20, 0, 20, 0.7);
  line020 -> SetLineColor(kGray+2);
  line020 -> SetLineWidth(2);
  line020 -> SetLineStyle(kDashed);

  TLine *line010 = new TLine(10, 0, 10, 0.7);
  line010 -> SetLineColor(kGray+2);
  line010 -> SetLineWidth(2);
  line010 -> SetLineStyle(kDashed);

  TF1 *func020 = new TF1("func020", "pol0", 0, 20);
  func020 -> SetLineColor(kRed+1);
  histResolution -> Fit(func020, "R0Q");

  TF1 *func010 = new TF1("func010", "pol0", 0, 10);
  func010 -> SetLineColor(kOrange+7);
  histResolution -> Fit(func010, "R0Q");

  TCanvas *canvasResolution = new TCanvas("canvasResolution", "", 800, 600);
  histResolution -> Draw("EP");
  func020 -> Draw("SAME");
  func010 -> Draw("SAME");
  line020 -> Draw();
  line010 -> Draw();

  TLatex *latexTitle = new TLatex();
  latexTitle -> SetTextSize(0.050);
  latexTitle -> SetNDC();
  latexTitle -> SetTextFont(42);
  latexTitle -> DrawLatex(0.50, 0.80, Form("#color[633]{<#it{R}_{2}^{SP}> [0-20%%] = %4.3f}", func020 -> GetParameter(0)));
  latexTitle -> DrawLatex(0.50, 0.70, Form("#color[807]{<#it{R}_{2}^{SP}> [0-10%%] = %4.3f}", func010 -> GetParameter(0)));

  canvasResolution -> SaveAs("figures/sp_resolution.pdf");
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetLegend(TLegend *legend){
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.045);
}