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

void BackGround_MixedSpect()
{
  TFile *file = new TFile("AnalysisResults_New.root");
  TList *l2 = (TList*) file->Get("analysis-event-mixing/output");
  TList *l_MEPM = (TList*)l2->FindObject("PairsMuonMEPM_matchedMchMid");
  TList *l_MEPP = (TList*)l2->FindObject("PairsMuonMEPP_matchedMchMid");
  TList *l_MEMM = (TList*)l2->FindObject("PairsMuonMEMM_matchedMchMid");
  TH1D *h_Mass_MEPM = (TH1D*)l_MEPM->FindObject("Mass");
  TH1D *h_Mass_MEPP = (TH1D*)l_MEPP->FindObject("Mass");
  TH1D *h_Mass_MEMM = (TH1D*)l_MEMM->FindObject("Mass");
  
  TH3D *h_R2SP1_3d = (TH3D*)l_MEPM->FindObject("R2SP1_CentFT0C");
  TH3D *h_R2SP2_3d = (TH3D*)l_MEPM->FindObject("R2SP2_CentFT0C");
  TH3D *h_U2Q2_CentFT0C_ev1 = (TH3D*)l_MEPM->FindObject("U2Q2_CentFT0C_ev1");
  TH3D *h_U2Q2_CentFT0C_ev2 = (TH3D*)l_MEPM->FindObject("U2Q2_CentFT0C_ev2");
  TH2D *h_cos2DeltaPhiMu1 = (TH2D*)l_MEPM->FindObject("Mass_cos2DeltaPhiMu1");
  TH2D *h_cos2DeltaPhiMu2 = (TH2D*)l_MEPM->FindObject("Mass_cos2DeltaPhiMu2");

  TList *l3 = (TList*) file->Get("analysis-same-event-pairing/output");
  TList *l_SEPM = (TList*)l3->FindObject("PairsMuonSEPM_matchedMchMid");
  TList *l_SEPP = (TList*)l3->FindObject("PairsMuonSEPP_matchedMchMid");
  TList *l_SEMM = (TList*)l3->FindObject("PairsMuonSEMM_matchedMchMid");
  TH1D *h_Mass_SEPM = (TH1D*)l_SEPM->FindObject("Mass");
  TH1D *h_Mass_SEPP = (TH1D*)l_SEPP->FindObject("Mass");
  TH1D *h_Mass_SEMM = (TH1D*)l_SEMM->FindObject("Mass");

  
  //U2Q2
  //h_R2SP1->Draw();
  //h_U2Q2_CentFT0C_ev1->Draw();
  //h_U2Q2_CentFT0C_ev1->GetYaxis()->SetRange(0,10);      // Select the centrality [0,1] = 0-10 %, [1,2] = 10-20 % etc.
  TH2D *h_U2Q21FT0C = (TH2D*)h_U2Q2_CentFT0C_ev1->Project3D("zx"); // Project zx to plot mass vs. U2Q2
  TProfile* prof_U2Q21 = (TProfile*)h_U2Q21FT0C->ProfileX();  // Plot average of <u2q2> vs. mass

  //h_U2Q2_CentFT0C_ev2->GetYaxis()->SetRange(0,10);      // Select the centrality [0,1] = 0-10 %, [1,2] = 10-20 % etc.
  TH2D *h_U2Q22FT0C = (TH2D*)h_U2Q2_CentFT0C_ev2->Project3D("zx"); // Project zx to plot mass vs. U2Q2
  TProfile* prof_U2Q22 = (TProfile*)h_U2Q22FT0C->ProfileX();  // Plot average of <u2q2> vs. mass

  //R2SP
  //h_R2SP1_3d->GetYaxis()->SetRange(0,1);
  TH2D *h_R2SP1 = (TH2D*)h_R2SP1_3d->Project3D("xy");
  TH2D *h_R2SP2 = (TH2D*)h_R2SP2_3d->Project3D("xy");
  TProfile* prof_R2SP1 = (TProfile*)h_R2SP1->ProfileX();  //Resolution as a function of centrality
  TProfile* prof_R2SP2 = (TProfile*)h_R2SP2->ProfileX();  //Resolution as a function of centrality

  TCanvas *c2 = new TCanvas("c2","Canvas Example",200,10,600,480);
  prof_R2SP2->Draw();
  
  //cos2(Deltaphi)
  TProfile* prof_cos2DeltaPhiMu1 = (TProfile*)h_cos2DeltaPhiMu1->ProfileX();
  TProfile* prof_cos2DeltaPhiMu2 = (TProfile*)h_cos2DeltaPhiMu2->ProfileX();
  
  //prof_cos2DeltaPhiMu2->Draw();
  //h_U2Q21FT0C->Draw("colz");
  
  //prof_U2Q22->Draw();
   TFile *file2 = new TFile("v2_Bkg.root", "RECREATE");
   Double_t u1q1, u2q2, r1, r2, cosdphi1, cosdphi2;
   Double_t erru1q1, erru2q2, errr1, errr2, errcosdphi1, errcosdphi2;
   double v2_bkg, v2_bkg_er1,v2_bkg_er2,v2_bkg_er;
   TH1F* prof_v2bkg[10];
   for(int i = 0;i<10;i++) 
   {
   prof_v2bkg[i] = new TH1F("v2Bkg_vs_mass"," v_{2}^{bkg} vs. mass",125,0,5);
   prof_v2bkg[i]->SetTitle(" v_{2}^{bkg} vs. mass");
   //prof_v2bkg[i]->SetMarkerStyle(20);
   //prof_v2bkg[i]->SetMarkerColor(kBlue);
   //prof_v2bkg[i]->SetMarkerSize(0.5);
   }

   
   //prof_R2SP1->Rebin(9);
   //prof_R2SP2->Rebin(9);
   //cout<<prof_R2SP1->GetXaxis()->GetNbins()<<endl;
  for(int k = 1; k<=prof_R2SP1->GetXaxis()->GetNbins(); k++)
    {
   r1 = prof_R2SP1->GetBinContent(k);
   r2 = prof_R2SP2->GetBinContent(k);

   errr1 = prof_R2SP1->GetBinError(k);
   errr2 = prof_R2SP2->GetBinError(k);

   //r1 = 3.09529;
   //r2 = 3.09529;
   cout<<" R1 "<<r1<<endl;
   cout<<" R2 "<<r2<<endl;
   cout<<" erR1 "<<errr1<<endl;
   cout<<" erR2 "<<errr2<<endl;
   cout<<" Centrality bin "<<k<<endl;
    }

  //r1 = 1;
  //r2 = 1;
  //errr1 = 0;
  //errr2 = 0;
  
   for (int i = 1; i<=prof_U2Q21->GetXaxis()->GetNbins(); i++) {

    cosdphi1 = prof_cos2DeltaPhiMu1->GetBinContent(i);
    cosdphi2 = prof_cos2DeltaPhiMu2->GetBinContent(i);
    
    errcosdphi1 = prof_cos2DeltaPhiMu1->GetBinError(i);
    errcosdphi2 = prof_cos2DeltaPhiMu2->GetBinError(i);

    //r1 = prof_R2SP1->GetBinContent(i);
    //r2 = prof_R2SP2->GetBinContent(i);

    //errr1 = prof_R2SP1->GetBinError(i);
    //errr2 = prof_R2SP2->GetBinError(i);

    //prof_U2Q21->Scale(1/0.0801749);
    //prof_U2Q22->Scale(1/0.0801749);

    //prof_U2Q21->Multiply(prof_cos2DeltaPhiMu1);
    //prof_U2Q22->Multiply(prof_cos2DeltaPhiMu2);

   
    u1q1 = prof_U2Q21->GetBinContent(i);
    u2q2 = prof_U2Q21->GetBinContent(i);

    erru1q1 = prof_U2Q21->GetBinError(i);
    erru2q2 = prof_U2Q21->GetBinError(i);
    
    //v2_bkg = (u1q1*cosdphi1) + (u2q2*cosdphi2);
    if (r2 == 0 || r1 ==0) {
            prof_v2bkg[1]->SetBinContent(i, -999.);
            prof_v2bkg[1]->SetBinError(i, -999.);
        } else  {
      v2_bkg = u1q1 + u2q2 ;
      //v2_bkg = (u1q1/r1)*(cosdphi1) + (u2q2/r2)*(cosdphi2);
    v2_bkg_er1 = (u1q1/r1)*(cosdphi1) * TMath::Sqrt((erru1q1/u1q1)*(erru1q1/u1q1) + (errr1/r1)*(errr1/r1) + (errcosdphi1/cosdphi1)*(errcosdphi1/cosdphi1));
    v2_bkg_er2 = (u2q2/r2)*(cosdphi2) * TMath::Sqrt((erru2q2/u2q2)*(erru2q2/u2q2) + (errr2/r2)*(errr2/r2) + (errcosdphi2/cosdphi2)*(errcosdphi2/cosdphi2));
    v2_bkg_er = v2_bkg_er1 + v2_bkg_er2;
        }
    //if(isnan(v2_bkg) || isinf(v2_bkg)) continue;
    //cout<<" mass bin "<<i<< " v2_bkg "<<v2_bkg<<endl;
    prof_v2bkg[1]->SetBinContent(i,v2_bkg);
    prof_v2bkg[1]->SetBinError(i,v2_bkg_er);

 }

   
  //prof_cos2DeltaPhiMu1->Draw();
  prof_v2bkg[1]->Write();
  //prof_R2SP2->Draw();
  //prof_v2bkg->Rebin(2);
  //prof_v2bkg[k]->Draw();
 TH1D *h1 = prof_cos2DeltaPhiMu1->ProjectionX();
 //h1->Draw();
 //prof_cos2DeltaPhiMu1->Draw();
 //prof_v2bkg[1]->Draw("same");
  double SE = h_Mass_SEPM->Integral(h_Mass_SEPM->GetXaxis()->FindBin(1.0),h_Mass_SEPM->GetXaxis()->FindBin(2.5)) + h_Mass_SEPM->Integral(h_Mass_SEPM->GetXaxis()->FindBin(3.72),h_Mass_SEPM->GetXaxis()->FindBin(5.0));

  double ME = h_Mass_MEPM->Integral(h_Mass_MEPM->GetXaxis()->FindBin(1.0),h_Mass_MEPM->GetXaxis()->FindBin(2.5)) + h_Mass_MEPM->Integral(h_Mass_MEPM->GetXaxis()->FindBin(3.72),h_Mass_MEPM->GetXaxis()->FindBin(5.0));

    cout<<" Normalization factor "<<ME/SE <<endl;
    h_Mass_SEPM->SetMarkerStyle(20);
    h_Mass_SEPM->SetMarkerColor(kBlue);
    h_Mass_MEPM->SetMarkerStyle(22);
    h_Mass_MEPM->SetMarkerColor(kRed);

    TCanvas *c0 = new TCanvas("c0","Canvas Example",200,10,600,480);
    h_Mass_SEPM->Draw();
    h_Mass_MEPM->Scale(SE/ME);
    h_Mass_MEPM->Draw("same");

    TLegend *legend2_1 = new TLegend(0.6125205,0.5806881,0.8948445,0.7828001,NULL,"brNDC");
    legend2_1->SetTextFont(42);
    legend2_1->SetTextSize(0.04);
    legend2_1->SetLineColor(0);
    legend2_1->SetFillColor(0);
    legend2_1->AddEntry(h_Mass_SEPM,"Same event #mu^{-}#mu^{-}","p");
    legend2_1->AddEntry(h_Mass_MEPM,"Mixed event #mu^{-}#mu^{-}","p");
    legend2_1->Draw();

    TCanvas *c1 = new TCanvas("c1","Canvas Example",200,10,600,480);
    c1->Divide(2,2);
    c1->cd(1);
    prof_R2SP2->Draw();
    c1->cd(2);
    //prof_R2SP2->Draw();
    prof_U2Q22->Draw();
    c1->cd(3);
    prof_cos2DeltaPhiMu2->Draw();
    c1->cd(4);
    prof_v2bkg[1]->Draw();

    file2->Close();
  
}
