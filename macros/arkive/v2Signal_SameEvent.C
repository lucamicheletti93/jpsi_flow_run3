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

void v2Signal_SameEvent()
{
  const int nMassBins = 125;
  const double minMassRange = 0;
  const double maxMassRange = 5;
  
  TFile *file = new TFile("AnalysisResults_New.root");
  TList *l3 = (TList*) file->Get("analysis-same-event-pairing/output");
  TList *l_SEPM = (TList*)l3->FindObject("PairsMuonSEPM_matchedMchMid");
  TList *l_SEPP = (TList*)l3->FindObject("PairsMuonSEPP_matchedMchMid");
  TList *l_SEMM = (TList*)l3->FindObject("PairsMuonSEMM_matchedMchMid");
  TH1D *h_Mass_SEPM = (TH1D*)l_SEPM->FindObject("Mass");
  TH1D *h_Mass_SEPP = (TH1D*)l_SEPP->FindObject("Mass");
  TH1D *h_Mass_SEMM = (TH1D*)l_SEMM->FindObject("Mass");

  THnSparseF* hs = (THnSparseF*)l_SEPM->FindObject("Mass_Pt_centrFT0C_V2"); // 0: mass; 1: pT ; 2:FT0C_Centrality; 3:u2q2_A; 4:cos2(phi-psi)
  hs->GetAxis(2)->SetRange(1,7); //10-70% centrality
  TH3D *h_R2SP_3d = (TH3D*)hs->Projection(0,1,2);

  TH3D *h_cos2DeltaPhi_3d = (TH3D*)hs->Projection(0,1,4); //Mass_centFT0C_cos2DeltaPhi
  TH3D *h_u2q2_3d = (TH3D*)hs->Projection(0,2,3);   //Mass_centFT0C_u2q2

  //cout<<h_cos2DeltaPhi_3d->GetYaxis()->GetNbins()<<endl;
  h_cos2DeltaPhi_3d->GetYaxis()->SetRange(h_cos2DeltaPhi_3d->GetYaxis()->FindBin(1.0),h_cos2DeltaPhi_3d->GetYaxis()->FindBin(2.0)); //pT: 1.0-2.0 GeV
  TH2D *h_cos2DeltaPhi = (TH2D*)h_cos2DeltaPhi_3d->Project3D("zx");
  TProfile* prof_cos2DeltaPhi = (TProfile*)h_cos2DeltaPhi->ProfileX();

  TH2D *h_u2q2 = (TH2D*)h_u2q2_3d->Project3D("zx");
  TProfile* prof_U2Q2 = (TProfile*)h_u2q2->ProfileX();

  prof_cos2DeltaPhi->Scale(1/0.87);
  prof_cos2DeltaPhi->Rebin(4);
  //prof_U2Q2->Draw();

  prof_U2Q2->Rebin(4);
  double R1 = 0.876378;
  double R2 = 0.876378;

  TFile *fOut = new TFile("v2_results.root", "RECREATE");
  prof_cos2DeltaPhi->Write("histV2EP",TObject::kWriteDelete);
  prof_U2Q2->Write("histV2SP",TObject::kWriteDelete);
  fOut -> Close();


  /*
  TCanvas* c0 = new TCanvas("c0", "c0", 500, 500);
  hs->Projection(0,1,2)->Draw();

  TCanvas* c1 = new TCanvas("c1", "c1", 500, 500);
  h_R2SP_3d->Draw();

  TCanvas* c2 = new TCanvas("c2", "c2", 500, 500);
  hs->Projection(2)->Draw();
  */

  /*
  //TH3D *h_R2SP_3d = (TH3D*)l_SEPM->FindObject("Mass_centFT0C_R2SP");
  //TH3D *h_R2EP_3d = (TH3D*)l_SEPM->FindObject("Mass_centFT0C_R2EP");
  
  TH3D *h_cos2DeltaPhi_3d = (TH3D*)l_SEPM->FindObject("Mass_centFT0C_cos2DeltaPhi");
  TH3D *h_u2q2_3d = (TH3D*)l_SEPM->FindObject("Mass_centFT0C_u2q2");


  //TH2D *h_R2EP = (TH2D*)h_R2EP_3d->Project3D("zx");
  //TProfile* prof_R2EP = (TProfile*)h_R2EP->ProfileX();  //Resolution as a function of centrality
  //prof_R2EP->Draw();

  TH2D *h_u2q2 = (TH2D*)h_u2q2_3d->Project3D("zx");
  TProfile* prof_U2Q2 = (TProfile*)h_u2q2->ProfileX();
  //prof_cos2DeltaPhi->Draw();

  TH2D *h_cos2DeltaPhi = (TH2D*)h_cos2DeltaPhi_3d->Project3D("zx");
  TProfile* prof_cos2DeltaPhi = (TProfile*)h_cos2DeltaPhi->ProfileX();

  TH1::AddDirectory(kFALSE);
  TH1D *histV2EP = h_cos2DeltaPhi->ProjectionX();
  TH1D *histV2SP = prof_U2Q2->ProjectionX();

  TCanvas *c0 = new TCanvas("c0","Canvas Example",200,10,600,480);
  histV2SP->Draw();
  TFile *fOut = new TFile("v2_results.root", "RECREATE");
  histV2SP->Write("histV2SP",TObject::kWriteDelete);
  histV2EP->Write("histV2EP",TObject::kWriteDelete);
  fOut -> Close();

  TCanvas *c1 = new TCanvas("c1","Canvas Example",200,10,600,480);
  c1->Divide(2,2);
  c1->cd(1);
  prof_R2EP->Draw();
  c1->cd(2);
  prof_cos2DeltaPhi->Draw();
  c1->cd(3);
  prof_U2Q2->Draw();
  */

}
