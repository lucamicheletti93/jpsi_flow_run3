import sys
import argparse
import yaml
import ROOT
import numpy as np
import array
import pandas as pd

sys.path.append('../utils')
from plot_library import LoadStyle, SetHistStat, SetHistSyst, SetGraStat, SetGraSyst, SetLegend
from fast_analysis_utils import ExtractFromYaml


def main():
    """
    function to compute the J/psi v2 in OO collisions
    """
    LoadStyle()
    ROOT.gStyle.SetOptStat(False)

    latexTitle = ROOT.TLatex()
    latexTitle.SetNDC()
    latexTitle.SetTextSize(0.05)
    latexTitle.SetTextFont(42)

    ptCentersOO = np.array([1.00, 3.00, 5.00, 7.00, 11.50])
    ptWidthsOO = np.array([1.00, 1.00, 1.00, 1.00, 3.50])
    ptSystWidthsOO = np.array([0.50, 0.50, 0.50, 0.50, 0.50])
    jpsiV2VsPtOO = np.array([0.01780, 0.00496, 0.04908, 0.07658, 0.20033])
    jpsiStatV2VsPtOO = np.array([0.01749, 0.01799, 0.03030, 0.04919, 0.07440])
    jpsiSystV2VsPtOO = np.array([0.00205, 0.00202, 0.00530, 0.01767, 0.04062])

    graStatV2VsPtOO = ROOT.TGraphErrors(len(ptCentersOO), ptCentersOO, jpsiV2VsPtOO, ptWidthsOO, jpsiStatV2VsPtOO)
    graSystV2VsPtOO = ROOT.TGraphErrors(len(ptCentersOO), ptCentersOO, jpsiV2VsPtOO, ptSystWidthsOO, jpsiSystV2VsPtOO)

    SetGraStat(graStatV2VsPtOO, 20, ROOT.kRed+1)
    SetGraSyst(graSystV2VsPtOO, 20, ROOT.kRed+1)

    ptCentersPbPb, ptWidthsPbPb, jpsiV2VsPtPbPb, jpsiStatV2VsPtPbPb, jpsiSystV2VsPtPbPb = ExtractFromYaml("HEP_Data/jpsi_v2_10_30_PbPb.yaml")
    jpsiPtSystWidthsPbPb = np.repeat(0.2, len(ptCentersPbPb))

    graStatV2VsPtPbPb = ROOT.TGraphErrors(len(ptCentersPbPb), np.array(ptCentersPbPb), np.array(jpsiV2VsPtPbPb), np.array(ptWidthsPbPb), np.array(jpsiStatV2VsPtPbPb))
    graSystV2VsPtPbPb = ROOT.TGraphErrors(len(ptCentersPbPb), np.array(ptCentersPbPb), np.array(jpsiV2VsPtPbPb), jpsiPtSystWidthsPbPb, np.array(jpsiSystV2VsPtPbPb))

    SetGraStat(graStatV2VsPtPbPb, 20, ROOT.kBlack, 0.7)
    SetGraSyst(graSystV2VsPtPbPb, 20, ROOT.kBlack, 0.7)



    ptCenterspPb, ptWidthspPb, jpsiV2VsPtpPb, jpsiStatV2VsPtpPb, jpsiSystV2VsPtpPb = ExtractFromYaml("HEP_Data/jpsi_v2_pPb.yaml")
    jpsiPtSystWidthspPb = np.repeat(0.2, len(ptCenterspPb))

    graStatV2VsPtpPb = ROOT.TGraphErrors(len(ptCenterspPb), np.array(ptCenterspPb), np.array(jpsiV2VsPtpPb), np.array(ptWidthspPb), np.array(jpsiStatV2VsPtpPb))
    graSystV2VsPtpPb = ROOT.TGraphErrors(len(ptCenterspPb), np.array(ptCenterspPb), np.array(jpsiV2VsPtpPb), jpsiPtSystWidthspPb, np.array(jpsiSystV2VsPtpPb))

    SetGraStat(graStatV2VsPtpPb, 20, ROOT.kAzure+4, 0.7)
    SetGraSyst(graSystV2VsPtpPb, 20, ROOT.kAzure+4, 0.7)

    ptCentersPbp, ptWidthsPbp, jpsiV2VsPtPbp, jpsiStatV2VsPtPbp, jpsiSystV2VsPtPbp = ExtractFromYaml("HEP_Data/jpsi_v2_Pbp.yaml")
    jpsiPtSystWidthsPbp = np.repeat(0.2, len(ptCentersPbp))

    graStatV2VsPtPbp = ROOT.TGraphErrors(len(ptCentersPbp), np.array(ptCentersPbp), np.array(jpsiV2VsPtPbp), np.array(ptWidthsPbp), np.array(jpsiStatV2VsPtPbp))
    graSystV2VsPtPbp = ROOT.TGraphErrors(len(ptCentersPbp), np.array(ptCentersPbp), np.array(jpsiV2VsPtPbp), jpsiPtSystWidthsPbp, np.array(jpsiSystV2VsPtPbp))

    SetGraStat(graStatV2VsPtPbp, 20, ROOT.kGreen+1, 0.7)
    SetGraSyst(graSystV2VsPtPbp, 20, ROOT.kGreen+1, 0.7)



    lineUnityVsPt = ROOT.TLine(0, 0, 15, 0)
    lineUnityVsPt.SetLineColor(ROOT.kGray+2)
    lineUnityVsPt.SetLineWidth(2)
    lineUnityVsPt.SetLineStyle(ROOT.kDashed)

    canvasV2VsPt = ROOT.TCanvas("canvasV2VsPt", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    histGridV2VsPt = ROOT.TH2D("histGridV2VsPt", ";#it{p}_{T} (GeV/#it{c});#it{v}_{2}^{J/#psi}", 100, 0, 15, 100, -0.10, 0.50)
    histGridV2VsPt.Draw()
    graStatV2VsPtOO.Draw("EP SAME")
    graSystV2VsPtOO.Draw("E2P SAME")
    graStatV2VsPtPbPb.Draw("EP SAME")
    graSystV2VsPtPbPb.Draw("E2P SAME")
    graStatV2VsPtpPb.Draw("EP SAME")
    graSystV2VsPtpPb.Draw("E2P SAME")
    graStatV2VsPtPbp.Draw("EP SAME")
    graSystV2VsPtPbp.Draw("E2P SAME")
    
    lineUnityVsPt.Draw("SAME")

    legendV2VsPt = ROOT.TLegend(0.20, 0.60, 0.40, 0.85, " ", "brNDC")
    SetLegend(legendV2VsPt)
    legendV2VsPt.SetTextSize(0.035)
    legendV2VsPt.AddEntry(graStatV2VsPtOO, "#sqrt{#it{s}} = 5.36 TeV, OO, 0#minus70%", "FP")
    legendV2VsPt.AddEntry(graSystV2VsPtPbPb, "#sqrt{#it{s}} = 5.02 TeV, Pb-Pb, 10#minus30%", "FP")
    legendV2VsPt.AddEntry(graSystV2VsPtpPb, "#sqrt{#it{s}} = 8.16 TeV, pPb, 2.03 < #it{y} < 3.53", "FP")
    legendV2VsPt.AddEntry(graSystV2VsPtPbp, "#sqrt{#it{s}} = 8.16 TeV, pPb, -4.46 < #it{y} < -2.96", "FP")
    legendV2VsPt.Draw("SAME")

    latexTitle.DrawLatex(0.22, 0.87, "ALICE Preliminary, J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")

    canvasV2VsPt.Update()
    

    input()

    canvasV2VsPt.SaveAs("flow_different_systems.pdf")
    exit()

if __name__ == '__main__':
    main()



