import sys
import argparse
import yaml
import ROOT
import numpy as np
import array
import pandas as pd

sys.path.append('../../utils')
from plot_library import LoadStyle, SetHistStat, SetHistSyst, SetGraStat, SetGraSyst, SetLegend
from fast_analysis_utils import ExtractFromYaml


def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--run", help="Produce QA plots", action="store_true")
    parser.add_argument("--compare", help="Compare results", action="store_true")
    args = parser.parse_args()
    

    with open(args.cfgFileName, 'r') as yml_cfg:
        inputCfg = yaml.load(yml_cfg, yaml.FullLoader)

    if args.run:
        qa(inputCfg)

    if args.compare:
        compare()

def qa(config):
    """
    function to produce QA plots
    """
    LoadStyle()
    ROOT.gStyle.SetOptStat(False)
    BrJpsiToMuMu =  0.05961
    errBrJpsiToMuMu = 0.033e-2
    BrPsi2sToMuMu =  0.008
    errBrPsi2sToMuMu = 0.6e-3

    deltaRap = config["inputs"]["deltaRap"]

    print("***** Extract raw yield *****")
    dfJpsiRawYieldVsPt = pd.read_csv(config["inputs"]["fInRawYieldVsPt"], sep=' ')
    dfJpsiWidthVsPt = pd.read_csv(config["inputs"]["fInWidthVsPt"], sep=' ')
    dfJpsiMeanVsPt = pd.read_csv(config["inputs"]["fInMeanVsPt"], sep=' ')

    ptMin = dfJpsiRawYieldVsPt["x_min"]
    ptMax = dfJpsiRawYieldVsPt["x_max"]
    ptCenters = (ptMax + ptMin) / 2.
    ptWidths = (ptMax - ptMin) / 2.
    ptSystWidths = sqrtsSystWidths = np.repeat(0.2, len(ptWidths))
    ptEdges = np.append(ptMin.to_numpy(), ptMax.to_numpy()[len(ptMax)-1])

    jpsiRawYieldVsPt = dfJpsiRawYieldVsPt["val"] / (2 * ptWidths)
    jpsiStatRawYieldVsPt = dfJpsiRawYieldVsPt["stat"] / (2 * ptWidths)
    jpsiSystRawYieldVsPt = dfJpsiRawYieldVsPt["syst"] / (2 * ptWidths)

    jpsiWidthVsPt = dfJpsiWidthVsPt["val"]
    jpsiStatWidthVsPt = dfJpsiWidthVsPt["stat"]
    jpsiSystWidthVsPt = dfJpsiWidthVsPt["syst"]

    jpsiMeanVsPt = dfJpsiMeanVsPt["val"]
    jpsiStatMeanVsPt = dfJpsiMeanVsPt["stat"]
    jpsiSystMeanVsPt = dfJpsiMeanVsPt["syst"]


    histStatRawYieldVsPt = ROOT.TH1D("histStatRawYieldVsPt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{p}_{T}", len(ptEdges)-1, ptEdges)
    histSystRawYieldVsPt = ROOT.TH1D("histSystRawYieldVsPt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{p}_{T}", len(ptEdges)-1, ptEdges)

    histStatWidthVsPt = ROOT.TH1D("histStatWidthVsPt", ";#it{p}_{T} (GeV/#it{c});#sigma_{J/#psi}", len(ptEdges)-1, ptEdges)
    histSystWidthVsPt = ROOT.TH1D("histSystWidthVsPt", ";#it{p}_{T} (GeV/#it{c});#sigma_{J/#psi}", len(ptEdges)-1, ptEdges)

    histStatMeanVsPt = ROOT.TH1D("histStatMeanVsPt", ";#it{p}_{T} (GeV/#it{c});#mu_{J/#psi}", len(ptEdges)-1, ptEdges)
    histSystMeanVsPt = ROOT.TH1D("histSystMeanVsPt", ";#it{p}_{T} (GeV/#it{c});#mu_{J/#psi}", len(ptEdges)-1, ptEdges)

    for iBin in range(0, len(ptMin)):
        histStatRawYieldVsPt.SetBinContent(iBin+1, jpsiRawYieldVsPt[iBin])
        histStatRawYieldVsPt.SetBinError(iBin+1, jpsiStatRawYieldVsPt[iBin])
        histSystRawYieldVsPt.SetBinContent(iBin+1, jpsiRawYieldVsPt[iBin])
        histSystRawYieldVsPt.SetBinError(iBin+1, jpsiSystRawYieldVsPt[iBin])

        histStatWidthVsPt.SetBinContent(iBin+1, jpsiWidthVsPt[iBin])
        histStatWidthVsPt.SetBinError(iBin+1, jpsiStatWidthVsPt[iBin])
        histSystWidthVsPt.SetBinContent(iBin+1, jpsiWidthVsPt[iBin])
        histSystWidthVsPt.SetBinError(iBin+1, jpsiSystWidthVsPt[iBin])

        histStatMeanVsPt.SetBinContent(iBin+1, jpsiMeanVsPt[iBin])
        histStatMeanVsPt.SetBinError(iBin+1, jpsiStatMeanVsPt[iBin])
        histSystMeanVsPt.SetBinContent(iBin+1, jpsiMeanVsPt[iBin])
        histSystMeanVsPt.SetBinError(iBin+1, jpsiSystMeanVsPt[iBin])

    canvasRawYieldVsPt = ROOT.TCanvas("canvasRawYieldVsPt", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histSystRawYieldVsPt.Draw("E2P SAME")
    histStatRawYieldVsPt.Draw("EP SAME")
    canvasRawYieldVsPt.Update()

    canvasWidthVsPt = ROOT.TCanvas("canvasWidthVsPt", "", 800, 600)
    histSystWidthVsPt.Draw("E2P SAME")
    histStatWidthVsPt.Draw("EP SAME")
    canvasWidthVsPt.Update()

    canvasMeanVsPt = ROOT.TCanvas("canvasMeanVsPt", "", 800, 600)
    histSystMeanVsPt.Draw("E2P SAME")
    histStatMeanVsPt.Draw("EP SAME")
    canvasMeanVsPt.Update()

    input()
    fOutName = config["outputs"]["fOut"]
    fOut = ROOT.TFile(fOutName, "RECREATE")
    histSystRawYieldVsPt.Write()
    histStatRawYieldVsPt.Write()
    histSystWidthVsPt.Write()
    histStatWidthVsPt.Write()
    histSystMeanVsPt.Write()
    histStatMeanVsPt.Write()
    fOut.Close()

def compare():
    """
    function to compare results
    """
    LoadStyle()
    ROOT.gStyle.SetOptStat(False)

    fInCentr010 = ROOT.TFile("output/jpsi_qa_centrality_0_10.root")
    histStatRawYieldCentr010VsPt = fInCentr010.Get("histStatRawYieldVsPt")
    histSystRawYieldCentr010VsPt = fInCentr010.Get("histSystRawYieldVsPt")
    histStatWidthCentr010VsPt = fInCentr010.Get("histStatWidthVsPt")
    histStatMeanCentr010VsPt = fInCentr010.Get("histStatMeanVsPt")

    fInCentr020 = ROOT.TFile("output/jpsi_qa_centrality_0_20.root")
    histStatRawYieldCentr020VsPt = fInCentr020.Get("histStatRawYieldVsPt")
    histSystRawYieldCentr020VsPt = fInCentr020.Get("histSystRawYieldVsPt")
    histStatWidthCentr020VsPt = fInCentr020.Get("histStatWidthVsPt")
    histStatMeanCentr020VsPt = fInCentr020.Get("histStatMeanVsPt")

    SetHistStat(histStatRawYieldCentr010VsPt, 20, ROOT.kOrange+7)
    SetHistSyst(histSystRawYieldCentr010VsPt, 20, ROOT.kOrange+7)
    SetHistStat(histStatWidthCentr010VsPt, 20, ROOT.kOrange+7)
    SetHistStat(histStatMeanCentr010VsPt, 20, ROOT.kOrange+7)

    SetHistStat(histStatRawYieldCentr020VsPt, 20, ROOT.kRed+1)
    SetHistSyst(histSystRawYieldCentr020VsPt, 20, ROOT.kRed+1)
    SetHistStat(histStatWidthCentr020VsPt, 20, ROOT.kRed+1)
    SetHistStat(histStatMeanCentr020VsPt, 20, ROOT.kRed+1)



    canvasRawYieldVsPt = ROOT.TCanvas("canvasRawYieldVsPt", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histStatRawYieldCentr020VsPt.Draw("EP")
    histSystRawYieldCentr020VsPt.Draw("E2P SAME")
    histStatRawYieldCentr010VsPt.Draw("EP SAME")
    histSystRawYieldCentr010VsPt.Draw("E2P SAME")

    legendRawYieldVsPt = ROOT.TLegend(0.70, 0.70, 0.85, 0.85, " ", "brNDC")
    SetLegend(legendRawYieldVsPt)
    legendRawYieldVsPt.SetTextSize(0.045)
    legendRawYieldVsPt.AddEntry(histStatRawYieldCentr020VsPt, "0 #minus 20%", "P")
    legendRawYieldVsPt.AddEntry(histStatRawYieldCentr010VsPt, "0 #minus 10%", "P")
    legendRawYieldVsPt.Draw("SAME")

    canvasRawYieldVsPt.Update()

    canvasWidthVsPt = ROOT.TCanvas("canvasWidthVsPt", "", 800, 600)
    histStatWidthCentr020VsPt.GetYaxis().SetRangeUser(0, 0.150)
    histStatWidthCentr020VsPt.Draw("EP")
    histStatWidthCentr010VsPt.Draw("EP SAME")

    legendWidthVsPt = ROOT.TLegend(0.70, 0.25, 0.85, 0.40, " ", "brNDC")
    SetLegend(legendWidthVsPt)
    legendWidthVsPt.SetTextSize(0.045)
    legendWidthVsPt.AddEntry(histStatWidthCentr020VsPt, "0 #minus 20%", "P")
    legendWidthVsPt.AddEntry(histStatWidthCentr010VsPt, "0 #minus 10%", "P")
    legendWidthVsPt.Draw("SAME")

    canvasWidthVsPt.Update()

    canvasMeanVsPt = ROOT.TCanvas("canvasMeanVsPt", "", 800, 600)
    histStatMeanCentr020VsPt.GetYaxis().SetRangeUser(3.00, 3.15)
    histStatMeanCentr020VsPt.Draw("EP")
    histStatMeanCentr010VsPt.Draw("EP SAME")

    legendMeanVsPt = ROOT.TLegend(0.70, 0.70, 0.85, 0.85, " ", "brNDC")
    SetLegend(legendMeanVsPt)
    legendMeanVsPt.SetTextSize(0.045)
    legendMeanVsPt.AddEntry(histStatMeanCentr020VsPt, "0 #minus 20%", "P")
    legendMeanVsPt.AddEntry(histStatMeanCentr010VsPt, "0 #minus 10%", "P")
    legendMeanVsPt.Draw("SAME")

    canvasMeanVsPt.Update()

    input()

    canvasRawYieldVsPt.SaveAs("figures/qa_raw_yield.pdf")
    canvasWidthVsPt.SaveAs("figures/qa_width.pdf")
    canvasMeanVsPt.SaveAs("figures/qa_mean.pdf")


if __name__ == '__main__':
    main()