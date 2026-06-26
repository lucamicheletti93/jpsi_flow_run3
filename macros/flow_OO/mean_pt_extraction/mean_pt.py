import sys
import argparse
import yaml
import ROOT
import numpy as np
import os

def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--do_prefilter", help="Apply selections", action="store_true")
    parser.add_argument("--do_mean_pt", help="Apply selections", action="store_true")
    args = parser.parse_args()
    

    with open(args.cfgFileName, 'r') as yml_cfg:
        inputCfg = yaml.load(yml_cfg, yaml.FullLoader)

    if args.do_prefilter:
        prefilter(inputCfg)

    if args.do_mean_pt:
        get_mean_pt(inputCfg)


def prefilter(config):
    """
    function to apply pre-filters on data trees (see config_fit.yml)
    """
    fIn = ROOT.TFile(config["prefilter"]["data"], "READ")
    treeIn = fIn.Get(config["prefilter"]["treeIn"])
    rDataFrame = ROOT.RDataFrame(treeIn).Filter(config["prefilter"]["cuts"])
    fOutName = config["prefilter"]["data"].replace(".root", config["prefilter"]["suffix"] + ".root")
    rDataFrame.Snapshot(config["prefilter"]["treeOut"], fOutName)
    fIn.Close()

def get_mean_pt(config):
    """
    function get the mean pt from the unbinned fits
    """
    fInList = config["input"]["fInList"]
    print(fInList)

    for fInName in fInList:
        fIn = ROOT.TFile(fInName, "READ")
        hWeight = fIn.Get("sig_Jpsi_sw")
        print(f'{fInName}: {hWeight.GetMean():.2f}')

if __name__ == '__main__':
    main()