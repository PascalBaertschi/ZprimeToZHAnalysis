#! /usr/bin/env python

import os, multiprocessing
import copy
import math
import numpy as np
from array import array
from ROOT import ROOT, gROOT, gStyle, gRandom, TSystemDirectory
from ROOT import TFile, TChain, TTree, TCut, TF1, TH1F, TH2F, THStack
from ROOT import TGraph, TGraphErrors, TGraphAsymmErrors, TVirtualFitter
from ROOT import TStyle, TCanvas, TPad, TCanvas
from ROOT import TLegend, TLatex, TText, TLine, TPaveText
from xsections import xsection

from samples import sample_2016, sample_2017, sample_2018, sample
from variables import variable
from selections import selection
from utils import *

########## SETTINGS ##########

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-v", "--variable", action="store", type="string", dest="variable", default="")
parser.add_option("-c", "--cut", action="store", type="string", dest="cut", default="")
parser.add_option("-y", "--year",action="store", type="string", dest="year", default="combined")
parser.add_option("-w", "--weighted",action="store_false", default=True, dest="weighted")
parser.add_option("-n", "--norm", action="store_true", default=False, dest="norm")
parser.add_option("-s", "--syst", action="store_true", default=False, dest="syst")
parser.add_option("-l", "--checkcut", action="store_true", default=False, dest="checkcut")
parser.add_option("-t", "--top", action="store_true", default=False, dest="top")
parser.add_option("-a", "--all", action="store_true", default=False, dest="all")
parser.add_option("-x", "--tagger", action="store_true", default=False, dest="tagger")
parser.add_option("-b", "--bash", action="store_true", default=False, dest="bash")
parser.add_option("-B", "--blind", action="store_true", default=False, dest="blind")
parser.add_option("-f", "--final", action="store_true", default=False, dest="final")
(options, args) = parser.parse_args()
if options.bash: gROOT.SetBatch(True)

########## SETTINGS ##########

#gROOT.SetBatch(True)
#gROOT.ProcessLine("TSystemDirectory::SetDirectory(0)")
#gROOT.ProcessLine("TH1::AddDirectory(kFALSE);")
gStyle.SetOptStat(0)
#TSystemDirectory.SetDirectory(0)


year        = options.year
cut         = options.cut

if year in ['2016','2017','2018']:
    NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/Ntuples%s/"%(year)
else:
    NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/"
OUTPUTDIR = "/work/pbaertsc/heavy_resonance/ZprimeToZHAnalysis/plots_%s/"%(year)

SIGNAL      = 1 # Signal magnification factor
RATIO       = 4 # 0: No ratio plot; !=0: ratio between the top and bottom pads
NORM        = options.norm
PARALLELIZE = False
BLIND = True

if year=='2016':
    LUMI=35920.
    sample=sample_2016
elif year=='2017':
    LUMI=41530.
    sample=sample_2017
elif year=='2018':
    LUMI=59740.
    sample=sample_2018
elif year=='combined':
    LUMI=137190
########## SAMPLES ##########
data = ['data_obs']
back = ["VV", "ST", "TTbarSL", "WJetsToLNu_HT", "DYJetsToNuNu_HT", "DYJetsToLL_HT"] 
sign = ['XZH_M800','XZH_M1000','XZH_M1200','XZH_M1400','XZH_M1600', 'XZH_M1800', 'XZH_M2000', 'XZH_M2500','XZH_M3000','XZH_M3500', 'XZH_M4000','XZH_M4500','XZH_M5000']
#
########## ######## ##########

jobs = []

def plot(var, cut):
    ### Preliminary Operations ###
    treeRead = not 'cutflow' in var
    #treeRead = False
    channel = cut
    if 'inc' in cut:
        eventWeightLuminame = 'eventWeightLumi_nobtag'
    else:
        eventWeightLuminame = 'eventWeightLumi'
    if "SB" in cut or "SR" in cut:
        channel_name = cut[:-2]
    else:
        channel_name = cut
    cut = 'mmincSB_noisocut'
    isBlind = BLIND and 'SR' in channel
    if year in ['2016','2017','2018']:
        NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/Ntuples%s/"%(year)
    else:
        NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/"
    stype = "HVT model B"
  
    
    if treeRead:
        for k in sorted(selection.keys(), key=len, reverse=True):
            if k in cut: cut = cut.replace(k, selection[k])
  
    # Determine Primary Dataset
    pd = []
    if any(w in cut for w in ['mn', 'mm', 'isZtoMM', 'isWtoMN']): pd += [x for x in sample['data_obs']['files'] if 'SingleMuon' in x]
    if any(w in cut for w in ['en', 'ee', 'isWtoEN', 'isZtoEE', 'emqq', 'isTtoEM']): pd += [x for x in sample['data_obs']['files'] if ('SingleElectron' in x or 'EGamma' in x)]
    if any(w in cut for w in ['nn', 'isZtoNN']): pd += [x for x in sample['data_obs']['files'] if 'MET' in x]
    if len(pd)==0: raw_input("Warning: Primary Dataset not recognized, continue?")
    if var == "PrefireWeight":
        pd = []

    print "Plotting from", ("tree" if treeRead else "file"), var, "in", channel, "channel with:"
    print "  dataset:", pd
    print "  cut    :", cut
    
    ### Create and fill MC histograms ###
    # Create dict
    file = {}
    tree_noisocut = {}
    tree_isocut = {}
    hist_noisocut = {}
    hist_isocut = {}
    cutstring_noisocut = "%s" % eventWeightLuminame + ("*(Mu1_highPtId==2 &&"+cut+")") 
    cutstring_isocut = "%s" % eventWeightLuminame + ("*(Mu1_relIso!=-1. && Mu1_relIso < 0.1 && Mu1_highPtId==2 &&"+cut+")") 
    #cutstring_noisocut = "%s" % eventWeightLuminame + ("*(Mu1_highPtId==1 &&"+cut+")") 
    #cutstring_isocut = "%s" % eventWeightLuminame + ("*(Mu1_relIso!=-1. && Mu1_relIso < 0.1 && Mu1_highPtId==1 &&"+cut+")")
    """
    if var == 'Mu1_pt':
        cutstring_noisocut = "%s" % eventWeightLuminame + ("*(Mu1_highPtId==2 &&"+cut+")") 
        cutstring_isocut = "%s" % eventWeightLuminame + ("*(Mu1_relIso!=-1. && Mu1_relIso < 0.1 && Mu1_highPtId==2 &&"+cut+")") 
    elif var == 'Mu2_pt':
        cutstring_noisocut = "%s" % eventWeightLuminame + ("*(Mu2_highPtId==1 &&"+cut+")") 
        cutstring_isocut = "%s" % eventWeightLuminame + ("*(Mu2_relIso!=-1. && Mu2_relIso < 0.1 && Mu2_highPtId==1 &&"+cut+")") 
    """
    ### Create and fill MC histograms ###
    for i, s in enumerate(sign):
        if treeRead: # Project from tree
            tree_noisocut[s] = TChain("tree")
            for j, ss in enumerate(sample[s]['files']):
                if not 'data' in s or ('data' in s and ss in pd):
                    tree_noisocut[s].Add(NTUPLEDIR + ss + ".root")
            if variable[var]['nbins']>0: 
                min_value = variable[var]['min']
                max_value =variable[var]['max']
                title = variable[var]['title']
                hist_noisocut[s] = TH1F(s, ";"+title+";Events;"+('log' if variable[var]['log'] else ''), variable[var]['nbins'], min_value, max_value)
            else: 
                hist_noisocut[s] = TH1F(s, ";"+variable[var]['title'], len(variable[var]['bins'])-1, array('f', variable[var]['bins']))
            hist_noisocut[s].Sumw2()
            tree_noisocut[s].Project(s, var, cutstring_noisocut)
        hist_noisocut[s].Scale(sample[s]['weight'] if hist_noisocut[s].Integral() >= 0 else 0)
   
    for i, s in enumerate(sign):
        if treeRead: # Project from tree
            tree_isocut[s] = TChain("tree")
            for j, ss in enumerate(sample[s]['files']):
                if not 'data' in s or ('data' in s and ss in pd):
                    tree_isocut[s].Add(NTUPLEDIR + ss + ".root")
            if variable[var]['nbins']>0: 
                min_value = variable[var]['min']
                max_value =variable[var]['max']
                title = variable[var]['title']
                hist_isocut[s] = TH1F(s, ";"+title+";Events;"+('log' if variable[var]['log'] else ''), variable[var]['nbins'], min_value, max_value)
            else: 
                hist_isocut[s] = TH1F(s, ";"+variable[var]['title'], len(variable[var]['bins'])-1, array('f', variable[var]['bins']))
            hist_isocut[s].Sumw2()
            tree_isocut[s].Project(s,var,cutstring_isocut)
        hist_isocut[s].Scale(sample[s]['weight'] if hist_isocut[s].Integral() >= 0 else 0)

        
    
    for i, s in enumerate(sign): addOverflow(hist_noisocut[s], False) # Add overflow
    for i, s in enumerate(sign): addOverflow(hist_isocut[s], False) # Add overflow
    for i, s in enumerate(sign): hist_noisocut[s].SetLineWidth(3)
    for i, s in enumerate(sign): hist_isocut[s].SetLineWidth(3)
    for i, s in enumerate(sign): sample[s]['plot'] = True#sample[s]['plot'] and s.startswith(channel[:2])

    
    # Create stack
    SigSum_noisocut = hist_noisocut[sign[0]].Clone("SigSum_noisocut")
    SigSum_noisocut.Reset("MICES")
    for i, s in enumerate(sign): SigSum_noisocut.Add(hist_noisocut[s], 1)
    SigSum_isocut = hist_isocut[sign[0]].Clone("SigSum_isocut")
    SigSum_isocut.Reset("MICES")
    for i, s in enumerate(sign): SigSum_isocut.Add(hist_isocut[s], 1)


    SigSum_noisocut.Divide(SigSum_isocut)

    SigSum_noisocut.GetYaxis().SetRangeUser(0.98,1.05)
    c1 = TCanvas("c1", "Efficiency Mu Iso", 800,800)
    SigSum_noisocut.Draw("HIST")
    drawCMS(LUMI,year,"Preliminary")
    drawRegion('XVH'+channel, True)
    drawAnalysis(channel)
    c1.Print("plotsSF/MuIsoEff_highPtId.png")
    c1.Print("plotsSF/MuIsoEff_highPtId.pdf")
    #c1.Print("plotsSF/MuIsoEff_trackerhighPtId.png")
    #c1.Print("plotsSF/MuIsoEff_trackerhighPtId.pdf")


plot(options.variable, options.cut)
