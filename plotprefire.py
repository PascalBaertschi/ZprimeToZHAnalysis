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

def ensureTFile(filename,option='READ'):
  """Open TFile, checking if the file in the given path exists."""
  if not os.path.isfile(filename):
    print '>>> ERROR! plotseff.ensureTFile: File in path "%s" does not exist!!'%(filename)
    exit(1)
  file = TFile(filename,option)
  if not file or file.IsZombie():
    print '>>> ERROR! plotseff.ensureTFile Could not open file by name "%s"'%(filename)
    exit(1)
  return file

########## SETTINGS ##########

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-v", "--variable", action="store", type="string", dest="variable", default="")
parser.add_option("-c", "--cut", action="store", type="string", dest="cut", default="")
parser.add_option("-y", "--year",action="store", type="string", dest="year", default="2017")
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
var         = options.variable

if year in ['2016','2017','2018']:
    NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/Ntuples%s/"%(year)
else:
    NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/"
OUTPUTDIR = "/work/pbaertsc/heavy_resonance/Analysis/plots_%s/"%(year)

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
#sign = ['XZH_M800','XZH_M1000','XZH_M1200','XZH_M1400','XZH_M1600', 'XZH_M1800', 'XZH_M2000', 'XZH_M2500','XZH_M3000','XZH_M3500', 'XZH_M4000','XZH_M4500','XZH_M5000']
sign_list = ["XZH_M1000","XZH_M2000","XZH_M3000"]
sign_list_vbf = ["XZHVBF_M1000","XZHVBF_M2000","XZHVBF_M3000"]



def plot(var, cut, signal):
    sign = signal
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
    if year in ['2016','2017','2018']:
        NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/Ntuples%s/"%(year)
    else:
        NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/"
    #showSignal = False if 'SB' in cut or 'TR' in cut else True #'SR' in channel or channel=='qqqq'#or len(channel)==5
    if var in ['dijet_VBF_mass','deltaR_VBF','deltaR_HVBFjet1','deltaR_HVBFjet2']:
        showSignal = False
    else:
        showSignal = True
    if treeRead:
        for k in sorted(selection.keys(), key=len, reverse=True):
            if k in cut: cut = cut.replace(k, selection[k])

    print "Plotting from", ("tree" if treeRead else "file"), var, "in", channel, "channel with:"
    print "  cut    :", cut
    
    ### Create and fill MC histograms ###
    # Create dict
    file = {}
    tree = {}
    hist = {}
    ### Create and fill MC histograms ###
    tree[sign] = TChain("tree")
    for j, ss in enumerate(sample[sign]['files']):
        tree[sign].Add(NTUPLEDIR + ss + ".root")
    if variable[var]['nbins']>0: 
      min_value = variable[var]['min']
      max_value =variable[var]['max']
      title = variable[var]['title']
      hist[sign] = TH1F(sign, ";"+title+";Events;"+('log' if variable[var]['log'] else ''), variable[var]['nbins'], min_value, max_value)
    else: hist[s] = TH1F(sign, ";"+variable[var]['title'], len(variable[var]['bins'])-1, array('f', variable[var]['bins']))
    hist[sign].Sumw2()
    cutstring = "%s" % eventWeightLuminame + ("*("+cut+")") 
    tree[sign].Project(sign, var, cutstring)
    if not tree[sign].GetTree()==None: hist[sign].SetOption("%s" % tree[sign].GetTree().GetEntriesFast())
    hist[sign].SetFillColor(sample[sign]['fillcolor'])
    hist[sign].SetFillStyle(sample[sign]['fillstyle'])
    hist[sign].SetLineColor(sample[sign]['linecolor'])
    hist[sign].SetLineStyle(sample[sign]['linestyle'])

    hist[sign].SetLineWidth(3)
    sample[sign]['plot'] = True

    # Legend
    leg = TLegend(0.6, 0.6, 0.9, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0) #1001
    leg.SetFillColor(0)
    leg.AddEntry(hist[sign], sample[sign]['label'], "fl")
        
    leg.SetY1(0.9-leg.GetNRows()*0.05)
    
    
    # --- Display ---
    c1 = TCanvas("c1", hist.values()[0].GetXaxis().GetTitle(), 800, 800)

    c1.SetTopMargin(0.06)
    c1.SetRightMargin(0.05)
    c1.SetTicks(1, 1)
    hist[sign].Draw("HIST")
    
    leg.Draw()
    drawCMS(LUMI,year,"Preliminary")
    drawRegion('XVH'+channel, True)
    drawAnalysis(channel)
    
    c1.Update()
    
    if gROOT.IsBatch():
        varname = var.replace('.', '_').replace('()', '')
        if not os.path.exists(OUTPUTDIR+channel_name): os.makedirs(OUTPUTDIR+channel_name)
        c1.Print("%s"% OUTPUTDIR +channel_name+"/"+varname+sign+".png")
        c1.Print("%s"% OUTPUTDIR +channel_name+"/"+varname+sign+".pdf")

    

var = 'PrefireWeight'
for cut in ['nninc','eeinc','mminc','nnincVBF','eeincVBF','mmincVBF']:
  if not 'VBF' in cut:
    for signal in sign_list:
      plot(var,cut,signal)
  else:
    for signal in sign_list_vbf:
      plot(var,cut,signal)

  



