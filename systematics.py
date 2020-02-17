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
parser.add_option("-s", "--syst", action="store_true", default=False, dest="syst")
parser.add_option("-e", "--effb", action="store_true", default=False, dest="effb")
parser.add_option("-j", "--jec", action="store_true", default=False, dest="jec")
parser.add_option("-m", "--jmc", action="store_true", default=False, dest="jmc")
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
sign_VBF = ['XZHVBF_M800','XZHVBF_M1000','XZHVBF_M1200','XZHVBF_M1400','XZHVBF_M1600', 'XZHVBF_M1800', 'XZHVBF_M2000', 'XZHVBF_M2500','XZHVBF_M3000','XZHVBF_M3500', 'XZHVBF_M4000','XZHVBF_M4500','XZHVBF_M5000']
#sign = ['XZH_M2000']
#sign = []


########## ######## ##########

jobs = []

def plot(var, cut, nm1=False):
    ### Preliminary Operations ###
    treeRead = True
    channel = cut
    BTagAK4deepup = False
    BTagAK4deepdown = False
    BTagAK4deep = False
    BTagAK8deepup = False
    BTagAK8deepdown = False
    BTagAK8deep = False
    if 'inc' in cut:
        eventWeightLuminame = 'eventWeightLumi_nobtag'
    else:
        eventWeightLuminame = 'eventWeightLumi'
    if "SB" in cut or "SR" in cut:
        channel_name = cut[:-2]
    else:
        channel_name = cut
    isBlind = BLIND and 'SR' in channel

    if year in ['2016','2017','2018']:
        NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/Ntuples%s/"%(year)
    else:
        NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/"
    #showSignal = False if 'SB' in cut or 'TR' in cut else True #'SR' in channel or channel=='qqqq'#or len(channel)==5
    if var in ['dijet_VBF_mass','deltaR_VBF','deltaR_HVBFjet1','deltaR_HVBFjet2']:
        showSignal = False
    else:
        showSignal = True
    if "cutflow" in var:
        showSignal = False
    if var in ['X_mass_jesUp','X_mass_jesDown','X_mass_jerUp','X_mass_jerDown','X_mass_nom','X_mass_MET_jesUp','X_mass_MET_jesDown','X_mass_MET_jerUp','X_mass_MET_jerDown','X_mass_MET_nom','H_mass_jmsUp','H_mass_jmsDown','H_mass_jmrUp','H_mass_jmrDown','H_mass_nom']:
        back = ["VV"]
    else:
        back = ["VV", "ST", "TTbarSL", "WJetsToLNu_HT", "DYJetsToNuNu_HT", "DYJetsToLL_HT"] 
    stype = "HVT model B"
    if channel.endswith('up'):
        TopBTagAK4deepup = True
        channel = channel[:-2]
        cut = channel
    if channel.endswith('down'):
        TopBTagAK4deepdown = True
        channel = channel[:-4]
        cut = channel
    if treeRead:
        for k in sorted(selection.keys(), key=len, reverse=True):
            if k in cut: cut = cut.replace(k, selection[k])
  
    # Determine Primary Dataset
    pd = []

    if var=='BTagAK4Weightdeep_up':
        var = 'X_mass'
        BTagAK4deepup = True
    elif var=='BTagAK4Weight_deep_down':
        var = 'X_mass'
        BTagAK4deepdown = True
    elif var=='BTagAK4Weight_deep':
        var = 'X_mass'
        BTagAK4deep = True
    elif var=='BTagAK8Weight_deep_up':
        var = 'X_mass'
        BTagAK8deepup = True
    elif var=='BTagAK8Weight_deep_down':
        var = 'X_mass'
        BTagAK8deepdown = True
    elif var=='BTagAK8Weight_deep':
        var = 'X_mass'
        BTagAK8deep = True   

    print "Plotting from", ("tree" if treeRead else "file"), var, "in", channel, "channel with:"
    print "  dataset:", pd
    print "  cut    :", cut
    
    ### Create and fill MC histograms ###
    # Create dict
    file = {}
    tree = {}
    hist = {}
    ### Create and fill MC histograms ###
    for i, s in enumerate(back+sign+sign_VBF):
        if treeRead: # Project from tree
            tree[s] = TChain("tree")
            for j, ss in enumerate(sample[s]['files']):
                if not 'data' in s or ('data' in s and ss in pd):
                    tree[s].Add(NTUPLEDIR + ss + ".root")
            if variable[var]['nbins']>0: 
                if 'X_mass' in var:
                    min_value = 0.0
                    max_value = 6000.
                else:
                    min_value = variable[var]['min']
                    max_value =variable[var]['max']
                title = variable[var]['title']
                hist[s] = TH1F(s, ";"+title+";Events;"+('log' if variable[var]['log'] else ''), variable[var]['nbins'], min_value, max_value)
            else: hist[s] = TH1F(s, ";"+variable[var]['title'], len(variable[var]['bins'])-1, array('f', variable[var]['bins']))
            hist[s].Sumw2()
            cutstring = "%s" % eventWeightLuminame + ("*("+cut+")") 
            if var=='LeptonWeightUp':
                cutstring = "(%s * LeptonWeightUp/LeptonWeight)" % eventWeightLuminame  + ("*("+cut+")") 
            elif var=='LeptonWeightDown':
                cutstring = "(%s * LeptonWeightDown/LeptonWeight)" % eventWeightLuminame  + ("*("+cut+")") 
            elif var=='TriggerWeightUp':
                cutstring = "(%s * TriggerWeightUp/TriggerWeight)" % eventWeightLuminame + ("*("+cut+")") 
            elif var=='TriggerWeightDown':
                cutstring = "(%s * TriggerWeightDown/TriggerWeight)" % eventWeightLuminame  + ("*("+cut+")") 
            #division by BTagAk8Weight_deep is because the weighted samples are used
            elif BTagAK4deep:
                cutstring = "(%s * BTagAK4Weight_deep)" % eventWeightLuminame + ("*("+cut+")")
            elif BTagAK4deepup:
                cutstring = "(%s * BTagAK4Weight_deep_up)" % eventWeightLuminame + ("*("+cut+")") 
            elif BTagAK4deepdown:
                cutstring = "(%s * BTagAK4Weight_deep_down)" % eventWeightLuminame + ("*("+cut+")") 
            elif BTagAK8deepup:
                cutstring = "(%s * BTagAK8Weight_deep_up/BTagAK8Weight_deep)" % eventWeightLuminame  + ("*("+cut+")") 
            elif BTagAK8deepdown:
                cutstring = "(%s * BTagAK8Weight_deep_down/BTagAK8Weight_deep)" % eventWeightLuminame  + ("*("+cut+")")    
            tree[s].Project(s, var, cutstring)
            if not tree[s].GetTree()==None: hist[s].SetOption("%s" % tree[s].GetTree().GetEntriesFast())
            #print tree[s].GetTree().GetEntriesFast()
        hist[s].Scale(sample[s]['weight'] if hist[s].Integral() >= 0 else 0)
        hist[s].SetFillColor(sample[s]['fillcolor'])
        hist[s].SetFillStyle(sample[s]['fillstyle'])
        hist[s].SetLineColor(sample[s]['linecolor'])
        hist[s].SetLineStyle(sample[s]['linestyle'])
    
    hist['BkgSum'] = hist['data_obs'].Clone("BkgSum") if 'data_obs' in hist else hist[back[0]].Clone("BkgSum")
    hist['BkgSum'].Reset("MICES")
    hist['BkgSum'].SetFillStyle(3003)
    hist['BkgSum'].SetFillColor(1)
    for i, s in enumerate(back): 
        hist['BkgSum'].Add(hist[s])
       # Set histogram style
    
    for i, s in enumerate(back+sign+sign_VBF+['BkgSum']): addOverflow(hist[s], False) # Add overflow
    for i, s in enumerate(sign+sign_VBF): hist[s].SetLineWidth(3)
    for i, s in enumerate(sign+sign_VBF): sample[s]['plot'] = True#sample[s]['plot'] and s.startswith(channel[:2])

    
    # Create stack
    bkg = THStack("Bkg", ";"+hist['BkgSum'].GetXaxis().GetTitle()+";Events")
    for i, s in enumerate(back): bkg.Add(hist[s])
    return hist

def calc_syst():
    Integral = 0.
    Integral_up = 0.
    Integral_down = 0.
    syst_trig = []
    syst_elec = []
    syst_muon = []
    for cut in ['nnbbSR','eebbSR','mmbbSR','nn0bSR','ee0bSR','mm0bSR','nnbbVBFSR','eebbVBFSR','mmbbVBFSR','nn0bVBFSR','ee0bVBFSR','mm0bVBFSR']:
        for var in ['LeptonWeightUp','LeptonWeightDown','LeptonWeight']:
            if 'Up' in var:
                Integral_up = plot(var,cut)['BkgSum'].Integral()
            elif 'Down' in var:
                Integral_down = plot(var,cut)['BkgSum'].Integral()
            else:
                Integral = plot(var,cut)['BkgSum'].Integral()
                syst_up = Integral_up/Integral-1.
                syst_down = Integral/Integral_down-1.
                syst_total = (syst_up+syst_down)/2
                print "up syst for %s:" % cut,syst_up
                print "down syst for %s:" % cut,syst_down
                print "total syst for %s" % cut,syst_total
                if 'ee' in cut:
                    syst_elec.append('%s : %.3f' % (cut[:-2],syst_total))
                    syst_muon.append('%s : %.3f' % (cut[:-2],0.000))
                elif 'mm' in cut:
                    syst_elec.append('%s : %.3f' % (cut[:-2],0.000))
                    syst_muon.append('%s : %.3f' % (cut[:-2],syst_total))
                else:
                    syst_elec.append('%s : %.3f' % (cut[:-2],0.00))
                    syst_muon.append('%s : %.3f' % (cut[:-2],0.00))
        Integral = 0.
        Integral_up = 0.
        Integral_down = 0.
        for var in ['TriggerWeightUp','TriggerWeightDown','TriggerWeight']:
            if 'Up' in var:
                Integral_up = plot(var,cut)['BkgSum'].Integral()
            elif 'Down' in var:
                Integral_down = plot(var,cut)['BkgSum'].Integral()
            else:
                Integral = plot(var,cut)['BkgSum'].Integral()
                syst_up = Integral_up/Integral-1.
                syst_down = Integral/Integral_down-1.
                syst_total = (syst_up+syst_down)/2
                print "up trig_syst for %s:" % cut,syst_up
                print "down trig_syst for %s:" % cut,syst_down
                print "total trig_syst for %s" % cut,syst_total
                syst_trig.append('%s : %.3f' % (cut[:-2],syst_total))
 
    print 'syst_trig = {',syst_trig,'}'
    print 'syst_elec = {',syst_elec,'}'
    print 'syst_muon = {',syst_muon,'}'


def calc_jec():
    sign_mean_up_list = []
    sign_mean_down_list = []
    sign_sigma_up_list = []
    sign_sigma_down_list = []
    sign_mean_list = []
    sign_sigma_list = []
    for cut in ['nnbbSR_nocut','eebbSR_nocut','nnbbVBFSR_nocut','eebbVBFSR_nocut']:
        if 'VBF' in cut:
            sign_list = sign_VBF
        else:
            sign_list = sign
        if 'nn' in cut: 
            var_list = ['X_mass_MET_jesUp','X_mass_MET_jesDown','X_mass_MET_jerUp','X_mass_MET_jerDown','X_mass_MET_nom']
        else:
            var_list = ['X_mass_jesUp','X_mass_jesDown','X_mass_jerUp','X_mass_jerDown','X_mass_nom']
        for var in var_list:
            if 'jesUp' in var:
                vv_mean_up = plot(var,cut)['VV'].GetMean()
                for signal in sign_list:
                    sign_mean_up_list.append(plot(var,cut)[signal].GetMean())
            elif 'jesDown' in var:
                vv_mean_down = plot(var,cut)['VV'].GetMean()
                for signal in sign_list:
                    sign_mean_down_list.append(plot(var,cut)[signal].GetMean())
            elif 'jerUp' in var:
                vv_sigma_up = plot(var,cut)['VV'].GetRMS()
                for signal in sign_list:
                    sign_sigma_up_list.append(plot(var,cut)[signal].GetRMS())
            elif 'jerDown' in var:
                vv_sigma_down = plot(var,cut)['VV'].GetRMS()
                for signal in sign_list:
                    sign_sigma_down_list.append(plot(var,cut)[signal].GetRMS())
            elif 'nom' in var:
                vv_mean= plot(var,cut)['VV'].GetMean()
                vv_sigma = plot(var,cut)['VV'].GetRMS()
                for signal in sign_list:
                    sign_mean_list.append(plot(var,cut)[signal].GetMean())
                    sign_sigma_list.append(plot(var,cut)[signal].GetRMS())
        print '%s VV jes: %.3f,%.3f'% (cut,vv_mean_down/vv_mean,vv_mean_up/vv_mean)
        print '%s VV jer: %.3f,%.3f'% (cut,vv_sigma_down/vv_sigma,vv_sigma_up/vv_sigma)
        print '%s Signal jes:'%cut
        for i in range(len(sign_list)):
            print "%s : [%.3f, %.3f]," %(sign_list[i][5:],sign_mean_down_list[i]/sign_mean_list[i],sign_mean_up_list[i]/sign_mean_list[i])
        print '%s Signal jer:'%cut
        for i in range(len(sign_list)):
            print "%s : [%.3f, %.3f]," %(sign_list[i][5:],sign_sigma_down_list[i]/sign_sigma_list[i],sign_sigma_up_list[i]/sign_sigma_list[i])
    
        

def calc_jmc():
    sign_mean_up_list = []
    sign_mean_down_list = []
    sign_sigma_up_list = []
    sign_sigma_down_list = []
    sign_mean_list = []
    sign_sigma_list = []
    for cut in ['nnbbSR_nocut','eebbSR_nocut','nnbbVBFSR_nocut','eebbVBFSR_nocut']:
        if 'VBF' in cut:
            sign_list = sign_VBF
        else:
            sign_list = sign
        for var in ['H_mass_jmsUp','H_mass_jmsDown','H_mass_jmrUp','H_mass_jmrDown','H_mass_nom']:
            if var=='H_mass_jmsUp':
                vv_mean_up = plot(var,cut)['VV'].GetMean()
                for signal in sign_list:
                    sign_mean_up_list.append(plot(var,cut)[signal].GetMean())
            elif var=='H_mass_jmsDown':
                vv_mean_down = plot(var,cut)['VV'].GetMean()
                for signal in sign_list:
                    sign_mean_down_list.append(plot(var,cut)[signal].GetMean())
            elif var=='H_mass_jmrUp':
                vv_sigma_up = plot(var,cut)['VV'].GetRMS()
                for signal in sign_list:
                    sign_sigma_up_list.append(plot(var,cut)[signal].GetRMS())
            elif var=='H_mass_jmrDown':
                vv_sigma_down = plot(var,cut)['VV'].GetRMS()
                for signal in sign_list:
                    sign_sigma_down_list.append(plot(var,cut)[signal].GetRMS())
            elif var=='H_mass_nom':
                vv_mean= plot(var,cut)['VV'].GetMean()
                vv_sigma = plot(var,cut)['VV'].GetRMS()
                for signal in sign_list:
                    sign_mean_list.append(plot(var,cut)[signal].GetMean())
                    sign_sigma_list.append(plot(var,cut)[signal].GetRMS())
        print '%s VV jms: %.3f,%.3f'% (cut,vv_mean_down/vv_mean,vv_mean_up/vv_mean)
        print '%s VV jmr: %.3f,%.3f'% (cut,vv_sigma_down/vv_sigma,vv_sigma_up/vv_sigma)
        print '%s Signal jms:'%cut
        for i in range(len(sign_list)):
            print "%s : [%.3f, %.3f]," %(sign_list[i][5:],sign_mean_down_list[i]/sign_mean_list[i],sign_mean_up_list[i]/sign_mean_list[i])
        print '%s Signal jmr:'%cut
        for i in range(len(sign_list)):
            print "%s : [%.3f, %.3f]," %(sign_list[i][5:],sign_sigma_down_list[i]/sign_sigma_list[i],sign_sigma_up_list[i]/sign_sigma_list[i])
        

def calc_effb():
    for cut in ['nnbbSR','nn0bSR']:
        Integral_up = []
        Integral_down = []
        Integral = []
        for var in ['BTagAK8Weight_deep_up','BTagAK8Weight_deep_down','BTagAK8Weight_deep']:
            hist = plot(var,cut)
            if 'up' in var:
                for signal in sign:
                    Integral_up.append(hist[signal].Integral())
            elif 'down' in var:
                for signal in sign:
                    Integral_down.append(hist[signal].Integral())
            else:
                for signal in sign:
                    Integral.append(hist[signal].Integral())
        print "Uncertainty for BTagAK8Weight_deep in channel %s" %cut
        for i,value in enumerate(Integral):
            print "%s : [%.3f, %.3f]," %(sign[i][5:],Integral_down[i]/Integral[i],Integral_up[i]/Integral[i])
            
    
    for cut in ['nnbbTR','nn0bTR']:
        Integral_up = 0.
        Integral_down = 0.
        Integral = 0.
        for var in ['BTagAK4Weight_deep_up','BTagAK4Weight_deep_down','BTagAK4Weight_deep']:
            hist = plot(var,cut)
            if 'up' in var:
                Integral_up = hist['BkgSum'].Integral()
            elif 'down' in var:
                Integral_down = hist['BkgSum'].Integral()
            else:
                Integral = hist['BkgSum'].Integral()
        print "Uncertainty for BTagAK4Weight_deep in channel %s" %cut
        print " [%.3f, %.3f]," %(Integral_down/Integral,Integral_up/Integral)
    





if options.syst: calc_syst()
elif options.jec: calc_jec()
elif options.jmc: calc_jmc()
elif options.effb: calc_effb()

