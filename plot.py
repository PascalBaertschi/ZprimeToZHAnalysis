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
#sign_list = ['XZH_M800','XZH_M1000','XZH_M1200','XZH_M1400','XZH_M1600', 'XZH_M1800', 'XZH_M2000', 'XZH_M2500','XZH_M3000','XZH_M3500', 'XZH_M4000','XZH_M4500','XZH_M5000']
sign_list = ['XZH_M1000','XZH_M2000','XZH_M3000']
sign_VBF_list = ['XZHVBF_M1000','XZHVBF_M2000','XZHVBF_M3000']
#sign_VBF_list = ['XZHVBF_M800','XZHVBF_M1000','XZHVBF_M1200','XZHVBF_M1400','XZHVBF_M1600','XZHVBF_M1800','XZHVBF_M2000','XZHVBF_M2500','XZHVBF_M3000','XZHVBF_3500','XZHVBF_M4000','XZHVBF_M4500','XZHVBF_M5000']
#sign_list = []
#
########## ######## ##########

topSF = {
    'nnbb'  : [1.102, 0.065],
    'embb'  : [1.004, 0.036],
    'nn0b'  : [1.230, 0.031],
    'em0b'  : [0.982, 0.009],
    'nnbbVBF'  : [1.102, 0.065],
    'embbVBF'  : [1.004, 0.036],
    'nn0bVBF'  : [1.230, 0.031],
    'em0bVBF'  : [0.982, 0.009],
    'eebb'  : [1.004, 0.036],
    'mmbb'  : [1.004, 0.036],
    'ee0b'  : [0.982, 0.009],
    'mm0b'  : [0.982, 0.009],
    'eebbVBF'  : [1.004, 0.036],
    'mmbbVBF'  : [1.004, 0.036],
    'ee0bVBF'  : [0.982, 0.009],
    'mm0bVBF'  : [0.982, 0.009]
}


########## ######## ##########

jobs = []

def plot(var, cut, norm=False, nm1=False):
    ### Preliminary Operations ###
    treeRead = not 'cutflow' in var
    #treeRead = False
    channel = cut
    TopBTagAK4deep = False
    TopBTagAK4deepup = False
    TopBTagAK4deepdown = False
    if 'inc' in cut:
        eventWeightLuminame = 'eventWeightLumi_nobtag'
    else:
        eventWeightLuminame = 'eventWeightLumi'
    if 'VBF' in cut:
        sign = sign_VBF_list
    else:
        sign = sign_list
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
    stype = "HVT model B"
    if channel.endswith('TR'):
        TopBTagAK4deep = True
        cut = channel
    elif channel.endswith('up'):
        TopBTagAK4deepup = True
        channel = channel[:-2]
        cut = channel
    elif channel.endswith('down'):
        TopBTagAK4deepdown = True
        channel = channel[:-4]
        cut = channel
    
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
    tree = {}
    hist = {}
    ### Create and fill MC histograms ###
    for i, s in enumerate(data+back+sign):
        if treeRead: # Project from tree
            tree[s] = TChain("tree")
            for j, ss in enumerate(sample[s]['files']):
                if not 'data' in s or ('data' in s and ss in pd):
                    tree[s].Add(NTUPLEDIR + ss + ".root")
            if variable[var]['nbins']>0: 
                min_value = variable[var]['min']
                max_value =variable[var]['max']
                title = variable[var]['title']
                if 'isZtoNN' in cut:
                    if var == 'MET':
                        min_value = 200
                        max_value = 2000
                    elif var == 'DPhi':
                        title = "#Delta #varphi (AK8 jet-#slash{E}_{T})"
                    elif var == 'VH_deltaR':
                        title = "#Delta R (#slash{E}_{T}, AK8 jet)"
                hist[s] = TH1F(s, ";"+title+";Events;"+('log' if variable[var]['log'] else ''), variable[var]['nbins'], min_value, max_value)
            else: hist[s] = TH1F(s, ";"+variable[var]['title'], len(variable[var]['bins'])-1, array('f', variable[var]['bins']))
            hist[s].Sumw2()
            cutstring = "%s" % eventWeightLuminame + ("*("+cut+")")
            if TopBTagAK4deep:
                cutstring = "(%s * BTagAK4Weight_deep)" % eventWeightLuminame + ("*("+cut+")")
            elif TopBTagAK4deepup:
                cutstring = "(%s * BTagAK4Weight_deep_up)" % eventWeightLuminame + ("*("+cut+")") 
            elif TopBTagAK4deepdown:
                cutstring = "(%s * BTagAK4Weight_deep_down)" % eventWeightLuminame + ("*("+cut+")") 
            tree[s].Project(s, var, cutstring)
            if not tree[s].GetTree()==None: hist[s].SetOption("%s" % tree[s].GetTree().GetEntriesFast())
            #print tree[s].GetTree().GetEntriesFast()
        else: # Histogram written to file
            for j, ss in enumerate(sample[s]['files']):
                file[ss] = TFile(NTUPLEDIR + ss + ".root", "R")
                if file[ss].IsZombie():
                    print "WARNING: file", NTUPLEDIR + ss + ".root", "does not exist"
                    continue
                tmphist = file[ss].Get(var)
                if tmphist==None: 
                    continue
                if s!="data_obs":
                    if 'VV' in ss: sample_name = ss.replace('VV/','')
                    elif 'ST' in ss: sample_name = ss.replace('ST/','')
                    elif 'TT' in ss: sample_name = ss.replace('TT/','')
                    elif 'WJ' in ss: sample_name = ss.replace('WJ/','')
                    elif 'ZJ' in ss: sample_name = ss.replace('ZJ/','')
                    elif 'DY' in ss: sample_name = ss.replace('DY/','')
                    elif 'XZH' in ss: sample_name = ss.replace('XZH/','')
                    elif 'XZH_VBF' in ss: sample_name = ss.replace('XZH_VBF/','')
                    ref_hist = file[ss].Get('Events')
                    totalEntries = ref_hist.GetBinContent(1)
                    #totalEntries = tmphist.GetBinContent(1)
                    #XS = xsection[sample_name]['xsec']*xsection[sample_name]['kfactor']*xsection[sample_name]['br']
                    #Leq = LUMI*XS/totalEntries if totalEntries > 0 else 0.
                    #tmphist.Scale(Leq)
                if not s in hist.keys(): hist[s] = tmphist
                else: hist[s].Add(tmphist)        
        hist[s].Scale(sample[s]['weight'] if hist[s].Integral() >= 0 else 0)
        hist[s].SetFillColor(sample[s]['fillcolor'])
        hist[s].SetFillStyle(sample[s]['fillstyle'])
        hist[s].SetLineColor(sample[s]['linecolor'])
        hist[s].SetLineStyle(sample[s]['linestyle'])
        
    #if channel.endswith('TR') and channel.replace('TR', '') in topSF:
    #    hist['TTbarSL'].Scale(topSF[channel.replace('TR', '')][0])
    #    hist['ST'].Scale(topSF[channel.replace('TR', '')][0])
    
    hist['BkgSum'] = hist['data_obs'].Clone("BkgSum") if 'data_obs' in hist else hist[back[0]].Clone("BkgSum")
    hist['BkgSum'].Reset("MICES")
    hist['BkgSum'].SetFillStyle(3003)
    hist['BkgSum'].SetFillColor(1)
    for i, s in enumerate(back): 
        hist['BkgSum'].Add(hist[s])
    if options.norm:
        for i, s in enumerate(back + ['BkgSum']): hist[s].Scale(hist[data[0]].Integral()/hist['BkgSum'].Integral())
    # Create data and Bkg sum histograms
#    if options.blind or 'SR' in channel:
#        hist['data_obs'] = hist['BkgSum'].Clone("data_obs")
#        hist['data_obs'].Reset("MICES")
    # Set histogram style
    hist['data_obs'].SetMarkerStyle(20)
    hist['data_obs'].SetMarkerSize(1.25)
    
    for i, s in enumerate(data+back+sign+['BkgSum']): addOverflow(hist[s], False) # Add overflow
    for i, s in enumerate(sign): hist[s].SetLineWidth(3)
    for i, s in enumerate(sign): sample[s]['plot'] = True#sample[s]['plot'] and s.startswith(channel[:2])

    
    # Create stack
    bkg = THStack("Bkg", ";"+hist['BkgSum'].GetXaxis().GetTitle()+";Events")
    for i, s in enumerate(back): bkg.Add(hist[s])

    
    # Legend
    leg = TLegend(0.65, 0.6, 0.95, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0) #1001
    leg.SetFillColor(0)
    if len(data) > 0:
        leg.AddEntry(hist[data[0]], sample[data[0]]['label'], "pe")
    for i, s in reversed(list(enumerate(['BkgSum']+back))):
        leg.AddEntry(hist[s], sample[s]['label'], "f")
    if showSignal:
        for i, s in enumerate(sign):
            if sample[s]['plot']: leg.AddEntry(hist[s], sample[s]['label'], "fl")
        
    leg.SetY1(0.9-leg.GetNRows()*0.05)
    
    
    # --- Display ---
    c1 = TCanvas("c1", hist.values()[0].GetXaxis().GetTitle(), 800, 800 if RATIO else 600)

    if RATIO:
        c1.Divide(1, 2)
        setTopPad(c1.GetPad(1), RATIO)
        setBotPad(c1.GetPad(2), RATIO)
    c1.cd(1)
    c1.GetPad(bool(RATIO)).SetTopMargin(0.06)
    c1.GetPad(bool(RATIO)).SetRightMargin(0.05)
    c1.GetPad(bool(RATIO)).SetTicks(1, 1)
    
    log = "log" in hist['BkgSum'].GetZaxis().GetTitle()
    if log: c1.GetPad(bool(RATIO)).SetLogy()

    # Draw
    bkg.Draw("HIST") # stack
    hist['BkgSum'].Draw("SAME, E2") # sum of bkg
    if not isBlind and len(data) > 0 and not 'Weight' in var: hist['data_obs'].Draw("SAME, PE") # data
    #data_graph.Draw("SAME, PE")
    if showSignal:
        smagn = 1. #if treeRead else 1.e2 #if log else 1.e2
        for i, s in enumerate(sign):
    #        if sample[s]['plot']:
                hist[s].Scale(smagn)
                hist[s].Draw("SAME, HIST") # signals Normalized, hist[s].Integral()*sample[s]['weight']
        textS = drawText(0.80, 0.9-leg.GetNRows()*0.05 - 0.02, stype+" (x%d)" % smagn, True)

    bkg.GetYaxis().SetTitleOffset(bkg.GetYaxis().GetTitleOffset()*1.075)
    bkg.SetMaximum((5. if log else 1.25)*max(bkg.GetMaximum(), hist['data_obs'].GetBinContent(hist['data_obs'].GetMaximumBin())+hist['data_obs'].GetBinError(hist['data_obs'].GetMaximumBin())))
    #if bkg.GetMaximum() < max(hist[sign[0]].GetMaximum(), hist[sign[-1]].GetMaximum()): bkg.SetMaximum(max(hist[sign[0]].GetMaximum(), hist[sign[-1]].GetMaximum())*1.25)
    bkg.SetMinimum(max(min(hist['BkgSum'].GetBinContent(hist['BkgSum'].GetMinimumBin()), hist['data_obs'].GetMinimum()), 5.e-1)  if log else 0.)
    if log:
        bkg.GetYaxis().SetNoExponent(bkg.GetMaximum() < 1.e4)
        bkg.GetYaxis().SetMoreLogLabels(True)
    
    #if log: bkg.SetMinimum(1)
    leg.Draw()
    drawCMS(LUMI,year,"Preliminary")
    drawRegion('XVH'+channel, True)
    drawAnalysis(channel)
    
    #if nm1 and not cutValue is None: drawCut(cutValue, bkg.GetMinimum(), bkg.GetMaximum()) #FIXME
    #if len(sign) > 0:
    #    if channel.startswith('X') and len(sign)>0: drawNorm(0.9-0.05*(leg.GetNRows()+1), "#sigma(X) = %.1f pb" % 1.)
    
    setHistStyle(bkg, 1.2 if RATIO else 1.1)
    setHistStyle(hist['BkgSum'], 1.2 if RATIO else 1.1)
       
    if RATIO:
        c1.cd(2)
        err = hist['BkgSum'].Clone("BkgErr;")
        err.SetTitle("")
        err.GetYaxis().SetTitle("Data / Bkg")
        for i in range(1, err.GetNbinsX()+1):
            err.SetBinContent(i, 1)
            if hist['BkgSum'].GetBinContent(i) > 0:
                err.SetBinError(i, hist['BkgSum'].GetBinError(i)/hist['BkgSum'].GetBinContent(i))
        setBotStyle(err)
        errLine = err.Clone("errLine")
        errLine.SetLineWidth(1)
        errLine.SetFillStyle(0)
        res = hist['data_obs'].Clone("Residues")
        for i in range(0, res.GetNbinsX()+1):
            if hist['BkgSum'].GetBinContent(i) > 0: 
                res.SetBinContent(i, res.GetBinContent(i)/hist['BkgSum'].GetBinContent(i))
                res.SetBinError(i, res.GetBinError(i)/hist['BkgSum'].GetBinContent(i))
        setBotStyle(res)
        if var=="eecutflow_inc":
            err.GetXaxis().SetBinLabel(1,"All")
            err.GetXaxis().SetBinLabel(2,"Z(ee) candidate")
            err.GetXaxis().SetBinLabel(3,"Ele p_{T}")
            err.GetXaxis().SetBinLabel(4,"Ele Id+Iso")
            err.GetXaxis().SetBinLabel(5,"Z boost")
            err.GetXaxis().SetBinLabel(6,"Trigger")
            err.GetXaxis().SetBinLabel(7,"Doubles")
        if var=="mmcutflow_inc":
            err.GetXaxis().SetBinLabel(1,"All")
            err.GetXaxis().SetBinLabel(2,"Z(#mu#mu) candidate")
            err.GetXaxis().SetBinLabel(3,"Muon p_{T}")
            err.GetXaxis().SetBinLabel(4,"Muon Id")
            err.GetXaxis().SetBinLabel(5,"Z boost")
            err.GetXaxis().SetBinLabel(6,"Trigger")
            err.GetXaxis().SetBinLabel(7,"Doubles")
        if var=="nncutflow_inc":
            err.GetXaxis().SetBinLabel(1,"All")
            err.GetXaxis().SetBinLabel(2,"MET > 250")
            err.GetXaxis().SetBinLabel(3,"Lepton veto")
            err.GetXaxis().SetBinLabel(4,"Trigger")
            err.GetXaxis().SetBinLabel(5,"Filter")
        #err.GetXaxis().SetLabelOffset(err.GetXaxis().GetLabelOffset()*5)
        #err.GetXaxis().SetTitleOffset(err.GetXaxis().GetTitleOffset()*2)
        err.Draw("E2")
        errLine.Draw("SAME, HIST")
        if not isBlind and len(data) > 0:
            res.Draw("SAME, PE0")
            #res_graph.Draw("SAME, PE0")
            if len(err.GetXaxis().GetBinLabel(1))==0: # Bin labels: not a ordinary plot
                drawRatio(hist['data_obs'], hist['BkgSum'])
                drawStat(hist['data_obs'], hist['BkgSum'])
    
    c1.Update()
    
    if var=="Events":
        if BLIND:
            removeData = False
            for i in range(hist[back[0]].GetNbinsX()):
                if any(hist[back[0]].GetXaxis().GetBinLabel(i+1) in x for x in ["H mass", "V mass"]): removeData = True
                if removeData:
                    hist['data_obs'].SetBinContent(i+1, -1)
                    res.SetBinContent(i+1, -1.e6)
        if True:
            bg, sn = 'TTbarSL', 'XVH_M2000'
            pres, mass, tag1, tag2 = -1, -1, -1, -1
            for i in range(hist[bg].GetNbinsX()):
                if hist[bg].GetXaxis().GetBinLabel(i+1)=='X mass': pres = i+1
                if hist[bg].GetXaxis().GetBinLabel(i+1)=='H mass': mass = i+1
                if hist[bg].GetXaxis().GetBinLabel(i+1)=='1 b-tag': tag1 = i+1
                if hist[bg].GetXaxis().GetBinLabel(i+1)=='2 b-tag': tag2 = i+1
            Hmis1 = hist[bg].GetBinContent(tag1)/hist[bg].GetBinContent(pres)
            Hmis2 = hist[bg].GetBinContent(tag2)/hist[bg].GetBinContent(pres)
            Heff1 = hist[sn].GetBinContent(tag1)/hist[sn].GetBinContent(pres)
            Heff2 = hist[sn].GetBinContent(tag2)/hist[sn].GetBinContent(pres)
            print "H tagging efficiency 1 b-tag: %.3f, background mistag: %.4f" % (Heff1-Heff2, Hmis1, )
            print "H tagging efficiency 2 b-tag: %.3f, background mistag: %.4f" % (Heff2, Hmis2, )
            
    #    if var=="X_mass" and cut=="XVHqq":
    #        for i, s in enumerate(data+back+['BkgSum']+sign): hist[s].GetXaxis().SetRangeUser(900, 6000)
    #        bkg.GetXaxis().SetRangeUser(900, 6000)
    #        err.GetXaxis().SetRangeUser(900, 6000)
    #        res.GetXaxis().SetRangeUser(900, 6000)
    
    if var=="MET_pt" and cut=="XZHnn":
        bkg.GetXaxis().SetRangeUser(250, 1500)
        err.GetXaxis().SetRangeUser(250, 1500)
        res.GetXaxis().SetRangeUser(250, 1500)
        
    if True and var=="Z_pt":
        ibin, fbin = hist['data_obs'].FindBin(200), hist['data_obs'].GetNbinsX()+1
        ndata, nmc = hist['data_obs'].Integral(ibin, fbin), hist['BkgSum'].Integral(ibin, fbin)
        print "Data/MC [200-Inf]:", ndata, "/", nmc, "=", ndata/nmc
    
    if False:
        ibin, fbin = hist['data_obs'].FindBin(900), hist['data_obs'].GetNbinsX()+1
        ndata, nmc = hist['data_obs'].Integral(ibin, fbin), hist['BkgSum'].Integral(ibin, fbin)
        print "QCD SF:", ndata/nmc
    
        
    #for s in sign:
    #    print s, 100*hist[s].GetRMS()/hist[s].GetMean()
    c1.Update()
    
    if gROOT.IsBatch():
        varname = var.replace('.', '_').replace('()', '')
        if not os.path.exists(OUTPUTDIR+channel_name): os.makedirs(OUTPUTDIR+channel_name)
        c1.Print("%s"% OUTPUTDIR +channel_name+"/"+varname+".png")
        c1.Print("%s"% OUTPUTDIR +channel_name+"/"+varname+".pdf")

    # Print table
    printTable(hist, sign)
    
    # Top scale factors
    if channel.endswith('TR'):
        ndata, nbkg, ntop = hist['data_obs'].Integral(), hist['BkgSum'].Integral(), hist['TTbarSL'].Integral()+hist['ST'].Integral()
        edata, nsub = math.sqrt(ndata), ndata-(nbkg-ntop)
        print "Top scale factor:\t%.3f +- %.3f" % (nsub/ntop, edata/ntop,) 
        if options.top: return [nsub/ntop, edata/ntop]
    
    if len(sign)>0 and 'monoH' in sign[0]:
        if var=='X_mass':
            for s in sign:
                print s,"| Bkg | S/sqrt(B+1)", "\t%-10.2f" % hist[s].Integral(), "\t| %-8.2f" % hist['BkgSum'].Integral(), "\t| %-8.2f" % (hist[s].Integral()/math.sqrt(hist['BkgSum'].Integral()+1), )
        if var=="Events":
            bin = 1
            print "monoH(bb) counts for <", hist[sign[0]].GetXaxis().GetBinLabel(bin), ">"
            for s in sign:
                print s, "\t%-10.2f" % hist[s].GetBinContent(bin), "\t+- %-8.2f" % hist[s].GetBinError(bin)
    #        hRes = res.Clone("hRes")
    #        fRes = TF1("fRes", "pol1", hRes.GetXaxis().GetXmin(), hRes.GetXaxis().GetXmax())
    ##        if 'mnbTR' in channel: hRes.Fit(fRes, "WEQ0", "", 950., 2200)
    ##        elif 'mnbbTR' in channel: hRes.Fit(fRes, "WEQ0", "", 750., 2900)
    ##        elif 'enbTR' in channel: hRes.Fit(fRes, "WEQ0", "", 750., 2600)
    ##        elif 'enbbTR' in channel: hRes.Fit(fRes, "WEQ0", "", 750., 1600)
    #        hRes.Fit(fRes, "CEQ0", "", 750., 3000)
    #        TVirtualFitter.GetFitter().GetConfidenceIntervals(hRes)#, 0.995)
    #        hRes.SetMarkerStyle(0)
    #        hRes.SetFillStyle(3002)
    #        hRes.SetFillColor(2)
    #        hRes.Draw("SAME, E3")
        
        #########
        
        
#        hMC = hist['BkgSum'].Clone("hMC")
#        #if 'QCD' in back: hMC.Add(hist['QCD'], -1)
#        hData = hist['data_obs'].Clone("hData")
#        
#        # Fit
#        c2 = TCanvas("c2", "Top control region", 800, 800)
#        c2.cd()
#        c2.GetPad(0).SetLogy()
#        #hData.Draw("PE")
#        
#        
#        fData = TF1("fData", "[0]*exp([1]*x+[2]/x)", hData.GetXaxis().GetXmin(), hData.GetXaxis().GetXmax())
##        fData.SetParLimits(2, 0., 1.e5)
#        fData.SetLineColor(1)
#        hData.Fit(fData, "EQ0", "") #, 750, 3900 if channel.count('b')>1 else 2500)
#        TVirtualFitter.GetFitter().GetConfidenceIntervals(hData) #, 0.683)
#        hData.SetMarkerStyle(0)
#        hData.SetFillStyle(3002)
#        hData.SetFillColor(1)
#        hData.Draw("CE3")
#        hData.GetXaxis().SetRangeUser(750., 3000)
#        
#        fMC = TF1("fMC", "[0]*exp([1]*x+[2]/x)", hMC.GetXaxis().GetXmin(), hMC.GetXaxis().GetXmax())
##        fMC.SetParLimits(2, 0., 1.e5)
#        fMC.SetLineColor(922)
#        fMC.SetLineStyle(7)
#        #if channel=='enbbTR': hMC.SetBinContent(20, 5)
#        hMC.Fit(fMC, "EQ0", "")#, 750, 4500 if 'm' in channel else 3500) #(4500 if channel.count('b')>1 else 3200))
#        #+("W" if channel.count('b')>1 else "")
##        fMC.SetRange(750, 4500)
##        gMC = TGraphErrors(hist['BkgSum'])
##        gMC.SetFillStyle(3003)
##        gMC.SetFillColor(922)
##        if gMC.GetN() > 0:
#        TVirtualFitter.GetFitter().GetConfidenceIntervals(hMC) #, 0.683)
#        hMC.SetFillStyle(3003)
#        hMC.SetFillColor(2)
#        hMC.Draw("CE3, SAME")
#        
#        hist['data_obs'].Draw("PE, SAME")
#        drawRegion('XVH'+channel, True)
#        if gROOT.IsBatch():
#            c2.Print("plots/"+channel+"/Shape.png")
#            c2.Print("plots/"+channel+"/Shape.pdf")
    if not gROOT.IsBatch(): raw_input("Press Enter to continue...")

    
########## ######## ##########


def plotNorm(var, cut, isPropaganda=False):
    if isPropaganda and not var in ["mass", "tau21", "dbt"]: return
    treeRead = not isPropaganda  #if cut=="" or cut in ["ZtoNN", "WtoEN", "WtoMN", "ZtoEE", "ZtoMM", "VtoQQ"] else True # Read from tree
    channel = cut
    if 'inc' in cut:
        eventWeightLuminame = 'eventWeightLumi_nobtag'
    else:
        eventWeightLuminame = 'eventWeightLumi'
    if "SB" in cut or "SR" in cut:
        channel_name = cut[:-2]
    else:
        channel_name = cut
    plotdir = cut if not treeRead else "tree"
    NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/Ntuples%s/"%(year)
    back = ["VV", "ST", "TTbarSL", "WJetsToLNu_HT", "DYJetsToNuNu_HT", "DYJetsToLL_HT"] 
    sign = ['XZH_M1200','XZH_M3000', 'XZH_M4000']
    colo = {'XZH_M1200' : 856, 'XZH_M3000':856, 'XZH_M4000' : 856}
    line = {'XZH_M1200' : 1, 'XZH_M3000': 2, 'XZH_M4000' : 3}
    
    if treeRead:
        for k in sorted(selection.keys(), key=len, reverse=True):
            if k in cut: cut = cut.replace(k, selection[k])

    file = {}
    tree = {}
    hist = {}
    
    pd = []
    if any(w in cut for w in ['mm', 'isZtoMM']): pd += [x for x in sample['data_obs']['files'] if 'SingleMuon' in x]
    if any(w in cut for w in ['ee', 'isWtoEN', 'isZtoEE', 'emqq', 'isTtoEM']): pd += [x for x in sample['data_obs']['files'] if ('SingleElectron' in x or 'EGamma' in x or 'SinglePhoton' in x)]
    if any(w in cut for w in ['nn', 'isZtoNN']): pd += [x for x in sample['data_obs']['files'] if 'MET' in x]
    if len(pd)==0: raw_input("Warning: Primary Dataset not recognized, continue?")
    
    ### Create and fill MC histograms ###
    for i, s in enumerate(data+back+sign):
        if treeRead: # Project from tree
            tree[s] = TChain("tree")
            for j, ss in enumerate(sample[s]['files']):
                if not 'data' in s or ('data' in s and ss in pd):
                    tree[s].Add(NTUPLEDIR + ss + ".root")
            if variable[var]['nbins']>0: hist[s] = TH1F(s, ";"+variable[var]['title'], variable[var]['nbins'], variable[var]['min'], variable[var]['max'])
            else: hist[s] = TH1F(s, ";"+variable[var]['title'], len(variable[var]['bins'])-1, array('f', variable[var]['bins']))
            hist[s].Sumw2()
            cutstring = "(%s)" % eventWeightLuminame + ("*("+cut+")" if len(cut)>0 else "")
            tree[s].Project(s, var, cutstring)
            hist[s].SetOption("%s" % tree[s].GetTree().GetEntriesFast())
        else: # Histogram written to file
            for j, ss in enumerate(sample[s]['files']):
                if not 'data' in s or ('data' in s and ss in pd):
                    file[ss] = TFile(ntupledir + ss + ".root", "R")
                    if file[ss].IsZombie():
                        print "WARNING: file", ntupledir + ss + ".root", "does not exist"
                        continue
                    varn, vart = var, ""
                    if isPropaganda:
                        if s in data: varn = 'Data_'+var
                        elif s in back: varn = 'Bkg_'+var
                        elif s in sign and 'XVH' in s and 'had' in ss: varn, vart = 'H_'+var, 'H(b#bar{b})'
                        elif s in sign and 'XWH' in s and 'had' in ss: varn, vart = 'V_'+var, 'W(q#bar{q})'
                        elif s in sign and 'XZH' in s and 'had' in ss: varn, vart = 'V_'+var, 'Z(q#bar{q})'
                    tmphist = file[ss].Get(cut+"/"+varn)
                    if tmphist==None: continue
                    tmphist.SetDirectory(0) # !!@@#$#@@@!
                    tmphist.SetName(varn+"_"+ss)
                    tmphist.SetTitle(vart)
                    if not s in hist.keys(): hist[s] = tmphist
                    else: hist[s].Add(tmphist)
                    
        hist[s].Scale(sample[s]['weight'] if hist[s].Integral() >= 0 else 0)
#        hist[s].SetFillColor(sample[s]['fillcolor'])
#        hist[s].SetFillStyle(sample[s]['fillstyle'])
#        hist[s].SetLineColor(sample[s]['linecolor'])
#        hist[s].SetLineStyle(sample[s]['linestyle'])
        
        if 'tau21' in var:
            hist[s].Rebin(4)
        if 'dbt' in var:
            hist[s].Rebin(2)
            if s[1]=='W' or s[1]=='Z': hist[s].Smooth()
        if 'mass' in var: 
            hist[s].Rebin(2)
    
    
    hist['BkgSum'] = hist['data_obs'].Clone("BkgSum") if 'data_obs' in hist else hist[back[0]].Clone("BkgSum")
    hist['BkgSum'].Reset("MICES")
#    hist['BkgSum'].SetFillStyle(-1)
#    hist['BkgSum'].SetFillColor(1)
#    hist['BkgSum'].SetLineStyle(-1)
#    hist['BkgSum'].SetMarkerStyle(20)
    hist['BkgSum'].SetLineColor(922)
    hist['BkgSum'].SetLineWidth(3)
    hist['BkgSum'].SetFillColor(920)
    hist['BkgSum'].SetFillStyle(3005)
    for i, s in enumerate(back): hist['BkgSum'].Add(hist[s])
    hist['BkgSum'].SetMaximum(1.25*max(hist['BkgSum'].GetMaximum(), hist[sign[-1]].GetMaximum()))
    hist['BkgSum'].Scale(hist['data_obs'].Integral()/hist['BkgSum'].Integral())
    
    hist['data_obs'].SetMarkerStyle(20)
    hist['data_obs'].SetMarkerSize(1.)
    hist['data_obs'].SetLineColor(1)
    
    maxY = hist['BkgSum'].GetMaximum()
    for i, s in enumerate(sign):
        hist[s].Scale(hist['BkgSum'].Integral()/hist[s].Integral() * hist[s].GetNbinsX()/hist['BkgSum'].GetNbinsX())
        hist[s].SetLineColor(colo[s])
        hist[s].SetLineStyle(line[s])
        hist[s].SetLineWidth(3)
        if hist[s].GetMaximum() > maxY: maxY =  hist[s].GetMaximum()
    
    if False:#BLIND:
        if 'tau21' in var:
            for i, s in enumerate(data):
                first, last = hist[s].FindBin(0), hist[s].FindBin(0.75)
                for j in range(first, last): hist[s].SetBinContent(j, -1.e-4)
            hist['BkgSum'].SetMaximum(2*hist['BkgSum'].GetMaximum())
        if 'mass' in var:
            for i, s in enumerate(data):
                first, last = hist[s].FindBin(65), hist[s].FindBin(135)
                for j in range(first, last): hist[s].SetBinContent(j, -1.e-4)
    
    lh = 0.045
    
    # Legend
    leg1 = TLegend(0.35, 0.55, 1.15, 0.90)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0) #1001
    leg1.SetFillColor(0)
    leg1.SetNColumns(2)
    leg1.AddEntry(hist['data_obs'], 'Data', "pe")
    leg1.AddEntry(hist['BkgSum'], 'Background simulation', "fl")
    leg1.SetY1(0.90-leg1.GetNRows()*lh)
    
    leg2 = TLegend(0.375, 0.55, 0.975, 0.90-leg1.GetNRows()*lh)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0) #1001
    leg2.SetFillColor(0)
    leg2.SetNColumns(2)
    for i, s in enumerate(sign):
        leg2.AddEntry(hist[s], sample[s]['label'].replace(" = ", "="), "l")
    leg2.SetY1(0.90-(leg1.GetNRows()+leg2.GetNRows())*lh)
    
    leg3 = TLegend(0.25, 0.55, 0.35, 0.90-leg1.GetNRows()*lh)
    leg3.SetBorderSize(0)
    leg3.SetFillStyle(0) #1001
    leg3.SetFillColor(0)
    for i, s in enumerate(sign):
        if i%2==0: leg3.AddEntry(None, hist[s].GetTitle(), "")
    leg3.SetY1(0.90-(leg1.GetNRows()+leg2.GetNRows())*lh)
    
    log = "log" in hist['BkgSum'].GetZaxis().GetTitle()
    
    # --- Display ---
    c1 = TCanvas("c1", hist.values()[0].GetXaxis().GetTitle(), 800, 600)
    c1.cd()
    c1.GetPad(0).SetTopMargin(0.06)
    c1.GetPad(0).SetRightMargin(0.05)
    c1.GetPad(0).SetTicks(1, 1)
    if log: c1.GetPad(0).SetLogy()
    setHistStyle(hist['BkgSum'], 1.1)
    
    # Magic trick: https://root-forum.cern.ch/t/changing-position-of-the-y-axis-exponent/12025/3
    from ROOT import TGaxis
    #TGaxis.SetMaxDigits(2)
    TGaxis.SetExponentOffset(-0.05, 0., 'y')
    
    # Draw
    hist['BkgSum'].Draw("HIST, C")
    for i, s in enumerate(sign):
        hist[s].Draw("SAME, HIST, C") # signals Normalized, hist[s].Integral()*sample[s]['weight']
#    bkg.GetYaxis().SetTitleOffset(bkg.GetYaxis().GetTitleOffset()*1.075)
#    bkg.SetMaximum((5. if log else 1.25)*max(bkg.GetMaximum(), hist['data_obs'].GetBinContent(hist['data_obs'].GetMaximumBin())+hist['data_obs'].GetBinError(hist['data_obs'].GetMaximumBin())))
#    if bkg.GetMaximum() < max(hist[sign[0]].GetMaximum(), hist[sign[-1]].GetMaximum()): bkg.SetMaximum(max(hist[sign[0]].GetMaximum(), hist[sign[-1]].GetMaximum())*1.25)
#    bkg.SetMinimum(max(min(hist['BkgSum'].GetBinContent(hist['BkgSum'].GetMinimumBin()), hist['data_obs'].GetMinimum()), 5.e-1)  if log else 0.)
#    if log:
#        bkg.GetYaxis().SetNoExponent(bkg.GetMaximum() < 1.e4)
#        bkg.GetYaxis().SetMoreLogLabels(True)
    hist['data_obs'].Draw("SAME, PEX0")
    #if log: bkg.SetMinimum(1)
    leg1.Draw()
    leg2.Draw()
    leg3.Draw()
    drawCMS(LUMI,year,"Preliminary") #Preliminary
    #drawRegion(channel, True)
    drawAnalysis(channel)
    
    if 'tau21' in var:
        hist['BkgSum'].GetXaxis().SetTitle("#it{N}-subjettiness  #tau_{21}")
        hist['BkgSum'].GetXaxis().SetRangeUser(0, 1)
        hist['BkgSum'].SetMaximum(maxY*1.6)
        lineT = drawLine(0.35, 0., 0.35, 0.7*hist['BkgSum'].GetMaximum())
        lineL = drawLine(0.75, 0., 0.75, 0.7*hist['BkgSum'].GetMaximum())
        textT = drawText(0.2, 0.7*hist['BkgSum'].GetMaximum(), "high purity")
        textL = drawText(0.55, 0.7*hist['BkgSum'].GetMaximum(), "low purity")
        #textA = drawText(0.85, 0.6*hist['BkgSum'].GetMaximum(), "anti-tag")
        for s in sign+['BkgSum']:
            print s, "high purity: %.2f, low purity: %.2f " % (hist[s].Integral(0, hist[s].FindBin(0.35)-1)/hist[s].Integral(), hist[s].Integral(hist[s].FindBin(0.35), hist[s].FindBin(0.75)-1)/hist[s].Integral())
            #for s in sign+['BkgSum']: for i in range(1, hist[s].GetNbinsX()): print i, hist[s].GetBinContent(i), hist[s].GetXaxis().GetBinUpEdge(i), hist[s].Integral(0, i)/hist[s].Integral()
    if 'dbt' in var:
        hist['BkgSum'].GetXaxis().SetTitle("b tagging discriminator")
        hist['BkgSum'].GetXaxis().SetRangeUser(-1.0, 1.10)
        hist['BkgSum'].SetMaximum(maxY*1.6)
        lineT = drawLine(0.9, 0., 0.9, 0.7*hist['BkgSum'].GetMaximum())
        lineL = drawLine(0.3, 0., 0.3, 0.7*hist['BkgSum'].GetMaximum())
        textT = drawText(0.99, 0.7*hist['BkgSum'].GetMaximum(), "tight")
        textL = drawText(0.6, 0.7*hist['BkgSum'].GetMaximum(), "loose")
        for s in sign+['BkgSum']:
            print s, "high purity: %.2f, low purity: %.2f, total: %.2f" % (hist[s].Integral(hist[s].FindBin(0.9)-1, hist[s].GetNbinsX()+1)/hist[s].Integral(), hist[s].Integral(hist[s].FindBin(0.3)-1, hist[s].FindBin(0.9))/hist[s].Integral(), hist[s].Integral(hist[s].FindBin(0.3)-1, hist[s].GetNbinsX()+1)/hist[s].Integral())
            #for i in range(1, hist[s].GetNbinsX()): print i, hist[s].GetBinContent(i), hist[s].GetXaxis().GetBinUpEdge(i), hist[s].Integral(0, i)/hist[s].Integral()
    if 'mass' in var:
        hist['BkgSum'].GetXaxis().SetTitle("Soft-drop PUPPI jet mass (GeV)")
        hist['BkgSum'].GetXaxis().SetRangeUser(0, 200.)
        hist['BkgSum'].SetMaximum(maxY*1.3)
        lineM = drawLine(85, 0., 85, 0.7*hist['BkgSum'].GetMaximum())
        lineZ = drawLine(105, 0., 105, 0.7*hist['BkgSum'].GetMaximum())
        lineH = drawLine(135, 0., 135, 0.7*hist['BkgSum'].GetMaximum())
        textL = drawText(95, 0.7*hist['BkgSum'].GetMaximum(), "Z")
        textA = drawText(120, 0.7*hist['BkgSum'].GetMaximum(), "H")
        #print hist['XVH_M1200'].Integral(hist['XVH_M1200'].FindBin(105), hist['XVH_M1200'].FindBin(135))/hist['XVH_M1200'].Integral()
    
    c1.Update()
    if gROOT.IsBatch():
        c1.Print("plots/Norm/"+var+".png")
        c1.Print("plots/Norm/"+var+".pdf")
#    setHistStyle(hist['BkgSum'], 1.2 if RATIO else 1.1)
    
    if not gROOT.IsBatch(): raw_input("Press Enter to continue...")



########## ######## ##########  
def plotTop():
    gROOT.SetBatch(True)
    #back.remove('QCD')
    CR = ['nnbb','embb','nn0b','em0b']
    SF = {}
    for tcr in CR:
        for v in ['H_mass', 'X_mass','H_pt', 'nJets', 'MinDPhi']:
            SF[tcr] = plot(v, tcr+'TR')
        
    print "topSF = {"
    for tcr in CR:
        print '    %-7s : [%.3f, %.3f],' % ("'"+tcr+"'", SF[tcr][0], SF[tcr][1])
    for tcr in ['eebb', 'mmbb','ee0b','mm0b']:
        region = tcr.replace('ee', 'em').replace('mm', 'em')
        print "    %-7s : [%.3f, %.3f]," % ("'"+tcr+"'", SF[region][0], SF[region][1])
    print "}"

def plotTopError():
    gROOT.SetBatch(True)
    #back.remove('QCD')
    CR = ['nnbb','embb','nn0b','em0b']
    SF = {}
    SF_up = {}
    SF_down = {}
    syst_trig = {'nnbb' : 0.000, 'eebb' : 0.004, 'mmbb' : 0.018, 'nn0b' : 0.000, 'ee0b' : 0.004, 'mm0b' : 0.016, 'nnbbVBF' : 0.000, 'eebbVBF' : 0.005, 'mmbbVBF' : 0.016, 'nn0bVBF' : 0.000, 'ee0bVBF' : 0.005, 'mm0bVBF' : 0.018}
    syst_elec = {'nnbb' : 0.000, 'eebb' : 0.067, 'mmbb' : 0.000, 'nn0b' : 0.000, 'ee0b' : 0.089, 'mm0b' : 0.000, 'nnbbVBF' : 0.000, 'eebbVBF' : 0.096, 'mmbbVBF' : 0.000, 'nn0bVBF' : 0.000, 'ee0bVBF' : 0.099, 'mm0bVBF' : 0.000}
    syst_muon = {'nnbb' : 0.000, 'eebb' : 0.000, 'mmbb' : 0.005, 'nn0b' : 0.000, 'ee0b' : 0.000, 'mm0b' : 0.005, 'nnbbVBF' : 0.000, 'eebbVBF' : 0.000, 'mmbbVBF' : 0.006, 'nn0bVBF' : 0.000, 'ee0bVBF' : 0.000, 'mm0bVBF' : 0.006}
    for tcr in CR:
        for v in ['H_mass', 'X_mass','H_pt', 'nJets', 'MinDPhi']:
            SF[tcr] = plot(v, tcr+'TR')
            SF_up[tcr] = plot(v, tcr+'TRup')
            SF_down[tcr] = plot(v, tcr+'TRdown')
    
    print "topSF = {"
    for tcr in CR:
        print '    %-7s : [%.3f, %.3f],' % ("'"+tcr+"'", SF[tcr][0], SF[tcr][1])
    for tcr in ['eebb', 'mmbb','ee0b','mm0b']: 
        region = tcr.replace('ee', 'em').replace('mm', 'em')
        print "    %-7s : [%.3f, %.3f]," % ("'"+tcr+"'", SF[region][0], SF[region][1])
    print "}"
    print "topSFerror = {"
    for tcr in CR:
        #print '    %-7s : [%.3f, %.3f],' % ("'"+tcr+'up'"'", SF_up[tcr][0], SF_up[tcr][1])
        #print '    %-7s : [%.3f, %.3f],' % ("'"+tcr+'down'"'", SF_down[tcr][0], SF_down[tcr][1])
        print tcr," systematic error:",round(abs(SF_up[tcr][0]-SF_down[tcr][0]),3)
    for tcr in ['eebb','mmbb','ee0b','mm0b']: 
        region = tcr.replace('ee', 'em').replace('mm', 'em')
        #print "    %-7s : [%.3f, %.3f]," % ("'"+tcr+'up'"'", SF_up[region][0], SF_up[region][1])
        #print "    %-7s : [%.3f, %.3f]," % ("'"+tcr+'down'"'", SF_down[region][0], SF_down[region][1])
        if 'ee' in tcr:
            tcr_mm = tcr.replace('ee', 'mm')
            print tcr," systematic error:",round(np.sqrt(syst_elec[tcr]**2+syst_muon[tcr_mm]**2),3)
        elif 'mm' in tcr:
            tcr_ee = tcr.replace('mm','ee')
            print tcr," systematic error:",round(np.sqrt((syst_elec[tcr_ee]+syst_trig[tcr_ee])**2+(syst_muon[tcr]+syst_trig[tcr])**2),3)
    print "}"

def checktagger(tagger, cut):
    c1 = TCanvas("c1", "Jetmass", 800, 600)
    c1.GetPad(0).SetTopMargin(0.06)
    c1.GetPad(0).SetRightMargin(0.06)
    c1.GetPad(0).SetTicky(2)
    gStyle.SetOptStat(0)
    #gStyle.SetOptTitle(0)
    background = back
    file = {}
    var = 'H_mass'
    channel = cut
    if "SB" in cut or "SR" in cut:
        channel_name = cut[:-2]
    else:
        channel_name = cut
    if 'inc' in cut:
        eventWeightLuminame = 'eventWeightLumi_nobtag'
    else:
        eventWeightLuminame = 'eventWeightLumi'
    NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/Ntuples%s/"%(year)
    if tagger=='H_dbt':
        wp_list = [0,0.6,0.7,0.8]
    elif tagger=='H_deepcsv2':
        wp_list = [0,0.05,0.1,0.2]
    elif tagger=='DeepTagMD_ZHbbvsQCD':
        wp_list = [0,0.4,0.6,0.8]
    color_list = [1,2,3,4,7,]
    leg = TLegend(0.7, 0.55, 0.9, 0.8)
    leg.SetBorderSize(0)
    for wp_order,wp in enumerate(wp_list):
        #for wp in [0.4,0.5,0.6,0.7]:
        if wp==0:
            cut_mod = "%s" %selection[cut]
        else:
            cut_tag = " && %s > %s" %(tagger,wp)
            cut_mod = "%s%s" %(selection[cut],cut_tag)
        varname = var.replace('.', '_').replace('()', '')
        tree = {}
        hist = {}
        nbins = variable[var]['nbins']*100
    
        ### Create and fill MC histograms ###
        for i, s in enumerate(background+sign):
            tree[s] = TChain("tree")
            for j, ss in enumerate(sample[s]['files']):
                tree[s].Add(NTUPLEDIR + ss + ".root")
            hist[s] = TH1F(s, ";"+variable[var]['title']+";Events;"+('log' if variable[var]['log'] else ''), nbins, variable[var]['min'], variable[var]['max'])
            hist[s].Sumw2()
            tree[s].Project(s, var, "(%s)*(" % eventWeightLuminame +cut_mod+")")
        BkgSum = hist[background[0]].Clone("BkgSum")
        BkgSum.Reset("MICES")
        for i, s in enumerate(background): BkgSum.Add(hist[s], 1)
        #SigSum = hist[sign[0]].Clone("SigSum")
        #SigSum.Reset("MICES")
        #for i, s in enumerate(sign): SigSum.Add(hist[s], 1)
        BkgSum_rebined = BkgSum.Clone("BkgSum_rebined")
        #SigSum_rebined = SigSum.Clone("SigSum_rebined")
        BkgSum_rebined.Rebin(200)
        #SigSum_rebined.Rebin(200)
        BkgSum_rebined.Scale(1./BkgSum_rebined.Integral())
        #SigSum_rebined.Scale(1./SigSum_rebined.Integral())
        BkgSum_rebined.SetTitle("%s" %tagger)
        if wp_order==0:
            BkgSum_rebined.SetLineColor(color_list[wp_order])
            BkgSum_rebined.Draw("HIST")
            leg.AddEntry(BkgSum_rebined,'no tagger','L')
        elif wp_order==len(wp_list)-1:
            BkgSum_rebined.SetLineColor(color_list[wp_order])
            BkgSum_rebined.Draw("SAME HIST")
            leg.AddEntry(BkgSum_rebined,'tagger>%s'%wp,'L')
            leg.Draw("SAME")
        else:
            BkgSum_rebined.SetLineColor(color_list[wp_order])
            BkgSum_rebined.Draw("SAME HIST")
            leg.AddEntry(BkgSum_rebined,'tagger>%s'%wp,'L')


    c1.Print("%s%s/%s_masssculpting.png"%(OUTPUTDIR,channel_name,tagger))
    c1.Print("%s%s/%s_masssculpting.pdf"%(OUTPUTDIR,channel_name,tagger))
     
def getwp(tagger, cut):
    sign_list = ['XZH_M1000','XZH_M2000']
    c1 = TCanvas("c1", "X_mass", 800, 600)
    c1.GetPad(0).SetTopMargin(0.06)
    c1.GetPad(0).SetRightMargin(0.06)
    c1.GetPad(0).SetTicky(2)
    gStyle.SetOptStat(0)
    #gStyle.SetOptTitle(0)
    leg = TLegend(0.7, 0.55, 0.9, 0.8)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    background = back
    channel = cut
    if "SB" in cut or "SR" in cut:
        channel_name = cut[:-2]
    else:
        channel_name = cut
    if 'inc' in cut:
        eventWeightLuminame = 'eventWeightLumi_nobtag'
    else:
        eventWeightLuminame = 'eventWeightLumi'
    NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/Ntuples%s/"%(year)
    file = {}
    var = 'X_mass'
    for sign in sign_list:
        wp_list = []
        fom_list = []
        for wp in np.arange(0.1,1.0,0.05):
            wp_list.append(wp)
            cut_tag = " && %s > %s" %(tagger,wp)
            cut_mod = "%s%s" %(selection[cut],cut_tag)
            varname = var.replace('.', '_').replace('()', '')
            tree = {}
            hist = {}
            nbins = variable[var]['nbins']*100
    
        ### Create and fill MC histograms ###
            for i, s in enumerate(background+[sign]):
                tree[s] = TChain("tree")
                for j, ss in enumerate(sample[s]['files']):
                    tree[s].Add(NTUPLEDIR + ss + ".root")
                hist[s] = TH1F(s, ";"+variable[var]['title']+";Events;"+('log' if variable[var]['log'] else ''), nbins, variable[var]['min'], variable[var]['max'])
                hist[s].Sumw2()
                tree[s].Project(s, var, "(%s)*(" % eventWeightLuminame +cut_mod+")")
            BkgSum = hist[background[0]].Clone("BkgSum")
            BkgSum.Reset("MICES")
            for i, s in enumerate(background): BkgSum.Add(hist[s], 1)
            SigSum = hist[sign].Clone("SigSum")
            SigSum.Reset("MICES")
            for i, s in enumerate([sign]): SigSum.Add(hist[s], 1)
            BkgSum_rebined = BkgSum.Clone("BkgSum_rebined")
            SigSum_rebined = SigSum.Clone("SigSum_rebined")
            BkgSum_rebined.Rebin(50)
            SigSum_rebined.Rebin(50)
            BkgSum_rebined.SetTitle("%s" %tagger)
            BkgSum_rebined.SetLineColor(4)
            SigSum_rebined.SetLineColor(2)
            #BkgSum_rebined.Draw("HIST")
            #SigSum_rebined.Draw("SAME HIST")
            nbins = BkgSum_rebined.GetNbinsX()
            max_value = SigSum_rebined.GetMaximum()
            binmax = SigSum_rebined.GetMaximumBin()
            for i in range(20):
                if SigSum_rebined.GetBinContent(binmax+i) > max_value/2.:
                    continue
                else:
                    right_border = binmax+i
                    break
            for i in range(20):
                if SigSum_rebined.GetBinContent(binmax-i) > max_value/2.:
                    continue
                else:
                    left_border = binmax-i
                    break
            S = SigSum_rebined.Integral(left_border,right_border)
            B_left = BkgSum.Integral(0, left_border)
            B_right = BkgSum.Integral(right_border,nbins)
            B = B_left + B_right
            fom = S/(np.sqrt(B)+1)
            if '2000' in sign:
                fom_list.append(fom*10)
            else:
                fom_list.append(fom)
            #print "fom:",fom,"for tagger:",tagger,"and wp:",wp
            x_left_border = SigSum_rebined.GetXaxis().GetBinCenter(left_border)
            x_right_border = SigSum_rebined.GetXaxis().GetBinCenter(right_border)
            line_left_border = TLine(x_left_border,0.,x_left_border,130.)
            line_right_border = TLine(x_right_border,0.,x_right_border,130.)
            #line_left_border.Draw("SAME")
            #line_right_border.Draw("SAME")
    
            #c1.Print("%s%s/%s_wp%s.png"%(OUTPUTDIR,channel,tagger,wp))
            #c1.Print("%s%s/%s_wp%s.pdf"%(OUTPUTDIR,channel,tagger,wp))
        if '1000' in sign:
            line = TGraph(len(wp_list),np.array(wp_list),np.array(fom_list))
            line.SetMarkerStyle(2)
            line.SetMarkerColor(2)
            line.SetTitle("FoM for different wp")
            line.GetXaxis().SetTitle("wp")
            line.GetYaxis().SetTitle("S/(sqrt(B)+1)")
            line.GetYaxis().SetRangeUser(0.,16.)
            leg.AddEntry(line,'1 TeV signal',"p")
            line.Draw("PA")
        else:
            line1 = TGraph(len(wp_list),np.array(wp_list),np.array(fom_list))
            line1.SetMarkerStyle(2)
            line1.SetMarkerColor(4)
            leg.AddEntry(line1,'2 TeV signal #times 10',"p")
            line1.Draw("P SAME")
            leg.Draw("SAME")
        
    c1.Print("%s%s/%s_wp.png"%(OUTPUTDIR,channel_name,tagger))
    c1.Print("%s%s/%s_wp.pdf"%(OUTPUTDIR,channel_name,tagger))
            
def plotROC_combined(var_list, cut,cut_sign):
    c1 = TCanvas("c1", "Signals", 800, 600)
    c1.GetPad(0).SetTopMargin(0.06)
    c1.GetPad(0).SetRightMargin(0.06)
    c1.GetPad(0).SetTicky(2)
    gStyle.SetOptStat(0)
    #gStyle.SetOptTitle(0)
    c1.GetPad(0).SetLogy()
    #leg = TLegend(0.6, 0.15, 0.9, 0.5)
    leg = TLegend(0.6, 0.15, 0.9, 0.4)
    #leg = TLegend(0.52, 0.15, 0.92, 0.4)
    leg.SetBorderSize(0)
    channel = cut
    if "SB" in cut or "SR" in cut:
        channel_name = cut[:-2]
    else:
        channel_name = cut
    if 'inc' in cut:
        eventWeightLuminame = 'eventWeightLumi_nobtag'
    else:
        eventWeightLuminame = 'eventWeightLumi'
    NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/Ntuples%s/"%(year)
    background = back
    file = {}
    H_csv1 = []
    H_deepcsv1 = []
    for var in var_list:
        varname = var.replace('.', '_').replace('()', '')
        wp = [0.5, ]
        tree = {}
        hist = {}
        nbins = variable[var]['nbins']*100
        ### Create and fill MC histograms ###
        for i, s in enumerate(background+sign):
            tree[s] = TChain("tree")
            for j, ss in enumerate(sample[s]['files']):
                tree[s].Add(NTUPLEDIR + ss + ".root")
            hist[s] = TH1F(s, ";"+variable[var]['title']+";Events;"+('log' if variable[var]['log'] else ''), nbins, variable[var]['min'], variable[var]['max'])
            hist[s].Sumw2()
            if 'XZH_M' in s:
                tree[s].Project(s, var, "(%s)*(" % eventWeightLuminame +selection[cut_sign]+")")
            else:
                tree[s].Project(s, var, "(%s)*(" % eventWeightLumiName +selection[cut]+")")
        BkgSum = hist[background[0]].Clone("BkgSum")
        BkgSum.Reset("MICES")
        for i, s in enumerate(background): BkgSum.Add(hist[s], 1)
        SigSum = hist[sign[0]].Clone("SigSum")
        SigSum.Reset("MICES")
        for i, s in enumerate(sign): SigSum.Add(hist[s], 1)

        if var == 'H_tau21':
            g = TGraph()
            g.SetLineWidth(3)
            g.SetTitle("tagger comparison;signal efficiency;background efficiency")
            g.SetLineColor(629)
            g.GetXaxis().SetRangeUser(0., 1.)
            g.GetYaxis().SetRangeUser(1.e-3, 1.)
            leg.AddEntry(g,'tau21','L')
            #print "bkgSum integral:",BkgSum.Integral(1,nbins)
            #print "sigSum integral:",SigSum.Integral(1,nbins)
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                g.SetPoint(g.GetN(), s, b)
            g.Draw("AL")
            drawCMS(-1, year, "Simulation", False)
        elif var == 'H_tau31':
            g1 = TGraph()
            g1.SetLineWidth(3)
            g1.SetLineColor(1)
            g1.SetLineStyle(7)
            leg.AddEntry(g1,'tau31','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                g1.SetPoint(g1.GetN(), s, b)
            g1.Draw("SAME L") 
        elif var == 'H_tau32':
            g2 = TGraph()
            g2.SetLineWidth(3)
            g2.SetLineColor(3)
            leg.AddEntry(g2,'tau32','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                g2.SetPoint(g2.GetN(), s, b)
            g2.Draw("SAME L") 
        elif var == 'H_tau41':
            g3 = TGraph()
            g3.SetLineWidth(3)
            g3.SetLineColor(4)
            g3.SetLineStyle(7)
            leg.AddEntry(g3,'tau41','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                g3.SetPoint(g3.GetN(), s, b)
            g3.Draw("SAME L") 
        elif var == 'H_tau42':
            g4 = TGraph()
            g4.SetLineWidth(3)
            g4.SetLineColor(5)
            leg.AddEntry(g4,'tau42','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                g4.SetPoint(g4.GetN(), s, b)
            g4.Draw("SAME L") 
        elif var == 'H_ddt':
            g5 = TGraph()
            g5.SetLineWidth(3)
            g5.SetLineColor(6)
            g5.SetLineStyle(7)
            leg.AddEntry(g5,'tau21_ddt','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                g5.SetPoint(g5.GetN(), s, b)
            g5.Draw("SAME L") 
        elif var == 'DeepTagMD_H4qvsQCD':
            g6 = TGraph()
            g6.SetLineWidth(3)
            g6.SetLineColor(7)
            leg.AddEntry(g6,'DeepTagMD_H4qvsQCD','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                g6.SetPoint(g6.GetN(), s, b)
            g6.Draw("SAME L") 
        elif var == 'DeepTagMD_WvsQCD':
            g7 = TGraph()
            g7.SetLineWidth(3)
            g7.SetLineColor(8)
            g7.SetLineStyle(7)
            leg.AddEntry(g7,'DeepTagMD_WvsQCD','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                g7.SetPoint(g7.GetN(), s, b)
            g7.Draw("SAME L")
        elif var == 'DeepTagMD_ZvsQCD':
            g8 = TGraph()
            g8.SetLineWidth(3)
            g8.SetLineColor(9)
            leg.AddEntry(g8,'DeepTagMD_ZvsQCD','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                g8.SetPoint(g8.GetN(), s, b)
            g8.Draw("SAME L")
            leg.Draw("SAME")
        elif var == 'H_csv1':
            g = TGraph()
            g.SetLineWidth(3)
            g.SetLineColor(44)
            g.SetLineStyle(7)
            g.SetTitle("b-tagger comparison;signal efficiency;background efficiency")
            g.GetXaxis().SetRangeUser(0., 1.)
            g.GetYaxis().SetRangeUser(1.e-3, 1.)
            leg.AddEntry(g,'csv1','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                H_csv1.append([b,s])
                g.SetPoint(g.GetN(), s, b)
            g.Draw("AL")
            drawCMS(-1, year, "Simulation", False)
        elif var == 'H_csv2':
            g1 = TGraph()
            g1.SetLineWidth(3)
            g1.SetLineColor(629)
            leg.AddEntry(g1,'csv','L')
            g2 = TGraph()
            g2.SetLineWidth(3)
            g2.SetLineColor(1)
            g2.SetLineStyle(7)
            leg.AddEntry(g2,'csv2','L')
            for i in range(1,nbins+1):
                b_1 = H_csv1[i-1][0]
                s_1 = H_csv1[i-1][1]
                b_2 = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s_2 = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                s = (s_1+s_2)/2
                b = (b_1+b_2)/2
                g1.SetPoint(g1.GetN(), s, b)
                g2.SetPoint(g2.GetN(), s_2, b_2)
            g1.Draw("SAME L")
            g2.Draw("SAME L")
        elif var == 'H_deepcsv1':
            g3 = TGraph()
            g3.SetLineWidth(3)
            g3.SetLineColor(2)
            g3.SetLineStyle(7)
            leg.AddEntry(g3,'deepcsv1','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                H_deepcsv1.append([b,s])
                g3.SetPoint(g3.GetN(),s,b)
            g3.Draw("SAME L")
        elif var == 'H_deepcsv2':
            g4 = TGraph()
            g4.SetLineWidth(3)
            g4.SetLineColor(3)
            g4.SetLineStyle(7)
            leg.AddEntry(g4,'deepcsv2','L')
            g5 = TGraph()
            g5.SetLineWidth(3)
            g5.SetLineColor(12)
            leg.AddEntry(g5, 'deepcsv', 'L')
            for i in range(1,nbins+1):
                b_1 = H_deepcsv1[i-1][0]
                s_1 = H_deepcsv1[i-1][1]
                b_2 = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s_2 = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                s = (s_1+s_2)/2
                b = (b_1+b_2)/2
                g4.SetPoint(g4.GetN(), s_2, b_2)
                g5.SetPoint(g5.GetN(), s, b)
            g4.Draw("SAME L")  
            g5.Draw("SAME L")
        elif var == 'BtagDeepB':
            g6 = TGraph()
            g6.SetLineWidth(3)
            g6.SetLineColor(5)
            leg.AddEntry(g6,'BtagDeepB','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                g6.SetPoint(g6.GetN(), s, b)
            g6.Draw("SAME L")
        elif var == 'H_dbt':
            g7 = TGraph()
            g7.SetLineWidth(3)
            g7.SetLineColor(6)
            g7.SetLineStyle(7)
            leg.AddEntry(g7,'DoubleB','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                g7.SetPoint(g7.GetN(), s, b)
            g7.Draw("SAME L")
        elif var == 'DeepTagMD_HbbvsQCD':
            g8 = TGraph()
            g8.SetLineWidth(3)
            g8.SetLineColor(7)
            leg.AddEntry(g8,'DeepTagMD_HbbvsQCD','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                g8.SetPoint(g8.GetN(), s, b)
            g8.Draw("SAME L")
        elif var == 'DeepTagMD_ZHbbvsQCD':
            g9 = TGraph()
            g9.SetLineWidth(3)
            g9.SetLineColor(8)
            leg.AddEntry(g9,'DeepTagMD_ZHbbvsQCD','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                g9.SetPoint(g9.GetN(), s, b)
            g9.Draw("SAME L")
        elif var == 'DeepTagMD_ZbbvsQCD':
            g10 = TGraph()
            g10.SetLineWidth(3)
            g10.SetLineColor(9)
            leg.AddEntry(g10,'DeepTagMD_ZbbvsQCD','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                g10.SetPoint(g10.GetN(), s, b)
            g10.Draw("SAME L")
        elif var == 'DeepTagMD_bbvsLight':
            g11 = TGraph()
            g11.SetLineWidth(3)
            g11.SetLineColor(33)
            g11.SetLineStyle(7)
            leg.AddEntry(g11,'DeepTagMD_bbvsLight','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                g11.SetPoint(g11.GetN(), s, b)
            g11.Draw("SAME L")
            leg.Draw("SAME")

    c1.Print("%s"% OUTPUTDIR +channel_name+"/btagger_combined_ROC.png")
    c1.Print("%s"% OUTPUTDIR +channel_name+"/btagger_combined_ROC.pdf")
    #c1.Print("%s"% OUTPUTDIR +channel_name+"/tagger_combined_ROC.png")
    #c1.Print("%s"% OUTPUTDIR +channel_name+"/tagger_combined_ROC.pdf")        
    if not gROOT.IsBatch(): raw_input("Press Enter to continue...")

def plotROC_wp(var_list, cut,cut_sign):
    c1 = TCanvas("c1", "Signals", 800, 600)
    c1.GetPad(0).SetTopMargin(0.06)
    c1.GetPad(0).SetRightMargin(0.06)
    c1.GetPad(0).SetTicky(2)
    gStyle.SetOptStat(0)
    #gStyle.SetOptTitle(0)
    c1.GetPad(0).SetLogy()
    leg = TLegend(0.6, 0.15, 0.9, 0.5)
    leg.SetBorderSize(0)
    channel = cut
    if "SB" in cut or "SR" in cut:
        channel_name = cut[:-2]
    else:
        channel_name = cut
    NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/Ntuples%s/"%(year)
    if 'inc' in cut:
        eventWeightLuminame = 'eventWeightLumi_nobtag'
    else:
        eventWeightLuminame = 'eventWeightLumi'
    background = back
    file = {}
    H_csv1 = []
    H_deepcsv1 = []
    H_deepcsv1_wp = []

    wp = [0.1241, 0.4184]

    for var in var_list:
        varname = var.replace('.', '_').replace('()', '')
        tree = {}
        hist = {}
        nbins = variable[var]['nbins']*100
    
        ### Create and fill MC histograms ###
        for i, s in enumerate(background+sign):
            tree[s] = TChain("tree")
            for j, ss in enumerate(sample[s]['files']):
                tree[s].Add(NTUPLEDIR + ss + ".root")
            hist[s] = TH1F(s, ";"+variable[var]['title']+";Events;"+('log' if variable[var]['log'] else ''), nbins, variable[var]['min'], variable[var]['max'])
            hist[s].Sumw2()
            if 'XZH_M' in s:
                tree[s].Project(s, var, "(%s)*(" % eventWeightLuminame +selection[cut_sign]+")")
            else:
                tree[s].Project(s, var, "(%s)*(" % eventWeightLuminame +selection[cut]+")")
        BkgSum = hist[background[0]].Clone("BkgSum")
        BkgSum.Reset("MICES")
        for i, s in enumerate(background): BkgSum.Add(hist[s], 1)
        SigSum = hist[sign[0]].Clone("SigSum")
        SigSum.Reset("MICES")
        for i, s in enumerate(sign): SigSum.Add(hist[s], 1)
        if var == 'H_deepcsv1':
            g = TGraph()
            g.SetLineWidth(3)
            g.SetLineColor(2)
            g.SetTitle("b-tagger comparison;signal efficiency;background efficiency")
            g.GetXaxis().SetRangeUser(0., 1.)
            g.GetYaxis().SetRangeUser(1.e-3, 1.)
            leg.AddEntry(g,'deepcsv1','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                H_deepcsv1.append([b,s])
                g.SetPoint(g.GetN(),s,b)
            g.Draw("AL")
            drawCMS(-1, year, "Simulation", False)
        elif var == 'H_deepcsv2':
            g1 = TGraph()
            g1.SetLineWidth(3)
            g1.SetLineColor(3)
            leg.AddEntry(g1,'deepcsv2','L')
            g2 = TGraph()
            g2.SetLineWidth(3)
            g2.SetLineColor(12)
            leg.AddEntry(g2, 'deepcsv', 'L')
            for i in range(1,nbins+1):
                b_1 = H_deepcsv1[i-1][0]
                s_1 = H_deepcsv1[i-1][1]
                b_2 = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s_2 = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                s = (s_1+s_2)/2
                b = (b_1+b_2)/2
                g1.SetPoint(g1.GetN(), s_2, b_2)
                g2.SetPoint(g2.GetN(), s, b)
            g1.Draw("SAME L")
            g2.Draw("SAME L")
        elif var == 'H_dbt':
            g3 = TGraph()
            g3.SetLineWidth(3)
            g3.SetLineColor(5)
            leg.AddEntry(g3,'DoubleB','L')
            for i in range(1,nbins+1):
                b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
                s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
                g3.SetPoint(g3.GetN(), s, b)
            g3.Draw("SAME L")
        #elif var == 'DeepTagMD_ZHbbvsQCD':
        #    g4 = TGraph()
        #    g4.SetLineWidth(3)
        #    g4.SetLineColor(6)
        #    leg.AddEntry(g4,'DeepTagMD_ZHbbvsQCD','L')
        #    for i in range(1,nbins+1):
        #        b = BkgSum.Integral(i, nbins)/BkgSum.Integral(1, nbins)
        #        s = SigSum.Integral(i, nbins)/SigSum.Integral(1, nbins)
        #        g4.SetPoint(g4.GetN(), s, b)
        #    g4.Draw("SAME L")

    #tree = {}
    #hist = {}

    #nbins = variable[var]['nbins']*100
    #for i, s in enumerate(background+sign):
    #    tree[s] = TChain("tree")
    #    for j, ss in enumerate(sample[s]['files']): 
    #        tree[s].Add(NTUPLEDIR + ss + ".root")
    #    hist[s] = TH1F(s, "DeepTagMD", nbins, 0, 1.)
    #    hist[s].Sumw2()
    #    if 'XZH_M' in s:
    #        tree[s].Project(s, "DeepTagMD_ZHbbvsQCD", "(eventWeightLumi)*("+selection[cut_sign]+")")
    #    else:
    #        tree[s].Project(s, "DeepTagMD_ZHbbvsQCD", "(eventWeightLumi)*("+selection[cut]+")")
    #
    #hist['BkgSum'] = hist[back[0]].Clone("BkgSum")
    #hist['BkgSum'].Reset("MICES")
    #hist['SigSum'] = hist['BkgSum'].Clone("SigSum")
    #for i, s in enumerate(back): hist['BkgSum'].Add(hist[s])
    #for i, s in enumerate(sign): hist['SigSum'].Add(hist[s])    
    
    #binL = hist['BkgSum'].GetXaxis().FindBin(0.6)
    #binM = hist['BkgSum'].GetXaxis().FindBin(0.7)
    #binT = hist['BkgSum'].GetXaxis().FindBin(0.9)
    
    #gWP1 = TGraph()
    #gWP1.SetMarkerStyle(20)
    #gWP1.SetMarkerColor(7)
    #gWP1.SetMarkerSize(1.5)
    #leg.AddEntry(gWP1,'deeptag loose','P')
    #gWP1.SetPoint(0, hist['SigSum'].Integral(binL, nbins)/hist['SigSum'].Integral(0, nbins), hist['BkgSum'].Integral(binL, nbins)/hist['BkgSum'].Integral(0, nbins))
    #gWP2 = TGraph()
    #gWP2.SetMarkerStyle(34)
    #gWP2.SetMarkerColor(7)
    #gWP2.SetMarkerSize(1.5)
    #leg.AddEntry(gWP2,'deeptag medium','P')
    #gWP2.SetPoint(0, hist['SigSum'].Integral(binM, nbins)/hist['SigSum'].Integral(0, nbins), hist['BkgSum'].Integral(binM, nbins)/hist['BkgSum'].Integral(0, nbins))
    #gWP3 = TGraph()
    #gWP3.SetMarkerStyle(21)
    #gWP3.SetMarkerColor(7)
    #gWP3.SetMarkerSize(1.5)
    #leg.AddEntry(gWP3,'deeptag tight','P')
    #gWP3.SetPoint(0, hist['SigSum'].Integral(binT, nbins)/hist['SigSum'].Integral(0, nbins), hist['BkgSum'].Integral(binT, nbins)/hist['BkgSum'].Integral(0, nbins))

    tree = {}
    hist = {}

    nbins = variable[var]['nbins']*100
    for i, s in enumerate(background+sign):
        tree[s] = TChain("tree")
        for j, ss in enumerate(sample[s]['files']): 
            tree[s].Add(NTUPLEDIR + ss + ".root")
        hist[s] = TH2F(s, ";CSV subjet 1;CSV subjet2;", nbins, 0, 1., nbins, 0, 1.)
        hist[s].Sumw2()
        if 'XZH_M' in s:
            
            tree[s].Project(s, "H_deepcsv2:H_deepcsv1", "(%s)*(" % eventWeightLuminame +selection[cut_sign]+")")
        else:
            tree[s].Project(s, "H_deepcsv2:H_deepcsv1", "(%s)*(" % eventWeightLuminame +selection[cut]+")")

    hist['BkgSum'] = hist[back[0]].Clone("BkgSum")
    hist['BkgSum'].Reset("MICES")
    hist['SigSum'] = hist['BkgSum'].Clone("SigSum")
    for i, s in enumerate(back): hist['BkgSum'].Add(hist[s])
    for i, s in enumerate(sign): hist['SigSum'].Add(hist[s]) 
    
    binL = hist['BkgSum'].GetXaxis().FindBin(0.1241)
    binM = hist['BkgSum'].GetXaxis().FindBin(0.4184)

    gWP4 = TGraph()
    gWP4.SetMarkerStyle(20)
    gWP4.SetMarkerColor(4)
    gWP4.SetMarkerSize(1.5)
    leg.AddEntry(gWP4,'deepcsv loose,loose','P')
    gWP4.SetPoint(0, hist['SigSum'].Integral(binL, nbins, binL, nbins)/hist['SigSum'].Integral(0, nbins, 0, nbins), hist['BkgSum'].Integral(binL, nbins, binL, nbins)/hist['BkgSum'].Integral(0, nbins, 0, nbins))
    gWP5 = TGraph()
    gWP5.SetMarkerStyle(34)
    gWP5.SetMarkerColor(4)
    gWP5.SetMarkerSize(1.5)
    leg.AddEntry(gWP5,'deepcsv medium,loose','P')
    gWP5.SetPoint(0, hist['SigSum'].Integral(binM, nbins, binL, nbins)/hist['SigSum'].Integral(0, nbins, 0, nbins), hist['BkgSum'].Integral(binM, nbins, binL, nbins)/hist['BkgSum'].Integral(0, nbins, 0, nbins))
    gWP6 = TGraph()
    gWP6.SetMarkerStyle(21)
    gWP6.SetMarkerColor(4)
    gWP6.SetMarkerSize(1.5)
    leg.AddEntry(gWP6,'deepcsv medium,medium','P')
    gWP6.SetPoint(0, hist['SigSum'].Integral(binM, nbins, binM, nbins)/hist['SigSum'].Integral(0, nbins, 0, nbins), hist['BkgSum'].Integral(binM, nbins, binM, nbins)/hist['BkgSum'].Integral(0, nbins, 0, nbins))
    #gWP7 = TGraph()
    #gWP7.SetMarkerStyle(34)
    #gWP7.SetMarkerColor(5)
    #gWP7.SetMarkerSize(1.4)
    #leg.AddEntry(gWP7,'loose,medium','P')
    #gWP7.SetPoint(0, hist['SigSum'].Integral(binL, nbins, binM, nbins)/hist['SigSum'].Integral(0, nbins, 0, nbins), hist['BkgSum'].Integral(binL, nbins, binM, nbins)/hist['BkgSum'].Integral(0, nbins, 0, nbins)) 
   
    #gWP1.Draw("SAME P")
    #gWP2.Draw("SAME P")
    #gWP3.Draw("SAME P")
    gWP4.Draw("SAME P")
    gWP5.Draw("SAME P")
    gWP6.Draw("SAME P")
    leg.Draw("SAME")

    c1.Print("%s"% OUTPUTDIR +channel_name+"/btagger_wp_ROC.png")
    c1.Print("%s"% OUTPUTDIR +channel_name+"/btagger_wp_ROC.pdf")
    if not gROOT.IsBatch(): raw_input("Press Enter to continue...")

def checkcut(checkvar,cut):
    background = back
    var = 'X_mass'
    channel = cut
    sign = ['XZH_M1400']
    lower_bound = 1300.
    upper_bound = 1500.
    if "SB" in cut or "SR" in cut:
        channel_name = cut[:-2]
    else:
        channel_name = cut
    NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/Ntuples%s/"%(year)
    if 'inc' in cut:
        eventWeightLuminame = 'eventWeightLumi_nobtag'
    else:
        eventWeightLuminame = 'eventWeightLumi'
    for cut_var in np.arange(0.1,3.0,0.1):
        cut_tag = " && %s < %s && X_mass > %f && X_mass < %f" %(checkvar,cut_var,lower_bound,upper_bound)
        cut_mod = "%s%s" %(selection[cut],cut_tag)
        varname = var.replace('.', '_').replace('()', '')
        tree = {}
        hist = {}
        nbins = variable[var]['nbins']
    
        ### Create and fill MC histograms ###
        for i, s in enumerate(background+sign):
            tree[s] = TChain("tree")
            for j, ss in enumerate(sample[s]['files']):
                tree[s].Add(NTUPLEDIR + ss + ".root")
            hist[s] = TH1F(s, ";"+variable[var]['title']+";Events;"+('log' if variable[var]['log'] else ''), nbins, variable[var]['min'], variable[var]['max'])
            hist[s].Sumw2()
            tree[s].Project(s, var, "(%s)*(" % eventWeightLuminame +cut_mod+")")
        BkgSum = hist[background[0]].Clone("BkgSum")
        BkgSum.Reset("MICES")
        for i, s in enumerate(background): BkgSum.Add(hist[s], 1)
        SigSum = hist[sign[0]].Clone("SigSum")
        SigSum.Reset("MICES")
        for i, s in enumerate(sign): SigSum.Add(hist[s], 1)
        background_integral = BkgSum.Integral()
        signal_integral = SigSum.Integral()
        #print "cut %f signal integral:"%cut_var, signal_integral
        #print "cut %f background integral:" %cut_var,background_integral
        print "cut %f signal/background:"%(cut_var), signal_integral/np.sqrt(background_integral+1)
        


def plotAll():
    gROOT.SetBatch(True)
    for c in ['mmincSB','eeincSB','nnincSB']:
        for v in ['eecutflow_inc','mmcutflow_inc','nncutflow_inc','nPV','H_mass','MET','nJets','nFatJets','MinDPhi','DPhi','nElectrons','nMuons','nTaus','H_dbt','H_deepcsv1','H_deepcsv2','H_pt','X_mass', 'V_pt','VH_deltaR','Mu1_pt','Mu1_eta','Mu2_pt','Mu2_eta','Ele1_pt','Ele1_eta','Ele2_pt','Ele2_eta','V_mass','X_pt','DEta','PrefireWeight','H_nhf','H_chf']:
            plot(v, c)
    for c in ['mmbbSB','mm0bSB','eebbSB','ee0bSB','nnbbSB','nn0bSB']:
        for v in ['nPV','H_mass','MET','nJets','nFatJets','MinDPhi','DPhi','nElectrons','nMuons','nTaus','H_dbt','H_deepcsv1','H_deepcsv2','H_pt','X_mass', 'V_pt','VH_deltaR','V_mass','X_pt','DEta']:
            plot(v, c)
    for c in ['mmincVBFSB','eeincVBFSB','nnincVBFSB']:
        for v in ['nPV','H_mass','MET','nJets','nFatJets','MinDPhi','DPhi','nElectrons','nMuons','nTaus','H_dbt','H_deepcsv1','H_deepcsv2','H_pt','X_mass', 'V_pt','VH_deltaR','Mu1_pt','Mu1_eta','Mu2_pt','Mu2_eta','Ele1_pt','Ele1_eta','Ele2_pt','Ele2_eta','V_mass','X_pt','DEta','PrefireWeight','deltaR_VBF','deltaR_HVBFjet1','deltaR_HVBFjet2','H_nhf','H_chf']:
            plot(v,c)
    for c in ['mmbbVBFSB','mm0bVBFSB','eebbVBFSB','ee0bVBFSB','nnbbVBFSB','nn0bVBFSB']:
        for v in ['nPV','H_mass','MET','nJets','nFatJets','MinDPhi','DPhi','nElectrons','nMuons','nTaus','H_dbt','H_deepcsv1','H_deepcsv2','H_pt','X_mass', 'V_pt','VH_deltaR','V_mass','X_pt','DEta']:
            plot(v, c)
   
"""
def plottagger():
    gROOT.SetBatch(True)
    for c in ['nnincSR','mmincSR','eeincSR']:
        c_sign = '%s_sign' %c
        v_list = ['H_csv1','H_csv2','H_deepcsv1','H_deepcsv2','BtagDeepB','H_dbt','DeepTagMD_HbbvsQCD','DeepTagMD_ZHbbvsQCD','DeepTagMD_ZbbvsQCD','DeepTagMD_bbvsLight']
        #v_list = ['H_tau21','H_tau31','H_tau32','H_tau41','H_tau42','H_ddt','DeepTagMD_H4qvsQCD','DeepTagMD_WvsQCD','DeepTagMD_ZvsQCD']
        plotROC_combined(v_list,c,c_sign)
"""
def plottagger():
    gROOT.SetBatch(True)
    for c in ['nnincSR','mmincSR','eeincSR']:
        c_sign = '%s_sign' %c
        v_list = ['H_deepcsv1','H_deepcsv2','H_dbt','DeepTagMD_ZHbbvsQCD']
        plotROC_wp(v_list,c,c_sign)
"""
def plottagger():
    gROOT.SetBatch(True)
    for c in ['mminc','eeinc','nninc']:
        for t in ['H_deepcsv2','H_dbt','DeepTagMD_ZHbbvsQCD']:
            checktagger(t,c)

def plottagger():
    gROOT.SetBatch(True)
    for c in ['mminc']:
        for t in ['H_deepcsv2','DeepTagMD_ZHbbvsQCD']:
            getwp(t,c)
"""
if options.all: plotAll()
elif options.tagger: plottagger()
elif options.norm: plotNorm(options.variable, options.cut)
elif options.top: plotTop()
#elif options.top: plotTopError()
elif options.checkcut: checkcut(options.variable,options.cut)
else: plot(options.variable, options.cut)

