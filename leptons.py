#! /usr/bin/env python

import os, sys, getopt
import copy, math
from array import array
from ROOT import gROOT, gStyle, gRandom
from ROOT import TFile, TChain, TTree, TCut, TH1F, TH2F, THStack, TGraph, TGaxis, TEfficiency
from ROOT import TStyle, TCanvas, TPad, TLegend, TLatex, TText

from samples import sample
from utils import *

gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.05)

#gStyle.SetErrorX(1)
NTUPLEDIR   = "/scratch/zucchett/Ntuple/VH/"
RATIO       = 4

back = ["DYJetsToLL_HT"]
sign = ["XZH_M800", "XZH_M1000", "XZH_M1200", "XZH_M1400", "XZH_M1800", "XZH_M2000", "XZH_M2500", "XZH_M3000", "XZH_M3500", "XZH_M4000", "XZH_M4500"]
#sign = ["XZZ_M2000", "XZZ_M2500", "XZZ_M3000", "XZZ_M3500", "XZZ_M4000", "XZZ_M4500"]
sigp = ['XZh_M1000', 'XZh_M1400', 'XZh_M2000', 'XZh_M2500', 'XZh_M3000', 'XZh_M3500', 'XZh_M4000', 'XZh_M4500']

def histograms(s = "XZh_M2000"):
    
    effh = []
    gens = []
    
    inFile = TFile(NTUPLEDIR + sample[s]['files'][0] + ".root", "READ")
    inFile.cd()
    for key in inFile.GetListOfKeys():
        obj = key.ReadObj()
          
        # Copy and rescale histograms
        if obj.IsA().InheritsFrom("TH1"):
            plotHist([obj])
        elif obj.IsA().InheritsFrom("TTree"):
            print "Skip tree", obj.GetName()
            # do nothing
        # Directories
        elif obj.IsFolder():
            subdir = obj.GetName()
            print "In directory", subdir
            inFile.cd(subdir)
            for subkey in inFile.GetDirectory(subdir).GetListOfKeys():
                subobj = subkey.ReadObj()
                if subobj.IsA().InheritsFrom("TH1"):
                    if 'Eff' in subdir: effh.append(subobj.GetName())
                    elif 'Gen' in subdir: gens.append(subobj.GetName())
                    elif not 'Counter' in subdir: plotHist([subobj])
            inFile.cd("..")
    inFile.Close()
    
    print "Now plotting efficiency histograms", effh

    for j, h in enumerate(effh):
        name = h.split("_", 1)[0].replace("Eff", "") #if not 'HEEP_E' in h else h.replace("Eff", "")
        num = TH1F()
        den = TH1F()
        files = {}
        for i, s in enumerate(sign):
            files[s] = TFile(NTUPLEDIR + sample[s]['files'][0] + ".root", "READ")
            files[s].cd()
            tnum = files[s].Get("Eff/"+h)
            tden = files[s].Get("Leptons/"+name)
            tnum.Multiply(tden)
            if i==0:
                num = tnum.Clone("Numerator")
                den = tden.Clone("Denominator")
            else:
                num.Add(tnum)
                den.Add(tden)
            #files[s].Close()
        num.Divide(den)
        num.SetName(h)
        plotHist([num])
    
    
    print "Now plotting gen signal histograms", effh

    for j, h in enumerate(gens):
        files = {}
        sigh = []
        for i, s in enumerate(sigp):
            files[s] = TFile(NTUPLEDIR + sample[s]['files'][0] + ".root", "READ")
            files[s].cd()
            sigh.append(files[s].Get("Gen/"+h))
            sigh[-1].SetTitle(s)
        plotHist(sigh)
            

def compareHist(variables, obj):
    
    num = {}
    den = {}
    eff = {}
    files = {}
    for i, s in enumerate(sign):
        for j, ss in enumerate(sample[s]['files']):
            files[ss] = TFile(NTUPLEDIR + ss + ".root", "READ")
    
    for variable in variables:
        var = variable[0]
        denname = var.split("_")[0]
        for i, f in files.iteritems():
#            tnum = f.Get("Muons/"+var)
#            tden = f.Get("Muons/"+denname)
            tnum = f.Get(obj+"/"+var)
            tden = f.Get(obj+"/"+denname)
            if tnum==None or tden==None: continue
            if not var in num.keys():
                num[var] = tnum.Clone("Numerator")
                den[var] = tden.Clone("Denominator")
            else:
                num[var].Add(tnum)
                den[var].Add(tden)
            #files[s].Close()
        num[var].Divide(den[var]) #eff[var] = TEfficiency(num[var], den[var])
        num[var].SetName("Multi_"+var)
        num[var].SetTitle(variable[1])
        num[var].GetYaxis().SetRangeUser(0., 1.)
        #print num[var].Integral()
    
    h = num.values()
    h.sort(key=lambda x: x.Integral(), reverse=True)
    plotHist(h)



def plotHist(h):
    colors = [633, 416+1, 862, 800, 921, 922, 1]
    n = len(h)
    
    c1 = TCanvas("c1", "Signals", 800, 600)
    c1.cd()
    #c1.GetPad(1).SetPadTopMargin(0.06)
    #c1.GetPad(1).SetPadRightMargin(0.05)

    c1.GetPad(0).SetTicky(2)
    
    for i in range(n):
        col = sample[h[i].GetTitle()]['linecolor'] if n > len(colors) else colors[i]
        h[i].SetMarkerColor(col)
        h[i].SetLineColor(col)
    
    if "Counter" in h[0].GetName():
        h[0].GetYaxis().SetRangeUser(0., -1)
    else:
        for i in range(n):
            h[i].SetMarkerStyle(20)
            #h[i].SetMarkerSize(1)
            h[i].SetLineWidth(2)
    
    ymax = 0.
    for i in range(n):
        if h[i].GetMaximum() > ymax: ymax = h[i].GetMaximum()
    h[0].SetMaximum(ymax*1.2)
    if h[0].Integral()/h[0].GetNbinsX() < 1.:
        #h[0].GetYaxis().SetTitle("Efficiency")
        h[0].GetYaxis().SetRangeUser(0., 1.)
    
    style = "HIST" if "Counter" in h[0].GetName() or "Barrel" in h[0].GetName() or "Endcap" in h[0].GetName() or h[0].GetName().startswith("Gen") else "PE1"
    h[0].Draw(style)
    for i in range(1, n): h[i].Draw("SAME, "+style)
    
    leg = TLegend(0.50, 0.55-0.05*n, 0.85, 0.55)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0) #1001
    leg.SetFillColor(0)
    for i in range(n):
        if h[i].GetTitle()!="":
            if h[i].GetTitle() in sample.keys(): leg.AddEntry(h[i], sample[h[i].GetTitle()]['label'], "pl")
            else: leg.AddEntry(h[i], h[i].GetTitle(), "pl")
    leg.Draw()
    drawCMS(-1, "Simulation", True)
    c1.Print("plots/Leptons/" + h[-1].GetName() + ".png")
    c1.Print("plots/Leptons/" + h[-1].GetName() + ".pdf")
    #raw_input("Press Enter to continue...")
    c1.Close()
    return True


#histograms()

#compareHist(variables = {"EffMuonZdR_TightTight" : "2 Tight muons;#Delta R_{#mu#mu} at gen level;Muon Reco+Id efficiency", "EffMuonZdR_HighptHighpt" : "2 HighPt muons;#Delta R_{#mu#mu} at gen level;Muon Reco+Id efficiency", "EffMuonZdR_HighptLoose" : "1 HighPt, 1 Loose;#Delta R_{#mu#mu} at gen level;Muon Reco+Id efficiency", "EffMuonZdR_HighptCustomTracker" : "1 HighPt, 1 Tracker;#Delta R_{#mu#mu} at gen level;Muon Reco+Id efficiency"})
#compareHist(variables = {"EffMuonZdR_HighptCustomTracker" : "1 HighPt, 1 Tracker, no iso;#Delta R_{#mu#mu} at gen level;Muon Reco+Id+Iso efficiency", "EffMuonZdR_HighptCustomTracker_pfIso" : "1 HighPt, 1 Tracker, PF iso 0.4;#Delta R_{#mu#mu} at gen level;Muon Reco+Id+Iso efficiency", "EffMuonZdR_HighptCustomTracker_miniIso" : "1 HighPt, 1 Tracker, mini iso;#Delta R_{#mu#mu} at gen level;Muon Reco+Id+Iso efficiency", "EffMuonZdR_HighptCustomTracker_trkIso" : "1 HighPt, 1 Tracker, tracker iso;#Delta R_{#mu#mu} at gen level;Muon Reco+Id+Iso efficiency"})

#compareHist(variables = {"EffElecZdR_Veto" : "Veto;#Delta R_{ee} at gen level;Electron Reco+Id efficiency", "EffElecZdR_Loose" : "Loose;#Delta R_{ee} at gen level;Electron Reco+Id efficiency", "EffElecZdR_Tight" : "Tight;#Delta R_{ee} at gen level;Electron Reco+Id efficiency", "EffElecZdR_HEEP" : "HEEP no iso;#Delta R_{ee} at gen level;Electron Reco+Id efficiency"})
#compareHist(variables = {"EffElecZdR_HEEP_pfIso" : "HEEP no iso + PF iso;#Delta R_{ee} at gen level;Electron Reco+Id+Iso efficiency", "EffElecZdR_HEEP_miniIso" : "HEEP no iso + mini iso;#Delta R_{ee} at gen level;Electron Reco+Id+Iso efficiency"})

#compareHist(variables = {"m_dR_reco" : "reco", "m_dR_pt" : "p_{T} thresholds", "m_dR_id1" : "high-p_{T} id", "m_dR_id2" : "tracker id", "m_dR_iso" : "tracker iso"})


#compareHist(variables = [["dR_HighptHighptId", "2 High-p_{T} Id"], ["dR_TightTightId", "2 Tight Id"], ["dR_MediumMediumId", "2 Medium Id"], ["dR_LooseLooseId", "2 Loose Id"], ["dR_HighptTrackerId", "1 High-p_{T} + 1 Tracker Id"]], obj = "Muons")
#compareHist(variables = [["dR_HighptTrackerId", "1 High-p_{T} + 1 High-p_{T} trk Id"], ["dR_HighptTrackerIdTrackerIso", "Id + modified trk Iso"], ["dR_HighptTrackerIdPFIso", "Id + PF Iso"]], obj = "Muons")
#compareHist(variables = [["pT1_HighptId", "Muon 1, HighPt Id"], ["pT2_HighptId", "Muon 2, HighPt Id"], ["pT1_TrackerId", "Muon 1, Tracker Id"], ["pT2_TrackerId", "Muon 2, Tracker Id"]], obj = "Muons")
#compareHist(variables = [["dR_VetoVetoId", "2 Veto Id"], ["dR_LooseLooseId", "2 Loose Id"], ["dR_MediumMediumId", "2 Medium Id"], ["dR_TightTightId", "2 Tight Id"], ["dR_HeepHeepId", "2 HEEP Id"]], obj = "Electrons")
#compareHist(variables = [["pT1_LooseId", "Electron 1, Loose Id"], ["pT2_LooseId", "Electron 2, Loose Id"]], obj = "Electrons")
compareHist(variables = [["dR_HighptTrackerId", "1 High-p_{T} + 1 High-p_{T} trk Id"], ["dR_HighptTrackerIdTrackerIso", "Id + modified trk Iso"], ["dR_HighptTrackerIdTrackerNonCorrIso", "Id + non corr. trk Iso"]], obj = "Muons")
