#! /usr/bin/env python

import os, sys, getopt
import copy, math
from array import array
from ROOT import gROOT, gStyle, gRandom
from ROOT import TFile, TChain, TTree, TCut, TH1F, TH2F, TF1, THStack, TGraph, TGraphErrors, TGraphAsymmErrors, TGaxis, TVirtualFitter
from ROOT import TStyle, TCanvas, TPad, TLegend, TLatex, TText, TPoint

from samples import sample_2016, sample_2017, sample_2018, sample
from variables import variable
from selections import selection
from utils import *

gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.06)

#gStyle.SetErrorX(1)
NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/combined_weighted/"
RATIO       = 4

########## SAMPLES ##########
data = ['data_obs']
back = ["VV", "ST", "TTbarSL", "WJetsToLNu_HT", "DYJetsToNuNu_HT", "DYJetsToLL_HT"] 
sign = ['XZH_M800','XZH_M1000','XZH_M1200','XZH_M1400','XZH_M1600', 'XZH_M1800', 'XZH_M2000', 'XZH_M2500','XZH_M3000','XZH_M3500', 'XZH_M4000','XZH_M4500','XZH_M5000','XZH_M5500','XZH_M6000']

colors = [616+4, 632, 800+7, 800, 416+1, 860+10, 600, 616, 921, 922]

def systematics(sys, isShape=False, antiCorr=False):
    treeRead = True
    file = {}
    hist = {}
    tree = {}
    histUp = {}
    histDown = {}
    up = {}
    down = {}
    gUp = TGraph()
    gDown = TGraph()
    gUp.SetLineWidth(3)
    gDown.SetLineWidth(3)
    gUp.SetLineColor(632)
    gDown.SetLineColor(602)
    BTagAK4deepup = False
    BTagAK4deepdown = False
    BTagAK4deep = False
    BTagAK8deepup = False
    BTagAK8deepdown = False
    BTagAK8deep = False
#    g.SetMarkerStyle(20)
#    g.SetMarkerColor(418)
#    g.SetMarkerSize(1.25)
    var = sys
    cut = 'nnbbSR'
    if var=='BTagAK4Weight_deep_up':
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
    upAvg, downAvg, upMin, downMin, upMax, downMax = 0., 0., 2., 2., 0., 0.
    for k in sorted(selection.keys(), key=len, reverse=True):
        if k in cut: 
            cut = cut.replace(k, selection[k])
    for i, s in enumerate(sign):
        if '_MZ' in s: m = int((s.split('_MZ')[1]).split('_MA')[0])
        elif '_M' in s: m = int(s.split('_M')[1])
        else: m = 0
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
            cutstring = "(eventWeightLumi)" + ("*("+cut+")")
            if var=='LeptonWeightUp':
                cutstring = "(eventWeightLumi * LeptonWeightUp/LeptonWeight)" + ("*("+cut+")") 
            elif var=='LeptonWeightDown':
                cutstring = "(eventWeightLumi * LeptonWeightDown/LeptonWeight)" + ("*("+cut+")") 
            elif var=='TriggerWeightUp':
                cutstring = "(eventWeightLumi * TriggerWeightUp/TriggerWeight)" + ("*("+cut+")") 
            elif var=='TriggerWeightDown':
                cutstring = "(eventWeightLumi * TriggerWeightDown/TriggerWeight)" + ("*("+cut+")") 
            #division by BTagAk4Weight_deep is because the weighted samples are used
            elif BTagAK4deepup:
                cutstring = "(eventWeightLumi * BTagAK4Weight_deep_up/BTagAK4Weight_deep)" + ("*("+cut+")") 
            elif BTagAK4deepdown:
                cutstring = "(eventWeightLumi * BTagAK4Weight_deep_down/BTagAK4Weight_deep)" + ("*("+cut+")") 
            elif BTagAK8deep:
                cutstring = "(eventWeightLumi * BTagAK8Weight_deep)" + ("*("+cut+")")    
            elif BTagAK8deepup:
                cutstring = "(eventWeightLumi * BTagAK8Weight_deep_up)" + ("*("+cut+")") 
            elif BTagAK8deepdown:
                cutstring = "(eventWeightLumi * BTagAK8Weight_deep_down)" + ("*("+cut+")")     
            tree[s].Project(s, var, cutstring)
            if not tree[s].GetTree()==None: hist[s].SetOption("%s" % tree[s].GetTree().GetEntriesFast())
        """
        for j, ss in enumerate(sample[s]['files']):
            file[ss] = TFile(NTUPLEDIR + ss + ".root", "READ")
            tmp = file[ss].Get("Sys/"+sys)
            if tmp == None: continue
            if not s in hist.keys(): hist[s] = tmp
            else:
                if antiCorr:
                    for x in range(1, hist[s].GetNbinsX()+1): hist[s].SetBinContent(x, hist[s].GetBinContent(x)+tmp.GetBinContent(tmp.GetNbinsX()+1-x))
                else: hist[s].Add( tmp )
            
            if isShape:
                tmp = file[ss].Get("Sys/"+sys+"_up")
                print "Sys/"+sys+"_up"
                if tmp == None: continue
                if not s in histUp.keys(): histUp[s] = tmp
                else: histUp[s].Add( tmp )
                #
                tmp = file[ss].Get("Sys/"+sys+"_down")
                if tmp == None: continue
                if not s in histDown.keys(): histDown[s] = tmp
                else: histDown[s].Add( tmp )
            
        if 'accept' in sys:
            norm = None
            for j, ss in enumerate(sample[s]['files']):
                tmp = file[ss].Get("Sys/PDF_scale")
                if tmp == None: continue
                if norm == None: norm = tmp
                else: norm.Add( tmp )
            hist[s].Divide(norm)
        """
        if(isShape):
            shape = TF1("shape", "gaus", 0, 5000)
            shapeUp = TF1("shapeUp", "gaus", 0, 5000)
            shapeDown = TF1("shapeDown", "gaus", 0, 5000)
            hist[s].Fit(shape, "Q0", "")
            histUp[s].Fit(shapeUp, "Q0", "")
            histDown[s].Fit(shapeDown, "Q0", "")
            if 'scale' in sys or 'unc' in sys:
                up[s] = histUp[s].GetMean()/hist[s].GetMean()
                down[s] = histDown[s].GetMean()/hist[s].GetMean()
#                up[s] = shapeUp.GetParameter(1)/shape.GetParameter(1)
#                down[s] = shapeDown.GetParameter(1)/shape.GetParameter(1)
            elif 'res' in sys:
                up[s] = histUp[s].GetRMS()/hist[s].GetRMS()
                down[s] = histDown[s].GetRMS()/hist[s].GetRMS()
#                up[s] = shapeUp.GetParameter(2)/shape.GetParameter(2)
#                down[s] = shapeDown.GetParameter(2)/shape.GetParameter(2)
        else:
            up[s] = hist[s].GetBinContent(hist[s].FindBin(+1))/hist[s].GetBinContent(hist[s].FindBin(0))
            down[s] = hist[s].GetBinContent(hist[s].FindBin(-1))/hist[s].GetBinContent(hist[s].FindBin(0))
        gUp.SetPoint(i, m, up[s])
        gDown.SetPoint(i, m, down[s])
        #if mass < 1000: continue
        upAvg += up[s]
        downAvg += down[s]
        if abs(up[s]) > upMax: upMax = abs(up[s])
        if abs(up[s]) < upMin: upMin = abs(up[s])
        if abs(down[s]) > downMax: downMax = abs(down[s])
        if abs(down[s]) < downMin: downMin = abs(down[s])
    
    upAvg /= len(sign)
    downAvg /= len(sign)

    print " ---", sys, "--- | up: %.3f, down: %.3f, average: %.3f" % (upAvg, downAvg, abs(upAvg -1 + 1.-downAvg)/2.), "|", "^{%.1f-%.1f}_{%.1f-%.1f}" % (100.*(1.-upMin), 100.*(1.-upMax), 100.*(downMin-1.), 100.*(downMax-1.))
    
    c1 = TCanvas("c1", "Signals", 800, 600)
    c1.cd()
    c1.GetPad(0).SetTicky(2)
    
    gUp.Draw("AL")
    gDown.Draw("SAME, L")
    gUp.GetYaxis().SetRangeUser(0.65, 1.35)
    gUp.GetXaxis().SetTitle("m_{X} (GeV)")
    gUp.GetYaxis().SetTitle("Uncertainty")
    
    leg = TLegend(0.5, 0.90-0.20, 0.9, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0) #1001
    leg.SetFillColor(0)
    leg.SetHeader(sys.replace('_', ' '))
    leg.AddEntry(gUp, "+1 s. d. (%.1f%%)" % (100.*(upAvg-1.)), "l")
    leg.AddEntry(gDown, " -1 s. d. (%.1f%%)" % (100.*(1.-downAvg)), "l")
    leg.Draw()
    
    drawCMS(-1, "Simulation", False)
    c1.Update()
    
    filename = sys
    if sys.startswith('W_mass') or sys.startswith('Z_mass'): filename += "_" + sign[0][:3]
#    c1.Print("plots/Systematics/"+filename+".png")
#    c1.Print("plots/Systematics/"+filename+".pdf")
    if not gROOT.IsBatch(): raw_input("Press Enter to continue...")
    if 'doubleB' in sys or 'subjet' in sys or 'mass' in sys:
        print sys + " = {",
        for m in [800, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 2000, 2500, 3000, 3500, 4000, 4500]:
            print "%d : [%.3f, %.3f], " % (m, gUp.Eval(m), gDown.Eval(m)), 
        print "}"
    if 'extr' in sys or 'tagging' in sys:
        print sys + " = {",
        for m in [800, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500]:
            print "%d : [%.3f, %.3f], " % (m, gUp.Eval(m), gDown.Eval(m)), 
        print "}"
    if 'QCD_scale' in sys or 'PDF_scale' in sys:
        print sys + " = {",
        for m in range(800, 4500+1, 10):
            print "%d : [%.3f, %.3f], " % (m, gUp.Eval(m), gDown.Eval(m)), 
        print "}"
    

def systematicsBkg(sys, s):
    file = {}
    hist, tmp = None, None
    for j, ss in enumerate(sample[s]['files']):
        file[ss] = TFile(NTUPLEDIR + ss + ".root", "READ")
        tmp = file[ss].Get("Sys/"+sys)
        if tmp.GetBinContent(tmp.FindBin(+1)) < 1.e-6: continue
        if hist == None: hist = tmp
        else: hist.Add( tmp )
    up = hist.GetBinContent(hist.FindBin(+1))/hist.GetBinContent(hist.FindBin(0))
    down = hist.GetBinContent(hist.FindBin(-1))/hist.GetBinContent(hist.FindBin(0))
    print " ---", sys, "--- | up: %.3f, down: %.3f, average: %.3f" % (up, down, abs(up -1 + 1.-down)/2.)

#systematics("J_energy_scale")
#systematics("J_energy_res")
#systematics("H_mass_scale")
#systematics("H_mass_res")
#systematics("V_mass_scale")
#systematics("V_mass_res")
#systematics("W_mass_scale")
#systematics("Z_mass_scale")
#systematics("W_mass_res", antiCorr=True) ###
#systematics("Z_mass_res", antiCorr=True) ###
#systematics("H_subjet_notag")
systematics("BTagAK4Weight_deep")
#systematics("H_subjet_loose")
#systematics("H_subjet_tight")
#systematics("H_doubleB_notag")
#systematics("H_doubleB_loose")
#systematics("H_doubleB_tight")
#systematics("BTagSF_veto")
#systematics("V_tau21_extr_HP")
#systematics("V_tau21_extr_LP")
#systematics("H_tagging")
#systematics("ElecT_eff")
#systematics("Elec1_eff")
#systematics("Elec2_eff")
#systematics("MuonT_eff")
#systematics("Muon1_eff")
#systematics("Muon2_eff")
#systematics("MuonTrkIso_eff")
#systematics("Pileup")
#systematics("PDF_scale")
#systematics("PDF_accept")
#systematics("QCD_scale")
#systematics("X_mass_scale", True)
#systematics("X_mass_res", True)
#systematics("X_mass_unc", True)
#systematicsBkg("PDF_scale", "ST")
#systematicsBkg("PDF_accept", "ST")
#systematicsBkg("QCD_scale", "ST")
#systematicsBkg("PDF_scale", "VV")
#systematicsBkg("PDF_accept", "VV")
#systematicsBkg("QCD_scale", "VV")
