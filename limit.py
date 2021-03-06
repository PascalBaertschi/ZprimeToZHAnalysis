#! /usr/bin/env python

import os, sys, getopt
import copy
import time
import math
from array import array
from ROOT import gROOT, gRandom
from ROOT import TFile, TTree, TCut, TH1F, TH2F, TGraph, TGraph2D, TGraphErrors, TGraphAsymmErrors
from ROOT import TStyle, TCanvas, TPad
from ROOT import TLegend, TLatex, TText, TColor

#from DMPD.Heppy.tools.samples import *
from utils import *
from hvtXs import HVT
from thdmXs import THDM

# Combine output
# 0 - Observed Limit
# 1 - Expected  2.5%
# 2 - Expected 16.0%
# 3 - Expected 50.0%
# 4 - Expected 84.0%
# 5 - Expected 97.5%
# 6 - Significance
# 7 - p-value
# 8 - Best fit r
# 9 - Best fit r down
#10 - Best fit r up

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-c", "--category", action="store", type="string", dest="category", default="")
parser.add_option("-m", "--method", action="store", type="string", dest="method", default="")
parser.add_option("-a", "--all", action="store_true", default=False, dest="all")
parser.add_option("-b", "--bash", action="store_true", default=False, dest="bash")
parser.add_option("-B", "--blind", action="store_true", default=False, dest="blind")
(options, args) = parser.parse_args()
if options.bash: gROOT.SetBatch(True)
gStyle.SetOptStat(0)

ANALYSIS    = "ZH"
HTOBB       = 0.5824#0.577
LUMI        =  137190.
MAXIMUM     = {'XVHsl' : 6000., 'XWHsl' : 6000., 'XZHsl' : 6000.}
YEAR        = 'combined'

#signals = [800, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000]
signals = range(800, 5000+1, 100)

VBr = {"XZHll" : 0.05826, "XZHee" : 0.05826/3., "XZHmm" : 0.05826/3., "XZHnn" : 0.1154, "XWHln" : 0.18712, "XWHen" : 0.18712/3., "XWHmn" : 0.18712/3.}
XVHahExp = { 1000.0 : 69.8125 , 1100.0 : 35.5 , 1200.0 : 20.0625 , 1300.0 : 14.9375 , 1400.0 : 11.8438 , 1500.0 : 9.4062 , 1600.0 : 7.4688 , 1700.0 : 6.0156 , 1800.0 : 5.0 , 1900.0 : 4.2344 , 2000.0 : 3.6406 , 2100.0 : 3.2031 , 2200.0 : 2.8359 , 2300.0 : 2.5547 , 2400.0 : 2.3203 , 2500.0 : 2.1094 , 2600.0 : 1.9297 , 2700.0 : 1.7734 , 2800.0 : 1.6562 , 2900.0 : 1.5352 , 3000.0 : 1.4336 , 3100.0 : 1.3398 , 3200.0 : 1.2578 , 3300.0 : 1.1914 , 3400.0 : 1.1367 , 3500.0 : 1.082 , 3600.0 : 1.043 , 3700.0 : 1.0039 , 3800.0 : 0.9805 , 3900.0 : 0.957 , 4000.0 : 0.9336 , 4100.0 : 0.918 , 4200.0 : 0.8945 , 4300.0 : 0.8789 , 4400.0 : 0.8633 , 4500.0 : 0.8398 , }


decay = {'nn' : "#nu#nu", 'ln' : "l#nu", 'en' : "eu#nu", 'mn' : "#mu#nu", 'll' : "ll", 'ee' : "ee", 'mm' : "#mu#mu", 'qq' : "q#bar{q}", 'wr' : "q#bar{q}", 'zr' : "q#bar{q}", 'b' : "b#bar{q}", '' : "q#bar{q}"}

theoryLabel = {'B3' : "HVT model B", 'A1' : "HVT model A", 'T1' : "2HDM Type-I", 'T2' : "2HDM Type-II"}
theoryLineColor = {'B3' : 629, 'A1' : 616-3, 'T1' : 880-4, 'T2' : 602}
theoryFillColor = {'B3' : 625, 'A1' : 616-7, 'T1' : 880-9, 'T2' : 856}
theoryFillStyle = {'B3' : 3002, 'A1' : 3013, 'T1' : 3002, 'T2' : 3013}



def fillValues(filename):
    val = {}
    mass = []
    for i, s in enumerate(signals):
        try:
            file = open(filename % s, 'r')
            if file==None:
                print "Signal", filename % s, "does not exist"
                continue
            val[s] = file.read().splitlines()
            if len(val[s]) <= 1:
                signals.remove(s)
                print "Signal", filename % s, "has no values"
                continue
            for i, f in enumerate(val[s]): val[s][i] = float(val[s][i])
            if 'fullCls' in filename: val[1:5]=sorted(val[1:5])
            if not s in mass: mass.append(s)
        except:
            print "File", filename % s, "does not exist"
    return mass, val





def limit(method, channel):
    particle = channel[1:2]
    particleP = particle+"'" if channel.startswith('X') else channel[0]
    THEORY = ['A1', 'B3'] if channel.startswith('X') else []
    
    suffix = ""
    if method=="hvt": suffix="_HVT"
    if method=="cls": suffix="_CLs"
    if method=="monoH": suffix="_monoH"

    filename = "./combine/" + method + "/" + channel + "_M%d.txt"
    mass, val = fillValues(filename)
    
    Obs0s = TGraph()
    Exp0s = TGraph()
    Exp1s = TGraphAsymmErrors()
    Exp2s = TGraphAsymmErrors()
    Sign = TGraph()
    pVal = TGraph()
    Best = TGraphAsymmErrors()
    Theory = {}

    for i, m in enumerate(mass):
        if not m in val:
            print "Key Error:", m, "not in value map"
            continue
        n = Exp0s.GetN()
        Obs0s.SetPoint(n, m, val[m][0])
        Exp0s.SetPoint(n, m, val[m][3])
        Exp1s.SetPoint(n, m, val[m][3])
        Exp1s.SetPointError(n, 0., 0., val[m][3]-val[m][2], val[m][4]-val[m][3])
        Exp2s.SetPoint(n, m, val[m][3])
        Exp2s.SetPointError(n, 0., 0., val[m][3]-val[m][1], val[m][5]-val[m][3])
        if len(val[m]) > 6: Sign.SetPoint(n, m, val[m][6])
        if len(val[m]) > 7: pVal.SetPoint(n, m, val[m][7])
        if len(val[m]) > 8: Best.SetPoint(n, m, val[m][8])
        if len(val[m]) > 10: Best.SetPointError(n, 0., 0., abs(val[m][9]), val[m][10])


    for t in THEORY:
        Theory[t] = TGraphAsymmErrors()
        for m in sorted(HVT[t]['Z']['XS'].keys()):
            if m < mass[0] or m > mass[-1]: continue
            XsW, XsW_Up, XsW_Down, XsZ, XsZ_Up, XsZ_Down = 0., 0., 0., 0., 0., 0.
            XsZ = 1000.*HVT[t]['Z']['XS'][m]*HVT[t]['Z']['BR'][m]
            XsZ_Up = XsZ*(1.+math.hypot(HVT[t]['Z']['QCD'][m][0]-1., HVT[t]['Z']['PDF'][m][0]-1.))
            XsZ_Down = XsZ*(1.-math.hypot(1.-HVT[t]['Z']['QCD'][m][0], 1.-HVT[t]['Z']['PDF'][m][0]))

            n = Theory[t].GetN()
            Theory[t].SetPoint(n, m, XsW+XsZ)
            Theory[t].SetPointError(n, 0., 0., (XsW-XsW_Down)+(XsZ-XsZ_Down), (XsW_Up-XsW)+(XsZ_Up-XsZ))

            Theory[t].SetLineColor(theoryLineColor[t])
            Theory[t].SetFillColor(theoryFillColor[t])
            Theory[t].SetFillStyle(theoryFillStyle[t])
            Theory[t].SetLineWidth(2)
            #Theory[t].SetLineStyle(7)


    Exp2s.SetLineWidth(2)
    Exp2s.SetLineStyle(1)
    Obs0s.SetLineWidth(3)
    Obs0s.SetMarkerStyle(0)
    Obs0s.SetLineColor(1)
    Exp0s.SetLineStyle(2)
    Exp0s.SetLineWidth(3)
    Exp1s.SetFillColor(417) #kGreen+1
    Exp1s.SetLineColor(417) #kGreen+1
    Exp2s.SetFillColor(800) #kOrange
    Exp2s.SetLineColor(800) #kOrange
    Exp2s.GetXaxis().SetTitle("m_{"+particleP+"} (GeV)")
    Exp2s.GetXaxis().SetTitleSize(Exp2s.GetXaxis().GetTitleSize()*1.25)
    Exp2s.GetXaxis().SetNoExponent(True)
    Exp2s.GetXaxis().SetMoreLogLabels(True)
    Exp2s.GetYaxis().SetTitle("#sigma("+particleP+") #bf{#it{#Beta}}("+particleP+" #rightarrow "+particle+"h) (fb)")
    Exp2s.GetYaxis().SetTitleOffset(1.5)
    Exp2s.GetYaxis().SetNoExponent(True)
    Exp2s.GetYaxis().SetMoreLogLabels()

    Sign.SetLineWidth(2)
    Sign.SetLineColor(629)
    Sign.GetXaxis().SetTitle("m_{"+particleP+"} (GeV)")
    Sign.GetXaxis().SetTitleSize(Sign.GetXaxis().GetTitleSize()*1.1)
    Sign.GetYaxis().SetTitle("Significance")

    pVal.SetLineWidth(2)
    pVal.SetLineColor(629)
    pVal.GetXaxis().SetTitle("m_{"+particleP+"} (GeV)")
    pVal.GetXaxis().SetTitleSize(pVal.GetXaxis().GetTitleSize()*1.1)
    pVal.GetYaxis().SetTitle("local p-Value")

    Best.SetLineWidth(2)
    Best.SetLineColor(629)
    Best.SetFillColor(629)
    Best.SetFillStyle(3003)
    Best.GetXaxis().SetTitle("m_{"+particleP+"} (GeV)")
    Best.GetXaxis().SetTitleSize(Best.GetXaxis().GetTitleSize()*1.1)
    Best.GetYaxis().SetTitle("Best Fit (pb)")



    c1 = TCanvas("c1", "Exclusion Limits", 800, 600)
    c1.cd()
    #SetPad(c1.GetPad(0))
    c1.GetPad(0).SetTopMargin(0.06)
    c1.GetPad(0).SetRightMargin(0.05)
    c1.GetPad(0).SetLeftMargin(0.12)
    c1.GetPad(0).SetTicks(1, 1)
    #c1.GetPad(0).SetGridx()
    #c1.GetPad(0).SetGridy()
    c1.GetPad(0).SetLogy()
    Exp2s.Draw("A3")
    Exp1s.Draw("SAME, 3")
    for t in THEORY:
        Theory[t].Draw("SAME, L3")
        Theory[t].Draw("SAME, L3X0Y0")
    Exp0s.Draw("SAME, L")
    if not options.blind: Obs0s.Draw("SAME, L")
    #setHistStyle(Exp2s)
    Exp2s.GetXaxis().SetTitleSize(0.050)
    Exp2s.GetYaxis().SetTitleSize(0.050)
    Exp2s.GetXaxis().SetLabelSize(0.045)
    Exp2s.GetYaxis().SetLabelSize(0.045)
    Exp2s.GetXaxis().SetTitleOffset(0.90)
    Exp2s.GetYaxis().SetTitleOffset(1.25)
    Exp2s.GetYaxis().SetMoreLogLabels(True)
    Exp2s.GetYaxis().SetNoExponent(True)
    Exp2s.GetYaxis().SetRangeUser(0.01, 5.e3)
    #else: Exp2s.GetYaxis().SetRangeUser(0.1, 1.e2)
    Exp2s.GetXaxis().SetRangeUser(mass[0], min(mass[-1], MAXIMUM[channel] if channel in MAXIMUM else 1.e6))
    drawAnalysis(channel)
    drawRegion(channel, True)
    drawCMS(LUMI, YEAR, "Preliminary") #Preliminary

    # legend
    top = 0.9
    nitems = 4 + len(THEORY)

    leg = TLegend(0.55, top-nitems*0.3/5., 0.98, top)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0) #1001
    leg.SetFillColor(0)
    leg.SetHeader("95% CL upper limits")
    leg.AddEntry(Obs0s, "Observed", "l")
    leg.AddEntry(Exp0s, "Expected", "l")
    leg.AddEntry(Exp1s, "#pm 1 std. deviation", "f")
    leg.AddEntry(Exp2s, "#pm 2 std. deviation", "f")
    for t in THEORY: leg.AddEntry(Theory[t], theoryLabel[t], "fl")
    leg.Draw()
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)
    latex.SetTextFont(42)
    #latex.DrawLatex(0.66, leg.GetY1()-0.045, particleP+" #rightarrow "+particle+"h")

    leg2 = TLegend(0.12, 0.225-2*0.25/5., 0.65, 0.225)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0) #1001
    leg2.SetFillColor(0)
    c1.GetPad(0).RedrawAxis()
    """
    if True and channel.endswith('sl'):
        mass, val = fillValues("./combine/alpha/X"+particle+"Hbonly_M%d.txt")
        Exp, Obs = TGraphAsymmErrors(), TGraphAsymmErrors()
        for i, m in enumerate(mass):
            if not m in val: continue
            Exp.SetPoint(Exp.GetN(), m, val[m][3])
            Obs.SetPoint(Obs.GetN(), m, val[m][0])
        Exp.SetLineWidth(3)
        Exp.SetLineColor(602) #602
        Exp.SetLineStyle(5)
        Obs.SetLineWidth(3)
        Obs.SetLineColor(602)
        Exp.Draw("SAME, L")
        Obs.Draw("SAME, L")
        
        mass15, val15 = fillValues("./combine/Vh_2015/X"+particle+"h_M%d.txt")
        Exp15, Obs15 = TGraphAsymmErrors(), TGraphAsymmErrors()
        for i, m in enumerate(mass15):
            if not m in val: continue
            Exp15.SetPoint(Exp15.GetN(), m, val15[m][3]*multF*Theory['B3'].GetY()[i]*(0.625 if particle=='V' and m>3000 else 1.))
            Obs15.SetPoint(Obs15.GetN(), m, val15[m][0]*multF*Theory['B3'].GetY()[i]*(0.625 if particle=='V' and m>3000 else 1.))
        Exp15.SetLineWidth(3)
        Exp15.SetLineColor(856) #602
        Exp15.SetLineStyle(6)
        Obs15.SetLineWidth(3)
        Obs15.SetLineColor(856)
        Exp15.Draw("SAME, L")
        #Obs15.Draw("SAME, L")
        
        leg2.AddEntry(Exp, "1+2 b-tag", "l")
    """
    if True and channel=='AZh':
        massLL, valLL = fillValues("./combine/AZh/AZhll_M%d.txt")
        ExpLL, ObsLL = TGraphAsymmErrors(), TGraphAsymmErrors()
        for i, m in enumerate(massLL):
            if not m in val: continue
            ExpLL.SetPoint(ExpLL.GetN(), m, valLL[m][3]*multF)
            ObsLL.SetPoint(ObsLL.GetN(), m, valLL[m][0]*multF)
        ExpLL.SetLineWidth(3)
        ExpLL.SetLineColor(833) #602
        ExpLL.SetLineStyle(5)
        ObsLL.SetLineWidth(3)
        ObsLL.SetLineColor(833)
        ExpLL.Draw("SAME, L")
        #ObsLL.Draw("SAME, L")

        massNN, valNN = fillValues("./combine/AZh/AZhnn_M%d.txt")
        ExpNN, ObsNN = TGraphAsymmErrors(), TGraphAsymmErrors()
        for i, m in enumerate(massNN):
            if not m in val: continue
            ExpNN.SetPoint(ExpNN.GetN(), m, valNN[m][3]*multF)
            ObsNN.SetPoint(ObsNN.GetN(), m, valNN[m][0]*multF)
        ExpNN.SetLineWidth(3)
        ExpNN.SetLineColor(855) #602
        ExpNN.SetLineStyle(6)
        ObsNN.SetLineWidth(3)
        ObsNN.SetLineColor(855)
        ExpNN.Draw("SAME, L")
        #ObsNN.Draw("SAME, L")

        leg2.AddEntry(ExpLL, "Expected, A #rightarrow Zh #rightarrow llb#bar{b}", "l")
        leg2.AddEntry(ExpNN, "Expected, A #rightarrow Zh #rightarrow #nu#nub#bar{b}", "l")
    
    if method=='combo':
        massAH, valAH = fillValues("./combine/dijet/X"+particle+"Hah_M%d.txt")
        ExpAH = TGraphAsymmErrors()
        for i, m in enumerate(massAH):
            if not m in val: continue
            ExpAH.SetPoint(ExpAH.GetN(), m, valAH[m][3]*multF)
        ExpAH.SetLineWidth(2)
        ExpAH.SetLineColor(602) #602
        ExpAH.SetLineStyle(4)
        ExpAH.Draw("SAME, L")

        massSL, valSL = fillValues("./combine/alpha/X"+particle+"Hsl_M%d.txt")
        ExpSL = TGraphAsymmErrors()
        for i, m in enumerate(massSL):
            if not m in val: continue
            ExpSL.SetPoint(ExpSL.GetN(), m, valSL[m][3]*multF)
        ExpSL.SetLineWidth(3)
        ExpSL.SetLineColor(860-9) #602
        ExpSL.SetLineStyle(7)
        ExpSL.Draw("SAME, L")

        leg2.AddEntry(ExpAH, "B2G-17-002", "l")
        leg2.AddEntry(ExpSL, "B2G-17-004", "l")

    leg2.Draw()
    if not options.blind: Obs0s.Draw("SAME, L")
    c1.GetPad(0).Update()

    if not gROOT.IsBatch(): raw_input("Press Enter to continue...")

    c1.Print("plotsLimit/Exclusion/"+channel+suffix+".png")
    c1.Print("plotsLimit/Exclusion/"+channel+suffix+".pdf")
    if 'ah' in channel or 'sl' in channel:
        c1.Print("plotsLimit/Exclusion/"+channel+suffix+".C")
        c1.Print("plotsLimit/Exclusion/"+channel+suffix+".root")

    for t in THEORY:
        print "Model", t, ":",
        for m in range(mass[0], mass[-1], 1):
            if not (Theory[t].Eval(m) > Obs0s.Eval(m)) == (Theory[t].Eval(m+1) > Obs0s.Eval(m+1)): print m,
        print ""

    print "p1s[",
    for i in range(Exp0s.GetN()):
        print Exp0s.GetY()[i]+Exp1s.GetErrorYhigh(i), ",",
    print "],"
    print "m1s[",
    for i in range(Exp0s.GetN()):
        print Exp0s.GetY()[i]-Exp1s.GetErrorYlow(i), ",",
    print "],"
    print "[",
    for i in range(Exp0s.GetN()):
        print Exp0s.GetY()[i], ",",
    print "]"

    #if not 'ah' in channel and not 'sl' in channel: return

    # ---------- Significance ----------
    c2 = TCanvas("c2", "Significance", 800, 600)
    c2.cd()
    c2.GetPad(0).SetTopMargin(0.06)
    c2.GetPad(0).SetRightMargin(0.05)
    c2.GetPad(0).SetTicks(1, 1)
    c2.GetPad(0).SetGridx()
    c2.GetPad(0).SetGridy()
    Sign.GetYaxis().SetRangeUser(0., 5.)
    Sign.Draw("AL3")
    drawCMS(LUMI, YEAR, "Preliminary")
    drawAnalysis(channel[1:3])
    c2.Print("plotsLimit/Significance/"+channel+suffix+".png")
    c2.Print("plotsLimit/Significance/"+channel+suffix+".pdf")
#    c2.Print("plotsLimit/Significance/"+channel+suffix+".root")
#    c2.Print("plotsLimit/Significance/"+channel+suffix+".C")

    # ---------- p-Value ----------
    c3 = TCanvas("c3", "p-Value", 800, 600)
    c3.cd()
    c3.GetPad(0).SetTopMargin(0.06)
    c3.GetPad(0).SetRightMargin(0.05)
    c3.GetPad(0).SetTicks(1, 1)
    c3.GetPad(0).SetGridx()
    c3.GetPad(0).SetGridy()
    c3.GetPad(0).SetLogy()
    pVal.Draw("AL3")
    pVal.GetYaxis().SetRangeUser(2.e-7, 0.5)

    ci = [1., 0.317310508, 0.045500264, 0.002699796, 0.00006334, 0.000000573303, 0.000000001973]
    line = TLine()
    line.SetLineColor(922)
    line.SetLineStyle(7)
    text = TLatex()
    text.SetTextColor(922)
    text.SetTextSize(0.025)
    text.SetTextAlign(12)
    for i in range(1, len(ci)-1):
        line.DrawLine(pVal.GetXaxis().GetXmin(), ci[i]/2, pVal.GetXaxis().GetXmax(), ci[i]/2);
        text.DrawLatex(pVal.GetXaxis().GetXmax()*1.01, ci[i]/2, "%d #sigma" % i);

    drawCMS(LUMI, YEAR, "Preliminary")
    drawAnalysis(channel[1:3])
    c3.Print("plotsLimit/pValue/"+channel+suffix+".png")
    c3.Print("plotsLimit/pValue/"+channel+suffix+".pdf")
#    c3.Print("plotsLimit/pValue/"+channel+suffix+".root")
#    c3.Print("plotsLimit/pValue/"+channel+suffix+".C")

    # --------- Best Fit ----------
    c4 = TCanvas("c4", "Best Fit", 800, 600)
    c4.cd()
    c4.GetPad(0).SetTopMargin(0.06)
    c4.GetPad(0).SetRightMargin(0.05)
    c4.GetPad(0).SetTicks(1, 1)
    c4.GetPad(0).SetGridx()
    c4.GetPad(0).SetGridy()
    Best.Draw("AL3")
    drawCMS(LUMI, YEAR, "Preliminary")
    drawAnalysis(channel[1:3])
    c4.Print("plotsLimit/BestFit/"+channel+suffix+".png")
    c4.Print("plotsLimit/BestFit/"+channel+suffix+".pdf")
#    c4.Print("plotsLimit/BestFit/"+channel+suffix+".root")
#    c4.Print("plotsLimit/BestFit/"+channel+suffix+".C")

    if not gROOT.IsBatch(): raw_input("Press Enter to continue...")

    if 'ah' in channel:
        outFile = TFile("bands.root", "RECREATE")
        outFile.cd()
        pVal.Write("graph")
        Best.Write("best")
        outFile.Close()


def limit2HDM():
    global signals
    signals = range(800, 2000+1, 50)
    multF = HTOBB
    THEORY = ['T1', 'T2']
    
    mass, val = fillValues("./combine/AZh/AZh_M%d.txt")
    Obs0s = TGraph()
    Exp0s = TGraph()
    Exp1s = TGraphAsymmErrors()
    Exp2s = TGraphAsymmErrors()
    
    massB, valB = fillValues("./combine/BBAZh/BBAZh_M%d.txt")
    Obs0sB = TGraph()
    Exp0sB = TGraph()
    Exp1sB = TGraphAsymmErrors()
    Exp2sB = TGraphAsymmErrors()
    
    for i, m in enumerate(mass):
        if not m in val:
            print "Key Error:", m, "not in value map"
            continue

        n = Exp0s.GetN()
        Obs0s.SetPoint(n, m, val[m][0]*multF)
        Exp0s.SetPoint(n, m, val[m][3]*multF)
        Exp1s.SetPoint(n, m, val[m][3]*multF)
        Exp1s.SetPointError(n, 0., 0., val[m][3]*multF-val[m][2]*multF, val[m][4]*multF-val[m][3]*multF)
        Exp2s.SetPoint(n, m, val[m][3]*multF)
        Exp2s.SetPointError(n, 0., 0., val[m][3]*multF-val[m][1]*multF, val[m][5]*multF-val[m][3]*multF)
        
        Obs0sB.SetPoint(n, m, valB[m][0]*multF)
        Exp0sB.SetPoint(n, m, valB[m][3]*multF)
        Exp1sB.SetPoint(n, m, valB[m][3]*multF)
        Exp1sB.SetPointError(n, 0., 0., valB[m][3]*multF-valB[m][2]*multF, valB[m][4]*multF-valB[m][3]*multF)
        Exp2sB.SetPoint(n, m, valB[m][3]*multF)
        Exp2sB.SetPointError(n, 0., 0., valB[m][3]*multF-valB[m][1]*multF, valB[m][5]*multF-valB[m][3]*multF)

    col = 629
    Exp2s.SetLineWidth(2)
    Exp2s.SetLineStyle(1)
    Obs0s.SetLineWidth(3)
    Obs0s.SetMarkerStyle(0)
    Obs0s.SetLineColor(1)
    Exp0s.SetLineStyle(2)
    Exp0s.SetLineWidth(3)
    Exp0s.SetLineColor(1)
#    Exp1s.SetFillColorAlpha(col, 0.4) #kGreen+1
#    Exp1s.SetLineColorAlpha(col, 0.4)
#    Exp2s.SetFillColorAlpha(col, 0.2) #kOrange
#    Exp2s.SetLineColorAlpha(col, 0.2)
    Exp1s.SetFillColor(417)
    Exp1s.SetLineColor(417)
    Exp2s.SetFillColor(800)
    Exp2s.SetLineColor(800)
    
    colB = 922
    Exp2sB.SetLineWidth(2)
    Obs0sB.SetLineStyle(9)
    Obs0sB.SetLineWidth(3)
    Obs0sB.SetMarkerStyle(0)
    Obs0sB.SetLineColor(colB)
    Exp0sB.SetLineStyle(8)
    Exp0sB.SetLineWidth(3)
    Exp0sB.SetLineColor(colB)
    Exp1sB.SetFillColorAlpha(colB, 0.4) #kGreen+1
    Exp1sB.SetLineColorAlpha(colB, 0.4)
    Exp2sB.SetFillColorAlpha(colB, 0.2) #kOrange
    Exp2sB.SetLineColorAlpha(colB, 0.2)
    
    Exp2s.GetXaxis().SetTitle("m_{A} (GeV)")
    Exp2s.GetXaxis().SetTitleSize(Exp2s.GetXaxis().GetTitleSize()*1.25)
    Exp2s.GetXaxis().SetNoExponent(True)
    Exp2s.GetXaxis().SetMoreLogLabels(True)
    Exp2s.GetYaxis().SetTitle("#sigma(A) #bf{#it{#Beta}}(A #rightarrow Zh) (fb)")
    Exp2s.GetYaxis().SetTitleOffset(1.5)
    Exp2s.GetYaxis().SetNoExponent(True)
    Exp2s.GetYaxis().SetMoreLogLabels()
    
    Theory = {}
    for t in THEORY:
        Theory[t] = TGraphAsymmErrors()
        for m in sorted(THDM[t]['ggA'].keys()):
            if m < mass[0] or m > mass[-1]: continue
            Xs, Xs_Up, Xs_Down = 0., 0., 0.
            Xs = THDM[t]['ggA'][m]
            Xs_Up = Xs*(1.+math.sqrt((THDM['PDF']['ggA'][m][0]-1.)**2 + (THDM['QCD']['ggA'][m][0]-1.)**2))
            Xs_Down = Xs*(1.-math.sqrt((1.-THDM['PDF']['ggA'][m][1])**2 + (1.-THDM['QCD']['ggA'][m][1])**2))
            n = Theory[t].GetN()
            Theory[t].SetPoint(n, m, Xs)
            Theory[t].SetPointError(n, 0., 0., (Xs-Xs_Down), (Xs_Up-Xs))

        Theory[t].SetLineColor(theoryLineColor[t])
        Theory[t].SetFillColor(theoryFillColor[t])
        Theory[t].SetFillStyle(theoryFillStyle[t])
        Theory[t].SetLineWidth(2)
            #Theory[t].SetLineStyle(7)

    c1 = TCanvas("c1", "Exclusion Limits", 800, 600)
    c1.cd()
    #SetPad(c1.GetPad(0))
    c1.GetPad(0).SetTopMargin(0.06)
    c1.GetPad(0).SetRightMargin(0.05)
    c1.GetPad(0).SetLeftMargin(0.12)
    c1.GetPad(0).SetTicks(1, 1)
    c1.GetPad(0).SetLogy()
    Exp2s.Draw("A3")
    Exp1s.Draw("SAME, 3")
    Exp0s.Draw("SAME, L")
#    Exp2sB.Draw("SAME, 3")
#    Exp1sB.Draw("SAME, 3")
    Exp0sB.Draw("SAME, L")
    if not options.blind:
        Obs0s.Draw("SAME, L")
        Obs0sB.Draw("SAME, L")
    for t in THEORY:
        Theory[t].Draw("SAME, L3")
        Theory[t].Draw("SAME, L3X0Y0")
    #setHistStyle(Exp2s)
#    Exp2s.GetXaxis().SetTitleSize(0.045)
#    Exp2s.GetYaxis().SetTitleSize(0.04)
#    Exp2s.GetXaxis().SetLabelSize(0.04)
#    Exp2s.GetYaxis().SetLabelSize(0.04)
#    Exp2s.GetXaxis().SetTitleOffset(1)
#    Exp2s.GetYaxis().SetTitleOffset(1.25)
    Exp2s.GetXaxis().SetTitleSize(0.050)
    Exp2s.GetYaxis().SetTitleSize(0.050)
    Exp2s.GetXaxis().SetLabelSize(0.045)
    Exp2s.GetYaxis().SetLabelSize(0.045)
    Exp2s.GetXaxis().SetTitleOffset(0.90)
    Exp2s.GetYaxis().SetTitleOffset(1.25)
    Exp2s.GetYaxis().SetMoreLogLabels(True)
    Exp2s.GetYaxis().SetNoExponent(True)
    Exp2s.GetYaxis().SetRangeUser(0.5, 1.e3)
    Exp2s.GetXaxis().SetRangeUser(mass[0], mass[-1])
    drawAnalysis('AZh')
    drawRegion('AZHsl', True)
    drawCMS(LUMI, YEAR, "Preliminary") #Preliminary

    # legend
    leg = TLegend(0.6, 0.90, 0.99, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0) #1001
    leg.SetFillColor(0)
    leg.SetHeader("95% CL upper limits")
    leg.AddEntry(None, "gg #rightarrow A #rightarrow Zh", "") #"95% CL upper limits"
    leg.AddEntry(Obs0s,  "Observed", "l")
    leg.AddEntry(Exp0s,  "Expected", "l")
    leg.AddEntry(Exp1s, "#pm 1 std. deviation", "f")
    leg.AddEntry(Exp2s, "#pm 2 std. deviation", "f")
    leg.AddEntry(None, "", "")
    leg.AddEntry(None, "bbA #rightarrow Zh", "")
    leg.AddEntry(Obs0sB,  "Observed", "l")
    leg.AddEntry(Exp0sB,  "Expected", "l")
    leg.SetY1(leg.GetY2()-leg.GetNRows()*0.045)
    leg.Draw()
    
    
#    latex = TLatex()
#    latex.SetNDC()
#    latex.SetTextSize(0.040)
#    latex.SetTextFont(42)
#    latex.DrawLatex(0.65, leg.GetY1()-0.045, "cos(#beta-#alpha)=0.25, tan(#beta)=1")

#    legB = TLegend(0.12, 0.4-4*0.3/5., 0.65, 0.4)
    legB = TLegend(0.15, 0.27, 0.68, 0.27)
    legB.SetBorderSize(0)
    legB.SetFillStyle(0) #1001
    legB.SetFillColor(0)
    for t in THEORY: legB.AddEntry(Theory[t], theoryLabel[t], "fl")
    legB.AddEntry(None, "cos(#beta-#alpha)=0.25, tan(#beta)=1", "")
    legB.SetY1(legB.GetY2()-legB.GetNRows()*0.045)
    legB.Draw()
    
    c1.GetPad(0).RedrawAxis()
    leg.Draw()

    c1.Update()
    
    if not gROOT.IsBatch(): raw_input("Press Enter to continue...")

    c1.Print("plotsLimit/Exclusion/THDM.png")
    c1.Print("plotsLimit/Exclusion/THDM.pdf")



def limitCompare(method):
    signal, particle, particleP = "XZH", "Z", "Z'"

    channels = ["0bonly","bbonly"]
    colors = [1, 610+4]#, 632, 800+7, 800, 416+1, 860+10, 600, 616, 921, 922]
    masses, vals, graphs = {}, {}, {}
    for j, c in enumerate(channels):
        masses[c], vals[c] = fillValues("./combine/" + method + "/"+signal+c+"_M%d.txt")
        graphs[c] = TGraph()
        n = 0
        #print vals[c]
        for i, m in enumerate(masses[c]):
            #if not signals[i] >= 1000: continue
            graphs[c].SetPoint(n, m, vals[c][m][3])
            n = n + 1
        graphs[c].SetLineColor(colors[j])
        graphs[c].SetLineWidth(3)
        graphs[c].Draw("SAME, L")
    comp = TCanvas("comp", "Exclusion Limits", 800, 600)
    comp.cd()
    #c1.SetPad(c1.GetPad(0))
    comp.GetPad(0).SetTopMargin(0.06)
    comp.GetPad(0).SetRightMargin(0.05)
    comp.GetPad(0).SetTicks(1, 1)
    comp.GetPad(0).SetLogy()
    for c in channels:
        graphs[c].Draw("SAME, L")
    graphs[channels[0]].GetXaxis().SetRangeUser(800., 5000.)
    graphs[channels[0]].GetYaxis().SetRangeUser(0.25, 2.5e4)
    graphs[channels[0]].GetXaxis().SetTitle("m_{"+particleP+"} (GeV)")
    graphs[channels[0]].GetYaxis().SetTitle("#sigma("+particleP+") #bf{#it{#Beta}}("+particleP+" #rightarrow "+particle+"H) (fb)")
    drawAnalysis(signal)
    #drawRegion(signal, True)
    drawCMS(LUMI, YEAR, "Preliminary")

    # legend
    top = 0.9
    leg = TLegend(0.4, top-len(channels)*0.2/5., 0.99, top)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0) #1001
    leg.SetFillColor(0)
    leg.SetHeader("95% CL expected limits")
    for c in channels:
        leg.AddEntry(graphs[c],  getChannel(c), "l")
    leg.Draw()

    comp.Update()
    comp.Print("plotsLimit/Multi.png")
    comp.Print("plotsLimit/Multi.pdf")


if options.all:
    for c in ["XZHnnbb", "XZHnn0b", "XZHeebb", "XZHee0b", "XZHmmbb", "XZHmm0b", "XZHVBFnnbbVBF", "XZHVBFnn0bVBF", "XZHVBFeebbVBF", "XZHVBFee0bVBF", "XZHVBFmmbbVBF", "XZHVBFmm0bVBF", "XZHsl","XZHslbb", "XZHsl0b", "XZHslbbVBF", "XZHsl0bVBF", "XZHslcomb"]:
        limit(options.method, c)
else:
    if not options.method=="2HDM": 
        limit(options.method, options.category)
        #limitCompare(options.method)
    else: limit2HDM()
  
