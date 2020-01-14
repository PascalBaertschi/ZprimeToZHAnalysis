#! /usr/bin/env python

import os, multiprocessing
import copy
import math
from array import array
from ROOT import ROOT, gROOT, gStyle, gRandom, TSystemDirectory
from ROOT import TFile, TChain, TTree, TCut, TH1, TH1F, TH2F, THStack, TGraph
from ROOT import TStyle, TCanvas, TPad
from ROOT import TLegend, TLatex, TText, TLine

from rooUtils import getCrossSection, getColor
from utils import drawCMS, drawRegion, drawAnalysis, drawChi2, setTopPad, setBotPad, setBotStyle, setHistStyle, convertHistToGraph, makeResidHist
from hvtXs import HVT

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-a", "--all", action="store_true", default=False, dest="all")
parser.add_option("-b", "--bash", action="store_true", default=False, dest="bash")
parser.add_option("-d", "--data", action="store_true", default=False, dest="data")
parser.add_option("-c", "--category", action="store", type="string", dest="category", default="")
parser.add_option("-2", "--category2", action="store", type="string", dest="category2", default="")
parser.add_option("-s", "--signal", action="store", type="string", dest="signal", default="XVH")
parser.add_option("-m", "--mass", action="store", type=int, dest="mass", default=2000)
parser.add_option("-f", "--fileName", action="store", type="string", dest="fileName", default="combine/test/mlfit_XVHsl_M2000.root")
parser.add_option("", "--fileName2", action="store", type="string", dest="fileName2", default="") #combine/test/mlfit_monoHnn_MZ3000_MA300.root
parser.add_option("", "--fileName3", action="store", type="string", dest="fileName3", default="")
parser.add_option("-o", "--outName", action="store", type="string", dest="outName", default="prePostFit")
parser.add_option("-v", "--verbose", action="store_true", default=False, dest="verbose")
(options, args) = parser.parse_args()
if options.bash: gROOT.SetBatch(True)

fileName = options.fileName
outName  = options.outName
readData = options.data
signal   = options.signal
signals  = [x for x in [options.fileName2, options.fileName3] if len(x)>0]
mass     = options.mass
VERBOSE  = options.verbose
RATIO    = 4
LUMI     = 35867

channelList = ['wrhpb', 'wrhpbb', 'wrlpb', 'wrlpbb', 'zrhpb', 'zrhpbb', 'zrlpb', 'zrlpbb', 'nnb', 'nnbb', 'enb', 'enbb', 'mnb', 'mnbb', 'eeb', 'mmb', 'eebb', 'mmbb']
signalList = ['XWH', 'XZH', 'XVH']

# Format:
#0 : sys name
#1 : b-only DX
#2 : b-only dDX
#3 : s+b DX
#4 : s+b dDX
#5 : rho

gStyle.SetOptStat(0)
gStyle.SetErrorX(0)


def plotPrePost(category, category2):
    if len(category)==0:
        print "Please select a category with the -c option"
        exit()

    is2 = len(category2)>0
    if 'l' in category: category, category2, is2 = category.replace('l', 'e'), category.replace('l', 'm'), True
    isAH = len(category) >= 5
    X_name = "JJ_mass" if isAH else "VH_mass"
    signalName = signal + category + "_M%d" % mass
    data = ["data_obs"]
    back = ["Bkg_"+category] if isAH else ["VV_"+category, "Top_"+category, "Vjets_"+category]
    sign = [signalName] if len(signal)>0 else []
    
    fitr = ["total", "total_signal", "total_background"]

    labels = {"data_obs" : "Data", "VV_"+category : "VV, Vh", "Top_"+category : "t#bar{t}, t+X", "Vjets_"+category : "V+jets", "Bkg_"+category : "Bkg fit", signalName : "Signal"}
    if 'ee' in category or 'mm' in category: labels["Vjets_"+category] = "Z(ll)+jets"
    elif 'en' in category or 'mn' in category: labels["Vjets_"+category] = "W(l#nu)+jets"
    elif 'nn' in category: labels["Vjets_"+category] = "Z(#nu#nu),W(l#nu)+jets"

    lastBin = 3500. if 'ee' in category or 'mm' in category or 'nn' in category else 4500.
    
#    xs = 0.
#    if 'XWH' in signal or 'XVH' in signal: xs += HVT['B3']['W']['XS'][mass]*HVT['B3']['W']['BR'][mass]
#    if 'XZH' in signal or 'XVH' in signal: xs += HVT['B3']['Z']['XS'][mass]*HVT['B3']['Z']['BR'][mass]
    xs = getCrossSection(signal, category, mass)
    #if len(signal2)>0: xs[signal2] = getCrossSection(signal2, category, 0)
    
    histData, histPre, histPost, graphPre, graphPost = None, {}, {}, {}, {}

    if readData:
        dataFile = TFile("workspace/"+category+".root", "READ")
        workspace = dataFile.Get("VH_2016")
        variable = workspace.var(X_name)
        dataset = workspace.data(data[0])
        if is2:
            dataFile2 = TFile("workspace/"+category2+".root", "READ")
            workspace2 = dataFile2.Get("VH_2016")
            dataset2 = workspace2.data(data[0])
            dataset.append(dataset2)
        data2 = dataset.createHistogram(variable, variable, variable.getBinning().numBins()/10, 1)
        data1 = data2.ProjectionX()
        data1.SetMarkerStyle(20)
        data1.SetMarkerSize(1.25)
        data1.SetLineColor(1)
        histData = data1
        graphData = convertHistToGraph(histData, True)
        width = histData.GetXaxis().GetBinWidth(1)
        

    inFile = TFile(options.fileName, "READ")
    if inFile==None:
        print "File", options.fileName, "not found"
        return
    if not inFile.GetDirectory("shapes_prefit/"+category):
        print "Category", category, "not recognized"
        return
    
    for i, h in enumerate(back+sign+fitr):
        histPre[h] = inFile.Get("shapes_prefit/"+category+"/"+h)
        if is2: histPre[h].Add( inFile.Get("shapes_prefit/"+category2+"/"+h.replace(category, category2)) )
        histPre[h].SetName(h+'_pre')
        histPre[h].SetLineColor(getColor(h, category))
        histPre[h].SetFillColor(getColor(h, category))
        histPre[h].SetLineWidth(3)
        if h in back: histPre[h].SetFillStyle(1001)
        elif h in fitr: histPre[h].SetFillStyle(3002)
        elif h in sign:
            histPre[h].SetTitle("m_{%s'} = %d GeV" % (signal[1], mass))
            histPre[h].SetOption("HVT model B g_{V}=3")
            histPre[h].SetFillStyle(1)
            histPre[h].SetLineStyle(3)
            histPre[h].SetLineWidth(6)
            if xs>0.: histPre[h].Scale(xs*1000.)
    #    histPre[h].Rebin(10)
        if readData:
            histPre[h].Scale(width)
            histPre[h].GetYaxis().SetTitle("Events / ( %d GeV )" % width)
        if 'nn' in category:
            histPre[h].GetXaxis().SetTitle(histPre[h].GetXaxis().GetTitle().replace('m_{VH}', 'm^{T}_{VH}'))
        

    for i, h in enumerate(back+sign+fitr):
        histPost[h] = inFile.Get("shapes_fit_b/"+category+"/"+h)
        if is2: histPost[h].Add( inFile.Get("shapes_fit_b/"+category2+"/"+h.replace(category, category2)) )
        histPost[h].SetName(h+'_post')
        histPost[h].SetLineColor(getColor(h, category))
        histPost[h].SetFillColor(getColor(h, category))
        histPost[h].SetLineWidth(3)
        if h in back: histPost[h].SetFillStyle(1001)
        elif h in fitr: histPost[h].SetFillStyle(3002)
        elif h in sign:
            histPost[h].SetTitle("m_{%s'} = %d GeV" % (signal[1], mass))
            histPost[h].SetOption("HVT model B g_{V}=3")
            histPost[h].SetFillStyle(1)
            histPost[h].SetLineStyle(3)
            histPost[h].SetLineWidth(6)
            if xs>0.: histPost[h].Scale(xs*1000.)
    #    histPost[h].Rebin(10)
        if readData:
            histPost[h].Scale(width)
            histPost[h].GetYaxis().SetTitle("Events / ( %d GeV )" % width)
        if 'nn' in category:
            histPost[h].GetXaxis().SetTitle(histPost[h].GetXaxis().GetTitle().replace('m_{VH}', 'm^{T}_{VH}'))

    # Set errors ot zero to have smooth curves
    for i, h in enumerate(back+sign):
        for i in range(histPre[h].GetNbinsX()): histPre[h].SetBinError(i+1, 0.)
        for i in range(histPost[h].GetNbinsX()): histPost[h].SetBinError(i+1, 0.)

    stackPre = THStack("Pre", ";"+histPre['total'].GetXaxis().GetTitle()+";"+histPre['total'].GetYaxis().GetTitle())
    for i, s in enumerate(back): stackPre.Add(histPre[s])

    stackPost = THStack("Post", ";"+histPost['total'].GetXaxis().GetTitle()+";"+histPost['total'].GetYaxis().GetTitle())
    for i, s in enumerate(back): stackPost.Add(histPost[s])


    for i, h in enumerate(back):
        tmpPre = histPre[back[i]].Clone(back[i]+"_stack_pre")
        for j in range(i): tmpPre.Add(histPre[back[j]])
        graphPre[back[i]] = convertHistToGraph(tmpPre)
        tmpPost = histPost[back[i]].Clone(back[i]+"_stack_post")
        for j in range(i): tmpPost.Add(histPost[back[j]])
        graphPost[back[i]] = convertHistToGraph(tmpPost)

    # Additional signal, if present
    inFiles = {}
    for i, f in enumerate(signals): #combine/test/mlfit_monoHnn_MZ3000_MA300.root
        signalName2 = f.replace('combine/test/mlfit_', '').replace('.root', '')
        if 'AZh' in signalName2: signalName2 = signalName2.replace('AZh', 'AZh'+category)
        signalSmpl2 = signalName2.replace(category, '')
        signal2 = signalName2.split('_M')[0].replace(category, '')
        try:
            mass2, mass2A = int(signalName2.split('_M')[1]), 0
        except:
            mass2, mass2A = int(signalName2.split('_MZ')[1].split('_MA')[0]), int(signalName2.split('_MA')[1])
        #
        inFiles[signalName2] = TFile(f, "READ")
        histSign2 = inFiles[signalName2].Get("shapes_prefit/"+category+"/"+signalName2)
        if is2: histSign2.Add( inFiles[signalName2].Get("shapes_prefit/"+category2+"/"+signalName2.replace(category, category2)) )
        histSign2.SetName(signalName2+'_pre')
        if signal2.startswith('X'): histSign2.SetTitle("m_{V'} = %d GeV" % mass2)
        elif signal2.startswith('A'): histSign2.SetTitle("m_{A} = %d GeV" % mass2)
        elif not mass2==0: histSign2.SetTitle("m_{Z'} = %d GeV" % mass2)
        if i==0: histSign2.SetOption("Z'-2HDM\nm_{A}=300 GeV" if 'monoH' in signalName2 else "Type-II 2HDM\ncos(#beta-#alpha) = 0.25\ntan#beta = 1")
        histSign2.SetLineColor(getColor(signalName2, category))
        histSign2.SetFillColor(getColor(signalName2, category))
        histSign2.SetFillStyle(1)
        histSign2.SetLineStyle(5+i)
        histSign2.SetLineWidth(5)
        if readData: histSign2.Scale(width)
        xs2 = getCrossSection(signalName2, category, 0)
        if xs2>0.: histSign2.Scale(xs2*1000.)
        sign += [signalName2]
        histPre[signalName2] = histSign2
        histPost[signalName2] = histSign2
        

    leg = TLegend(0.6-0.005, 0.6, 0.925, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0) #1001
    leg.SetFillColor(0)
    if readData:
        leg.AddEntry(graphData, labels[data[0]], "pe")
    for i, s in reversed(list(enumerate(back))):
        leg.AddEntry(histPre[s], labels[s], "f")
    #for i, s in enumerate(sign):
    #    leg.AddEntry(histPre[s], labels[s], "fl")

    ###

    c1 = TCanvas("c1", "Pre-Post Fit", 800, 800)
    c1.Divide(1, 2)
    c1.cd(1)
    setTopPad(c1.GetPad(1), RATIO)
    setBotPad(c1.GetPad(2), RATIO)
    c1.GetPad(1).SetTopMargin(0.06)
    c1.GetPad(1).SetRightMargin(0.05)
    c1.GetPad(1).SetBottomMargin(0.01)
    c1.GetPad(1).SetTicks(1, 1)
    
    #histData.Draw("APE" if d==0 else "SAME, PE")
    #for i, h in enumerate(back): graphPost[h].Draw("ACL" if i==0 else "C")
    stackPost.Draw("C")
    setHistStyle(stackPost, 1.2)
    stackPost.SetMaximum(stackPost.GetMaximum()*5.)
    stackPost.SetMinimum(max(stackPost.GetMinimum(), 0.2))
    stackPost.GetXaxis().SetRangeUser(stackPost.GetXaxis().GetXmin(), lastBin)
    histPost['total'].SetLineWidth(1)
    histPost['total'].Draw("SAME, E3")
    histPre['total_background'].SetLineColor(921)
    histPre['total_background'].SetLineStyle(2)
    histPre['total_background'].SetLineWidth(3)
    histPre['total_background'].SetFillColor(1)
    histPre['total_background'].SetFillStyle(0)
    histPre['total_background'].Draw("SAME, HIST")
    for i, s in enumerate(sign):
        histPre[s].Draw("SAME, L")
    if readData:
        graphData.Draw("SAME, PE0")
    
    
    leg.AddEntry(histPost['total'], "Bkg. unc.", "f")
    leg.AddEntry(histPre['total_background'], "Pre-fit", "l")
    for i, s in enumerate(sign):
        for o in histPre[s].GetOption().split('\n'):
            if len(o)>0.: leg.AddEntry(None, o, "")
        leg.AddEntry(histPre[s], histPre[s].GetTitle(), "l")
    leg.SetY1(0.9-leg.GetNRows()*0.060)
    leg.Draw()
#    if len(sign)>0:
#        latex = TLatex()
#        latex.SetNDC()
#        latex.SetTextSize(0.045)
#        latex.SetTextFont(42)
#        latex.DrawLatex(0.67, leg.GetY1()-0.045, "HVT model B g_{V}=3")
    drawCMS(LUMI, "") #Preliminary
    drawRegion('XVH'+options.category, True)
    drawAnalysis(category)

    c1.GetPad(1).SetLogy()

    c1.cd(2)
    err = histPost['total'].Clone("BkgErr;")
    err.SetTitle("")
    err.Reset("MICES")
    err.GetYaxis().SetTitle("(N^{data}-N^{bkg})/#sigma")
    setBotStyle(err, 4+1)
    err.GetXaxis().SetTitleSize(0.16);
    #err.GetXaxis().SetTitleOffset(1.25);
    err.GetYaxis().SetTitleOffset(0.33);
    err.GetXaxis().SetRangeUser(err.GetXaxis().GetXmin(), lastBin)
    err.GetYaxis().SetRangeUser(-5., 5.)
    err.SetLineWidth(2)
    err.SetLineStyle(2)
    err.SetFillStyle(0)
    err.Draw("L") #"E2"
    if readData:
        pulls = makeResidHist(graphData, histPost['total'])
        #setBotStyle(pulls, RATIO, False)
        pulls.Draw("SAME, PE0")
        #drawRatio(hist['data_obs'], hist['BkgSum'])
        #drawStat(hist['data_obs'], hist['BkgSum'])
        chi2, nbins, npar = 0., 0, 0
        for i in range(0, pulls.GetN()):
            if graphData.GetY()[i] > 1.e-3:
                nbins = nbins + 1
                chi2 += pulls.GetY()[i]**2
        #drawChi2(chi2, nbins-npar, True)
        
    c1.Update()
    
    if not gROOT.IsBatch(): raw_input("Press Enter to continue...")

    
    c1.Print("plotsPrePost/BkgSR_"+options.category+".png")
    c1.Print("plotsPrePost/BkgSR_"+options.category+".pdf")

    ###

    if VERBOSE:
        c2 = TCanvas("c2", "Pre-Post Fit", 1200, 600)
        c2.Divide(2, 1)
        c2.cd(1)
        c2.GetPad(1).SetTopMargin(0.06)
        c2.GetPad(1).SetRightMargin(0.05)
        c2.GetPad(1).SetBottomMargin(0.10)
        c2.GetPad(1).SetTicks(1, 1)

        #histData.Draw("APE" if d==0 else "SAME, PE")
        #for i, h in enumerate(back): graphPre[h].Draw("FL" if i==0 else "SAME, FL")
        stackPre.Draw("HIST")
        histPre['total'].Draw("SAME, E2")
        stackPre.SetMaximum(stackPre.GetMaximum()*5.)
        stackPre.SetMinimum(max(stackPre.GetMinimum(), 0.01))
        histData.Draw("SAME, PE")

        leg.Draw()
        drawCMS(LUMI, "Preliminary")
        drawRegion('XVH'+category, True)
        drawAnalysis(category)

        c2.cd(2)
        c2.GetPad(2).SetTopMargin(0.06)
        c2.GetPad(2).SetRightMargin(0.05)
        c2.GetPad(2).SetBottomMargin(0.10)
        c2.GetPad(2).SetTicks(1, 1)

        #histData.Draw("APE" if d==0 else "SAME, PE")
        #for i, h in enumerate(back): graphPost[h].Draw("ACL" if i==0 else "C")
        stackPost.Draw("HIST")
        histPost['total'].Draw("SAME, E2")
        stackPost.SetMaximum(stackPost.GetMaximum()*5.)
        stackPost.SetMinimum(max(stackPost.GetMinimum(), 0.01))
        histData.Draw("SAME, PE")

        leg.Draw()
        drawCMS(LUMI, "Preliminary")
        drawRegion('XVH'+category, True)
        drawAnalysis(category)

        c2.GetPad(1).SetLogy()
        c2.GetPad(2).SetLogy()

        c2.Print("combine/test/"+signalName+"_prepost.png")
        c2.Print("combine/test/"+signalName+"_prepost.pdf")

        c2.Close()





if options.all:
    for c in channelList: plotPrePost(c)
else:
    plotPrePost(options.category, options.category2)

'''
source combineTest.sh datacards/alpha/XVHsl_M2000.txt
python plotPrePost.py -b -d -c nnb -s XVH -m 2000
python plotPrePost.py -b -d -c nnbb -s XVH -m 2000
python plotPrePost.py -b -d -c enb -s XVH -m 2000
python plotPrePost.py -b -d -c mnb -s XVH -m 2000
python plotPrePost.py -b -d -c lnb -s XVH -m 2000
python plotPrePost.py -b -d -c enbb -s XVH -m 2000
python plotPrePost.py -b -d -c mnbb -s XVH -m 2000
python plotPrePost.py -b -d -c lnbb -s XVH -m 2000
python plotPrePost.py -b -d -c eeb -s XVH -m 2000
python plotPrePost.py -b -d -c mmb -s XVH -m 2000
python plotPrePost.py -b -d -c llb -s XVH -m 2000
python plotPrePost.py -b -d -c eebb -s XVH -m 2000
python plotPrePost.py -b -d -c mmbb -s XVH -m 2000
python plotPrePost.py -b -d -c llbb -s XVH -m 2000

python plotPrePost.py -b -d -c nnb -s XVH -m 2000 --fileName2 combine/test/mlfit_monoHnnb_MZ1400_MA300.root --fileName3 combine/test/mlfit_monoHnnb_MZ3000_MA300.root
python plotPrePost.py -b -d -c nnbb -s XVH -m 2000 --fileName2 combine/test/mlfit_monoHnnbb_MZ1400_MA300.root --fileName3 combine/test/mlfit_monoHnnbb_MZ3000_MA300.root
python plotPrePost.py -b -d -c llb -s XVH -m 2000 --fileName2 combine/test/mlfit_AZh_M1000.root
python plotPrePost.py -b -d -c llbb -s XVH -m 2000 --fileName2 combine/test/mlfit_AZh_M1000.root
'''
