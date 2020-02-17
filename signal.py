#/usr/bin/env python

import os, sys, getopt, multiprocessing
import copy, math, pickle
from array import array
from ROOT import gROOT, gSystem, gStyle, gRandom
from ROOT import TMath, TFile, TChain, TTree, TCut, TH1F, TH2F, THStack, TGraph, TGaxis
from ROOT import TStyle, TCanvas, TPad, TLegend, TLatex, TText, TColor
from ROOT import TH1, TF1, TGraph, TGraphErrors, TGraphAsymmErrors, TVirtualFitter

gSystem.Load("PDFs/HWWLVJRooPdfs_cxx.so")
from ROOT import RooFit, RooRealVar, RooDataHist, RooDataSet, RooAbsData, RooAbsReal, RooAbsPdf, RooPlot, RooBinning, RooCategory, RooSimultaneous, RooArgList, RooArgSet, RooWorkspace, RooMsgService
from ROOT import RooFormulaVar, RooGenericPdf, RooGaussian, RooExponential, RooPolynomial, RooChebychev, RooBreitWigner, RooCBShape, RooDoubleCrystalBall, RooExtendPdf, RooAddPdf

from alpha import drawPlot
from rooUtils import *
from samples import sample,sample_2016,sample_2017,sample_2018
from selections import selection
from xsections import xsection

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-a", "--all", action="store_true", default=False, dest="all")
parser.add_option("-b", "--bash", action="store_true", default=False, dest="bash")
parser.add_option("-c", "--channel", action="store", type="string", dest="channel", default="")
parser.add_option("-s", "--signal", action="store", type="string", dest="signal", default="XZH")
parser.add_option("-e", "--efficiency", action="store_true", default=False, dest="efficiency")
parser.add_option("-p", "--parallelize", action="store_true", default=False, dest="parallelize")
parser.add_option("-v", "--verbose", action="store_true", default=False, dest="verbose")
(options, args) = parser.parse_args()
if options.bash: gROOT.SetBatch(True)

colour = [
    TColor(1001, 0., 0., 0., "black", 1.),
    TColor(1002, 230./255, 159./255, 0., "orange", 1.),
    TColor(1003, 86./255, 180./255, 233./255, "skyblue", 1.),
    TColor(1004, 0., 158./255, 115./255, "bluishgreen", 1.),
    TColor(1005, 0., 114./255, 178./255, "blue", 1.),
    TColor(1006, 213./255, 94./255, 0., "vermillion", 1.),
    TColor(1007, 204./255, 121./255, 167./255, "reddishpurple", 1.),
]

########## SETTINGS ##########

# Silent RooFit
RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)

#gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.05)
gStyle.SetErrorX(0.)

NTUPLEDIR   = "/work/pbaertsc/heavy_resonance/"
PLOTDIR     = "/work/pbaertsc/heavy_resonance/ZprimeToZHAnalysis/plotsSignal/"
CARDDIR     = "datacards/"
WORKDIR     = "workspace/"
RATIO       = 4
LUMI        = 137190.
YEAR        = 'combined'
VERBOSE     = options.verbose
PARALLELIZE = True
READTREE    = True
SIGNAL      = options.signal


channelList = ['nnbb', 'mmbb', 'eebb','nn0b','ee0b','mm0b','nnbbVBF', 'mmbbVBF', 'eebbVBF','nn0bVBF','ee0bVBF','mm0bVBF']
signalList = ['XZH','XZHVBF']

color = {'nnbb' : 634, 'eebb' : 418, 'mmbb' : 602,'nn0b' : 634, 'ee0b' : 418, 'mm0b' : 602,'nnbbVBF' : 634, 'eebbVBF' : 418, 'mmbbVBF' : 602,'nn0bVBF' : 634,'ee0bVBF' : 418,'mm0bVBF' : 602}

jobs = []

LOWMIN = 30.
LOWMAX = 65. 
LOWINT = 85.
SIGMIN = 105.
SIGMAX = 135.

HIGMIN = 135.
HIGMAX = 250.


XBINMIN=750.
XBINMAX= 6750.
XBINS  = 120

XTBINMIN = 1000.
XTBINMAX= 6750.
XTBINS  = 115

def signal(channel, stype):
    if 'VBF' in channel:
        stype = 'XZHVBF'
    else:
        stype = 'XZH'
    # HVT model
    if stype.startswith('X'):
        signalType = 'HVT'
        genPoints = [800, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
        massPoints = [x for x in range(800, 5000+1, 100)]
        interPar = True
    else:
        print "Signal type", stype, "not recognized"
        return
    
    n = len(genPoints)  
    
    category = channel
    cColor = color[category] if category in color else 1

    nElec = channel.count('e')
    nMuon = channel.count('m')
    nLept = nElec + nMuon
    nBtag = channel.count('b')
    if '0b' in channel:
        nBtag = 0

    X_name = "VH_mass"

    if not os.path.exists(PLOTDIR+stype+category): os.makedirs(PLOTDIR+stype+category)

    #*******************************************************#
    #                                                       #
    #              Variables and selections                 #
    #                                                       #
    #*******************************************************#
    X_mass = RooRealVar(  "X_mass",    "m_{ZH}",       XBINMIN, XBINMAX, "GeV")
    J_mass = RooRealVar(  "H_mass",   "jet mass",        LOWMIN, HIGMAX, "GeV")
    V_mass = RooRealVar(  "V_mass", "V jet mass",           -9.,  1.e6, "GeV")
    CSV1    = RooRealVar( "H_csv1",           "",         -999.,     2.     )
    CSV2    = RooRealVar( "H_csv2",           "",         -999.,     2.     )
    DeepCSV1= RooRealVar( "H_deepcsv1",       "",         -999.,     2.     )
    DeepCSV2= RooRealVar( "H_deepcsv2",       "",         -999.,     2.     )
    H_ntag  = RooRealVar( "H_ntag",           "",           -9.,     9.     )
    H_dbt   = RooRealVar( "H_dbt",            "",           -2.,     2.     )
    H_tau21 = RooRealVar( "H_tau21",          "",           -9.,     2.     )
    H_eta = RooRealVar( "H_eta",              "",           -9.,     9.     )
    H_tau21_ddt = RooRealVar( "H_ddt",  "",           -9.,     2.     )
    MaxBTag = RooRealVar( "MaxBTag",          "",          -10.,     2.     )
    H_chf   = RooRealVar( "H_chf",            "",           -1.,     2.     )
    MinDPhi = RooRealVar( "MinDPhi",          "",           -1.,    99.     )
    DPhi    = RooRealVar( "DPhi",             "",           -1.,    99.     )
    DEta    = RooRealVar( "DEta",             "",           -1.,    99.     )
    Mu1_relIso = RooRealVar( "Mu1_relIso",    "",           -1.,    99.     )
    Mu2_relIso = RooRealVar( "Mu2_relIso",    "",           -1.,    99.     )
    nTaus   = RooRealVar( "nTaus",            "",           -1.,    99.     )
    Vpt     = RooRealVar( "V.Pt()",           "",           -1.,   1.e6     )
    V_pt     = RooRealVar( "V_pt",            "",           -1.,   1.e6     )
    H_pt     = RooRealVar( "H_pt",            "",           -1.,   1.e6     )
    VH_deltaR=RooRealVar( "VH_deltaR",        "",           -1.,    99.     )
    isZtoNN = RooRealVar( "isZtoNN",          "",            0.,     2.     )
    isZtoEE = RooRealVar( "isZtoEE",          "",            0.,     2.     )
    isZtoMM = RooRealVar( "isZtoMM",          "",            0.,     2.     )
    isHtobb = RooRealVar( "isHtobb",          "",            0.,     2.     )
    isVBF   = RooRealVar( "isVBF",            "",            0.,     2.     )
    isMaxBTag_loose = RooRealVar( "isMaxBTag_loose", "",     0.,     2.     )
    weight  = RooRealVar( "eventWeightLumi",  "",         -1.e9,   1.e9     )

    Xmin = XBINMIN
    Xmax = XBINMAX

    # Define the RooArgSet which will include all the variables defined before
    # there is a maximum of 9 variables in the declaration, so the others need to be added with 'add'
    variables = RooArgSet(X_mass, J_mass, V_mass, CSV1, CSV2, H_ntag, H_dbt, H_tau21)
    variables.add(RooArgSet(DEta, DPhi, MaxBTag, MinDPhi, nTaus, Vpt))
    variables.add(RooArgSet(DeepCSV1, DeepCSV2,VH_deltaR, H_tau21_ddt))
    variables.add(RooArgSet(isZtoNN, isZtoEE, isZtoMM, isHtobb, isMaxBTag_loose, weight))
    variables.add(RooArgSet(isVBF, Mu1_relIso, Mu2_relIso, H_chf, H_pt, V_pt,H_eta))
    #X_mass.setRange("X_extended_range", X_mass.getMin(), X_mass.getMax())
    X_mass.setRange("X_reasonable_range", X_mass.getMin(), X_mass.getMax())
    X_mass.setRange("X_integration_range", Xmin, Xmax)
    X_mass.setBins(int((X_mass.getMax() - X_mass.getMin())/100))
    binsXmass = RooBinning(int((X_mass.getMax() - X_mass.getMin())/100), X_mass.getMin(), X_mass.getMax())
    X_mass.setBinning(binsXmass, "PLOT")
    massArg = RooArgSet(X_mass)

    # Cuts
    SRcut = selection[category]+selection['SR']
    print "  Cut:\t", SRcut
    #*******************************************************#
    #                                                       #
    #                    Signal fits                        #
    #                                                       #
    #*******************************************************#

    treeSign = {}
    setSignal = {}

    vmean  = {}
    vsigma = {}
    valpha1 = {}
    vslope1 = {}
    smean  = {}
    ssigma = {}
    salpha1 = {}
    sslope1 = {}
    salpha2 = {}
    sslope2 = {}
    a1 = {}
    a2 = {}
    sbrwig = {}
    signal = {}
    signalExt = {}
    signalYield = {}
    signalIntegral = {}
    signalNorm = {}
    signalXS = {}
    frSignal = {}
    frSignal1 = {}
    frSignal2 = {}
    frSignal3 = {}

    # Signal shape uncertainties (common amongst all mass points)
    xmean_fit = RooRealVar("sig_p1_fit", "Variation of the resonance position with the fit uncertainty", 0.005, -1., 1.)
    smean_fit = RooRealVar("CMSRunII_sig_p1_fit", "Change of the resonance position with the fit uncertainty", 0., -10, 10)
    xmean_jes = RooRealVar("sig_p1_scale_jes", "Variation of the resonance position with the jet energy scale", 0.010, -1., 1.) #0.001
    smean_jes = RooRealVar("CMSRunII_sig_p1_jes", "Change of the resonance position with the jet energy scale", 0., -10, 10)
    xmean_e = RooRealVar("sig_p1_scale_e", "Variation of the resonance position with the electron energy scale", 0.001, -1., 1.)
    smean_e = RooRealVar("CMSRunII_sig_p1_scale_e", "Change of the resonance position with the electron energy scale", 0., -10, 10)
    xmean_m = RooRealVar("sig_p1_scale_m", "Variation of the resonance position with the muon energy scale", 0.001, -1., 1.)
    smean_m = RooRealVar("CMSRunII_sig_p1_scale_m", "Change of the resonance position with the muon energy scale", 0., -10, 10)

    xsigma_fit = RooRealVar("sig_p2_fit", "Variation of the resonance width with the fit uncertainty", 0.02, -1., 1.)
    ssigma_fit = RooRealVar("CMSRunII_sig_p2_fit", "Change of the resonance width with the fit uncertainty", 0., -10, 10)
    xsigma_jes = RooRealVar("sig_p2_scale_jes", "Variation of the resonance width with the jet energy scale", 0.010, -1., 1.) #0.001
    ssigma_jes = RooRealVar("CMSRunII_sig_p2_jes", "Change of the resonance width with the jet energy scale", 0., -10, 10)
    xsigma_jer = RooRealVar("sig_p2_scale_jer", "Variation of the resonance width with the jet energy resolution", 0.020, -1., 1.)
    ssigma_jer = RooRealVar("CMSRunII_sig_p2_jer", "Change of the resonance width with the jet energy resolution", 0., -10, 10)
    xsigma_e = RooRealVar("sig_p2_scale_e", "Variation of the resonance width with the electron energy scale", 0.001, -1., 1.)
    ssigma_e = RooRealVar("CMSRunII_sig_p2_scale_e", "Change of the resonance width with the electron energy scale", 0., -10, 10)
    xsigma_m = RooRealVar("sig_p2_scale_m", "Variation of the resonance width with the muon energy scale", 0.040, -1., 1.)
    ssigma_m = RooRealVar("CMSRunII_sig_p2_scale_m", "Change of the resonance width with the muon energy scale", 0., -10, 10)
    
    xalpha1_fit = RooRealVar("sig_p3_fit", "Variation of the resonance alpha with the fit uncertainty", 0.03, -1., 1.)
    salpha1_fit = RooRealVar("CMSRunII_sig_p3_fit", "Change of the resonance alpha with the fit uncertainty", 0., -10, 10)
    
    xslope1_fit = RooRealVar("sig_p4_fit", "Variation of the resonance slope with the fit uncertainty", 0.10, -1., 1.)
    sslope1_fit = RooRealVar("CMSRunII_sig_p4_fit", "Change of the resonance slope with the fit uncertainty", 0., -10, 10)

    xmean_fit.setConstant(True)
    smean_fit.setConstant(True)
    xmean_jes.setConstant(True)
    smean_jes.setConstant(True)
    xmean_e.setConstant(True)
    smean_e.setConstant(True)
    xmean_m.setConstant(True)
    smean_m.setConstant(True)
    
    xsigma_fit.setConstant(True)
    ssigma_fit.setConstant(True)
    xsigma_jes.setConstant(True)
    ssigma_jes.setConstant(True)
    xsigma_jer.setConstant(True)
    ssigma_jer.setConstant(True)
    xsigma_e.setConstant(True)
    ssigma_e.setConstant(True)
    xsigma_m.setConstant(True)
    ssigma_m.setConstant(True)
    
    xalpha1_fit.setConstant(True)
    salpha1_fit.setConstant(True)
    xslope1_fit.setConstant(True)
    sslope1_fit.setConstant(True)

    # the alpha method is now done.
    for m in massPoints:
        signalString = "M%d" % m
        signalMass = "%s_M%d" % (stype, m)
        signalName = "%s%s_M%d" % (stype, category, m)
        signalColor = sample[signalMass]['linecolor'] if signalName in sample else 1

        # define the signal PDF
        vmean[m] = RooRealVar(signalName + "_vmean", "Crystal Ball mean", m, m*0.5, m*1.25)
        smean[m] = RooFormulaVar(signalName + "_mean", "@0*(1+@1*@2)*(1+@3*@4)*(1+@5*@6)*(1+@7*@8)", RooArgList(vmean[m], xmean_e, smean_e, xmean_m, smean_m, xmean_jes, smean_jes, xmean_fit, smean_fit))

        vsigma[m] = RooRealVar(signalName + "_vsigma", "Crystal Ball sigma", m*0.035, m*0.01, m*0.4)
        sigmaList = RooArgList(vsigma[m], xsigma_e, ssigma_e, xsigma_m, ssigma_m, xsigma_jes, ssigma_jes, xsigma_jer, ssigma_jer)
        sigmaList.add(RooArgList(xsigma_fit, ssigma_fit))
        ssigma[m] = RooFormulaVar(signalName + "_sigma", "@0*(1+@1*@2)*(1+@3*@4)*(1+@5*@6)*(1+@7*@8)*(1+@9*@10)", sigmaList)
        
        valpha1[m] = RooRealVar(signalName + "_valpha1", "Crystal Ball alpha", 1.,  0., 5.) # number of sigmas where the exp is attached to the gaussian core. >0 left, <0 right
        salpha1[m] = RooFormulaVar(signalName + "_alpha1", "@0*(1+@1*@2)", RooArgList(valpha1[m], xalpha1_fit, salpha1_fit))

        vslope1[m] = RooRealVar(signalName + "_vslope1", "Crystal Ball slope", 10., 1., 60.) # slope of the power tail   #10 1 60
        sslope1[m] = RooFormulaVar(signalName + "_slope1", "@0*(1+@1*@2)", RooArgList(vslope1[m], xslope1_fit, sslope1_fit))

        salpha2[m] = RooRealVar(signalName + "_alpha2", "Crystal Ball alpha", 2,  1., 5.) # number of sigmas where the exp is attached to the gaussian core. >0 left, <0 right
        sslope2[m] = RooRealVar(signalName + "_slope2", "Crystal Ball slope", 10, 1.e-1, 115.) # slope of the power tail
        #define polynomial
        #a1[m] = RooRealVar(signalName + "_a1", "par 1 for polynomial", m, 0.5*m, 2*m)
        a1[m] = RooRealVar(signalName + "_a1", "par 1 for polynomial", 0.001*m, 0.0005*m, 0.01*m)
        a2[m] = RooRealVar(signalName + "_a2", "par 2 for polynomial", 0.05, -1.,1.)
        #if channel=='nnbbVBF' or channel=='nn0bVBF':
        #    signal[m] = RooPolynomial(signalName,"m_{%s'} = %d GeV" % (stype[1], m) , X_mass, RooArgList(a1[m],a2[m]))
        #else:
        #    signal[m] = RooCBShape(signalName, "m_{%s'} = %d GeV" % (stype[1], m), X_mass, smean[m], ssigma[m], salpha1[m], sslope1[m]) # Signal name does not have the channel
        signal[m] = RooCBShape(signalName, "m_{%s'} = %d GeV" % (stype[1], m), X_mass, smean[m], ssigma[m], salpha1[m], sslope1[m]) # Signal name does not have the channel
        # extend the PDF with the yield to perform an extended likelihood fit
        signalYield[m] = RooRealVar(signalName+"_yield", "signalYield", 100, 0., 1.e6)
        signalNorm[m] = RooRealVar(signalName+"_norm", "signalNorm", 1., 0., 1.e6)
        signalXS[m] = RooRealVar(signalName+"_xs", "signalXS", 1., 0., 1.e6)
        signalExt[m] = RooExtendPdf(signalName+"_ext", "extended p.d.f", signal[m], signalYield[m])
        
        vslope1[m].setMax(50.)
        vslope1[m].setVal(20.)
        #valpha1[m].setVal(1.0)
        #valpha1[m].setConstant(True)
        
        if 'bb' in channel and 'VBF' not in channel:
            if 'nn' in channel:
                valpha1[m].setVal(0.5)
        elif '0b' in channel and 'VBF' not in channel:
            if 'nn' in channel:
                if m==800:
                    valpha1[m].setVal(2.)
                    vsigma[m].setVal(m*0.04)
            elif 'ee' in channel:
                valpha1[m].setVal(0.8)
                if m==800:
                    #valpha1[m].setVal(1.2)
                    valpha1[m].setVal(2.5)
                    vslope1[m].setVal(50.)
            elif 'mm' in channel:
                if m==800:
                    valpha1[m].setVal(2.)
                    vsigma[m].setVal(m*0.03)
                else:
                    vmean[m].setVal(m*0.9)
                    vsigma[m].setVal(m*0.08)
        elif 'bb' in channel and 'VBF' in channel:
            if 'nn' in channel:
                if m!=1800:
                    vmean[m].setVal(m*0.8)
                vsigma[m].setVal(m*0.08)
                valpha1[m].setMin(1.)
            elif 'ee' in channel:
                valpha1[m].setVal(0.7)
            elif 'mm' in channel:
                if m==800:
                    vslope1[m].setVal(50.)
                valpha1[m].setVal(0.7)
        elif '0b' in channel and 'VBF' in channel:
            if 'nn' in channel:
                valpha1[m].setVal(3.) 
                vmean[m].setVal(m*0.8)
                vsigma[m].setVal(m*0.08)
                valpha1[m].setMin(1.)
            elif 'ee' in channel:
                if m<2500:
                    valpha1[m].setVal(2.)
                if m==800:
                    vsigma[m].setVal(m*0.05)
                elif m==1000:
                    vsigma[m].setVal(m*0.03)
                elif m>1000 and m<1800:
                    vsigma[m].setVal(m*0.04)
            elif 'mm' in channel:
                if m<2000:
                    valpha1[m].setVal(2.)
                if m==1000 or m==1800:
                    vsigma[m].setVal(m*0.03)
                elif m==1200 or m==1600:
                    vsigma[m].setVal(m*0.04)

            
        #if m < 1000: vsigma[m].setVal(m*0.06)

        # If it's not the proper channel, make it a gaussian
        #if nLept==0 and 'VBF' in channel:
        #    valpha1[m].setVal(5)
        #    valpha1[m].setConstant(True)
        #    vslope1[m].setConstant(True)
        #    salpha2[m].setConstant(True)
        #    sslope2[m].setConstant(True)

        
        # ---------- if there is no simulated signal, skip this mass point ----------
        if m in genPoints:
            if VERBOSE: print " - Mass point", m

            # define the dataset for the signal applying the SR cuts
            treeSign[m] = TChain("tree")
            for j, ss in enumerate(sample[signalMass]['files']):
                treeSign[m].Add(NTUPLEDIR + ss + ".root")
            
            if treeSign[m].GetEntries() <= 0.:
                if VERBOSE: print " - 0 events available for mass", m, "skipping mass point..."
                signalNorm[m].setVal(-1)
                vmean[m].setConstant(True)
                vsigma[m].setConstant(True)
                salpha1[m].setConstant(True)
                sslope1[m].setConstant(True)
                salpha2[m].setConstant(True)
                sslope2[m].setConstant(True)
                signalNorm[m].setConstant(True)
                signalXS[m].setConstant(True)
                continue
            
            setSignal[m] = RooDataSet("setSignal_"+signalName, "setSignal", variables, RooFit.Cut(SRcut), RooFit.WeightVar(weight), RooFit.Import(treeSign[m]))
            if VERBOSE: print " - Dataset with", setSignal[m].sumEntries(), "events loaded"
            
            # FIT
            signalYield[m].setVal(setSignal[m].sumEntries())
            
            if treeSign[m].GetEntries(SRcut) > 5:
                if VERBOSE: print " - Running fit"
 
                frSignal[m] = signalExt[m].fitTo(setSignal[m], RooFit.Save(1), RooFit.Extended(True), RooFit.SumW2Error(True), RooFit.PrintLevel(-1))
                if VERBOSE: print "********** Fit result [", m, "] **", category, "*"*40, "\n", frSignal[m].Print(), "\n", "*"*80
                if VERBOSE: frSignal[m].correlationMatrix().Print()
                drawPlot(signalMass, stype+channel, X_mass, signal[m], setSignal[m], frSignal[m])
            
            else:
                print "  WARNING: signal", stype, "and mass point", m, "in channel", channel, "has 0 entries or does not exist"          
            # Remove HVT cross section (which is the same for Zlep and Zinv)
            if stype == "XZHVBF":
                sample_name = 'Zprime_VBF_Zh_Zlephinc_narrow_M-%d' % m
            else:
                sample_name = 'ZprimeToZHToZlepHinc_narrow_M%d' % m

            xs = xsection[sample_name]['xsec']
            
            signalXS[m].setVal(xs * 1000.)
            
            signalIntegral[m] = signalExt[m].createIntegral(massArg, RooFit.NormSet(massArg), RooFit.Range("X_integration_range"))
            boundaryFactor = signalIntegral[m].getVal()
            if VERBOSE: 
                print " - Fit normalization vs integral:", signalYield[m].getVal(), "/", boundaryFactor, "events"
            if channel=='nnbb' and m==5000:
                signalNorm[m].setVal(2.5)
            elif channel=='nn0b' and m==5000:
                signalNorm[m].setVal(6.7)
            else:
                signalNorm[m].setVal( boundaryFactor * signalYield[m].getVal() / signalXS[m].getVal()) # here normalize to sigma(X) x Br(X->VH) = 1 [fb]
            
            
        a1[m].setConstant(True)
        a2[m].setConstant(True)
        vmean[m].setConstant(True)
        vsigma[m].setConstant(True)
        valpha1[m].setConstant(True)
        vslope1[m].setConstant(True)
        salpha2[m].setConstant(True)
        sslope2[m].setConstant(True)
        signalNorm[m].setConstant(True)
        signalXS[m].setConstant(True)

    #*******************************************************#
    #                                                       #
    #                 Signal interpolation                  #
    #                                                       #
    #*******************************************************#


    # ====== CONTROL PLOT ======
    c_signal = TCanvas("c_signal", "c_signal", 800, 600)
    c_signal.cd()
    frame_signal = X_mass.frame()
    for m in genPoints[:-2]:
        if m in signalExt.keys():
            signal[m].plotOn(frame_signal, RooFit.LineColor(sample["%s_M%d" % (stype, m)]['linecolor']), RooFit.Normalization(signalNorm[m].getVal(), RooAbsReal.NumEvent), RooFit.Range("X_reasonable_range"))
    frame_signal.GetXaxis().SetRangeUser(0, 6500)
    frame_signal.Draw()
    drawCMS(-1, YEAR, "Simulation")
    drawAnalysis(channel)
    drawRegion(channel)
    c_signal.SaveAs(PLOTDIR+"/"+stype+category+"/"+stype+"_Signal.pdf")
    c_signal.SaveAs(PLOTDIR+"/"+stype+category+"/"+stype+"_Signal.png")
    #if VERBOSE: raw_input("Press Enter to continue...")
    # ====== CONTROL PLOT ======

    # Normalization
    gnorm = TGraphErrors()
    gnorm.SetTitle(";m_{X} (GeV);integral (GeV)")
    gnorm.SetMarkerStyle(20)
    gnorm.SetMarkerColor(1)
    gnorm.SetMaximum(0)
    inorm = TGraphErrors()
    inorm.SetMarkerStyle(24)
    fnorm = TF1("fnorm", "pol9", 800, 5000) #"pol5" if not channel=="XZHnnbb" else "pol6" #pol5*TMath::Floor(x-1800) + ([5]*x + [6]*x*x)*(1-TMath::Floor(x-1800))
    fnorm.SetLineColor(920)
    fnorm.SetLineStyle(7)
    fnorm.SetFillColor(2)
    fnorm.SetLineColor(cColor)

    # Mean
    gmean = TGraphErrors()
    gmean.SetTitle(";m_{X} (GeV);gaussian mean (GeV)")
    gmean.SetMarkerStyle(20)
    gmean.SetMarkerColor(cColor)
    gmean.SetLineColor(cColor)
    imean = TGraphErrors()
    imean.SetMarkerStyle(24)
    fmean = TF1("fmean", "pol1", 0, 5000)
    fmean.SetLineColor(2)
    fmean.SetFillColor(2)

    # Width
    gsigma = TGraphErrors()
    gsigma.SetTitle(";m_{X} (GeV);gaussian width (GeV)")
    gsigma.SetMarkerStyle(20)
    gsigma.SetMarkerColor(cColor)
    gsigma.SetLineColor(cColor)
    isigma = TGraphErrors()
    isigma.SetMarkerStyle(24)
    fsigma = TF1("fsigma", "pol1", 0, 5000)
    fsigma.SetLineColor(2)
    fsigma.SetFillColor(2)

    # Alpha1
    galpha1 = TGraphErrors()
    galpha1.SetTitle(";m_{X} (GeV);crystal ball lower alpha")
    galpha1.SetMarkerStyle(20)
    galpha1.SetMarkerColor(cColor)
    galpha1.SetLineColor(cColor)
    ialpha1 = TGraphErrors()
    ialpha1.SetMarkerStyle(24)
    falpha1 = TF1("falpha", "pol0", 0, 5000)
    falpha1.SetLineColor(2)
    falpha1.SetFillColor(2)

    # Slope1
    gslope1 = TGraphErrors()
    gslope1.SetTitle(";m_{X} (GeV);exponential lower slope (1/Gev)")
    gslope1.SetMarkerStyle(20)
    gslope1.SetMarkerColor(cColor)
    gslope1.SetLineColor(cColor)
    islope1 = TGraphErrors()
    islope1.SetMarkerStyle(24)
    fslope1 = TF1("fslope", "pol0", 0, 5000)
    fslope1.SetLineColor(2)
    fslope1.SetFillColor(2)

    # Alpha2
    galpha2 = TGraphErrors()
    galpha2.SetTitle(";m_{X} (GeV);crystal ball upper alpha")
    galpha2.SetMarkerStyle(20)
    galpha2.SetMarkerColor(cColor)
    galpha2.SetLineColor(cColor)
    ialpha2 = TGraphErrors()
    ialpha2.SetMarkerStyle(24)
    falpha2 = TF1("falpha", "pol0", 0, 5000)
    falpha2.SetLineColor(2)
    falpha2.SetFillColor(2)

    # Slope2
    gslope2 = TGraphErrors()
    gslope2.SetTitle(";m_{X} (GeV);exponential upper slope (1/Gev)")
    gslope2.SetMarkerStyle(20)
    gslope2.SetMarkerColor(cColor)
    gslope2.SetLineColor(cColor)
    islope2 = TGraphErrors()
    islope2.SetMarkerStyle(24)
    fslope2 = TF1("fslope", "pol0", 0, 5000)
    fslope2.SetLineColor(2)
    fslope2.SetFillColor(2)



    n = 0
    for i, m in enumerate(genPoints):
        if not m in signalNorm.keys(): continue
        if signalNorm[m].getVal() < 1.e-6: continue
        signalString = "M%d" % m
        signalName = "%s_M%d" % (stype, m)

        if gnorm.GetMaximum() < signalNorm[m].getVal(): gnorm.SetMaximum(signalNorm[m].getVal())
        gnorm.SetPoint(n, m, signalNorm[m].getVal())
        gmean.SetPoint(n, m, vmean[m].getVal())
        gmean.SetPointError(n, 0, min(vmean[m].getError(), vmean[m].getVal()*0.02))
        gsigma.SetPoint(n, m, vsigma[m].getVal())
        gsigma.SetPointError(n, 0, min(vsigma[m].getError(), vsigma[m].getVal()*0.05))
        galpha1.SetPoint(n, m, valpha1[m].getVal())
        galpha1.SetPointError(n, 0, min(valpha1[m].getError(), valpha1[m].getVal()*0.10))
        gslope1.SetPoint(n, m, vslope1[m].getVal())
        gslope1.SetPointError(n, 0, min(vslope1[m].getError(), vslope1[m].getVal()*0.10))
        galpha2.SetPoint(n, m, salpha2[m].getVal())
        galpha2.SetPointError(n, 0, min(salpha2[m].getError(), salpha2[m].getVal()*0.10))
        gslope2.SetPoint(n, m, sslope2[m].getVal())
        gslope2.SetPointError(n, 0, min(sslope2[m].getError(), sslope2[m].getVal()*0.10))
        n = n + 1
    print "fit on gmean:"
    gmean.Fit(fmean, "Q0", "SAME")
    print "fit on gsigma:"
    gsigma.Fit(fsigma, "Q0", "SAME")
    print "fit on galpha:"
    galpha1.Fit(falpha1, "Q0", "SAME")
    print "fit on gslope:"
    gslope1.Fit(fslope1, "Q0", "SAME")
    galpha2.Fit(falpha2, "Q0", "SAME")
    gslope2.Fit(fslope2, "Q0", "SAME")
    #for m in [5000, 5500]: gnorm.SetPoint(gnorm.GetN(), m, gnorm.Eval(m, 0, "S"))
    gnorm.Fit(fnorm, "Q", "SAME", 700, 5000)

    for m in massPoints:
        signalName = "%s_M%d" % (stype, m)
        
        if vsigma[m].getVal() < 10.: vsigma[m].setVal(10.)

        # Interpolation method
        syield = gnorm.Eval(m)
        spline = gnorm.Eval(m, 0, "S")
        sfunct = fnorm.Eval(m)
        
        #delta = min(abs(1.-spline/sfunct), abs(1.-spline/syield))
        delta = abs(1.-spline/sfunct) if sfunct > 0 else 0
        syield = spline
               
        if interPar:
            jmean = gmean.Eval(m)
            jsigma = gsigma.Eval(m)
            jalpha1 = galpha1.Eval(m)
            jslope1 = gslope1.Eval(m)
        else:
            jmean = fmean.GetParameter(0) + fmean.GetParameter(1)*m + fmean.GetParameter(2)*m*m
            jsigma = fsigma.GetParameter(0) + fsigma.GetParameter(1)*m + fsigma.GetParameter(2)*m*m
            jalpha1 = falpha1.GetParameter(0) + falpha1.GetParameter(1)*m + falpha1.GetParameter(2)*m*m
            jslope1 = fslope1.GetParameter(0) + fslope1.GetParameter(1)*m + fslope1.GetParameter(2)*m*m

        inorm.SetPoint(inorm.GetN(), m, syield)
        signalNorm[m].setVal(syield)

        imean.SetPoint(imean.GetN(), m, jmean)
        if jmean > 0: vmean[m].setVal(jmean)

        isigma.SetPoint(isigma.GetN(), m, jsigma)
        if jsigma > 0: vsigma[m].setVal(jsigma)

        ialpha1.SetPoint(ialpha1.GetN(), m, jalpha1)
        if not jalpha1==0: valpha1[m].setVal(jalpha1)

        islope1.SetPoint(islope1.GetN(), m, jslope1)
        if jslope1 > 0: vslope1[m].setVal(jslope1)
    

    c1 = TCanvas("c1", "Crystal Ball", 1200, 800)
    c1.Divide(2, 2)
    c1.cd(1)
    gmean.SetMinimum(0.)
    gmean.Draw("APL")
    imean.Draw("P, SAME")
    drawRegion(channel)
    c1.cd(2)
    gsigma.SetMinimum(0.)
    gsigma.Draw("APL")
    isigma.Draw("P, SAME")
    drawRegion(channel)
    c1.cd(3)
    galpha1.Draw("APL")
    ialpha1.Draw("P, SAME")
    drawRegion(channel)
    galpha1.GetYaxis().SetRangeUser(0., 5.)
    c1.cd(4)
    gslope1.Draw("APL")
    islope1.Draw("P, SAME")
    drawRegion(channel)
    gslope1.GetYaxis().SetRangeUser(0., 125.)
    if False:
        c1.cd(5)
        galpha2.Draw("APL")
        ialpha2.Draw("P, SAME")
        drawRegion(channel)
        c1.cd(6)
        gslope2.Draw("APL")
        islope2.Draw("P, SAME")
        drawRegion(channel)
        gslope2.GetYaxis().SetRangeUser(0., 10.)


    c1.Print(PLOTDIR+stype+category+"/"+stype+"_SignalShape.pdf")
    c1.Print(PLOTDIR+stype+category+"/"+stype+"_SignalShape.png")


    c2 = TCanvas("c2", "Signal Efficiency", 800, 600)
    c2.cd(1)
    gnorm.SetMarkerColor(cColor)
    gnorm.SetMarkerStyle(20)
    gnorm.SetLineColor(cColor)
    gnorm.SetLineWidth(2)
    gnorm.Draw("APL")
    inorm.Draw("P, SAME")
    gnorm.GetXaxis().SetRangeUser(genPoints[0]-100, genPoints[-1]+100)
    gnorm.GetYaxis().SetRangeUser(0., gnorm.GetMaximum()*1.25)
    drawCMS(-1,YEAR , "Simulation")
    drawAnalysis(channel)
    drawRegion(channel)
    c2.Print(PLOTDIR+stype+category+"/"+stype+"_SignalNorm.pdf")
    c2.Print(PLOTDIR+stype+category+"/"+stype+"_SignalNorm.png")





    #*******************************************************#
    #                                                       #
    #                   Generate workspace                  #
    #                                                       #
    #*******************************************************#

    # create workspace
    w = RooWorkspace("ZH_RunII", "workspace")
    for m in massPoints:
        getattr(w, "import")(signal[m], RooFit.Rename(signal[m].GetName()))
        getattr(w, "import")(signalNorm[m], RooFit.Rename(signalNorm[m].GetName()))
        getattr(w, "import")(signalXS[m], RooFit.Rename(signalXS[m].GetName()))
    w.writeToFile("%s%s.root" % (WORKDIR, stype+channel), True)
    print "Workspace", "%s%s.root" % (WORKDIR, stype+channel), "saved successfully"
    sys.exit()

def efficiency(stype, Zlep=True):
    genPoints = [800, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500]
    eff = {}
    
    channels = [x for x in channelList if len(x)<5]

    for channel in channels:
        treeSign = {}
        ngenSign = {}
        nevtSign = {}
        eff[channel] = TGraphErrors()

        for i, m in enumerate(genPoints):
            signName = "%s_M%d" % (stype, m) #"%s_M%d" % (channel[:3], m)
            ngenSign[m] = 0.
            nevtSign[m] = 0.
            for j, ss in enumerate(sample[signMass]['files']):
                if 'nn' in channel and not 'Zinv' in ss: continue
                if ('en' in channel or 'mn' in channel) and not 'Wlep' in ss: continue
                if ('ee' in channel or 'mm' in channel) and not 'Zlep' in ss: continue
                if Zlep and 'Zinv' in ss: continue
                if not Zlep and 'Zlep' in ss: continue

                sfile = TFile(NTUPLEDIR + ss + ".root", "READ")
                if not sfile.Get("Events")==None:
                    ngenSign[m] += sfile.Get("Events").GetEntries()
                    # From trees
                    treeSign[m] = sfile.Get("tree")
                    nevtSign[m] += treeSign[m].GetEntries(selection[channel] + selection['SR'])
                else:
                    ngenSign[m] = -1
                    print "Failed reading file", NTUPLEDIR + ss + ".root"
                sfile.Close()
            if nevtSign[m] == 0 or ngenSign[m] < 0: continue
            # Gen Br
            n = eff[channel].GetN()
            eff[channel].SetPoint(n, m, nevtSign[m]/ngenSign[m])
            eff[channel].SetPointError(n, 0, math.sqrt(nevtSign[m])/ngenSign[m])

        eff[channel].SetMarkerColor(color[channel])
        eff[channel].SetMarkerStyle(20)
        eff[channel].SetLineColor(color[channel])
        eff[channel].SetLineWidth(2)
        if channel.count('b')==1: eff[channel].SetLineStyle(3)

    n = max([eff[x].GetN() for x in channels])
    maxEff = 0.

    # Total efficiency
    eff["sum"] = TGraphErrors(n)
    eff["sum"].SetMarkerStyle(24)
    eff["sum"].SetMarkerColor(1)
    eff["sum"].SetLineWidth(2)
    for i in range(n):
        tot, mass = 0., 0.
        for channel in channels:
            if eff[channel].GetN() > i:
                tot += eff[channel].GetY()[i]
                mass = eff[channel].GetX()[i]
                if tot > maxEff: maxEff = tot
        eff["sum"].SetPoint(i, mass, tot)


    leg = TLegend(0.15, 0.60, 0.95, 0.8)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0) #1001
    leg.SetFillColor(0)
    leg.SetNColumns(len(channels)/4)
    for i, channel in enumerate(channels):
        if eff[channel].GetN() > 0: leg.AddEntry(eff[channel], getChannel(channel), "pl")
    leg.SetY1(leg.GetY2()-len([x for x in channels if eff[x].GetN() > 0])/2.*0.045)

    legS = TLegend(0.55, 0.85-0.045, 0.95, 0.85)
    legS.SetBorderSize(0)
    legS.SetFillStyle(0) #1001
    legS.SetFillColor(0)
    legS.AddEntry(eff['sum'], "Total efficiency", "pl")

    c1 = TCanvas("c1", "Signal Efficiency", 1200, 800)
    c1.cd(1)
    eff['sum'].Draw("APL")
    for i, channel in enumerate(channels): eff[channel].Draw("SAME, PL")
    leg.Draw()
    legS.Draw()
    setHistStyle(eff["sum"], 1.1)
    eff["sum"].SetTitle(";m_{"+stype[1]+"'} (GeV);Acceptance #times efficiency")
    eff["sum"].SetMinimum(0.)
    eff["sum"].SetMaximum(max(1., maxEff*1.5)) #0.65
    eff["sum"].GetXaxis().SetTitleSize(0.045)
    eff["sum"].GetYaxis().SetTitleSize(0.045)
    eff["sum"].GetYaxis().SetTitleOffset(1.1)
    eff["sum"].GetXaxis().SetTitleOffset(1.05)
    eff["sum"].GetXaxis().SetRangeUser(750, 5500)
    if stype=='XWH' or (stype=='XZH' and Zlep): line = drawLine(750, 2./3., 4500, 2./3.)
    drawCMS(-1,YEAR, "Simulation") #Preliminary
    drawAnalysis("ZH")

    suffix = ""
    if stype=='XZH' and Zlep: suffix = "ll"
    elif stype=='XZH' and not Zlep: suffix = "nn"
    elif stype=='XWH': suffix = "ln"

    c1.Print("plotsSignal/Efficiency/"+stype+suffix+".pdf")
    c1.Print("plotsSignal/Efficiency/"+stype+suffix+".png")

    # print
    print "category",
    for m in range(0, eff["sum"].GetN()):
        print " & %d" % int(eff["sum"].GetX()[m]),
    print "\\\\", "\n\\hline"
    for i, channel in enumerate(channels+["sum"]):
        if channel=='sum': print "\\hline"
        print getChannel(channel).replace("high ", "H").replace("low ", "L").replace("purity", "P").replace("b-tag", ""),
        for m in range(0, eff[channel].GetN()):
            print "& %.1f" % (100.*eff[channel].GetY()[m]),
        print "\\\\"



def efficiencyAll():
    #signals = {'XZHeebb':['eebb'],'XZHmmbb':['mmbb'],'XZHnnbb':['nnbb'],'XZHee0b':['ee0b'],'XZHmm0b':['mm0b'],'XZHnn0b':['nn0b'],'XZHVBFeebbVBF':['eebbVBF'],'XZHVBFmmbbVBF':['mmbbVBF'],'XZHVBFnnbbVBF':['nnbbVBF'],'XZHVBFee0bVBF':['ee0bVBF'],'XZHVBFmm0bVBF':['mm0bVBF'],'XZHVBFnn0bVBF':['nn0bVBF']}
    labels = {'XZHeebb' : "eeb#bar{b}",'XZHmmbb' : "#mu#mub#bar{b}",'XZHnnbb' : "#nu#nub#bar{b}",'XZHee0b' : "ee0b",'XZHmm0b' : "#mu#mu0b",'XZHnn0b' : "#nu#nu0b",'XZHVBFeebbVBF' : "eeb#bar{b}VBF",'XZHVBFmmbbVBF' : "#mu#mub#bar{b}VBF",'XZHVBFnnbbVBF' : "#nu#nub#bar{b}VBF",'XZHVBFee0bVBF' : "ee0bVBF",'XZHVBFmm0bVBF' : "#mu#mu0bVBF",'XZHVBFnn0bVBF' : "#nu#nu0bVBF"}
    colors = {'XZHeebb' : 2, 'XZHmmbb' : 4, 'XZHnnbb' : 2,'XZHee0b' : 3, 'XZHmm0b' : 6, 'XZHnn0b' : 4,'XZHVBFeebbVBF' : 2, 'XZHVBFmmbbVBF' : 4, 'XZHVBFnnbbVBF' : 2,'XZHVBFee0bVBF' : 3, 'XZHVBFmm0bVBF' : 6, 'XZHVBFnn0bVBF' : 4}
    styles = {'XZHeebb' : 1, 'XZHmmbb' : 1, 'XZHnnbb' : 1,'XZHee0b' : 1, 'XZHmm0b' : 1, 'XZHnn0b' : 1,'XZHVBFeebbVBF' : 1, 'XZHVBFmmbbVBF' : 1, 'XZHVBFnnbbVBF' : 1,'XZHVBFee0bVBF' : 1, 'XZHVBFmm0bVBF' : 1, 'XZHVBFnn0bVBF' : 1}
    marker = {'XZHeebb' : 22, 'XZHmmbb' : 20, 'XZHnnbb' : 22,'XZHee0b' : 22, 'XZHmm0b' : 20, 'XZHnn0b' : 20, 'XZHVBFeebbVBF' : 22, 'XZHVBFmmbbVBF' : 20, 'XZHVBFnnbbVBF' : 22,'XZHVBFee0bVBF' : 22, 'XZHVBFmm0bVBF' : 20, 'XZHVBFnn0bVBF' : 20}
    genPoints = [800, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
    eff = {}

    for signal_samples in ['ZlepHinc','ZinvHinc','ZinvHincVBF','ZlepHincVBF']:
        if signal_samples == 'ZinvHinc':
            signals = {'XZHnnbb':['nnbb'],'XZHnn0b':['nn0b']}
            sign_list = ['XZHnnbb','XZHnn0b']
        elif signal_samples == 'ZlepHinc':
            signals = {'XZHeebb':['eebb'],'XZHmmbb':['mmbb'],'XZHee0b':['ee0b'],'XZHmm0b':['mm0b']}
            sign_list = ['XZHeebb','XZHmmbb','XZHee0b', 'XZHmm0b']
        elif signal_samples == 'ZinvHincVBF':
            signals = {'XZHVBFnnbbVBF':['nnbbVBF'],'XZHVBFnn0bVBF':['nn0bVBF']}
            sign_list = ['XZHVBFnnbbVBF', 'XZHVBFnn0bVBF']
        elif signal_samples == 'ZlepHincVBF':
            signals = {'XZHVBFeebbVBF':['eebbVBF'],'XZHVBFmmbbVBF':['mmbbVBF'],'XZHVBFee0bVBF':['ee0bVBF'],'XZHVBFmm0bVBF':['mm0bVBF']}
            sign_list = ['XZHVBFeebbVBF', 'XZHVBFmmbbVBF','XZHVBFee0bVBF', 'XZHVBFmm0bVBF']
        for sign, channels in signals.iteritems():
      
            treeSign = {}
            ngenSign = {}
            nevtSign = {}
            eff[sign] = TGraphErrors()
            eff[sign].SetTitle(sign)
            eff[sign].SetMarkerColor(colors[sign])
            eff[sign].SetMarkerSize(1.25)
            eff[sign].SetLineColor(colors[sign])
            eff[sign].SetLineWidth(2)
            eff[sign].SetLineStyle(styles[sign])
            eff[sign].SetMarkerStyle(marker[sign])

            for i, m in enumerate(genPoints):
                neff = 0.
                for channel in channels:
                    if signal_samples == 'ZinvHinc':
                        file_list = ['Ntuples2016/XZH/ZprimeToZHToZinvHall_narrow_M%s'%m,'Ntuples2017/XZH/ZprimeToZHToZinvHall_narrow_M%s'%m,'Ntuples2018/XZH/ZprimeToZHToZinvHall_narrow_M%s'%m]
                    elif signal_samples == 'ZlepHinc':
                        file_list = ['Ntuples2016/XZH/ZprimeToZHToZlepHinc_narrow_M%s'%m,'Ntuples2017/XZH/ZprimeToZHToZlepHinc_narrow_M%s'%m,'Ntuples2018/XZH/ZprimeToZHToZlepHinc_narrow_M%s'%m]
                    elif signal_samples == 'ZinvHincVBF':
                        file_list = ['Ntuples2016/XZHVBF/Zprime_VBF_Zh_Zinvhinc_narrow_M-%s'%m,'Ntuples2017/XZHVBF/Zprime_VBF_Zh_Zinvhinc_narrow_M-%s'%m,'Ntuples2018/XZHVBF/Zprime_VBF_Zh_Zinvhinc_narrow_M-%s'%m]
                    elif signal_samples == 'ZlepHincVBF':
                        file_list = ['Ntuples2016/XZHVBF/Zprime_VBF_Zh_Zlephinc_narrow_M-%s'%m,'Ntuples2017/XZHVBF/Zprime_VBF_Zh_Zlephinc_narrow_M-%s'%m,'Ntuples2018/XZHVBF/Zprime_VBF_Zh_Zlephinc_narrow_M-%s'%m]
                    #if not 'VBF' in channel:
                    #    signMass = "XZH_M%d" % m
                    #else:
                    #    signMass = "XZHVBF_M%d" % m
                    ngenSign[m] = 0.
                    nevtSign[m] = 0.
                    #for j, ss in enumerate(sample[signMass]['files']):
                    for j, ss in enumerate(file_list):
                        sfile = TFile(NTUPLEDIR + ss + ".root", "READ")
                        if not sfile.Get("Events")==None:
                            ngenSign[m] += sfile.Get("Events").GetEntries()
                            # From trees
                            treeSign[m] = sfile.Get("tree")
                            nevtSign[m] += treeSign[m].GetEntries(selection[channel] + selection['SR'])
                        else:
                            ngenSign[m] = -1
                            print "Failed reading file", NTUPLEDIR + ss + ".root"
                        sfile.Close()
                    if nevtSign[m] == 0 or ngenSign[m] < 0: continue
                    # Gen Br
                    #print "m:",m
                    #print "nevtSign:",nevtSign[m]
                    #print "ngenSign:",ngenSign[m]
                    neff += nevtSign[m]/ngenSign[m]
                if 'ln' in sign or 'll' in sign: neff *= 1.5
                n = eff[sign].GetN()
                eff[sign].SetPoint(n, m, neff)
                eff[sign].SetPointError(n, 0, 0)

        n = 0. #max([eff[x].GetN() for x in channels])
        maxEff = 0.
        #sign_list = ['XZHeebb','XZHmmbb','XZHnnbb', 'XZHee0b', 'XZHmm0b', 'XZHnn0b','XZHVBFeebbVBF', 'XZHVBFmmbbVBF', 'XZHVBFnnbbVBF','XZHVBFee0bVBF', 'XZHVBFmm0bVBF', 'XZHVBFnn0bVBF']
        leg = TLegend(0.15, 0.15, 0.95, 0.35)
        #leg = TLegend(0.15, 0.7, 0.95, 0.8)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0) #1001
        leg.SetFillColor(0)
        for sign in sign_list:
            #if eff[sign].GetN() > 0:
            leg.AddEntry(eff[sign], labels[sign], "pl")
            n += 1
        leg.SetNColumns(int(n/3))
        leg.SetY1(leg.GetY2()-n*0.045/leg.GetNColumns())
    
        n_error = max([eff[x].GetN() for x in sign_list])
    
        # Total efficiency
        eff["sum"] = TGraphErrors(n_error)
        eff["sum"].SetMarkerStyle(24)
        eff["sum"].SetMarkerColor(1)
        eff["sum"].SetLineWidth(2)
        for i in range(n_error):
            tot, mass = 0., 0.
            for sign in sign_list:
                if eff[sign].GetN() > i:
                    tot += eff[sign].GetY()[i]
                    mass = eff[sign].GetX()[i]
            if tot > maxEff: maxEff = tot
            eff["sum"].SetPoint(i, mass, tot)

        #legS = TLegend(0.55, 0.85-0.045, 0.95, 0.85)
        legS = TLegend(0.55, 0.35-0.045, 0.95, 0.35)
        legS.SetBorderSize(0)
        legS.SetFillStyle(0) #1001
        legS.SetFillColor(0)
        legS.AddEntry(eff['sum'], "%s Total efficiency"%signal_samples, "pl")

        c1 = TCanvas("c1", "Signal Efficiency", 1200, 800)
        c1.cd(1)
        c1.GetPad(0).SetTicks(1, 1)
        c1.SetLogy()
        first = sign_list[0]
        #if eff['XZHeebb'].GetN()!=0:
        #    first = 'XZHeebb'
        #else:
        #    first = 'XZHnnbb'
        eff[first].Draw("APL")
        for sign, channels in signals.iteritems():
            eff[sign].Draw("APL" if i==0 else "SAME, PL")
        eff["sum"].Draw("SAME, PL")
        leg.Draw()
        legS.Draw()

        setHistStyle(eff[first], 1.1)
        eff[first].SetTitle(";m_{X} (GeV);Acceptance #times efficiency")
        eff[first].SetMinimum(0.)
        eff[first].SetMaximum(max(1., maxEff*1.5)) #0.65
        eff[first].GetXaxis().SetTitleSize(0.045)
        eff[first].GetYaxis().SetTitleSize(0.045)
        eff[first].GetXaxis().SetLabelSize(0.045)
        eff[first].GetYaxis().SetLabelSize(0.045)
        eff[first].GetYaxis().SetTitleOffset(1.1)
        eff[first].GetXaxis().SetTitleOffset(1.05)
        eff[first].GetXaxis().SetRangeUser(750,5500)
        eff[first].GetYaxis().SetRangeUser(0., 0.4)
        drawCMS(-1,YEAR, "Simulation")
        """
        latex = TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.05)
        latex.SetTextColor(1)
        latex.SetTextFont(42)
        latex.SetTextAlign(13)
        latex.DrawLatex(0.83, 0.99, "(13 TeV)")
        latex.SetTextFont(62)
        latex.SetTextSize(0.06)
        latex.DrawLatex(0.15, 0.90, "CMS")
        latex.SetTextSize(0.05)
        latex.SetTextFont(52)
        """
        #c1.Print("plotsSignal/Efficiency/Efficiency.pdf")
        #c1.Print("plotsSignal/Efficiency/Efficiency.png")
        c1.Print("plotsSignal/Efficiency/%s_Efficiency.pdf"%signal_samples)
        c1.Print("plotsSignal/Efficiency/%s_Efficiency.png"%signal_samples)


def plot(channel):
    # ====== CONTROL PLOT ======
    c_signal = TCanvas("c_signal", "c_signal", 800, 600)
    c_signal.cd()
    frame_signal = X_mass.frame()
    for m in genPoints[1:-1]: msignal[m].plotOn(frame_signal, RooFit.LineColor(sample["%s_M%d" % (channel[:3], m)]['linecolor']), RooFit.Normalization(1., RooAbsReal.NumEvent), RooFit.Name("M%d" % m))
    frame_signal.GetXaxis().SetRangeUser(XBINMIN, XBINMAX)
    frame_signal.GetYaxis().SetRangeUser(0, 1)
    frame_signal.GetYaxis().SetTitle("Arbitray units")
    frame_signal.Draw()
    drawCMS(-1,YEAR, "Simulation")
    drawAnalysis("ZH")
    drawRegion(channel)
    mleg = TLegend(0.75, 0.6, 0.98, 0.92)
    mleg.SetBorderSize(0)
    mleg.SetFillStyle(0) #1001
    mleg.SetFillColor(0)
    for m in genPoints[1:-1]: mleg.AddEntry("M%d" % m, "m_{X} = %d GeV" % m, "L")
    mleg.Draw()
    c_signal.SaveAs(PLOTDIR+"/"+channel+"/Signal.pdf")
    c_signal.SaveAs(PLOTDIR+"/"+channel+"/Signal.png")
    if VERBOSE: raw_input("Press Enter to continue...")


if __name__ == "__main__":
    if options.efficiency:
        efficiencyAll()

    elif options.all:
        gROOT.SetBatch(True)
        for c in channelList:
            if 'VBF' in c:
                s = 'XZHVBF'
            else:
                s = 'XZH'
                if PARALLELIZE:
                    p = multiprocessing.Process(target=signal, args=(c, s,))
                    jobs.append(p)
                    p.start()
                else:
                    signal(c, s)
    else:
        signal(options.channel, options.signal)



# python parametrization.py | grep @ > parameters_.py
