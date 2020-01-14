import os, re
from ROOT import TFile,TCanvas,gStyle, TLegend, TGraph
import numpy as np


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

path = '/work/pbaertsc/heavy_resonance/NanoTreeProducer/CorrectionTools/DYCorrection/'


def drawQCDNLOCorr():
  pt_list = []
  sf_list_z = []
  sf_list_w = []
  z_a = 1.423
  z_b = 2.257
  z_c = 0.451
  w_a = 1.024
  w_b = 3.072
  w_c = 0.749
  for pt in np.arange(0,1600,10):
    pt_list.append(float(pt))
    sf_z = z_a*np.exp(-z_b*pt/1000)+z_c
    sf_w = w_a*np.exp(-w_b*pt/1000)+w_c
    sf_list_z.append(sf_z)
    sf_list_w.append(sf_w)
  n = len(pt_list)
  pt_list = np.array(pt_list)
  sf_list_z = np.array(sf_list_z)
  sf_list_w = np.array(sf_list_w)
  line1 = TGraph(n,pt_list,sf_list_z)
  line2 = TGraph(n,pt_list,sf_list_w)
  c1 = TCanvas("c1","sf plot", 800,600)
  line1.SetTitle("QCD NLO SF")
  line1.GetXaxis().SetTitle("V p_{T} (GeV)")
  line1.GetYaxis().SetTitle("QCD NLO SF")
  line1.SetLineColor(4)
  line2.SetLineColor(2)
  line1.SetLineWidth(3)
  line2.SetLineWidth(3)
  leg = TLegend(0.7,0.65,0.85,0.8)
  leg.AddEntry(line1,"Z+Jets","L")
  leg.AddEntry(line2,"W+Jets","L")
  line1.Draw("AC")
  line2.Draw("C SAME")
  leg.Draw("SAME")
  c1.SaveAs("./plotsSF/QCDNLOCorr.pdf")
  c1.SaveAs("./plotsSF/QCDNLOCorr.png")
  c1.Close()


def drawQCDNNLOCorr():
    filename = path+"lindert_qcd_nnlo_sf.root"
    histname1 = "eej"
    histname2 = "evj"
    histname3 = "vvj"
    file1     = ensureTFile(filename)
    hist1     = file1.Get(histname1)
    hist2     = file1.Get(histname2)
    hist3     = file1.Get(histname3)
    hist1.SetDirectory(0)
    hist2.SetDirectory(0)
    hist3.SetDirectory(0)
    file1.Close()
    c1 = TCanvas("c1","sf plot", 800,600)
 
    hist1.SetTitle("QCD NNLO SF")
    hist1.GetXaxis().SetTitle("V p_{T} (GeV)")
    hist1.GetYaxis().SetTitle("QCD NNLO SF")
    hist1.GetYaxis().SetRangeUser(1.,1.18)
    hist1.SetLineColor(3)
    hist2.SetLineColor(2)
    hist3.SetLineColor(4)
    hist1.SetLineWidth(3)
    hist2.SetLineWidth(3)
    hist3.SetLineWidth(3)
    hist1.SetMarkerStyle(8)
    hist1.SetMarkerColor(3)
    hist2.SetMarkerStyle(8)
    hist2.SetMarkerColor(2)
    hist3.SetMarkerStyle(8)
    hist3.SetMarkerColor(4)
    leg = TLegend(0.15,0.6,0.3,0.8)
    leg.AddEntry(hist1,"DY+Jets")
    leg.AddEntry(hist2,"W+Jets")
    leg.AddEntry(hist3,"Z+Jets")
    c1.SetLogx()
    hist1.Draw("LP HIST")
    hist2.Draw("LP HIST SAME")
    hist3.Draw("LP HIST SAME")
    leg.Draw("SAME")
    gStyle.SetOptStat(0)
    c1.SaveAs("./plotsSF/QCDNNLOCorr.pdf")
    c1.SaveAs("./plotsSF/QCDNNLOCorr.png")
    c1.Close()



def drawEWKNLOCorr():
    filename1 = path+"merged_kfactors_zjets.root"
    filename2 = path+"merged_kfactors_wjets.root"
    histname = "kfactor_monojet_ewk"
    file1     = ensureTFile(filename1)
    file2     = ensureTFile(filename2)
    hist1     = file1.Get(histname)
    hist2     = file2.Get(histname)
    hist1.SetDirectory(0)
    hist2.SetDirectory(0)
    file1.Close()
    file2.Close()
    c1 = TCanvas("c1","sf plot", 800,600)
 
    hist1.SetTitle("EWK NLO SF")
    hist1.GetXaxis().SetTitle("V p_{T} (GeV)")
    hist1.GetYaxis().SetTitle("EWK NLO SF")
    hist1.GetYaxis().SetRangeUser(0.7,1.0)
    hist1.SetLineColor(4)
    hist2.SetLineColor(2)
    hist1.SetLineWidth(3)
    hist1.SetMarkerStyle(8)
    hist1.SetMarkerColor(4)
    hist2.SetMarkerStyle(8)
    hist2.SetMarkerColor(2)
    leg = TLegend(0.7,0.65,0.85,0.8)
    leg.AddEntry(hist1,"Z+Jets")
    leg.AddEntry(hist2,"W+Jets")
    hist1.Draw("LP HIST")
    hist2.Draw("LP HIST SAME")
    leg.Draw("SAME")
    gStyle.SetOptStat(0)
    c1.SaveAs("./plotsSF/EWKNLOCorr.pdf")
    c1.SaveAs("./plotsSF/EWKNLOCorr.png")
    c1.Close()
  
    
drawQCDNLOCorr()    
drawQCDNNLOCorr()
drawEWKNLOCorr()


