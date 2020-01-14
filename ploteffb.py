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

signal_list = [800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500,5000]


def draweffb(btag):
  upline = TGraph()
  downline = TGraph()
  if btag=='2b':
    eff_b = { 800 : [0.981, 1.019], 1000 : [0.975, 1.025], 1200 : [0.969, 1.031], 1400 : [0.957, 1.043], 1600 : [0.948, 1.052], 1800 : [0.944, 1.056], 2000 : [0.943, 1.056], 2500 : [0.945, 1.055], 3000 : [0.952, 1.048], 3500 : [0.961, 1.039], 4000 : [0.974, 1.026], 4500 : [0.988, 1.011], 5000 : [1.031, 0.969]}
  else:
    eff_b = {800 : [0.989, 1.011], 1000 : [0.991, 1.009], 1200 : [0.994, 1.006], 1400 : [0.995, 1.005], 1600 : [0.995, 1.005], 1800 : [0.997, 1.003], 2000 : [0.997, 1.003], 2500 : [1.004, 0.996], 3000 : [1.016, 0.984], 3500 : [1.026, 0.974], 4000 : [1.039, 0.961], 4500 : [1.054, 0.946], 5000 : [1.092, 0.907]}
  for i,signal in enumerate(signal_list):
    upline.SetPoint(i,signal,eff_b[signal][1])
    downline.SetPoint(i,signal,eff_b[signal][0])

  c1 = TCanvas("c1","sf plot", 800,600)
  upline.SetTitle("b-tagging uncertainty %s" %btag)
  upline.GetXaxis().SetTitle("m_{X} (GeV)")
  upline.GetYaxis().SetTitle("Uncertainty")
  upline.GetYaxis().SetRangeUser(0.7,1.3)
  upline.SetLineColor(4)
  downline.SetLineColor(2)
  upline.SetLineWidth(3)
  downline.SetLineWidth(3)
  leg = TLegend(0.7,0.65,0.85,0.8)
  leg.AddEntry(upline,"+ 1 s. d.","L")
  leg.AddEntry(downline,"- 1 s. d.","L")
  upline.Draw("AC")
  downline.Draw("C SAME")
  leg.Draw("SAME")
  c1.SaveAs("./plotsSF/effb_%s.pdf" % btag)
  c1.SaveAs("./plotsSF/effb_%s.png" % btag)
  c1.Close()

  
    
draweffb('0b')


