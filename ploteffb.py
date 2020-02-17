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
    eff_b = {800 : [0.961, 1.040], 1000 : [0.955, 1.046], 1200 : [0.944, 1.058], 1400 : [0.913, 1.093], 1600 : [0.889, 1.119], 1800 : [0.877, 1.133], 2000 : [0.873, 1.138], 2500 : [0.869, 1.141], 3000 : [0.868, 1.142], 3500 : [0.868, 1.143], 4000 : [0.868, 1.143], 4500 : [0.866, 1.145], 5000 : [0.862, 1.150]}
  else:
    eff_b = {800 : [1.007, 0.993], 1000 : [1.009, 0.991], 1200 : [1.013, 0.988], 1400 : [1.025, 0.976], 1600 : [1.033, 0.968], 1800 : [1.038, 0.963], 2000 : [1.041, 0.961], 2500 : [1.050, 0.952], 3000 : [1.054, 0.948], 3500 : [1.057, 0.945], 4000 : [1.059, 0.944], 4500 : [1.061, 0.941], 5000 : [1.063, 0.940]}
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
  leg = TLegend(0.7,0.75,0.85,0.85)
  leg.AddEntry(upline,"+ 1 s. d.","L")
  leg.AddEntry(downline,"- 1 s. d.","L")
  upline.Draw("AC")
  downline.Draw("C SAME")
  leg.Draw("SAME")
  c1.SaveAs("./plotsSF/effb_%s.pdf" % btag)
  c1.SaveAs("./plotsSF/effb_%s.png" % btag)
  c1.Close()

def drawjec(cut):
  scale_upline = TGraph()
  scale_downline = TGraph()
  res_upline = TGraph()
  res_downline = TGraph()
  if cut=='nnbb':
    eff_scale_vv = [0.993,1.012]
    eff_res_vv = [0.996,0.995]
    eff_scale = {800 : [0.991, 1.009], 1000 : [0.992, 1.008], 1200 : [0.993, 1.007], 1400 : [0.993, 1.007], 1600 : [0.993, 1.007], 1800 : [0.993, 1.007], 2000 : [0.993, 1.007], 2500 : [0.993, 1.007], 3000 : [0.992, 1.008], 3500 : [0.992, 1.008], 4000 : [0.992, 1.009], 4500 : [0.991, 1.009], 5000 : [0.991, 1.009]}
    eff_res = {800 : [1.007, 1.015], 1000 : [1.005, 1.010], 1200 : [1.004, 1.008], 1400 : [1.003, 1.005], 1600 : [1.002, 1.005], 1800 : [1.002, 1.005], 2000 : [1.001, 1.005], 2500 : [1.001, 1.008], 3000 : [1.001, 1.008], 3500 : [1.001, 1.008], 4000 : [1.001, 1.011], 4500 : [1.001, 1.009], 5000 : [1.001, 1.013]}
  elif cut=='eebb':
    eff_scale_vv = [0.997,1.003]
    eff_res_vv = [0.999,0.998]
  elif cut =='nnbbVBF':
    eff_scale_vv = [0.993,1.006]
    eff_res_vv = [1.000,1.005]
  elif cut == 'eebbVBF':
    eff_scale_vv = [0.982,0.999]
    eff_res_vv = [1.000,0.999]
  
  for i,signal in enumerate(signal_list):
    scale_upline.SetPoint(i,signal,eff_scale[signal][1])
    scale_downline.SetPoint(i,signal,eff_scale[signal][0])
  c1 = TCanvas("c1","sf plot", 800,600)
  scale_upline.SetTitle("Jet energy scale %s" %cut)
  scale_upline.GetXaxis().SetTitle("m_{X} (GeV)")
  scale_upline.GetYaxis().SetTitle("Uncertainty")
  scale_upline.GetYaxis().SetRangeUser(0.8,1.2)
  scale_upline.SetLineColor(4)
  scale_downline.SetLineColor(2)
  scale_upline.SetLineWidth(3)
  scale_downline.SetLineWidth(3)
  leg = TLegend(0.7,0.75,0.85,0.85)
  leg.AddEntry(scale_upline,"+ 1 s. d.","L")
  leg.AddEntry(scale_downline,"- 1 s. d.","L")
  scale_upline.Draw("AC")
  scale_downline.Draw("C SAME")
  leg.Draw("SAME")
  c1.SaveAs("./plotsSF/jes_%s.pdf" % cut)
  c1.SaveAs("./plotsSF/jes_%s.png" % cut)
  c1.Close() 
  for i,signal in enumerate(signal_list):
    res_upline.SetPoint(i,signal,eff_res[signal][1])
    res_downline.SetPoint(i,signal,eff_res[signal][0])
  c2 = TCanvas("c1","sf plot", 800,600)
  res_upline.SetTitle("Jet energy resolution %s" %cut)
  res_upline.GetXaxis().SetTitle("m_{X} (GeV)")
  res_upline.GetYaxis().SetTitle("Uncertainty")
  res_upline.GetYaxis().SetRangeUser(0.8,1.2)
  res_upline.SetLineColor(4)
  res_downline.SetLineColor(2)
  res_upline.SetLineWidth(3)
  res_downline.SetLineWidth(3)
  leg = TLegend(0.7,0.75,0.85,0.85)
  leg.AddEntry(res_upline,"+ 1 s. d.","L")
  leg.AddEntry(res_downline,"- 1 s. d.","L")
  res_upline.Draw("AC")
  res_downline.Draw("C SAME")
  leg.Draw("SAME")
  c2.SaveAs("./plotsSF/jer_%s.pdf" % cut)
  c2.SaveAs("./plotsSF/jer_%s.png" % cut)
  c2.Close()   

def drawjmc(cut):
  scale_upline = TGraph()
  scale_downline = TGraph()
  res_upline = TGraph()
  res_downline = TGraph()
  if cut=='nnbb':
    eff_scale_vv = [0.996,1.006]
    eff_res_vv = [0.958,1.049]
    eff_scale = {800 : [0.994, 1.005], 1000 : [0.994, 1.005], 1200 : [0.995, 1.006], 1400 : [0.994, 1.006], 1600 : [0.994, 1.006], 1800 : [0.994, 1.006], 2000 : [0.994, 1.006], 2500 : [0.994, 1.006], 3000 : [0.994, 1.006], 3500 : [0.994, 1.006], 4000 : [0.994, 1.006], 4500 : [0.994, 1.005], 5000 : [0.995, 1.006]}
    eff_res = {800 : [0.956, 1.075], 1000 : [0.960, 1.062], 1200 : [0.962, 1.064], 1400 : [0.964, 1.065], 1600 : [0.959, 1.066], 1800 : [0.964, 1.068], 2000 : [0.962, 1.074], 2500 : [0.964, 1.078], 3000 : [0.958, 1.083], 3500 : [0.957, 1.082], 4000 : [0.958, 1.082], 4500 : [0.980, 1.088], 5000 : [0.980, 1.089]}
  elif cut=='eebb':
    eff_scale_vv = [0.994,1.006]
    eff_res_vv = [0.936,1.095]
  elif cut=='nnbbVBF':
    eff_scale_vv = [0.996,1.004]
    eff_res_vv = [0.990,1.066]
  elif cut=='eebbVBF':
    eff_scale_vv = [0.997,1.008]
    eff_res_vv = [0.880,1.091]
  
  for i,signal in enumerate(signal_list):
    scale_upline.SetPoint(i,signal,eff_scale[signal][1])
    scale_downline.SetPoint(i,signal,eff_scale[signal][0])
  c1 = TCanvas("c1","sf plot", 800,600)
  scale_upline.SetTitle("H mass scale %s" %cut)
  scale_upline.GetXaxis().SetTitle("m_{X} (GeV)")
  scale_upline.GetYaxis().SetTitle("Uncertainty")
  scale_upline.GetYaxis().SetRangeUser(0.8,1.2)
  scale_upline.SetLineColor(4)
  scale_downline.SetLineColor(2)
  scale_upline.SetLineWidth(3)
  scale_downline.SetLineWidth(3)
  leg = TLegend(0.7,0.75,0.85,0.85)
  leg.AddEntry(scale_upline,"+ 1 s. d.","L")
  leg.AddEntry(scale_downline,"- 1 s. d.","L")
  scale_upline.Draw("AC")
  scale_downline.Draw("C SAME")
  leg.Draw("SAME")
  c1.SaveAs("./plotsSF/jms_%s.pdf" % cut)
  c1.SaveAs("./plotsSF/jms_%s.png" % cut)
  c1.Close() 
  for i,signal in enumerate(signal_list):
    res_upline.SetPoint(i,signal,eff_res[signal][1])
    res_downline.SetPoint(i,signal,eff_res[signal][0])
  c2 = TCanvas("c1","sf plot", 800,600)
  res_upline.SetTitle("H mass resolution %s" %cut)
  res_upline.GetXaxis().SetTitle("m_{X} (GeV)")
  res_upline.GetYaxis().SetTitle("Uncertainty")
  res_upline.GetYaxis().SetRangeUser(0.8,1.2)
  res_upline.SetLineColor(4)
  res_downline.SetLineColor(2)
  res_upline.SetLineWidth(3)
  res_downline.SetLineWidth(3)
  leg = TLegend(0.7,0.75,0.85,0.85)
  leg.AddEntry(res_upline,"+ 1 s. d.","L")
  leg.AddEntry(res_downline,"- 1 s. d.","L")
  res_upline.Draw("AC")
  res_downline.Draw("C SAME")
  leg.Draw("SAME")
  c2.SaveAs("./plotsSF/jmr_%s.pdf" % cut)
  c2.SaveAs("./plotsSF/jmr_%s.png" % cut)
  c2.Close()   


#draweffb('2b')
#draweffb('0b')
drawjec('nnbb')
drawjmc('nnbb')

