import os, re
from ROOT import TFile,TCanvas,gStyle
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

path = '/work/pbaertsc/heavy_resonance/NanoTreeProducer/CorrectionTools/leptonEfficiencies/'

def draweffplot(filename1, histname1, name1, ptvseta):
    if 'Electron' in filename1:
      lepton = "electron"
    else:
      lepton = "muon"
    file1     = ensureTFile(filename1)
    hist1     = file1.Get(histname1)
    hist = hist1
    hist.SetDirectory(0)
    file1.Close()
    if '2016' in filename1:
      year = '2016'
    elif '2017' in filename1:
      year = '2017'
    elif '2018' in filename1:
      year = '2018'
    c1 = TCanvas("c1","sf plot", 800,600)
    hist.SetOption("COLZ")
    hist.SetTitle("%s trigger scale factor %s" %(lepton,year))
    if ptvseta:
      hist.GetXaxis().SetTitle("#eta")
      hist.GetYaxis().SetTitle("p_{T}")
      if lepton =="muon":
        if year == '2016':
          hist.GetYaxis().SetRangeUser(0,800)
        else:
          hist.GetYaxis().SetRangeUser(0,1200)
      else:
        if year == '2017':
          hist.GetYaxis().SetRangeUser(0,1000)
        else:
          hist.GetYaxis().SetRangeUser(0,1000)
      hist.GetYaxis().SetTitleOffset(1.2)
    else:
      hist.GetXaxis().SetTitle("p_{T}")
      hist.GetYaxis().SetTitle("#eta")
    hist.GetZaxis().SetRangeUser(0.7,1.0)
    hist.Draw()
    gStyle.SetOptStat(0)
    c1.SaveAs("./plotsSF/%s_trigger_sf_%s.pdf" %(lepton,year))
    c1.SaveAs("./plotsSF/%s_trigger_sf_%s.png" %(lepton,year))
    c1.Close()


def draweffplot_2files(filename1, histname1, name1, filename2, histname2, name2, ptvseta):
    if 'Electron' in filename1:
      lepton = "electron"
    else:
      lepton = "muon"
    file1     = ensureTFile(filename1)
    hist1     = file1.Get(histname1)
    file2     = ensureTFile(filename2)
    hist2     = file2.Get(histname2)
    hist1.Divide(hist2)
    hist = hist1
    hist.SetDirectory(0)
    file1.Close()
    file2.Close()
    if '2016' in filename1:
      year = '2016'
    elif '2017' in filename1:
      year = '2017'
    elif '2018' in filename1:
      year = '2018' 
    c1 = TCanvas("c1","sf plot", 800,600)
    hist.SetOption("COLZ")
    hist.SetTitle("%s trigger scale factor %s" %(lepton,year))
    if ptvseta:
      hist.GetXaxis().SetTitle("#eta")
      hist.GetYaxis().SetTitle("p_{T}")
      hist.GetYaxis().SetTitleOffset(1.2)
    else:
      hist.GetXaxis().SetTitle("p_{T}")
      hist.GetYaxis().SetTitle("#eta")
    hist.GetZaxis().SetRangeUser(0.7,1.0)
    hist.Draw()
    gStyle.SetOptStat(0)
    c1.SaveAs("./plotsSF/%s_trigger_sf_%s.pdf" %(lepton,year))
    c1.SaveAs("./plotsSF/%s_trigger_sf_%s.png" %(lepton,year))
    c1.Close()

  
    
    

#draweffplot_2files(path+"ElectronPOG/Run2016/Ele115_passingTight_2016.root",'EGamma_EffData2D_ABSOLUTE','ele_trig_data', path+"ElectronPOG/Run2016/Ele115_passingTight_2016.root",'EGamma_EffMC2D_ABSOLUTE','ele_trig_mc',ptvseta=True)
#draweffplot_2files(path+"ElectronPOG/Run2017/Ele115orEle35_SF_2017.root",'ELE_DATA_ABSOLUTE','ele_trig_data', path+"ElectronPOG/Run2017/Ele115orEle35_SF_2017.root",'ELE_MC_ABSOLUTE','ele_trig_mc',ptvseta=False)
#draweffplot_2files(path+"ElectronPOG/Run2018/Ele115orEle35_SF_2018.root",'ELE_DATA_ABSOLUTE','ele_trig_data', path+"ElectronPOG/Run2018/Ele115orEle35_SF_2018.root",'ELE_MC_ABSOLUTE','ele_trig_mc',ptvseta=False)
draweffplot(path+"ElectronPOG/Run2016/Ele115orEleIso27orPho175_SF_2016.root","SF_TH2F","el_trig",ptvseta=True)
draweffplot(path+"ElectronPOG/Run2017/Ele115orEleIso35orPho200_SF_2017.root","SF_TH2F","el_trig",ptvseta=True)
draweffplot(path+"ElectronPOG/Run2018/Ele115orEleIso32orPho200_SF_2018.root","SF_TH2F","el_trig",ptvseta=True)
draweffplot(path+"MuonPOG/Run2016/SingleMuonTrigger_2016.root","Mu50_OR_TkMu50_PtEtaBins/abseta_pt_ratio","mu_trig",ptvseta=True)
draweffplot(path+"MuonPOG/Run2017/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root","Mu50_PtEtaBins/abseta_pt_ratio",'mu_trig',ptvseta=True)
draweffplot(path+"MuonPOG/Run2018/EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root","Mu50_OR_OldMu100_OR_TkMu100_PtEtaBins/abseta_pt_ratio",'mu_trig',ptvseta=True)

