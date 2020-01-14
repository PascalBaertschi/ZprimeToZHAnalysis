    
variable = {
    'H_pt': {
        'title' : "AK8 jet p_{T} (GeV)",
        'nbins' : 36,
        'min' : 200,
        'max' : 2000,
        'log' : True,
    },
    'H_eta': {
        'title' : "AK8 jet #eta",
        'nbins' : 36,
        'min' : 0.,
        'max' : 5.,
        'log' : False,
    },
    'H_phi': {
        'title' : "AK8 jet #phi",
        'nbins' : 36,
        'min' : -3.15,
        'max' : 3.15,
        'log' : False,
    },
    'H_M': {
        'title' : "AK8 mass (GeV)",
        'nbins' : 50,
        'min' : 0,
        'max' : 250,
        'log' : False,
    },
    'H_mass': {
        'title' : "jet soft drop mass (GeV)",
        'nbins' : 44,
        'min' : 30,
        'max' : 250,
        'log' : False,
    },
    'H_tau21': {
        'title' : "AK8 jet #tau_{21}",
        'nbins' : 50,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'H_tau31': {
        'title' : "AK8 jet #tau_{31}",
        'nbins' : 50,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'H_tau32': {
        'title' : "AK8 jet #tau_{32}",
        'nbins' : 50,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'H_tau41': {
        'title' : "AK8 jet #tau_{41}",
        'nbins' : 50,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'H_tau42': {
        'title' : "AK8 jet #tau_{42}",
        'nbins' : 50,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'H_ddt': {
        'title' : "AK8 jet DDT",
        'nbins' : 50,
        'min' : -0.2,
        'max' : 1.2,
        'log' : False,
    },
    'H_tau21_ddt': {
        'title' : "AK8 jet DDT",
        'nbins' : 50,
        'min' : -0.2,
        'max' : 1.2,
        'log' : False,    
    },
    'H_csv1': {
        'title' : "AK8 jet subjet1 CSV",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'H_csv2': {
        'title' : "AK8 jet subjet2 CSV",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'H_deepcsv1': {
        'title' : "AK8 jet subjet1 DeepCSV",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'H_deepcsv2': {
        'title' : "AK8 jet subjet2 DeepCSV",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'BtagDeepB': {
        'title' : "AK8 jet btagDeepB",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'BtagHbb': {
        'title' : "AK8 jet btagHbb",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'DeepTagMD_H4qvsQCD': {
        'title' : "AK8 jet btagH4qvsQCD",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'DeepTagMD_HbbvsQCD': {
        'title' : "AK8 jet btagHbbvsQCD",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'DeepTagMD_ZHbbvsQCD': {
        'title' : "AK8 jet btagZHbbvsQCD",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'DeepTagMD_ZbbvsQCD': {
        'title' : "AK8 jet btagZbbvsQCD",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'DeepTagMD_bbvsLight': {
        'title' : "AK8 jet btagbbvsLight",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'DeepTagMD_WvsQCD': {
        'title' : "AK8 jet btagWvsQCD",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'DeepTagMD_ZvsQCD': {
        'title' : "AK8 jet btagZvsQCD",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'QCDNLO_Corr': {
        'title' : "QCDNLO_Corr",
        'nbins' : 25,
        'min' : 0.5,
        'max' : 1.1,
        'log' : False,
    },
    'QCDNNLO_Corr': {
        'title' : "QCDNNLO_Corr",
        'nbins' : 25,
        'min' : 0.5,
        'max' : 1.1,
        'log' : False,
    },
    'EWKNLO_Corr': {
        'title' : "EWKNLO_Corr",
        'nbins' : 50,
        'min' : 0.75,
        'max' : 1.05,
        'log' : False,
    },
    'eecutflow': {
        'title' : "",
        'nbins' : 10,
        'min' : 0,
        'max' : 10,
        'log' : True,
    },
    'mmcutflow': {
        'title' : "",
        'nbins' : 10,
        'min' : 0,
        'max' : 10,
        'log' : True,
    },
    'nncutflow': {
        'title' : "",
        'nbins' : 10,
        'min' : 0,
        'max' : 10,
        'log' : True,
    },
    'H_ntag': {
        'title' : "AK8 jet number of b-tagged subjets",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : True,
    },
    'H_dbt': {
        'title' : "AK8 jet double b-tagger",
        'nbins' : 50,
        'min' : -1,
        'max' : 1,
        'log' : False,
    },
    'H_chf': {
        'title' : "AK8 jet charged hadron fraction",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'H_nhf': {
        'title' : "AK8 jet neutral hadron fraction",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'H_phf': {
        'title' : "AK8 jet photon hadron fraction",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'H_chm': {
        'title' : "AK8 jet charged hadron multiplicity",
        'nbins' : 100,
        'min' : 0,
        'max' : 100,
        'log' : False,
    },
    'H_ptB': {
        'title' : "AK8 jet subjet p_{T} balance",
        'nbins' : 50,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'H_hadronflavour': {
        'title' : "H hadron flavour",
        'nbins' : 7,
        'min' : 0,
        'max' : 6,
        'log' : False,
    },
    'Mu1_pt': {
        'title' : "muon 1 p_{T} (GeV)",
        'nbins' : 36,
        'min' : 0,
        'max' : 2000,
        'log' : True,
    },
    'Mu1_eta': {
        'title' : "muon 1 #eta",
        'nbins' : 36,
        'min' : -3,
        'max' : 3,
        'log' : False,
    },
    'Mu1_relIso': {
        'title' : "muon 1 relIso",
        'nbins' : 36,
        'min' : 0,
        'max' : 5,
        'log' : False,
    },
    'Mu2_pt': {
        'title' : "muon 2 p_{T} (GeV)",
        'nbins' : 36,
        'min' : 0,
        'max' : 1000,
        'log' : True,
    },
    'Mu2_eta': {
        'title' : "muon 2 #eta",
        'nbins' : 36,
        'min' : -3,
        'max' : 3,
        'log' : False,
    },
    'Mu2_relIso': {
        'title' : "muon 2 relIso",
        'nbins' : 36,
        'min' : 0,
        'max' : 5,
        'log' : False,
    },
    'Ele1_pt': {
        'title' : "electron 1 p_{T} (GeV)",
        'nbins' : 36,
        'min' : 0,
        'max' : 1000,
        'log' : True,
    },
    'Ele1_eta': {
        'title' : "electron 1 #eta",
        'nbins' : 36,
        'min' : -3,
        'max' : 3,
        'log' : False,
    },
    'Ele2_pt': {
        'title' : "electron 2 p_{T} (GeV)",
        'nbins' : 36,
        'min' : 0,
        'max' : 1000,
        'log' : True,
    },
    'Ele2_eta': {
        'title' : "electron 2 #eta",
        'nbins' : 36,
        'min' : -3,
        'max' : 3,
        'log' : False,
    },
    'V_pt': {
        'title' : "Z p_{T} (GeV)",
        'nbins' : 36,
        'min' : 200,
        'max' : 2000,
        'log' : True,
    },
    'V_eta': {
        'title' : "Z #eta",
        'nbins' : 36,
        'min' : 0.,
        'max' : 5.,
        'log' : False,
    },
    'V_phi': {
        'title' : "Z #phi",
        'nbins' : 36,
        'min' : -3.15,
        'max' : 3.15,
        'log' : False,
    },
    'V_mass': {
        'title' : "Z mass (GeV)",
        'nbins' : 75,
        'min' : 70,
        'max' : 110,
        'log' : False,
    },
    'MET': {
        'title' : "#slash{E}_{T} (GeV)",
        'nbins' : 36,
        'min' : 0,
        'max' : 400,
        'log' : True,
    },
    'V_tau21': {
        'title' : "V jet #tau_{21}",
        'nbins' : 50,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'V_ddt': {
        'title' : "V jet DDT",
        'nbins' : 50,
        'min' : -0.2,
        'max' : 1.2,
        'log' : False,
    },
    'V_csv1': {
        'title' : "V jet subjet1 CSV",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'V_csv2': {
        'title' : "V jet subjet2 CSV",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'V_dbt': {
        'title' : "V jet double b-tagger",
        'nbins' : 50,
        'min' : -1,
        'max' : 1,
        'log' : False,
    },
    'VH_deltaR': {
        'title' : "#Delta R (Z, AK8 jet)",
        'nbins' : 25,
        'min' : 0,
        'max' : 5,
        'log' : False,
    },
    'MaxBTag': {
        'title' : "max DeepCSV AK4 jets",
        'nbins' : 25,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'MinDPhi': {
        'title' : "min #Delta #varphi (jet-#slash{E}_{T})",
        'nbins' : 28,
        'min' : 0,
        'max' : 3.15,
        'log' : True,
    },
    'DEta': {
        'title' : "#Delta #eta (Z-AK8 jet)",
        'nbins' : 25,
        'min' : 0,
        'max' : 5.,
        'log' : False,
    },
    'DPhi': {
        'title' : "#Delta #varphi (Z-AK8 jet)",
        'nbins' : 28,
        'min' : 0,
        'max' : 3.15,
        'log' : False,
    },
     'nFatJets': {
        'title' : "number of AK8 jets",
        'nbins' : 6,
        'min' : -0.5,
        'max' : 6.5,
        'log' : True,
    },
    'nJets': {
        'title' : "number of AK4 jets",
        'nbins' : 6,
        'min' : -0.5,
        'max' : 6.5,
        'log' : True,
    },
    'nElectrons': {
        'title' : "number of e",
        'nbins' : 4,
        'min' : -0.5,
        'max' : 4.5,
        'log' : True,
    },
    'nMuons': {
        'title' : "number of #mu",
        'nbins' : 4,
        'min' : -0.5,
        'max' : 4.5,
        'log' : True,
    },
    'nTaus': {
        'title' : "number of #tau",
        'nbins' : 4,
        'min' : -0.5,
        'max' : 4.5,
        'log' : True,
    },
    'nPV': {
        'title' : "number of primary vertices",
        'nbins' : 51,
        'min' : -0.5,
        'max' : 50.5,
        'log' : True,
    },
    'X_pt': {
        'title' : "X p_{T} (GeV)",
        'nbins' : 36,
        'min' : 0,
        'max' : 2000,
        'log' : True,
    },
    'X_eta': {
        'title' : "X #eta",
        'nbins' : 36,
        'min' : 0.,
        'max' : 5.,
        'log' : False,
    },
    'X_phi': {
        'title' : "X #phi",
        'nbins' : 36,
        'min' : -3.15,
        'max' : 3.15,
        'log' : False,
    },
    'X.M()': {
        'title' : "m_{X} (GeV)",
        'nbins' : 80,
        'min' : 750,
        'max' : 5750,
        'log' : True,
    },
    'X_mass': {
        'title' : "m_{X} (GeV)",
        'nbins' : 40, #10 15
        'min' : 750,
        'max' : 5750, #2750 3750 6750
        'log' : True,
    },
    'CosThetaStar': {
        'title' : "cos(#theta*)",
        'nbins' : 50,
        'min' : -1,
        'max' : 1,
        'log' : False,
    },
    'CosTheta1': {
        'title' : "cos(#theta_{1})",
        'nbins' : 50,
        'min' : -1,
        'max' : 1,
        'log' : False,
    },
    'CosTheta2': {
        'title' : "cos(#theta_{2})",
        'nbins' : 50,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'V.Pt()/X_mass': {
        'title' : "p_{T}^{V}/m_{X}",
        'nbins' : 50,
        'min' : 0,
        'max' : 1,
        'log' : False,
    },
    'isBoostedTau': {
        'title' : "overlap with VH #rightarrow #tau#tau",
        'nbins' : 2,
        'min' : 0.5,
        'max' : 1.5,
        'log' : True,
    },
    'isBoosted4B': {
        'title' : "overlap with HH #rightarrow 4b",
        'nbins' : 2,
        'min' : 0.5,
        'max' : 1.5,
        'log' : True,
    },
    'dijet_VBF_mass':{
        'title' : "dijet mass VBF",
         'nbins' : 36,
         'min'  : 200,
         'max' : 2000,
         'log' : False,
    },
    'deltaR_VBF' :{
        'title' : "#Delta R(jet1,jet2)",
        'nbins' : 50,
        'min' : 0.,
        'max' : 10.,
        'log' : False,
    },
    'deltaR_HVBFjet1': {
        'title' : "#DeltaR(AK8 jet,VBF jet",
        'nbins' : 50,
        'min' : 0.,
        'max' : 10.,
        'log' : False,
    },
    'deltaR_HVBFjet2': {
        'title' : "#DeltaR(AK8 jet, VBF jet)",
        'nbins' : 50,
        'min' : 0.,
        'max' : 10.,
        'log' : False,
    },
    'PrefireWeight' : {
        'title' : "PrefireWeight",
        'nbins' : 50,
        'min' : 0.,
        'max' : 2.,
        'log' : False,
    },
    'LeptonWeight': {
        'title' : "LeptonWeight",
        'nbins' : 50,
        'min' : 0.,
        'max' : 5.,
        'log' : False,
    },
    'LeptonWeightUp': {
        'title' : "LeptonWeight up",
        'nbins' : 50,
        'min' : 0.,
        'max' : 5.,
        'log' : False,
    },
    'LeptonWeightDown': {
        'title' : "LeptonWeight down",
        'nbins' : 50,
        'min' : 0.,
        'max' : 5.,
        'log' : False,
    },
    'TriggerWeight': {
        'title' : "TriggerWeight",
        'nbins' : 50,
        'min' : 0.,
        'max' : 5.,
        'log' : False,
    },
    'TriggerWeightUp': {
        'title' : "TriggerWeight up",
        'nbins' : 50,
        'min' : 0.,
        'max' : 5.,
        'log' : False,
    },
    'TriggerWeightDown': {
        'title' : "TriggerWeight down",
        'nbins' : 50,
        'min' : 0.,
        'max' : 5.,
        'log' : False,
    },
    'BTagAK4Weight_deep': {
        'title' : "BTagAK4Weight deep",
        'nbins' : 50,
        'min' : 0.,
        'max' : 5.,
        'log' : False,
    },
    'BTagAK4Weight_deep_up': {
        'title' : "BTagAK4Weight deep up",
        'nbins' : 50,
        'min' : 0.,
        'max' : 5.,
        'log' : False,
    },
    'BTagAK4Weight_deep_down': {
        'title' : "BTagAK4Weight deep down",
        'nbins' : 50,
        'min' : 0.,
        'max' : 5.,
        'log' : False,
   },
    'BTagAK8Weight_deep': {
        'title' : "BTagAK8Weight deep",
        'nbins' : 50,
        'min' : 0.,
        'max' : 5.,
        'log' : False,
    },
    'BTagAK8Weight_deep_up': {
        'title' : "BTagAK8Weight deep up",
        'nbins' : 50,
        'min' : 0.,
        'max' : 5.,
        'log' : False,
    },
    'BTagAK8Weight_deep_down': {
        'title' : "BTagAK8Weight deep down",
        'nbins' : 50,
        'min' : 0.,
        'max' : 5.,
        'log' : False,
   }
}

