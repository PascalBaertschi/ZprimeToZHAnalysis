#! /usr/bin/env python

#(abs(lepton1_eta)<1.4442 || abs(lepton1_eta)>1.566) &&

selection = {
   #bb
   "nnbb" : "isZtoNN && !isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750 && !isMaxBTag_loose && isHtobb",
   "eebb" : "isZtoEE && !isVBF && DEta<1.7 && X_mass > 750 && VH_deltaR >2 && isHtobb",
   "mmbb" : "isZtoMM && !isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && DEta<1.7 && X_mass > 750 && VH_deltaR >2 && isHtobb",
   #bb in sideband
   "nnbbSB" : "isZtoNN && !isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750 && !isMaxBTag_loose && isHtobb && ((H_mass>30. && H_mass<105.) || H_mass>135.)",
   "eebbSB" : "isZtoEE && !isVBF && DEta<1.7 && X_mass > 750 && VH_deltaR > 2 && isHtobb && ((H_mass>30. && H_mass<105.) || H_mass>135.)",
   "mmbbSB" : "isZtoMM && !isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && DEta<1.7 && X_mass > 750 && VH_deltaR > 2 && isHtobb && ((H_mass>30. && H_mass<105.) || H_mass>135.)",
   #bb in signal range
   "nnbbSR" : "isZtoNN && !isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750  && !isMaxBTag_loose && isHtobb  && (H_mass>105. && H_mass<135.)",
   "eebbSR" : "isZtoEE && !isVBF && DEta<1.7 && X_mass > 750 && VH_deltaR > 2  && isHtobb && (H_mass>105. && H_mass<135.)",
   "mmbbSR" : "isZtoMM && !isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && DEta<1.7 && X_mass > 750 && VH_deltaR > 2 && isHtobb && (H_mass>105. && H_mass<135.)",

   ### no btag channel
   "nn0b" : "isZtoNN && !isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750 && !isMaxBTag_loose && !isHtobb",
   "ee0b" : "isZtoEE && !isVBF && DEta<1.7 && X_mass > 750 && VH_deltaR > 2 && !isHtobb",
   "mm0b" : "isZtoMM && !isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && DEta<1.7 && X_mass > 750 && VH_deltaR > 2 && !isHtobb",
   #0b in sideband
   "nn0bSB" : "isZtoNN && !isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750 &&  !isMaxBTag_loose && !isHtobb && ((H_mass>30. && H_mass<105.) || H_mass>135.)",
   "ee0bSB" : "isZtoEE && !isVBF && DEta<1.7 && X_mass > 750 && VH_deltaR > 2 && !isHtobb && ((H_mass>30. && H_mass<105.) || H_mass>135.)",
   "mm0bSB" : "isZtoMM && !isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && DEta<1.7 && X_mass > 750 && VH_deltaR > 2 && !isHtobb && ((H_mass>30. && H_mass<105.) || H_mass>135.) ",
   #0b in signal range
   "nn0bSR" : "isZtoNN && !isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750 && !isMaxBTag_loose && !isHtobb && (H_mass>105. && H_mass<135.)",
   "ee0bSR" : "isZtoEE && !isVBF && DEta<1.7 && X_mass > 750 && VH_deltaR > 2 && !isHtobb && (H_mass>105. && H_mass<135.)",
   "mm0bSR" : "isZtoMM && !isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && DEta<1.7 && X_mass > 750 && VH_deltaR > 2 && !isHtobb && (H_mass>105. && H_mass<135.)",

   #inclusive
   "nninc" : "isZtoNN && !isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750",
   "eeinc" : "isZtoEE && !isVBF && DEta<1.7 && X_mass > 750 && VH_deltaR > 2",
   "mminc" : "isZtoMM && !isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && DEta<1.7 && X_mass > 750 && VH_deltaR > 2",
   #inclusive sideband
   "nnincSB" : "isZtoNN && !isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750 && ((H_mass>30. && H_mass<105.) || H_mass>135.)",
   "eeincSB" : "isZtoEE && !isVBF && DEta<1.7 && X_mass > 750 && VH_deltaR >2 &&  ((H_mass>30. && H_mass<105.) || H_mass>135.)",
   "mmincSB" : "isZtoMM && !isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && DEta<1.7 && X_mass > 750 && VH_deltaR >2 && ((H_mass>30. && H_mass<105.) || H_mass>135.)",
   #inclusive signal range
   "nnincSR" : "isZtoNN && !isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750 && (H_mass>105. && H_mass<135.)",
   "eeincSR" : "isZtoEE && !isVBF && DEta<1.7 && X_mass > 750 && VH_deltaR >2 && (H_mass>105. && H_mass<135.)",
   "mmincSR" : "isZtoMM && !isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && DEta<1.7 && X_mass > 750 && VH_deltaR > 2 && (H_mass>105. && H_mass<135.)",

   #top control region
   "nnbbTR" : "isZtoNN && !isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750  && isMaxBTag_tight && isHtobb && (H_mass>30. && H_mass<250.)",
   "embbTR" : "isTtoEM && !isVBF && V_mass>110. && V_pt>120. && isHtobb && (H_mass>30. && H_mass<250.)",
   "nn0bTR" : "isZtoNN && !isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750 && isMaxBTag_tight && !isHtobb && (H_mass>30. && H_mass<250.)",
   "em0bTR" : "isTtoEM && !isVBF && V_mass>110. && V_pt>120. && !isHtobb && (H_mass>30. && H_mass<250.)",
   "nnbbIR" : "isZtoNN && !isVBF && MinDPhi>0.5 && DPhi>2 && X_mass > 750 && VH_deltaR >2 && isMaxBTag_medium && !isMaxBTag_tight && isHtobb && (H_mass>30. && H_mass<250.)",
   "nn0bIR" : "isZtoNN && !isVBF && MinDPhi>0.5 && DPhi>2 && X_mass > 750 && VH_deltaR >2 && isMaxBTag_medium && !isMaxBTag_tight && !isHtobb && (H_mass>30. && H_mass<250.)",


   ###  VBF  ###
   #bb
   "nnbbVBF" : "isZtoNN && isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750 && !isMaxBTag_loose && isHtobb",
   "eebbVBF" : "isZtoEE && isVBF && X_mass > 750 && VH_deltaR > 2 && isHtobb",
   "mmbbVBF" : "isZtoMM && isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && X_mass > 750 && VH_deltaR > 2 && isHtobb", 
   #bb in sideband
   "nnbbVBFSB" : "isZtoNN && isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750 && !isMaxBTag_loose && isHtobb && ((H_mass>30. && H_mass<105.) || H_mass>135.)",
   "eebbVBFSB" : "isZtoEE && isVBF && X_mass > 750 && VH_deltaR > 2 && isHtobb && ((H_mass>30. && H_mass<105.) || H_mass>135.)",
   "mmbbVBFSB" : "isZtoMM && isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && X_mass > 750 && VH_deltaR > 2 && isHtobb && ((H_mass>30. && H_mass<105.) || H_mass>135.)",
   #bb in signal range
   "nnbbVBFSR" : "isZtoNN && isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750 && !isMaxBTag_loose && isHtobb  && (H_mass>105. && H_mass<135.)",
   "eebbVBFSR" : "isZtoEE && isVBF && X_mass > 750 && VH_deltaR > 2  && isHtobb && (H_mass>105. && H_mass<135.)",
   "mmbbVBFSR" : "isZtoMM && isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && X_mass > 750 && VH_deltaR > 2 && isHtobb && (H_mass>105. && H_mass<135.)",

   ### no btag channel
   "nn0bVBF" : "isZtoNN && isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750 && !isMaxBTag_loose && !isHtobb",
   "ee0bVBF" : "isZtoEE && isVBF && X_mass > 750 && VH_deltaR > 2 && !isHtobb",
   "mm0bVBF" : "isZtoMM && isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && X_mass > 750 && VH_deltaR > 2 && !isHtobb",
   #0b in sideband
   "nn0bVBFSB" : "isZtoNN && isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750 && !isMaxBTag_loose && !isHtobb && ((H_mass>30. && H_mass<105.) || H_mass>135.)",
   "ee0bVBFSB" : "isZtoEE && isVBF && X_mass > 750 && VH_deltaR > 2 && !isHtobb && ((H_mass>30. && H_mass<105.) || H_mass>135.)",
   "mm0bVBFSB" : "isZtoMM && isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && X_mass > 750 && VH_deltaR > 2 && !isHtobb && ((H_mass>30. && H_mass<105.) || H_mass>135.) ",   #0b in signal range
   "nn0bVBFSR" : "isZtoNN && isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750  && !isMaxBTag_loose && !isHtobb && (H_mass>105. && H_mass<135.)",
   "ee0bVBFSR" : "isZtoEE && isVBF && X_mass > 750 && VH_deltaR > 2 && !isHtobb && (H_mass>105. && H_mass<135.)",
   "mm0bVBFSR" : "isZtoMM && isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && X_mass > 750 && VH_deltaR > 2 && !isHtobb && (H_mass>105. && H_mass<135.)",
 
   #inclusive
   #"nnincVBF" : "isZtoNN && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750",
   #"eeincVBF" : "isZtoEE && X_mass > 750 && VH_deltaR > 2",
   #"mmincVBF" : "isZtoMM && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && X_mass > 750 && VH_deltaR > 2",
   "nnincVBF" : "isZtoNN && isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750",
   "eeincVBF" : "isZtoEE && isVBF && X_mass > 750 && VH_deltaR > 2",
   "mmincVBF" : "isZtoMM && isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && X_mass > 750 && VH_deltaR > 2",
   #inclusive sideband
   "nnincVBFSB" : "isZtoNN && isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750 && ((H_mass>30. && H_mass<105.) || H_mass>135.)",
   "eeincVBFSB" : "isZtoEE && isVBF && X_mass > 750 && VH_deltaR >2 &&  ((H_mass>30. && H_mass<105.) || H_mass>135.)",
   "mmincVBFSB" : "isZtoMM && isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && X_mass > 750 && VH_deltaR >2 && ((H_mass>30. && H_mass<105.) || H_mass>135.)",
   #inclusive signal range
   "nnincVBFSR" : "isZtoNN && isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750  && (H_mass>105. && H_mass<135.)",
   "eeincVBFSR" : "isZtoEE && isVBF && X_mass > 750 && VH_deltaR >2 && (H_mass>105. && H_mass<135.)",
   "mmincVBFSR" : "isZtoMM && isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && X_mass > 750 && VH_deltaR > 2 && (H_mass>105. && H_mass<135.)",

   #top control region
   "nnbbVBFTR" : "isZtoNN && isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750 && isMaxBTag_tight && isHtobb && (H_mass>30. && H_mass<250.)",
   "embbVBFTR" : "isTtoEM && isVBF && V_mass>110. && V_pt>120. && isHtobb && (H_mass>30. && H_mass<250.)",
   "nn0bVBFTR" : "isZtoNN && isVBF && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750 && isMaxBTag_tight && !isHtobb && (H_mass>30. && H_mass<250.)",
   "em0bVBFTR" : "isTtoEM && isVBF && V_mass>110. && V_pt>120. && !isHtobb && (H_mass>30. && H_mass<250.)",

   #special
   "eebbDEtacheck" : "isZtoEE && !isVBF && X_mass > 750 && VH_deltaR > 2  && isHtobb",
   "mmbbDEtacheck" : "isZtoMM && !isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && X_mass > 750 && VH_deltaR > 2 && isHtobb",

   "eebbVBFDEtacheck" : "isZtoEE && isVBF && X_mass > 750 && VH_deltaR >2 && isHtobb",
   "mmbbVBFDEtacheck" : "isZtoMM && isVBF && Mu1_relIso!=-1. && Mu1_relIso < 0.15 && Mu2_relIso!=-1. && Mu2_relIso < 0.15 && X_mass > 750 && VH_deltaR >2 && isHtobb",

   #inclusive signal range only for q signal
   #"nnincSR_sign" : "isZtoNN && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750  && (H_mass>105. && H_mass<135.) && !isHtobb",
   #"eeincSR_sign" : "isZtoEE && DEta<1.7 && X_mass > 750 && VH_deltaR >2 && (H_mass>105. && H_mass<135.) && !isHtobb",
   #"mmincSR_sign" : "isZtoMM && DEta<1.7 && X_mass > 750 && VH_deltaR >2 && (H_mass>105. && H_mass<135.) && !isHtobb",
   #inclusive signal range only for bb signal
   "nnincSR_sign" : "isZtoNN && MinDPhi>0.5 && DPhi>2 && H_chf>0.1 && (V_pt/H_pt)>0.6 && X_mass > 750 && (H_mass>105. && H_mass<135.) && H_hadronflavour==5",
   "eeincSR_sign" : "isZtoEE && DEta<1.7 && X_mass > 750 && VH_deltaR >2 && (H_mass>105. && H_mass<135.) && H_hadronflavour==5",
   "mmincSR_sign" : "isZtoMM && DEta<1.7 && X_mass > 750 && VH_deltaR >2 && (H_mass>105. && H_mass<135.) && H_hadronflavour==5",
   #ranges
   "SR" : " && (H_mass>105. && H_mass<135.)",
   "VR" : " && (H_mass>65. && H_mass<105.)",
   "SB" : " && ((H_mass>30. && H_mass<65.) || H_mass>135.)",
}
