{
  //load
  gROOT->Reset();
  gROOT->ProcessLine("#include <utility>");
  string path("/afs/cern.ch/user/y/yohay/CMSSW_5_2_4_patch3/src/BoostedTauAnalysis/TauAnalyzer/");
  path+="test/";
  gSystem->Load((path + "STLDictionary.so").c_str());
//   gROOT->LoadMacro((path + "Plot.C++").c_str());
  gSystem->Load((path + "Plot_C.so").c_str());

  //unit strings
  string unitPTTau("Visible gen p_{T}^{#tau} (GeV)");
  string unitPTMu("Gen p_{T}^{#mu} (GeV)");
  string unitEtaTau("Visible gen #eta^{#tau}");
  string unitEtaMu("Gen #eta^{#mu}");
  string unitDR("#DeltaR(visible gen #tau, gen #mu)");
  string noUnit("");

  //map of bin labels for certain efficiency plots
  vector<string> binLabels;
  binLabels.push_back("kNull");
  binLabels.push_back("kOneProng0PiZero");
  binLabels.push_back("kOneProng1PiZero");
  binLabels.push_back("kOneProng2PiZero");
  binLabels.push_back("kOneProng3PiZero");
  binLabels.push_back("kOneProngNPiZero");
  binLabels.push_back("kTwoProng0PiZero");
  binLabels.push_back("kTwoProng1PiZero");
  binLabels.push_back("kTwoProng2PiZero");
  binLabels.push_back("kTwoProng3PiZero");
  binLabels.push_back("kTwoProngNPiZero");
  binLabels.push_back("kThreeProng0PiZero");
  binLabels.push_back("kThreeProng1PiZero");
  binLabels.push_back("kThreeProng2PiZero");
  binLabels.push_back("kThreeProng3PiZero");
  binLabels.push_back("kThreeProngNPiZero");
  binLabels.push_back("kRareDecayMode");
  map<string, vector<string> > binLabelMap;
  binLabelMap["muHadGenDecayMode"] = binLabels;
  binLabelMap["muHadCorrectRecoDecayModeGenDecayMode"] = binLabels;
  binLabelMap["muHadRecoDecayMode"] = binLabels;
  binLabelMap["muHadGen1ProngRecoDecayMode"] = binLabels;
  binLabelMap["muHadGen1Prong1Pi0RecoDecayMode"] = binLabels;
  binLabelMap["muHadGen3ProngRecoDecayMode"] = binLabels;

  //map of inputs to efficiency histograms
  map<string, pair<string, string> > effHistMap;
  effHistMap["muHadPFTauMatchVisibleGenPT"] = make_pair(string("muHadVisibleGenPT"), unitPTTau);
  effHistMap["muHadPFTauMatchVisibleGenEta"] = make_pair(string("muHadVisibleGenEta"), unitEtaTau);
  effHistMap["muHadMuMatchVisibleGenPT"] = make_pair(string("muHadVisibleGenPT"), unitPTTau);
  effHistMap["muHadMuMatchVisibleGenEta"] = make_pair(string("muHadVisibleGenEta"), unitEtaTau);
  effHistMap["muHadPFTauMatchMuMatchVisibleGenPT"] = 
    make_pair(string("muHadVisibleGenPT"), unitPTTau);
  effHistMap["muHadPFTauMatchMuMatchVisibleGenEta"] = 
    make_pair(string("muHadVisibleGenEta"), unitEtaTau);
  effHistMap["muHadPFTauMatchGenMuPT"] = make_pair(string("muHadGenMuPT"), unitPTMu);
  effHistMap["muHadPFTauMatchGenMuEta"] = make_pair(string("muHadGenMuEta"), unitEtaMu);
  effHistMap["muHadMuMatchGenMuPT"] = make_pair(string("muHadGenMuPT"), unitPTMu);
  effHistMap["muHadMuMatchGenMuEta"] = make_pair(string("muHadGenMuEta"), unitEtaMu);
  effHistMap["muHadPFTauMatchMuMatchGenMuPT"] = make_pair(string("muHadGenMuPT"), unitPTMu);
  effHistMap["muHadPFTauMatchMuMatchGenMuEta"] = make_pair(string("muHadGenMuEta"), unitEtaMu);
  effHistMap["muHadPFTauMatchGenDR"] = make_pair(string("muHadGenDR"), unitDR);
  effHistMap["muHadMuMatchGenDR"] = make_pair(string("muHadGenDR"), unitDR);
  effHistMap["muHadPFTauMatchMuMatchGenDR"] = make_pair(string("muHadGenDR"), unitDR);
  effHistMap["muHadPFTauGenMuMatchGenDR"] = make_pair(string("muHadGenDR"), unitDR);
  effHistMap["muHadPFChargedHadronMatchGenDR"] = make_pair(string("muHadGenDR"), unitDR);
  effHistMap["muHadSharedJetMethod1GenDR"] = make_pair(string("muHadGenDR"), unitDR);
  effHistMap["muHadSharedJetMethod2GenDR"] = make_pair(string("muHadGenDR"), unitDR);
  effHistMap["muHadCorrectRecoDecayModeGenDecayMode"] = 
    make_pair(string("muHadGenDecayMode"), noUnit);
  effHistMap["muHadSharedJetCorrectRecoDecayModeGenDecayMode"] = 
    make_pair(string("muHadGenDecayMode"), noUnit);

  //map of inputs to 1D histograms
  map<string, string> hist1DMap;
  hist1DMap["muHadGenDR"] = unitDR;
  hist1DMap["muHadGenDecayMode"] = noUnit;
  hist1DMap["muHadCorrectRecoDecayModeGenDecayMode"] = noUnit;
  hist1DMap["muHadRecoDecayMode"] = noUnit;
  hist1DMap["muHadGen1ProngRecoDecayMode"] = noUnit;
  hist1DMap["muHadGen1Prong1Pi0RecoDecayMode"] = noUnit;
  hist1DMap["muHadGen3ProngRecoDecayMode"] = noUnit;

//   vector<pair<pair<TFile*, Option_t*>, pair<Color_t, Style_t> > > histMap;
//   histMap.push_back(make_pair(make_pair(), make_pair()));

//   map<pair<string, string>, vector<pair<pair<TFile*, Option_t*>, pair<Color_t, Style_t> > > > 
//     canvasMap;
//   canvasMap[make_pair(string("muHadGen1ProngRecoDecayMode"), noUnit)] = ;
//   canvasMap[make_pair(string("muHadGen1Prong1Pi0RecoDecayMode"), noUnit)] = ;
//   canvasMap[make_pair(string("muHadGen3ProngRecoDecayMode"), noUnit)] = ;

  //plot
  string beginSignal("/data1/yohay/NMSSMHiggs_gg_");
  string endSignal("_genMuPTGt20GeV_genTauPTGt15GeV.root");
  string beginDY("/data1/yohay/Summer12_DYToTauTau_");
  string endDY("_genMuPTGt20GeV_genTauPTGt15GeV.root");
  plotNice(beginSignal + "analysis" + endSignal, effHistMap, binLabelMap, hist1DMap, 
	   beginSignal + "plots" + endSignal, "noPDF");
//   plotNice(beginDY + "analysis" + endDY, effHistMap, binLabelMap, hist1DMap, 
// 	   beginDY + "plots" + endDY, "noPDF");

// plotNiceDifferentFiles(map<pair<string, string>, vector<pair<pair<TFile*, Option_t*>, 
// 			    pair<Color_t, Style_t> > > >& canvasMap, 
// 			    const map<string, vector<string> >& binLabelMap, 
// 			    const string& outputFileName, const string& savePath)
}
