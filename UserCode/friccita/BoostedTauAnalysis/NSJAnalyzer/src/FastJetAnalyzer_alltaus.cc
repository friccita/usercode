// -*- C++ -*-
//
// Package:    FastJetAnalyzer
// Class:      FastJetAnalyzer
// 
/**\class FastJetAnalyzer_alltaus FastJetAnalyzer_alltaus.cc BoostedTauAnalysis/FastJetAnalyzer/src/FastJetAnalyzer_alltaus.cc

 Description: CMSSW Framework version of Thaler's FastJetExample.cc
              to calculate N-subjettiness for all tau decay modes

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesca Ricci-Tam,6 R-025,+41227672274,
//         Created:  Mon Sep 10 10:28:20 CEST 2012
// $Id: FastJetAnalyzer_alltaus.cc,v 1.1 2013/05/06 09:18:52 friccita Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>
#include <TH1.h>
#include <TTree.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "AnalysisDataFormats/TopObjects/interface/CATopJetTagInfo.h"
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
//#include <fastjet/ActiveAreaSpec.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "BoostedTauAnalysis/NSJAnalyzer/interface/Nsubjettiness.h"
#include "BoostedTauAnalysis/NSJAnalyzer/interface/Njettiness.hh"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "BoostedTauAnalysis/Common/interface/GenTauDecayID.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TauReco/interface/PFTauTagInfo.h"
#include "DataFormats/TauReco/interface/PFTauDecayMode.h"
#include "DataFormats/TauReco/interface/PFTauDecayModeAssociation.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "RecoTauTag/TauTagTools/interface/GeneratorTau.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/PseudoJet.hh"

using namespace std;
using namespace fastjet;
using namespace reco;
using namespace pat;
using namespace edm;


//
// class declaration
//

class FastJetAnalyzer_alltaus : public edm::EDAnalyzer {
   public:
      explicit FastJetAnalyzer_alltaus(const edm::ParameterSet&);
      ~FastJetAnalyzer_alltaus();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

  TFile* out_;
  TTree *jetTree;
  TH1F *jetPt_hist;
  TH1F *jetMass_hist;
  TH1I *intime_PU_vertices;
  TH1I* cands_beforeGrooming;
  TH1I* cands_afterGrooming;
  
  TH1F *pjmj_PU0;
  TH1F *pjmj_PUAll;
  TH2F *nSubj_vs_pjmj_PU0;
  TH2F *nSubj_vs_pjmj_PULo;
  TH2F *nSubj_vs_pjmj_PUHi;
  TH2F *nSubj_vs_pjmj_PUAll;

  TH1F *tau1_One;
  TH1I *npIT_One;
  TH1I *cands_One;
  TH1F *pT_One;
  TH1F *tau1_NotOne;
  TH1I *npIT_NotOne;
  TH1I *cands_NotOne;
  TH1F *pT_NotOne;
  TH1F *pT_Groomed;
  TH1F *pT_highest_One;
  TH1F *pT_secondhighest_One;
  TH1F *maxDelR_One;
  TH1F *pT_highest_NotOne;
  TH1F *pT_secondhighest_NotOne;
  TH1F *maxDelR_NotOne;
  TH1F *delR_sisters_One;
  TH1F *delR_sisters_NotOne;
  TH1I *decayModes_One;
  TH1I *decayModes_NotOne;

  TH2F *t3t1_vs_t4t1;
  TH2F *t3t1_vs_t2t1;
  TH2F *t2t1_vs_t4t1;

  TH1F *nSubj1_PU0_hist;
  TH1F *nSubj2_PU0_hist;
  TH1F *nSubj3_PU0_hist;
  TH1F *nSubj4_PU0_hist;

  TH1F *t3t1_PUAll_Ratio;
  TH1F *t3t1_ee;
  TH1F *t3t1_mm;
  TH1F *t3t1_em;
  TH1F *t3t1_eh;
  TH1F *t3t1_mh;
  TH1F *t3t1_hh;

  TH1F *t3t1_PU0_Ratio;
  TH1F *t2t1_PU0_Ratio;
  TH1F *t3t2_PU0_Ratio;
  TH1F *t4t1_PU0_Ratio;
  TH1F *t1t2_PU0_Ratio;
  TH1F *t2t3_PU0_Ratio;
  TH1F *t3t4_PU0_Ratio;

  TH1F *nSubj1_PULo_hist;
  TH1F *nSubj2_PULo_hist;
  TH1F *nSubj3_PULo_hist;
  TH1F *nSubj4_PULo_hist;
  TH1F *t3t1_PULo_Ratio;
  TH1F *t2t1_PULo_Ratio;
  TH1F *t3t2_PULo_Ratio;
  TH1F *t4t1_PULo_Ratio;
  TH1F *t1t2_PULo_Ratio;
  TH1F *t2t3_PULo_Ratio;
  TH1F *t3t4_PULo_Ratio;

  TH1F *nSubj1_PUHi_hist;
  TH1F *nSubj2_PUHi_hist;
  TH1F *nSubj3_PUHi_hist;
  TH1F *nSubj4_PUHi_hist;
  TH1F *t3t1_PUHi_Ratio;
  TH1F *t2t1_PUHi_Ratio;
  TH1F *t3t2_PUHi_Ratio;
  TH1F *t4t1_PUHi_Ratio;
  TH1F *t1t2_PUHi_Ratio;
  TH1F *t2t3_PUHi_Ratio;
  TH1F *t3t4_PUHi_Ratio;

  unsigned int momPDGID_;
  double genMuTauPTMin_;
  double genMuPTMin_;
  edm::ParameterSet* cfg_;
  edm::InputTag jetSrc_;
  std::string Sample_;
  std::string outFileName_;
  std::string outputTextFile_;
  int nFilt_;
  double rFilt_;
  double trimPtFracMin_;
  double zCut_;
  double RcutFactor_;

  ofstream SurvivingAxes;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FastJetAnalyzer_alltaus::FastJetAnalyzer_alltaus(const edm::ParameterSet& iConfig):out_(0), jetTree(0)
{
  //now do what ever initialization is needed
  Sample_ = iConfig.getParameter<std::string>("Sample");
  jetSrc_ = iConfig.getParameter<edm::InputTag>("jetSrc");
  outFileName_ = iConfig.getParameter<std::string>("outFileName");
  outputTextFile_ = iConfig.getParameter<std::string>("outputTextFile");
  nFilt_ = iConfig.getParameter<int>("nFilt");
  rFilt_ = iConfig.getParameter<double>("rFilt");
  trimPtFracMin_ = iConfig.getParameter<double>("trimPtFracMin");
  zCut_ = iConfig.getParameter<double>("zCut");
  RcutFactor_ = iConfig.getParameter<double>("RcutFactor");
  cfg_ = const_cast<edm::ParameterSet*>(&iConfig);
}


FastJetAnalyzer_alltaus::~FastJetAnalyzer_alltaus()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//


// ------------ method called for each event  ------------
void
FastJetAnalyzer_alltaus::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ // start analyze
  cout.precision(10);
  //// fastjet filter stuff
  
  fastjet::Filter trimmer( fastjet::Filter(fastjet::JetDefinition(fastjet::antikt_algorithm, rFilt_), fastjet::SelectorPtFractionMin(trimPtFracMin_)));
  fastjet::Filter filter( fastjet::Filter(fastjet::JetDefinition(fastjet::antikt_algorithm, rFilt_), fastjet::SelectorNHardest(nFilt_)));
  fastjet::Pruner pruner(fastjet::antikt_algorithm, zCut_, RcutFactor_);

  ////// pileup stuff /////////////

  edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
  std::vector<PileupSummaryInfo>::const_iterator PVI;

  float npT=-1.;
  float npIT=-1.;
  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
    int BX = PVI->getBunchCrossing();
    if(BX == 0) { // "0" is the in-time crossing, negative values are the early crossings, positive are late
      npT = PVI->getTrueNumInteractions();
      npIT = PVI->getPU_NumInteractions();
      intime_PU_vertices->Fill((int)npIT);
    }
  }
  
  float jetPt, jetMass;
  //float minMass, topMass;
  int nSubjets;
  float nSubj1, nSubj2, nSubj3, nSubj4;
  
  // Get gen particle collection
   Handle<GenParticleCollection> genParticles;
   iEvent.getByLabel("genParticles", genParticles);


   // Get jets
  Handle< std::vector<reco::PFJet> > jetH;
  iEvent.getByLabel(jetSrc_, jetH);
  vector<reco::PFJet> const * jets = jetH.product();
  

  if (Sample_ == "Wh1")
    {
      // Collect hadronically decaying gen taus
      // with muonically decaying sisters
      // from a decays

      std::vector<unsigned int> keysToIgnore;
      std::vector<GenTauDecayID> tauDecays; // vector of taus from a decays
      std::vector<reco::GenParticle*> genTausFromA;
      for (GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen)
	{ // loop over genparticles to get tau decays
	  
	  if((fabs(iGen->pdgId()) == 15) && (iGen->status() == 3))
	    { // if it is a status-3 tau
	      if((fabs(iGen->mother()->pdgId()) == 36))
		{ // if its mother was an a
		  //cout << "A status 3 tau from an a has been found" << endl;
		  //cout << "status of a = " << iGen->mother()->status() << endl;
		  genTausFromA.push_back(const_cast<reco::GenParticle*>(&*iGen));
		  try{
		    //cout << "Seeking tau decay" << endl;
		    GenTauDecayID tauDecay(*cfg_, genParticles, iGen-genParticles->begin());
		    //cout << "got tau decay" << endl;
		    tauDecays.push_back(tauDecay);
		    //cout << "Found tau decay" << endl;
		  }
		  catch (char * str) { cout << "Tau Not Found" << str << '\n'; }
		} // if its mother was an a
	    } // if it is a status-3 tau
	  
	} // loop over genparticles to get tau decays

      for (unsigned i = 0; i < jets->size(); i++)
	{ // loop over pf jets
	  if (jets->at(i).pt() > 30.)
	    { // if pfjet pt > 30
	      if (genTausFromA.size() != 0)
		{ // if there are gen taus from a
		  double mindelR = 100.;
		  double delR = 100.;
		  unsigned int tauIterator = 0;
		  unsigned int tauPointer = 0;
		  for (unsigned int j = 0; j < genTausFromA.size(); j++)
		    { // loop over gen taus from a
		      
		      delR = deltaR(jets->at(i).eta(), jets->at(i).phi(), genTausFromA.at(j)->eta(), genTausFromA.at(j)->phi());
		      if (delR < mindelR)
			{
			  mindelR = delR;
			  tauPointer = j;
			}
		      tauIterator += 1;
		    } // loop over gen taus from a
		  if (mindelR < 0.3)
		    { // if mindelR < 0.3, jet is matched to a gen tau from an a
		      double delR_sisters = 9999.;
		      int decayModeCode = 6;
		      bool mode_ee = false;
		      bool mode_em = false;
		      bool mode_mm = false;
		      bool mode_eh = false;
		      bool mode_mh = false;
		      bool mode_hh = false;
		      for (unsigned k = 0; k != genTausFromA.size(); k++)
			{ // loop over taus to find sister
			  if (k == tauPointer)
			    continue;
			  else
			    { // if k != taupointer
			      if (genTausFromA.at(tauPointer)->motherRef().key() == genTausFromA.at(k)->motherRef().key())
				{ // calculate `delR and print decay modes
				  delR_sisters = deltaR(genTausFromA.at(tauPointer)->eta(),genTausFromA.at(tauPointer)->phi(),genTausFromA.at(k)->eta(),genTausFromA.at(k)->phi());
				  
				  // decayModeCode dictionary //
				  // 0 = ee; 1 = mumu; 2 = emu; //
				  // 3 = ehad; 4 = muhad; 5 = hadhad; //
				  // 6 = none of the above; //
				  try
				    {//try
				      tauDecays.at(tauPointer).findSister();
				      const unsigned int iSister = tauDecays.at(tauPointer).getSisterIndex();
				      //cout << "Found sister? " << tauDecays.at(tauPointer).foundSister() << endl;
				      if (tauDecays.at(tauPointer).foundSister() == true) {
				      bool tauWasHad = false;
				      bool sisterWasHad = false;
				      int tauHadMode = 99;
				      int sistertauHadMode = 99;
				      int tauMode = tauDecays.at(tauPointer).tauDecayType();
				      int sistertauMode = tauDecays.at(tauPointer).sisterDecayType();
				      if (tauMode == GenTauDecayID::HAD)
					{
					  tauHadMode = tauDecays.at(tauPointer).tauHadronicDecayType();
					  tauWasHad = true;
					}
				      if (sistertauMode == GenTauDecayID::HAD)
					{
					  sistertauHadMode = tauDecays.at(tauPointer).sisterHadronicDecayType();
					  sisterWasHad = true;
					}
				      if (tauWasHad == false && sisterWasHad == false)
					{ // fully leptonic
					  if (tauMode == GenTauDecayID::MU)
					    {
					      if (sistertauMode == GenTauDecayID::MU)
						{ decayModeCode = 1; mode_mm = true; }
					      else if (sistertauMode == GenTauDecayID::E)
						{ decayModeCode = 2; mode_em = true; }
					    }
					  else if (tauMode == GenTauDecayID::E)
					    {
					      if (sistertauMode == GenTauDecayID::E)
						{ decayModeCode = 0; mode_ee = true; }
					      else if (sistertauMode == GenTauDecayID::MU)
						{ decayModeCode = 2; mode_em = true; }
					    }
					  continue;
					} // fully leptonic
				      else if (tauWasHad == true && sisterWasHad == false)
					{ // HadLep
					  if (sistertauMode == GenTauDecayID::E)
					    {
					      { decayModeCode = 3; mode_eh = true; }
					    }
					  else if (sistertauMode == GenTauDecayID::MU)
					    {
					      { decayModeCode = 4; mode_mh = true; }
					    }
					  continue;
					} // HadLep
				      else if (tauWasHad == false && sisterWasHad == true)
					{ // LepHad
					  if (tauMode == GenTauDecayID::E)
					    {
					      { decayModeCode = 3; mode_eh = true; }
					    }
					  else if (tauMode == GenTauDecayID::MU)
					    {
			        	      { decayModeCode = 4; mode_mh = true; }
					    }
					  continue;
					} // LepHad
				      else if (tauWasHad == true && sisterWasHad == true)
					{ // fully hadronic
					  decayModeCode = 5;
					  mode_hh = true;
					  /*
					  if (tauHadMode == 0 && sistertauHadMode == 0)
					    decayModeCode = 0; // (1P,1P)
					  else if ((tauHadMode == 0 && sistertauHadMode == 1) || (tauHadMode == 1 && sistertauHadMode == 0))
					    decayModeCode = 1; // (1P,1PS)
					  else if ((tauHadMode == 0 && sistertauHadMode == 2) || (tauHadMode == 2 && sistertauHadMode == 0))
					    decayModeCode = 2; // (1P,1P2S)
					  else if ((tauHadMode == 0 && sistertauHadMode == 10) || (tauHadMode == 10 && sistertauHadMode == 0))
					    decayModeCode = 3; // (1P,3P)
					  else if ((tauHadMode == 1 && sistertauHadMode == 1))
					    decayModeCode = 4; // (1PS,1PS)
					  else if ((tauHadMode == 1 && sistertauHadMode == 2) || (tauHadMode == 2 && sistertauHadMode == 1))
					    decayModeCode = 5; // (1PS,1P2S)
					  else if ((tauHadMode == 1 && sistertauHadMode == 10) || (tauHadMode == 10 && sistertauHadMode == 1))
					    decayModeCode = 6; // (1PS,3P)
					  else if ((tauHadMode == 2 && sistertauHadMode == 2))
					    decayModeCode = 7; // (1P2S,1P2S)
					  else if ((tauHadMode == 2 && sistertauHadMode == 10) || (tauHadMode == 10 && sistertauHadMode == 2))
					    decayModeCode = 8; // (1P2S,3P)
					  else if ((tauHadMode == 10 && sistertauHadMode == 10))
					    decayModeCode = 9; // (3P,3P)
					  */
					  //decayModeCode = tauHadMode;
					} // fully hadronic
				      else
					{
					  decayModeCode = 6;
					  //cout << "The combination was: " << tauHadMode << ", " << sistertauHadMode << endl;
					}
				      }
				    }//try
				  catch (char * str) { cout << "Decay mode not found" << str << '\n'; }
				} // calculate delR and print decay modes
			    } // if k != taupointer
			  
			} // loop over taus to find sister
		      //}
		    
		      genTausFromA.erase(genTausFromA.begin()+tauPointer);
		      tauDecays.erase(tauDecays.begin()+tauPointer);

		      jetPt = jets->at(i).pt();
		      jetPt_hist->Fill(jets->at(i).pt());
		      jetMass = jets->at(i).mass();
		      jetMass_hist->Fill(jets->at(i).mass());
		      nSubjets = jets->at(i).numberOfDaughters();
		      
		      pjmj_PUAll->Fill(jetPt/jetMass);
		      if (npIT == 0 && jetMass != 0)
			pjmj_PU0->Fill(jetPt/jetMass);
		      
		      float d_0 = 0;
		      float NsubVals[10] = { 0.0 } ;
		      float minDR = 999.;
		      float dPhi, dEta, dR;
		      
		      vector<reco::PFCandidatePtr> pfCands = jets->at(i).getPFConstituents();
		      vector<const reco::Candidate*> subjets = jets->at(i).getJetConstituentsQuick();
		      vector<const reco::PFCandidate*> all_particles;
		      for (unsigned j = 0; j < pfCands.size(); j++){
			const reco::PFCandidate *thisPF = pfCands.at(j).get(); 
			all_particles.push_back( thisPF );	
		      }
		      
		      vector<fastjet::PseudoJet> FJparticles;
		      
		      for (unsigned particle = 0; particle < all_particles.size(); particle++){
			const reco::PFCandidate *thisParticle = all_particles.at(particle);
			FJparticles.push_back( fastjet::PseudoJet( thisParticle->px(), thisParticle->py(), thisParticle->pz(), thisParticle->energy() ) );
		      }
		      //cout << "# of constituents before grooming: " << FJparticles.size() << endl;
		      cands_beforeGrooming->Fill(FJparticles.size());
		      // use groomed
		      fastjet::JetDefinition jet_def(fastjet::kt_algorithm, 0.5);
		      fastjet::ClusterSequence thisClustering(FJparticles, jet_def);
		      vector<fastjet::PseudoJet> FJjet = thisClustering.inclusive_jets();
		      fastjet::PseudoJet thisMainJet = FJjet.at(0);
		      fastjet::PseudoJet thisGroomedJet = pruner(thisMainJet);
		      vector<fastjet::PseudoJet> FJparticles2 = thisGroomedJet.constituents();
		      double groomedjetpt = thisGroomedJet.pt();
		      pT_Groomed->Fill(groomedjetpt);
		      //cout << "# of constituents after grooming: " << FJparticles2.size() << endl;
		      cands_afterGrooming->Fill(FJparticles2.size());
		      
		      //vector<fastjet::PseudoJet> FJparticles2 = thisClustering.constituents(thisGroomedJet); // thisMainJet changed to thisGroomedJet
		      
		      
		      int index = 0;
		 
		      NsubParameters paraNsub = NsubParameters(1.0, 0.8);
		      
		      //Nsubjettiness routine(1, Njettiness::antikt_0p2_axes, paraNsub);
		      Nsubjettiness routine1(1, Njettiness::kt_axes, 1.0, 0.8, 10000.0);

		      //double tau1 = routine.getTau(FJparticles2); // FJparticles -> FJparticles 2 
		      double tau1 = routine1.result(thisGroomedJet);
		      Nsubjettiness routine2(2, Njettiness::kt_axes, 1.0, 0.8, 10000.0);
		      Nsubjettiness routine3(3, Njettiness::kt_axes, 1.0, 0.8, 10000.0);
		      Nsubjettiness routine4(4, Njettiness::kt_axes, 1.0, 0.8, 10000.0);

		      double tau2 = routine2.result(thisGroomedJet); // ""
		      double tau3 = routine3.result(thisGroomedJet); // ""
		      double tau4 = routine4.result(thisGroomedJet); // ""
		      //std::vector<fastjet::PseudoJet> remainingSubaxes = routine.currentAxes();

		      //std::vector<fastjet::PseudoJet> remainingSubaxes = routine.getJets(FJparticles2);

		      // use groomed
		      
		      if(tau3/tau1 == 1.)
			{
			  cout << "# of constituents before grooming: " << FJparticles.size() << endl;
			  cout << "# of constituents after grooming: " << FJparticles2.size() << endl;
			  cout << "tau1 = " << tau1 << endl;
			  cout << "tau3 = " << tau3 << endl;
			  if (delR_sisters != 9999.)
			    delR_sisters_One->Fill(delR_sisters);
			}
		      else
			{
			  if (delR_sisters != 9999.)
			    delR_sisters_NotOne->Fill(delR_sisters);
			}

		      if(tau3/tau1 == 1.)
			{ 
			  decayModes_One->Fill(decayModeCode);
			}
		      else
			{ 
			  decayModes_NotOne->Fill(decayModeCode);
			}
			  		      
		      if(tau1 != 0.)
			{
			  nSubj_vs_pjmj_PUAll->Fill(jetPt/jetMass, tau3/tau1);
			  t3t1_vs_t4t1->Fill(tau4/tau1, tau3/tau1);
			  t3t1_vs_t2t1->Fill(tau2/tau1, tau3/tau1);
			  t2t1_vs_t4t1->Fill(tau4/tau1, tau2/tau1);
			  t3t1_PUAll_Ratio->Fill(tau3/tau1);
			  if (mode_ee == true)
			    t3t1_ee->Fill(tau3/tau1);
			  if (mode_mm == true)
			    t3t1_mm->Fill(tau3/tau1);
			  if (mode_em == true)
			    t3t1_em->Fill(tau3/tau1);
			  if (mode_eh == true)
			    t3t1_eh->Fill(tau3/tau1);
			  if (mode_mh == true)
			    t3t1_mh->Fill(tau3/tau1);
			  if (mode_hh == true)
			    t3t1_hh->Fill(tau3/tau1);

			  if (tau3/tau1 == 1.)
			    {
			      tau1_One->Fill(tau1);
			      npIT_One->Fill(npIT);
			      cands_One->Fill(FJparticles2.size());
			      pT_One->Fill(groomedjetpt);
			    }
			  else
			    {
			      tau1_NotOne->Fill(tau1);
			      npIT_NotOne->Fill(npIT);
			      cands_NotOne->Fill(FJparticles2.size());
			      pT_NotOne->Fill(groomedjetpt);
			    }
			  
			  // now look at the jet constituent pT
			  double pmax = -100000.;
			  double p = 0.;
			  unsigned highestpointer = 0;
			  for(unsigned particle = 0; particle != FJparticles2.size(); particle++)
			    { // find highest pT
			      p = FJparticles2.at(particle).pt();
			      if (p > pmax)
				{ pmax = p;
				  highestpointer = particle;
				}
			    } // find highest pT
			  double psecondmax = -100000.;
			  double maxDelR = -1.;
			  double etahigh = FJparticles2.at(highestpointer).eta();
			  double phihigh = FJparticles2.at(highestpointer).phi();
			  double etapart, phipart;
			  for(unsigned particle = 0; particle != FJparticles2.size(); particle++)
			    { // find second highest pT
			      p = FJparticles2.at(particle).pt();
			      etapart = FJparticles2.at(particle).eta();
			      phipart = FJparticles2.at(particle).phi();
			      if (p > psecondmax && p < pmax)
				{ psecondmax = p; }
			      double delR = deltaR(etahigh, phihigh, etapart, phipart);
			      if (delR > maxDelR)
				maxDelR = delR;	     
			    } // find second highest pT
			  if (tau3/tau1 == 1.)
			    {
			      pT_highest_One->Fill(pmax);
			      pT_secondhighest_One->Fill(psecondmax);
			      maxDelR_One->Fill(maxDelR);
			      cout << "Highest pT = " << pmax << endl;
			      cout << "Second-highest pT = " << psecondmax << endl;
			      cout << "Maximum delR from highest pT cand = " << maxDelR << endl;
			    }
			  else
			    {
			      pT_highest_NotOne->Fill(pmax);
			      pT_secondhighest_NotOne->Fill(psecondmax);
			      maxDelR_NotOne->Fill(maxDelR);
			    }
			  
			}
		      if (npIT == 0)
			{
			  nSubj1_PU0_hist->Fill(tau1);
			  nSubj2_PU0_hist->Fill(tau2);
			  nSubj3_PU0_hist->Fill(tau3);
			  nSubj4_PU0_hist->Fill(tau4);
			  if(tau1 != 0.)
			    {
			      t3t1_PU0_Ratio->Fill(tau3/tau1);
			      t2t1_PU0_Ratio->Fill(tau2/tau1);
			      nSubj_vs_pjmj_PU0->Fill(jetPt/jetMass, tau3/tau1);
			    }
			  if(tau2 != 0.)
			    {
			      t3t2_PU0_Ratio->Fill(tau3/tau2);
			      t1t2_PU0_Ratio->Fill(tau1/tau2);
			    }
			  if(tau1 != 0.)
			    t4t1_PU0_Ratio->Fill(tau4/tau1);
			  if(tau3 != 0.)
			    t2t3_PU0_Ratio->Fill(tau2/tau3);
			  if(tau4 != 0.)
			    t3t4_PU0_Ratio->Fill(tau3/tau4);
			}
		      else if (npIT < 20)
			{
			  nSubj1_PULo_hist->Fill(tau1);
			  nSubj2_PULo_hist->Fill(tau2);
			  nSubj3_PULo_hist->Fill(tau3);
			  nSubj4_PULo_hist->Fill(tau4);
			  if(tau1 != 0.)
			    {
			      t3t1_PULo_Ratio->Fill(tau3/tau1);
			      t2t1_PULo_Ratio->Fill(tau2/tau1);
			      nSubj_vs_pjmj_PULo->Fill(jetPt/jetMass, tau3/tau1);
			    }
			  if(tau2 != 0.)
			    {
			      t3t2_PULo_Ratio->Fill(tau3/tau2);
			      t1t2_PULo_Ratio->Fill(tau1/tau2);
			    }
			  if(tau1 != 0.)
			    t4t1_PULo_Ratio->Fill(tau4/tau1);
			  if(tau3 != 0.)
			    t2t3_PULo_Ratio->Fill(tau2/tau3);
			  if(tau4 != 0.)
			    t3t4_PULo_Ratio->Fill(tau3/tau4);
			}
		      else
			{
			  nSubj1_PUHi_hist->Fill(tau1);
			  nSubj2_PUHi_hist->Fill(tau2);
			  nSubj3_PUHi_hist->Fill(tau3);
			  nSubj4_PUHi_hist->Fill(tau4);
			  if(tau1 != 0.)
			    {
			      t3t1_PUHi_Ratio->Fill(tau3/tau1);
			      t2t1_PUHi_Ratio->Fill(tau2/tau1);
			      nSubj_vs_pjmj_PUHi->Fill(jetPt/jetMass, tau3/tau1);
			    }
			  if(tau2 != 0.)
			    {
			      t3t2_PUHi_Ratio->Fill(tau3/tau2);
			      t1t2_PUHi_Ratio->Fill(tau1/tau2);
			    }
			  if(tau1 != 0.)
			    t4t1_PUHi_Ratio->Fill(tau4/tau1);
			  if(tau3 != 0.)
			    t2t3_PUHi_Ratio->Fill(tau2/tau3);
			  if(tau4 != 0.)
			    t3t4_PUHi_Ratio->Fill(tau3/tau4);
			}
		      
		      nSubj1 = tau1; //NsubVals[0];	
		      nSubj2 = tau2;//NsubVals[1];	
		      nSubj3 = tau3;//NsubVals[2];	
		      nSubj4 = tau4;//NsubVals[3];	
		      jetTree->Fill();
		      
		    } // if jet is matched to a gen tau from an a
		  
		} // if there are gen taus from a
	    } // if pfjet pt > 30
	} // loop over pf jets
      
    } // if Sample = Wh1
  
  if (Sample_ == "WJets")
    { // if Sample = WJets
      
      std::vector<reco::GenParticle*> genMuonsFromW;
      // Loop over genparticles to find muons from W->mumu
      for (GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen)
	{ // loop over genparticles
	  if((fabs(iGen->pdgId()) == 13) && (fabs(iGen->mother()->pdgId()) != 13))
	    { // if it's a muon whose mother is not a muon
	      if((fabs(iGen->mother()->pdgId()) == 24))
		{ // if its mother is a W
		  genMuonsFromW.push_back(const_cast<reco::GenParticle*>(&*iGen));
		} // if its mother is a W
	    } // if it's a muon whose mother is not a muon
	} // loop over genparticles
      
      
      for (unsigned i = 0; i < jets->size(); i++)
	{ // loop over pf jets
	  
	  std::vector<reco::PFCandidatePtr> JetPFCands = jets->at(i).getPFConstituents();
	  
	  bool matchedToMuonFromW = false;
	  double mindelR = 100.;
	  double delR = 100;
	  unsigned int muonIterator = 0;
	  unsigned int muonPointer = 0;
	  
	  for(std::vector<edm::Ptr<reco::PFCandidate> >::iterator it = JetPFCands.begin(); it != JetPFCands.end(); ++it)
	    { // loop over PF candidates
	      if (matchedToMuonFromW == true)
		break;
	      reco::PFCandidate pfCand = *it;
	      if (pfCand.particleId() == 3)
		{// if it's a PF muon
		  reco::MuonRef theRecoMuon = pfCand.muonRef();
		  muonPointer = 0;
		  //cout << genMuonsFromW.size() << endl;
		  for (unsigned int j = 0; j < genMuonsFromW.size(); j++)
		    { // matching to gen W->mumu
		      delR = deltaR(theRecoMuon->eta(), theRecoMuon->phi(),genMuonsFromW.at(j)->eta(), genMuonsFromW.at(j)->phi());
		      if (delR < mindelR)
			{
			  mindelR = delR;
			  muonPointer = j;
			}
		      muonIterator += 1;
		    } // matching to gen W->mumu
		  if (mindelR < 0.3)
		    {
		      matchedToMuonFromW = true;
		      genMuonsFromW.erase(genMuonsFromW.begin()+muonPointer);
		    }
		}// if it's a PF muon 
	    } // loop over PF candidates
	  
	  if((matchedToMuonFromW == false) && (jets->at(i).pt() > 30.))
	    {//if(matchedToMuonFromW == false)
	      
	      jetPt = jets->at(i).pt();
	      jetPt_hist->Fill(jets->at(i).pt());
	      jetMass = jets->at(i).mass();
	      jetMass_hist->Fill(jets->at(i).mass());
	      nSubjets = jets->at(i).numberOfDaughters();
	      
	      pjmj_PUAll->Fill(jetPt/jetMass);
	      if (npIT == 0 && jetMass != 0)
		pjmj_PU0->Fill(jetPt/jetMass);
	      
	      float d_0 = 0;
	      float NsubVals[10] = { 0.0 } ;
	      float minDR = 999.;
	      float dPhi, dEta, dR;
	      
	      vector<reco::PFCandidatePtr> pfCands = jets->at(i).getPFConstituents();
	      vector<const reco::Candidate*> subjets = jets->at(i).getJetConstituentsQuick();
	      vector<const reco::PFCandidate*> all_particles;
	      for (unsigned j = 0; j < pfCands.size(); j++){
		const reco::PFCandidate *thisPF = pfCands.at(j).get(); 
		all_particles.push_back( thisPF );	
	      }
	      
	      vector<fastjet::PseudoJet> FJparticles;
	      
	      for (unsigned particle = 0; particle < all_particles.size(); particle++){
		const reco::PFCandidate *thisParticle = all_particles.at(particle);
		FJparticles.push_back( fastjet::PseudoJet( thisParticle->px(), thisParticle->py(), thisParticle->pz(), thisParticle->energy() ) );
	      }
	      //cout << "# of constituents before grooming: " << FJparticles.size() << endl;
	      cands_beforeGrooming->Fill(FJparticles.size());
	      // use groomed
	      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.5);
	      fastjet::ClusterSequence thisClustering(FJparticles, jet_def);
	      vector<fastjet::PseudoJet> FJjet = thisClustering.inclusive_jets();
	      fastjet::PseudoJet thisMainJet = FJjet.at(0);
	      fastjet::PseudoJet thisGroomedJet = pruner(thisMainJet);
	      vector<fastjet::PseudoJet> FJparticles2 = thisGroomedJet.constituents();
	      double groomedjetpt = thisGroomedJet.pt();
	      pT_Groomed->Fill(groomedjetpt);
	      //cout << "# of constituents after grooming: " << FJparticles2.size() << endl;
	      cands_afterGrooming->Fill(FJparticles2.size());
	      
	      //vector<fastjet::PseudoJet> FJparticles2 = thisClustering.constituents(thisGroomedJet); // thisMainJet changed to thisGroomedJet
	      
	      
	      int index = 0;
	      
	      NsubParameters paraNsub = NsubParameters(1.0, 0.8);
	      
	      //Nsubjettiness routine(1, Njettiness::antikt_0p2_axes, paraNsub);
	      Nsubjettiness routine1(1, Njettiness::kt_axes, 1.0, 0.8, 10000.0);
	      
	      //double tau1 = routine.getTau(FJparticles2); // FJparticles -> FJparticles 2 
	      double tau1 = routine1.result(thisGroomedJet);
	      Nsubjettiness routine2(2, Njettiness::kt_axes, 1.0, 0.8, 10000.0);
	      Nsubjettiness routine3(3, Njettiness::kt_axes, 1.0, 0.8, 10000.0);
	      Nsubjettiness routine4(4, Njettiness::kt_axes, 1.0, 0.8, 10000.0);
	      
	      double tau2 = routine2.result(thisGroomedJet); // ""
	      double tau3 = routine3.result(thisGroomedJet); // ""
	      double tau4 = routine4.result(thisGroomedJet); // ""
	      
	      // use groomed
	      
	      if(tau3/tau1 == 1.)
		{
		  cout << "# of constituents before grooming: " << FJparticles.size() << endl;
		  cout << "# of constituents after grooming: " << FJparticles2.size() << endl;
		  cout << "tau1 = " << tau1 << endl;
		  cout << "tau3 = " << tau3 << endl;
		}
	      	      
	      if(tau1 != 0.)
		{
		  nSubj_vs_pjmj_PUAll->Fill(jetPt/jetMass, tau3/tau1);
		  t3t1_vs_t4t1->Fill(tau4/tau1, tau3/tau1);
		  t3t1_vs_t2t1->Fill(tau2/tau1, tau3/tau1);
		  t2t1_vs_t4t1->Fill(tau4/tau1, tau2/tau1);
		  t3t1_PUAll_Ratio->Fill(tau3/tau1);

		  if (tau3/tau1 == 1.)
		    {
		      tau1_One->Fill(tau1);
		      npIT_One->Fill(npIT);
		      cands_One->Fill(FJparticles2.size());
		      pT_One->Fill(groomedjetpt);
		    }
		  else
		    {
		      tau1_NotOne->Fill(tau1);
		      npIT_NotOne->Fill(npIT);
		      cands_NotOne->Fill(FJparticles2.size());
		      pT_NotOne->Fill(groomedjetpt);
		    }
		  
		  // now look at the jet constituent pT
		  double pmax = -100000.;
		  double p = 0.;
		  unsigned highestpointer = 0;
		  for(unsigned particle = 0; particle != FJparticles2.size(); particle++)
		    { // find highest pT
		      p = FJparticles2.at(particle).pt();
		      if (p > pmax)
			{ pmax = p;
			  highestpointer = particle;
			}
		    } // find highest pT
		  double psecondmax = -100000.;
		  double maxDelR = -1.;
		  double etahigh = FJparticles2.at(highestpointer).eta();
		  double phihigh = FJparticles2.at(highestpointer).phi();
		  double etapart, phipart;
		  for(unsigned particle = 0; particle != FJparticles2.size(); particle++)
		    { // find second highest pT
		      p = FJparticles2.at(particle).pt();
		      etapart = FJparticles2.at(particle).eta();
		      phipart = FJparticles2.at(particle).phi();
		      if (p > psecondmax && p < pmax)
			{ psecondmax = p; }
		      double delR = deltaR(etahigh, phihigh, etapart, phipart);
		      if (delR > maxDelR)
			maxDelR = delR;	     
		    } // find second highest pT
		  if (tau3/tau1 == 1.)
		    {
		      pT_highest_One->Fill(pmax);
		      pT_secondhighest_One->Fill(psecondmax);
		      maxDelR_One->Fill(maxDelR);
		      cout << "Highest pT = " << pmax << endl;
		      cout << "Second-highest pT = " << psecondmax << endl;
		      cout << "Maximum delR from highest pT cand = " << maxDelR << endl;
		    }
		  else
		    {
		      pT_highest_NotOne->Fill(pmax);
		      pT_secondhighest_NotOne->Fill(psecondmax);
		      maxDelR_NotOne->Fill(maxDelR);
		    }
		  
			}
	      if (npIT == 0)
		{
		  nSubj1_PU0_hist->Fill(tau1);
		  nSubj2_PU0_hist->Fill(tau2);
		  nSubj3_PU0_hist->Fill(tau3);
		  nSubj4_PU0_hist->Fill(tau4);
		  if(tau1 != 0.)
		    {
		      t3t1_PU0_Ratio->Fill(tau3/tau1);
		      t2t1_PU0_Ratio->Fill(tau2/tau1);
		      nSubj_vs_pjmj_PU0->Fill(jetPt/jetMass, tau3/tau1);
		    }
		  if(tau2 != 0.)
		    {
		      t3t2_PU0_Ratio->Fill(tau3/tau2);
		      t1t2_PU0_Ratio->Fill(tau1/tau2);
		    }
		  if(tau1 != 0.)
		    t4t1_PU0_Ratio->Fill(tau4/tau1);
		  if(tau3 != 0.)
		    t2t3_PU0_Ratio->Fill(tau2/tau3);
		  if(tau4 != 0.)
		    t3t4_PU0_Ratio->Fill(tau3/tau4);
		}
	      else if (npIT < 20)
		{
		  nSubj1_PULo_hist->Fill(tau1);
		  nSubj2_PULo_hist->Fill(tau2);
		  nSubj3_PULo_hist->Fill(tau3);
		  nSubj4_PULo_hist->Fill(tau4);
		  if(tau1 != 0.)
		    {
		      t3t1_PULo_Ratio->Fill(tau3/tau1);
		      t2t1_PULo_Ratio->Fill(tau2/tau1);
		      nSubj_vs_pjmj_PULo->Fill(jetPt/jetMass, tau3/tau1);
		    }
		  if(tau2 != 0.)
		    {
		      t3t2_PULo_Ratio->Fill(tau3/tau2);
		      t1t2_PULo_Ratio->Fill(tau1/tau2);
		    }
		  if(tau1 != 0.)
		    t4t1_PULo_Ratio->Fill(tau4/tau1);
		  if(tau3 != 0.)
		    t2t3_PULo_Ratio->Fill(tau2/tau3);
		  if(tau4 != 0.)
		    t3t4_PULo_Ratio->Fill(tau3/tau4);
		}
	      else
		{
		  nSubj1_PUHi_hist->Fill(tau1);
		  nSubj2_PUHi_hist->Fill(tau2);
		  nSubj3_PUHi_hist->Fill(tau3);
		  nSubj4_PUHi_hist->Fill(tau4);
		  if(tau1 != 0.)
		    {
		      t3t1_PUHi_Ratio->Fill(tau3/tau1);
		      t2t1_PUHi_Ratio->Fill(tau2/tau1);
		      nSubj_vs_pjmj_PUHi->Fill(jetPt/jetMass, tau3/tau1);
		    }
		  if(tau2 != 0.)
		    {
		      t3t2_PUHi_Ratio->Fill(tau3/tau2);
		      t1t2_PUHi_Ratio->Fill(tau1/tau2);
		    }
		  if(tau1 != 0.)
		    t4t1_PUHi_Ratio->Fill(tau4/tau1);
		  if(tau3 != 0.)
		    t2t3_PUHi_Ratio->Fill(tau2/tau3);
		  if(tau4 != 0.)
		    t3t4_PUHi_Ratio->Fill(tau3/tau4);
		}
	      
	      nSubj1 = tau1; //NsubVals[0];	
	      nSubj2 = tau2;//NsubVals[1];	
	      nSubj3 = tau3;//NsubVals[2];	
	      nSubj4 = tau4;//NsubVals[3];	
	      jetTree->Fill();
	      
	    } // if jet is not matched to a mu from a W
	} // loop over pf jets
      
    } // if Sample = WJets

} // end analyze


// ------------ method called once each job just before starting event loop  ------------
void 
FastJetAnalyzer_alltaus::beginJob()
{
  //  out_ = new TFile("/data1/friccita/NSJdatasets_10242012/ak5_groomed/FastJetAnalysis_signalWH_alltaus_pileup_pruned_z0p05_D0p3.root", "RECREATE");
  out_ = new TFile(outFileName_.c_str(), "RECREATE");
  SurvivingAxes.open(outputTextFile_.c_str());
  //TTree *jetTree = new TTree("jets", "jets");
  jetTree = new TTree("jets", "jets");
  TBranch *jetPt = jetTree->Branch("jetPt", &jetPt, "jetPt/F");
  TBranch *jetMass = jetTree->Branch("jetMass", &jetMass, "jetMass/F");
  //TBranch *minMass = jetTree->Branch("minMass", &minMass, "minMass/F");
  // TBranch *topMass = jetTree->Branch("topMass", &topMass, "topMass/F");
  TBranch *nSubjets = jetTree->Branch("nSubjets", &nSubjets, "nSubjets/I");
  TBranch *nSubj1 = jetTree->Branch("nSubj1",  &nSubj1, "nSubj1/F");
  TBranch *nSubj2 = jetTree->Branch("nSubj2", &nSubj2, "nSubj2/F");
  TBranch *nSubj3 = jetTree->Branch("nSubj3", &nSubj3, "nSubj3/F");
  TBranch *nSubj4 = jetTree->Branch("nSubj4", &nSubj4, "nSubj4/F");

  jetPt_hist = new TH1F("jetPt_hist", "pT of jets", 100., 0., 200.);
  jetMass_hist = new TH1F("jetMass_hist", "mass of jets", 100., 0., 200.);
  intime_PU_vertices = new TH1I("intime_PU_vertices", "Number of in-time pileup vertices", 100, 0, 100);
  cands_beforeGrooming = new TH1I("cands_beforeGrooming", "Number of jet constituents before grooming", 50, 0, 50);
  cands_afterGrooming = new TH1I("cands_afterGrooming", "Number of jet constituents after grooming", 50, 0, 50);

  tau1_One = new TH1F("tau1_One", "tau1 when tau3/tau1 = 1", 110., 0., 1.1);
  cands_One = new TH1I("cands_One", "jetcands when tau3/tau1 = 1", 50, 0, 50);
  npIT_One = new TH1I("npIT_One", "npIT when tau3/tau1 = 1", 100, 0, 100);
  pT_One = new TH1F("pT_One", "groomed jet pT when tau3/tau1 = 1", 100., 0., 200.);
  tau1_NotOne = new TH1F("tau1_NotOne", "tau1 when tau3/tau1 != 1", 110., 0., 1.1);
  cands_NotOne = new TH1I("cands_NotOne", "jetcands when tau3/tau1 != 1", 50, 0, 50);
  npIT_NotOne = new TH1I("npIT_NotOne", "npIT when tau3/tau1 1= 1", 100, 0, 100);
  pT_NotOne = new TH1F("pT_NotOne", "groomed jet pT when tau3/tau1 != 1", 100., 0., 200.);
  pT_Groomed = new TH1F("pT_Groomed", "groomed jet pT", 100., 0., 200.);
  pT_highest_One = new TH1F("pT_highest_One", "highest pT cand when tau3/tau1 = 1", 100., 0., 200.);
  pT_secondhighest_One = new TH1F("pT_secondhighest_One", "second highest pT cand when tau3/tau1 = 1", 100., 0., 200.);
  maxDelR_One = new TH1F("maxDelR_One", "max delR from highest pT cand when tau3/tau1 = 1", 100, 0., 2.);
  pT_highest_NotOne = new TH1F("pT_highest_NotOne", "highest pT cand when tau3/tau1 != 1", 100., 0., 200.);
  pT_secondhighest_NotOne = new TH1F("pT_secondhighest_NotOne", "second highest pT cand when tau3/tau1 != 1", 100., 0., 200.);
  maxDelR_NotOne = new TH1F("maxDelR_NotOne", "max delR from highest pT cand when tau3/tau1 != 1", 100, 0., 2.);
  delR_sisters_One = new TH1F("delR_sisters_One", "delR between sister taus from a decay", 100, 0., 2.);
  delR_sisters_NotOne = new TH1F("delR_sisters_NotOne", "delR between sister taus from a decay", 100, 0., 2.);
  decayModes_One = new TH1I("decayModes_One", "decayModes when tau3/tau1 = 1", 11, 0, 11);
  decayModes_NotOne = new TH1I("decayModes_NotOne", "decayModes when tau3/tau1 != 1", 11, 0, 11);

  pjmj_PU0 = new TH1F("pjmj_PU0", "pT(jet)/mass(jet) for PU0", 30, 0., 60.);
  pjmj_PUAll = new TH1F("pjmj_PUAll", "pT(jet)/mass(jet) for all PU", 30, 0., 60.);
  nSubj_vs_pjmj_PU0 = new TH2F("nSubj_vs_pjmj_PU0", "tau_3/tau_1 vs pT(jet)/mass(jet) for PU0", 60, 0., 60., 22, 0., 1.1);
  nSubj_vs_pjmj_PULo = new TH2F("nSubj_vs_pjmj_PULo", "tau_3/tau_1 vs pT(jet)/mass(jet) for PULo", 60, 0., 60., 22, 0., 1.1);
  nSubj_vs_pjmj_PUHi = new TH2F("nSubj_vs_pjmj_PUHi", "tau_3/tau_1 vs pT(jet)/mass(jet) for PUHi", 60, 0., 60., 22, 0., 1.1);
  nSubj_vs_pjmj_PUAll = new TH2F("nSubj_vs_pjmj_PUAll", "tau_3/tau_1 vs pT(jet)/mass(jet) for all PU", 60, 0., 60., 22, 0., 1.1);

  t3t1_vs_t4t1 = new TH2F("t3t1_vs_t4t1", "tau_3/tau_1 vs tau_4/tau_1 for all PU", 22, 0., 1.1, 22, 0., 1.1);
  t3t1_vs_t2t1 = new TH2F("t3t1_vs_t2t1", "tau_3/tau_1 vs tau_2/tau_1 for all PU", 22, 0., 1.1, 22, 0., 1.1);
  t2t1_vs_t4t1 = new TH2F("t2t1_vs_t4t1", "tau_2/tau_1 vs tau_4/tau_1 for all PU", 22, 0., 1.1, 22, 0., 1.1);

  t3t1_PUAll_Ratio = new TH1F("t3t1_PUAll_Ratio", "Ratio of tau_3 to tau_1 for all PU", 15, 0., 1.5);
  t3t1_ee = new TH1F("t3t1_ee", "tau_3/tau_1 for ee decay mode", 15, 0., 1.5);
  t3t1_mm = new TH1F("t3t1_mm", "tau_3/tau_1 for mumu decay mode", 15, 0., 1.5);
  t3t1_em = new TH1F("t3t1_em", "tau_3/tau_1 for emu decay mode", 15, 0., 1.5);
  t3t1_eh = new TH1F("t3t1_eh", "tau_3/tau_1 for ehad decay mode", 15, 0., 1.5);
  t3t1_mh = new TH1F("t3t1_mh", "tau_3/tau_1 for muhad decay mode", 15, 0., 1.5);
  t3t1_hh = new TH1F("t3t1_hh", "tau_3/tau_1 for hadhad decay mode", 15, 0., 1.5);

  // Zero pileup
  nSubj1_PU0_hist = new TH1F("nSubj1_PU0_hist", "1-subjettiness, zero pileup", 110., 0., 1.1);
  nSubj2_PU0_hist = new TH1F("nSubj2_PU0_hist", "2-subjettiness, zero pileup", 110., 0., 1.1);
  nSubj3_PU0_hist = new TH1F("nSubj3_PU0_hist", "3-subjettiness, zero pileup", 110., 0., 1.1);
  nSubj4_PU0_hist = new TH1F("nSubj4_PU0_hist", "4-subjettiness, zero pileup", 110., 0., 1.1);
  t3t1_PU0_Ratio = new TH1F("t3t1_PU0_Ratio", "Ratio of tau_3 to tau_1, zero pileup", 15., 0., 1.5);
  t2t1_PU0_Ratio = new TH1F("t2t1_PU0_Ratio", "Ratio of tau_2 to tau_1, zero pileup", 15., 0., 1.5);
  t3t2_PU0_Ratio = new TH1F("t3t2_PU0_Ratio", "Ratio of tau_3 to tau_2, zero pileup", 15., 0., 1.5);
  t4t1_PU0_Ratio = new TH1F("t4t1_PU0_Ratio", "Ratio of tau_4 to tau_1, zero pileup", 15., 0., 1.5);
  t1t2_PU0_Ratio = new TH1F("t1t2_PU0_Ratio", "Ratio of tau_1 to tau_2, zero pileup", 15., 0., 1.5);
  t2t3_PU0_Ratio = new TH1F("t2t3_PU0_Ratio", "Ratio of tau_2 to tau_3, zero pileup", 15., 0., 1.5);
  t3t4_PU0_Ratio = new TH1F("t3t4_PU0_Ratio", "Ratio of tau_3 to tau_4, zero pileup", 15., 0., 1.5);
  
  // Low pileup (< 20)
  nSubj1_PULo_hist = new TH1F("nSubj1_PULo_hist", "1-subjettiness, low pileup", 110., 0., 1.1);
  nSubj2_PULo_hist = new TH1F("nSubj2_PULo_hist", "2-subjettiness, low pileup", 110., 0., 1.1);
  nSubj3_PULo_hist = new TH1F("nSubj3_PULo_hist", "3-subjettiness, low pileup", 110., 0., 1.1);
  nSubj4_PULo_hist = new TH1F("nSubj4_PULo_hist", "4-subjettiness, low pileup", 110., 0., 1.1);
  t3t1_PULo_Ratio = new TH1F("t3t1_PULo_Ratio", "Ratio of tau_3 to tau_1, low pileup", 15., 0., 1.5);
  t2t1_PULo_Ratio = new TH1F("t2t1_PULo_Ratio", "Ratio of tau_2 to tau_1, low pileup", 15., 0., 1.5);
  t3t2_PULo_Ratio = new TH1F("t3t2_PULo_Ratio", "Ratio of tau_3 to tau_2, low pileup", 15., 0., 1.5);
  t4t1_PULo_Ratio = new TH1F("t4t1_PULo_Ratio", "Ratio of tau_4 to tau_1, low pileup", 15., 0., 1.5);
  t1t2_PULo_Ratio = new TH1F("t1t2_PULo_Ratio", "Ratio of tau_1 to tau_2, low pileup", 15., 0., 1.5);
  t2t3_PULo_Ratio = new TH1F("t2t3_PULo_Ratio", "Ratio of tau_2 to tau_3, low pileup", 15., 0., 1.5);
  t3t4_PULo_Ratio = new TH1F("t3t4_PULo_Ratio", "Ratio of tau_3 to tau_4, low pileup", 15., 0., 1.5);

  // High pileup (=> 20)
  nSubj1_PUHi_hist = new TH1F("nSubj1_PUHi_hist", "1-subjettiness, high pileup", 110., 0., 1.1);
  nSubj2_PUHi_hist = new TH1F("nSubj2_PUHi_hist", "2-subjettiness, high pileup", 110., 0., 1.1);
  nSubj3_PUHi_hist = new TH1F("nSubj3_PUHi_hist", "3-subjettiness, high pileup", 110., 0., 1.1);
  nSubj4_PUHi_hist = new TH1F("nSubj4_PUHi_hist", "4-subjettiness, high pileup", 110., 0., 1.1);
  t3t1_PUHi_Ratio = new TH1F("t3t1_PUHi_Ratio", "Ratio of tau_3 to tau_1, high pileup", 15., 0., 1.5);
  t2t1_PUHi_Ratio = new TH1F("t2t1_PUHi_Ratio", "Ratio of tau_2 to tau_1, high pileup", 15., 0., 1.5);
  t3t2_PUHi_Ratio = new TH1F("t3t2_PUHi_Ratio", "Ratio of tau_3 to tau_2, high pileup", 15., 0., 1.5);
  t4t1_PUHi_Ratio = new TH1F("t4t1_PUHi_Ratio", "Ratio of tau_4 to tau_1, high pileup", 15., 0., 1.5);
  t1t2_PUHi_Ratio = new TH1F("t1t2_PUHi_Ratio", "Ratio of tau_1 to tau_2, high pileup", 15., 0., 1.5);
  t2t3_PUHi_Ratio = new TH1F("t2t3_PUHi_Ratio", "Ratio of tau_2 to tau_3, high pileup", 15., 0., 1.5);
  t3t4_PUHi_Ratio = new TH1F("t3t4_PUHi_Ratio", "Ratio of tau_3 to tau_4, high pileup", 15., 0., 1.5);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FastJetAnalyzer_alltaus::endJob() 
{
  out_->cd();
  // jetTree->Write();
  jetPt_hist->Write();
  jetMass_hist->Write();
  intime_PU_vertices->Write();
  cands_beforeGrooming->Write();
  cands_afterGrooming->Write();

  if (pjmj_PU0->Integral() != 0.)
    pjmj_PU0->Scale(1./pjmj_PU0->Integral());
  pjmj_PU0->Write();
  if (pjmj_PUAll->Integral() != 0.)
    pjmj_PUAll->Scale(1./pjmj_PUAll->Integral());
  pjmj_PUAll->Write();

  tau1_One->Write();
  npIT_One->Write();
  cands_One->Write();
  pT_One->Write();
  tau1_NotOne->Write();
  npIT_NotOne->Write();
  cands_NotOne->Write();
  pT_NotOne->Write();
  pT_Groomed->Write();
  pT_highest_One->Write();
  pT_secondhighest_One->Write();
  maxDelR_One->Write();
  pT_highest_NotOne->Write();
  pT_secondhighest_NotOne->Write();
  maxDelR_NotOne->Write();
  delR_sisters_One->Write();
  delR_sisters_NotOne->Write();
  decayModes_One->Write();
  decayModes_NotOne->Write();

  //TCanvas *TwoD_PU0 = new TCanvas("TwoD_PU0", "nSubj_vs_pjmj_PU0", 750, 750);
  //TwoD_PU0->cd();
  nSubj_vs_pjmj_PU0->Write();
  //out_->cd();
  //TwoD_PU0->Write();

  //TCanvas *TwoD_PULo = new TCanvas("TwoD_PULo", "nSubj_vs_pjmj_PULo", 750, 750);
  //TwoD_PULo->cd();
  nSubj_vs_pjmj_PULo->Write();
  //out_->cd();
  //TwoD_PULo->Write();

  //TCanvas *TwoD_PUHi = new TCanvas("TwoD_PUHi", "nSubj_vs_pjmj_PUHi", 750, 750);
  //TwoD_PUHi->cd();
  nSubj_vs_pjmj_PUHi->Write();
  //out_->cd();
  //TwoD_PUHi->Write();

  //TCanvas *TwoD_PUAll = new TCanvas("TwoD_PUAll", "nSubj_vs_pjmj_PUAll", 750, 750);
  //TwoD_PUAll->cd();
  nSubj_vs_pjmj_PUAll->Write();
  //out_->cd();
  //TwoD_PUAll->Write();
  t3t1_vs_t4t1->Write();
  t3t1_vs_t2t1->Write();
  t2t1_vs_t4t1->Write();

  if (t3t1_PUAll_Ratio->Integral()!= 0.)
    t3t1_PUAll_Ratio->Scale(1./t3t1_PUAll_Ratio->Integral());
  t3t1_PUAll_Ratio->Write();
  t3t1_ee->Write();
  t3t1_mm->Write();
  t3t1_em->Write();
  t3t1_eh->Write();
  t3t1_mh->Write();
  t3t1_hh->Write();

  // zero pileup
  nSubj1_PU0_hist->Write();
  nSubj2_PU0_hist->Write();
  nSubj3_PU0_hist->Write();
  nSubj4_PU0_hist->Write();
  if (t3t1_PU0_Ratio->Integral()!=0.)
    t3t1_PU0_Ratio->Scale(1./t3t1_PU0_Ratio->Integral());
  t3t1_PU0_Ratio->Write();
  if (t2t1_PU0_Ratio->Integral()!=0.)
    t2t1_PU0_Ratio->Scale(1./t2t1_PU0_Ratio->Integral());
  t2t1_PU0_Ratio->Write();
  if (t3t2_PU0_Ratio->Integral()!=0.)
    t3t2_PU0_Ratio->Scale(1./t3t2_PU0_Ratio->Integral());
  t3t2_PU0_Ratio->Write();
  if (t4t1_PU0_Ratio->Integral()!=0.)
    t4t1_PU0_Ratio->Scale(1./t4t1_PU0_Ratio->Integral());
  t4t1_PU0_Ratio->Write();
  if (t1t2_PU0_Ratio->Integral()!=0.)
    t1t2_PU0_Ratio->Scale(1./t1t2_PU0_Ratio->Integral());
  t1t2_PU0_Ratio->Write();
  if (t2t3_PU0_Ratio->Integral()!=0.)
    t2t3_PU0_Ratio->Scale(1./t2t3_PU0_Ratio->Integral());
  t2t3_PU0_Ratio->Write();
  if (t3t4_PU0_Ratio->Integral()!=0.)
    t3t4_PU0_Ratio->Scale(1./t3t4_PU0_Ratio->Integral());
  t3t4_PU0_Ratio->Write();

  // low pileup
  nSubj1_PULo_hist->Write();
  nSubj2_PULo_hist->Write();
  nSubj3_PULo_hist->Write();
  nSubj4_PULo_hist->Write();
  if (t3t1_PULo_Ratio->Integral()!=0.)
    t3t1_PULo_Ratio->Scale(1./t3t1_PULo_Ratio->Integral());
  t3t1_PULo_Ratio->Write();
  if (t2t1_PULo_Ratio->Integral()!=0.)
    t2t1_PULo_Ratio->Scale(1./t2t1_PULo_Ratio->Integral());
  t2t1_PULo_Ratio->Write();
  if (t3t2_PULo_Ratio->Integral()!=0.)
    t3t2_PULo_Ratio->Scale(1./t3t2_PULo_Ratio->Integral());
  t3t2_PULo_Ratio->Write();
  if (t4t1_PULo_Ratio->Integral()!=0.)
    t4t1_PULo_Ratio->Scale(1./t4t1_PULo_Ratio->Integral());
  t4t1_PULo_Ratio->Write();
  if (t1t2_PULo_Ratio->Integral()!=0.)
    t1t2_PULo_Ratio->Scale(1./t1t2_PULo_Ratio->Integral());
  t1t2_PULo_Ratio->Write();
  if (t2t3_PULo_Ratio->Integral()!=0.)
    t2t3_PULo_Ratio->Scale(1./t2t3_PULo_Ratio->Integral());
  t2t3_PULo_Ratio->Write();
  if (t3t4_PULo_Ratio->Integral()!=0.)
    t3t4_PULo_Ratio->Scale(1./t3t4_PULo_Ratio->Integral());
  t3t4_PULo_Ratio->Write();

  // high pileup
  nSubj1_PUHi_hist->Write();
  nSubj2_PUHi_hist->Write();
  nSubj3_PUHi_hist->Write();
  nSubj4_PUHi_hist->Write();
  if (t3t1_PUHi_Ratio->Integral()!=0.)
    t3t1_PUHi_Ratio->Scale(1./t3t1_PUHi_Ratio->Integral());
  t3t1_PUHi_Ratio->Write();
  if (t2t1_PUHi_Ratio->Integral()!=0.)
    t2t1_PUHi_Ratio->Scale(1./t2t1_PUHi_Ratio->Integral());
  t2t1_PUHi_Ratio->Write();
  if (t3t2_PUHi_Ratio->Integral()!=0.)
    t3t2_PUHi_Ratio->Scale(1./t3t2_PUHi_Ratio->Integral());
  t3t2_PUHi_Ratio->Write();
  if (t4t1_PUHi_Ratio->Integral()!=0.)
    t4t1_PUHi_Ratio->Scale(1./t4t1_PUHi_Ratio->Integral());
  t4t1_PUHi_Ratio->Write();
  if (t1t2_PUHi_Ratio->Integral()!=0.)
    t1t2_PUHi_Ratio->Scale(1./t1t2_PUHi_Ratio->Integral());
  t1t2_PUHi_Ratio->Write();
  if (t2t3_PUHi_Ratio->Integral()!=0.)
    t2t3_PUHi_Ratio->Scale(1./t2t3_PUHi_Ratio->Integral());
  t2t3_PUHi_Ratio->Write();
  if (t3t4_PUHi_Ratio->Integral()!=0.)
    t3t4_PUHi_Ratio->Scale(1./t3t4_PUHi_Ratio->Integral());
  t3t4_PUHi_Ratio->Write();


  out_->Write();
  out_->Close();
  SurvivingAxes.close();
}

// ------------ method called when starting to processes a run  ------------
void 
FastJetAnalyzer_alltaus::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
FastJetAnalyzer_alltaus::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
FastJetAnalyzer_alltaus::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
FastJetAnalyzer_alltaus::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FastJetAnalyzer_alltaus::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FastJetAnalyzer_alltaus);
