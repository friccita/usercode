// -*- C++ -*-
//
// Package:    TauAnalyzer
// Class:      TauAnalyzer
// 
/**\class TauAnalyzer TauAnalyzer.cc BoostedTauAnalysis/TauAnalyzer/src/TauAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 16:40:51 CEST 2012
// $Id: TauAnalyzer.cc,v 1.2 2012/08/10 13:50:38 yohay Exp $
//
//

/* understand triggerability of signal
   - ZH: fraction of events passing Z-->ee/mumu triggers
   - gg fusion: fraction of events passing mu + had trigger
   - first plot raw number, then vs. offline trigger object pT/eta cut and NPV/rho

   nice triggers in use:
   - HLT_Mu17_Mu8
   - HLT_Mu8_DiJet30
   - HLT_IsoMu12_DoubleCentralJet65
   - HLT_IsoMu17_eta2p1_CentralPFNoPUJet30
   - HLT_IsoMu17_eta2p1_DiCentralPFNoPUJet30
   - HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30
   - HLT_IsoMu24
   - HLT_Mu40
   - HLT_IsoMu17_eta2p1_LooseIsoPFTau20
   - HLT_LooseIsoPFTau35_Trk20_Prong1
   - HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Prong1

   make nice plotting code

   identify a-->tautau-->muhadX decays, and plot:
   - gen mu pT

   identify a-->tautau-->muhadX decays and plot:
   do same for matched HPS taus
   do same for matched HPS taus that were rerun with pT sorting, not isolation sorting
   do same for matched HPS taus that were rerun with isolation sorting excluding high energy muons 
   (how high is high?  get it from the mu pT plot)
   do same for matched HPS taus that were rerun with sorting that prefers high energy muons in 
   isolation cone (how high is high?  get it from the mu pT plot)

   identify a-->tautau-->muhadX decays, and for PF jets matched to the had part, plot:
   - fraction of time it's also matched to the mu part
   - visible pT
   - eta
   - MET
   - pT of the muon matched to the gen mu in the decay
   - energy fractions (esp. of muons and charged tracks)
   - multiplicities (esp. of muons and charged tracks)
   - dR to the PF jet matched to the other a
   - dR between jet and highest energy constituent muon
   - had+mu invariant mass
   - had+mu pT
   do same for matched PF taus
   do same for matched HPS taus
   do same for matched HPS taus that were rerun with pT sorting, not isolation sorting
   do same for matched HPS taus that were rerun with isolation sorting excluding high energy muons 
   (how high is high?  get it from the mu pT plot)
   do same for matched HPS taus that were rerun with sorting that prefers high energy muons in 
   isolation cone (how high is high?  get it from the mu pT plot)
   compare to same for PF jets matched to Z-->tau-->had decays

   backgrounds to the gg fusion search:
   - big one is mu-enriched QCD

   repeat for jets in mu-enriched QCD that have a real muon from heavy flavor decay
   repeat for jets in regular QCD

   try to develop an ID that can at least cut out some QCD
 */

// system include files
#include <memory>
#include <string>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
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
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
//
// class declaration
//

class TauAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TauAnalyzer(const edm::ParameterSet&);
      ~TauAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      void reset();
      double dR(const std::vector<double>& dR2, const unsigned int pos = 0) const;
      void fillDRHistogramNearFar(const unsigned int, const bool, 
				  const std::vector<double>&) const;
      void fillDRHistogramLeadTrail(const bool, const unsigned int, 
				    const std::vector<double>&) const;
      void fillDRHistogramHighLowPT(const bool, const std::vector<double>&) const;
      template<typename T, typename U>
        const unsigned int numMatchingObjects(const std::vector<T*>& objectsToMatch, 
					      const U& referenceObject) const
        {
          std::vector<double> dR2Min = 
          Common::sortByProximity(objectsToMatch, referenceObject);
          std::vector<double>::const_iterator iDR2 = dR2Min.begin();
          while ((iDR2 != dR2Min.end()) && (dR(dR2Min, iDR2 - dR2Min.begin()) < 0.3)) {
	    ++iDR2;
          }
          return (iDR2 - dR2Min.begin());
        }

      // ----------member data ---------------------------

      //output
      TFile* out_;

      //input
      std::string outFileName_;
      edm::InputTag genParticleTag_;
      edm::InputTag tauTag_;
      edm::InputTag muonTag_;
      edm::InputTag vtxTag_;
      std::vector<edm::InputTag> HPSDiscriminatorTags_;
      std::map<std::string, edm::Handle<reco::PFTauDiscriminator> > HPSDiscriminators_;
      unsigned int momPDGID_;
      double genMuTauPTMin_;
      double genMuPTMin_;
      double effVsEtaPTMin_;
      edm::ParameterSet* cfg_;

      //histograms
      TH1F* status3TauMultiplicity_;
      TH1F* HPSTauMultiplicity_;
      TH1F* HPSTausPassingEtaCutMultiplicity_;
      TH1F* HPSIsoTausPassingEtaCutMultiplicity_;
      TH1F* dRToHPSTau0_;
      TH1F* dRToHPSTau1_;
      TH1F* dRToHPSIsoTau0_;
      TH1F* dRToHPSIsoTau1_;
      TH2F* HPSTauMultiplicityVsDRToHPSTau0_;
      TH2F* HPSIsoTauMultiplicityVsDRToHPSIsoTau0_;
      TH1F* dRLeadingGenTauToHPSTau0_;
      TH1F* dRTrailingGenTauToHPSTau0_;
      TH1F* dRLeadingGenTauToHPSIsoTau0_;
      TH1F* dRTrailingGenTauToHPSIsoTau0_;
      TH1F* PFTauPT_;
      TH1F* HPSTauPT_;
      TH1F* HPSIsoTauPT_;
      TH1F* HPSIsoTauPTPassingEtaCut_;
      TH1F* genTauPT_;
      TH1F* dRToHPSTau0GenPTGe30_;
      TH1F* dRToHPSIsoTau0GenPTGe30_;
      TH1F* HPSTauEta_;
      TH1F* genTauEta_;
      TH1F* muHadPFTauMatchMultiplicity_;
      TH1F* muHadPFTauIsoMatchMultiplicity_;
      TH1F* muHadMuonMatchMultiplicity_;
      TH1F* muHadPFTauMatchMuonMatchMultiplicity_;
      TH1F* muHadPFTauIsoMatchMuonMatchMultiplicity_;
      TH1F* muHadVisibleGenPT_;
      TH1F* muHadPFTauMatchVisibleGenPT_;
      TH1F* muHadPFTauIsoMatchVisibleGenPT_;
      TH1F* muHadMuMatchVisibleGenPT_;
      TH1F* muHadPFTauMatchMuMatchVisibleGenPT_;
      TH1F* muHadPFTauIsoMatchMuMatchVisibleGenPT_;
      TH1F* muHadVisibleGenEta_;
      TH1F* muHadPFTauMatchVisibleGenEta_;
      TH1F* muHadPFTauIsoMatchVisibleGenEta_;
      TH1F* muHadMuMatchVisibleGenEta_;
      TH1F* muHadPFTauMatchMuMatchVisibleGenEta_;
      TH1F* muHadPFTauIsoMatchMuMatchVisibleGenEta_;
      TH1F* muHadGenMuPT_;
      TH1F* muHadPFTauMatchGenMuPT_;
      TH1F* muHadPFTauIsoMatchGenMuPT_;
      TH1F* muHadMuMatchGenMuPT_;
      TH1F* muHadPFTauMatchMuMatchGenMuPT_;
      TH1F* muHadPFTauIsoMatchMuMatchGenMuPT_;
      TH1F* muHadGenMuEta_;
      TH1F* muHadPFTauMatchGenMuEta_;
      TH1F* muHadPFTauIsoMatchGenMuEta_;
      TH1F* muHadMuMatchGenMuEta_;
      TH1F* muHadPFTauMatchMuMatchGenMuEta_;
      TH1F* muHadPFTauIsoMatchMuMatchGenMuEta_;
      TH1F* muHadGenDR_;
      TH1F* muHadPFTauMatchGenDR_;
      TH1F* muHadPFTauIsoMatchGenDR_;
      TH1F* muHadMuMatchGenDR_;
      TH1F* muHadPFTauMatchMuMatchGenDR_;
      TH1F* muHadPFTauIsoMatchMuMatchGenDR_;

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
TauAnalyzer::TauAnalyzer(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
  tauTag_(iConfig.getParameter<edm::InputTag>("tauTag")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonTag")),
  vtxTag_(iConfig.getParameter<edm::InputTag>("vtxTag")),
  HPSDiscriminatorTags_(iConfig.getParameter<std::vector<edm::InputTag> >("HPSDiscriminatorTags")),
  momPDGID_(iConfig.getParameter<unsigned int>("momPDGID")),
  genMuTauPTMin_(iConfig.getParameter<double>("genMuTauPTMin")),
  genMuPTMin_(iConfig.getParameter<double>("genMuPTMin")),
  effVsEtaPTMin_(iConfig.getParameter<double>("effVsEtaPTMin")),
  cfg_(const_cast<edm::ParameterSet*>(&iConfig))
{
   //now do what ever initialization is needed

}


TauAnalyzer::~TauAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   reset();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
TauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get GEN particles
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByLabel(genParticleTag_, pGenParticles);

  //get HPS PF taus
  edm::Handle<reco::PFTauCollection> pTaus;
  iEvent.getByLabel(tauTag_, pTaus);

  //get HPS discriminators
  for (std::vector<edm::InputTag>::const_iterator 
	 iHPSDiscriminatorTag = HPSDiscriminatorTags_.begin(); 
       iHPSDiscriminatorTag != HPSDiscriminatorTags_.end(); ++iHPSDiscriminatorTag) {
    edm::Handle<reco::PFTauDiscriminator> pHPSDiscriminator;
    iEvent.getByLabel(*iHPSDiscriminatorTag, pHPSDiscriminator);
    HPSDiscriminators_[iHPSDiscriminatorTag->label()] = pHPSDiscriminator;
  }

  //get muons
  edm::Handle<reco::MuonCollection> pMuons;
  iEvent.getByLabel(muonTag_, pMuons);

  //get vertices
  edm::Handle<reco::VertexCollection> pVertices;
  iEvent.getByLabel(vtxTag_, pVertices);

  //HPS tau multiplicity
  unsigned int nHPSTaus = 0;
  unsigned int nHPSTausPassingEtaCut = 0;
  unsigned int nHPSIsoTausPassingEtaCut = 0;
  std::vector<reco::PFTau*> PFTaus;
  std::vector<reco::PFTau*> HPSIsoTaus;
  std::vector<reco::PFTau*> HPSTaus;
  for (reco::PFTauCollection::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
       ++iTau) {
    PFTaus.push_back(const_cast<reco::PFTau*>(&*iTau));
    reco::PFTauRef tauRef(pTaus, iTau - pTaus->begin());
    std::map<std::string, edm::Handle<reco::PFTauDiscriminator> >::const_iterator 
      decayModeFinding = HPSDiscriminators_.find("hpsPFTauDiscriminationByDecayModeFinding");
    if ((decayModeFinding != HPSDiscriminators_.end()) && 
	((*(decayModeFinding->second))[tauRef] == 1.0)) {
      ++nHPSTaus;
      std::map<std::string, edm::Handle<reco::PFTauDiscriminator> >::const_iterator 
	looseIsoMVA = HPSDiscriminators_.find("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr");
      if (fabs(iTau->eta()) < 2.3) {
	++nHPSTausPassingEtaCut;
	if ((looseIsoMVA != HPSDiscriminators_.end()) && 
	    ((*(looseIsoMVA->second))[tauRef] == 1.0)) {
	  ++nHPSIsoTausPassingEtaCut;

	  //fill pT distribution of all isolated HPS taus passing eta cut
	  HPSIsoTauPTPassingEtaCut_->Fill(iTau->pt());
	}
      }

      //save isolated HPS tau
      if ((looseIsoMVA != HPSDiscriminators_.end()) && 
	  ((*(looseIsoMVA->second))[tauRef] == 1.0)) {
	HPSIsoTaus.push_back(const_cast<reco::PFTau*>(&*iTau));

	//fill pT distribution of all isolated HPS taus
	HPSIsoTauPT_->Fill(iTau->pt());
      }

      //save non-isolated HPS tau
      HPSTaus.push_back(const_cast<reco::PFTau*>(&*iTau));

      //fill kinematic distributions of all HPS taus
      HPSTauPT_->Fill(iTau->pt());
      HPSTauEta_->Fill(iTau->eta());
    }

    //fill pT distribution of all PF taus
    PFTauPT_->Fill(iTau->pt());
  }

  //fill multiplicity plots
  HPSTauMultiplicity_->Fill(nHPSTaus);
  HPSTausPassingEtaCutMultiplicity_->Fill(nHPSTausPassingEtaCut);
  HPSIsoTausPassingEtaCutMultiplicity_->Fill(nHPSIsoTausPassingEtaCut);

  //get multiplicity of status 3 taus and Z-->tautau decays
  unsigned int nStatus3TausPerEvt = 0;
  std::vector<reco::Candidate*> ZDecayTaus;
  std::vector<GenTauDecayID> tauDecays;
  for (reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin(); 
       iGenParticle != pGenParticles->end(); ++iGenParticle) {

    //instantiate GenTauDecayID object
    try {
      GenTauDecayID tauDecay(*cfg_, pGenParticles, iGenParticle - pGenParticles->begin());

      //look for a status 3 tau from boson decay
      if (tauDecay.tauIsStatus3DecayProduct()) {
	++nStatus3TausPerEvt;
	const reco::Candidate* pZDecayTau = dynamic_cast<const reco::Candidate*>(&*iGenParticle);
	ZDecayTaus.push_back(const_cast<reco::Candidate*>(pZDecayTau));

	//sort non-isolated HPS taus by proximity to GEN tau
	std::vector<double> dR2MinHPSTau = Common::sortByProximity(HPSTaus, *iGenParticle);
	fillDRHistogramNearFar(0, false, dR2MinHPSTau);
	fillDRHistogramNearFar(1, false, dR2MinHPSTau);

	//sort isolated HPS taus by proximity to GEN tau
	std::vector<double> dR2MinHPSIsoTau = Common::sortByProximity(HPSIsoTaus, *iGenParticle);
	fillDRHistogramNearFar(0, true, dR2MinHPSIsoTau);
	fillDRHistogramNearFar(1, true, dR2MinHPSIsoTau);

	//min dR for gen tau pT > 30 GeV
	if (iGenParticle->pt() > 30.0) {
	  fillDRHistogramHighLowPT(false, dR2MinHPSTau);
	  fillDRHistogramHighLowPT(true, dR2MinHPSIsoTau);
	}

	//fill kinematic distributions of all GEN taus from desired resonance decay (Z or a)
	genTauPT_->Fill(iGenParticle->pt());
	genTauEta_->Fill(iGenParticle->eta());

	//save the boson-matched taus for later use
	tauDecays.push_back(tauDecay);
      }
    }
    catch (std::string& ex) { throw cms::Exception("TauAnalyzer") << ex; }
  }
  status3TauMultiplicity_->Fill(nStatus3TausPerEvt);

  //dR MC matching plots split by leading and trailing GEN tau
  Common::sortByPT(ZDecayTaus);
  if (ZDecayTaus.size() >= 2) {
    std::vector<double> dR2MinLeadingGenTauHPSTau = 
      Common::sortByProximity(HPSTaus, *(ZDecayTaus[0]));
    std::vector<double> dR2MinTrailingGenTauHPSTau = 
      Common::sortByProximity(HPSTaus, *(ZDecayTaus[1]));
    std::vector<double> dR2MinLeadingGenTauHPSIsoTau = 
      Common::sortByProximity(HPSIsoTaus, *(ZDecayTaus[0]));
    std::vector<double> dR2MinTrailingGenTauHPSIsoTau = 
      Common::sortByProximity(HPSIsoTaus, *(ZDecayTaus[1]));
    fillDRHistogramLeadTrail(false, 1, dR2MinLeadingGenTauHPSTau);
    fillDRHistogramLeadTrail(false, 2, dR2MinTrailingGenTauHPSTau);
    fillDRHistogramLeadTrail(true, 1, dR2MinLeadingGenTauHPSIsoTau);
    fillDRHistogramLeadTrail(true, 2, dR2MinTrailingGenTauHPSIsoTau);
  }

  //identify the first good vertex (the "primary" (?))
  reco::VertexCollection::const_iterator iVtx = pVertices->begin();
  reco::Vertex* pPV = NULL;
  while ((iVtx != pVertices->end()) && (pPV == NULL)) {
    if (!iVtx->isFake() && 
	(iVtx->ndof() > 4) && 
	(fabs(iVtx->x()) <= 24.0/*cm*/) && 
	(fabs(iVtx->position().Rho()) <= 2.0/*cm*/)) pPV = const_cast<reco::Vertex*>(&*iVtx);
    ++iVtx;
  }

  /*create a collection of muons passing the 2012 tight selection (cf. 
    https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId)*/
  std::vector<reco::Muon*> tightMuons;
  for (reco::MuonCollection::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); 
       ++iMuon) {
    if (muon::isTightMuon(*iMuon, *pPV) && 
	iMuon->isPFMuon() && 
	(fabs(iMuon->innerTrack()->dz(pPV->position())) < 0.5) && 
	(iMuon->track()->hitPattern().trackerLayersWithMeasurement() > 5)) {
      tightMuons.push_back(const_cast<reco::Muon*>(&*iMuon));
    }
  }

  //loop over taus from boson decay
  std::vector<unsigned int> keysToIgnore;
  for (std::vector<GenTauDecayID>::iterator iTau = tauDecays.begin(); 
       iTau != tauDecays.end(); ++iTau) {

    //try for exceptions
    try {

      //look for a hadronic tau decay from the signal sample
      if ((iTau->tauDecayType() == GenTauDecayID::HAD) && (momPDGID_ == GenTauDecayID::APDGID)) { // for a
      //  if (iTau->tauDecayType() == GenTauDecayID::HAD) { // for Z

	//look for the other tau decay product of the a: is it a tau-->mu decay?
	iTau->findSister();
	const unsigned int iSister = iTau->getSisterIndex();
	if ((std::find(keysToIgnore.begin(), keysToIgnore.end(), iSister) == 
	     keysToIgnore.end()) && (iTau->sisterDecayType() == GenTauDecayID::MU)) {

	  //does the sister decay mu pass the pT threshold?
	  reco::GenParticleRef sisterRef(pGenParticles, iSister);
	  if (sisterRef->pt() > genMuTauPTMin_) {

	    //get the visible had 4-vector
	    reco::GenParticleRef tauRef(pGenParticles, iTau->getTauIndex());
	    reco::GenParticleRef status2HadTauRef = tauRef->daughterRef(0);
	    reco::LeafCandidate::LorentzVector visibleHadTauP4;
	    for (unsigned int iDaughter = 0; iDaughter < status2HadTauRef->numberOfDaughters(); 
		 ++iDaughter) {
	      reco::GenParticleRef hadTauDaughterRef = status2HadTauRef->daughterRef(iDaughter);
	      if ((fabs(hadTauDaughterRef->pdgId()) != GenTauDecayID::ENEUTRINOPDGID) && 
		  (fabs(hadTauDaughterRef->pdgId()) != GenTauDecayID::MUNEUTRINOPDGID) && 
		  (fabs(hadTauDaughterRef->pdgId()) != GenTauDecayID::TAUNEUTRINOPDGID)) {
		visibleHadTauP4+=hadTauDaughterRef->p4();
	      }
	    }
	    const double visibleHadTauPT = visibleHadTauP4.Pt();
	    const double visibleHadTauEta = visibleHadTauP4.Eta();

	    //get the muon 4-vector
	    int iGenMu = -1;
	    unsigned int iDaughter = 0;
	    reco::GenParticleRef status2MuTauRef = sisterRef->daughterRef(0);
	    while ((iDaughter < status2MuTauRef->numberOfDaughters()) && (iGenMu == -1)) {
	      reco::GenParticleRef muTauDaughter = status2MuTauRef->daughterRef(iDaughter);
	      if (fabs(muTauDaughter->pdgId()) == GenTauDecayID::MUPDGID) {
		iGenMu = muTauDaughter.key();
	      }
	      ++iDaughter;
	    }
	    if (iGenMu == -1) throw cms::Exception("TauAnalyzer") << "Muon not found.\n";
	    reco::GenParticleRef muRef(pGenParticles, iGenMu);
	    const double genMuPT = muRef->pt();
	    const double genMuEta = muRef->eta();

	    //get dR between visible parts of hadronic tau and muonic tau
	    const double muHadDR = 
	      reco::deltaR(visibleHadTauEta, visibleHadTauP4.Phi(), genMuEta, muRef->phi());

	    //fill denominator histograms
	    if (genMuPT > genMuPTMin_) {
	      muHadVisibleGenPT_->Fill(visibleHadTauPT);
	      if (visibleHadTauPT > effVsEtaPTMin_) muHadVisibleGenEta_->Fill(visibleHadTauEta);
	      muHadGenMuPT_->Fill(genMuPT);
	      if (genMuPT > effVsEtaPTMin_) muHadGenMuEta_->Fill(genMuEta);
	      muHadGenDR_->Fill(muHadDR);
	    }

	    //how many matching PF taus?
	    const unsigned int nMatchingTaus = numMatchingObjects(PFTaus, *tauRef);
	    muHadPFTauMatchMultiplicity_->Fill(nMatchingTaus);

	    //how many matching PF muons?
	    const unsigned int nMatchingMuons = numMatchingObjects(tightMuons, *muRef);
	    muHadMuonMatchMultiplicity_->Fill(nMatchingMuons);

	    //how many matching PF taus (LooseIsoMVA)?
	    const unsigned int nMatchingIsoTaus = numMatchingObjects(HPSIsoTaus, *tauRef);
	    muHadPFTauIsoMatchMultiplicity_->Fill(nMatchingIsoTaus);

	    //how many with >1 PF tau match and >1 muon match?
	    if (nMatchingTaus > 0) muHadPFTauMatchMuonMatchMultiplicity_->Fill(nMatchingMuons);

	    //how many with >1 PF tau (LooseIsoMVA) match and >1 muon match?
	    if (nMatchingIsoTaus > 0) muHadPFTauIsoMatchMuonMatchMultiplicity_->Fill(nMatchingMuons);

	    //fill numerator histograms
	    if (genMuPT > genMuPTMin_) {
	      if (nMatchingTaus > 0) {
		muHadPFTauMatchVisibleGenPT_->Fill(visibleHadTauPT);
		if (visibleHadTauPT > effVsEtaPTMin_) {
		  muHadPFTauMatchVisibleGenEta_->Fill(visibleHadTauEta);
		}
		muHadPFTauMatchGenMuPT_->Fill(genMuPT);
		if (genMuPT > effVsEtaPTMin_) muHadPFTauMatchGenMuEta_->Fill(genMuEta);
		muHadPFTauMatchGenDR_->Fill(muHadDR);
	      }
	      if (nMatchingMuons > 0) {
		muHadMuMatchVisibleGenPT_->Fill(visibleHadTauPT);
		if (visibleHadTauPT > effVsEtaPTMin_) {
		  muHadMuMatchVisibleGenEta_->Fill(visibleHadTauEta);
		}
		muHadMuMatchGenMuPT_->Fill(genMuPT);
		if (genMuPT > effVsEtaPTMin_) muHadMuMatchGenMuEta_->Fill(genMuEta);
		muHadMuMatchGenDR_->Fill(muHadDR);
	      }
	      if (nMatchingIsoTaus > 0) {
		muHadPFTauIsoMatchVisibleGenPT_->Fill(visibleHadTauPT);
		if (visibleHadTauPT > effVsEtaPTMin_) {
		  muHadPFTauIsoMatchVisibleGenEta_->Fill(visibleHadTauEta);
		}
		muHadPFTauIsoMatchGenMuPT_->Fill(genMuPT);
		if (genMuPT > effVsEtaPTMin_) muHadPFTauIsoMatchGenMuEta_->Fill(genMuEta);
		muHadPFTauIsoMatchGenDR_->Fill(muHadDR);
	      }
	      if ((nMatchingTaus > 0) && (nMatchingMuons > 0)) {
		muHadPFTauMatchMuMatchVisibleGenPT_->Fill(visibleHadTauPT);
		if (visibleHadTauPT > effVsEtaPTMin_) {
		  muHadPFTauMatchMuMatchVisibleGenEta_->Fill(visibleHadTauEta);
		}
		muHadPFTauMatchMuMatchGenMuPT_->Fill(genMuPT);
		if (genMuPT > effVsEtaPTMin_) muHadPFTauMatchMuMatchGenMuEta_->Fill(genMuEta);
		muHadPFTauMatchMuMatchGenDR_->Fill(muHadDR);
	      }
	      if ((nMatchingIsoTaus > 0) && (nMatchingMuons > 0)) {
		muHadPFTauIsoMatchMuMatchVisibleGenPT_->Fill(visibleHadTauPT);
		if (visibleHadTauPT > effVsEtaPTMin_) {
		  muHadPFTauIsoMatchMuMatchVisibleGenEta_->Fill(visibleHadTauEta);
		}
		muHadPFTauIsoMatchMuMatchGenMuPT_->Fill(genMuPT);
		if (genMuPT > effVsEtaPTMin_) muHadPFTauIsoMatchMuMatchGenMuEta_->Fill(genMuEta);
		muHadPFTauIsoMatchMuMatchGenDR_->Fill(muHadDR);
	      }
	    }
	  }

	  /*add this tau's key to the list of keys to ignore when we get to its sister so the 
	    histograms don't get filled twice*/
	  keysToIgnore.push_back(iTau->getTauIndex());
	}
      }
    }
    catch (std::string& ex) { throw cms::Exception("TauAnalyzer") << ex; }
  }

  //check distribution of min dR between gen jet and nearest reco tau
  //check distribution of min dR between gen jet and 2nd nearest reco tau

}


// ------------ method called once each job just before starting event loop  ------------
void 
TauAnalyzer::beginJob()
{
  //open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book histograms
  status3TauMultiplicity_ = new TH1F("status3TauMultiplicity", "", 7, -0.5, 6.5);
  HPSTauMultiplicity_ = new TH1F("HPSTauMultiplicity", "", 30, -0.5, 29.5);
  HPSTausPassingEtaCutMultiplicity_ = 
    new TH1F("HPSTausPassingEtaCutMultiplicity", "", 30, -0.5, 29.5);
  HPSIsoTausPassingEtaCutMultiplicity_ = 
    new TH1F("HPSIsoTausPassingEtaCutMultiplicity", "", 30, -0.5, 29.5);
  dRToHPSTau0_ = new TH1F("dRToHPSTau0", "", 500, 0.0, 10.0);
  dRToHPSTau1_ = new TH1F("dRToHPSTau1", "", 500, 0.0, 10.0);
  dRToHPSIsoTau0_ = new TH1F("dRToHPSIsoTau0", "", 500, 0.0, 10.0);
  dRToHPSIsoTau1_ = new TH1F("dRToHPSIsoTau1", "", 500, 0.0, 10.0);
  HPSTauMultiplicityVsDRToHPSTau0_ = 
    new TH2F("HPSTauMultiplicityVsDRToHPSTau0", "", 500, 0.0, 10.0, 30, -0.5, 29.5);
  HPSIsoTauMultiplicityVsDRToHPSIsoTau0_ = new TH2F("HPSIsoTauMultiplicityVsDRToHPSIsoTau0", "", 500, 0.0, 10.0, 30, -0.5, 29.5);
  dRLeadingGenTauToHPSTau0_ = new TH1F("dRLeadingGenTauToHPSTau0", "", 500, 0.0, 10.0);
  dRTrailingGenTauToHPSTau0_ = new TH1F("dRTrailingGenTauToHPSTau0", "", 500, 0.0, 10.0);
  dRLeadingGenTauToHPSIsoTau0_ = new TH1F("dRLeadingGenTauToHPSIsoTau0", "", 500, 0.0, 10.0);
  dRTrailingGenTauToHPSIsoTau0_ = new TH1F("dRTrailingGenTauToHPSIsoTau0", "", 500, 0.0, 10.0);
  PFTauPT_ = new TH1F("PFTauPT", "", 50, 0.0, 100.0);
  HPSTauPT_ = new TH1F("HPSTauPT", "", 50, 0.0, 100.0);
  HPSIsoTauPT_ = new TH1F("HPSIsoTauPT", "", 50, 0.0, 100.0);
  HPSIsoTauPTPassingEtaCut_ = new TH1F("HPSIsoTauPTPassingEtaCut", "", 50, 0.0, 100.0);
  genTauPT_ = new TH1F("genTauPT", "", 50, 0.0, 100.0);
  dRToHPSTau0GenPTGe30_ = new TH1F("dRToHPSTau0GenPTGe30", "", 500, 0.0, 10.0);
  dRToHPSIsoTau0GenPTGe30_ = new TH1F("dRToHPSIsoTau0GenPTGe30", "", 500, 0.0, 10.0);
  HPSTauEta_ = new TH1F("HPSTauEta", "", 100, -5.0, 5.0);
  genTauEta_ = new TH1F("genTauEta", "", 100, -5.0, 5.0);
  muHadPFTauMatchMultiplicity_ = new TH1F("muHadPFTauMatchMultiplicity", "", 4, -0.5, 3.5);
  muHadPFTauIsoMatchMultiplicity_ = new TH1F("muHadPFTauIsoMatchMultiplicity", "", 4, -0.5, 3.5);
  muHadMuonMatchMultiplicity_ = new TH1F("muHadMuonMatchMultiplicity", "", 4, -0.5, 3.5);
  muHadPFTauMatchMuonMatchMultiplicity_ = new TH1F("muHadPFTauMatchMuonMatchMultiplicity", "", 4, -0.5, 3.5);
  muHadPFTauIsoMatchMuonMatchMultiplicity_ = new TH1F("muHadPFTauIsoMatchMuonMatchMultiplicity", "", 4, -0.5, 3.5);
  muHadVisibleGenPT_ = new TH1F("muHadVisibleGenPT", "", 10, 0.0, 100.0);
  muHadPFTauMatchVisibleGenPT_ = new TH1F("muHadPFTauMatchVisibleGenPT", "", 10, 0.0, 100.0);
  muHadPFTauIsoMatchVisibleGenPT_ = new TH1F("muHadPFTauIsoMatchVisibleGenPT", "", 10, 0.0, 100.0);
  muHadMuMatchVisibleGenPT_ = new TH1F("muHadMuMatchVisibleGenPT", "", 10, 0.0, 100.0);
  muHadPFTauMatchMuMatchVisibleGenPT_ = new TH1F("muHadPFTauMatchMuMatchVisibleGenPT", "", 10, 0.0, 100.0);
  muHadPFTauIsoMatchMuMatchVisibleGenPT_ = new TH1F("muHadPFTauIsoMatchMuMatchVisibleGenPT", "", 10, 0.0, 100.0);
  muHadVisibleGenEta_ = new TH1F("muHadVisibleGenEta", "", 10, -2.5, 2.5);
  muHadPFTauMatchVisibleGenEta_ = new TH1F("muHadPFTauMatchVisibleGenEta", "", 10, -2.5, 2.5);
  muHadPFTauIsoMatchVisibleGenEta_ = new TH1F("muHadPFTauIsoMatchVisibleGenEta", "", 10, -2.5, 2.5);
  muHadMuMatchVisibleGenEta_ = new TH1F("muHadMuMatchVisibleGenEta", "", 10, -2.5, 2.5);
  muHadPFTauMatchMuMatchVisibleGenEta_ = new TH1F("muHadPFTauMatchMuMatchVisibleGenEta", "", 10, -2.5, 2.5);
  muHadPFTauIsoMatchMuMatchVisibleGenEta_ = new TH1F("muHadPFTauIsoMatchMuMatchVisibleGenEta", "", 10, -2.5, 2.5);
  muHadGenMuPT_ = new TH1F("muHadGenMuPT", "", 10, 0.0, 100.0);
  muHadPFTauMatchGenMuPT_ = new TH1F("muHadPFTauMatchGenMuPT", "", 10, 0.0, 100.0);
  muHadPFTauIsoMatchGenMuPT_ = new TH1F("muHadPFTauIsoMatchGenMuPT", "", 10, 0.0, 100.0);
  muHadMuMatchGenMuPT_ = new TH1F("muHadMuMatchGenMuPT", "", 10, 0.0, 100.0);
  muHadPFTauMatchMuMatchGenMuPT_ = new TH1F("muHadPFTauMatchMuMatchGenMuPT", "", 10, 0.0, 100.0);
  muHadPFTauIsoMatchMuMatchGenMuPT_ = new TH1F("muHadPFTauIsoMatchMuMatchGenMuPT", "", 10, 0.0, 100.0);
  muHadGenMuEta_ = new TH1F("muHadGenMuEta", "", 10, -2.5, 2.5);
  muHadPFTauMatchGenMuEta_ = new TH1F("muHadPFTauMatchGenMuEta", "", 10, -2.5, 2.5);
  muHadPFTauIsoMatchGenMuEta_ = new TH1F("muHadPFTauIsoMatchGenMuEta", "", 10, -2.5, 2.5);
  muHadMuMatchGenMuEta_ = new TH1F("muHadMuMatchGenMuEta", "", 10, -2.5, 2.5);
  muHadPFTauMatchMuMatchGenMuEta_ = new TH1F("muHadPFTauMatchMuMatchGenMuEta", "", 10, -2.5, 2.5);
  muHadPFTauIsoMatchMuMatchGenMuEta_ = new TH1F("muHadPFTauIsoMatchMuMatchGenMuEta", "", 10, -2.5, 2.5);
  muHadGenDR_ = new TH1F("muHadGenDR", "", 15, 0.0, 3.0);
  muHadPFTauMatchGenDR_ = new TH1F("muHadPFTauMatchGenDR", "", 15, 0.0, 3.0);
  muHadPFTauIsoMatchGenDR_ = new TH1F("muHadPFTauIsoMatchGenDR", "", 15, 0.0, 3.0);
  muHadMuMatchGenDR_ = new TH1F("muHadMuMatchGenDR", "", 15, 0.0, 3.0);
  muHadPFTauMatchMuMatchGenDR_ = new TH1F("muHadPFTauMatchMuMatchGenDR", "", 15, 0.0, 3.0);
  muHadPFTauIsoMatchMuMatchGenDR_ = new TH1F("muHadPFTauIsoMatchMuMatchGenDR", "", 15, 0.0, 3.0);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TauAnalyzer::endJob() 
{
  //make efficiency plots
  TGraphAsymmErrors effPFTauMatchVsVisibleGenPT(muHadPFTauMatchVisibleGenPT_, muHadVisibleGenPT_);
  TGraphAsymmErrors effPFTauMatchVsVisibleGenEta(muHadPFTauMatchVisibleGenEta_, muHadVisibleGenEta_);
  TGraphAsymmErrors effPFTauIsoMatchVsVisibleGenPT(muHadPFTauIsoMatchVisibleGenPT_, muHadVisibleGenPT_);
  TGraphAsymmErrors effPFTauIsoMatchVsVisibleGenEta(muHadPFTauIsoMatchVisibleGenEta_, muHadVisibleGenEta_);
  TGraphAsymmErrors effMuMatchVsVisibleGenPT(muHadMuMatchVisibleGenPT_, muHadVisibleGenPT_);
  TGraphAsymmErrors effMuMatchVsVisibleGenEta(muHadMuMatchVisibleGenEta_, muHadVisibleGenEta_);
  TGraphAsymmErrors effPFTauMatchMuMatchVsVisibleGenPT(muHadPFTauMatchMuMatchVisibleGenPT_, muHadVisibleGenPT_);
  // everything up to here is ok
  TGraphAsymmErrors effPFTauMatchMuMatchVsVisibleGenEta(muHadPFTauMatchMuMatchVisibleGenEta_, muHadVisibleGenEta_);
  TGraphAsymmErrors effPFTauIsoMatchMuMatchVsVisibleGenPT(muHadPFTauIsoMatchMuMatchVisibleGenPT_, muHadVisibleGenPT_);
  TGraphAsymmErrors effPFTauIsoMatchMuMatchVsVisibleGenEta(muHadPFTauIsoMatchMuMatchVisibleGenEta_, muHadVisibleGenEta_);
  TGraphAsymmErrors effPFTauMatchVsGenMuPT(muHadPFTauMatchGenMuPT_, muHadGenMuPT_);
  TGraphAsymmErrors effPFTauMatchVsGenMuEta(muHadPFTauMatchGenMuEta_, muHadGenMuEta_);
  TGraphAsymmErrors effPFTauIsoMatchVsGenMuPT(muHadPFTauIsoMatchGenMuPT_, muHadGenMuPT_);
  TGraphAsymmErrors effPFTauIsoMatchVsGenMuEta(muHadPFTauIsoMatchGenMuEta_, muHadGenMuEta_);
  TGraphAsymmErrors effMuMatchVsGenMuPT(muHadMuMatchGenMuPT_, muHadGenMuPT_);
  TGraphAsymmErrors effMuMatchVsGenMuEta(muHadMuMatchGenMuEta_, muHadGenMuEta_);
  TGraphAsymmErrors effPFTauMatchMuMatchVsGenMuPT(muHadPFTauMatchMuMatchGenMuPT_, muHadGenMuPT_);
  TGraphAsymmErrors effPFTauMatchMuMatchVsGenMuEta(muHadPFTauMatchMuMatchGenMuEta_, muHadGenMuEta_);
  TGraphAsymmErrors effPFTauMatchVsGenMuHadDR(muHadPFTauMatchGenDR_, muHadGenDR_);
  TGraphAsymmErrors effPFTauIsoMatchMuMatchVsGenMuPT(muHadPFTauIsoMatchMuMatchGenMuPT_, muHadGenMuPT_);
  TGraphAsymmErrors effPFTauIsoMatchMuMatchVsGenMuEta(muHadPFTauIsoMatchMuMatchGenMuEta_, muHadGenMuEta_);
  TGraphAsymmErrors effPFTauIsoMatchVsGenMuHadDR(muHadPFTauIsoMatchGenDR_, muHadGenDR_);
  TGraphAsymmErrors effMuMatchVsGenMuHadDR(muHadMuMatchGenDR_, muHadGenDR_);
  TGraphAsymmErrors effPFTauMatchMuMatchVsGenMuHadDR(muHadPFTauMatchMuMatchGenDR_, muHadGenDR_);
  TGraphAsymmErrors effPFTauIsoMatchMuMatchVsGenMuHadDR(muHadPFTauIsoMatchMuMatchGenDR_, muHadGenDR_);

  //make efficiency canvases
  TCanvas effPFTauMatchVsVisibleGenPTCanvas("effPFTauMatchVsVisibleGenPTCanvas", "", 600, 600);
  TCanvas effPFTauMatchVsVisibleGenEtaCanvas("effPFTauMatchVsVisibleGenEtaCanvas", "", 600, 600);

  TCanvas effPFTauIsoMatchVsVisibleGenPTCanvas("effPFTauIsoMatchVsVisibleGenPTCanvas", "", 600, 600);
  TCanvas effPFTauIsoMatchVsVisibleGenEtaCanvas("effPFTauIsoMatchVsVisibleGenEtaCanvas", "", 600, 600);

  TCanvas effMuMatchVsVisibleGenPTCanvas("effMuMatchVsVisibleGenPTCanvas", "", 600, 600);
  TCanvas effMuMatchVsVisibleGenEtaCanvas("effMuMatchVsVisibleGenEtaCanvas", "", 600, 600);
  TCanvas effPFTauMatchMuMatchVsVisibleGenPTCanvas("effPFTauMatchMuMatchVsVisibleGenPTCanvas", 
						   "", 600, 600);
  TCanvas effPFTauMatchMuMatchVsVisibleGenEtaCanvas("effPFTauMatchMuMatchVsVisibleGenEtaCanvas", 
						    "", 600, 600);
  TCanvas effPFTauMatchVsGenMuPTCanvas("effPFTauMatchVsGenMuPTCanvas", "", 600, 600);
  TCanvas effPFTauMatchVsGenMuEtaCanvas("effPFTauMatchVsGenMuEtaCanvas", "", 600, 600);


  TCanvas effPFTauIsoMatchMuMatchVsVisibleGenPTCanvas("effPFTauIsoMatchMuMatchVsVisibleGenPTCanvas", 
						   "", 600, 600);
  TCanvas effPFTauIsoMatchMuMatchVsVisibleGenEtaCanvas("effPFTauIsoMatchMuMatchVsVisibleGenEtaCanvas", 
						    "", 600, 600);
  TCanvas effPFTauIsoMatchVsGenMuPTCanvas("effPFTauIsoMatchVsGenMuPTCanvas", "", 600, 600);
  TCanvas effPFTauIsoMatchVsGenMuEtaCanvas("effPFTauIsoMatchVsGenMuEtaCanvas", "", 600, 600);


  TCanvas effMuMatchVsGenMuPTCanvas("effMuMatchVsGenMuPTCanvas", "", 600, 600);
  TCanvas effMuMatchVsGenMuEtaCanvas("effMuMatchVsGenMuEtaCanvas", "", 600, 600);
  TCanvas effPFTauMatchMuMatchVsGenMuPTCanvas("effPFTauMatchMuMatchVsGenMuPTCanvas", "", 600, 600);
  TCanvas effPFTauMatchMuMatchVsGenMuEtaCanvas("effPFTauMatchMuMatchVsGenMuEtaCanvas", "", 600, 600);
  TCanvas effPFTauMatchVsGenMuHadDRCanvas("effPFTauMatchVsGenMuHadDRCanvas", "", 600, 600);


  TCanvas effPFTauIsoMatchMuMatchVsGenMuPTCanvas("effPFTauIsoMatchMuMatchVsGenMuPTCanvas", "", 600, 600);
  TCanvas effPFTauIsoMatchMuMatchVsGenMuEtaCanvas("effPFTauIsoMatchMuMatchVsGenMuEtaCanvas", "", 600, 600);
  TCanvas effPFTauIsoMatchVsGenMuHadDRCanvas("effPFTauIsoMatchVsGenMuHadDRCanvas", "", 600, 600);


  TCanvas effMuMatchVsGenMuHadDRCanvas("effMuMatchVsGenMuHadDRCanvas", "", 600, 600);
  TCanvas effPFTauMatchMuMatchVsGenMuHadDRCanvas("effPFTauMatchMuMatchVsGenMuHadDRCanvas", "", 600, 600);
  TCanvas effPFTauIsoMatchMuMatchVsGenMuHadDRCanvas("effPFTauIsoMatchMuMatchVsGenMuHadDRCanvas", "", 600, 600);

  //draw efficiency plots
  effPFTauMatchVsVisibleGenPTCanvas.cd();
  effPFTauMatchVsVisibleGenPT.Draw("AP");
  effPFTauMatchVsVisibleGenEtaCanvas.cd();
  effPFTauMatchVsVisibleGenEta.Draw("AP");
  effPFTauIsoMatchVsVisibleGenPTCanvas.cd();
  effPFTauIsoMatchVsVisibleGenPT.Draw("AP");
  effPFTauIsoMatchVsVisibleGenEtaCanvas.cd();
  effPFTauIsoMatchVsVisibleGenEta.Draw("AP");
  effMuMatchVsVisibleGenPTCanvas.cd();
  effMuMatchVsVisibleGenPT.Draw("AP");
  effMuMatchVsVisibleGenEtaCanvas.cd();
  effMuMatchVsVisibleGenEta.Draw("AP");
  effPFTauMatchMuMatchVsVisibleGenPTCanvas.cd();
  effPFTauMatchMuMatchVsVisibleGenPT.Draw("AP");
  effPFTauMatchMuMatchVsVisibleGenEtaCanvas.cd();
  effPFTauMatchMuMatchVsVisibleGenEta.Draw("AP");
  effPFTauMatchVsGenMuPTCanvas.cd();
  effPFTauMatchVsGenMuPT.Draw("AP");
  effPFTauMatchVsGenMuEtaCanvas.cd();
  effPFTauMatchVsGenMuEta.Draw("AP");
  effPFTauIsoMatchMuMatchVsVisibleGenPTCanvas.cd();
  effPFTauIsoMatchMuMatchVsVisibleGenPT.Draw("AP");
  effPFTauIsoMatchMuMatchVsVisibleGenEtaCanvas.cd();
  effPFTauIsoMatchMuMatchVsVisibleGenEta.Draw("AP");
  effPFTauIsoMatchVsGenMuPTCanvas.cd();
  effPFTauIsoMatchVsGenMuPT.Draw("AP");
  effPFTauIsoMatchVsGenMuEtaCanvas.cd();
  effPFTauIsoMatchVsGenMuEta.Draw("AP");
  effMuMatchVsGenMuPTCanvas.cd();
  effMuMatchVsGenMuPT.Draw("AP");
  effMuMatchVsGenMuEtaCanvas.cd();
  effMuMatchVsGenMuEta.Draw("AP");
  effPFTauMatchMuMatchVsGenMuPTCanvas.cd();
  effPFTauMatchMuMatchVsGenMuPT.Draw("AP");
  effPFTauMatchMuMatchVsGenMuEtaCanvas.cd();
  effPFTauMatchMuMatchVsGenMuEta.Draw("AP");
  effPFTauMatchVsGenMuHadDRCanvas.cd();
  effPFTauMatchVsGenMuHadDR.Draw("AP");
  effPFTauIsoMatchMuMatchVsGenMuPTCanvas.cd();
  effPFTauIsoMatchMuMatchVsGenMuPT.Draw("AP");
  effPFTauIsoMatchMuMatchVsGenMuEtaCanvas.cd();
  effPFTauIsoMatchMuMatchVsGenMuEta.Draw("AP");
  effPFTauIsoMatchVsGenMuHadDRCanvas.cd();
  effPFTauIsoMatchVsGenMuHadDR.Draw("AP");
  effMuMatchVsGenMuHadDRCanvas.cd();
  effMuMatchVsGenMuHadDR.Draw("AP");
  effPFTauMatchMuMatchVsGenMuHadDRCanvas.cd();
  effPFTauMatchMuMatchVsGenMuHadDR.Draw("AP");
  effPFTauIsoMatchMuMatchVsGenMuHadDRCanvas.cd();
  effPFTauIsoMatchMuMatchVsGenMuHadDR.Draw("AP");

  //write output file
  out_->cd();
  status3TauMultiplicity_->Write();
  HPSTauMultiplicity_->Write();
  HPSTausPassingEtaCutMultiplicity_->Write();
  HPSIsoTausPassingEtaCutMultiplicity_->Write();
  dRToHPSTau0_->Write();
  dRToHPSTau1_->Write();
  dRToHPSIsoTau0_->Write();
  dRToHPSIsoTau1_->Write();
  HPSTauMultiplicityVsDRToHPSTau0_->Write();
  HPSIsoTauMultiplicityVsDRToHPSIsoTau0_->Write();
//   dRLeadingGenTauToHPSTau0_->Write();
//   dRTrailingGenTauToHPSTau0_->Write();
//   dRLeadingGenTauToHPSIsoTau0_->Write();
//   dRTrailingGenTauToHPSIsoTau0_->Write();
  PFTauPT_->Write();
  HPSTauPT_->Write();
  HPSIsoTauPT_->Write();
  HPSIsoTauPTPassingEtaCut_->Write();
  genTauPT_->Write();
  dRToHPSTau0GenPTGe30_->Write();
  dRToHPSIsoTau0GenPTGe30_->Write();
  HPSTauEta_->Write();
  genTauEta_->Write();

  muHadPFTauMatchMultiplicity_->Write();
  muHadPFTauIsoMatchMultiplicity_->Write();
  muHadMuonMatchMultiplicity_->Write();
  muHadPFTauMatchMuonMatchMultiplicity_->Write();
  muHadPFTauIsoMatchMuonMatchMultiplicity_->Write();
  muHadVisibleGenPT_->Write();
  muHadPFTauMatchVisibleGenPT_->Write();
  muHadPFTauIsoMatchVisibleGenPT_->Write();
  muHadMuMatchVisibleGenPT_->Write();
  muHadPFTauMatchMuMatchVisibleGenPT_->Write();
  muHadPFTauIsoMatchMuMatchVisibleGenPT_->Write();
  muHadVisibleGenEta_->Write();
  muHadPFTauMatchVisibleGenEta_->Write();
  muHadPFTauIsoMatchVisibleGenEta_->Write();
  muHadMuMatchVisibleGenEta_->Write();
  muHadPFTauMatchMuMatchVisibleGenEta_->Write();
  muHadPFTauIsoMatchMuMatchVisibleGenEta_->Write();
  muHadGenMuPT_->Write();
  muHadPFTauMatchGenMuPT_->Write();
  muHadPFTauIsoMatchGenMuPT_->Write();
  muHadMuMatchGenMuPT_->Write();
  muHadPFTauMatchMuMatchGenMuPT_->Write();
  muHadPFTauIsoMatchMuMatchGenMuPT_->Write();
  muHadGenMuEta_->Write();
  muHadPFTauMatchGenMuEta_->Write();
  muHadPFTauIsoMatchGenMuEta_->Write();
  muHadMuMatchGenMuEta_->Write();
  muHadPFTauMatchMuMatchGenMuEta_->Write();
  muHadPFTauIsoMatchMuMatchGenMuEta_->Write();
  muHadGenDR_->Write();
  muHadPFTauMatchGenDR_->Write();
  muHadMuMatchGenDR_->Write();
  muHadPFTauMatchMuMatchGenDR_->Write();
  muHadPFTauIsoMatchMuMatchGenDR_->Write();
  effPFTauMatchVsVisibleGenPT.Write();
  effPFTauIsoMatchVsVisibleGenPT.Write();
  effPFTauMatchVsVisibleGenEta.Write();
  effPFTauIsoMatchVsVisibleGenEta.Write();
  effMuMatchVsVisibleGenPT.Write();
  effMuMatchVsVisibleGenEta.Write();
  effPFTauMatchMuMatchVsVisibleGenPT.Write();
  effPFTauMatchMuMatchVsVisibleGenEta.Write();
  effPFTauMatchVsGenMuPT.Write();
  effPFTauMatchVsGenMuEta.Write();
  effPFTauIsoMatchMuMatchVsVisibleGenPT.Write();
  effPFTauIsoMatchMuMatchVsVisibleGenEta.Write();
  effPFTauIsoMatchVsGenMuPT.Write();
  effPFTauIsoMatchVsGenMuEta.Write();
  effMuMatchVsGenMuPT.Write();
  effMuMatchVsGenMuEta.Write();
  effPFTauMatchMuMatchVsGenMuPT.Write();
  effPFTauMatchMuMatchVsGenMuEta.Write();
  effPFTauMatchVsGenMuHadDR.Write();
  effPFTauIsoMatchMuMatchVsGenMuPT.Write();
  effPFTauIsoMatchMuMatchVsGenMuEta.Write();
  effPFTauIsoMatchVsGenMuHadDR.Write();
  effMuMatchVsGenMuHadDR.Write();
  effPFTauMatchMuMatchVsGenMuHadDR.Write();
  effPFTauMatchVsVisibleGenPTCanvas.Write();
  effPFTauMatchVsVisibleGenEtaCanvas.Write();
  effPFTauIsoMatchMuMatchVsGenMuHadDR.Write();
  effPFTauIsoMatchVsVisibleGenPTCanvas.Write();
  effPFTauIsoMatchVsVisibleGenEtaCanvas.Write();
  effMuMatchVsVisibleGenPTCanvas.Write();
  effMuMatchVsVisibleGenEtaCanvas.Write();
  effPFTauMatchMuMatchVsVisibleGenPTCanvas.Write();
  effPFTauMatchMuMatchVsVisibleGenEtaCanvas.Write();
  effPFTauMatchVsGenMuPTCanvas.Write();
  effPFTauMatchVsGenMuEtaCanvas.Write();
  effPFTauIsoMatchMuMatchVsVisibleGenPTCanvas.Write();
  effPFTauIsoMatchMuMatchVsVisibleGenEtaCanvas.Write();
  effPFTauIsoMatchVsGenMuPTCanvas.Write();
  effPFTauIsoMatchVsGenMuEtaCanvas.Write();
  effMuMatchVsGenMuPTCanvas.Write();
  effMuMatchVsGenMuEtaCanvas.Write();
  effPFTauMatchMuMatchVsGenMuPTCanvas.Write();
  effPFTauMatchMuMatchVsGenMuEtaCanvas.Write();
  effPFTauMatchVsGenMuHadDRCanvas.Write();
  effPFTauIsoMatchMuMatchVsGenMuPTCanvas.Write();
  effPFTauIsoMatchMuMatchVsGenMuEtaCanvas.Write();
  effPFTauIsoMatchVsGenMuHadDRCanvas.Write();
  effMuMatchVsGenMuHadDRCanvas.Write();
  effPFTauMatchMuMatchVsGenMuHadDRCanvas.Write();
  effPFTauIsoMatchMuMatchVsGenMuHadDRCanvas.Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
void 
TauAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TauAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TauAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TauAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void TauAnalyzer::reset()
{
  if (out_ != NULL) delete out_;
  out_ = NULL;
  outFileName_.clear();
  HPSDiscriminatorTags_.clear();
  HPSDiscriminators_.clear();
  momPDGID_ = 0;
  genMuTauPTMin_ = 0.0;
  genMuPTMin_ = 0.0;
  effVsEtaPTMin_ = 0.0;
  cfg_ = NULL;
  if (status3TauMultiplicity_ != NULL) delete status3TauMultiplicity_;
  status3TauMultiplicity_ = NULL;
  if (HPSTauMultiplicity_ != NULL) delete HPSTauMultiplicity_;
  HPSTauMultiplicity_ = NULL;
  if (HPSTausPassingEtaCutMultiplicity_ != NULL) delete HPSTausPassingEtaCutMultiplicity_;
  HPSTausPassingEtaCutMultiplicity_ = NULL;
  if (HPSIsoTausPassingEtaCutMultiplicity_ != NULL) delete HPSIsoTausPassingEtaCutMultiplicity_;
  HPSIsoTausPassingEtaCutMultiplicity_ = NULL;
  if (dRToHPSTau0_ != NULL) delete dRToHPSTau0_;
  dRToHPSTau0_ = NULL;
  if (dRToHPSTau1_ != NULL) delete dRToHPSTau1_;
  dRToHPSTau1_ = NULL;
  if (dRToHPSIsoTau0_ != NULL) delete dRToHPSIsoTau0_;
  dRToHPSIsoTau0_ = NULL;
  if (dRToHPSIsoTau1_ != NULL) delete dRToHPSIsoTau1_;
  dRToHPSIsoTau1_ = NULL;
  if (HPSTauMultiplicityVsDRToHPSTau0_ != NULL) delete HPSTauMultiplicityVsDRToHPSTau0_;
  HPSTauMultiplicityVsDRToHPSTau0_ = NULL;
  if (HPSIsoTauMultiplicityVsDRToHPSIsoTau0_ != NULL) {
    delete HPSIsoTauMultiplicityVsDRToHPSIsoTau0_;
  }
  HPSIsoTauMultiplicityVsDRToHPSIsoTau0_ = NULL;
  if (dRLeadingGenTauToHPSTau0_ != NULL) delete dRLeadingGenTauToHPSTau0_;
  dRLeadingGenTauToHPSTau0_ = NULL;
  if (dRTrailingGenTauToHPSTau0_ != NULL) delete dRTrailingGenTauToHPSTau0_;
  dRTrailingGenTauToHPSTau0_ = NULL;
  if (dRLeadingGenTauToHPSIsoTau0_ != NULL) delete dRLeadingGenTauToHPSIsoTau0_;
  dRLeadingGenTauToHPSIsoTau0_ = NULL;
  if (dRTrailingGenTauToHPSIsoTau0_ != NULL) delete dRTrailingGenTauToHPSIsoTau0_;
  dRTrailingGenTauToHPSIsoTau0_ = NULL;
  if (PFTauPT_ != NULL) delete PFTauPT_;
  PFTauPT_ = NULL;
  if (HPSTauPT_ != NULL) delete HPSTauPT_;
  HPSTauPT_ = NULL;
  if (HPSIsoTauPT_ != NULL) delete HPSIsoTauPT_;
  HPSIsoTauPT_ = NULL;
  if (HPSIsoTauPTPassingEtaCut_ != NULL) delete HPSIsoTauPTPassingEtaCut_;
  HPSIsoTauPTPassingEtaCut_ = NULL;
  if (genTauPT_ != NULL) delete genTauPT_;
  genTauPT_ = NULL;
  if (dRToHPSTau0GenPTGe30_ != NULL) delete dRToHPSTau0GenPTGe30_;
  dRToHPSTau0GenPTGe30_ = NULL;
  if (dRToHPSIsoTau0GenPTGe30_ != NULL) delete dRToHPSIsoTau0GenPTGe30_;
  dRToHPSIsoTau0GenPTGe30_ = NULL;
  if (HPSTauEta_ != NULL) delete HPSTauEta_;
  HPSTauEta_ = NULL;
  if (genTauEta_ != NULL) delete genTauEta_;
  genTauEta_ = NULL;
  //oof
  if (muHadPFTauMatchMultiplicity_ != NULL) delete muHadPFTauMatchMultiplicity_;
  muHadPFTauMatchMultiplicity_ = NULL;
 if (muHadPFTauIsoMatchMultiplicity_ != NULL) delete muHadPFTauIsoMatchMultiplicity_;
  muHadPFTauIsoMatchMultiplicity_ = NULL;
  if (muHadMuonMatchMultiplicity_ != NULL) delete muHadMuonMatchMultiplicity_;
  muHadMuonMatchMultiplicity_ = NULL;
  if (muHadPFTauMatchMuonMatchMultiplicity_ != NULL) delete muHadPFTauMatchMuonMatchMultiplicity_;
  muHadPFTauMatchMuonMatchMultiplicity_ = NULL;
  if (muHadPFTauIsoMatchMuonMatchMultiplicity_ != NULL) delete muHadPFTauIsoMatchMuonMatchMultiplicity_;
  muHadPFTauIsoMatchMuonMatchMultiplicity_ = NULL;
  if (muHadVisibleGenPT_ != NULL) delete muHadVisibleGenPT_;
  muHadVisibleGenPT_ = NULL;
  if (muHadPFTauMatchVisibleGenPT_ != NULL) delete muHadPFTauMatchVisibleGenPT_;
  muHadPFTauMatchVisibleGenPT_ = NULL;
 if (muHadPFTauIsoMatchVisibleGenPT_ != NULL) delete muHadPFTauIsoMatchVisibleGenPT_;
  muHadPFTauIsoMatchVisibleGenPT_ = NULL;
  if (muHadMuMatchVisibleGenPT_ != NULL) delete muHadMuMatchVisibleGenPT_;
  muHadMuMatchVisibleGenPT_ = NULL;
  if (muHadPFTauMatchMuMatchVisibleGenPT_ != NULL) delete muHadPFTauMatchMuMatchVisibleGenPT_;
  muHadPFTauMatchMuMatchVisibleGenPT_ = NULL;
  if (muHadPFTauIsoMatchMuMatchVisibleGenPT_ != NULL) delete muHadPFTauIsoMatchMuMatchVisibleGenPT_;
  muHadPFTauIsoMatchMuMatchVisibleGenPT_ = NULL;
  if (muHadVisibleGenEta_ != NULL) delete muHadVisibleGenEta_;
  muHadVisibleGenEta_ = NULL;
  if (muHadPFTauMatchVisibleGenEta_ != NULL) delete muHadPFTauMatchVisibleGenEta_;
  muHadPFTauMatchVisibleGenEta_ = NULL;
  if (muHadPFTauIsoMatchVisibleGenEta_ != NULL) delete muHadPFTauIsoMatchVisibleGenEta_;
  muHadPFTauIsoMatchVisibleGenEta_ = NULL;
  if (muHadMuMatchVisibleGenEta_ != NULL) delete muHadMuMatchVisibleGenEta_;
  muHadMuMatchVisibleGenEta_ = NULL;
  if (muHadPFTauMatchMuMatchVisibleGenEta_ != NULL) delete muHadPFTauMatchMuMatchVisibleGenEta_;
  muHadPFTauMatchMuMatchVisibleGenEta_ = NULL;
  if (muHadPFTauIsoMatchMuMatchVisibleGenEta_ != NULL) delete muHadPFTauIsoMatchMuMatchVisibleGenEta_;
  muHadPFTauIsoMatchMuMatchVisibleGenEta_ = NULL;
  if (muHadGenMuPT_ != NULL) delete muHadGenMuPT_;
  muHadGenMuPT_ = NULL;
  if (muHadPFTauMatchGenMuPT_ != NULL) delete muHadPFTauMatchGenMuPT_;
  muHadPFTauMatchGenMuPT_ = NULL;
  if (muHadPFTauIsoMatchGenMuPT_ != NULL) delete muHadPFTauIsoMatchGenMuPT_;
  muHadPFTauIsoMatchGenMuPT_ = NULL;
  if (muHadMuMatchGenMuPT_ != NULL) delete muHadMuMatchGenMuPT_;
  muHadMuMatchGenMuPT_ = NULL;
  if (muHadPFTauMatchMuMatchGenMuPT_ != NULL) delete muHadPFTauMatchMuMatchGenMuPT_;
  muHadPFTauMatchMuMatchGenMuPT_ = NULL;
  if (muHadPFTauIsoMatchMuMatchGenMuPT_ != NULL) delete muHadPFTauIsoMatchMuMatchGenMuPT_;
  muHadPFTauIsoMatchMuMatchGenMuPT_ = NULL;
  if (muHadGenMuEta_ != NULL) delete muHadGenMuEta_;
  muHadGenMuEta_ = NULL;
  if (muHadPFTauMatchGenMuEta_ != NULL) delete muHadPFTauMatchGenMuEta_;
  muHadPFTauMatchGenMuEta_ = NULL;
  if (muHadPFTauIsoMatchGenMuEta_ != NULL) delete muHadPFTauIsoMatchGenMuEta_;
  muHadPFTauIsoMatchGenMuEta_ = NULL;
  if (muHadMuMatchGenMuEta_ != NULL) delete muHadMuMatchGenMuEta_;
  muHadMuMatchGenMuEta_ = NULL;
  if (muHadPFTauMatchMuMatchGenMuEta_ != NULL) delete muHadPFTauMatchMuMatchGenMuEta_;
  muHadPFTauMatchMuMatchGenMuEta_ = NULL;
  if (muHadPFTauIsoMatchMuMatchGenMuEta_ != NULL) delete muHadPFTauIsoMatchMuMatchGenMuEta_;
  muHadPFTauIsoMatchMuMatchGenMuEta_ = NULL;
  if (muHadGenDR_ != NULL) delete muHadGenDR_;
  muHadGenDR_ = NULL;
  if (muHadPFTauMatchGenDR_ != NULL) delete muHadPFTauMatchGenDR_;
  muHadPFTauMatchGenDR_ = NULL;
  if (muHadPFTauIsoMatchGenDR_ != NULL) delete muHadPFTauIsoMatchGenDR_;
  muHadPFTauIsoMatchGenDR_ = NULL;
  if (muHadMuMatchGenDR_ != NULL) delete muHadMuMatchGenDR_;
  muHadMuMatchGenDR_ = NULL;
  if (muHadPFTauMatchMuMatchGenDR_ != NULL) delete muHadPFTauMatchMuMatchGenDR_;
  muHadPFTauMatchMuMatchGenDR_ = NULL;
  if (muHadPFTauIsoMatchMuMatchGenDR_ != NULL) delete muHadPFTauIsoMatchMuMatchGenDR_;
  muHadPFTauIsoMatchMuMatchGenDR_ = NULL;
}

double TauAnalyzer::dR(const std::vector<double>& dR2, const unsigned int pos) const
{
  double dR = -1.0;
  if (dR2.size() == 0) dR = 100.0;
  else if (pos >= dR2.size()) dR = sqrt(dR2[dR2.size() - 1]);
  else dR = sqrt(dR2[pos]);
  return dR;
}

void TauAnalyzer::fillDRHistogramNearFar(const unsigned int pos, const bool iso, 
					 const std::vector<double>& dR2) const
{
  double deltaR = dR(dR2, pos);
  if (!iso) {
    if (pos == 0) {
      dRToHPSTau0_->Fill(deltaR);
      HPSTauMultiplicityVsDRToHPSTau0_->Fill(deltaR, dR2.size());
    }
    else if (pos == 1) dRToHPSTau1_->Fill(deltaR);
    else {
      std::stringstream err;
      err << "Unsupported position " << pos << " in void TauAnalyzer::fillDRHistogram(const ";
      err << "unsigned int pos, const unsigned int iso, const std::vector<double> dR2) const.\n";
      edm::LogError("TauAnalyzer") << err.str();
    }
  }
  else {
    if (pos == 0) {
      dRToHPSIsoTau0_->Fill(deltaR);
      HPSIsoTauMultiplicityVsDRToHPSIsoTau0_->Fill(deltaR, dR2.size());
    }
    else if (pos == 1) dRToHPSIsoTau1_->Fill(deltaR);
    else {
      std::stringstream err;
      err << "Unsupported position " << pos << " in void TauAnalyzer::fillDRHistogram(const ";
      err << "unsigned int pos, const unsigned int iso, const std::vector<double> dR2) const.\n";
      edm::LogError("TauAnalyzer") << err.str();
    }
  }
}

void TauAnalyzer::fillDRHistogramLeadTrail(const bool iso, const unsigned int lead, 
					   const std::vector<double>& dR2) const
{
  double deltaR = dR(dR2);
  if (!iso) {
    if (lead == 1) dRLeadingGenTauToHPSTau0_->Fill(deltaR);
    else if (lead == 2) dRTrailingGenTauToHPSTau0_->Fill(deltaR);
  }
  else {
    if (lead == 1) dRLeadingGenTauToHPSIsoTau0_->Fill(deltaR);
    else if (lead == 2) dRTrailingGenTauToHPSIsoTau0_->Fill(deltaR);
  }
}

void TauAnalyzer::fillDRHistogramHighLowPT(const bool iso, const std::vector<double>& dR2) const
{
  double deltaR = dR(dR2);
  if (!iso) dRToHPSTau0GenPTGe30_->Fill(deltaR);
  else dRToHPSIsoTau0GenPTGe30_->Fill(deltaR);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TauAnalyzer);
