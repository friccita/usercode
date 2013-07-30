// -*- C++ -*-
//
// Package:    JetChecker
// Class:      JetChecker
// 
/**\class JetChecker JetChecker.cc BoostedTauAnalysis/JetChecker/src/JetChecker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesca Ricci-Tam,6 R-025,+41227672274,
//         Created:  Mon Sep 10 10:28:20 CEST 2012
// $Id: JetChecker.cc,v 1.1 2013/05/06 09:18:44 friccita Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>

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
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

using namespace edm;
using namespace reco;
using namespace std;

//
// class declaration
//

class JetChecker : public edm::EDAnalyzer {
   public:
      explicit JetChecker(const edm::ParameterSet&);
      ~JetChecker();

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
  std::string outFileName_;
  edm::InputTag JetCollectionBefore_;
  edm::InputTag JetCollectionAfter_;
  TH1F *pt_difference;
  TH1F *eta_difference;
  TH1F *phi_difference;
  TH1I *trackno_difference;
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
JetChecker::JetChecker(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  outFileName_ = iConfig.getParameter<std::string>("outFileName");
  JetCollectionBefore_ = iConfig.getParameter<edm::InputTag>("JetCollectionBefore");
  JetCollectionAfter_ = iConfig.getParameter<edm::InputTag>("JetCollectionAfter");

}


JetChecker::~JetChecker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
JetChecker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  Handle< std::vector<reco::PFJet> > jetBeforeHandle;
  iEvent.getByLabel(JetCollectionBefore_, jetBeforeHandle);
  vector<reco::PFJet> const * jets_before = jetBeforeHandle.product();

  Handle< std::vector<reco::PFJet> > jetAfterHandle;
  iEvent.getByLabel(JetCollectionAfter_, jetAfterHandle);
  vector<reco::PFJet> const * jets_after = jetAfterHandle.product();

  cout << "Size of original PFJet collection = " << jets_before->size() << endl;
  cout << "Size of post-PF2PAT PFJet collection = " << jets_after->size() << endl;

  double mindelR = 100.;
  double delR = 100.;
  int jPosition = 0;
  bool shouldIgnore;

  std::vector<unsigned int> ignorePositions;
  for (unsigned int i = 0; i < jets_before->size(); i++)
    { // loop over original PFJets
      for (unsigned int j = 0; j < jets_after->size(); j++)
	{ // loop over PF2PAT PFJets
	  shouldIgnore = false;
	  mindelR = 100.;
	  if (ignorePositions.size() != 0)
	    { // if there are things to ignore
	      for (unsigned int k = 0; k < ignorePositions.size(); k++)
		{
		  if (j == ignorePositions.at(k))
		    {
		      shouldIgnore = true;
		      break;
		    }
		}
	    } // if there are things to ignore
	  if (shouldIgnore == true)
	    break;
	  else
	    {
	      delR = deltaR(jets_before->at(i).eta(), jets_before->at(i).phi(), jets_after->at(j).eta(), jets_after->at(j).phi());
	      if (delR < mindelR)
		{
		  mindelR = delR;
		  jPosition = j;
		}
	    }
	} // loop over PF2PAT PFJets
      if (mindelR < 0.3)
	{ // if matched
	  ignorePositions.push_back(jPosition);
	  double ptdiff = jets_before->at(i).pt() - jets_after->at(jPosition).pt();
	  double etadiff = fabs(jets_before->at(i).eta() - jets_after->at(jPosition).eta());
	  double phidiff = fabs(jets_before->at(i).phi() - jets_after->at(jPosition).phi());
	  if (phidiff > 3.14159265)
	    phidiff -= 3.14159265;
	  pt_difference->Fill(ptdiff);
	  eta_difference->Fill(etadiff);
	  phi_difference->Fill(phidiff);

	  reco::TrackRefVector tracksbefore = jets_before->at(i).getTrackRefs();
	  reco::TrackRefVector tracksafter = jets_after->at(jPosition).getTrackRefs();
	  int trackno_before = tracksbefore.size();
	  int trackno_after = tracksafter.size();
	  cout << "Number of tracks in original PFJet = " << trackno_before << endl;
	  cout << "Number of tracks after PF2PAT procedure = " << trackno_after << endl;
	  trackno_difference->Fill(trackno_before - trackno_after);;
	  
	} // if matched
      
    } // loop over original PFJets

}


// ------------ method called once each job just before starting event loop  ------------
void 
JetChecker::beginJob()
{
  out_ = new TFile(outFileName_.c_str(), "RECREATE");
  pt_difference = new TH1F("pt_difference", "Difference in pT between matched jets", 201, -100., 100.);
  eta_difference = new TH1F("eta_difference", "Difference in eta between matched jets", 80, -4., 4.);
  phi_difference = new TH1F("phi_difference", "Difference in phi between matched jets", 24, 0., 6.);
  trackno_difference = new TH1I("trackno_difference", "Difference in the number of tracks between matched jets", 201, -100, 100);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetChecker::endJob() 
{
  out_->cd();
  pt_difference->Write();
  eta_difference->Write();
  phi_difference->Write();
  trackno_difference->Write();
  out_->Write();
  out_->Close();
  
}

// ------------ method called when starting to processes a run  ------------
void 
JetChecker::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
JetChecker::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
JetChecker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
JetChecker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetChecker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetChecker);
