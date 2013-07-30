// -*- C++ -*-
//
// Package:    GenMatchedRecoObjectProducer
// Class:      CustomTauSelector
// 
/**\class CustomTauSelector CustomTauSelector.cc 
   BoostedTauAnalysis/GenMatchedRecoObjectProducer/src/CustomTauSelector.cc

 Description: create a collection of custom selected hadronic taus to put in the event, and stop 
 processing if no taus are selected

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Fri Aug 24 17:10:12 CEST 2012
// $Id: CustomTauSelector.cc,v 1.1 2012/09/19 09:42:49 yohay Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"

//
// class declaration
//

class CustomTauSelector : public edm::EDFilter {
public:
  explicit CustomTauSelector(const edm::ParameterSet&);
  ~CustomTauSelector();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  virtual bool beginRun(edm::Run&, edm::EventSetup const&);
  virtual bool endRun(edm::Run&, edm::EventSetup const&);
  virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  // ----------member data ---------------------------

  //input tag for reco tau collection
  edm::InputTag tauTag_;

  //vector of input tags, 1 for each discriminator the tau should pass
  std::vector<edm::InputTag> tauDiscriminatorTags_;

  //|eta| cut
  double etaMax_;

  //minimum number of objects that must be found to pass the filter
  unsigned int minNumObjsToPassFilter_;
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
CustomTauSelector::CustomTauSelector(const edm::ParameterSet& iConfig) :
  tauTag_(iConfig.getParameter<edm::InputTag>("tauTag")),
  tauDiscriminatorTags_(iConfig.getParameter<std::vector<edm::InputTag> >("tauDiscriminatorTags")),
  etaMax_(iConfig.getParameter<double>("etaMax")),
  minNumObjsToPassFilter_(iConfig.getParameter<unsigned int>("minNumObjsToPassFilter"))
{
  produces<reco::PFTauCollection>();
}


CustomTauSelector::~CustomTauSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool CustomTauSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //create pointer to output collection
  std::auto_ptr<reco::PFTauCollection> tauColl(new reco::PFTauCollection);

  //get taus
  edm::Handle<reco::PFTauCollection> pTaus;
  iEvent.getByLabel(tauTag_, pTaus);

  //get tau discriminators
  std::vector<edm::Handle<reco::PFTauDiscriminator> > 
    pTauDiscriminators(tauDiscriminatorTags_.size(), edm::Handle<reco::PFTauDiscriminator>());
  for (std::vector<edm::InputTag>::const_iterator iTag = tauDiscriminatorTags_.begin(); 
       iTag != tauDiscriminatorTags_.end(); ++iTag) {
    iEvent.getByLabel(*iTag, pTauDiscriminators[iTag - tauDiscriminatorTags_.begin()]);
  }

  //fill STL container with taus passing specified discriminators in specified eta range
  std::vector<reco::PFTau*> taus = Common::getRecoTaus(pTaus, pTauDiscriminators, etaMax_);

  //fill output collection
  for (std::vector<reco::PFTau*>::const_iterator iTau = taus.begin(); iTau != taus.end(); 
       ++iTau) { tauColl->push_back(**iTau); }
  iEvent.put(tauColl);

  //if not enough taus passing cuts were found in this event, stop processing
  return (taus.size() >= minNumObjsToPassFilter_);
}

// ------------ method called once each job just before starting event loop  ------------
void 
CustomTauSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CustomTauSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
CustomTauSelector::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
CustomTauSelector::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
CustomTauSelector::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
CustomTauSelector::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CustomTauSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(CustomTauSelector);
