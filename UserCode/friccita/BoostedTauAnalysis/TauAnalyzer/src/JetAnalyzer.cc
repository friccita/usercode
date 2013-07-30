// -*- C++ -*-
//
// Package:    TauAnalyzer
// Class:      JetAnalyzer
// 
/**\class JetAnalyzer JetAnalyzer.cc 
   BoostedTauAnalysis/TauAnalyzer/src/JetAnalyzer.cc

   Description: analyze jet activity related to taus

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 16:40:51 CEST 2012
// $Id: JetAnalyzer.cc,v 1.1 2012/09/19 10:57:14 yohay Exp $
//
//


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
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "BoostedTauAnalysis/Common/interface/GenTauDecayID.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

//
// class declaration
//

class JetAnalyzer : public edm::EDAnalyzer {
public:
  explicit JetAnalyzer(const edm::ParameterSet&);
  ~JetAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  //delete memory
  void reset(const bool);

  //fill pT histogram with object of arbitrary type
  template<typename U>
  void fillPTHistogramArbitrary(edm::Handle<edm::View<U> >& pView, TH1F* hist)
  {
    for (unsigned int i = 0; i < pView->size(); ++i) hist->Fill(pView->refAt(i)->pt());
  }

  //fill ET histogram
  template<typename U>
  void fillETHistogram(edm::Handle<edm::View<U> >& pView, TH1F* hist)
  {
    for (unsigned int i = 0; i < pView->size(); ++i) hist->Fill(pView->refAt(i)->et());
  }

  //fill HT histogram
  template<typename U>
  void fillHTHistogram(std::map<TH1F*, edm::Handle<edm::View<U> > >& map, 
		       const unsigned int numObjs)
  {
    double HT = 0.0;
    unsigned int nObjs = 0;
    for (typename std::map<TH1F*, edm::Handle<edm::View<U> > >::iterator i = map.begin(); 
	 i != map.end(); ++i) {
      nObjs+=i->second->size();
      for (unsigned int j = 0; j < i->second->size(); ++j) HT+=(i->second->refAt(j)->et());
    }
    if (nObjs == 2) HT_->Fill(HT);
  }

  //make visible gen tau pT canvas for 1 decay mode, multiple pT ranks
  void makePTRankCanvas(TCanvas&, TLegend&, const std::string&, std::vector<TH1F*>&);

  //format and draw multiple pT histograms on one canvas
  void drawMultiplePTHistograms(TCanvas&, std::vector<TH1F*>&, const std::vector<unsigned int>&, 
				const std::vector<unsigned int>&, TLegend&, 
				const std::vector<std::string>&, const std::string&);

  // ----------member data ---------------------------

  //pointer to output file object
  TFile* out_;

  //name of output file
  std::string outFileName_;

  //vector of jet tags
  std::vector<edm::InputTag> jetTags_;

  //MET tag
  edm::InputTag METTag_;

  //marker colors for histograms with different pT rank
  std::vector<unsigned int> pTRankColors_;

  //marker styles for histograms with different pT rank
  std::vector<unsigned int> pTRankStyles_;

  //legend entries for histograms with different pT rank
  std::vector<std::string> pTRankEntries_;

  //histogram of eta for the single collection
  TH1F* eta_;

  //vector of jet pT histograms, highest pT rank first
  std::vector<TH1F*> jetPTHists_;

  //histogram of MET
  TH1F* MET_;

  //histogram of HT
  TH1F* HT_;
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
JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  jetTags_(iConfig.getParameter<std::vector<edm::InputTag> >("jetTags")),
  METTag_(iConfig.getParameter<edm::InputTag>("METTag")),
  pTRankColors_(iConfig.getParameter<std::vector<unsigned int> >("pTRankColors")),
  pTRankStyles_(iConfig.getParameter<std::vector<unsigned int> >("pTRankStyles")),
  pTRankEntries_(iConfig.getParameter<std::vector<std::string> >("pTRankEntries"))
{
  //now do what ever initialization is needed
  reset(false);
}

JetAnalyzer::~JetAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void JetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get jet tags
  std::map<TH1F*, edm::Handle<edm::View<reco::PFJet> > > jetCollMap;
  for (std::vector<edm::InputTag>::const_iterator iTag = jetTags_.begin(); 
       iTag != jetTags_.end(); ++iTag) {
    edm::Handle<edm::View<reco::PFJet> > pView;
    iEvent.getByLabel(*iTag, pView);
    jetCollMap[jetPTHists_[iTag - jetTags_.begin()]] = pView;
  }

  //get MET tag
  edm::Handle<edm::View<reco::PFMET> > pMET;
  iEvent.getByLabel(METTag_, pMET);

  //fill pT histograms for jets, 1 per pT rank
  for (std::map<TH1F*, edm::Handle<edm::View<reco::PFJet> > >::iterator i = jetCollMap.begin(); 
       i != jetCollMap.end(); ++i) { fillPTHistogramArbitrary(i->second, i->first); }

  //plot MET distribution
  fillETHistogram(pMET, MET_);

  /*plot HT distribution (sum ET of jets matched to di-tau objects in the case where there are 2 
    di-tau objects with jet matches)*/
  fillHTHistogram(jetCollMap, 2);
}


// ------------ method called once each job just before starting event loop  ------------
void JetAnalyzer::beginJob()
{
  //open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book jet pT histograms split by pT rank
  for (std::vector<edm::InputTag>::const_iterator iTag = jetTags_.begin(); 
       iTag != jetTags_.end(); ++iTag) {
    jetPTHists_.push_back(new TH1F((iTag->label() + "_" + iTag->instance() + "_pT")
				   .c_str(), "", 20, 0.0, 100.0));
  }

  //book reco MET histogram
  MET_ = new TH1F("MET", "", 20, 0.0, 100.0);

  //book reco HT histogram
  HT_ = new TH1F("HT", "", 100, 0.0, 500.0);
}

// ------------ method called once each job just after ending the event loop  ------------
void JetAnalyzer::endJob() 
{
  //make the jet pT canvas
  TCanvas jetPTRankCanvas("jetPTRankCanvas", "", 600, 600);
  TLegend jetPTRankLegend(0.4, 0.6, 0.8, 0.8);
  makePTRankCanvas(jetPTRankCanvas, jetPTRankLegend, 
		   "gg fusion NMSSM Higgs-matched AK5 jets", jetPTHists_);

  //write output file
  out_->cd();
  for (std::vector<TH1F*>::iterator iHist = jetPTHists_.begin(); 
       iHist != jetPTHists_.end(); ++iHist) { (*iHist)->Write(); }
  MET_->Write();
  HT_->Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
void JetAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void JetAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void JetAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, 
				       edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void JetAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, 
				     edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void JetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void JetAnalyzer::makePTRankCanvas(TCanvas& canvas, TLegend& legend, 
				   const std::string& header, 
				   std::vector<TH1F*>& hists)
{
  drawMultiplePTHistograms(canvas, hists, pTRankColors_, pTRankStyles_, legend, 
			   pTRankEntries_, header);
}

void JetAnalyzer::drawMultiplePTHistograms(TCanvas& canvas, 
					   std::vector<TH1F*>& hists, 
					   const std::vector<unsigned int>& colors, 
					   const std::vector<unsigned int>& styles, 
					   TLegend& legend, 
					   const std::vector<std::string>& entries, 
					   const std::string& header)
{
  Common::setLegendOptions(legend, header.c_str());
  canvas.cd();
  Common::setCanvasOptions(canvas, 1, 0, 0);
  TH1F* pHistWithMaxMaxBin = NULL;
  Double_t maxBinContent = 0.0;
  for (std::vector<TH1F*>::iterator iHist = hists.begin(); 
       iHist != hists.end(); ++iHist) {
    const unsigned int i = iHist - hists.begin();
    Common::setHistogramOptions(*iHist, colors[i], 0.7, styles[i], 1.0, 
				"p_{T} (GeV)", "", 0.04);
    legend.AddEntry(*iHist, entries[i].c_str(), "l");
    Int_t histMaxBin = (*iHist)->GetMaximumBin();
    Double_t histMaxBinContent = (*iHist)->GetBinContent(histMaxBin);
    if (histMaxBinContent > maxBinContent) {
      maxBinContent = histMaxBinContent;
      pHistWithMaxMaxBin = *iHist;
    }
  }
  pHistWithMaxMaxBin->Draw();
  for (std::vector<TH1F*>::iterator iHist = hists.begin(); 
       iHist != hists.end(); ++iHist) {
    (*iHist)->Draw("SAME");
  }
  legend.Draw();
  canvas.Write();
}

void JetAnalyzer::reset(const bool doDelete)
{
  if ((doDelete) && (out_ != NULL)) delete out_;
  out_ = NULL;
  for (std::vector<TH1F*>::iterator iHist = jetPTHists_.begin(); 
       iHist != jetPTHists_.end(); ++iHist) {
    if ((doDelete) && (*iHist != NULL)) delete *iHist;
    *iHist = NULL;
  }
  if ((doDelete) && (MET_ != NULL)) delete MET_;
  MET_ = NULL;
  if ((doDelete) && (HT_ != NULL)) delete HT_;
  HT_ = NULL;
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetAnalyzer);
