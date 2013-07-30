#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "TH1F.h"
#include "TLegend.h"

void minDR::setCandidate(const reco::Candidate* cand)
{
  cand_ = const_cast<reco::Candidate*>(cand);
}

void minDR::deleteCandidate() { cand_ = NULL; }

reco::Candidate* minDR::getCandidate() const { return cand_; }

void Common::sortByPT(std::vector<reco::Candidate*>& objects)
{
  std::sort(objects.begin(), objects.end(), compareCandidatePT);
}

void Common::sortByPT(std::vector<GenTauDecayID>& objects)
{
  std::sort(objects.begin(), objects.end(), compareGenTauDecayIDPT);
}

void Common::sortByPT(std::vector<reco::PFJet*>& objects)
{
  std::sort(objects.begin(), objects.end(), comparePFJetPT);
}

void Common::sortByPT(std::vector<reco::PFTau*>& objects)
{
  std::sort(objects.begin(), objects.end(), comparePFTauPT);
}

void Common::sortByPT(std::vector<reco::Muon*>& objects)
{
  std::sort(objects.begin(), objects.end(), compareMuonPT);
}

void Common::sortByPT(std::vector<reco::MuonRef>& objects)
{
  std::sort(objects.begin(), objects.end(), compareMuonRefPT);
}

void Common::sortByPT(std::vector<reco::PFJetRef>& objects)
{
  std::sort(objects.begin(), objects.end(), comparePFJetRefPT);
}

void Common::sortByPT(std::vector<reco::PFTauRef>& objects)
{
  std::sort(objects.begin(), objects.end(), comparePFTauRefPT);
}

bool Common::compareCandidatePT(reco::Candidate* object1, reco::Candidate* object2)
{
  return (object1->pt() < object2->pt());
}

bool Common::compareGenTauDecayIDPT(GenTauDecayID object1, GenTauDecayID object2)
{
  return (object1.getVisibleTauP4().Pt() < object2.getVisibleTauP4().Pt());
}

bool Common::comparePFJetPT(reco::PFJet* object1, reco::PFJet* object2)
{
  return (object1->pt() < object2->pt());
}

bool Common::comparePFTauPT(reco::PFTau* object1, reco::PFTau* object2)
{
  return (object1->pt() < object2->pt());
}

bool Common::compareMuonPT(reco::Muon* object1, reco::Muon* object2)
{
  return (object1->pt() < object2->pt());
}

bool Common::compareMuonRefPT(reco::MuonRef object1, reco::MuonRef object2)
{
  return (object1->pt() < object2->pt());
}

bool Common::comparePFJetRefPT(reco::PFJetRef object1, reco::PFJetRef object2)
{
  return (object1->pt() < object2->pt());
}

bool Common::comparePFTauRefPT(reco::PFTauRef object1, reco::PFTauRef object2)
{
  return (object1->pt() < object2->pt());
}

reco::Vertex* Common::getPrimaryVertex(edm::Handle<reco::VertexCollection>& pVertices)
{
  reco::VertexCollection::const_iterator iVtx = pVertices->begin();
  reco::Vertex* pPV = NULL;
  while ((iVtx != pVertices->end()) && (pPV == NULL)) {
    if (!iVtx->isFake() && 
	(iVtx->ndof() > 4) && 
	(fabs(iVtx->x()) <= 24.0/*cm*/) && 
	(fabs(iVtx->position().Rho()) <= 2.0/*cm*/)) pPV = const_cast<reco::Vertex*>(&*iVtx);
    ++iVtx;
  }
  return pPV;
}

float Common::getMuonCombPFIso(const reco::Muon& muon, const double PUSubtractionCoeff)
{
  const reco::MuonPFIsolation isoBlock = muon.pfIsolationR04();
  return (isoBlock.sumChargedHadronPt + 
	  std::max(0.0, (double)(isoBlock.sumNeutralHadronEt + isoBlock.sumPhotonEt - 
				 PUSubtractionCoeff*isoBlock.sumPUPt)));
}

std::vector<reco::Muon*>
Common::getTightPFIsolatedRecoMuons(edm::Handle<reco::MuonCollection>& pMuons, 
				    const reco::Vertex* pPV, const double PUSubtractionCoeff, 
				    const double PFIsoMax, const double etaMax, const bool passIso)
{
  return getTightIsolatedRecoMuons(pMuons, pPV, true, PUSubtractionCoeff, PFIsoMax, etaMax, 
				   passIso);
}

std::vector<reco::Muon*>
Common::getTightDetectorIsolatedRecoMuons(edm::Handle<reco::MuonCollection>& pMuons, 
					  const reco::Vertex* pPV, const double detectorIsoMax, 
					  const double etaMax, const bool passIso)
{
  return getTightIsolatedRecoMuons(pMuons, pPV, false, 0.0, detectorIsoMax, etaMax, passIso);
}

std::vector<reco::PFTau*>
Common::getRecoTaus(const edm::Handle<reco::PFTauCollection>& pTaus, 
		    const std::vector<edm::Handle<reco::PFTauDiscriminator> >& pTauDiscriminators, 
		    const double etaMax)
{
  std::vector<reco::PFTau*> taus;
  for (reco::PFTauCollection::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
       ++iTau) {
    reco::PFTauRef tauRef(pTaus, iTau - pTaus->begin());
    bool passTauDiscriminators = true;
    std::vector<edm::Handle<reco::PFTauDiscriminator> >::const_iterator iDiscriminator = 
      pTauDiscriminators.begin();
    while ((iDiscriminator != pTauDiscriminators.end()) && passTauDiscriminators) {
      if ((**iDiscriminator)[tauRef] != 1.0) passTauDiscriminators = false;
      ++iDiscriminator;
    }
    if (passTauDiscriminators && ((etaMax == -1.0) || (fabs(iTau->eta()) < etaMax))) {
      taus.push_back(const_cast<reco::PFTau*>(&*iTau));
    }
  }
  return taus;
}

//set canvas drawing options
void Common::setCanvasOptions(TCanvas& canvas, const Int_t grid, const Int_t logY, 
			      const Int_t logZ)
{
  canvas.SetFillStyle(0);
  canvas.SetFillColor(0);
  canvas.SetGrid(grid, grid);
  canvas.SetLogy(logY);
  canvas.SetLogz(logZ);
}

//set canvas margins
void Common::setCanvasMargins(TCanvas& canvas, const float left, const float top, 
			      const float right, const float bottom)
{
  canvas.cd()->SetLeftMargin(left);
  canvas.cd()->SetTopMargin(top);
  canvas.cd()->SetRightMargin(right);
  canvas.cd()->SetBottomMargin(bottom);
}

//set axis drawing options
void Common::setAxisOptions(TAxis* axis, const Float_t labelSize, const Float_t titleSize, 
			    const Float_t titleOffset, const char* title)
{
  axis->SetLabelFont(42);
  axis->SetLabelOffset(0.007);
  axis->SetLabelSize(labelSize);
  axis->SetTitleFont(42);
  axis->SetTitleSize(titleSize);
  axis->SetTitleOffset(titleOffset);
  axis->SetTitle(title);
}

//set axis labels
void Common::setAxisLabels(TAxis* axis, const std::vector<std::string>& binLabels)
{
  for (Int_t iBin = 1; iBin <= axis->GetNbins(); ++iBin) {
    if (iBin <= (int)binLabels.size()) axis->SetBinLabel(iBin, binLabels[iBin - 1].c_str());
  }
}

//set graph drawing options
void Common::setGraphOptions(TGraphAsymmErrors& graph, const Color_t color, const Size_t size, 
			     const Style_t style, const char* xAxisTitle, const char* yAxisTitle)
{
  graph.SetMarkerColor(color);
  graph.SetMarkerSize(size);
  graph.SetMarkerStyle(style);
  graph.SetLineColor(color);
  graph.SetTitle("");
  setAxisOptions(graph.GetXaxis(), 0.05, 0.05, 0.9, xAxisTitle);
  setAxisOptions(graph.GetYaxis(), 0.05, 0.05, 1.05, yAxisTitle);
}

//set histogram drawing options
void Common::setHistogramOptions(TH1F* histogram, const Color_t color, const Size_t size, 
				 const Style_t style, const Double_t scale, 
				 const char* xAxisTitle, const char* yAxisTitle, 
				 const Double_t xAxisLabelSize)
{
  histogram->SetMarkerColor(color);
  histogram->SetMarkerSize(size);
  histogram->SetMarkerStyle(style);
  histogram->SetLineColor(color);
  histogram->SetLineWidth(1);
  histogram->SetFillStyle(0);
  setAxisOptions(histogram->GetXaxis(), xAxisLabelSize, 0.05, 0.9, xAxisTitle);
  setAxisOptions(histogram->GetYaxis(), 0.05, 0.05, 1.05, yAxisTitle);
  histogram->Scale(scale);
}

void Common::setHistogramOptions(TH2F* histogram, const Color_t color, const Size_t size, 
				 const Style_t style, const Float_t yAxisTitleOffset, 
				 const Double_t scale, const char* xAxisTitle, 
				 const char* yAxisTitle)
{
  histogram->SetMarkerColor(color);
  histogram->SetMarkerSize(size);
  histogram->SetMarkerStyle(style);
  histogram->SetLineColor(color);
  histogram->SetLineWidth(1);
  histogram->SetFillStyle(0);
  setAxisOptions(histogram->GetXaxis(), 0.04, 0.04, 1.1, xAxisTitle);
  setAxisOptions(histogram->GetYaxis(), 0.04, 0.04, yAxisTitleOffset, yAxisTitle);
  setAxisOptions(histogram->GetZaxis(), 0.04, 0.04, histogram->GetZaxis()->GetTitleOffset(), "");
  histogram->Scale(scale);
}

void Common::setLegendOptions(TLegend& legend, const char* header)
{
  legend.SetFillColor(0);
  legend.SetTextFont(42);
  legend.SetHeader(header);
}

unsigned int
Common::getStatus3Key(const edm::Handle<reco::GenParticleRefVector>& pSelectedGenParticles, 
		      const edm::Handle<reco::GenParticleCollection>& pGenParticles, 
		      const unsigned int i)
{
  reco::GenParticleRef genParticleRef = pSelectedGenParticles->at(i);
  unsigned int key = genParticleRef.key();
  if (genParticleRef->status() == 1) {
    key = reco::GenParticleRef(pGenParticles, key)->motherRef()->motherRef().key();
  }
  return key;
}

std::vector<reco::Muon*>
Common::getTightIsolatedRecoMuons(edm::Handle<reco::MuonCollection>& pMuons, 
				  const reco::Vertex* pPV, const bool usePFIso, 
				  const double PUSubtractionCoeff, const double isoMax, 
				  const double etaMax, const bool passIso)
{
  std::vector<reco::Muon*> tightMuons;
  for (reco::MuonCollection::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); 
       ++iMuon) {
    if (muon::isTightMuon(*iMuon, *pPV) && 
	iMuon->isPFMuon() && 
	(fabs(iMuon->innerTrack()->dz(pPV->position())) < 0.5) && 
	(iMuon->track()->hitPattern().trackerLayersWithMeasurement() > 5) && 
	((etaMax == -1.0) || (fabs(iMuon->eta()) < etaMax))) {
      float iso = 0.0;
      if (usePFIso) {
	iso = getMuonCombPFIso(*iMuon, PUSubtractionCoeff)/iMuon->pt();
      }
      else {
	const reco::MuonIsolation isoBlock = iMuon->isolationR03();
	iso = (isoBlock.sumPt + isoBlock.emEt + isoBlock.hadEt)/iMuon->pt();
      }
      if ((isoMax == -1.0) || (passIso && (iso < isoMax)) || (!passIso && (iso >= isoMax))) {
	tightMuons.push_back(const_cast<reco::Muon*>(&*iMuon));
      }
    }
  }
  return tightMuons;
}
