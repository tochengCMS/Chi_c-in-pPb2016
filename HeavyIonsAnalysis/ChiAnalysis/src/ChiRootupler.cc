// -*- C++ -*-
//
// Package:    ChiRootupler
// Class:      ChiRootupler
// 
/**
 Description: Saves the muon, dimuon and chi candidate information

 Implementation:  Ota Kukral based on work of Andre Stahl and Stefano Argiro, Alessandro Degano  and the Torino team, Alberto Sanchez
*/

// system include files
#include <memory>

// user include files
#include <HeavyIonsAnalysis/ChiAnalysis/interface/ChiRootupler.h>
//#include <HeavyIonsAnalysis/ChiAnalysis/interface/mydict.cxx>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <DataFormats/PatCandidates/interface/UserData.h> 
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"
#include <vector>
#include <sstream>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include <boost/foreach.hpp>



//constructor
ChiRootupler::ChiRootupler(const edm::ParameterSet & iConfig) :
	//muon_label(consumes<edm::View <pat::Muon> >(iConfig.getParameter < edm::InputTag >("muon_cand"))),
	muon_label(consumes< pat::MuonCollection >(iConfig.getParameter < edm::InputTag >("muon_cand"))),
	dimuon_label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag >("dimuon_cand"))),
	photon_label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag >("photon_cand"))),
	conversion_label(consumes<reco::ConversionCollection>(iConfig.getParameter < edm::InputTag >("conversions_ch"))),
	//chi_label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag >("chi_cand"))),
	primaryVertices_label(consumes<reco::VertexCollection>(iConfig.getParameter < edm::InputTag >("primaryVertices"))),
	triggerResults_label(consumes<edm::TriggerResults>(iConfig.getParameter < edm::InputTag >("TriggerResults"))),
	genParticles_label(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag >("genParticlesTag"))),
	flag_doMC(iConfig.getParameter < bool >("isMC"))
{

	// initialize arrays
	muon_p4 = new TClonesArray("TLorentzVector", 1000);
	muonGen_p4 = new TClonesArray("TLorentzVector", 1000);

	dimuon_p4 = new TClonesArray("TLorentzVector", 100);
	dimuon_vtx = new TClonesArray("TVector3", 100);
	
	conv_vtx = new TClonesArray("TVector3", 1000);
	conv_p4 = new TClonesArray("TLorentzVector", 1000);
	convGen_p4 = new TClonesArray("TLorentzVector", 1000);

	gen_chic_p4 = new TClonesArray("TLorentzVector", 100);
	gen_Jpsi_p4 = new TClonesArray("TLorentzVector", 100);
	gen_muon_p4 = new TClonesArray("TLorentzVector", 100);
	gen_phot_p4 = new TClonesArray("TLorentzVector", 100);

	chi_p4 = new TClonesArray("TLorentzVector", 100);


	edm::Service < TFileService > fs;
	event_tree = fs->make < TTree >("event_tree", "Tree with event");

	//general information
	event_tree->Branch("run", &runNumber, "run/L");
	event_tree->Branch("event", &eventNumber, "event/L");
	event_tree->Branch("nPrimVertices", &nPrimVertices, "nPrimVertices/L");
	event_tree->Branch("muonPerEvent_noCuts", &muonPerEvent, "muonPerEvent/I");
	event_tree->Branch("convPerEvent_noCuts", &convPerTriggeredEvent, "convPerEvent/I");
	event_tree->Branch("dimuonPerEvent", &dimuonPerEvent, "convPerEvent/I");

	//	PV
	event_tree->Branch("pvtx_z", &pvtx_z);
	event_tree->Branch("pvtx_x", &pvtx_x);
	event_tree->Branch("pvtx_y", &pvtx_y);
	event_tree->Branch("pvtx_nTracks", &pvtx_nTracks);
	event_tree->Branch("pvtx_isFake", &pvtx_isFake);



	// muon
	event_tree->Branch("muonIsGlobal", &muonIsGlobal);
	event_tree->Branch("muonIsTracker", &muonIsTracker);
	event_tree->Branch("muonIsPF", &muonIsPF);
	event_tree->Branch("muonIsSoft", &muonIsSoft);
	event_tree->Branch("muonIsTight", &muonIsTight);
	event_tree->Branch("muonIsNotGlobalNorTracker", &muonIsNotGlobalNorTracker);
	event_tree->Branch("muonIDHas_TMOneStationTight", &muonIDHas_TMOneStationTight);
	event_tree->Branch("muonInnerTrack_dxy", &muonInnerTrack_dxy);
	event_tree->Branch("muonInnerTrack_dz", &muonInnerTrack_dz);
	event_tree->Branch("muonTrackerLayersWithMeasurement", &muonTrackerLayersWithMeasurement);
	event_tree->Branch("muonPixelLayersWithMeasurement", &muonPixelLayersWithMeasurement);
	event_tree->Branch("muonQuality_isHighPurity", &muonQuality_isHighPurity);
	event_tree->Branch("muon_charge", &muon_charge);
	event_tree->Branch("muon_eta", &muon_eta);
	event_tree->Branch("muon_pt", &muon_pt);
	event_tree->Branch("muon_p4", "TClonesArray", &muon_p4, 32000, 0);
	if (flag_doMC) {
		event_tree->Branch("muon_isMatchedMC", &muon_isMatchedMC);
		event_tree->Branch("muonGen_eta", &muonGen_eta);
		event_tree->Branch("muonGen_pt", &muonGen_pt);
		event_tree->Branch("muonGen_p4", &muonGen_p4, 32000, 0);
		event_tree->Branch("muonGen_rDelta", &muonGen_rDelta);
		event_tree->Branch("muonGen_ptDelta", &muonGen_ptDelta);
		event_tree->Branch("muonGen_ptDeltaRel", &muonGen_ptDeltaRel);
	}
	if (flag_saveExtraThings) {
		event_tree->Branch("patMuon", "std::vector <pat::Muon>", &patMuonStored);
	}

	// dimuon - TBU

	event_tree->Branch("dimuon_p4", "TClonesArray", &dimuon_p4, 32000, 0);
	event_tree->Branch("dimuon_eta", &dimuon_eta);
	event_tree->Branch("dimuon_pt", &dimuon_pt);
	event_tree->Branch("dimuon_charge", &dimuon_charge);
	event_tree->Branch("dimuon_vtx", "TClonesArray", &dimuon_vtx, 32000, 0);
	event_tree->Branch("dimuon_dz_dimuonvtx_pvtx", &dimuon_dz_dimuonvtx_pvtx);

	if (flag_saveExtraThings) {
		event_tree->Branch("dimuonStored", "std::vector <pat::CompositeCandidate>", &dimuonStored);
	}
	event_tree->Branch("dimuon_muon1_position", &dimuon_muon1_position);
	event_tree->Branch("dimuon_muon2_position", &dimuon_muon2_position);
	event_tree->Branch("dimuon_ctpv", &dimuon_ctpv);
	event_tree->Branch("dimuon_ctpvError", &dimuon_ctpvError);





	// conversions
	event_tree->Branch("convQuality_isHighPurity", &convQuality_isHighPurity);
	event_tree->Branch("convQuality_isGeneralTracksOnly", &convQuality_isGeneralTracksOnly);
	event_tree->Branch("conv_vtx", "TClonesArray", &conv_vtx, 32000, 0);
	event_tree->Branch("conv_vertexPositionRho", &conv_vertexPositionRho);
	event_tree->Branch("conv_sigmaTkVtx1", &conv_sigmaTkVtx1);
	event_tree->Branch("conv_sigmaTkVtx2", &conv_sigmaTkVtx2);
	event_tree->Branch("conv_tkVtxCompatibilityOK", &conv_tkVtxCompatibilityOK);
	event_tree->Branch("conv_tkVtxCompatible_bestVertex", &conv_tkVtxCompatible_bestVertex);
	event_tree->Branch("conv_tkVtxCompatible_secondBestVertexA", &conv_tkVtxCompatible_secondBestVertexA);
	event_tree->Branch("conv_tkVtxCompatible_secondBestVertexB", &conv_tkVtxCompatible_secondBestVertexB);
	event_tree->Branch("conv_compatibleInnerHitsOK", &conv_compatibleInnerHitsOK);
	if (flag_saveExtraThings) {
		event_tree->Branch("conv_hitPat1", "reco::HitPattern", &conv_hitPat1);
		event_tree->Branch("conv_hitPat2", "reco::HitPattern", &conv_hitPat2);
		event_tree->Branch("conv_tkVtxCompatibilityOK_test", &conv_tkVtxCompatibilityOK_test);
		event_tree->Branch("conv_tkVtxCompatible_bestVertex_test", &conv_tkVtxCompatible_bestVertex_test);
		event_tree->Branch("conv_tkVtxCompatible_secondBestVertexA_test", &conv_tkVtxCompatible_secondBestVertexA_test);
		event_tree->Branch("conv_tkVtxCompatible_secondBestVertexB_test", &conv_tkVtxCompatible_secondBestVertexB_test);
	}
	event_tree->Branch("conv_vertexChi2Prob", &conv_vertexChi2Prob);
	event_tree->Branch("conv_zOfPriVtx", &conv_zOfPriVtx);
	event_tree->Branch("conv_zOfPriVtxFromTracks", &conv_zOfPriVtxFromTracks);
	event_tree->Branch("conv_dzToClosestPriVtx", &conv_dzToClosestPriVtx);
	event_tree->Branch("conv_dxyPriVtx_Tr1", &conv_dxyPriVtx_Tr1);
	event_tree->Branch("conv_dxyPriVtx_Tr2", &conv_dxyPriVtx_Tr2);
	event_tree->Branch("conv_dxyPriVtxTimesCharge_Tr1", &conv_dxyPriVtxTimesCharge_Tr1);
	event_tree->Branch("conv_dxyPriVtxTimesCharge_Tr2", &conv_dxyPriVtxTimesCharge_Tr2);
	event_tree->Branch("conv_dxyError_Tr1", &conv_dxyError_Tr1);
	event_tree->Branch("conv_dxyError_Tr2", &conv_dxyError_Tr2);

	event_tree->Branch("conv_tk1NumOfDOF", &conv_tk1NumOfDOF);
	event_tree->Branch("conv_tk2NumOfDOF", &conv_tk2NumOfDOF);
	event_tree->Branch("conv_track1Chi2", &conv_track1Chi2);
	event_tree->Branch("conv_track2Chi2", &conv_track2Chi2);
	event_tree->Branch("conv_minDistanceOfApproach", &conv_minDistanceOfApproach);

	event_tree->Branch("conv_p4", "TClonesArray", &conv_p4, 32000, 0);
	event_tree->Branch("conv_eta", &conv_eta);
	event_tree->Branch("conv_pt", &conv_pt);



	if (flag_doMC) {
		event_tree->Branch("conv_isMatchedMC", &conv_isMatchedMC);
		event_tree->Branch("convGen_eta", &convGen_eta);
		event_tree->Branch("convGen_pt", &convGen_pt);
		event_tree->Branch("convGen_p4", "TClonesArray", &convGen_p4, 32000, 0);
		event_tree->Branch("convGen_rDelta", &convGen_rDelta);
		event_tree->Branch("convGen_ptDelta", &convGen_ptDelta);
		event_tree->Branch("convGen_ptDeltaRel", &convGen_ptDeltaRel);
		event_tree->Branch("convGen_motherCode", &convGen_motherCode);
	}

	
	// MC general
	if (flag_doMC) {
		event_tree->Branch("gen_pdgId", &gen_pdgId);

		event_tree->Branch("gen_chic_pt", &gen_chic_pt);
		event_tree->Branch("gen_chic_eta", &gen_chic_eta);
		event_tree->Branch("gen_chic_p4", "TClonesArray", &gen_chic_p4, 32000, 0);

		event_tree->Branch("gen_Jpsi_pt", &gen_Jpsi_pt);
		event_tree->Branch("gen_Jpsi_eta", &gen_Jpsi_eta);
		event_tree->Branch("gen_Jpsi_p4", "TClonesArray", &gen_Jpsi_p4, 32000, 0);
		event_tree->Branch("gen_Jpsi_matchPosition", &gen_Jpsi_matchPosition);
		event_tree->Branch("gen_Jpsi_nMatches", &gen_Jpsi_nMatches);
		event_tree->Branch("gen_Jpsi_rDelta", &gen_Jpsi_rDelta);
		event_tree->Branch("gen_Jpsi_ptDeltaRel", &gen_Jpsi_ptDeltaRel);

		event_tree->Branch("gen_muon_charge", &gen_muon_charge);
		event_tree->Branch("gen_muon_pt", &gen_muon_pt);
		event_tree->Branch("gen_muon_eta", &gen_muon_eta);
		event_tree->Branch("gen_muon_p4", "TClonesArray", &gen_muon_p4, 32000, 0);
		event_tree->Branch("gen_muon_matchPosition", &gen_muon_matchPosition);
		event_tree->Branch("gen_muon_nMatches", &gen_muon_nMatches);
		event_tree->Branch("gen_muon_rDelta", &gen_muon_rDelta);
		event_tree->Branch("gen_muon_ptDeltaRel", &gen_muon_ptDeltaRel);


		event_tree->Branch("gen_phot_pt", &gen_phot_pt);
		event_tree->Branch("gen_phot_eta", &gen_phot_eta);
		event_tree->Branch("gen_phot_p4", "TClonesArray", &gen_phot_p4, 32000, 0);
		event_tree->Branch("gen_conv_matchPosition", &gen_conv_matchPosition);
		event_tree->Branch("gen_conv_nMatches", &gen_conv_nMatches);
		event_tree->Branch("gen_conv_rDelta", &gen_conv_rDelta);
		event_tree->Branch("gen_conv_ptDeltaRel", &gen_conv_ptDeltaRel);
	}

	//chi

	event_tree->Branch("chi_p4", "TClonesArray", &chi_p4, 32000, 0);
	event_tree->Branch("chi_eta", &chi_eta);
	event_tree->Branch("chi_pt", &chi_pt);
	event_tree->Branch("chi_daughterJpsi_position", &chi_daughterJpsi_position);
	event_tree->Branch("chi_daughterConv_position", &chi_daughterConv_position);
	event_tree->Branch("chi_dzPhotToDimuonVtx", &chi_dzPhotToDimuonVtx);
	event_tree->Branch("chi_dxyPhotToDimuonVtx", &chi_dxyPhotToDimuonVtx);
	event_tree->Branch("chiStored", "std::vector <pat::CompositeCandidate>", &chiStored);
	//event_tree->Branch("chi_pt", &chi_pt);


	// chi - TBD
	chi_tree = fs->make < TTree >("chiTree", "Tree of chi");

	chi_tree->Branch("run", &runNumber, "run/L");
	chi_tree->Branch("event", &eventNumber, "event/L");
	chi_tree->Branch("chi_cand", "pat::CompositeCandidate", &chi_cand);
	chi_tree->Branch("chi_dzPhotToDimuonVtx", &chi_dzPhotToDimuonVtx, "chi_dzPhotToDimuonVtx/D");

	chi_tree->Branch("primary_v", "TVector3", &primary_v);
	chi_tree->Branch("secondary_v", "TVector3", &secondary_v);
	chi_tree->Branch("dimuon_v", "TVector3", &dimuon_v);

	//chi_tree->Branch("chi_p4", "TLorentzVector", &chi_p4);
	chi_tree->Branch("chi_dimuon_p4", "TLorentzVector", &chi_dimuon_p4);
	chi_tree->Branch("muonP_p4", "TLorentzVector", &muonP_p4);
	chi_tree->Branch("muonN_p4", "TLorentzVector", &muonN_p4);
	chi_tree->Branch("photon_p4", "TLorentzVector", &photon_p4);

	chi_tree->Branch("ele_lowerPt_pt", &ele_lowerPt_pt, "ele_lowerPt_pt/D");
	chi_tree->Branch("ele_higherPt_pt", &ele_higherPt_pt, "ele_higherPt_pt/D");

	chi_tree->Branch("trigger", &trigger, "trigger/I");

	/*if (flag_doMC) {
		chi_tree->Branch("gen_chic_p4", "TLorentzVector", &gen_chic_p4);
		chi_tree->Branch("chic_pdgId", &chic_pdgId, "chic_pdgId/I");
		chi_tree->Branch("gen_jpsi_p4", "TLorentzVector", &gen_jpsi_p4);
		chi_tree->Branch("gen_photon_p4", "TLorentzVector", &gen_photon_p4);
		chi_tree->Branch("gen_muonP_p4", "TLorentzVector", &gen_muonP_p4);
		chi_tree->Branch("gen_muonM_p4", "TLorentzVector", &gen_muonM_p4);
	}
	genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
	packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");*/

}

ChiRootupler::~ChiRootupler() {}

//
// member functions
//

////Check recursively if any ancestor of particle is the given one
//bool ChiRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
//   if (ancestor == particle ) return true;
//   if (particle->numberOfMothers() && isAncestor(ancestor,particle->mother(0))) return true;
//   return false;
//}

// ------------ method called for each event  ------------
void ChiRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {
	using namespace edm;
	using namespace std;

	//edm::Handle  < edm::View <pat::Muon> > muon_handle;
	edm::Handle  < std::vector <pat::Muon> > muon_handle;
	iEvent.getByToken(muon_label, muon_handle);

	edm::Handle < std::vector < pat::CompositeCandidate > >dimuon_handle;
	iEvent.getByToken(dimuon_label, dimuon_handle);

	edm::Handle < std::vector < pat::CompositeCandidate > >photon_handle;
	iEvent.getByToken(photon_label, photon_handle);

	edm::Handle< std::vector <reco::Conversion>> conversion_handle;
	iEvent.getByToken(conversion_label, conversion_handle);

	//edm::Handle < std::vector < pat::CompositeCandidate >>chi_handle;
	//iEvent.getByToken(chi_label, chi_handle);

	edm::Handle  < reco::VertexCollection> primaryVertices_handle;
	iEvent.getByToken(primaryVertices_label, primaryVertices_handle);

	edm::Handle < edm::TriggerResults > triggerResults_handle;
	iEvent.getByToken(triggerResults_label, triggerResults_handle);

	edm::Handle <reco::GenParticleCollection> genParticles_handle;
	iEvent.getByToken(genParticles_label, genParticles_handle);

	

	//edm::Handle<reco::VertexCollection> priVtxs;
	//event.getByToken(thePVsToken_, priVtxs);



	//edm::Handle<reco::PFCandidateCollection> pfcandidates;
	//event.getByToken(pfCandidateCollectionToken_, pfcandidates);

	//const reco::PFCandidateCollection pfphotons = selectPFPhotons(*pfcandidates);






	/////////////////////
	//// S T A R T //////
	////////////////////

	//general info
	runNumber = iEvent.id().run();
	eventNumber = iEvent.id().event();
	nPrimVertices = 0;
	if (primaryVertices_handle.isValid()) {
		nPrimVertices = primaryVertices_handle->size();
	}
	muonPerEvent = 0;
	if (muon_handle.isValid()) {
		muonPerEvent = muon_handle->size(); //all muons without any cuts
	}
	convPerTriggeredEvent = 0;
	if (conversion_handle.isValid())
	{
		convPerTriggeredEvent = conversion_handle->size(); //all conversions without any cuts
	}
	dimuonPerEvent = 0;
	if (dimuon_handle.isValid())
	{
		dimuonPerEvent = dimuon_handle->size(); //all conversions without any cuts
	}

	//PV
	if (primaryVertices_handle.isValid()) {
		for (uint i = 0; i < primaryVertices_handle->size(); i++) {
			const reco::Vertex& pvtx = primaryVertices_handle->at(i);
			pvtx_z.push_back(pvtx.z());
			pvtx_x.push_back(pvtx.x());
			pvtx_y.push_back(pvtx.y());
			pvtx_nTracks.push_back(pvtx.nTracks());
			pvtx_isFake.push_back(pvtx.isFake());
		}

	//BOOST_FOREACH(const reco::Vertex& vtx, *primaryVertices_handle.product()) {
	//	(*primaryVertices_handle.product())[closest_pv_index].z()
	//	i++;
	}
	else cout << "Problem with PV handle" << endl;

	int pvtx_index = 0; //top primary vertex used for now

	//muons
	if (muon_handle.isValid()) {
		for (uint i = 0; i < muon_handle->size(); i++) {
			const pat::Muon& patMuon = muon_handle->at(i);
			const reco::Vertex& pvtx = primaryVertices_handle->at(pvtx_index); //use our selected PV
			muonIsGlobal.push_back(patMuon.isGlobalMuon());
			muonIsTracker.push_back(patMuon.isTrackerMuon());
			muonIsPF.push_back(patMuon.isPFMuon());
			muonIsSoft.push_back(patMuon.isSoftMuon(pvtx));
			muonIsTight.push_back(patMuon.isTightMuon(pvtx));
			if (!patMuon.isGlobalMuon() && !patMuon.isTrackerMuon()) { muonIsNotGlobalNorTracker.push_back(true); }
			else muonIsNotGlobalNorTracker.push_back(false); // just for convenience

			muonIDHas_TMOneStationTight.push_back(patMuon.muonID("TMOneStationTight"));
			if (patMuon.isGlobalMuon() || patMuon.isTrackerMuon()) {
				muonInnerTrack_dxy.push_back(patMuon.innerTrack()->dxy());
				muonInnerTrack_dz.push_back(patMuon.innerTrack()->dz());
				muonTrackerLayersWithMeasurement.push_back(patMuon.innerTrack()->hitPattern().trackerLayersWithMeasurement());
				muonPixelLayersWithMeasurement.push_back(patMuon.innerTrack()->hitPattern().pixelLayersWithMeasurement());
				reco::TrackBase::TrackQuality tq = reco::TrackBase::qualityByName("highPurity");//high purity=2 //see DataFormats/TrackReco/interface/TrackBase.h
				muonQuality_isHighPurity.push_back(patMuon.innerTrack()->quality(tq));
			}
			else
			{
				muonInnerTrack_dxy.push_back(-100);
				muonInnerTrack_dz.push_back(-100);
				muonTrackerLayersWithMeasurement.push_back(-1);
				muonPixelLayersWithMeasurement.push_back(-1);
				muonQuality_isHighPurity.push_back(0);
			}
			muon_charge.push_back(patMuon.charge());
			muon_eta.push_back(patMuon.eta());
			muon_pt.push_back(patMuon.pt());
			TLorentzVector muon_p4_aux;
			muon_p4_aux.SetPtEtaPhiM(patMuon.pt(), patMuon.eta(), patMuon.phi(), patMuon.mass());
			new ((*muon_p4)[i]) TLorentzVector(muon_p4_aux);
			// MC generated information
			if (flag_doMC) {
				bool muon_isMatchedMC_flag = patMuon.genLepton(); //false for null pointer, true if match exists, from pat::Lepton.h
				muon_isMatchedMC.push_back(muon_isMatchedMC_flag);
				if (muon_isMatchedMC_flag) {
					const reco::GenParticle genMuon = *patMuon.genLepton();
					muonGen_eta.push_back(genMuon.eta());
					muonGen_pt.push_back(genMuon.pt());
					TLorentzVector muonGen_p4_aux;
					muonGen_p4_aux.SetPtEtaPhiM(genMuon.pt(), genMuon.eta(), genMuon.phi(), genMuon.mass());
					new ((*muonGen_p4)[i]) TLorentzVector(muonGen_p4_aux);
					muonGen_rDelta.push_back(reco::deltaR(patMuon, genMuon)); //sqrt(phi^2+eta^2)
					muonGen_ptDelta.push_back(patMuon.pt() - genMuon.pt());
					muonGen_ptDeltaRel.push_back((patMuon.pt() - genMuon.pt()) / genMuon.pt()); //defined in the matcher to be divided by genPt //MCTruthMatchers.cc
				}
				else { //default values to store if no match - in principle can be ommitted, just for direct looking at branches
					muonGen_eta.push_back(-5);
					muonGen_pt.push_back(0);
					new ((*muonGen_p4)[i]) TLorentzVector(TLorentzVector());
					muonGen_rDelta.push_back(-5);
					muonGen_ptDelta.push_back(-5);
					muonGen_ptDeltaRel.push_back(-5);
				}
			}


			if (flag_saveExtraThings)
			{
				patMuonStored.push_back(patMuon);
			}

		}
	}
	else cout << "Problem with muon handle" << endl;
	
	//
	//if (genParticles_handle.isValid()) {
	//	//cout << "Gen: "<<genParticles_handle->size() << endl;
	//	/*for (uint i = 0; i < genParticles_handle->size(); i++) {
	//		const reco::GenParticle& genParticle = genParticles_handle->at(i);// genParticle(size_t idx = 0) const;
	//		if (i == 100) { cout << "GEN pT " << genParticle.pt() << endl; }
	//	}*/
	//}
	//else cout << "Problem with gen handle" << endl;
	

	if (dimuon_handle.isValid()) 
	{
		for (uint i = 0; i < dimuon_handle->size(); i++) {
			const pat::CompositeCandidate& dimuon = dimuon_handle->at(i);
			TLorentzVector dimuon_p4_aux;
			dimuon_p4_aux.SetPtEtaPhiM(dimuon.pt(), dimuon.eta(), dimuon.phi(), dimuon.mass());
			new ((*dimuon_p4)[i]) TLorentzVector(dimuon_p4_aux);

			dimuon_charge.push_back(dimuon.charge());
			dimuon_eta.push_back(dimuon.eta());
			dimuon_pt.push_back(dimuon.pt());

			if (flag_saveExtraThings)
			{
				dimuonStored.push_back(dimuon);
			}

			const reco::Vertex* dimuon_recovtx = dimuon.userData<reco::Vertex>("commonVertex");
			TVector3 dimuon_vtx_aux;
			dimuon_vtx_aux.SetXYZ(dimuon_recovtx->x(), dimuon_recovtx->y(), dimuon_recovtx->z());
			new ((*dimuon_vtx)[i]) TVector3(dimuon_vtx_aux);
			dimuon_dz_dimuonvtx_pvtx.push_back(dimuon_recovtx->z() - (*primaryVertices_handle.product())[pvtx_index].z());
			int muonPos1 = dimuon.userInt("muonPosition1");
			int muonPos2 = dimuon.userInt("muonPosition2");
			if (fabs(muon_handle->at(muonPos2).pt() - dimuon.daughter("muon2")->pt()) > 0.01) { cout << "SOMETHING WRONG WITH THE MATCHING FROM DIMUON TO MUON - possibly different muon collections used"<<endl; }//muon matching assumes that muon collection going to dimuon producer and here is the same
			dimuon_muon1_position.push_back(muonPos1);
			dimuon_muon2_position.push_back(muonPos2);
			dimuon_ctpv.push_back(dimuon.userFloat("ppdlPV"));
			dimuon_ctpvError.push_back(dimuon.userFloat("ppdlErrPV"));
			//cout << "Muon positions  " << muonPos1 << "   " << muonPos2 << endl;

			//double pokus2 = (dynamic_cast <pat::CompositeCandidate *>(dimuon.daughter("dimuon")))->userFloat("ppdlErrPV");
			//double ptTest = dimuon.daughter("muon1")->pt();
			//double etaTest = dimuon.daughter("muon1")->eta();
			//double ptTestMu = muon_handle->at(muonPos1).pt();
			//double etaTestMu = muon_handle->at(muonPos1).eta();
			//cout << "pt from daughter and muon collection  " << ptTest << "  " << ptTestMu << " and eta " <<etaTest<<"   " << etaTestMu<<endl << endl;
			////dimuon.userData<reco::Vertex>("PVwithmuons");
			//const reco::Vertex* opv = dimuon.userData<reco::Vertex>("PVwithmuons");
			////const reco::Vertex* opv = dimuon.userData<reco::Vertex>("commonVertex");
			////const reco::Vertex* opv = (dynamic_cast < pat::CompositeCandidate*> (&dimuon))->userData<reco::Vertex>("PVwithmuons");
			////cout << opv << endl;
			//primary_v.SetXYZ(opv->x(), opv->y(), opv->z());
			////cout << "Ntracks in vertex  " << opv->tracksSize() << endl;
			////cout << eventNumber << "   " << i << " " << dimuon.userFloat("DCA")<< endl;
			////cout << opv->x() << endl;
			//dimuonStored = dimuon;
			//dimuon_tree->Fill();
		}
	}
	else cout << "Problem with dimuon handle" << endl;



	// conversions

	if (conversion_handle.isValid())
	{
		for (uint i = 0; i < conversion_handle->size(); i++) {
			const reco::Conversion& candPhoton = conversion_handle->at(i);
			TLorentzVector conv_p4_aux;
			conv_p4_aux.SetXYZT(candPhoton.refittedPair4Momentum().x(), candPhoton.refittedPair4Momentum().y(), candPhoton.refittedPair4Momentum().z(), candPhoton.refittedPair4Momentum().t());
			new ((*conv_p4)[i]) TLorentzVector(conv_p4_aux);
			conv_eta.push_back(conv_p4_aux.Eta());
			conv_pt.push_back(conv_p4_aux.Pt());

			convQuality_isHighPurity.push_back(candPhoton.quality((reco::Conversion::ConversionQuality)(8))); //8 is high purity, see reco::Conversion Class Reference
			convQuality_isGeneralTracksOnly.push_back(candPhoton.quality((reco::Conversion::ConversionQuality)(0))); //0 is general tracks only, see reco::Conversion Class Reference
			const reco::Vertex conv_recovtx = candPhoton.conversionVertex();
			TVector3 conv_vtx_aux;
			conv_vtx_aux.SetXYZ(conv_recovtx.x(), conv_recovtx.y(), conv_recovtx.z());
			new ((*conv_vtx)[i]) TVector3(conv_vtx_aux);
			conv_vertexPositionRho.push_back(candPhoton.conversionVertex().position().rho());
			bool conv_tkVtxCompatible_bestVertex_aux, conv_tkVtxCompatible_secondBestVertexA_aux, conv_tkVtxCompatible_secondBestVertexB_aux;
			double conv_sigmaTkVtx1_aux, conv_sigmaTkVtx2_aux;
			if (flag_saveExtraThings)
			{
				conv_tkVtxCompatibilityOK_test.push_back(Conv_checkTkVtxCompatibility(candPhoton, *primaryVertices_handle.product(), 20, conv_tkVtxCompatible_bestVertex_aux, conv_tkVtxCompatible_secondBestVertexA_aux, conv_tkVtxCompatible_secondBestVertexB_aux, conv_sigmaTkVtx1_aux, conv_sigmaTkVtx2_aux));
				conv_tkVtxCompatible_bestVertex_test.push_back(conv_tkVtxCompatible_bestVertex_aux);
				conv_tkVtxCompatible_secondBestVertexA_test.push_back(conv_tkVtxCompatible_secondBestVertexA_aux);
				conv_tkVtxCompatible_secondBestVertexB_test.push_back(conv_tkVtxCompatible_secondBestVertexB_aux);
			}
			conv_tkVtxCompatibilityOK.push_back(Conv_checkTkVtxCompatibility(candPhoton, *primaryVertices_handle.product(), conv_TkVtxCompSigmaCut, conv_tkVtxCompatible_bestVertex_aux, conv_tkVtxCompatible_secondBestVertexA_aux, conv_tkVtxCompatible_secondBestVertexB_aux, conv_sigmaTkVtx1_aux, conv_sigmaTkVtx2_aux));
			conv_tkVtxCompatible_bestVertex.push_back(conv_tkVtxCompatible_bestVertex_aux);
			conv_tkVtxCompatible_secondBestVertexA.push_back(conv_tkVtxCompatible_secondBestVertexA_aux);
			conv_tkVtxCompatible_secondBestVertexB.push_back(conv_tkVtxCompatible_secondBestVertexB_aux);
			conv_sigmaTkVtx1.push_back(conv_sigmaTkVtx1_aux);
			conv_sigmaTkVtx2.push_back(conv_sigmaTkVtx2_aux);
			
			if (candPhoton.tracks().size() == 2) {
				const edm::RefToBase<reco::Track> conv_tk1 = candPhoton.tracks().at(0);
				const edm::RefToBase<reco::Track> conv_tk2 = candPhoton.tracks().at(1);

				reco::HitPattern hitPatA = conv_tk1->hitPattern();
				reco::HitPattern hitPatB = conv_tk2->hitPattern();
				conv_hitPat1.push_back(hitPatA);
				conv_hitPat2.push_back(hitPatB);
				conv_compatibleInnerHitsOK.push_back((Conv_foundCompatibleInnerHits(hitPatA, hitPatB) && Conv_foundCompatibleInnerHits(hitPatB, hitPatA)));

				//find vertex that points closest - // right now turned off, the 0 vertex used
				int closest_pv_index = 0;
				//int i = 0;
				//BOOST_FOREACH(const reco::Vertex& vtx, *primaryVertices_handle.product()) {
					//if (fabs(candPhoton.zOfPrimaryVertexFromTracks(vtx.position()) - vtx.z()) < fabs(candPhoton.zOfPrimaryVertexFromTracks((*primaryVertices_handle.product())[closest_pv_index].position()) - (*primaryVertices_handle.product())[closest_pv_index].z())) { closest_pv_index = i; }
					//i++;
				//}
				conv_zOfPriVtx.push_back((*primaryVertices_handle.product())[closest_pv_index].z());
				conv_zOfPriVtxFromTracks.push_back(candPhoton.zOfPrimaryVertexFromTracks((*primaryVertices_handle.product())[closest_pv_index].position()));
				conv_dzToClosestPriVtx.push_back(candPhoton.zOfPrimaryVertexFromTracks((*primaryVertices_handle.product())[closest_pv_index].position()) - (*primaryVertices_handle.product())[closest_pv_index].z());
				// Now check impact parameter wtr with the just found closest primary vertex
				conv_dxyPriVtx_Tr1.push_back(conv_tk1->dxy((*primaryVertices_handle.product())[closest_pv_index].position()));
				conv_dxyPriVtx_Tr2.push_back(conv_tk2->dxy((*primaryVertices_handle.product())[closest_pv_index].position()));
				conv_dxyPriVtxTimesCharge_Tr1.push_back(conv_tk1->dxy((*primaryVertices_handle.product())[closest_pv_index].position())*conv_tk1->charge());
				conv_dxyPriVtxTimesCharge_Tr2.push_back(conv_tk2->dxy((*primaryVertices_handle.product())[closest_pv_index].position())*conv_tk2->charge());
				conv_dxyError_Tr1.push_back(conv_tk1->dxyError());
				conv_dxyError_Tr2.push_back(conv_tk2->dxyError());

				conv_tk1NumOfDOF.push_back(conv_tk1->ndof());
				conv_tk2NumOfDOF.push_back(conv_tk2->ndof());
				conv_track1Chi2.push_back(candPhoton.tracks().at(0)->normalizedChi2());
				conv_track2Chi2.push_back(candPhoton.tracks().at(1)->normalizedChi2());

			}
			else conv_compatibleInnerHitsOK.push_back(-1);

			conv_vertexChi2Prob.push_back(ChiSquaredProbability(candPhoton.conversionVertex().chi2(), candPhoton.conversionVertex().ndof()));
			conv_minDistanceOfApproach.push_back(candPhoton.distOfMinimumApproach());
			//if (candPhoton.distOfMinimumApproach() > -10 && candPhoton.distOfMinimumApproach() < 10) { conv_minDistanceOfApproach = candPhoton.distOfMinimumApproach(); }
			//else conv_minDistanceOfApproach = 0;

			//MC for conversions
			if (flag_doMC)
			{
				if (genParticles_handle.isValid()) {
					reco::GenParticle genConv_best = reco::GenParticle();
					bool conv_isMatchedMC_aux = false;
					for (uint i = 0; i < genParticles_handle->size(); i++) {
						const reco::GenParticle& genParticle = genParticles_handle->at(i);

						int pdgId = genParticle.pdgId();
						if (pdgId != 22) { continue; } //if not photon, don't bother
						if (genParticle.status() != 1) {
							//cout << "notStable" << endl;
							//cout << genParticle.daughter(0)->pdgId() << endl;
							//cout << "ptdif " << genParticle.daughter(0)->pt() - genParticle.pt() << endl;
							continue;
						} //if not stable, don't bother

						bool genParticleMatched = false; //this particular - is it matched?
						genParticleMatched = Conv_isMatched(candPhoton.refittedPair4Momentum(), genParticle, conv_maxDeltaR, conv_maxDPtRel);
						if (genParticleMatched == true) {
							if (conv_isMatchedMC_aux == false) {//first one found
								conv_isMatchedMC_aux = true;
								genConv_best = genParticle;
							}
							else { //check whether the second match is better by deltaR than the first, save the better one
								if (reco::deltaR(candPhoton.refittedPair4Momentum(), genConv_best) > reco::deltaR(candPhoton.refittedPair4Momentum(), genParticle)) { genConv_best = genParticle; }
							}
						}
					}
					if (conv_isMatchedMC_aux) {
						conv_isMatchedMC.push_back(true);
						convGen_eta.push_back(genConv_best.eta());
						convGen_pt.push_back(genConv_best.pt());
						TLorentzVector convGen_p4_aux;
						convGen_p4_aux.SetPtEtaPhiM(genConv_best.pt(), genConv_best.eta(), genConv_best.phi(), genConv_best.mass());
						new ((*convGen_p4)[i]) TLorentzVector(convGen_p4_aux);

						convGen_rDelta.push_back(reco::deltaR(candPhoton.refittedPair4Momentum(), genConv_best)); //sqrt(phi^2+eta^2)
						convGen_ptDelta.push_back(candPhoton.refittedPair4Momentum().pt() - genConv_best.pt());
						convGen_ptDeltaRel.push_back((candPhoton.refittedPair4Momentum().pt() - genConv_best.pt()) / genConv_best.pt()); //defined in the matcher to be divided by genPt //MCTruthMatchers.cc
						convGen_motherCode.push_back(genConv_best.mother()->pdgId());
					}
					else { //default values to store if no match - in principle can be ommitted, just for direct looking at branches
						conv_isMatchedMC.push_back(false);
						convGen_eta.push_back(-5);
						convGen_pt.push_back(0);
						new ((*convGen_p4)[i]) TLorentzVector(TLorentzVector());
						convGen_rDelta.push_back(-5);
						convGen_ptDelta.push_back(-5);
						convGen_ptDeltaRel.push_back(-5);
						convGen_motherCode.push_back(-1);
					}
				}
				else cout << "Problem with gen handle" << endl;
			}
		}
	}
	else { cout << "Conversions handle problem" << endl; }



//////////////////////
//    C H I       ////
/////////////////////


	/// create the collection 

	int n_chic = 0;  // counter for n chic in event
	pat::CompositeCandidateCollection* chiCandColl = new pat::CompositeCandidateCollection;
	for (pat::CompositeCandidateCollection::const_iterator dimuonCand = dimuon_handle->begin(); dimuonCand != dimuon_handle->end(); ++dimuonCand)
	{
		// use only trigger-matched Jpsi or Upsilon if so requested 
		//if (triggerMatch_) {
		//if (!dimuonCand->userInt("isTriggerMatched")) continue;
		//}

		// loop on conversion candidates, make chi cand
		for (pat::CompositeCandidateCollection::const_iterator photCand = photon_handle->begin(); photCand != photon_handle->end(); ++photCand) {

			chi_cand = makeChiCandidate(*dimuonCand, *photCand);

			//if (!cutDeltaMass(chiCand, *dimuonCand)) {
			//	delta_mass_fail++;
			//	continue;
			//}



																				  //if (!cutdz(dz)) {
																				  //	dz_cut_fail++;
																				  //	continue;
																				  //}


			TLorentzVector chi_p4_aux;
			chi_p4_aux.SetPtEtaPhiM(chi_cand.pt(), chi_cand.eta(), chi_cand.phi(), chi_cand.mass());
			new ((*chi_p4)[n_chic]) TLorentzVector(chi_p4_aux);
			++n_chic;
			chi_eta.push_back(chi_p4_aux.Eta());
			chi_pt.push_back(chi_p4_aux.Pt());
			chi_daughterJpsi_position.push_back(std::distance(dimuon_handle->begin(), dimuonCand));
			chi_daughterConv_position.push_back(std::distance(photon_handle->begin(), photCand));
			chi_dzPhotToDimuonVtx.push_back(Getdz(*photCand, dimuonCand->vertex())); 
			chi_dxyPhotToDimuonVtx.push_back(Getdxy(*photCand, dimuonCand->vertex()));

			chiCandColl->push_back(chi_cand);
		}
	}

	// old stuff, to be moved to the loop above and updated

	if (chiCandColl) {
		unsigned int csize = chiCandColl->size();
		for (unsigned int i = 0; i < csize; i++) {
			chi_cand = chiCandColl->at(i);
			chi_dimuon_p4.SetPtEtaPhiM(chi_cand.daughter("dimuon")->pt(), chi_cand.daughter("dimuon")->eta(),
				chi_cand.daughter("dimuon")->phi(), chi_cand.daughter("dimuon")->mass());
			photon_p4.SetPtEtaPhiM(chi_cand.daughter("photon")->pt(), chi_cand.daughter("photon")->eta(),
				chi_cand.daughter("photon")->phi(), chi_cand.daughter("photon")->mass());
			reco::Candidate::LorentzVector vP = chi_cand.daughter("dimuon")->daughter("muon1")->p4();
			reco::Candidate::LorentzVector vM = chi_cand.daughter("dimuon")->daughter("muon2")->p4();

			if (chi_cand.daughter("dimuon")->daughter("muon1")->charge() < 0) {
				vP = chi_cand.daughter("dimuon")->daughter("muon2")->p4();
				vM = chi_cand.daughter("dimuon")->daughter("muon1")->p4();
			}

			muonP_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
			muonN_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());

			Double_t ele1_pt = (dynamic_cast<const pat::CompositeCandidate *>(chi_cand.daughter("photon"))->userData<reco::Track>("track0"))->pt();
			Double_t ele2_pt = (dynamic_cast<const pat::CompositeCandidate *>(chi_cand.daughter("photon"))->userData<reco::Track>("track1"))->pt();

			if (ele1_pt > ele2_pt) {
				ele_higherPt_pt = ele1_pt;
				ele_lowerPt_pt = ele2_pt;
			}
			else {
				ele_higherPt_pt = ele2_pt;
				ele_lowerPt_pt = ele1_pt;
			}

			pi0_abs_mass = 1.;//pi0_abs_values[0];


			reco::Candidate::Point dimuon_vtx = chi_cand.daughter("dimuon")->vertex();
			dimuon_v.SetXYZ(dimuon_vtx.X(), dimuon_vtx.Y(), dimuon_vtx.Z());

			const reco::Vertex *opv = (dynamic_cast <pat::CompositeCandidate *>(chi_cand.daughter("dimuon")))->userData<reco::Vertex>("PVwithmuons");
			//const reco::Vertex *opv = chi_cand.daughter("dimuon")->userData<reco::Vertex>("PVwithmuons");
			primary_v.SetXYZ(opv->x(), opv->y(), opv->z());

			// 2012 parameterization
			//double sigma = Y_sig_par_A + Y_sig_par_B * pow(fabs(dimuon_p4.Rapidity()), 2) +
				//Y_sig_par_C * pow(fabs(dimuon_p4.Rapidity()), 3);

			//psi1S_nsigma = fabs(dimuon_p4.M() - psi1SMass) / sigma;
			//psi2S_nsigma = fabs(dimuon_p4.M() - psi2SMass) / sigma;

			runNumber = iEvent.id().run();
			eventNumber = iEvent.id().event();
			chi_tree->Fill();
		}

	}
	else { std::cout << "no valid chi handle" << std::endl; }


	////////////////
	////// GEN ////
	///////////////

	if (flag_doMC) {
		if (genParticles_handle.isValid()) {
			//bool foundChicInEvent=false;
			for (uint i = 0; i < genParticles_handle->size(); i++) {
				const reco::GenParticle& genParticle = genParticles_handle->at(i);
				int gen_pdgId_aux = genParticle.pdgId();
				if (gen_pdgId_aux != PythCode_chic0 && gen_pdgId_aux != PythCode_chic1 && gen_pdgId_aux != PythCode_chic2) { continue; } //if not chi, skip
				if (genParticle.isLastCopy() == false) { continue; } // if not last copy of a given particle, skip
				//foundChicInEvent
				// it is the chic
				gen_pdgId.push_back(gen_pdgId_aux);
				gen_chic_pt.push_back(genParticle.pt());
				gen_chic_eta.push_back(genParticle.eta());
				TLorentzVector gen_chic_p4_aux;
				gen_chic_p4_aux.SetXYZT(genParticle.p4().x(), genParticle.p4().y(), genParticle.p4().z(), genParticle.p4().t());
				new ((*gen_chic_p4)[i]) TLorentzVector(gen_chic_p4_aux);

				int nDaughters = genParticle.numberOfDaughters();
				if (nDaughters != 2) { cout << "weird gen decay" << endl; } //all of them should be J/psi + gamma
				for (int j = 0; j < nDaughters; j++) {
					const reco::Candidate& gen_chiDaughter = *genParticle.daughter(j);
					//if (gen_chiDaughter->isLastCopy() == false) { cout << "Warning: daughter is not a last copy!" << endl; } //doesn't work
					int dauId = gen_chiDaughter.pdgId();
					if (dauId == 443) {//Jpsi
						gen_Jpsi_pt.push_back(gen_chiDaughter.pt());
						gen_Jpsi_eta.push_back(gen_chiDaughter.eta());
						TLorentzVector gen_Jpsi_p4_aux;
						gen_Jpsi_p4_aux.SetXYZT(gen_chiDaughter.p4().x(), gen_chiDaughter.p4().y(), gen_chiDaughter.p4().z(), gen_chiDaughter.p4().t());
						new ((*gen_Jpsi_p4)[i]) TLorentzVector(gen_Jpsi_p4_aux);
						int nMatches_aux;
						double jpsi_rDelta_aux, jpsi_ptDeltaRel_aux;
						gen_Jpsi_matchPosition.push_back(MatchGen(gen_chiDaughter, dimuon_handle, jpsi_maxDeltaR, jpsi_maxDPtRel, nMatches_aux, jpsi_rDelta_aux, jpsi_ptDeltaRel_aux));
						gen_Jpsi_nMatches.push_back(nMatches_aux);
						gen_Jpsi_rDelta.push_back(jpsi_rDelta_aux);
						gen_Jpsi_ptDeltaRel.push_back(jpsi_ptDeltaRel_aux);

						//check the decay muons
						int nDaughtJpsi = gen_chiDaughter.numberOfDaughters();
						if (nDaughtJpsi != 2) { cout << "weird Jpsi decay" << endl; } //all of them should be mu+mu-
						for (int k = 0; k < nDaughters; k++) { // loop over jpsi daughters
							const reco::Candidate& gen_chiMuon = *gen_chiDaughter.daughter(k);
							//cout << gen_chiMuon.pdgId() << endl;
							if (abs(gen_chiMuon.pdgId()) != 13) { cout << "Warning, jpsi decay daughter is not a muon" << endl; }
							gen_muon_charge.push_back(gen_chiMuon.charge());
							gen_muon_pt.push_back(gen_chiMuon.pt());
							gen_muon_eta.push_back(gen_chiMuon.eta());
							TLorentzVector gen_muon_p4_aux;
							gen_muon_p4_aux.SetXYZT(gen_chiMuon.p4().x(), gen_chiMuon.p4().y(), gen_chiMuon.p4().z(), gen_chiMuon.p4().t());
							new ((*gen_muon_p4)[i]) TLorentzVector(gen_muon_p4_aux);
							int nMatchesMuon_aux;
							double muon_rDelta_aux, muon_ptDeltaRel_aux;
							gen_muon_matchPosition.push_back(MatchGen(gen_chiMuon, muon_handle, muon_maxDeltaR, muon_maxDPtRel, nMatchesMuon_aux, muon_rDelta_aux, muon_ptDeltaRel_aux));
							gen_muon_nMatches.push_back(nMatchesMuon_aux);
							gen_muon_rDelta.push_back(muon_rDelta_aux);
							gen_muon_ptDeltaRel.push_back(muon_ptDeltaRel_aux);
						}
					}
					else if (dauId == 22) {//photon
						gen_phot_pt.push_back(gen_chiDaughter.pt());
						gen_phot_eta.push_back(gen_chiDaughter.eta());
						TLorentzVector gen_phot_p4_aux;
						gen_phot_p4_aux.SetXYZT(gen_chiDaughter.p4().x(), gen_chiDaughter.p4().y(), gen_chiDaughter.p4().z(), gen_chiDaughter.p4().t());
						new ((*gen_phot_p4)[i]) TLorentzVector(gen_phot_p4_aux);
						int nMatches_aux;
						double conv_rDelta_aux, conv_ptDeltaRel_aux;
						gen_conv_matchPosition.push_back(MatchGen(gen_chiDaughter, conversion_handle, conv_maxDeltaR, conv_maxDPtRel, nMatches_aux, conv_rDelta_aux, conv_ptDeltaRel_aux));
						gen_conv_nMatches.push_back(nMatches_aux);
						gen_conv_rDelta.push_back(conv_rDelta_aux);
						gen_conv_ptDeltaRel.push_back(conv_ptDeltaRel_aux);

						//MatchGen(gen_chiDaughter, conv_maxDeltaR, conv_maxDPtRel);

						//int photGen_motherCode = gen_chiDaughter->mother()->pdgId();
						//int photGen_status = gen_chiDaughter->status();
						//cout << "Photon " << photGen_motherCode << " " << photGen_status << endl;
					}
					else {
						cout << "Different decay: " << dauId << endl;
					}
				}
				
				//{ cout << " event and GEN ID, status and n daughters, first id " << eventNumber << "   " << genParticle.pdgId()<< "   " <<st << "    " <<n << "   " << genParticle.daughter(0)->pdgId()<< "  pt:" << genParticle.pt()<< endl; }
				//if (nChicPerEvent > 1) {
					//cout << endl << endl << endl;
				//}

			break; //We found chic, no need to go through the rest of the collection. However, very few events (<1%) exist with multiple chic in them - this takes only the first one into account
			}
		}
		else cout << "Problem with gen handle" << endl;
	}




	//if (flag_doMC) {
	//    // Pruned particles are the one containing "important" stuff
	//    edm::Handle<std::vector<reco::GenParticle> > pruned;
	//    iEvent.getByToken(genCands_, pruned);
	//    // Packed particles are all the status 1, so usable to remake jets
	//    // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
	//    edm::Handle<std::vector<pat::PackedGenParticle> > packed;
	//    iEvent.getByToken(packCands_, packed);

	//    //let's try to find all status1 originating directly from a Chi_c1 meson decay

	//    chic_pdgId = 0;
	//    int foundit = 0;
	//    for (size_t i = 0; i < pruned->size(); i++) {
	//	    int p_id = abs((*pruned)[i].pdgId());
	//	    if (p_id == 20443 || p_id == 445 || p_id == 10443) {
	//	  	  chic_pdgId = p_id;
	//	  	  foundit++;
	//	  	  const reco::Candidate * chic = &(*pruned)[i];
	//	  	  MC_chic_p4.SetPtEtaPhiM(chic->pt(), chic->eta(), chic->phi(), chic->mass());

	//	  	  const reco::Candidate * j1 = nullptr;
	//	  	  const reco::Candidate * j2 = nullptr;

	//	  	  for (size_t j = 0; j < packed->size(); j++) {
	//	  		  //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
	//	  		  const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
	//	  		  if (motherInPrunedCollection != nullptr && isAncestor(chic, motherInPrunedCollection)) {
	//	  			  const reco::Candidate * d = &(*packed)[j];
	//	  			  int dauId = d->pdgId();

	//	  			  if (dauId == 13 && motherInPrunedCollection->pdgId() == 443) {
	//	  				  gen_muonM_p4.SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
	//	  				  foundit++;
	//	  				  j1 = motherInPrunedCollection;
	//	  			  }
	//	  			  if (dauId == -13 && motherInPrunedCollection->pdgId() == 443) {
	//	  				  gen_muonP_p4.SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
	//	  				  foundit++;
	//	  				  j2 = motherInPrunedCollection;
	//	  			  }
	//	  			  if (dauId == 22) {
	//	  				  gen_photon_p4.SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
	//	  				  foundit++;
	//	  			  }
	//	  		  }
	//	  		  if (foundit == 4) break;
	//	  	  }
	//	  	  if (foundit == 4) {
	//	  		  if (j1 != nullptr && j2 != nullptr && j1 == j2) {
	//	  			  gen_jpsi_p4.SetPtEtaPhiM(j1->pt(), j1->eta(), j1->phi(), j1->mass());
	//	  			  break;
	//	  		  }
	//	  		  else {
	//	  			  std::cout << "Mother of muons does not match (" << j1->pdgId() << "," << j2->pdgId() << ")" << std::endl;
	//	  			  foundit = 0;
	//	  			  chic_pdgId = 0;
	//	  		  }
	//	  	  }
	//	  	  else {
	//	  		  std::cout << "Found just " << foundit << " out of 4 particles" << std::endl;
	//	  		  foundit = 0;
	//	  		  chic_pdgId = 0;
	//	  	  }
	//	    }  // if ( p_id
	//    } // for (size

	//    if (!chic_pdgId) { // sanity check
	//	    std::cout << "Rootupler does not found the given decay " <<
	//	  	  iEvent.id().run() << "," << iEvent.id().event() << std::endl;
	//    }
	//}


	//grab Trigger informations
	// save it in variable trigger, trigger is an int between 0 and 15, in binary it is:
	// (pass 11)(pass 8)(pass 7)(pass 5)
	// es. 11 = pass 5, 7 and 11
	// es. 4 = pass only 8

	trigger = 0;
	if (triggerResults_handle.isValid()) {

		const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);

		unsigned int NTRIGGERS = 9;
		//  string TriggersToTest[NTRIGGERS] = {"HLT_Dimuon10_Jpsi_Barrel","HLT_Dimuon16_Jpsi","HLT_Dimuon20_Jpsi","HLT_Dimuon8_Upsilon_Barrel","HLT_Dimuon13_Upsilon","HLT_Dimuon8_PsiPrime_Barrel","HLT_Dimuon13_PsiPrime","HLT_Mu16_TkMu0_dEta18_Onia","HLT_Mu25_TkMu0_dEta18_Onia"};
		string TriggersToTest[NTRIGGERS] = { "HLT_PAL1DoubleMuOpen","HLT_Dimuon16_Jpsi","HLT_PAL1DoubleMuOpen_OS","HLT_Dimuon8_Upsilon_Barrel","HLT_Dimuon13_Upsilon","HLT_PAL1DoubleMuOpen_SS","HLT_Dimuon13_PsiPrime","HLT_Mu16_TkMu0_dEta18_Onia","HLT_Mu25_TkMu0_dEta18_Onia" };

		for (unsigned int i = 0; i < NTRIGGERS; i++) {
			for (int version = 1; version < 5; version++) {
				stringstream ss;
				ss << TriggersToTest[i] << "_v" << version;
				unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label().c_str());
				if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
					trigger += (1 << i);
					break;
				}
			}
		}

	}
	else {
		std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;
	} // if (trigger...


	event_tree->Fill();
	Clear();
}

/////////////////////////////////////////////
///// G E N E R A L   F U N C T I O N S ///// 
/////////////////////////////////////////////

// ------------ method called once each job just before starting event loop  ------------
void ChiRootupler::beginJob() 
{

}

// ------------ method called once each job just after ending the event loop  ------------
void ChiRootupler::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void ChiRootupler::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void ChiRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void ChiRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void ChiRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ clear the vectors after each event
void ChiRootupler::Clear()
{
	//PV

	pvtx_z.clear();
	pvtx_x.clear();
	pvtx_y.clear();
	pvtx_nTracks.clear();
	pvtx_isFake.clear();


	//muon
	muonIsGlobal.clear();
	muonIsTracker.clear();
	muonIsPF.clear();
	muonIsSoft.clear();
	muonIsTight.clear();
	muonIsNotGlobalNorTracker.clear();
	muonIDHas_TMOneStationTight.clear();
	muonInnerTrack_dxy.clear();
	muonInnerTrack_dz.clear();
	muonTrackerLayersWithMeasurement.clear();
	muonPixelLayersWithMeasurement.clear();
	muonQuality_isHighPurity.clear();
	muon_charge.clear();
	muon_eta.clear();
	muon_pt.clear();
	muon_p4->Clear(); //TLorentzVector
	patMuonStored.clear();

	//muon MC
	muon_isMatchedMC.clear();
	muonGen_eta.clear();
	muonGen_pt.clear();
	muonGen_p4->Clear(); //TLorentzVector
	muonGen_rDelta.clear();
	muonGen_ptDelta.clear();
	muonGen_ptDeltaRel.clear();
	
	//dimuon 
	dimuon_p4->Clear();
	dimuon_eta.clear();
	dimuon_pt.clear();
	dimuon_charge.clear(); 
	dimuon_vtx->Clear();
	dimuon_dz_dimuonvtx_pvtx.clear();
	dimuonStored.clear();
	dimuon_muon1_position.clear();
	dimuon_muon2_position.clear();
	dimuon_ctpv.clear();
	dimuon_ctpvError.clear();

	//conv
	convQuality_isHighPurity.clear();
	convQuality_isGeneralTracksOnly.clear();
	conv_vtx->Clear();
	conv_vertexPositionRho.clear();
	conv_sigmaTkVtx1.clear();
	conv_sigmaTkVtx2.clear();
	conv_tkVtxCompatibilityOK.clear();
	conv_tkVtxCompatible_bestVertex.clear();
	conv_tkVtxCompatible_secondBestVertexA.clear();
	conv_tkVtxCompatible_secondBestVertexB.clear();
	conv_tkVtxCompatibilityOK_test.clear(); 
	conv_tkVtxCompatible_bestVertex_test.clear(); 
	conv_tkVtxCompatible_secondBestVertexA_test.clear(); 
	conv_tkVtxCompatible_secondBestVertexB_test.clear(); 

	conv_compatibleInnerHitsOK.clear();
	conv_hitPat1.clear();
	conv_hitPat2.clear();
	conv_isCustomHighPurity.clear();
	conv_zOfPriVtx.clear();
	conv_zOfPriVtxFromTracks.clear();
	conv_dzToClosestPriVtx.clear();
	conv_dxyPriVtx_Tr1.clear();
	conv_dxyPriVtx_Tr2.clear();
	conv_dxyPriVtxTimesCharge_Tr1.clear();
	conv_dxyPriVtxTimesCharge_Tr2.clear();
	conv_dxyError_Tr1.clear();
	conv_dxyError_Tr2.clear();

	conv_tk1NumOfDOF.clear();
	conv_tk2NumOfDOF.clear();
	conv_track1Chi2.clear();
	conv_track2Chi2.clear();
	conv_vertexChi2Prob.clear();
	conv_minDistanceOfApproach.clear();
	conv_p4->Clear(); 
	conv_eta.clear();
	conv_pt.clear();


	//conv MC
	conv_isMatchedMC.clear();
	convGen_eta.clear();
	convGen_pt.clear();
	convGen_p4->Clear();
	convGen_rDelta.clear();
	convGen_ptDelta.clear();
	convGen_ptDeltaRel.clear();
	convGen_motherCode.clear();

	//MC
	gen_pdgId.clear();
	gen_chic_pt.clear();
	gen_chic_eta.clear();
	gen_chic_p4->Clear();
	gen_Jpsi_pt.clear();
	gen_Jpsi_eta.clear();
	gen_Jpsi_p4->Clear();
	gen_Jpsi_matchPosition.clear();
	gen_Jpsi_nMatches.clear();
	gen_Jpsi_rDelta.clear();
	gen_Jpsi_ptDeltaRel.clear();
	gen_muon_charge.clear();
	gen_muon_pt.clear();
	gen_muon_eta.clear();
	gen_muon_p4->Clear();
	gen_muon_matchPosition.clear();
	gen_muon_nMatches.clear();
	gen_muon_rDelta.clear();
	gen_muon_ptDeltaRel.clear();
	gen_phot_pt.clear();
	gen_phot_eta.clear();
	gen_phot_p4->Clear();
	gen_conv_matchPosition.clear();
	gen_conv_nMatches.clear();
	gen_conv_rDelta.clear();
	gen_conv_ptDeltaRel.clear();

	//chi
	chi_p4->Clear();
	chi_eta.clear();
	chi_pt.clear();
	chi_daughterJpsi_position.clear();
	chi_daughterConv_position.clear();
	chi_dzPhotToDimuonVtx.clear();
	chi_dxyPhotToDimuonVtx.clear();
	chiStored.clear();

}


//////////////////////////////////////////////
///// P H O T O N   C H E C K S    //////////  based on HeavyFlavorAnalysis.Onia2MuMu.OniaPhotonConversionProducer
////////////////////////////////////////////

bool ChiRootupler::lt_comparator(std::pair<double, short> a, std::pair<double, short> b) {
	return a.first < b.first;
}

bool ChiRootupler::Conv_checkTkVtxCompatibility(const reco::Conversion& conv, const reco::VertexCollection& priVtxs, double sigmaTkVtxComp_, bool& Flag_Best_Out, bool& Flag_SecondBestA_Out, bool& Flag_SecondBestB_Out, double& sigmaMinValue1Out, double& sigmaMinValue2Out)
{
	sigmaMinValue1Out = 100; //starting values
	sigmaMinValue2Out = 100;
	Flag_Best_Out = false;
	Flag_SecondBestA_Out = false;
	Flag_SecondBestB_Out = false;

	std::vector< std::pair< double, short> > idx[2];  //dz distance, vertex number
	short ik = -1;
	BOOST_FOREACH(edm::RefToBase<reco::Track> tk, conv.tracks()) {  //ik - going over conversion tracks
		ik++;
		short count = -1;
		BOOST_FOREACH(const reco::Vertex& vtx, priVtxs) {  //count - going over vertices
			count++;
			double dz_ = tk->dz(vtx.position());
			double dzError_ = tk->dzError();
			dzError_ = sqrt(dzError_*dzError_ + vtx.covariance(2, 2));

			if ((ik==0)&&((fabs(dz_) / dzError_) < sigmaMinValue1Out)) {//save lowest value of sigma
				sigmaMinValue1Out = fabs(dz_) / dzError_;
			}
			if ((ik == 1) && ((fabs(dz_) / dzError_) < sigmaMinValue2Out)) {
				sigmaMinValue2Out = fabs(dz_) / dzError_;
			}

			if (fabs(dz_) / dzError_ > sigmaTkVtxComp_) continue;

			idx[ik].push_back(std::pair<double, short>(fabs(dz_), count));
		}
		if (idx[ik].size() == 0) { return false; }

		std::stable_sort(idx[ik].begin(), idx[ik].end(), lt_comparator);
	}

	Flag_Best_Out = (idx[0][0].second == idx[1][0].second);
	Flag_SecondBestA_Out = (idx[0][1].second == idx[1][0].second); //second best vertex is OK
	Flag_SecondBestB_Out = (idx[0][0].second == idx[1][1].second);

	if (idx[0][0].second == idx[1][0].second || idx[0][1].second == idx[1][0].second || idx[0][0].second == idx[1][1].second) return true;  //original return
	return false;
}

bool ChiRootupler::Conv_foundCompatibleInnerHits(const reco::HitPattern& hitPatA, const reco::HitPattern& hitPatB) { //directly copied from HeavyFlavorAnalysis.Onia2MuMu.OniaPhotonConversionProducer
	size_t count = 0;
	uint32_t oldSubStr = 0;
	for (int i = 0; i<hitPatA.numberOfHits(reco::HitPattern::HitCategory::TRACK_HITS) && count<2; i++) {
		uint32_t hitA = hitPatA.getHitPattern(reco::HitPattern::HitCategory::TRACK_HITS, i);
		if (!hitPatA.validHitFilter(hitA) || !hitPatA.trackerHitFilter(hitA)) continue;

		if (hitPatA.getSubStructure(hitA) == oldSubStr && hitPatA.getLayer(hitA) == oldSubStr)
			continue;

		if (hitPatB.getTrackerMonoStereo(reco::HitPattern::HitCategory::TRACK_HITS, hitPatA.getSubStructure(hitA), hitPatA.getLayer(hitA)) != 0)
			return true;

		oldSubStr = hitPatA.getSubStructure(hitA);
		count++;
	}
	return false;
}

//bool OniaPhotonConversionProducer::HighpuritySubset(const reco::Conversion& conv, const reco::VertexCollection& priVtxs) {	// select high purity conversions same way as OniaPhotonConversionProducer:
//	// vertex chi2 cut
//	if (ChiSquaredProbability(conv.conversionVertex().chi2(), conv.conversionVertex().ndof())< _vertexChi2ProbCut) return false;
//
//	// d0 cut
//	// Find closest primary vertex
//	int closest_pv_index = 0;
//	int i = 0;
//	BOOST_FOREACH(const reco::Vertex& vtx, priVtxs) {
//		if (conv.zOfPrimaryVertexFromTracks(vtx.position()) < conv.zOfPrimaryVertexFromTracks(priVtxs[closest_pv_index].position())) closest_pv_index = i;
//		i++;
//	}
//	// Now check impact parameter wtr with the just found closest primary vertex
//	BOOST_FOREACH(const edm::RefToBase<reco::Track> tk, conv.tracks()) if (-tk->dxy(priVtxs[closest_pv_index].position())*tk->charge() / tk->dxyError()<0) return false;
//
//	// chi2 of single tracks
//	BOOST_FOREACH(const edm::RefToBase<reco::Track> tk, conv.tracks()) if (tk->normalizedChi2() > _trackchi2Cut) return false;
//
//	// dof for each track  
//	BOOST_FOREACH(const edm::RefToBase<reco::Track> tk, conv.tracks()) if (tk->ndof()< TkMinNumOfDOF_) return false;
//
//	// distance of approach cut
//	if (conv.distOfMinimumApproach() < _minDistanceOfApproachMinCut || conv.distOfMinimumApproach() > _minDistanceOfApproachMaxCut) return false;
//
//	return true;
//}

bool ChiRootupler::Conv_isMatched(const math::XYZTLorentzVectorF& reco_conv, const reco::GenParticle& gen_phot, double maxDeltaR, double maxDPtRel)
{
	return reco::deltaR(reco_conv, gen_phot) < maxDeltaR && ((reco_conv.pt() - gen_phot.pt()) / (gen_phot.pt()+ 1E-9)) < maxDPtRel;
}



//////////////////////////////////
////  C H I   P A R T S   ////
//////////////////////////////


const pat::CompositeCandidate ChiRootupler::makeChiCandidate(const pat::CompositeCandidate& dimuon,	const pat::CompositeCandidate& photon) {
	pat::CompositeCandidate chiCand;
	chiCand.addDaughter(dimuon, "dimuon");
	chiCand.addDaughter(photon, "photon");
	chiCand.setVertex(dimuon.vertex());
	reco::Candidate::LorentzVector vChic = dimuon.p4() + photon.p4();
	chiCand.setP4(vChic);
	return chiCand;
}

double ChiRootupler::Getdz(const pat::CompositeCandidate& c, const reco::Candidate::Point &p) {
	reco::Candidate::LorentzVector mom = c.p4();
	reco::Candidate::Point vtx = c.vertex();
	double dz = (vtx.Z() - p.Z()) - ((vtx.X() - p.X())*mom.X() + (vtx.Y() - p.Y())*mom.Y()) / mom.Rho() * mom.Z() / mom.Rho();
	return dz;
}

double ChiRootupler::Getdxy(const pat::CompositeCandidate& c, const reco::Candidate::Point &p) {
	reco::Candidate::LorentzVector mom = c.p4();
	reco::Candidate::Point vtx = c.vertex();
	double dx = (vtx.X() - p.X()) - ((vtx.Z() - p.Z()) * mom.X() / mom.Z());
	double dy = (vtx.Y() - p.Y()) - ((vtx.Z() - p.Z()) * mom.Y() / mom.Z());
	double dxy = sqrt(dx*dx+dy*dy);
	return dxy;
}


///////////////////////////////////
////  G E N   P A R T S  //////////
//////////////////////////////////////

template <typename particle_in, typename particle_type>
int ChiRootupler::MatchGen(particle_in& myGenParticle, edm::Handle< std::vector <particle_type>>& collToBeMatched, double maxDeltaR, double maxDPtRel, int& nMatches_Out, double& rDelta_Out, double& ptDeltaRel_Out) //returns position of the best match
{
	nMatches_Out = 0;
	rDelta_Out = -5;
	ptDeltaRel_Out = -5;
	bool isMatched = false;
	int bestMatchPosition = -1;

	if (collToBeMatched.isValid())
	{
		for (uint j = 0; j < collToBeMatched->size(); j++) { //loop over the collection, find the best match
			const particle_type& myParticle = collToBeMatched->at(j);
			const reco::Candidate& recoParticle = GetRecoCandidate(myParticle); //get reco candidate of whatever type went in
			const reco::Candidate& genParticle = GetRecoCandidate(myGenParticle);
			double deltaR = reco::deltaR(recoParticle, genParticle);// reco::deltaR(inCand.eta(), inCand.phi(), mCand.eta(), mCand.phi());
			double dPtRel = std::fabs(recoParticle.pt() - genParticle.pt()) / (double)(genParticle.pt() + 1E-9); //defined in the matcher to be divided by genPt //MCTruthMatchers.cc
			//std::cout << "pt reco "<< recoParticle.pt()<<" pt gen "<< genParticle.pt()<<"  and relative "<<dPtRel << "   difference "<< std::fabs(recoParticle.pt() - genParticle.pt())<<"  and denom "<< (genParticle.pt() + 1E-9) <<std::endl;
			
			if ((deltaR < maxDeltaR) && (dPtRel < maxDPtRel) && (recoParticle.charge() == genParticle.charge())) //this particular is matched
			{
				if (isMatched == false)//first match
				{
					isMatched = true;
					nMatches_Out++;
					bestMatchPosition = j;
					rDelta_Out = deltaR;
					ptDeltaRel_Out = dPtRel;
				}
				else //it has been matched already
				{
					nMatches_Out++;
					//compare with the best match by R, save the better one
					const reco::Candidate& prevMatch = GetRecoCandidate(collToBeMatched->at(bestMatchPosition));
					if (deltaR < reco::deltaR(prevMatch, genParticle)) {
						bestMatchPosition = j;
						rDelta_Out = deltaR;
						ptDeltaRel_Out = dPtRel;
					}
				}

			}

		}
		return bestMatchPosition;
	}
	else {
		std::cout << "Passed handle invalid" << std::endl;
		return -2; //-2 is a problem
	}
}


template<typename T>  reco::LeafCandidate ChiRootupler::GetRecoCandidate(const T& a) { return (reco::LeafCandidate) a; }

reco::LeafCandidate ChiRootupler::GetRecoCandidate(const reco::Conversion& conv) 
{
	reco::Candidate::LorentzVector p4;
	if (conv.conversionVertex().isValid()) {
		p4.SetPxPyPzE(conv.refittedPair4Momentum().Px(), conv.refittedPair4Momentum().Py(), conv.refittedPair4Momentum().Pz(), conv.refittedPair4Momentum().E());
	}
	else { 
		p4.SetPxPyPzE(conv.pairMomentum().X(), conv.pairMomentum().Y(), conv.pairMomentum().Z(), conv.pairMomentum().R()); //assumes no mass
	}
	reco::LeafCandidate out;
	out.setP4(p4);
	out.setCharge(0.0);
	return out;
}

reco::LeafCandidate ChiRootupler::GetRecoCandidate(const pat::CompositeCandidate& comp)
{
	//const reco::Candidate cand = (PATObject<reco::CompositeCandidate>(comp)).originalObject();
	//const reco::Candidate cand2 = dynamic_cast<const reco::Candidate> *cand;
	//return (reco::LeafCandidate) cand;

	reco::Candidate::LorentzVector p4 = comp.p4();
	reco::LeafCandidate out;
	out.setP4(p4);
	out.setCharge(comp.charge());
	return out;
}


reco::LeafCandidate ChiRootupler::GetRecoCandidate(const pat::Muon& comp)
{
	//const reco::Candidate cand = (PATObject<reco::CompositeCandidate>(comp)).originalObject();
	//const reco::Candidate cand2 = dynamic_cast<const reco::Candidate> *cand;
	//return (reco::LeafCandidate) cand;

	reco::Candidate::LorentzVector p4 = comp.p4();
	reco::LeafCandidate out;
	out.setP4(p4);
	out.setCharge(comp.charge());
	return out;
}











////////////////////////////////////////////////

//define this as a plug-in
DEFINE_FWK_MODULE(ChiRootupler);
