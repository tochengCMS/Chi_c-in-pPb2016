#ifndef __ChiRootupler_h_
#define __ChiRootupler_h_

  // Declaration of ChiRootupler
  // Description: Saves the muon, dimuon and chi candidate information
  // Implementation:  Ota Kukral based on work of Andre Stahl and Stefano Argiro, Alessandro Degano  and the Torino team, Alberto Sanchez

//O: includes should be cleaned up and moved to cc
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"
#include <TClonesArray.h>
#include <vector>
#include <sstream>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"


//
// class declaration
//

class ChiRootupler :public edm::EDAnalyzer {
public:
	explicit ChiRootupler(const edm::ParameterSet &);
	~ChiRootupler();
	//static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
	//bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);

private:
	
	virtual void analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup);						// ------------ method called for each event  ------------
		
	virtual void beginJob();																				// ------------ method called once each job just before starting event loop  ------------
	virtual void endJob();																					// ------------ method called once each job just after ending the event loop  ------------
	virtual void beginRun(edm::Run const &, edm::EventSetup const &);										// ------------ method called when starting to processes a run  ------------
	virtual void endRun(edm::Run const &, edm::EventSetup const &);											// ------------ method called when ending the processing of a run  ------------
	virtual void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);				// ------------ method called when starting to processes a luminosity block  ------------
	virtual void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);					// ------------ method called when ending the processing of a luminosity block  ------------

	// general functions
	void Clear();

	// photon checks
	static bool lt_comparator(std::pair<double, short> a, std::pair<double, short> b); // comparator for checking conversions
	bool Conv_checkTkVtxCompatibility(const reco::Conversion& conv, const reco::VertexCollection&  priVtxs, double sigmaTkVtxComp_, bool& Flag_Best_Out, bool& Flag_SecondBestA_Out, bool& Flag_SecondBestB_Out, double& sigmaMinValue1Out, double& sigmaMinValue2Out);
	bool Conv_foundCompatibleInnerHits(const reco::HitPattern& hitPatA, const reco::HitPattern& hitPatB);
	//photon MC
	bool Conv_isMatched(const math::XYZTLorentzVectorF& reco_conv, const reco::GenParticle& gen_phot, double maxDeltaR, double maxDPtRel);

	// chi functions
	const pat::CompositeCandidate makeChiCandidate(const pat::CompositeCandidate&, const pat::CompositeCandidate&);
	double Getdz(const pat::CompositeCandidate&, const reco::Candidate::Point &);

	// gen functions
	template <typename particle_in, typename particle_type>
	int MatchGen(particle_in& genParticle, edm::Handle< std::vector <particle_type>>& collToBeMatched, double maxDeltaR, double maxDPtRel, int& nMatches_Out, double& rDelta_Out, double& ptDeltaRel_Out);
	template <typename T>  reco::LeafCandidate GetRecoCandidate(const T& a); //default, however momentum is stored differently for conversion and pats
	reco::LeafCandidate GetRecoCandidate(const reco::Conversion& conv); //momentum is stored differently for conversion
	reco::LeafCandidate GetRecoCandidate(const pat::CompositeCandidate& comp);
	reco::LeafCandidate GetRecoCandidate(const pat::Muon& comp);

	//template <typename particle_type>
	//int MatchGen(particle_type& genParticle, double maxDeltaR, double maxDPtRel) { return 0; }

	//int MatchGen(const reco::Candidate& genParticle, double maxDeltaR, double maxDPtRel) { return 0; }

	// ----------member data ---------------------------
	std::string file_name;

	//edm::EDGetTokenT< edm::View <pat::Muon> > muon_label; //is a muon collection
	edm::EDGetTokenT<pat::MuonCollection>  muon_label;
	edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_label;
	edm::EDGetTokenT<pat::CompositeCandidateCollection> photon_label;
	edm::EDGetTokenT<reco::ConversionCollection> conversion_label;
	//edm::EDGetTokenT<pat::CompositeCandidateCollection> chi_label;
	edm::EDGetTokenT<reco::VertexCollection> primaryVertices_label;
	edm::EDGetTokenT<edm::TriggerResults> triggerResults_label;
	edm::EDGetTokenT<reco::GenParticleCollection> genParticles_label;

	bool flag_doMC;
	bool flag_saveExtraThings = false;  //saves all muons, and other not so necessary data
	
	const int PythCode_chic0 = 10441; //Pythia codes
	const int PythCode_chic1 = 20443;
	const int PythCode_chic2 = 445;

	// constants for muon cuts
	const double muon_maxDeltaR = 0.5;
	const double muon_maxDPtRel = 0.5;

	// constants for dimuon cuts
	const double jpsi_maxDeltaR = 0.5;
	const double jpsi_maxDPtRel = 0.5;

	// constants for conversion cuts
	const double conv_TkVtxCompSigmaCut = 50.0;
	const double conv_maxDeltaR = 0.2;
	const double conv_maxDPtRel = 1;



	TTree* dimuon_tree;
	TTree* chi_tree;

	TTree* event_tree;

	//general
	long runNumber;
	long eventNumber;
	long nPrimVertices;
	int muonPerEvent;
	int convPerTriggeredEvent;
	

	//muon info
	std::vector <bool> muonIsGlobal;
	std::vector <bool> muonIsTracker;
	std::vector <bool> muonIsPF;
	std::vector <bool> muonIsNotGlobalNorTracker;
	std::vector <bool> muonIDHas_TMOneStationTight;
	std::vector <double> muonInnerTrack_dxy;
	std::vector <double> muonInnerTrack_dz;
	std::vector <int> muonTrackerLayersWithMeasurement;
	std::vector <int> muonPixelLayersWithMeasurement;
	std::vector <bool> muonQuality_isHighPurity;
	std::vector <int> muon_charge;
	std::vector <double> muon_eta;
	std::vector <double> muon_pt;
	TClonesArray* muon_p4; //TLorentzVector
	std::vector <pat::Muon> patMuonStored;

	//muon MC
	std::vector <bool> muon_isMatchedMC;
	std::vector <double> muonGen_eta;
	std::vector <double> muonGen_pt;
	TClonesArray* muonGen_p4; //TLorentzVector
	std::vector <double> muonGen_rDelta;
	std::vector <double> muonGen_ptDelta;
	std::vector <double> muonGen_ptDeltaRel;


	//dimuon info

	TClonesArray*  dimuon_p4; //TLorentzVector
	std::vector <double> dimuon_eta;
	std::vector <double> dimuon_pt;
	std::vector <double> dimuon_charge; //crosscheck
	TClonesArray* dimuon_vtx; //TVector3
	std::vector <pat::CompositeCandidate> dimuonStored;
	std::vector <double> dimuon_pmuon_position; //stores position of positive muon in muon collection
	std::vector <double> dimuon_nmuon_position; //stores position of negative muon in muon collection




	//conversion info

	std::vector <bool> convQuality_isHighPurity;
	std::vector <bool> convQuality_isGeneralTracksOnly;
	TClonesArray* conv_vtx; //TVector3
	std::vector <double> conv_vertexPositionRho;
	std::vector <double> conv_sigmaTkVtx1;
	std::vector <double> conv_sigmaTkVtx2;
	std::vector <bool> conv_tkVtxCompatibilityOK;
	std::vector <bool> conv_tkVtxCompatible_bestVertex;
	std::vector <bool> conv_tkVtxCompatible_secondBestVertexA;
	std::vector <bool> conv_tkVtxCompatible_secondBestVertexB;
	std::vector <bool> conv_tkVtxCompatibilityOK_test; //test ones
	std::vector <bool> conv_tkVtxCompatible_bestVertex_test; //test ones
	std::vector <bool> conv_tkVtxCompatible_secondBestVertexA_test; //test
	std::vector <bool> conv_tkVtxCompatible_secondBestVertexB_test; //test

	std::vector <int> conv_compatibleInnerHitsOK; //-1: less than 2 tracks, 0: not compatible, 1: yes
	std::vector <reco::HitPattern> conv_hitPat1;
	std::vector <reco::HitPattern> conv_hitPat2;
	std::vector <bool> conv_isCustomHighPurity;//tbd - is just a sum of some other cuts, not creating at the time
	std::vector <double> conv_vertexChi2Prob;
	std::vector <double> conv_zOfPriVtx;
	std::vector <double> conv_zOfPriVtxFromTracks;
	std::vector <double> conv_dzToClosestPriVtx;
	std::vector <double> conv_dxyPriVtx_Tr1;
	std::vector <double> conv_dxyPriVtx_Tr2;
	std::vector <double> conv_dxyPriVtxTimesCharge_Tr1;
	std::vector <double> conv_dxyPriVtxTimesCharge_Tr2;
	std::vector <double> conv_dxyError_Tr1;
	std::vector <double> conv_dxyError_Tr2;

	std::vector <int> conv_tk1NumOfDOF;
	std::vector <int> conv_tk2NumOfDOF;
	std::vector <double> conv_track1Chi2;
	std::vector <double> conv_track2Chi2;

	std::vector <double> conv_minDistanceOfApproach;
	TClonesArray*  conv_p4; //TLorentzVector
	std::vector <double> conv_eta;
	std::vector <double> conv_pt;


	//conv MC
	std::vector <bool> conv_isMatchedMC;
	std::vector <double> convGen_eta;
	std::vector <double> convGen_pt;
	TClonesArray* convGen_p4; //TLorentzVector
	std::vector <double> convGen_rDelta;
	std::vector <double> convGen_ptDelta;
	std::vector <double> convGen_ptDeltaRel;
	std::vector <int> convGen_motherCode;
	

	// MC general

	std::vector <int> gen_pdgId;
	std::vector <double> gen_chic_pt;
	std::vector <double> gen_chic_eta;
	TClonesArray* gen_chic_p4; //TLorentzVector
	std::vector <double> gen_Jpsi_pt;
	std::vector <double> gen_Jpsi_eta;
	std::vector <int> gen_Jpsi_matchPosition;
	std::vector <int> gen_Jpsi_nMatches;
	std::vector <double> gen_Jpsi_rDelta; //in principle duplicates information
	std::vector <double> gen_Jpsi_ptDeltaRel;//in principle duplicates information
	TClonesArray* gen_Jpsi_p4; //TLorentzVector

	std::vector <int> gen_muon_charge;
	std::vector <double> gen_muon_pt;
	std::vector <double> gen_muon_eta;
	std::vector <int> gen_muon_matchPosition;
	std::vector <int> gen_muon_nMatches;
	std::vector <double> gen_muon_rDelta; //in principle duplicates information
	std::vector <double> gen_muon_ptDeltaRel;//in principle duplicates information
	TClonesArray* gen_muon_p4; //TLorentzVector

	std::vector <double> gen_phot_pt;
	std::vector <double> gen_phot_eta;
	TClonesArray* gen_phot_p4; //TLorentzVector
	std::vector <int> gen_conv_matchPosition;
	std::vector <int> gen_conv_nMatches;
	std::vector <double> gen_conv_rDelta; //in principle duplicates information
	std::vector <double> gen_conv_ptDeltaRel;//in principle duplicates information


	// chi

	TClonesArray*  chi_p4; //TLorentzVector
	std::vector <double> chi_eta;
	std::vector <double> chi_pt;
	pat::CompositeCandidate chi_cand;
	std::vector <double> chi_daughterJpsi_position; //stores position of daughter Jpsi in dimuon collection
	std::vector <double> chi_daughterConv_position; //stores position of daughter photon (conversion)
	std::vector <double> chi_dzPhotToDimuonVtx;


	// run and vertex info
	//TBU
	TVector3 primary_v;
	TVector3 secondary_v;
	TVector3 dimuon_v;

	// chi info


	//double chi_dzPhotToDimuonVtx;

	// Lorentz vectors
	//TLorentzVector chi_p4;
	TLorentzVector chi_dimuon_p4;
	TLorentzVector muonP_p4;
	TLorentzVector muonN_p4;
	TLorentzVector photon_p4;

	//Various
	Double_t ele_lowerPt_pt;
	Double_t ele_higherPt_pt;
	Double_t ctpv;
	Double_t ctpv_error;
	Double_t pi0_abs_mass;
	Double_t psi1S_nsigma;
	Double_t psi2S_nsigma;
	Double_t dz;
	Int_t trigger;


	//MC
	TLorentzVector MC_chic_p4;
	Int_t chic_pdgId;
	TLorentzVector gen_jpsi_p4;
	TLorentzVector gen_photon_p4;
	TLorentzVector gen_muonP_p4;
	TLorentzVector gen_muonM_p4;
	edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
	edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;

	// static data member definitions
	//
	const double pi0_mass = 0.1349766;
	const Double_t psi1SMass = 3.09691;
	const Double_t psi2SMass = 3.68610;

	/*
	// 2011 par
	static const double Y_sig_par_A = 0.058;
	static const double Y_sig_par_B = 0.047;
	static const double Y_sig_par_C = 0.22;
	*/

	// 2012 par
	const double Y_sig_par_A = 62.62;
	const double Y_sig_par_B = 56.3;
	const double Y_sig_par_C = -20.77;


};

#endif 

