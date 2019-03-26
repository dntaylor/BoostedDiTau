// -*- C++ -*-
//
// Package:    MuonCleanedMiniAODJetProducer
// Class:      MuonCleanedMiniAODJetProducer
// 
/**\class MuonCleanedMiniAODJetProducer MuonCleanedMiniAODJetProducer.cc

 Description: Removes PF muons from PFJet candidates and reconstructs the jets
	          Associates those muons to the jets from which they were removed

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesca Ricci-Tam,6 R-025,+41227672274,
//     Contributer:  Devin Taylor
//         Created:  Fri Aug 31 13:01:48 CEST 2012
//
//


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "TLorentzVector.h"
#include "TMath.h"

//
// class declaration
//

class MuonCleanedMiniAODJetProducer : public edm::stream::EDProducer<>
{
   public:
      explicit MuonCleanedMiniAODJetProducer(const edm::ParameterSet&);
      ~MuonCleanedMiniAODJetProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------

      // source of the jets to be cleaned of muons
      edm::EDGetTokenT<reco::PFJetCollection> jetSrc_;

      // source of muons that, if found within jet, should be removed
      edm::EDGetTokenT<edm::RefVector<pat::MuonCollection> > muonSrc_;

      // source of PF candidates
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandSrc_;

      edm::ParameterSet* cfg_;

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
MuonCleanedMiniAODJetProducer::MuonCleanedMiniAODJetProducer(const edm::ParameterSet& iConfig):
  jetSrc_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetSrc"))),
  muonSrc_(consumes<edm::RefVector<pat::MuonCollection> >(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  pfCandSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandSrc")))
{
  cfg_ = const_cast<edm::ParameterSet*>(&iConfig);

  //register your products
  produces<reco::PFJetCollection>();
  produces<edm::ValueMap<bool> >( "jetCleanedValueMap" );
  produces<pat::PackedCandidateCollection >( "particleFlowMuonCleaned" );

}


MuonCleanedMiniAODJetProducer::~MuonCleanedMiniAODJetProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
MuonCleanedMiniAODJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::PFJetCollection> pfJets;
  iEvent.getByToken(jetSrc_, pfJets);
  std::unique_ptr<reco::PFJetCollection> SetOfJets( new reco::PFJetCollection );

  edm::Handle<edm::RefVector<pat::MuonCollection> > muons;
  iEvent.getByToken(muonSrc_, muons);

  edm::Handle<pat::PackedCandidateCollection> pfCands;
  iEvent.getByToken(pfCandSrc_, pfCands);
  std::unique_ptr<pat::PackedCandidateCollection > pfCandsExcludingMuons(new pat::PackedCandidateCollection);

  std::vector<reco::CandidatePtr> muonPFs;
  if (muons.isValid()) 
  {
    for (edm::RefVector<pat::MuonCollection>::const_iterator iMuon = muons->begin(); iMuon != muons->end(); ++iMuon)
    {
      for (unsigned int j = 0; j<(*iMuon)->numberOfSourceCandidatePtrs(); j++) {
        muonPFs.push_back((*iMuon)->sourceCandidatePtr(j));
      }
    }
  }

  // fill the vector of particleFlow without muons
  // TODO: error when active
  // Inconsistency RefCore::pushBackItem: Ref or Ptr is inconsistent with RefVector (PtrVector)id = (4:2072) should be (6:28)
  for (unsigned int i = 0; i != pfCands->size(); i++) {
    const pat::PackedCandidate &pf = (*pfCands)[i];
    //pat::PackedCandidateRef pfRef(pfCands,i);
    if (std::find(muonPFs.begin(), muonPFs.end(), reco::CandidatePtr(pfCands,i)) == muonPFs.end()) {
      pfCandsExcludingMuons->push_back(pf);
    }
  }

  //vector of bools holding the signal muon tag decision for each jet
  std::vector<bool> muonTagDecisions;

  // Do cleaning
  for (reco::PFJetCollection::const_iterator iJet = pfJets->begin(); iJet != pfJets->end(); ++iJet)
  {
    std::vector<reco::CandidatePtr> jetPFCands = iJet->daughterPtrVector();
    reco::PFJet::Specific specs = iJet->getSpecific();
    math::XYZTLorentzVector pfmomentum;
    std::vector<edm::Ptr<reco::Candidate> > jetConstituents;
    jetConstituents.clear();

    //flag indicating whether >=0 muons were tagged for removal
    bool taggedMuonForRemoval = false;

    for (std::vector<reco::CandidatePtr>::iterator i = jetPFCands.begin(); i != jetPFCands.end(); ++i)
    {
      reco::CandidatePtr pfCand = *i;
      
      //does this muon pass the desired muon ID?
      if (std::find(muonPFs.begin(), muonPFs.end(), pfCand) != muonPFs.end())
   	  {
        specs.mMuonEnergy -= pfCand->p4().e();
        specs.mMuonMultiplicity -= 1;
        specs.mChargedMuEnergy -= pfCand->p4().e();
        specs.mChargedMultiplicity -= 1;
        //save tag decision for this muon
        taggedMuonForRemoval = true;

      }
      else
      {
        pfmomentum += pfCand->p4(); // total p4()
        jetConstituents.push_back(pfCand);
      }
    } // loop over PF candidates

    // Build a new jet without the muon
    reco::PFJet muonfreePFJet(pfmomentum, specs, jetConstituents);
    SetOfJets->push_back( muonfreePFJet );

    //if at least 1 muon was tagged for removal, save a positive muon tag decision for this jet
    muonTagDecisions.push_back(taggedMuonForRemoval);


  } // loop over jets
  
  const edm::OrphanHandle<reco::PFJetCollection> cleanedJetsRefProd = iEvent.put(std::move(SetOfJets));

  //fill the value map of muon tag decision for each cleaned jet
  std::unique_ptr<edm::ValueMap<bool> > valMap(new edm::ValueMap<bool>());
  edm::ValueMap<bool>::Filler filler(*valMap);
  filler.insert(cleanedJetsRefProd, muonTagDecisions.begin(), muonTagDecisions.end());
  filler.fill();
  iEvent.put(std::move(valMap), "jetCleanedValueMap" );

  //put the soft-muon-free PF cands into the event
  iEvent.put(std::move(pfCandsExcludingMuons), "particleFlowMuonCleaned");

}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonCleanedMiniAODJetProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonCleanedMiniAODJetProducer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonCleanedMiniAODJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonCleanedMiniAODJetProducer);
