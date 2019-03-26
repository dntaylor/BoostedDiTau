// -*- C++ -*-
//
// Package:    ElectronCleanedMiniAODJetProducer
// Class:      ElectronCleanedMiniAODJetProducer
// 
/**\class ElectronCleanedMiniAODJetProducer ElectronCleanedMiniAODJetProducer.cc

 Description: Removes PF electrons from PFJet candidates and reconstructs the jets
	          Associates those electrons to the jets from which they were removed

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
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "TLorentzVector.h"
#include "TMath.h"

//
// class declaration
//

class ElectronCleanedMiniAODJetProducer : public edm::stream::EDProducer<>
{
   public:
      explicit ElectronCleanedMiniAODJetProducer(const edm::ParameterSet&);
      ~ElectronCleanedMiniAODJetProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------

      // source of the jets to be cleaned of electrons
      edm::EDGetTokenT<reco::PFJetCollection> jetSrc_;

      // source of electrons that, if found within jet, should be removed
      edm::EDGetTokenT<pat::ElectronCollection> electronSrc_;

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
ElectronCleanedMiniAODJetProducer::ElectronCleanedMiniAODJetProducer(const edm::ParameterSet& iConfig):
  jetSrc_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetSrc"))),
  electronSrc_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSrc"))),
  pfCandSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandSrc")))
{
  cfg_ = const_cast<edm::ParameterSet*>(&iConfig);

  //register your products
  produces<reco::PFJetCollection>();
  produces<edm::ValueMap<bool> >( "jetCleanedValueMap" );
  produces<pat::PackedCandidateCollection >( "particleFlowElectronCleaned" );

}


ElectronCleanedMiniAODJetProducer::~ElectronCleanedMiniAODJetProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ElectronCleanedMiniAODJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::PFJetCollection> pfJets;
  iEvent.getByToken(jetSrc_, pfJets);
  std::unique_ptr<reco::PFJetCollection> SetOfJets( new reco::PFJetCollection );

  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronSrc_, electrons);

  edm::Handle<pat::PackedCandidateCollection> pfCands;
  iEvent.getByToken(pfCandSrc_, pfCands);
  std::unique_ptr<pat::PackedCandidateCollection > pfCandsExcludingElectrons(new pat::PackedCandidateCollection);

  std::vector<reco::CandidatePtr> electronPFs;
  if (electrons.isValid()) 
  {
    for (pat::ElectronCollection::const_iterator iElectron = electrons->begin(); iElectron != electrons->end(); ++iElectron)
    {
      for (unsigned int j = 0; j<iElectron->numberOfSourceCandidatePtrs(); j++) {
        electronPFs.push_back(iElectron->sourceCandidatePtr(j));
      }
    }
  }

  // fill the vector of particleFlow without electrons
  // TODO: error when active
  // Inconsistency RefCore::pushBackItem: Ref or Ptr is inconsistent with RefVector (PtrVector)id = (4:2072) should be (6:28)
  for (unsigned int i = 0; i != pfCands->size(); i++) {
    const pat::PackedCandidate &pf = (*pfCands)[i];
    //pat::PackedCandidateRef pfRef(pfCands,i);
    if (std::find(electronPFs.begin(), electronPFs.end(), reco::CandidatePtr(pfCands,i)) == electronPFs.end()) {
      pfCandsExcludingElectrons->push_back(pf);
    }
  }

  //vector of bools holding the signal electron tag decision for each jet
  std::vector<bool> electronTagDecisions;

  // Do cleaning
  for (reco::PFJetCollection::const_iterator iJet = pfJets->begin(); iJet != pfJets->end(); ++iJet)
  {
    std::vector<reco::CandidatePtr> jetPFCands = iJet->daughterPtrVector();
    reco::PFJet::Specific specs = iJet->getSpecific();
    math::XYZTLorentzVector pfmomentum;
    std::vector<edm::Ptr<reco::Candidate> > jetConstituents;
    jetConstituents.clear();

    //flag indicating whether >=0 electrons were tagged for removal
    bool taggedElectronForRemoval = false;

    for (std::vector<reco::CandidatePtr>::iterator i = jetPFCands.begin(); i != jetPFCands.end(); ++i)
    {
      reco::CandidatePtr pfCand = *i;
      
      //does this electron pass the desired electron ID?
      if (std::find(electronPFs.begin(), electronPFs.end(), pfCand) != electronPFs.end())
   	  {
        specs.mElectronEnergy -= pfCand->p4().e();
        specs.mElectronMultiplicity -= 1;
        specs.mChargedMuEnergy -= pfCand->p4().e();
        specs.mChargedMultiplicity -= 1;
        //save tag decision for this electron
        taggedElectronForRemoval = true;

      }
      else
      {
        pfmomentum += pfCand->p4(); // total p4()
        jetConstituents.push_back(pfCand);
      }
    } // loop over PF candidates

    // Build a new jet without the electron
    reco::PFJet electronfreePFJet(pfmomentum, specs, jetConstituents);
    SetOfJets->push_back( electronfreePFJet );

    //if at least 1 electron was tagged for removal, save a positive electron tag decision for this jet
    electronTagDecisions.push_back(taggedElectronForRemoval);


  } // loop over jets
  
  const edm::OrphanHandle<reco::PFJetCollection> cleanedJetsRefProd = iEvent.put(std::move(SetOfJets));

  //fill the value map of electron tag decision for each cleaned jet
  std::unique_ptr<edm::ValueMap<bool> > valMap(new edm::ValueMap<bool>());
  edm::ValueMap<bool>::Filler filler(*valMap);
  filler.insert(cleanedJetsRefProd, electronTagDecisions.begin(), electronTagDecisions.end());
  filler.fill();
  iEvent.put(std::move(valMap), "jetCleanedValueMap" );

  //put the soft-electron-free PF cands into the event
  iEvent.put(std::move(pfCandsExcludingElectrons), "particleFlowElectronCleaned");

}

// ------------ method called once each job just before starting event loop  ------------
void 
ElectronCleanedMiniAODJetProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronCleanedMiniAODJetProducer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronCleanedMiniAODJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronCleanedMiniAODJetProducer);
