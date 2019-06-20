import FWCore.ParameterSet.Config as cms

maxEvents = 10000
doLower = True
doMuTau = True
doETau = True
doTauTau = False
nCore = 8

from Configuration.StandardSequences.Eras import eras
era = eras.Run2_2018
process = cms.Process("BoostedTauReco", era)
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source(
    "PoolSource", 
    fileNames=readFiles, 
    secondaryFileNames=secFiles,
)

process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(maxEvents)
)

readFiles.extend([
    #'file:patMiniAOD_standard.root',
    #'/store/relval/CMSSW_10_5_0_pre1/RelValZTT_13/MINIAODSIM/PU25ns_103X_upgrade2018_realistic_v8-v1/20000/EA29017F-9967-3F41-BB8A-22C44A454235.root',
    '/store/mc/RunIIFall17MiniAODv2/SUSYGluGluToHToAA_AToMuMu_AToTauTau_M-125_M-5_TuneCUETP8M1_13TeV_madgraph_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/B6992772-6143-E811-9E47-002481D24972.root',
])

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')


import Configuration.EventContent.EventContent_cff as evtContent
process.output = cms.OutputModule(
    'PoolOutputModule',
    fileName=cms.untracked.string('miniAOD_BoostedDiTauReco.root'),
    fastCloning=cms.untracked.bool(False),
    dataset=cms.untracked.PSet(
        dataTier=cms.untracked.string('MINIAODSIM'),
        filterName=cms.untracked.string('')
    ),
    outputCommands = evtContent.MINIAODSIMEventContent.outputCommands,
    SelectEvents=cms.untracked.PSet(
        SelectEvents=cms.vstring('*',)
        )
    )

process.out = cms.EndPath(process.output)

######################
### rerun tau reco ###
######################
#import RecoTauTag.Configuration.tools.adaptToRunAtMiniAOD as tauAtMiniTools
import BoostedDiTau.DiTauProducers.adaptToRunAtMiniAOD as tauAtMiniTools
tauAtMiniTools.addTauReReco(process)

# lower the pt threshold
def lowerTauPt(process,postfix='',tauPt=8, jetPt=5):
    from FWCore.ParameterSet.MassReplace import massSearchReplaceParam
    massSearchReplaceParam(getattr(process,'miniAODTausTask'+postfix),'minJetPt',14,jetPt)
    getattr(process,'selectedPatTaus'+postfix).cut = cms.string("pt > {} && tauID(\'decayModeFindingNewDMs\')> 0.5".format(tauPt))

###############################
### lower pt threshold taus ###
###############################

if doLower:
    # TODO: use the Jet reclustering tools rather than by hand
    pfCandSrc = 'packedPFCandidates'
    jetSrc = 'patJetsLowerPt'
    recoJetSrc = 'ak4PFJetsLowerPt'
    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
    ak4PFJetsNew = ak4PFJets.clone(
        src=cms.InputTag(pfCandSrc)
    )
    setattr(process,recoJetSrc,ak4PFJetsNew)
    from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import _patJets
    patJetsNew = _patJets.clone(
        jetSource            = cms.InputTag(recoJetSrc),
        addJetCorrFactors    = cms.bool(False),
        jetCorrFactorsSource = cms.VInputTag(),
        addBTagInfo          = cms.bool(False),
        addDiscriminators    = cms.bool(False),
        discriminatorSources = cms.VInputTag(),
        addAssociatedTracks  = cms.bool(False),
        addJetCharge         = cms.bool(False),
        addGenPartonMatch    = cms.bool(False),
        embedGenPartonMatch  = cms.bool(False),
        addGenJetMatch       = cms.bool(False),
        getJetMCFlavour      = cms.bool(False),
        addJetFlavourInfo    = cms.bool(False),
    )
    setattr(process,jetSrc,patJetsNew)
    
    process.lowerPtTausTask = cms.Task(
        getattr(process,recoJetSrc),
        getattr(process,jetSrc),
    )
    
    tauAtMiniTools.adaptTauToMiniAODReReco(process, jetSrc=jetSrc, postfix='LowerPt')
    process.miniAODTausTaskLowerPt.add(process.lowerPtTausTask)
    lowerTauPt(process,postfix='LowerPt')
    process.output.outputCommands.append('keep *_selectedPatTausLowerPt_*_*')

#########################
### muon cleaned taus ###
#########################

if doMuTau:
    #jetSrc = 'patJetsMuonCleaned'
    #recoJetSrc = 'ak4PFJetsMuonCleaned'
    #process.recoMuonsForJetCleaning = cms.EDFilter('PATMuonRefSelector',
    #    src = cms.InputTag('slimmedMuons'),
    #    cut = cms.string('pt > 3.0 && isPFMuon && (isGlobalMuon || isTrackerMuon)'),
    #)
    #
    #from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
    #process.ak4PFJetsForMuonCleaning = ak4PFJets.clone(
    #    src=cms.InputTag('packedPFCandidates')
    #)
    #
    #process.ak4PFJetsMuonCleaned = cms.EDProducer(
    #    'MuonCleanedMiniAODJetProducer',
    #    jetSrc = cms.InputTag("ak4PFJetsForMuonCleaning"),
    #    muonSrc = cms.InputTag("recoMuonsForJetCleaning"),
    #    pfCandSrc = cms.InputTag("packedPFCandidates"),
    #)
    #from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import _patJets
    #patJetsNew = _patJets.clone(
    #    jetSource            = cms.InputTag(recoJetSrc),
    #    addJetCorrFactors    = cms.bool(False),
    #    jetCorrFactorsSource = cms.VInputTag(),
    #    addBTagInfo          = cms.bool(False),
    #    addDiscriminators    = cms.bool(False),
    #    discriminatorSources = cms.VInputTag(),
    #    addAssociatedTracks  = cms.bool(False),
    #    addJetCharge         = cms.bool(False),
    #    addGenPartonMatch    = cms.bool(False),
    #    embedGenPartonMatch  = cms.bool(False),
    #    addGenJetMatch       = cms.bool(False),
    #    getJetMCFlavour      = cms.bool(False),
    #    addJetFlavourInfo    = cms.bool(False),
    #)
    #setattr(process,jetSrc,patJetsNew)
    #
    #
    #process.muonCleanedHPSPFTausTask = cms.Task(
    #    process.recoMuonsForJetCleaning,
    #    process.ak4PFJetsForMuonCleaning,
    #    process.ak4PFJetsMuonCleaned,
    #    process.patJetsMuonCleaned,
    #)
    #
    #pfCandSrc = cms.InputTag(recoJetSrc,'particleFlowMuonCleaned')
    ##pfCandSrc = "packedPFCandidates"
    #
    #tauAtMiniTools.adaptTauToMiniAODReReco(process, jetSrc=jetSrc, pfCandSrc=pfCandSrc, postfix='MuonCleaned')
    #process.miniAODTausTaskMuonCleaned.add(process.muonCleanedHPSPFTausTask)
    #lowerTauPt(process,postfix='MuonCleaned')
    #process.output.outputCommands.append('keep *_selectedPatTausMuonCleaned_*_*')

    # alternative approach, full recluster jets without muons
    jetSrc = 'patJetsMuonCleaned'
    recoJetSrc = 'ak4PFJetsMuonCleaned'
    pfCandSrc = "muonCleanedPackedPFCandidates"
    
    process.recoMuonsForJetCleaning = cms.EDFilter('PATMuonRefSelector',
        src = cms.InputTag('slimmedMuons'),
        cut = cms.string('pt > 3.0 && isPFMuon && (isGlobalMuon || isTrackerMuon)'),
    )

    process.muonCleanedPackedPFCandidates = cms.EDProducer("MuonCleanedPackedCandidateProducer",
        src = cms.InputTag("recoMuonsForJetCleaning"),
        pfCandSrc = cms.InputTag("packedPFCandidates"),
    )
    
    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
    process.ak4PFJetsMuonCleaned = ak4PFJets.clone(
        src=cms.InputTag('muonCleanedPackedPFCandidates')
    )
    
    from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import _patJets
    patJetsNew = _patJets.clone(
        jetSource            = cms.InputTag(recoJetSrc),
        addJetCorrFactors    = cms.bool(False),
        jetCorrFactorsSource = cms.VInputTag(),
        addBTagInfo          = cms.bool(False),
        addDiscriminators    = cms.bool(False),
        discriminatorSources = cms.VInputTag(),
        addAssociatedTracks  = cms.bool(False),
        addJetCharge         = cms.bool(False),
        addGenPartonMatch    = cms.bool(False),
        embedGenPartonMatch  = cms.bool(False),
        addGenJetMatch       = cms.bool(False),
        getJetMCFlavour      = cms.bool(False),
        addJetFlavourInfo    = cms.bool(False),
    )
    setattr(process,jetSrc,patJetsNew)
    
    
    process.muonCleanedHPSPFTausTask = cms.Task(
        process.recoMuonsForJetCleaning,
        process.muonCleanedPackedPFCandidates,
        process.ak4PFJetsMuonCleaned,
        process.patJetsMuonCleaned,
    )
    
    tauAtMiniTools.adaptTauToMiniAODReReco(process, jetSrc=jetSrc, pfCandSrc=pfCandSrc, postfix='MuonCleaned')
    process.miniAODTausTaskMuonCleaned.add(process.muonCleanedHPSPFTausTask)
    lowerTauPt(process,postfix='MuonCleaned')
    process.output.outputCommands.append('keep *_selectedPatTausMuonCleaned_*_*')


#############################
### electron cleaned taus ###
#############################

if doETau:
    from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
    # 2018
    #setupEgammaPostRecoSeq(process,
    #                       runEnergyCorrections=False, #as energy corrections are not yet availible for 2018
    #                       era='2018-Prompt')
    # 2017
    setupEgammaPostRecoSeq(process,
                           runVID=False, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                           era='2017-Nov17ReReco')  
    # 2016
    #setupEgammaPostRecoSeq(process,
    #                       runEnergyCorrections=False, #corrections by default are fine so no need to re-run
    #                       era='2016-Legacy')

    # alternative approach, full recluster jets without electrons
    jetSrc = 'patJetsElectronCleaned'
    recoJetSrc = 'ak4PFJetsElectronCleaned'
    pfCandSrc = "electronCleanedPackedPFCandidates"
    
    process.recoElectronsForJetCleaning = cms.EDFilter('PATElectronRefSelector',
        src = cms.InputTag('slimmedElectrons'),
        cut = cms.string('electronID("mvaEleID-Fall17-noIso-V1-wp90")'),
    )

    process.electronCleanedPackedPFCandidates = cms.EDProducer("ElectronCleanedPackedCandidateProducer",
        src = cms.InputTag("recoElectronsForJetCleaning"),
        pfCandSrc = cms.InputTag("packedPFCandidates"),
    )
    
    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
    process.ak4PFJetsElectronCleaned = ak4PFJets.clone(
        src=cms.InputTag('electronCleanedPackedPFCandidates')
    )
    
    from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import _patJets
    patJetsNew = _patJets.clone(
        jetSource            = cms.InputTag(recoJetSrc),
        addJetCorrFactors    = cms.bool(False),
        jetCorrFactorsSource = cms.VInputTag(),
        addBTagInfo          = cms.bool(False),
        addDiscriminators    = cms.bool(False),
        discriminatorSources = cms.VInputTag(),
        addAssociatedTracks  = cms.bool(False),
        addJetCharge         = cms.bool(False),
        addGenPartonMatch    = cms.bool(False),
        embedGenPartonMatch  = cms.bool(False),
        addGenJetMatch       = cms.bool(False),
        getJetMCFlavour      = cms.bool(False),
        addJetFlavourInfo    = cms.bool(False),
    )
    setattr(process,jetSrc,patJetsNew)
    
    
    process.electronCleanedHPSPFTausTask = cms.Task(
        process.recoElectronsForJetCleaning,
        process.electronCleanedPackedPFCandidates,
        process.ak4PFJetsElectronCleaned,
        process.patJetsElectronCleaned,
    )
    
    
    tauAtMiniTools.adaptTauToMiniAODReReco(process, jetSrc=jetSrc, pfCandSrc=pfCandSrc, postfix='ElectronCleaned')
    process.miniAODTausTaskElectronCleaned.add(process.electronCleanedHPSPFTausTask)
    lowerTauPt(process,postfix='ElectronCleaned')
    process.output.outputCommands.append('keep *_selectedPatTausElectronCleaned_*_*')

    process.TauRecoElectronCleaned += process.egammaPostRecoSeq

####################
### boosted taus ###
####################

if doTauTau:
    pfCandSrc = 'packedPFCandidates'
    jetSrc = 'patJetsBoosted'
    recoJetSrc = 'boostedTauSeeds'
    
    import CommonTools.ParticleFlow.pfNoPileUp_cff as boostedTaus
    process.pfPileUpForBoostedTaus = boostedTaus.pfPileUp.clone(
        PFCandidates = cms.InputTag('packedPFCandidates'),
        checkClosestZVertex = cms.bool(False)
    )
    process.pfNoPileUpForBoostedTaus = boostedTaus.pfNoPileUp.clone(
        topCollection = cms.InputTag('pfPileUpForBoostedTaus'),
        bottomCollection = cms.InputTag('packedPFCandidates')
    )
    
    import RecoJets.JetProducers.ak4PFJets_cfi as boostedTaus2
    import RecoJets.JetProducers.CMSBoostedTauSeedingParameters_cfi as boostedTaus3
    process.ca8PFJetsCHSprunedForBoostedTaus = boostedTaus2.ak4PFJets.clone(
        boostedTaus3.CMSBoostedTauSeedingParameters,
        #src = cms.InputTag('pfNoPileUpForBoostedTaus'),
        jetPtMin = cms.double(10.0),
        doAreaFastjet = cms.bool(True),
        nFilt = cms.int32(100),
        rParam = cms.double(0.8),
        jetAlgorithm = cms.string("CambridgeAachen"),
        writeCompound = cms.bool(True),
        jetCollInstanceName = cms.string('subJetsForSeedingBoostedTaus')
    )
    
    process.boostedTauSeeds = cms.EDProducer("BoostedTauSeedsProducer",
        subjetSrc = cms.InputTag('ca8PFJetsCHSprunedForBoostedTaus', 'subJetsForSeedingBoostedTaus'),
        pfCandidateSrc = cms.InputTag('packedPFCandidates'),
        verbosity = cms.int32(0)
    )
    
    from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import _patJets
    patJetsNew = _patJets.clone(
        jetSource            = cms.InputTag(recoJetSrc),
        addJetCorrFactors    = cms.bool(False),
        jetCorrFactorsSource = cms.VInputTag(),
        addBTagInfo          = cms.bool(False),
        addDiscriminators    = cms.bool(False),
        discriminatorSources = cms.VInputTag(),
        addAssociatedTracks  = cms.bool(False),
        addJetCharge         = cms.bool(False),
        addGenPartonMatch    = cms.bool(False),
        embedGenPartonMatch  = cms.bool(False),
        addGenJetMatch       = cms.bool(False),
        getJetMCFlavour      = cms.bool(False),
        addJetFlavourInfo    = cms.bool(False),
    )
    setattr(process,jetSrc,patJetsNew)
    
    process.boostedHPSPFTausTask = cms.Task(
        process.pfPileUpForBoostedTaus,
        process.pfNoPileUpForBoostedTaus,
        process.ca8PFJetsCHSprunedForBoostedTaus,
        process.boostedTauSeeds,
        getattr(process,jetSrc),
    )
    
    tauAtMiniTools.adaptTauToMiniAODReReco(process, jetSrc=jetSrc, postfix='Boosted')
    process.miniAODTausTaskBoosted.add(process.boostedHPSPFTausTask)
    lowerTauPt(process,postfix='Boosted')
    process.recoTauAK4Jets08RegionBoosted.pfCandAssocMapSrc = cms.InputTag('boostedTauSeeds', 'pfCandAssocMapForIsolation')
    process.ak4PFJetsRecoTauChargedHadronsBoosted.builders[1].dRcone = cms.double(0.3)
    process.ak4PFJetsRecoTauChargedHadronsBoosted.builders[1].dRconeLimitedToJetArea = cms.bool(True)
    process.hpsPFTauDiscriminationByLooseMuonRejectionSimpleBoosted.dRmuonMatch = 0.1
    process.hpsPFTauDiscriminationByTightMuonRejectionSimpleBoosted.dRmuonMatch = 0.1
    process.output.outputCommands.append('keep *_selectedPatTausBoosted_*_*')


###############
### options ###
###############
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
)
process.options.numberOfThreads = cms.untracked.uint32(nCore)
process.options.numberOfStreams = cms.untracked.uint32(0)

process.options = cms.untracked.PSet(
    process.options,
    wantSummary=cms.untracked.bool(False)
)
