import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import CandVars,Var


def customCalo(process):

    process.MessageLogger.cerr.FwkReport.reportEvery=100

    if not process.nanogenSequence.contains(process.finalGenParticles):
        raise ValueError("Can't find finalGenParticles in in nanogenSequence")

        
    process.finalGenParticles.select=cms.vstring('keep status==1')
    
    #add calo particles and sim clusters
    process.load('DPGAnalysis.CaloNanoAOD.caloParticles_cff')
    process.caloParticleTable.cut = cms.string("abs(momentum().eta) > 1.5 || abs(momentum().eta) < 3")
    process.caloParticleTables = cms.Sequence(process.caloParticleTable)
    process.nanogenSequence.insert(0, process.caloParticleTables)

    process.load('DPGAnalysis.CaloNanoAOD.simClusters_cff')
    process.simClusterTable.cut = cms.string("abs(momentum().eta) > 1.5 || abs(momentum().eta) < 3")
    process.simClusterTables = cms.Sequence(process.simClusterTable)
    process.nanogenSequence.insert(0, process.simClusterTables)

    from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
    randSvc=RandomNumberServiceHelper(process.RandomNumberGeneratorService)
    randSvc.populate()


    return process
