import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.StandardSequences.Eras import eras

options = VarParsing()
options.register ("doseMap", "SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-6.2.0.1.txt",  VarParsing.multiplicity.singleton, VarParsing.varType.string)
options.register ("geometry", "GeometryExtended2026D92Reco",  VarParsing.multiplicity.singleton, VarParsing.varType.string)
options.register ("voltage", 600, VarParsing.multiplicity.singleton, VarParsing.varType.int)
options.register ("annealing", 90, VarParsing.multiplicity.singleton, VarParsing.varType.int)
options.register ("fscale", 1.0, VarParsing.multiplicity.singleton, VarParsing.varType.float)
options.parseArguments()

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
process = cms.Process('demo',Phase2C11I13M9)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.Geometry.{}_cff'.format(options.geometry))
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.source = cms.Source("EmptySource")

from SimCalorimetry.HGCalSimAlgos.hgcSensorOpParams_cfi import hgcSiSensorIleakRadDam,hgcSiSensorCCE 
#from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import HGCAL_ileakParam_toUse, HGCAL_cceParams_toUse


ileakconds=f'{options.voltage}V_{options.annealing}m'
HGCAL_ileakParam_toUse    = cms.PSet(
    ileakParam = cms.double( hgcSiSensorIleakRadDam(ileakconds) )
)

process.analyses=cms.Sequence()
for var in ['nom','up','dn']:
    cceconds=f'{options.voltage}V_{var}_{options.annealing}m'
    HGCAL_cceParams_toUse = cms.PSet(
        cceParamFine  = cms.vdouble(hgcSiSensorCCE(120,cceconds)),
        cceParamThin  = cms.vdouble(hgcSiSensorCCE(200,cceconds)),
        cceParamThick = cms.vdouble(hgcSiSensorCCE(300,cceconds)),
    )

    ana_name=f'conds_{var}'
    setattr(process,
            ana_name,
            cms.EDAnalyzer("HGCalConditionsByAlgoAnalyzer",
                           scaleByDoseFactor  = cms.double(options.fscale),
                           doseMap            = cms.string( options.doseMap ),
                           confAlgo           = cms.vuint32(0,0,0),
                           ileakParam         = HGCAL_ileakParam_toUse,
                           cceParams          = HGCAL_cceParams_toUse,
                           referenceIdark = cms.double(0.5),
                           aimMIPtoADC        = cms.int32(10),
                           ignoreGainSettings = cms.bool(False)
            )
    )
    process.analyses.insert(0,getattr(process,ana_name))
    

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(f"condsbyalgo_output_{options.geometry}_{options.voltage}_{options.annealing}_{options.fscale:3.2f}.root") )

process.p = cms.Path(process.analyses)

