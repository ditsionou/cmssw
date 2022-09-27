import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.StandardSequences.Eras import eras

options = VarParsing()
options.register ("doseMap",    "SimCalorimetry/HGCalSimProducers/data/doseParams_3000fb_fluka-6.2.0.1.txt",  VarParsing.multiplicity.singleton, VarParsing.varType.string)
options.register ("geometry",   "GeometryExtended2026D92Reco",  VarParsing.multiplicity.singleton, VarParsing.varType.string)
options.register ("conditions", "CERN21_600V_120m", VarParsing.multiplicity.singleton, VarParsing.varType.string)
options.parseArguments()

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
process = cms.Process('demo',Phase2C11I13M9)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.Geometry.{}_cff'.format(options.geometry))
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.source = cms.Source("EmptySource")

from SimCalorimetry.HGCalSimAlgos.hgcSensorOpParams_cfi import hgcSiSensorIleak,hgcSiSensorCCE 
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import HGCAL_ileakParam_toUse, HGCAL_cceParams_toUse

HGCAL_ileakParam_toUse    = cms.PSet(
    ileakParam = cms.vdouble( hgcSiSensorIleak(options.conditions) )
)

HGCAL_cceParams_toUse = cms.PSet(
    cceParamFine  = cms.vdouble(hgcSiSensorCCE(120,options.conditions)),
    cceParamThin  = cms.vdouble(hgcSiSensorCCE(200,options.conditions)),
    cceParamThick = cms.vdouble(hgcSiSensorCCE(300,options.conditions)),
)

process.conds_3iab = cms.EDAnalyzer("HGCalConditionsByAlgoAnalyzer",
                                    scaleByDoseFactor  = cms.double(1.),
                                    doseMap            = cms.string( options.doseMap ),
                                    confAlgo           = cms.vuint32(0,0,0),
                                    ileakParam         = HGCAL_ileakParam_toUse,
                                    cceParams          = HGCAL_cceParams_toUse,
                                    referenceIdark = cms.double(0.5),
                                    aimMIPtoADC        = cms.int32(10),
                                    ignoreGainSettings = cms.bool(False)
)

#process.plotter_eol_nogain = process.plotter_eol.clone( ignoreGainSettings = cms.bool(True) )

#process.plotter_start = process.plotter_eol.clone( doseMapAlgo=cms.uint32(3) )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("condsbyalgo_output_{}_{}.root".format(options.geometry,options.conditions))
                               )

process.p = cms.Path(process.conds_3iab
                     #*process.plotter_eol_nogain
                     #*process.plotter_start
)
