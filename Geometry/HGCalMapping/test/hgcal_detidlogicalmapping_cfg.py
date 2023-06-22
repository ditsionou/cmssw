import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.StandardSequences.Eras import eras
from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

options = VarParsing()
options.register ("geom", "v17",  VarParsing.multiplicity.singleton, VarParsing.varType.string)
options.parseArguments()

process = cms.Process("demo",eras.Phase2C17I13M9)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if options.geom == 'v16':
    process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
elif options.geom == 'v17':
    process.load('Configuration.Geometry.GeometryExtended2026D92Reco_cff')
else:
    raise Exception('UNKNOWN GEOMETRY!')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.source = cms.Source("EmptySource")
process.MessageLogger = cms.Service("MessageLogger",
    cerr = cms.untracked.PSet(
        enable = cms.untracked.bool(False)
    ),
    cout = cms.untracked.PSet(
        enable = cms.untracked.bool(True),
        threshold = cms.untracked.string('INFO')
    )
)

process.analyzer = cms.EDAnalyzer("HGCalDetIdLogicalMappingTester",
                                  modFile = cms.string('Geometry/HGCalMapping/data/modulelocator.txt'),
                                  siFile = cms.string('Geometry/HGCalMapping/data/WaferCellMapTraces.txt'),
                                  sipmFile = cms.string('Geometry/HGCalMapping/data/channels_sipmontile.hgcal.txt'),
                                  label = cms.string(''))

process.p = cms.Path(
    process.analyzer
)
