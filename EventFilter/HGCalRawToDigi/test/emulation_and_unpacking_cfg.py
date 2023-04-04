import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("TEST")

options = VarParsing.VarParsing('standard')
options.register('mode', 'trivial', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,
                 'type of emulation')
options.register('fedId', 0, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,
                 'emulated FED id')
options.register('debug', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,
                 'debugging mode')
options.register('dumpFedRawData', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,
                 'also dump the FEDRawData content')
options.register('numChannelsPerERx', 37, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,
                 'number of channels enabled per ERx')
options.register('numERxsPerECOND', 12, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,
                 'number of ERxs enabled per ECON-D')
options.register('activeECONDs', [], VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.int,
                 'list of ECON-Ds enabled')
options.register('ECONDsInPassthrough', [], VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.int,
                 'list of ECON-Ds in passthrough mode')
options.register('ECONDsInCharacterisation', [], VarParsing.VarParsing.multiplicity.list, VarParsing.VarParsing.varType.int,
                 'list of ECON-Ds in characterisation mode')
options.register('ECONDToTStatus', 3, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,
                 'default ToT status bits (aka TcTp bits) value')
options.register('numCaptureBlocks', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,
                 'number of capture blocks to emulate')
options.register('storeOutput', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,
                 'also store the output into an EDM file')
options.register('storeRAWOutput', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,
                 'also store the RAW output into a streamer file')
options.register('storeEmulatorInfo', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int,
                 'also store the emulator metadata')
options.register('inputFiles',
                 'file:/eos/cms/store/group/dpg_hgcal/tb_hgcal/2023/labtest/module822/pedestal_run0.root',
                 VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string,
                 'input TB file')
options.maxEvents = 1  # number of events to emulate
options.output = 'output.root'  # output EDM file
options.secondaryOutput = 'output.raw'  # output streamer file
options.parseArguments()

process.load('EventFilter.HGCalRawToDigi.hgcalEmulatedSlinkRawData_cfi')
process.load('EventFilter.HGCalRawToDigi.hgcalDigis_cfi')
#process.load('Validation.HGCalValidation.hgcalEmulValidation_cfi')
process.load('EventFilter.HGCalRawToDigi.hgcalEmulatorTest_cfi')  #FIXME

process.load("FWCore.MessageService.MessageLogger_cfi")
if options.debug:
    process.MessageLogger.cerr.threshold = "DEBUG"
    process.MessageLogger.debugModules = ["hgcalEmulatedSlinkRawData", "hgcalDigis"]

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    hgcalEmulatedSlinkRawData = cms.PSet(initialSeed = cms.untracked.uint32(42))
)

process.source = cms.Source("EmptySource")

process.hgcalEmulatedSlinkRawData.emulatorType = options.mode
if process.hgcalEmulatedSlinkRawData.emulatorType == 'hgcmodule':
    process.hgcalEmulatedSlinkRawData.inputs = cms.untracked.vstring(options.inputFiles)

process.hgcalEmulatedSlinkRawData.storeEmulatorInfo = bool(options.storeEmulatorInfo)
process.hgcalEmulatedSlinkRawData.slinkParams.numCaptureBlocks = options.numCaptureBlocks
econd_id = 0
print('S-link: number of capture blocks: {}'.format(
    process.hgcalEmulatedSlinkRawData.slinkParams.numCaptureBlocks.value()))
for econd in process.hgcalEmulatedSlinkRawData.slinkParams.ECONDs:
    # must use 'cms.' python configuration types
    econd.active = cms.bool((econd_id in options.activeECONDs))
    econd.passthroughMode = cms.bool((econd_id in options.ECONDsInPassthrough))
    econd.characterisationMode = cms.bool((econd_id in options.ECONDsInCharacterisation))
    econd.enabledERxs = cms.vuint32([i for i in range(options.numERxsPerECOND)])
    econd.numChannelsPerERx = cms.uint32(options.numChannelsPerERx)
    econd.defaultToTStatus = cms.uint32(options.ECONDToTStatus)
    print('ECON-D {}: active? {}, TcTp: {}, enabled eRxs: {}, number of channels/eRx: {}, passthrough? {}, characterisation? {}'.format(
        econd_id, bool(econd.active), econd.defaultToTStatus.value(),
        [i for i in econd.enabledERxs], econd.numChannelsPerERx.value(),
        bool(econd.passthroughMode), bool(econd.characterisationMode)))
    econd_id += 1

process.hgcalDigis.src = cms.InputTag('hgcalEmulatedSlinkRawData')
process.hgcalDigis.fedIds = cms.vuint32(options.fedId)
process.hgcalDigis.numERxsInECOND = options.numERxsPerECOND
process.hgcalDigis.maxCaptureBlock = options.numCaptureBlocks

process.p = cms.Path(process.hgcalEmulatedSlinkRawData * process.hgcalDigis)

if options.dumpFedRawData:
    process.dump = cms.EDAnalyzer("DumpFEDRawDataProduct",
        label = cms.untracked.InputTag('hgcalEmulatedSlinkRawData'),
        feds = cms.untracked.vint32(options.fedId),
        dumpPayload = cms.untracked.bool(True)
    )
    process.p *= process.dump

process.outpath = cms.EndPath()

if options.storeOutput:
    process.output = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string(options.output),
        outputCommands = cms.untracked.vstring(
            'drop *',
            'keep *_hgcalEmulatedSlinkRawData_*_*',
            'keep *_hgcalDigis_*_*',
        )
    )
    process.outpath += process.output

if options.storeRAWOutput:
    process.outputRAW = cms.OutputModule("FRDOutputModule",
        source = cms.InputTag('hgcalEmulatedSlinkRawData'),
        frdVersion = cms.untracked.uint32(6),
        frdFileVersion = cms.untracked.uint32(1),
        fileName = cms.untracked.string(options.secondaryOutput)
    )
    process.outpath += process.outputRAW
