import FWCore.ParameterSet.Config as cms

"""
This parameter set defines the baseline parameters for the HGCROC emulation.
The main parameters are the following
* opMode - operation mode (0=default, 1=characterization)
* adcFSC - full scale charge for the ADC in fC
* adcPulse - pre-amp pulse shape (used for charge leakage)
* totFSC - full scale charge for the time over threshold in fC
* totOnset - the threshold at which the ToT mode is triggered in fC
* totBxUndershoot - number of additional bunches in which the channel is busy due to baseline undershoot (ToT recovery)
* totChargeDrainParam - the parameterization for the number of bunches in which the ToT will keep the channel in busy state
* totChargeDrainJitterParam - the parameterization of the jitter on the charge drain in ToT mode
* toaOnset - threshold in fC for the toa onset
* toaFSC - "full scale charge" for the toa measurement in ns
* toaJitter - jitter for the stochastic component of the time-of-arrival (ns)
* toaClockOffset - clock offset used to generate the constant term (ns)
* pedestal - baseline pedestal in fC
* noiseJitter - the stochastic component of the noise used for smearing
* commonNoise - common mode in fC (in practice will be added to the pedestal but it's useful have as separate component)

"""
hgcrocEmulator = cms.PSet( 
    opMode = cms.uint32(0),
    adcFSC = cms.double(160.),
    adcPulse = cms.vdouble(0.00, 0.017,   0.817,   0.163,  0.003,  0.000),
    totFSC = cms.double(10000.),
    totOnset = cms.double(160.),
    totBxUndershoot = cms.uint32(2),
    totChargeDrainParam = cms.vdouble(0.00077481,-0.001128,0.09327157),
    totChargeDrainJitterParam = cms.vdouble(1.37510708e-05,3.49998849e-03),
    toaOnset = cms.double(12.),
    toaFSC = cms.double(25.),
    toaJitter = cms.double(25.),
    toaClockOffset = cms.double(0.02),
    pulseAvgT = cms.vdouble(0., 0., 1.0, 0.0963/0.9037, 0., 0.),
    tdcResolutionInPs = cms.double( 0.001 ),
    pedestal = cms.double(0.),
    noiseJitter = cms.double(0.4),
    commonNoise = cms.double(0.)
)
