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
* toaOnset - threshold in fC for the toa onset
* toaFSC - "full scale charge" for the toa measurement in ns
* toaJitter - jitter for the stochastic component of the time-of-arrival (ns)
* toaClockOffset - clock offset used to generate the constant term (ns)
"""
hgcrocEmulator = cms.PSet( 
    opMode = cms.uint32(0),
    adcFSC = cms.double(160.),
    adcPulse = cms.vdouble(0.00, 0.017,   0.817,   0.163,  0.003,  0.000),
    totFSC = cms.double(10000.),
    totOnset = cms.double(160.),
    totBxUndershoot = cms.uint32(2),
    totChargeDrainParam = cms.vdouble(-919.13, 365.36, -14.10, 0.2,
                                      -21.85, 49.39,  22.21,  0.8,
                                      -0.28,   27.14,  43.95,
                                      3.89048 ),
    toaOnset = cms.double(12.),
    toaFSC = cms.double(25.),
    toaJitter = cms.double(25.),
    toaClockOffset = cms.double(0.02),

    pulseAvgT         = cms.vdouble(0.00, 23.42298,13.16733,6.41062,5.03946,4.5320),
    tdcResolutionInPs = cms.double( 0.001 )
)
