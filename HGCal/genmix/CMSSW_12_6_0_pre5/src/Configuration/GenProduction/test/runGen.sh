#!/bin/bash

cmssw=/afs/cern.ch/user/p/psilva/work/HGCal/genmix/CMSSW_12_6_0_pre5/src/
work_dir=`pwd`
cd $cmssw
eval `scram r -sh`
cd $work_dir

jobid=$1
cfg=$2
nevts=$3
outdir=/eos/cms/store/cmst3/group/hgcal/CMG_studies/psilva/NANOGEN

cmsDriver.py Configuration/GenProduction/python/${cfg}_cff.py \
    --python_filename  ${cfg}_cfg.py \
    --eventcontent NANOAODGEN \
    --datatier GEN-SIM-DIGI-NANOAOD \
    --fileout file:${cfg}.root \
    --step GEN,SIM,DIGI:pdigi_valid,NANOGEN \
    --mc \
    -n ${nevts} \
    --geometry Extended2026D94 \
    --beamspot HLLHC14TeV  \
    --conditions auto:phase2_realistic_T21 \
    --era Phase2C18I13M9 \
    --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_3000 \
    --customise Configuration/DataProcessing/Utils.addMonitoring \
    --customise PhysicsTools/NanoAOD/nanogen_cff.pruneGenParticlesNano \
    --customise Configuration/GenProduction/nanogencalo_custom_cff.customCalo

mv -v ${cfg}.root ${outdir}/${cfg}_${jobid}.root

