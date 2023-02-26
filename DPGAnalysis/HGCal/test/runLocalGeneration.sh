#!/bin/bash

work_dir=`pwd`
echo "Working directory is: $work_dir"

#source environment
cmssw_release=$1
cd ${cmssw_release}/src
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
echo "CMSSW_BASE=$CMSSW_BASE"
cd ${work_dir}

#configure script
n=$2
cfg=$3
outf=$4
scenario=$5
digi_custom=""
if [ $scenario == "aged_3iab" ]; then
    digi_custom="--customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_3000"
    echo "Will apply aging for 3/ab"
fi

outdir=`dirname ${outf}`
echo "Will generate $n events from cfg=${cfg} with output @ ${outf}"

#change if needed for other conditions
common="--conditions auto:phase2_realistic_T25 --geometry Extended2026D99 --era Phase2C17I13M9 --eventcontent FEVTDEBUG"
custom="from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper; randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService); randSvc.populate();"

cmsDriver.py ${cfg} -s GEN,SIM -n ${n} ${common} --beamspot HLLHC --datatier GEN-SIM --customise_command "${custom}" --fileout file:step1.root

cmsDriver.py step2  -s DIGI:pdigi_valid,L1TrackTrigger,L1,DIGI2RAW,HLT:@fake2 ${common} ${digi_custom} --datatier GEN-SIM-DIGI-RAW -n -1 --filein file:step1.root  --fileout file:step2.root 

cmsDriver.py step3  -s RAW2DIGI,RECO,RECOSIM ${common} --datatier GEN-SIM-RECO -n -1 --filein  file:step2.root  --fileout file:step3.root

mkdir -p ${outdir}
cp -v step3.root ${outf}
rm step*.root
echo "All done!"
