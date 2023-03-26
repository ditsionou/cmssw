#!/bin/bash

GEOMETRY=Extended2026D99
ERA=Phase2C17I13M9
CONDITIONS=auto:phase2_realistic_T25
nevts=1000
simonly=0
avg_pileup=0

#parse command line options
SHORT=c:,m:,n:,s:,o:,a:,p:,l:,h
LONG=cmssw:,mcm:,nevts:,simonly:,output:,aged:,pileup_input:,avg_pileup:,help
OPTS=$(getopt -a -n weather --options $SHORT --longoptions $LONG -- "$@")
eval set -- "$OPTS"
while :
do
  case "$1" in
    -c | --cmssw )
      cmssw="$2"
      shift 2
      ;;
    -m | --mcm )
      mcmfragment="$2"
      shift 2
      ;;
    -n | --nevts )
      nevts="$2"
      shift 2
      ;;
    -s | --simonly )
      simonly="$2"
      shift 2
      ;;
    -o | --output )
      output="$2"
      shift 2
      ;;
    -a | --aged )
      aged="$2"
      shift 2
      ;;
    -p | --pileup_input )
      pileup_input="$2"
      shift 2
      ;;
    -l | --avg_pileup )
      avg_pileup="$2"
      shift 2
      ;;
    -h | --help)
      echo ""
      echo "runLocalGeneration.sh -c cmssw_dir -m mcm_fragment -n nevts -s simonlyflag -o output -a aged"
      echo "Some examples are the following"
      echo "  mcm_fragment = min.bias EGM-Run3Summer19GS-00020 "
      echo "                 ttbar L1T-PhaseIITDRSpring19GS-00005"
      echo "  simonlyflag = 0/1, if 1 will stop after SIM and copy it to the output"
      echo "  aged - if not given it will do vanilla"
      echo "       - startup SimCalorimetry/HGCalSimProducers/hgcalDigitizer_cfi.HGCal_setRealisticStartupNoise"
      echo "       - 3/ab SimCalorimetry/HGCalSimProducers/hgcalDigitizer_cfi.HGCal_setEndOfLifeNoise"
      echo ""
      exit 2
      ;;
    --)
      shift;
      break
      ;;
    *)
      echo "Unexpected option: $1"
      ;;
  esac
done

#prepare the output
outdir=`dirname ${output}`
mkdir -p ${outdir}
echo "Output will be available in ${outdir}"

#setup CMSSW
work=`pwd`
cd ${cmssw}
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
# Download fragment from McM
cfg=Configuration/GenProduction/python/${mcmfragment}-fragment.py
curl -s -k https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_fragment/${mcmfragment} --retry 3 --create-dirs -o ${cfg}
[ -s ${cfg} ] || exit $?;
ls ${cfg}
scram b
cd ../..
cd ${work}

#run generation
cmsDriver.py ${cfg} -s GEN,SIM -n ${nevts} \
             --conditions ${CONDITIONS} --geometry ${GEOMETRY} --era ${ERA} \
             --eventcontent FEVTDEBUG --beamspot HLLHC --datatier GEN-SIM --fileout file:step1.root \
             --customise_command "from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper; randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService); randSvc.populate();"

if [ $simonly == "1" ]; then
    echo "This is a simonly production - moving to output"
    mv -v step1.root ${output}
    exit -1
fi

#run digitization
if [ -z $aged ]; then
    echo "No customisation will be applied - vanilla digis will run"
else
    aging_customise="--customise ${aged}"
    echo "Aging will be customised with ${aging_customise}"
fi
if [ -z $avg_pileup ] || [ -z $pileup_input ]; then
    echo "No pileup scenario"
else
    echo "Will generate average pileup of ${avg_pileup} with files from ${pileup_input}"
    pileup_input=`find ${pileup_input} -iname "*.root" -printf "file:%h/%f,"`
    pileup_input=${pileup_input::-1}
    pileup_costumise="--pileup AVE_${avg_pileup}_BX_25ns --pileup_input ${pileup_input}"
fi
cmsDriver.py step2 -s DIGI:pdigi_valid,L1TrackTrigger,L1,DIGI2RAW,HLT:@fake2 \
             --conditions ${CONDITIONS} --geometry ${GEOMETRY} --era ${ERA} \
             --eventcontent FEVTDEBUG --datatier GEN-SIM-DIGI-RAW -n -1 --filein file:step1.root --fileout file:step2.root \
             ${aging_customise} ${pileup_costumise}

#run reconstruction step
cmsDriver.py step3 -s RAW2DIGI,RECO,RECOSIM \
             --conditions ${CONDITIONS} --geometry ${GEOMETRY} --era ${ERA} \
             --eventcontent FEVTDEBUG --datatier GEN-SIM-RECO -n -1 --filein file:step2.root --fileout file:step3.root

#move output and cleanup
mv -v step3.root ${output}
rm *.root
