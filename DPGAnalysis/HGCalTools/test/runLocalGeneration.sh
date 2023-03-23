#!/bin/bash

nevts=1000
simonly=0
SHORT=c:,m:,n:,s:,o:,h
LONG=cmssw:,mcm:,nevts:,simonly:,output:,help
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
    -h | --help)
      echo ""
      echo "runLocalGeneration.sh -c cmssw_dir -m mcm_fragment -n nevts -s simonlyflag -o output"
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
cd ${cmssw}/src
eval `scram r -sh`
# Download fragment from McM
cfg=Configuration/GenProduction/python/${mcmfragment}-fragment.py
curl -s -k https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_fragment/${mcmfragment} --retry 3 --create-dirs -o ${cfg}
[ -s ${cfg} ] || exit $?;
ls ${cfg}
scram b
cd ../..
cd ${work}

#run generation
cmsDriver.py ${cfg} -s GEN,SIM -n ${nevts} --conditions auto:phase2_realistic_T25 \
             --geometry Extended2026D99 --era Phase2C17I13M9 --eventcontent FEVTDEBUG \
             --beamspot HLLHC --datatier GEN-SIM --fileout file:step1.root \
             --customise_command "from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper; randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService); randSvc.populate();"

if [ $simonly == "1" ]; then
    echo "This is a simonly production - moving to output"
    mv -v step1.root ${output}
    exit -1
fi

