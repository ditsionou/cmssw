# DPGAnalysis/HGCalTools

## Local production scripts

These scripts run with HTCondor at CERN. The base script can be run as

```
sh runLocalGeneration.sh -c $CMSSW_BASE -m mcm_fragment -n #events -o output_file
```

It will fetch the mcm_fragment and run the cmsDriver commands to produce the GEN-SIM-DIGI-RECO file.
Some extra options are:

* `-s/--simonly` stops the generation at the GEN-SIM tier
* `-p/--pu` the average pileup to generate
* `--minbias` the directory with min-bias files to use
* `-a/--aged` the aged scenario to apply
* `-l/--avg_pileup` average pileup (140, 200, not all values are possible)
* `-p/--pileup_input` directory with min. bias files for the mixing

Geometry, conditions era are hardcoded at the start of the script for the moment. Edit as needed.
The script can be run with condor by using:

```
condor_submit condor.sub
```

A condor submission file with pileup is also available `condor_pu.sub`.
In both cases edit before submitting :)