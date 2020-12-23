#!/bin/sh
cd /storage/user/idutta/CMSSW_9_4_2/src/181120/HHLooper/TaggerFits/condor/condor_output/condor_logs
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`
#cd $TMP
/storage/user/idutta/CMSSW_9_4_2/src/181120/HHLooper/TaggerFits/analyzeTaggerFitsAna $1 $5/$2 $3 $4 $6

