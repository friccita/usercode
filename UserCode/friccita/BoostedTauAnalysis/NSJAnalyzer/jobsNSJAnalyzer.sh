#!/bin/bash

export SCRAM_ARCH=slc5_amd64_gcc462
pwd=$PWD
echo ${pwd}

#------------------------SET VARIABLES---------------------------------
# NFiles is the number of jobs you want to run concurrently
sourceDir=/afs/cern.ch/user/f/friccita/CMSSW_5_2_5/src/BoostedTauAnalysis/NSJAnalyzer
NFiles=100
#inputFiles is the number of files you are running over
inputFilesF=`cmsLs /store/user/friccita/WJetsToMuNu_skim/ | wc -l`
inputFilesR=`cmsLs /store/user/yohay/Summer12_WToMuNu_skim/ | wc -l`
inputFiles=$((${inputFilesF}+${inputFilesR}))
#inputFiles=10
# final_destination is the address to which you want to write the output
final_destination=/afs/cern.ch/user/f/friccita/CMSSW_5_2_5/src/BoostedTauAnalysis/NSJAnalyzer/BkgOutput_11202012
#----------------------------------------------------------------------

#fileList=`cmsLs /store/user/friccita/WJetsToMuNu_skim/ | grep root | awk '{ print $5 }'`
cmsLs /store/user/friccita/WJetsToMuNu_skim/ | grep root | awk '{ print $5 }' > fileList_fran.txt
cmsLs /store/user/yohay/Summer12_WToMuNu_skim/ | grep root | awk '{ print $5 }' > fileList_rachel.txt

cat fileList_fran.txt fileList_rachel.txt > fileList.txt
rm fileList_fran.txt
rm fileList_rachel.txt

p=1
prefix="root://eoscms//eos/cms"

# run below for N files
while [ $p -le $inputFiles ]
  do
  beginning=$p
  if [ $(($p + 4)) -ge $inputFiles ]; then
      end=$inputFiles
  else
      end=$(($p + 4))
  fi
  input_file_list=""
echo $beginning $end

#  if [ ${end} -ge $inputFiles ]; then
#      fileSubList=`head -n${inputFiles} wjets.txt | tail -n${beginning}`
#  else
#      fileSubList=`head -n${end} wjets.txt | tail -n${beginning}`
#  fi

  for linenumber in `seq ${beginning} ${end}`
    do
    
  #get linenumber-th filename
    file=`sed -n "${linenumber} p" fileList.txt `
    echo $file
  #append prefix
    file_name="${prefix}${file}"
    echo $file_name
  #append newline if necessary
    if [ $linenumber -le $end ]; then
	file_name="${file_name}\n"
	input_file_list="${input_file_list}${file_name}"
   else
	file_name="${file_name}"
	input_file_list="${input_file_list}${file_name}"
    fi
  done

echo -e $input_file_list > input_file_list_${p}.txt

cat>${final_destination}/nsjanalyzer_cfg_${p}.py<<EOF
import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("Nsubjettiness")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("BoostedTauAnalysis.NSJAnalyzer.nsjanalyzer_cfi")
process.load("BoostedTauAnalysis.WmunuFilter.WmunuFilter_cff")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.GlobalTag.globaltag = "START52_V9::All"

mylist = FileUtils.loadListFromFile ('${sourceDir}/input_file_list_${p}.txt')
readFiles = cms.untracked.vstring(*mylist)
process.source = cms.Source("PoolSource",
                            fileNames = readFiles
                           )

process.FastJetAnalyzer_alltaus.outFileName = cms.string('FastJetAnalysis_WJets_bkg_pileup_pruned_${p}.root')
process.FastJetAnalyzer_alltaus.outputTextFile = cms.string('JetConstituentInfo_bkg_${p}.txt')
process.FastJetAnalyzer_alltaus.jetSrc = cms.InputTag("ak5PFJets");
process.FastJetAnalyzer_alltaus.Sample = cms.string('WJets')
process.p = cms.Path(process.WmunuFilter*process.FastJetAnalyzer_alltaus)

EOF
  
# ...and you have just built your python file!
  
cat>JobStep_${p}.sh<<EOF
#!/bin/bash

EDMFile="FastJetAnalysis_WJets_bkg_pileup_pruned_${p}.root"
TXTFile="JetConstituentInfo_bkg_${p}.txt"
jobDir="/afs/cern.ch/user/f/friccita/CMSSW_5_2_5/src/BoostedTauAnalysis/NSJAnalyzer"
fileNamePrefix="nsjanalyzer_cfg_${p}"
cfgfilefolder="${final_destination}"

cd \$jobDir
eval \`scramv1 runtime -sh\`
cd -
cp \${cfgfilefolder}/\${fileNamePrefix}.py .
#cmsRun \${fileNamePrefix}.py >& outputlog_${p}.txt
cmsRun \${fileNamePrefix}.py
cmsStage -f \$EDMFile /store/user/friccita/WJetsNSJ
cmsStage -f \$TXTFile /store/user/friccita/WJetsNSJ
rm \$EDMFile \${fileNamePrefix}.py
exit 0

EOF

# ...and you have just built the job to run your python file!

chmod a+x JobStep_${p}.sh
bsub -q 1nh JobStep_${p}.sh
p=$(($p + 5))
  
done

exit 0