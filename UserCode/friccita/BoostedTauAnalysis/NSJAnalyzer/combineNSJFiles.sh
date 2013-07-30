#!/bin/bash

export SCRAM_ARCH=slc5_amd64_gcc462
pwd=$PWD
echo ${pwd}
eval `scramv1 runtime -sh`
cmsLs /store/user/friccita/WJetsNSJ/ | grep root | awk '{ print $5 }' > NSJfileList.txt

copy_destination=/data1/friccita/NSJdatasets_10242012/WJetsTemp

#inputFiles is the number of files you are running over
inputFiles=`cmsLs /store/user/friccita/WJetsNSJ/ | grep root | wc -l`
p=1
endp=(${inputFiles} + 1)
startp=1
prefix="root://eoscms//eos/cms"

for linenumber in `seq ${startp} ${endp}`
  do
  file=""
  file=`sed -n "${linenumber} p" NSJfileList.txt `
  if [ $linenumber -le $endp ]; then  
      echo $file
      cmsStage -f $file ${copy_destination}
  fi

done

exit 0
