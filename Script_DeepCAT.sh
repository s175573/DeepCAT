#!/bin/bash -f
# If raw TCR repertoire sequencing data are available please place data 
# in the "Input" folder
var1="Output"
# If raw TCR repertoire sequencing data are not available 
# please use our Sample Data as an example 
#python DeepCAT_modif.py $var3  $var4   
var2="iSMART_results"
var3="DeepCAT_CHKP"
#
args=("$@")
if [ ${args[0]} == '-r' ]; then
  if [ "$(ls -A ${args[1]})" ]; then
     if [ ! -d $var2 ]; then
        mkdir $var2
     fi
     python PrepareAdaptiveFile.py ${args[1]} $var1
     python iSMARTm.py -d $var1 -o $var2  
  else
     echo "Error! The" ${args[1]} "directory is empty"
     exit 1
  fi   
elif [ ${args[0]} == '-t' ]; then
    echo "The files in ${args[1]} will be processed"
    if [ ! -d $var2 ]; then
      mkdir $var2
    fi
    python iSMARTm.py -d ${args[1]} -o $var2
fi  
python DeepCAT.py $var2  $var3
