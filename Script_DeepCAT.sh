#!/bin/bash -f
# If raw TCR repertoire sequencing data are available please place data 
# in the "Input" folder
var1="Input"
var2="Output"
# If raw TCR repertoire sequencing data are not available 
# please use our Sample Data as an example 
#python DeepCAT_modif.py $var3  $var4   
var3="iSMART_results"
var4="DeepCAT_CHKP"
#
if [ ! -d $var1 ]; then
  mkdir $var1
fi
if [ "$(ls -A $var1)" ]; then
  if [ ! -d $var3 ]; then
    mkdir $var3
  fi
  python PrepareAdaptiveFile.py $var1 $var2
  python iSMARTm.py -d $var2 -o $var3  
else
    echo "Input directory is Empty,the files in SampleData will be processed"
    var2="SampleData"
    if [ ! -d $var3 ]; then
      mkdir $var3
    fi
    python iSMARTm.py -d $var2 -o $var3
fi  
python DeepCAT.py $var3  $var4
