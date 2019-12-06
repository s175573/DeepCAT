#! usr/bin/python

import os
from os.path import exists
import numpy as np
import csv
from csv import reader
import sys

indir=sys.argv[1]
outdir=sys.argv[2]
thr=10000

#def PrepareAdaptiveFile(indir,outdir,thr=10000):
ffs=os.listdir(indir)
for ff in ffs:
    if('.tsv' not in ff):  
       continue
    ff0=ff
    if not os.path.exists(outdir):
       os.makedirs(outdir)   
    str1='TestReal-'   
    newff=outdir+'/'+str1+ff0
#    if exists(newff)==False:
#       continue
    csv_reader = reader(open(indir+'/'+ff,"r"), delimiter='\t', quotechar="\"") 
    ddnew=[] 
    for row in csv_reader:
      if '*' not in row[1]:
        if 'X' not in row[1]:
#         if '^C.+F$' not in row[1]:
          if (len(row[1])>=10) and (len(row[1])<=24):
            if 'unresolved' not in row[5]:
              if (row[1][0]=='C') and (row[1][-1]=='F'):
                  ddnew.append(row)
    ddnew=np.array(ddnew)      
    sorted_array = ddnew[ddnew[:,3].argsort()]   
    reverse_array = sorted_array[::-1]
    if len(reverse_array)>thr:
       col1=reverse_array[0:thr,1]
       col2=reverse_array[0:thr,5]
       col3=reverse_array[0:thr,3]
    else:
       col1=reverse_array[:,1]
       col2=reverse_array[:,5]
       col3=reverse_array[:,3]
    c=zip(col1,col2,col3)         
    with open(newff, 'w') as f:
       writer = csv.writer(f, delimiter='\t')
       writer.writerows(c)
