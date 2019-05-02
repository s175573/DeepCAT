# DeepCAT
Deep Learning Method to Identify Cancer Associated TCRs

DeepCAT is a computational method based on convolutional neural network to exclusively identify cancer-associated beta chain TCR hypervariable CDR3 sequences. The input data were generated from tumor RNA-seq data and TCR repertoire sequencing data of healthy donors. Users do not need to perform training or evaluation. Instead, users can directly apply the PredictCancer function in the package, after downloading the CHKP folder. 
Standard pipeline of using DeepCAT:

Input data is beta CDR3 amino acid sequence filtered out the non-productive and non-amino-acid letters, with the first column being the CDR3 sequences, and optional other columns in the file. It is recommended to document variale gene symbols and clonal frequencies for each CDR3 sequence in the second and third columns. Please see the data/ folder for an example input.

If starting from raw TCR repertoire sequencing data produced from AdaptiveBiotech immuneAnalyzer, one can use the `PrepareAdaptiveInput.R` function to pre-process the data. Then apply `iSMARTm.py` to perform clustering:

In R console:
```
source('PrepareAdaptiveInput.R')
PrepareAdaptiveInput(InputDataFolder, OutputDataFolder)
```
In Mac or Linux terminal:
```
python iSMARTm.py -d OutputDataFolder -o iSMART_results
```

After running iSMART, use DeepCAT classifiers to perform cancer score prediction:

```
from DeepCAT import *
dir_prefix='YOUR_DIR_To_CHKP_Folder/'
ffs=os.listdir('iSMART_results/')
ff=ffs[0]
ff='YOUR_INPUT_FILE_OF_CDR3s.txt'
_, XX,_=PredictCancer(DIR+ff, dir_prefix=dir_prefix)
np.mean(XX)  ## cancer score
```
