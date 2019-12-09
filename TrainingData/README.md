The training data contains two files, each is a list of CDR3s coming from either cancer or healthy individuals. 

To train DeepCAT from scratch, please download all the files, put them in the same directory as the DeepCAT source code, and run the following command in the python console:

```python
from DeepCAT import *
batchTrain(ftumor='TumorCDR3.txt',n=10, STEPs=20000, rate=0.2, fnormal='NormalCDR3.txt')
```

This function performs n (=10 here) times cross-validation by subsampling 1-rate (80%) of the data for training, and the remaining 20% for validation. 

It will create a subdirectory under the current path, `/tmp/`, which stores all the checkpoint folders and ROC curves for each run. 
