#! usr/bin/python
## CNN model for tumor-specific CDR3 sequence prediction

import sys,os,re,csv,pathlib
import tensorflow as tf
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
#sys.path.insert(0,'/Users/bo/Library/Python/2.7/bin')
#import skimage

tf.logging.set_verbosity(tf.logging.ERROR)
AAs=np.array(list('WFGAVILMPYSTNQCKRHDE'))
curPath=os.getcwd()
##AAidx_file='AAindexNormalized.txt' ## AA index reached AUC about 61% for L=14. Worse than AdaBoost
##AAidx_file='AtchleyFactors.txt'  ## Atchley factors work worse than using 544 AA index
AAidx_file='AAidx_PCA.txt' ## works like a charm!!!
gg=open(AAidx_file)
AAidx_Names=gg.readline().strip().split('\t')
AAidx_Dict={}
for ll in gg.readlines():
    ll=ll.strip().split('\t')
    AA=ll[0]
    tag=0
    vv=[]
    for xx in ll[1:]:
        vv.append(float(xx))
    if tag==1:
        continue
    AAidx_Dict[AA]=vv

Nf=len(AAidx_Dict['C'])

pat=re.compile('[\\*_XB]')  ## non-productive CDR3 patterns

def OneHotEncoding(Seq):
    Seq_aa=list(Seq)
    Ns=len(Seq_aa)
    OHE=np.zeros([20,Ns])
    for ii in range(Ns):
        aa=Seq_aa[ii]
        vv=np.where(AAs==aa)
        OHE[vv,ii]=1
    OHE=OHE.astype(np.float32)
    return OHE

def AAindexEncoding(Seq):
    Ns=len(Seq)
    AAE=np.zeros([Ns, Nf])
    for kk in range(Ns):
        ss=Seq[kk]
        AAE[kk,]=AAidx_Dict[ss]
    AAE=np.transpose(AAE.astype(np.float32))
    return AAE

def GetFeatureLabels(TumorCDR3s, NonTumorCDR3s):
    nt=len(TumorCDR3s)
    nc=len(NonTumorCDR3s)
    LLt=[len(ss) for ss in TumorCDR3s]
    LLt=np.array(LLt)
    LLc=[len(ss) for ss in NonTumorCDR3s]
    LLc=np.array(LLc)
    NL=range(12,17)
    FeatureDict={}
    LabelDict={}
    for LL in NL:
        vvt=np.where(LLt==LL)[0]
        vvc=np.where(LLc==LL)[0]
        Labels=[1]*len(vvt)+[0]*len(vvc)
        Labels=np.array(Labels)
        Labels=Labels.astype(np.int32)
        data=[]
        for ss in TumorCDR3s[vvt]:
            if len(pat.findall(ss))>0:
                continue
            data.append(AAindexEncoding(ss))
#            data.append(OneHotEncoding(ss))
        for ss in NonTumorCDR3s[vvc]:
            if len(pat.findall(ss))>0:
                continue
            data.append(AAindexEncoding(ss))
#            data.append(OneHotEncoding(ss))
        data=np.array(data)
        features={'x':data,'LL':LL}
        FeatureDict[LL]=features
        LabelDict[LL]=Labels
    return FeatureDict, LabelDict

def cnn_model_CDR3_LL12(features, labels, mode):
  """Model function for CNN."""
  # Input Layer
  data=features['x']
#  LL=features['LL']
  input_layer = tf.reshape(data, [-1, Nf, 12, 1])
  # Convolutional Layer #1
  conv1 = tf.layers.conv2d(
      inputs=input_layer,
      filters=8,
      kernel_size=[Nf,2],
      padding="valid",
      activation=tf.nn.relu)
  # Pooling Layer #1
  pool1 = tf.layers.max_pooling2d(inputs=conv1, pool_size=[1, 2], strides=[1,1])  ## stride used to be 2
  # Convolutional Layer #2 and Pooling Layer #2
  conv2 = tf.layers.conv2d(
      inputs=pool1,
      filters=16,
      kernel_size=[1,2],
      padding="valid",
      activation=tf.nn.relu)
  pool2 = tf.layers.max_pooling2d(inputs=conv2, pool_size=[1, 2], strides=[1,1])
  SHAPE=pool2.shape
  pool2_flat = tf.reshape(pool2, [-1, int(SHAPE[1]*SHAPE[2]*SHAPE[3])])
  dense = tf.layers.dense(inputs=pool2_flat, units=10, activation=tf.nn.relu)
  dropout = tf.layers.dropout(
      inputs=dense, rate=0.4, training=mode == tf.estimator.ModeKeys.TRAIN)
  # Logits Layer
  logits = tf.layers.dense(inputs=dropout, units=2)
  predictions = {
      # Generate predictions (for PREDICT and EVAL mode)
      "classes": tf.argmax(input=logits, axis=1),
      # Add `softmax_tensor` to the graph. It is used for PREDICT and by the
      # `logging_hook`.
      "probabilities": tf.nn.softmax(logits, name="softmax_tensor")
  }
  if mode == tf.estimator.ModeKeys.PREDICT:
    return tf.estimator.EstimatorSpec(mode=mode, predictions=predictions)
  # Calculate Loss (for both TRAIN and EVAL modes)
  loss = tf.losses.sparse_softmax_cross_entropy(labels=labels, logits=logits)
  # Configure the Training Op (for TRAIN mode)
  if mode == tf.estimator.ModeKeys.TRAIN:
    optimizer = tf.train.GradientDescentOptimizer(learning_rate=0.001)
    train_op = optimizer.minimize(
        loss=loss,
        global_step=tf.train.get_global_step())
    return tf.estimator.EstimatorSpec(mode=mode, loss=loss, train_op=train_op)
  # Add evaluation metrics (for EVAL mode)
  eval_metric_ops = {
      "accuracy": tf.metrics.accuracy(
          labels=labels, predictions=predictions["classes"])
  }
  return tf.estimator.EstimatorSpec(
      mode=mode, loss=loss, eval_metric_ops=eval_metric_ops)

def cnn_model_CDR3_LL13(features, labels, mode):
  """Model function for CNN."""
  # Input Layer
  data=features['x']
#  LL=features['LL']
  input_layer = tf.reshape(data, [-1, Nf, 13, 1])
  # Convolutional Layer #1
  conv1 = tf.layers.conv2d(
      inputs=input_layer,
      filters=8,
      kernel_size=[Nf,2],
      padding="valid",
      activation=tf.nn.relu)
  # Pooling Layer #1
  pool1 = tf.layers.max_pooling2d(inputs=conv1, pool_size=[1, 2], strides=[1,1])  ## stride used to be 2
  # Convolutional Layer #2 and Pooling Layer #2
  conv2 = tf.layers.conv2d(
      inputs=pool1,
      filters=16,
      kernel_size=[1,2],
      padding="valid",
      activation=tf.nn.relu)
  pool2 = tf.layers.max_pooling2d(inputs=conv2, pool_size=[1, 2], strides=[1,1])
  SHAPE=pool2.shape
  pool2_flat = tf.reshape(pool2, [-1, int(SHAPE[1]*SHAPE[2]*SHAPE[3])])
  dense = tf.layers.dense(inputs=pool2_flat, units=10, activation=tf.nn.relu)
  dropout = tf.layers.dropout(
      inputs=dense, rate=0.4, training=mode == tf.estimator.ModeKeys.TRAIN)
  # Logits Layer
  logits = tf.layers.dense(inputs=dropout, units=2)
  predictions = {
      # Generate predictions (for PREDICT and EVAL mode)
      "classes": tf.argmax(input=logits, axis=1),
      # Add `softmax_tensor` to the graph. It is used for PREDICT and by the
      # `logging_hook`.
      "probabilities": tf.nn.softmax(logits, name="softmax_tensor")
  }
  if mode == tf.estimator.ModeKeys.PREDICT:
    return tf.estimator.EstimatorSpec(mode=mode, predictions=predictions)
  # Calculate Loss (for both TRAIN and EVAL modes)
  loss = tf.losses.sparse_softmax_cross_entropy(labels=labels, logits=logits)
  # Configure the Training Op (for TRAIN mode)
  if mode == tf.estimator.ModeKeys.TRAIN:
    optimizer = tf.train.GradientDescentOptimizer(learning_rate=0.001)
    train_op = optimizer.minimize(
        loss=loss,
        global_step=tf.train.get_global_step())
    return tf.estimator.EstimatorSpec(mode=mode, loss=loss, train_op=train_op)
  # Add evaluation metrics (for EVAL mode)
  eval_metric_ops = {
      "accuracy": tf.metrics.accuracy(
          labels=labels, predictions=predictions["classes"])
  }
  return tf.estimator.EstimatorSpec(
      mode=mode, loss=loss, eval_metric_ops=eval_metric_ops)

def cnn_model_CDR3_LL14(features, labels, mode):
  """Model function for CNN."""
  # Input Layer
  data=features['x']
#  LL=features['LL']
  input_layer = tf.reshape(data, [-1, Nf, 14, 1])
  # Convolutional Layer #1
  conv1 = tf.layers.conv2d(
      inputs=input_layer,
      filters=8,
      kernel_size=[Nf,2],
      padding="valid",
      activation=tf.nn.relu)
  # Pooling Layer #1
  pool1 = tf.layers.max_pooling2d(inputs=conv1, pool_size=[1, 2], strides=[1,1])  ## stride used to be 2
  # Convolutional Layer #2 and Pooling Layer #2
  conv2 = tf.layers.conv2d(
      inputs=pool1,
      filters=16,
      kernel_size=[1,2],
      padding="valid",
      activation=tf.nn.relu)
  pool2 = tf.layers.max_pooling2d(inputs=conv2, pool_size=[1, 2], strides=[1,1])
  SHAPE=pool2.shape
  pool2_flat = tf.reshape(pool2, [-1, int(SHAPE[1]*SHAPE[2]*SHAPE[3])])
  dense = tf.layers.dense(inputs=pool2_flat, units=10, activation=tf.nn.relu)
  dropout = tf.layers.dropout(
      inputs=dense, rate=0.4, training=mode == tf.estimator.ModeKeys.TRAIN)
  # Logits Layer
  logits = tf.layers.dense(inputs=dropout, units=2)
  predictions = {
      # Generate predictions (for PREDICT and EVAL mode)
      "classes": tf.argmax(input=logits, axis=1),
      # Add `softmax_tensor` to the graph. It is used for PREDICT and by the
      # `logging_hook`.
      "probabilities": tf.nn.softmax(logits, name="softmax_tensor")
  }
  if mode == tf.estimator.ModeKeys.PREDICT:
    return tf.estimator.EstimatorSpec(mode=mode, predictions=predictions)
  # Calculate Loss (for both TRAIN and EVAL modes)
  loss = tf.losses.sparse_softmax_cross_entropy(labels=labels, logits=logits)
  # Configure the Training Op (for TRAIN mode)
  if mode == tf.estimator.ModeKeys.TRAIN:
    optimizer = tf.train.GradientDescentOptimizer(learning_rate=0.001)
    train_op = optimizer.minimize(
        loss=loss,
        global_step=tf.train.get_global_step())
    return tf.estimator.EstimatorSpec(mode=mode, loss=loss, train_op=train_op)
  # Add evaluation metrics (for EVAL mode)
  eval_metric_ops = {
      "accuracy": tf.metrics.accuracy(
          labels=labels, predictions=predictions["classes"])
  }
  return tf.estimator.EstimatorSpec(
      mode=mode, loss=loss, eval_metric_ops=eval_metric_ops)


def cnn_model_CDR3_LL15(features, labels, mode):
  """Model function for CNN."""
  # Input Layer
  data=features['x']
#  LL=features['LL']
  input_layer = tf.reshape(data, [-1, Nf, 15, 1])
  # Convolutional Layer #1
  conv1 = tf.layers.conv2d(
      inputs=input_layer,
      filters=8,
      kernel_size=[Nf,2],
      padding="valid",
      activation=tf.nn.relu)
  # Pooling Layer #1
  pool1 = tf.layers.max_pooling2d(inputs=conv1, pool_size=[1, 2], strides=[1,1])  ## stride used to be 2
  # Convolutional Layer #2 and Pooling Layer #2
  conv2 = tf.layers.conv2d(
      inputs=pool1,
      filters=16,
      kernel_size=[1,2],
      padding="valid",
      activation=tf.nn.relu)
  pool2 = tf.layers.max_pooling2d(inputs=conv2, pool_size=[1, 2], strides=[1,1])
  SHAPE=pool2.shape
  pool2_flat = tf.reshape(pool2, [-1, int(SHAPE[1]*SHAPE[2]*SHAPE[3])])
  dense = tf.layers.dense(inputs=pool2_flat, units=10, activation=tf.nn.relu)
  dropout = tf.layers.dropout(
      inputs=dense, rate=0.4, training=mode == tf.estimator.ModeKeys.TRAIN)
  # Logits Layer
  logits = tf.layers.dense(inputs=dropout, units=2)
  predictions = {
      # Generate predictions (for PREDICT and EVAL mode)
      "classes": tf.argmax(input=logits, axis=1),
      # Add `softmax_tensor` to the graph. It is used for PREDICT and by the
      # `logging_hook`.
      "probabilities": tf.nn.softmax(logits, name="softmax_tensor")
  }
  if mode == tf.estimator.ModeKeys.PREDICT:
    return tf.estimator.EstimatorSpec(mode=mode, predictions=predictions)
  # Calculate Loss (for both TRAIN and EVAL modes)
  loss = tf.losses.sparse_softmax_cross_entropy(labels=labels, logits=logits)
  # Configure the Training Op (for TRAIN mode)
  if mode == tf.estimator.ModeKeys.TRAIN:
    optimizer = tf.train.GradientDescentOptimizer(learning_rate=0.001)
    train_op = optimizer.minimize(
        loss=loss,
        global_step=tf.train.get_global_step())
    return tf.estimator.EstimatorSpec(mode=mode, loss=loss, train_op=train_op)
  # Add evaluation metrics (for EVAL mode)
  eval_metric_ops = {
      "accuracy": tf.metrics.accuracy(
          labels=labels, predictions=predictions["classes"])
  }
  return tf.estimator.EstimatorSpec(
      mode=mode, loss=loss, eval_metric_ops=eval_metric_ops)

def cnn_model_CDR3_LL16(features, labels, mode):
  """Model function for CNN."""
  # Input Layer
  data=features['x']
#  LL=features['LL']
  input_layer = tf.reshape(data, [-1, Nf, 16, 1])
  # Convolutional Layer #1
  conv1 = tf.layers.conv2d(
      inputs=input_layer,
      filters=8,
      kernel_size=[Nf,2],
      padding="valid",
      activation=tf.nn.relu)
  # Pooling Layer #1
  pool1 = tf.layers.max_pooling2d(inputs=conv1, pool_size=[1, 2], strides=[1,1])  ## stride used to be 2
  # Convolutional Layer #2 and Pooling Layer #2
  conv2 = tf.layers.conv2d(
      inputs=pool1,
      filters=16,
      kernel_size=[1,2],
      padding="valid",
      activation=tf.nn.relu)
  pool2 = tf.layers.max_pooling2d(inputs=conv2, pool_size=[1, 2], strides=[1,1])
  SHAPE=pool2.shape
  pool2_flat = tf.reshape(pool2, [-1, int(SHAPE[1]*SHAPE[2]*SHAPE[3])])
  dense = tf.layers.dense(inputs=pool2_flat, units=10, activation=tf.nn.relu)
  dropout = tf.layers.dropout(
      inputs=dense, rate=0.4, training=mode == tf.estimator.ModeKeys.TRAIN)
  # Logits Layer
  logits = tf.layers.dense(inputs=dropout, units=2)
  predictions = {
      # Generate predictions (for PREDICT and EVAL mode)
      "classes": tf.argmax(input=logits, axis=1),
      # Add `softmax_tensor` to the graph. It is used for PREDICT and by the
      # `logging_hook`.
      "probabilities": tf.nn.softmax(logits, name="softmax_tensor")
  }
  if mode == tf.estimator.ModeKeys.PREDICT:
    return tf.estimator.EstimatorSpec(mode=mode, predictions=predictions)
  # Calculate Loss (for both TRAIN and EVAL modes)
  loss = tf.losses.sparse_softmax_cross_entropy(labels=labels, logits=logits)
  # Configure the Training Op (for TRAIN mode)
  if mode == tf.estimator.ModeKeys.TRAIN:
    optimizer = tf.train.GradientDescentOptimizer(learning_rate=0.001)
    train_op = optimizer.minimize(
        loss=loss,
        global_step=tf.train.get_global_step())
    return tf.estimator.EstimatorSpec(mode=mode, loss=loss, train_op=train_op)
  # Add evaluation metrics (for EVAL mode)
  eval_metric_ops = {
      "accuracy": tf.metrics.accuracy(
          labels=labels, predictions=predictions["classes"])
  }
  return tf.estimator.EstimatorSpec(
      mode=mode, loss=loss, eval_metric_ops=eval_metric_ops)

ModelDict={12:cnn_model_CDR3_LL12,
            13:cnn_model_CDR3_LL13,
            14:cnn_model_CDR3_LL14,
            15:cnn_model_CDR3_LL15,
            16:cnn_model_CDR3_LL16}

def TrainAndEvaluate(TrainFeature, TrainLabels, EvalFeature, EvalLabels, STEPs=10000, ID='', dir_prefix='/tmp/'):
    ## Train CNN model:
    for LL in range(12,17):
        CDR3_classifier=tf.estimator.Estimator(model_fn=ModelDict[LL],model_dir=dir_prefix+'/CDR3_classifier_PCA_LL'+str(LL)+'_L2_k2f8d10_'+ID+'/')
        train_input_fn=tf.estimator.inputs.numpy_input_fn(
            x={'x':TrainFeature[LL]['x']},        
            y=TrainLabels[LL],
            batch_size=100,
            num_epochs=None,
            shuffle=True)
        CDR3_classifier.train(input_fn=train_input_fn,steps=STEPs)
        eval_input_fn=tf.estimator.inputs.numpy_input_fn(
            x={'x':EvalFeature[LL]['x']},        
            y=EvalLabels[LL],
            num_epochs=1,
            shuffle=False)
        eval_results=CDR3_classifier.evaluate(input_fn=eval_input_fn)
        print(eval_results)
        
def PredictEvaluation(EvalFeature,EvalLabels=None,makePlot=False,ID='',dir_prefix=curPath+'/tmp/'):
    PredictClass={}
    PredictLabels={}
    AUCDict={}
    for LL in range(12,17):
        CDR3_classifier=tf.estimator.Estimator(model_fn=ModelDict[LL],model_dir=dir_prefix+'/CDR3_classifier_PCA_LL'+str(LL)+'_L2_k2f8d10_'+ID+'/')
        if EvalLabels is None:
            eval_input_fn=tf.estimator.inputs.numpy_input_fn(
                x={'x':EvalFeature[LL]['x']},        
                num_epochs=1,
                shuffle=False)
            eval_results=CDR3_classifier.predict(input_fn=eval_input_fn)
            xx=[]
            for x in eval_results:
                xx.append(x['probabilities'][1])
            AUC=None
            YY=None
        else:
            eval_input_fn=tf.estimator.inputs.numpy_input_fn(
                x={'x':EvalFeature[LL]['x']},        
                y=EvalLabels[LL],
                num_epochs=1,
                shuffle=False)
            eval_results=CDR3_classifier.predict(input_fn=eval_input_fn)
            xx=[]
            for x in eval_results:
                xx.append(x['probabilities'][1])
            YY=EvalLabels[LL]
            xy=list(zip(xx,YY))
            xy.sort()
            xs=[x for x,y in xy]
            ys=[y for x,y in xy]
            AUC=roc_auc_score(ys,xs)
        PredictClass[LL]=xx
        AUCDict[LL]=AUC
        PredictLabels[LL]=YY
    if makePlot:
        LLcolors=['b','g','r','c','m']
        LegendLabels=[]
        plt.figure(figsize=(7,7))
        font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 22}
        mpl.rc('font', **font)
        hhList=[]
        for LL in range(12,17):
            xx=PredictClass[LL]
            yy=PredictLabels[LL]
            ycurve=roc_curve(yy,xx)
            hh,=plt.plot(ycurve[0],ycurve[1],LLcolors[LL-12],lw=2)
            hhList.append(hh)
            LegendLabels.append(str(LL)+' ('+str(np.round(AUCDict[LL],2))+')')
        plt.plot([0,1],[0,1],ls='dashed',lw=2)
        plt.xlabel('False Positive Rate',fontsize=22)
        plt.ylabel('True Positive Rate',fontsize=22)
        legend=plt.legend(hhList,LegendLabels,fontsize=22,title='Length (AUC)')
        #plt.show()
        plt.savefig(dir_prefix+'/ROC-'+ID+'.png',dpi=300)
    return PredictClass, PredictLabels, AUCDict

def batchTrain(ftumor, fnormal,feval_tumor,feval_normal, rate=0.33,n=100,STEPs=10000,dir_prefix=curPath+'/tmp'):
    ## rate: cross validation ratio: 0.2 means 80% samples will be used for training
    ## n: number of subsamplings
    pathlib.Path(dir_prefix).mkdir(parents=True, exist_ok=True) 
    tumorCDR3s=[]
    g=open(ftumor)
    for ll in g.readlines():
        rr=ll.strip()
        if not rr.startswith('C') or not rr.endswith('F'):
            print("Non-standard CDR3s. Skipping.")
            continue
        tumorCDR3s.append(rr)
    normalCDR3s=[]
    g=open(fnormal)
    for ll in g.readlines():
        rr=ll.strip()
        if not rr.startswith('C') or not rr.endswith('F'):
            print("Non-standard CDR3s. Skipping.")
            continue
        normalCDR3s.append(rr)
    count=0
    nt=len(tumorCDR3s)
    nn=len(normalCDR3s)
    vt_idx=range(0,nt)
    vn_idx=range(0,nn)
    nt_s=int(np.ceil(nt*(1-rate)))
    nn_s=int(np.ceil(nn*(1-rate)))
    PredictClassList=[]
    PredictLabelList=[]
    AUCDictList=[]
    while count<n:
        print("==============Training cycle %d.=============" %(count))
        ID=str(count)
        vt_train=np.random.choice(vt_idx,nt_s,replace=False)
        vt_test=[x for x in vt_idx if x not in vt_train]
        vn_train=np.random.choice(vn_idx,nn_s,replace=False)
        vn_test=[x for x in vn_idx if x not in vn_train]
        sTumorTrain=np.array(tumorCDR3s)[vt_train]
        sNormalTrain=np.array(normalCDR3s)[vn_train]
        sTumorTest=np.array(tumorCDR3s)[vt_test]
        sNormalTest=np.array(normalCDR3s)[vn_test]
        ftrain_tumor=dir_prefix+'/sTumorTrain-'+str(ID)+'.txt'
        ftrain_normal=dir_prefix+'/sNormalTrain-'+str(ID)+'.txt'
        feval_tumor=dir_prefix+'/sTumorTest-'+str(ID)+'.txt'
        feval_normal=dir_prefix+'/sNormalTest-'+str(ID)+'.txt'
        h=open(ftrain_tumor,'w')
        _=[h.write(x+'\n') for x in sTumorTrain]
        h.close()
        h=open(ftrain_normal,'w')
        _=[h.write(x+'\n') for x in sNormalTrain]
        h.close()
        h=open(feval_tumor,'w')
        _=[h.write(x+'\n') for x in sTumorTest]
        h.close()
        h=open(feval_normal,'w')
        _=[h.write(x+'\n') for x in sNormalTest]
        h.close()
        g=open(ftrain_tumor)
        Train_Tumor=[]
        for line in g.readlines():
            Train_Tumor.append(line.strip())
        Train_Tumor=np.array(Train_Tumor)
        g=open(ftrain_normal)
        Train_Normal=[]
        for line in g.readlines():
            Train_Normal.append(line.strip())
        Train_Normal=np.array(Train_Normal)
        TrainFeature, TrainLabels=GetFeatureLabels(Train_Tumor,Train_Normal)
        g=open(feval_tumor)
        Eval_Tumor=[]
        for line in g.readlines():
            Eval_Tumor.append(line.strip())
        Eval_Tumor=np.array(Eval_Tumor)
        g=open(feval_normal)
        Eval_Normal=[]
        for line in g.readlines():
            Eval_Normal.append(line.strip())
        Eval_Normal=np.array(Eval_Normal)
        EvalFeature, EvalLabels=GetFeatureLabels(Eval_Tumor,Eval_Normal)
        TrainAndEvaluate(TrainFeature, TrainLabels, EvalFeature, EvalLabels,STEPs=STEPs,ID=ID,dir_prefix=dir_prefix)
        PC,PL,AD=PredictEvaluation(EvalFeature,EvalLabels=EvalLabels,makePlot=False,ID=ID,dir_prefix=dir_prefix)
        PredictClassList.append(PC)
        PredictLabelList.append(PL)
        AUCDictList.append(AD)
        count+=1
    return PredictClassList, PredictLabelList, AUCDictList

def PredictCancer(f,dir_prefix):
    ## f: input iSMART result file
    ## N: top N most frequent CDR3s will be included in the analysis
    gf=open(f)
    CDR3s=[]
    for ll in gf.readlines():
        cc=ll.strip().split('\t')[0]
        if not cc.startswith('C') or not cc.endswith('F'):
            continue
        CDR3s.append(cc)
    CDR3sDict={}
    for cc in CDR3s:
        if len(pat.findall(cc))>0:
            continue
        ll=len(cc)
        ccF=AAindexEncoding(cc)
        if ll not in CDR3sDict:
            CDR3sDict[ll]=[ccF]
        else:
            CDR3sDict[ll].append(ccF)
    ScoreDict={}
    XX=[]
    for LL in range(12,17):
        CDR3_classifier=tf.estimator.Estimator(model_fn=ModelDict[LL],model_dir=dir_prefix+'/CDR3_classifier_PCA_LL'+str(LL)+'_L2_k2f8d10_tCi01'+'/')
        if LL in CDR3sDict:
            eval_input_fn=tf.estimator.inputs.numpy_input_fn(
                x={'x':np.array(CDR3sDict[LL])},        
                num_epochs=1,
                shuffle=False)
        else:
            continue
        eval_results=CDR3_classifier.predict(input_fn=eval_input_fn)
        xx=[]
        for x in eval_results:
            xx.append(x['probabilities'][1])
        ScoreDict[LL]=xx
        XX+=xx
    mms=[]
    for kk in ScoreDict:
        mms.append((kk,np.mean(ScoreDict[kk])))
    CancerScore=np.mean(XX)
    return CancerScore,XX

#    return mms, XX, ScoreDict, CancerScore

def PredictBatch(DIR, dir_prefix=curPath+'/tmp/'):
    ffs=os.listdir(DIR)
    mmsList=[]
    SDList=[]
    XXList=[]
    for ff in ffs:
        mms, XX, SD=PredictCancer(DIR+ff, dir_prefix=dir_prefix)
        mmsList.append(mms)
        SDList.append(SD)
        XXList.append(XX)
    return ffs, mmsList, SDList, XXList
    
if len(sys.argv) > 1:   
 DIR=sys.argv[1]
 DIR1=os.path.basename(DIR)
 ffs=os.listdir(DIR)
 dir_prefix=sys.argv[2]
 CC=[]
 ffss=[]
 for ff in ffs:
   if ff == 'README.md':
     continue
   else: 
     score,XX1 = PredictCancer(DIR+'/'+ff, dir_prefix+'/tmp/')
     CC.append(score)
     ffss.append(ff)  
 CC=np.array(CC)
 ffss=np.array(ffss)
 if sys.argv[3] == '-t':
   with open('Cancer_score_'+DIR1+'.txt', 'w') as f:
      writer = csv.writer(f, delimiter='\t')
      writer.writerows(zip(ffss,CC))  
 elif sys.argv[3] == '-r':       
   with open('Cancer_score.txt', 'w') as f:
      writer = csv.writer(f, delimiter='\t')
      writer.writerows(zip(ffss,CC))  
  
  