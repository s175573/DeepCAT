## TCRboost
## Method for computational diagnosis of cancer using blood TCR repertoire data
## Aug 15th, 2018
## run this script using an R console in the same directory where the script is located. Put data/ folder in the same directory.

library(JOUSBoost)
library(pROC)
#library(seqinr)

#load('data/TumorSpecificCDR3.Rdata')  ## TCR beta chain CDR3s from TCGA
#load('data/VDJdbAll.Rdata')
#load('data/PvalMats_r.Rdata')
#load('data/tvModels_r.Rdata')
#load('data/AAindex.Rdata')

PrepareAdaptiveFile <- function(indir,outdir,thr=10000){
  ffs=dir(indir,full.names=T)
  for(ff in ffs){
    if(length(grep('\\.tsv',ff))==0)next
    ff0=unlist(strsplit(ff,'\\/'))
    ff0=ff0[length(ff0)]
    if(!file.exists(outdir))dir.create(outdir)
    newff=paste(outdir,'TestReal-',ff0,sep='')
    if(exists(newff))next
    ddnew=read.table(ff,header=T,sep='\t',stringsAsFactors=F)
    tmp.vv.rm=grep('\\*',ddnew[,2])
    if(length(tmp.vv.rm)>0)ddnew=ddnew[-tmp.vv.rm,]
    tmp.vv.rm=grep('X',ddnew[,2])
    if(length(tmp.vv.rm)>0)ddnew=ddnew[-tmp.vv.rm,]
    tmp.nn=nchar(ddnew[,2])
    tmp.vv=which(tmp.nn>=10 & tmp.nn<=24)
    #tmp.vv=which(tmp.nn>=1)
    ddnew=ddnew[tmp.vv,]
    ddnew=ddnew[which(ddnew[,6]!='unresolved'),]
    ddnew=ddnew[grep('^C.+F$',ddnew[,2]),]
    #ss=sum(ddnew[,3])
    #ddnew[,4]=ddnew[,3]/ss
    ddnew=ddnew[order(ddnew[,4],decreasing=T),]
    #q75=quantile(ddnew[,3],0.75)
    #THR=min(ddnew[,3]+1)
    #tmp.vv=which(ddnew[,3]>max(q75,THR))
    #if(length(tmp.vv)<=thr)tmp.vv=1:thr
    if(nrow(ddnew)>thr)ddnew=ddnew[1:thr,c(2,6,4)] else ddnew=ddnew[,c(2,6,4)]
    #ddnew=ddnew[1:5000,]
    ddnew=ddnew[1:min(thr,nrow(ddnew)),]
    write.table(ddnew,file=newff,quote=F,sep='\t',row.names=F)
  }
}

Prepare.iRepertoire.old <- function(indir, outdir, thr=10000){
  ffs=dir(indir,full.names=T)
  for(ff in ffs){
    if(length(grep('\\.csv',ff))==0)next
    ff0=unlist(strsplit(ff,'\\/'))
    ff0=ff0[length(ff0)]
    if(!file.exists(outdir))dir.create(outdir)
    newff=paste(outdir,'TestReal-',ff0,sep='')
    if(exists(newff))next
    ddnew=read.table(ff,header=T,sep=',',stringsAsFactors=F)
    ddnew=ddnew[order(ddnew[,4],decreasing=T),]
    if(nrow(ddnew)>thr)ddnew=ddnew[1:thr,c(1,2,4)] else ddnew=ddnew[,c(1,2,4)]
    ddnew[,2]=gsub('^h','',ddnew[,2])
    write.table(ddnew,file=newff,quote=F,sep='\t',row.names=F)
  }  
}

Prepare.iRepertoire <- function(indir, outdir, thr=5000){
  ## For RNA data, naive T cells are more likely to fluctuate with higher abundance
  ffs=dir(indir,full.names=T)
  for(ff in ffs){
    if(length(grep('\\.csv',ff))==0)next
    ff0=unlist(strsplit(ff,'\\/'))
    ff0=ff0[length(ff0)]
    if(!file.exists(outdir))dir.create(outdir)
    newff=paste(outdir,'TestReal-',ff0,sep='')
    if(exists(newff))next
    ddnew=read.table(ff,header=T,sep=',',stringsAsFactors=F)
    ddnew=ddnew[order(ddnew[,'copy'],decreasing=T),]
    tmp.vv=grep('\\*',ddnew[,1])
    if(length(tmp.vv)>0)ddnew=ddnew[-tmp.vv,]
    if(nrow(ddnew)>thr)ddnew=ddnew[1:thr,c('CDR3.pep.','V','copy')] else ddnew=ddnew[,c('CDR3.pep.','V','copy')]
    ddnew[,2]=gsub('^h','',ddnew[,2])
    ddnew[,1]=paste('C',ddnew[,1],'F',sep='')
    write.table(ddnew,file=newff,quote=F,sep='\t',row.names=F)
  }
}

PrepareAdaptiveFileEff <- function(indir,outdir,min.thr=4){
  ## Mar 25th, 2019 update: estimate effector cell using the min.thr*(minimum frequency)
  ffs=dir(indir,full.names=T)
  for(ff in ffs){
    if(length(grep('\\.tsv',ff))==0)next
    ff0=unlist(strsplit(ff,'\\/'))
    ff0=ff0[length(ff0)]
    if(!file.exists(outdir))dir.create(outdir)
    newff=paste(outdir,'TestReal-',ff0,sep='')
    if(exists(newff))next
    ddnew=read.table(ff,header=T,sep='\t',stringsAsFactors=F)
    tmp.vv.rm=grep('\\*',ddnew[,2])
    if(length(tmp.vv.rm)>0)ddnew=ddnew[-tmp.vv.rm,]
    tmp.vv.rm=grep('X',ddnew[,2])
    if(length(tmp.vv.rm)>0)ddnew=ddnew[-tmp.vv.rm,]
    tmp.nn=nchar(ddnew[,2])
    tmp.vv=which(tmp.nn>=10 & tmp.nn<=24)
    #tmp.vv=which(tmp.nn>=1)
    ddnew=ddnew[tmp.vv,]
    ddnew=ddnew[which(ddnew[,6]!='unresolved'),]
    ddnew=ddnew[grep('^C.+F$',ddnew[,2]),]
    #ss=sum(ddnew[,3])
    #ddnew[,4]=ddnew[,3]/ss
    ddnew=ddnew[order(ddnew[,4],decreasing=T),]
    #q75=quantile(ddnew[,3],0.75)
    #THR=min(ddnew[,3]+1)
    #tmp.vv=which(ddnew[,3]>max(q75,THR))
    #if(length(tmp.vv)<=thr)tmp.vv=1:thr
    MM=min(ddnew[,3])
    thr=min.thr*MM
    ddnew=ddnew[which(ddnew[,3]>thr),c(2,6,4)]
    if(nrow(ddnew)>=10000)ddnew=ddnew[1:10000,]
#    if(nrow(ddnew)>thr)ddnew=ddnew[1:thr,c(2,6,4)] else ddnew=ddnew[,c(2,6,4)]
    #ddnew=ddnew[1:5000,]
#    ddnew=ddnew[1:min(thr,nrow(ddnew)),]
    write.table(ddnew,file=newff,quote=F,sep='\t',row.names=F)
  }
}

run_iSMART <- function(indir,outdir){
  ## Place iSMART python source code, and related data files in the same directory as TCRboost
  ffs=dir(indir,full.names=T)
  if(!file.exists(outdir))dir.create(outdir)
  #cmd=paste('python iSMARTm.py -d ',indir,' -v -o ',outdir,sep='')
  cmd=paste('python iSMARTm.py -d ',indir,' -o ',outdir,sep='') ## Changed on Oct 11th
  #cmd=paste('python iSMARTm.py -d ',indir,' -v -o ',outdir,sep='') ## Changed on Mar 26, 2019
  system(cmd)
}

run_iSMARTeff <- function(indir,outdir){
  ## Place iSMART python source code, and related data files in the same directory as TCRboost
  ffs=dir(indir,full.names=T)
  if(!file.exists(outdir))dir.create(outdir)
  #cmd=paste('python iSMARTm.py -d ',indir,' -v -o ',outdir,sep='')
  cmd=paste('python iSMARTm.py -d ',indir,' -o ',outdir,sep='') ## Changed on Oct 11th, 2018
  system(cmd)
}

PrepareInputFiles <- function(indir,is.adaptive=TRUE,is.irepertoire=FALSE){
  iSMART_input_dir=paste(indir,'/for_iSMART_Mar15/',sep='')
  iSMART_output_dir=paste(indir,'/iSMART_results_Mar15/',sep='')
  if(is.adaptive){
    cat('Converting input file for iSMART.')
    PrepareAdaptiveFile(indir,iSMART_input_dir)
    cat('done! \n')
    cat('Running iSMART. ')
    run_iSMART(iSMART_input_dir,iSMART_output_dir)
    cat('done! \n')
  }
  if(is.irepertoire){
    cat('Converting input file for iSMART.')
    Prepare.iRepertoire(indir,iSMART_input_dir)
    cat('done! \n')
    cat('Running iSMART. ')
    run_iSMART(iSMART_input_dir,iSMART_output_dir)
    cat('done! \n')    
  }
}

AnalyzeShugay <- function(indir){
  iSMART_input_dir=paste(indir,'/for_iSMART_Mar15/',sep='')
  if(!file.exists(iSMART_input_dir))dir.create(iSMART_input_dir)
  iSMART_output_dir=paste(indir,'/iSMART_results_Mar15/',sep='')
  if(!file.exists(iSMART_output_dir))dir.create(iSMART_output_dir)
  ffs=dir(indir,full.names=T)
  for(ff in ffs){
    if(length(grep('txt',ff))==0)next
    tmp=unlist(strsplit(ff,'/'))
    ff0=gsub('\\.txt','',tmp[length(tmp)])
    ffnew=paste(iSMART_input_dir,ff0,'-top10000.txt',sep='')
    dd=read.table(ff,header=T,sep='\t',stringsAsFactors=F)
    tmp.vv=grep('[_XB]{1,3}',dd[,4])
    if(length(tmp.vv)>0)dd=dd[-tmp.vv,]
    write.table(dd[1:10000,c(4,5,1,2)],file=ffnew,sep='\t',quote=F,row.names=F)
  }
  run_iSMART(iSMART_input_dir,iSMART_output_dir)
}

PrepareInputFilesEff <- function(indir,is.adaptive=TRUE){
  iSMART_input_dir=paste(indir,'/for_iSMART_eff/',sep='')
  iSMART_output_dir=paste(indir,'/iSMART_results_eff/',sep='')
  if(is.adaptive){
    cat('Converting input file for iSMART.')
    PrepareAdaptiveFileEff(indir,iSMART_input_dir)
    cat('done! \n')
    cat('Running iSMART. ')
    run_iSMARTeff(iSMART_input_dir,iSMART_output_dir)
    cat('done! \n')
  }
}


TrainModel <- function(nn,sample.rate=0.5){
  ## The final output model has randomness due to the randomly sampled input, and the stochastic process of model training.
  ## This function has been pre-run to select the best predictive model. No need to run in applications.
  tumor.CDR3.all.nn=tumor.CDR3.all[which(nchar(tumor.CDR3.all)==nn)]
  NT=length(tumor.CDR3.all.nn)
  xx1=sample(tumor.CDR3.all.nn,NT*sample.rate)
  xx1.mat=c()
  for(kk in 3:(nn-3)){
    xx1.mat=cbind(xx1.mat,substr(xx1,kk,kk))
  }
  #control.CDR3s=VDJdb.all[which(VDJdb.all[,'Species']=='HomoSapiens'),'CDR3']
  control.CDR3s.nn=control.CDR3s[which(nchar(control.CDR3s)==nn)]
  NN=length(control.CDR3s.nn)
  yy=sample(control.CDR3s.nn,NN*sample.rate)
  yy.mat=c()
  for(kk in 3:(nn-3)){
    yy.mat=cbind(yy.mat,substr(yy,kk,kk))
  }
  
  tmp.PvalMat=c()
  for(kk in 1:ncol(yy.mat)){
    print(kk)
    for(ii in 1:nrow(AAindex.n)){
      x1=AAindex.n[ii,xx1.mat[,kk]]
      y1=AAindex.n[ii,yy.mat[,kk]]
      tmp.fc1=mean(x1,na.rm=T)/mean(y1,na.rm=T)
      tmp1=t.test(x1,y1)
      tmp.line=c(ii,kk+2,tmp.fc1,tmp1$p.value)
      tmp.PvalMat=rbind(tmp.PvalMat,tmp.line)
    }
  }
  
  XXa=c()
  
  for(kk in 3:(nn-4)){
    tmp.dd=tmp.PvalMat[which(tmp.PvalMat[,2]==kk),]
    tmp.vv=which(tmp.dd[,3]>=1.1)
    colnn=paste(kk,'_',rownames(AAindex.n)[tmp.vv],sep='')
    XX=c()
    for(ii in tmp.vv){
      XX=cbind(XX,AAindex.n[ii,xx1.mat[,kk-2]])
    }
    
    XX0=c()
    for(ii in tmp.vv){
      XX0=cbind(XX0,AAindex.n[ii,yy.mat[,kk-2]])
    }
    
    XX=as.matrix(XX)
    XX0=as.matrix(XX0)
    colnames(XX)=colnames(XX0)=colnn
    XXa=cbind(XXa,rbind(XX,XX0))
  }
  
  YY=c(rep(1,nrow(XX)),rep(-1,nrow(XX0)))
  
  adaboost(XXa,YY,n_rounds=50,verbose=T,tree_depth = 10)->tmp 
  
  test.CDR3s=tumor.CDR3.all.nn[which(!tumor.CDR3.all.nn %in% xx1)]
  yy2=control.CDR3s.nn[which(!control.CDR3s.nn %in% yy)]
  
  tt.mat=c()
  for(kk in 3:(nn-3)){
    tt.mat=cbind(tt.mat,substr(test.CDR3s,kk,kk))
  }
  yy2.mat=c()
  for(kk in 3:(nn-3)){
    yy2.mat=cbind(yy2.mat,substr(yy2,kk,kk))
  }
  
  XXa.test=c()
  for(kk in 3:(nn-4)){
    tmp.dd=tmp.PvalMat[which(tmp.PvalMat[,2]==kk),]
    tmp.vv=which(tmp.dd[,3]>=1.1)
    XX=c()
    for(ii in tmp.vv){
      XX=cbind(XX,AAindex.n[ii,tt.mat[,kk-2]])
    }
    
    XX0=c()
    for(ii in tmp.vv){
      XX0=cbind(XX0,AAindex.n[ii,yy2.mat[,kk-2]])
    }
    
    XX=as.matrix(XX)
    XX0=as.matrix(XX0)
    XXa.test=cbind(XXa.test,rbind(XX,XX0))
  }
  colnames(XXa.test)=colnames(XXa)
  YY.test=c(rep(1,nrow(tt.mat)),rep(-1,nrow(yy2.mat)))
  
  predict(tmp,XXa.test)->tmp.p
  
  return(list(Model=tmp,CVmat=table(tmp.p,YY.test),PvalMat=tmp.PvalMat))
}

TrainModel.PCA <- function(nn,sample.rate=0.5){
  ## The final output model has randomness due to the randomly sampled input, and the stochastic process of model training.
  ## This function has been pre-run to select the best predictive model. No need to run in applications.
  tumor.CDR3.all.nn=tumor.CDR3.all[which(nchar(tumor.CDR3.all)==nn)]
  NT=length(tumor.CDR3.all.nn)
  xx1=sample(tumor.CDR3.all.nn,NT*sample.rate)
  xx1.mat=c()
  for(kk in 3:(nn-3)){
    xx1.mat=cbind(xx1.mat,substr(xx1,kk,kk))
  }
  #control.CDR3s=VDJdb.all[which(VDJdb.all[,'Species']=='HomoSapiens'),'CDR3']
  control.CDR3s.nn=control.CDR3s[which(nchar(control.CDR3s)==nn)]
  NN=length(control.CDR3s.nn)
  yy=sample(control.CDR3s.nn,NN*sample.rate)
  yy.mat=c()
  for(kk in 3:(nn-3)){
    yy.mat=cbind(yy.mat,substr(yy,kk,kk))
  }

  AAmat=na.omit(AAindex.n)
  pca=prcomp(t(AAmat))$x
  pca=pca[,1:15]  ## select top 15 PCs
  XX1=c()
  for(ii in 1:nrow(xx1.mat)){
    tmp.y=xx1.mat[ii,]
    tmp.pca=c()
    for(y in tmp.y){
      tmp.pca=c(tmp.pca,pca[y,])
    }
    XX1=rbind(XX1,tmp.pca)
  }

  XX0=c()
  for(ii in 1:nrow(yy.mat)){
    tmp.y=yy.mat[ii,]
    tmp.pca=c()
    for(y in tmp.y){
      tmp.pca=c(tmp.pca,pca[y,])
    }
    XX0=rbind(XX0,tmp.pca)
  }  
  YY=c(rep(1,nrow(XX1)),rep(-1,nrow(XX0)))
  XXa=rbind(XX1,XX0)
  
  adaboost(XXa,YY,n_rounds=50,verbose=T,tree_depth = 10)->tmp 
  
  test.CDR3s=tumor.CDR3.all.nn[which(!tumor.CDR3.all.nn %in% xx1)]
  yy2=control.CDR3s.nn[which(!control.CDR3s.nn %in% yy)]
  
  tt.mat=c()
  for(kk in 3:(nn-3)){
    tt.mat=cbind(tt.mat,substr(test.CDR3s,kk,kk))
  }
  yy2.mat=c()
  for(kk in 3:(nn-3)){
    yy2.mat=cbind(yy2.mat,substr(yy2,kk,kk))
  }
  
  XX1.test=c()
  for(ii in 1:nrow(tt.mat)){
    tmp.y=tt.mat[ii,]
    tmp.pca=c()
    for(y in tmp.y){
      tmp.pca=c(tmp.pca,pca[y,])
    }
    XX1.test=rbind(XX1.test,tmp.pca)
  }
  
  XX0.test=c()
  for(ii in 1:nrow(yy2.mat)){
    tmp.y=yy2.mat[ii,]
    tmp.pca=c()
    for(y in tmp.y){
      tmp.pca=c(tmp.pca,pca[y,])
    }
    XX0.test=rbind(XX0.test,tmp.pca)
  } 
  XXa.test=rbind(XX1.test, XX0.test)
  colnames(XXa.test)=colnames(XXa)
  YY.test=c(rep(1,nrow(tt.mat)),rep(-1,nrow(yy2.mat)))
  
  predict(tmp,XXa.test)->tmp.p
  
  return(list(Model=tmp,CVmat=table(tmp.p,YY.test)))
}

PredictCancer <- function(filename, nn, tvModel, PvalMat, min.n=5){
  COLnames=c()
  for(kk in 3:(nn-4)){
    tmp.dd=PvalMat[which(PvalMat[,2]==kk),]
    tmp.vv=which(tmp.dd[,3]>=1.1)
    colnn=paste(kk,'_',rownames(AAindex.n)[tmp.vv],sep='')
    COLnames=c(COLnames,colnn)
  }
  ddnew=read.table(filename,stringsAsFactors=F,header=T,sep='\t')
  ddnew=as.matrix(ddnew)
  rownames(ddnew)=ddnew[,1]
  tmp.xx=ddnew[,1]
  tmp.xx=tmp.xx[which(nchar(tmp.xx)==nn)]
  if(length(tmp.xx)<min.n)return()
  tmp.xx.mat=c()
  for(kk in 3:(nn-3)){
    tmp.xx.mat=cbind(tmp.xx.mat,substr(tmp.xx,kk,kk))
  }
  
  tmp.XXa=c()
  for(kk in 3:(nn-4)){
    tmp.dd=PvalMat[which(PvalMat[,2]==kk),]
    tmp.vv=which(tmp.dd[,3]>=1.1)
    XX=c()
    for(ii in tmp.vv){
      XX=cbind(XX,AAindex.n[ii,tmp.xx.mat[,kk-2]])
    }
    XX=as.matrix(XX)
    tmp.XXa=cbind(tmp.XXa,XX)
  }
  colnames(tmp.XXa)=COLnames
  rownames(tmp.XXa)=NULL
  tree_seq = seq_along(tvModel$alphas)
  
  f = 0
  for (i in tree_seq) {
    tree = tvModel$trees[[i]]
    tree$terms = tvModel$terms
    pred = as.integer(as.character(stats::predict(tree, data.frame(tmp.XXa), type = "class")))
    f = f + (tvModel$alphas[i] * pred)
  }
  tmp.xx.prob=1/(1 + exp(-2 * f))
  #tmp.xx.prob=JOUSBoost::predict.adaboost(tvModel,tmp.XXa,type='prob')
  names(tmp.xx.prob)=tmp.xx
  if('frequencyCount....' %in% colnames(ddnew))tmp.freq=as.numeric(ddnew[tmp.xx,'frequencyCount....']) else tmp.freq=NULL
  return(list(prob=tmp.xx.prob,freq=tmp.freq))
}

PredictCancer.PCA <- function(filename, nn, tvModel, min.n=5){
  ddnew=read.table(filename,stringsAsFactors=F,header=T,sep='\t')
  ddnew=as.matrix(ddnew)
  rownames(ddnew)=ddnew[,1]
  tmp.xx=ddnew[,1]
  tmp.xx=tmp.xx[which(nchar(tmp.xx)==nn)]
  if(length(tmp.xx)<min.n)return()
  tmp.xx.mat=c()
  for(kk in 3:(nn-3)){
    tmp.xx.mat=cbind(tmp.xx.mat,substr(tmp.xx,kk,kk))
  }
  
  tmp.XXa=c()
  for(ii in 1:nrow(tmp.xx.mat)){
    tmp.y=tmp.xx.mat[ii,]
    tmp.pca=c()
    for(y in tmp.y){
      tmp.pca=c(tmp.pca,pca[y,])
    }
    tmp.XXa=rbind(tmp.XXa,tmp.pca)
  }
  colnames(tmp.XXa)=rep(paste('PC',1:15,sep=''),ncol(tmp.xx.mat))
  rownames(tmp.XXa)=NULL
  tree_seq = seq_along(tvModel$alphas)
  
  f = 0
  for (i in tree_seq) {
    tree = tvModel$trees[[i]]
    tree$terms = tvModel$terms
    pred = as.integer(as.character(stats::predict(tree, data.frame(tmp.XXa), type = "class")))
    f = f + (tvModel$alphas[i] * pred)
  }
  tmp.xx.prob=1/(1 + exp(-2 * f))
  #tmp.xx.prob=JOUSBoost::predict.adaboost(tvModel,tmp.XXa,type='prob')
  names(tmp.xx.prob)=tmp.xx
  if('frequencyCount....' %in% colnames(ddnew))tmp.freq=as.numeric(ddnew[tmp.xx,'frequencyCount....']) else tmp.freq=NULL
  return(list(prob=tmp.xx.prob,freq=tmp.freq))
}

PredictCombined <- function(filename,min.L=12,max.L=16,min.n=5){
  xx.prob=0
  xx.prob.f=0
  xx.prob.list=c()
  xx.freq.list=c()
  count=0
  short.list=c()
  for(nn in min.L:max.L){
    tvModel=eval(parse(text=paste('tvModel',nn,sep='')))
    PvalMat=eval(parse(text=paste('PvalMat',nn,sep='')))
    #tvModel=tvModels_r[[nn-11]]
    #PvalMat=PvalMats_r[[nn-11]]
    tmp.xx=PredictCancer(filename,nn,tvModel,PvalMat,min.n)
    tmp.xx.prob=tmp.xx$prob
    tmp.xx.freq=tmp.xx$freq
    if(length(tmp.xx.prob)>0){
      xx.prob=xx.prob+mean(tmp.xx.prob)
      xx.prob.f=xx.prob.f+mean(tmp.xx.prob*tmp.xx.freq,na.rm=T)
      tmp=list(tmp.xx.prob)
      names(tmp)=nn
      xx.prob.list=c(xx.prob.list,tmp)
      tmp=list(tmp.xx.freq)
      names(tmp)=nn
      xx.freq.list=c(xx.freq.list,tmp)
      count=count+1
    }else{
      short.list=c(short.list,nn)
    }
  }
  xx.prob=xx.prob/count
  xx.prob.f=xx.prob.f/count
  return(list(score=xx.prob,score_f=xx.prob.f,data=xx.prob.list,freq=xx.freq.list,sl=short.list))
}

PredictCombined.PCA <- function(filename,min.L=12,max.L=16,min.n=5){
  xx.prob=0
  xx.prob.f=0
  xx.prob.list=c()
  xx.freq.list=c()
  count=0
  short.list=c()
  for(nn in min.L:max.L){
    tvModel=eval(parse(text=paste('tvModel',nn,'p',sep='')))
    #tvModel=tvModels_r[[nn-11]]
    #PvalMat=PvalMats_r[[nn-11]]
    tmp.xx=PredictCancer.PCA(filename,nn,tvModel,min.n)
    tmp.xx.prob=tmp.xx$prob
    tmp.xx.freq=tmp.xx$freq
    if(length(tmp.xx.prob)>0){
      xx.prob=xx.prob+mean(tmp.xx.prob)
      xx.prob.f=xx.prob.f+mean(tmp.xx.prob*tmp.xx.freq,na.rm=T)
      tmp=list(tmp.xx.prob)
      names(tmp)=nn
      xx.prob.list=c(xx.prob.list,tmp)
      tmp=list(tmp.xx.freq)
      names(tmp)=nn
      xx.freq.list=c(xx.freq.list,tmp)
      count=count+1
    }else{
      short.list=c(short.list,nn)
    }
  }
  xx.prob=xx.prob/count
  xx.prob.f=xx.prob.f/count
  return(list(score=xx.prob,score_f=xx.prob.f,data=xx.prob.list,freq=xx.freq.list,sl=short.list))
}

batchPrediction <- function(DIR,min.L=12,max.L=16,min.n=5,type=c('tv')){

  ## run these lines to let the predict function be happy  
  set.seed(111)
  dat = circle_data(n = 500)
  train_index = sample(1:500, 400)
  ada = adaboost(dat$X[train_index,], dat$y[train_index], tree_depth = 2, n_rounds = 10, verbose = TRUE)
  yhat = predict(ada, dat$X[-train_index, ])
  ## end
  
  ffs=dir(DIR,full.names=T)
  if(length(ffs)==0)return()
  scores=c()
  DataList=c()
  FreqList=c()
  for(ff in ffs){
    print(ff)
    if(type=='tv') tmp=PredictCombined(ff,min.L,max.L,min.n)
    if(type=='pca') tmp=PredictCombined.PCA(ff,min.L,max.L,min.n)
    scores=c(scores,tmp$score)
    DataList=c(DataList,list(tmp$data))
    FreqList=c(FreqList,list(tmp$freq))
  }
  names(scores)=ffs
  return(list(S=scores,DL=DataList,FREQ=FreqList))
}

SubsampleEntropy <- function(DIR, nsub=10000,ratio=NA,outDIR=NA,PAT='\\.tsv'){
  ffs=dir(DIR,full.names=T)
  require(entropy)
  EE=c()
  ff0=c()
  for(ff in ffs){
    if(length(grep(PAT,ff))==0)next
    print(ff)
    tmp.dd=read.table(ff,header=T,sep='\t',stringsAsFactors=F)
    if(nrow(tmp.dd)==0)next
    tmp.vv=grep('\\*',tmp.dd[,2])
    if(length(tmp.vv)>0)tmp.dd=tmp.dd[-tmp.vv,]
    tmp.vv=grep('X',tmp.dd[,2])
    if(length(tmp.vv)>0)tmp.dd=tmp.dd[-tmp.vv,]
    tmp.nn=nchar(tmp.dd[,2])
#    tmp.dd=tmp.dd[which(tmp.nn>1),]
    tmp.vv=which(tmp.nn>=10 & tmp.nn<=24)
    tmp.dd=tmp.dd[tmp.vv,]
    tmp.dd=tmp.dd[which(tmp.dd[,6]!='unresolved'),]
    tmp.dd=tmp.dd[grep('^C.+F$',tmp.dd[,2]),]
    tmp.dd=tmp.dd[order(tmp.dd[,4],decreasing=T),]
    tmp.dd[,2]=make.names(tmp.dd[,2],unique=T)
    tmp.sim=c()
    nrep=table(as.numeric(tmp.dd[,3]))
    nrep=nrep[which(as.numeric(names(nrep))>=1)]
    for(ii in 1:length(nrep)){
      tmp.vv=which(as.numeric(tmp.dd[,3])==as.numeric(names(nrep)[ii]))
      if(names(nrep)[ii]=="1"){
        tmp.sim=c(tmp.sim,tmp.vv) 
      }else{
        tmp.sim=c(tmp.sim,rep(tmp.vv,each=as.numeric(names(nrep)[ii])))
      } 
    }
    tmp.sim=tmp.dd[tmp.sim,c(2,6)]
    Ne=nsub
    if(!is.na(ratio))Ne=floor(ratio*nrow(tmp.sim))
    if(nrow(tmp.sim)<Ne)next
    tmp.sim1=tmp.sim[sample(1:nrow(tmp.sim),Ne),]
    tmp.tt=sort(table(paste(tmp.sim1[,1],tmp.sim1[,2],sep=':')),decreasing=T)
    tmp.nns=strsplit(names(tmp.tt),':')
    tmp.cc=sapply(tmp.nns,function(x)x[1])
    tmp.Vg=sapply(tmp.nns,function(x)x[2])
    tmp.dd.ss=cbind(tmp.cc,tmp.Vg,tmp.tt)
    tmp.dd.ss[,1]=gsub('\\..+','',tmp.dd.ss[,1])
    if(is.na(ratio))ratio=round(nsub/nrow(tmp.sim),2)
    if(!is.na(outDIR)){
      f1=sapply(strsplit(ff,'\\/'),function(x)x[length(x)])
      f1=paste(outDIR,f1,'-SubSample-',ratio,'.tsv',sep='')
      write.table(tmp.dd.ss,file=f1,row.names=F,col.names=T,quote=F,sep='\t')
    }
    new.ff=as.numeric(tmp.dd.ss[,3])
    new.ff=new.ff/sum(new.ff)
    EE=c(EE,entropy:::entropy(new.ff))
    ff0=c(ff0,ff)
  }
  if(length(ff0)==0)return(c())
  ff0=sapply(strsplit(ff0,'\\/'),function(x)x[length(x)])
  ff0=gsub('\\.tsv','',ff0)
  names(EE)=ff0
  return(EE)
}

PrepareSubsample.iSMART <- function(indir){
  ffs=dir(indir,full.names=T)
  ffs=ffs[grep('\\.tsv',ffs)]
  outdir=paste(indir,'/for_iSMART/',sep='')
  if(!file.exists(outdir))dir.create(outdir)
  for(ff in ffs){
    ff0=sapply(strsplit(ff,'\\/'),function(x)x[length(x)])
    dd=read.table(ff,header=T,sep='\t',stringsAsFactors=F)
    dd0=dd[1:10000,]
    write.table(dd0,file=paste(outdir,ff0,sep=''),row.names=F,quote=F,sep='\t')
  }
  outdir2=paste(indir,'/iSMART_results/',sep='')
  if(!file.exists(outdir2))dir.create(outdir2)
  #cmd=paste('python iSMARTm.py -d ',indir,' -v -o ',outdir,sep='')
  cmd=paste('python iSMARTm.py -d ',outdir,' -o ',outdir2,sep='') ## Changed on Oct 11th
  system(cmd)
}


MovingAverage <- function(x, step=10, window=100){
  nx=length(x)
  ii=0
  MA=c()
  while(ii<=nx){
    ii1=ii+window
    if(ii1>nx)ii1=nx
    xx=x[ii:ii1]
    MA=c(MA,mean(xx))
    ii=ii+step
  }
  return(MA)
}

