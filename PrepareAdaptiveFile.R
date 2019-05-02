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
