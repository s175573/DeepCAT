#!usr/bin/python
## Pairwise distance estimation for homologous and non-homologous TCR CDR3 sequences
## immuno-Similarity Measurement by Aligning Receptors of T cell (iSMART)
## Apr 14th, 2018, first version

## Ultra-fast alignment version for iSMART, based on motif finding, July 15th, 2018
## Motif-guided iSMART (MiSMART) is suitable for analyzing large scale (500K-1M) sequence datasets.

import numpy as np
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
from itertools import chain
from optparse import OptionParser
import random, time
from functools import partial
from multiprocessing import Pool
#from statsmodels.sandbox.stats.multicomp import multipletests as mlt
import sys,os,re,resource

t0=time.time()
sys.setrecursionlimit(1000000)

def SplitMotif(seq,m=6,gap=1,strip=True):
    ## gap is either 0 or 1
    if strip:
        ns=len(seq)
        if ns>=10:
            seq=seq[2:(ns-2)]
        else:
            return []
    ns=len(seq)
    if ns<=6:
        return []
    motifList=[seq[xx:(xx+m)] for xx in range(0,ns-m+1)]
    if gap==1:
        for ii in range(1,m):
            motifList+=[seq[xx:(xx+ii)]+'.'+seq[(xx+ii+1):(xx+m+1)] for xx in range(0,ns-m)]
    return motifList           

def IndexSeqByMotif(seqs,seqIDs,m=6,gap=1):
    Ns=len(seqs)
    seqDict={}
    for ii in range(0,Ns):
        ss=seqs[ii]
        MM=SplitMotif(ss,m=m,gap=gap)
        seqDict[seqIDs[ii]]=MM
    motifDict={}
    for kk in seqDict.keys():
        vv=seqDict[kk]
        for mm in vv:
            if mm in motifDict:
                motifDict[mm].append(kk)
            else:
                motifDict[mm]=[kk]
    motifDictNew={}
    for kk in motifDict:
        if len(motifDict[kk])==1:
            continue
        motifDictNew[kk]=motifDict[kk]
    return motifDictNew

def GenerateMotifGraph(mD,seqs,seqIDs):
    SeqShareGraph={}
    mDL={}
    for kk in mD:
        vv=mD[kk]
        LL=[]
        for v in vv:
            LL.append(len(seqs[v]))
        mDL[kk]=LL
    for kk in mD:
        vv=mD[kk]
        LL=mDL[kk]
        nv=len(vv)
        for ii in range(0,nv):
            id_1=vv[ii]
            L1=LL[ii]
            for jj in range(ii,nv):
                if jj==ii:
                    continue
                id_2=vv[jj]
                L2=LL[jj]
                if L2 != L1:
                    continue
                if id_1 not in SeqShareGraph:
                    SeqShareGraph[id_1]=[id_2]
                elif id_2 not in SeqShareGraph[id_1]:
                    SeqShareGraph[id_1].append(id_2)
                if id_2 not in SeqShareGraph:
                    SeqShareGraph[id_2]=[id_1]
                elif id_1 not in SeqShareGraph[id_2]:
                    SeqShareGraph[id_2].append(id_1)
    return SeqShareGraph

def dfs(graph, start):
    '''
    Non-resursive depth first search
    '''
    visited = set()
    stack = [start]
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            stack.extend(set(graph[vertex]) - visited)
    
    return visited

def IdentifyMotifCluster(SSG):
    ## Input SeqShareGraph dictionary representation of sparse matrix
    POS=list(SSG.keys())
    NP=len(POS)
    ClusterList=[]
    tmpL=list(chain(*ClusterList))
    count=0
##    def LoadComm(STACK,cur_ii):
##        if cur_ii in STACK:
##                return
##        else:
##                STACK.append(cur_ii)
##                vv=SSG[cur_ii]
##                for v in vv:
##                    #v_idx=POS.index(v)
##                    if v not in STACK:
##                        LoadComm(STACK,v)
##        return STACK
    for ii in POS:
        if ii not in tmpL:
#            STACK=LoadComm([],ii)
            STACK=dfs(SSG,ii)
            ClusterList.append(list(STACK))
            tmpL=list(chain(*ClusterList))
            count+=1
            if count % 200 ==0:
                print("    Solved %d clusters" %(count))
##    ClusterList_ss=[]
##    for cc in ClusterList:
##        CL=[]
##        for pp in cc:
##            CL.append(POS[pp])
##        ClusterList_ss.append(CL)
    return ClusterList

def ParseFa(fname):
    InputStr=open(fname).readlines()
    FaDict={}
    seq=''
    for line in InputStr:
        if line.startswith('>'):
            if len(seq)>0:
                FaDict[seqHead]=seq
                seq=''
            seqHead=line.strip()
        else:
            seq+=line.strip()
    if seqHead not in FaDict:
        FaDict[seqHead]=seq
    return FaDict

cur_dir=os.path.dirname(os.path.realpath(__file__))+'/'
print(cur_dir)

def PreCalculateVgeneDist(VgeneFa="Imgt_Human_TRBV.fasta"):
    ## Only run one time if needed
    FaDict=ParseFa(cur_dir+VgeneFa)
    VScore={}
    CDR1Dict={}
    CDR2Dict={}
    for kk in FaDict:
        if '|' in kk:
            VV=kk.split('|')[1]
        else:
            VV=kk[1:]
        CDR1Dict[VV]=FaDict[kk][26:37]  ## Imgt CDR1: 27 - 38
        CDR2Dict[VV]=FaDict[kk][55:64]  ## Imgt CDR2: 56 - 65
    Vkeys=[]    
    Vkeys=CDR1Dict.keys()
    nn=len(Vkeys)
    for ii in range(0,nn):
        V1=list(Vkeys)[ii]
        s1_CDR1=CDR1Dict[V1]
        s1_CDR2=CDR2Dict[V1]
        for jj in range(ii,nn):
            V2=list(Vkeys)[jj]
            s2_CDR1=CDR1Dict[V2]
            s2_CDR2=CDR2Dict[V2]
            score1=SeqComparison(s1_CDR1,s2_CDR1)
            score2=SeqComparison(s2_CDR2,s2_CDR2)
            #print score1+score2
            VScore[(V1,V2)]=score1+score2
    gg=open('VgeneScores.txt','w')
    for kk in VScore:
        vv=VScore[kk]
        line=kk[0]+'\t'+kk[1]+'\t'+str(vv)+'\n'
        gg.write(line)
    gg.close()

def GetCor(m1,m2):
    ## Given di amino acid motifs m1 and m2, get correlation
    if '-' in m1 or '-' in m2:
        return (0.1,0.9)
    if '*' in m1 or '*' in m2:
        return (0.05,0.95)  ## Arbitrary low score for gap
    if '.' in m1 or '.' in m2:
        return (0.05,0.95)  ## Arbitrary low score for gap  
    COR=DDcor[(m1,m2)] 
    return (COR,1)

def InsertGap(Seq,n):
    ## Insert n gaps to Seq; n<=2
    if n==0:
        return [Seq]
    ns=len(Seq)
    SeqList=[]
    if(n==1):
        for kk in range(0,ns+1):
            SeqNew=Seq[0:kk]+'-'+Seq[kk:]
            SeqList.append(SeqNew)
    if(n==2):
        for kk in range(0,ns+1):
            SeqNew=Seq[0:kk]+'-'+Seq[kk:]
            for jj in range(0,ns+2):
                SeqNew0=SeqNew[0:jj]+'-'+SeqNew[jj:]
                SeqList.append(SeqNew0)
    return SeqList

def SeqComparison(s1,s2,gap=-6):
    n=len(s1)
    CorList=[]
    score=0
    for kk in range(0,n):
        aa=s1[kk]
        bb=s2[kk]
        if aa in ['.','-','*'] or bb in ['.','-','*']:
            if aa!=bb:
                score += gap
            continue
        if aa==bb:
            score += min(4,blosum62[(aa,aa)])
            continue
        KEY=(aa,bb)
        if KEY not in blosum62:
            KEY=(bb,aa)
        if KEY not in blosum62:
            print(KEY)
            raise "Non-standard amino acid coding!"
        score+=blosum62[KEY]
    return score

def SeqComparison_Di(s1,s2,gap=-6):
    ## Older version that allows di amino acid replacement.
    n=len(s1)
    CorList=[]
    score=0
    for kk in range(0,n-1):
        m1=s1[kk:(kk+2)]
        m2=s2[kk:(kk+2)]
        (Cor,PP)=GetCor(m1,m2)
        CorList.append(Cor)
        aa=s1[kk]
        if kk==0:
            if Cor>=COR_0005:
                bb=s1[kk]
            else:
                bb=s2[kk]
        else:
            Cor1=CorList[kk-1]
            if Cor1>=COR_0005 or Cor>=COR_0005:
                bb=s1[kk]
            else:
                bb=s2[kk]
        if aa=='-' or bb=='-':
            score+= gap
            continue
        if aa=='*' or bb=='*':
            score+= gap
            continue
        if aa=='.' or bb=='.':
            if aa=='.' and bb=='.':
                continue
            else:
                score+=gap
                continue
        if aa==bb:
            score+= min(4,blosum62[(aa,aa)])
            continue
        KEY=(aa,bb)
        if KEY not in blosum62:
            KEY=(bb,aa)
        if KEY not in blosum62:
            print(KEY)
            raise "Non-standard amino acid coding!"
        score+=blosum62[KEY]
    aa=s1[n-1]
    bb=s2[n-1]
    if aa in ['.','-','*'] or bb in ['.','-','*']:
        if not (aa=='.' and bb=='.'):
            score+= gap
        else:
            score+=0
    else:
        if aa==bb:
#            score+= min(4,blosum62[(aa,aa)])
            score+=blosum62[(aa,aa)]
        else:
            KEY=(aa,bb)
            if KEY not in blosum62:
                KEY=(bb,aa)
            if KEY not in blosum62:
                print(KEY)
                raise "Non-standard amino acid coding!"
            score+=blosum62[KEY]
    return score

def NHLocalAlignment(Seq1,Seq2,gap_thr=1,gap=-6,Di=False):
    n1=len(Seq1)
    n2=len(Seq2)
    if n1<n2:
        Seq=Seq1
        Seq1=Seq2
        Seq2=Seq
        nn=n2-n1
    else:
        nn=n1-n2
    if nn>gap_thr:
        return -1
    #alns=pairwise2.align.localms(Seq1,Seq2,m,s,g,ge)
    SeqList1=[Seq1]
    SeqList2=InsertGap(Seq2,nn)
    alns=[]
    SCOREList=[]
    for s1 in SeqList1:
        for s2 in SeqList2:
            if Di:
                SCOREList.append(SeqComparison_Di(s1,s2,gap))
            else:
                SCOREList.append(SeqComparison(s1,s2,gap))
#            alns.append((s1,s2))
##    SCOREList=[]
##    for seq in SeqList:
##        s1=aln[0]
##        s2=aln[1]
##        SCORE=SeqComparison(s1,s2)
##        SCOREList.append(SCORE)
    maxS=max(SCOREList)
#    ALN=alns[np.where(np.array(SCOREList)==maxS)[0][0]]
    return maxS

def fun_map(p,f):
    ## Fake function for passing multiple arguments to Pool.map()
    return f(*p)

def falign(xx,st,VScore={}, Seqs=[], Vgene=[], UseV=True, gapn=1, gap=-6):
    ii=xx[0]
    jj=xx[1]
    V1=Vgene[ii]
    V2=Vgene[jj]
    mid1=Seqs[ii][st:-2]
    mid2=Seqs[jj][st:-2]
    if UseV:
        if V2==V1:
            V_score=4
        else:
            Vkey=(V1,V2)
            if Vkey not in VScore:
                Vkey=(V2,V1)
            if Vkey not in VScore:
                #print("V gene not found!")
                VScore=0
            else:
                V_score=VScore[Vkey]/20.0
    else:
        V_score=4.0
    aln=NHLocalAlignment(mid1,mid2,gapn,gap)
    score=aln/float(max(len(mid1),len(mid2)))+V_score
    return score

def PWalign(Seqs,Vgene=[],ID=[],VScore={}, gap=-6,gapn=1,UseV=True,cutoff=7,Nthread=1,Di=False):
    ## Wrapper function
    ns=len(Seqs)
    if ns != len(Vgene):
        if len(Vgene)==0:
            Vgene=['']*ns
            ID=range(0,ns)
        else:
            raise "Incompatible variable gene number!"
    z=sorted(zip(Seqs,Vgene,ID),key=lambda pair:len(pair[0]))
    Seqs=[x for x,y,t in z]
    Vgene=[x for y,x,t in z]
    ID=[x for t,y,x in z]
    del z
    PWscore={}
    st=4
    if not UseV:
        st=2
    t1=time.time()
    if Nthread==1:
        for ii in range(0,ns):
            V1=Vgene[ii]
            if ii % 100 ==0:
                t2=time.time()
#                print('%d: Time elapsed %f' %(ii, t2-t1))
            for jj in range(ii,ns):
                if ii==jj:
                    continue
                V2=Vgene[jj]
                mid1=Seqs[ii][st:-2]
                mid2=Seqs[jj][st:-2]
                if UseV:
                    if V2==V1:
                        V_score=4
                    else:
                        Vkey=(V1,V2)
                        if Vkey not in VScore:
                            Vkey=(V2,V1)
                        if Vkey not in VScore:
                            #print("V gene not found!")
                            continue
                        else:
                            V_score=VScore[Vkey]/20.0  ## Take the floor of the float number
                else:
                    V_score=4.0
                aln=NHLocalAlignment(mid1,mid2,gapn,gap,Di)
                #print aln
    #            J_score=NHLocalAlignment(Jend1,Jend2,gap=False)[0]
                score=aln/float(max(len(mid1),len(mid2)))+V_score
                if score>=cutoff:
                    PWscore[(ii,jj)]=1
    else:
        # Multi-thread processing
        p=Pool(Nthread)
        XX=[]
        for ii in range(0,ns):
            for jj in range(ii,ns):
                if ii==jj:
                    continue
                else:
                    XX.append([ii,jj])
        para= []
        for xx in XX:
            para.append((xx,st,VScore,Seqs, Vgene, UseV, gapn, gap))
        pl_out=p.map(partial(fun_map,f=falign),para)
        p.close()
        p.join()
        ## End multiple processing
        for kk in range(0,len(XX)):
            score=pl_out[kk]
            if score>=cutoff:
                PWscore[(XX[kk][0],XX[kk][1])]=1

    return (PWscore,Seqs,Vgene,ID)  

def IdentifyCDR3Clusters(PWscore,cutoff=7):
    POS=np.array(list(PWscore.keys()))[np.where(np.array(list(PWscore.values()))==1)]
    if len(POS)<=0:
       # print("Too few clustered CDR3s! Please check your repertoire data.")
        return []
    POS=list(POS)
    POS=np.array([list(map(lambda x:x[0],POS)),list(map(lambda x:x[1],POS))])
    uniquePos=list(set(list(POS[0])+list(POS[1])))
    ClusterList=[]
    tmpL=list(chain(*ClusterList))
    def LoadComm(STACK,cur_ii):
        if cur_ii in STACK:
                return
        else:
                STACK.append(cur_ii)
                vv=list(POS[1][np.where(POS[0]==cur_ii)])+list(POS[0][np.where(POS[1]==cur_ii)])
                for v in vv:
                        LoadComm(STACK,v)
        return STACK
    for ii in uniquePos:
        if ii in tmpL:
            continue
        else:
            STACK=LoadComm([],ii)
            ClusterList.append(STACK)
            tmpL=list(chain(*ClusterList))
    return ClusterList

def CompareClusters(CLinfo1, CLinfo2, VScore,gapn=1, gap=-6, cutoff=7,UseV=True,Di=False):
    CL1=CLinfo1[0]
    CL2=CLinfo2[0]
    Seqs1=CLinfo1[1]
    Seqs2=CLinfo2[1]
    Vgene1=CLinfo1[2]
    Vgene2=CLinfo2[2]
    n1=len(CL1)
    n2=len(CL2)
    #print "Processing %d * %d clusters" %(n1,n2)
    MergedCL=[]
    MergedSeq=[]
    MergedVgene=[]
    for ii in range(0,n1):
        #print ii
        seqs1=list(np.array(Seqs1)[CL1[ii]])
        VG1=list(np.array(Vgene1)[CL1[ii]])
        L1=np.median(list(map(lambda x:len(x), seqs1)))
        for jj in range(0,n2):
            seqs2=list(np.array(Seqs2)[CL2[jj]])
            VG2=list(np.array(Vgene2)[CL2[jj]])
            L2=np.median(list(map(lambda x:len(x), seqs2)))
            if L2<=L1-1 or L2>=L1+1:
                continue
            Scores=[]
            st=4
            if not UseV:
                st=2
            for tt1 in zip(seqs1,VG1):
                ss1=tt1[0]
                vv1=tt1[1]
                mid1=ss1[st:-2]
                for tt2 in zip(seqs2,VG2):
                    ss2=tt2[0]
                    vv2=tt2[1]
                    mid2=ss2[st:-2]
                    Score=NHLocalAlignment(mid1,mid2,gapn,gap,Di)/float(max(len(mid1),len(mid2)))
                    if UseV:
                        if vv1==vv2:
                            V_score=4
                        else:
                            Vkey=(vv1,vv2)
                            if Vkey not in VScore:
                                Vkey=(vv2,vv1)
                            if Vkey not in VScore:
                                #print("V gene not found!")
                                V_score=0
                            else:
                                V_score=VScore[Vkey]/20.0  ## Take the floor of the float number
                    else:
                        V_score=4.0
                    Score+=V_score
                    #print ss1, ss2, Score
                    Scores.append(Score)
            Scores_sorted=sorted(Scores,reverse=True)
            if Scores_sorted[0]>=cutoff and Scores_sorted[1]>=cutoff:
                #print [ii,jj]
                MergedCL.append([ii,jj])
                MergedSeq.append((seqs1,seqs2))
                MergedVgene.append((VG1,VG2))
    return (MergedCL,MergedSeq,MergedVgene)

def ObtainCL(InputFile,VScore, gap, gapn, cutoff=7, UseV=True, outDir='./',Nthread=1,Di=False):
    ff=open(InputFile)
    OutF=outDir+'/'+re.sub('.txt','',InputFile.split('/')[-1])
    Seqs0=[]
    Vgene0=[]
    ID=[]
    count=0
    ALLLines=ff.readlines()
    CDR3Dict={}
    for line in ALLLines[1:]: ## get rid of header line
        ww=line.strip().split('\t')
        Seqs0.append(ww[0])
        if UseV:
            Vgene0.append(ww[1])
        else:
            Vgene0.append('')
        CDR3Dict[count]=ww[1:]
        ID.append(count)
        count+=1
    print("Building motif index")
    mD=IndexSeqByMotif(Seqs0,ID)
    print("Generating motif sharing graph")
    SSG=GenerateMotifGraph(mD,Seqs0,ID)
    print("Dividing motif sharing graph")
    mClusters=IdentifyMotifCluster(SSG)
    g=open(OutF+'_ClusteredCDR3s_'+str(cutoff)+'.txt','w')
    g.write(ALLLines[0].strip()+'\t'+'Group'+'\n')
    gr=0
    CL=[]
    Seqs=[]
    Vgene=[]
    PWscore={}
    for mID in mClusters:
        mSeqs=[]
        mVgene=[]
        for mm in mID:
            mSeqs.append(Seqs0[mm])
            mVgene.append(Vgene0[mm])
#        print("  Processing %d sequences." %(len(mSeqs)))
        TMP=PWalign(mSeqs,mVgene,mID,VScore, gap, gapn, UseV=UseV,cutoff=cutoff,Nthread=Nthread,Di=Di)
        PWscore=TMP[0]
        Seqs=TMP[1]
        Vgene=TMP[2]
        ID=TMP[3]
        CL=IdentifyCDR3Clusters(PWscore,cutoff=cutoff)
        for cl in CL:
            gr+=1
            for ss in cl:
                cdr3=Seqs[ss]
                tmpID=ID[ss]
                Line='\t'.join([cdr3]+CDR3Dict[tmpID]+[str(gr)])+'\n'  ## Must add statistical significance estimation
                g.write(Line)
    g.close()
    return (CL,Seqs,Vgene,PWscore,CDR3Dict)

def ParseCLFile(CLfile):
    ff=open(CLfile)
    ALLLines=ff.readlines()
    Seqs=[]
    Vgene=[]
    CL=[]
    gp=1
    groups=[]
    count=0
    for line in ALLLines[1:]:
        ww=line.strip().split('\t')
        cdr3=ww[0]
        vv=ww[1]
        Seqs.append(cdr3)
        Vgene.append(vv)
        ID=int(ww[-1])
        if ID>gp:
            CL.append(groups)
            groups=[]
            gp+=1
        groups.append(count)
        count+=1
    CL.append(groups)    
    return (CL,Seqs,Vgene)

def CommandLineParser():
    parser=OptionParser()
    print('''
iSMART is a highly specific tools for dividing TCR beta chain repertoire sequencing
data into antigen-specific groups. Similarity between different repertoires is also
compared through commonly shared CDR3 groups. iSMART is developed by Li lab at UTSW.
All rights reserved.
Input columns:
1. CDR3 amino acid sequence (Starting from C, ending with the first F/L in motif [FL]G.G)
2. Variable gene name in Imgt format: TRBVXX-XX*XX
3. Joining gene name (optional)
4. Frequency (optional)
5. Other information (optional)
''')
    parser.add_option("-d","--directory",dest="Directory",help="Input repertoire sequencing file directory. Please make sure that all the files in the directory are input files.",default="")
    parser.add_option("-f","--file",dest="File",default='',help="Input single file of CDR3 sequences for grouping")
    parser.add_option("-F","--fileList",dest="files",default='',help='Alternative input: a file containing the full path to all the files. If given, overwrite -d and -f option')
    parser.add_option("-t","--threshold",dest="thr",default=7.5,help="Threshold for calling similar CDR3 groups. The higher the more specific.")
    parser.add_option("-o","--output",dest="OutDir",default='./',help="Output directory for intermediate and final outputs.")
    parser.add_option("-g","--GapPenalty",dest="Gap",default= -6,help="Gap penalty,default= -6")
    parser.add_option("-n","--GapNumber",dest="GapN",default=1,help="Maximum number of gaps allowed when performing alignment. Max=1, default=1")
    parser.add_option("-V","--VariableGeneFa",dest="VFa",default="Imgt_Human_TRBV.fasta",help="IMGT Human beta variable gene sequences")
    parser.add_option("-W","--KeepPairwiseMatrix",dest="PW",default=False,action="store_true",help="If true, iSMART will keep the pairwise alignment score matrix. Make sure you have enough disk space when dealing with large samples. Default: False")    
    parser.add_option("-I","--CrossInteraction",dest='I',default=False,action='store_true',help="If true, iSMART takes the clonal group files to compute sharing between individuals.")
    parser.add_option("-C","--CrossComparison",dest='C',default=False,action='store_true',help="If true, iSMART compares all the CDR3 clusters in the input to the directory specified in -r.")
    parser.add_option("-r","--referenceCohort",dest='R',default='',help="See -C option")
    parser.add_option("-v","--VariableGene",dest="V",default=True,action="store_false",help="If False, iSMART will omit variable gene information and use CDR3 sequences only. This will yield reduced specificity. The cut-off will automatically become the current value-4.0")
    parser.add_option("-N","--NumberOfThreads",dest="NN",default=1,help="Number of threads for multiple processing. Not working so well.")
#    parser.add_option("-D","--UseDiAAmat",dest="Di",default=False,action="store_true",help="If True, iSMART will use a predefined di-amino acid substitution matrix in sequence comparison.")
    return parser.parse_args()

def main():
    (opt,_)=CommandLineParser()
    FileDir=opt.Directory
    if len(FileDir)>0:
            files=os.listdir(FileDir)
            files0=[]
            for ff in files:
                if os.path.splitext(ff)[1] == ".tsv":
                    ff=FileDir+'/'+ff
                    files0.append(ff)    
            files=files0
    else:
            files=[]
    File=opt.File
    if len(File)>0:
            files=[File]
    FileList=opt.files
    if len(FileList)>0:
            files=[]
            fL=open(FileList)
            for ff in fL.readlines():
                    files.append(ff.strip())
    VFa=opt.VFa
    PreCalculateVgeneDist(VFa)
    vf=open('./VgeneScores.txt')  ## Use tcrDist's Vgene 80-score calculation
    VScore={}
    while 1:
        line=vf.readline()
        if len(line)==0:
            break
        ww=line.strip().split('\t')
        VScore[(ww[0],ww[1])]=int(ww[2])
    Gap=int(opt.Gap)
    Gapn=int(opt.GapN)
    cutoff=float(opt.thr)
    OutDir=opt.OutDir
    PW=opt.PW
    II=opt.I
    CC=opt.C
    RR=opt.R
    VV=opt.V
    NN=int(opt.NN)
    Di=False
    DataDict={}
    if CC:
        print("Compare input file with reference data")
        RefFiles=os.listdir(RR)
        gg=open(OutDir+"CrossReference.txt",'w')
        gg.write("CDR3\tVgene\tIndividualGroupID\tCrossGroupID\tSampleID\n")
        for f1 in files:
            TMP1=ParseCLFile(ff)
            print("Processing %s" %(f1))
            for fr in RefFiles:
                TMPr=ParseCLFile(RR+fr)
                sys.stdout.write('.')
                sys.stdout.flush()
                MC=CompareClusters(TMP1,TMPr,VScore,Gapn,Gap,cutoff,VV,Di)
                n=len(MC[0])
                for kk in range(0,n):
                    gg1=MC[0][kk][0]
                    gg2=MC[0][kk][1]
                    ww1=MC[1][kk][0]
                    ww2=MC[1][kk][1]
                    vv1=MC[2][kk][0]
                    vv2=MC[2][kk][1]
                    nw1=len(ww1)
                    nw2=len(ww2)
                    for ss in range(0,nw1):
                        line=ww1[ss]+'\t'+vv1[ss]+'\t'+str(gg1)+'\t'+str(groupID)+'\t'+files[ii]+'\n'
                        gg.write(line)
                    for ss in range(0,nw2):
                        line=ww2[ss]+'\t'+vv2[ss]+'\t'+str(gg2)+'\t'+str(groupID)+'\t'+files[jj]+'\n'
                        gg.write(line)
                    groupID+=1
            print('')
        gg.close()
        return 
    for ff in files:
        print("Processing %s" %(ff))
        if II:
            TMP=ParseCLFile(ff)
        else:
            TMP=ObtainCL(ff,VScore,Gap,Gapn,cutoff,VV, OutDir, NN,Di)
        if PW:
            PWscore=TMP[0]
            OutF=OutDir+ff+"_PWscores.txt"
            gg=open(OutF,'w')
            line='\t'.join(TMP[1])+'\n'
            gg.write(line)
            line='\t'.join(TMP[2])+'\n'
            gg.write(line)
            for ii in range(0,len(PWscore)):
                line='\t'.join(map(str,list(PWscore[ii])))+'\n'
                gg.write(line)
            gg.close()
        DataDict[ff]=TMP
##    gg=open(OutDir+"CrossComparison.txt",'w') ## Must add statistical significance estimation later
##    gg.write("CDR3\tVgene\tIndividualGroupID\tCrossGroupID\tSampleID\n")
##    if len(files)<= -1:
##        nn=len(files)
##        groupID=0
##        print "Pairwise comparison of %d repertoires" %(nn)
##        for ii in xrange(0,nn):
##            print ii
##            CLinfo1=DataDict[files[ii]]
##            for jj in xrange(ii,nn):
##                if jj==ii:
##                    continue
##                CLinfo2=DataDict[files[jj]]
##                MC=CompareClusters(CLinfo1,CLinfo2,VScore,Gapn,Gap,cutoff)
##                n=len(MC[0])
##                for kk in xrange(0,n):
##                    gg1=MC[0][kk][0]
##                    gg2=MC[0][kk][1]
##                    ww1=MC[1][kk][0]
##                    ww2=MC[1][kk][1]
##                    vv1=MC[2][kk][0]
##                    vv2=MC[2][kk][1]
##                    nw1=len(ww1)
##                    nw2=len(ww2)
##                    for ss in xrange(0,nw1):
##                        line=ww1[ss]+'\t'+vv1[ss]+'\t'+str(gg1)+'\t'+str(groupID)+'\t'+files[ii]+'\n'
##                        gg.write(line)
##                    for ss in xrange(0,nw2):
##                        line=ww2[ss]+'\t'+vv2[ss]+'\t'+str(gg2)+'\t'+str(groupID)+'\t'+files[jj]+'\n'
##                        gg.write(line)
##                    groupID+=1
##    gg.close()

if __name__=='__main__':
    main()
    print("Total time elapsed: %f" %(time.time()-t0))
    print("Maximum memory usage: %f MB" %(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000000))