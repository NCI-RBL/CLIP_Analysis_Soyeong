### random Intron sequence location

RandomIntron_chuck=function(){
  args <- commandArgs(trailingOnly = TRUE)
  nrepeat <- args[1] #100
  window <- args[2] #100
  meanBPP <- args[3] #.2
  maxBPP<- args[4] #.9
  Nlocations <- args[5] #nrow(Peaksdata4_chunk2_out) exon=138 intron=421
  UseCores <- args[6]
  folder <- args[7]
  CLIPtype <- as.character(args[8]) #Introns or #Exons
  RunMode <- as.character(args[9]) #Genomic or Transcriptomic
  TrasncLoc <- as.character(args[10]) # all, 3UTR, 5UTR, exons_CDS, introns
  
  getseq=T 
  
  # Rscript RandomExon_chunk.R 1 100 .2 .9 138 1
  
  window=as.numeric(window)
  nrepeat=as.numeric(nrepeat)
  meanBPP=as.numeric(meanBPP)
  maxBPP=as.numeric(maxBPP)
  Nlocations=as.numeric(Nlocations)
  UseCores=as.numeric(UseCores)

  #   print(window)
  # print(nrepeat)
  # print(maxBPP)
  # print(Nlocations)
  
  library(pacman)
  # p_install('GeneStructureTools')
  library(GeneStructureTools)
  library(VariantAnnotation,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library(GenomicRanges,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  #library(WhopGenome)
  #library(trackViewer)
  #library(vcfR)
  library(seqinr,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library(ggplot2,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library("viridis",quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library(edgeR,quietly = T,verbose = F)
  library('GenomicFeatures',quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library('rtracklayer',quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  #library(GeneStructureTools)
  library(matrixStats,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library(plyr,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library(tidyr,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library(fitdistrplus,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library(stringr,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
  library(data.table)
  library(reshape)
  library(GeneStructureTools)
  library(parallel)
  
  
  fasta=paste0('/data/RBL_NCI/Phil/Reference/fasta/mm10/GRCm38.primary_assembly.gencodeM23.genome.fa')
  intronbed=paste0('/data/RBL_NCI/Phil/Reference/gtf/mouse/mm10/Gencode_VM23/KnownGene/KnownGene_GRCm38_exons.bed')
  scratch=paste0('/data/RBL_NCI/Wolin/Phil/mESC_clip/structure/Random',folder,'/',CLIPtype,'/scr')
  scratch_UP=paste0('/data/RBL_NCI/Wolin/Phil/mESC_clip/structure/Random',folder,'/',CLIPtype,'/',TrasncLoc,'/',RunMode,'/upstream/scr')
  scratch_DN=paste0('/data/RBL_NCI/Wolin/Phil/mESC_clip/structure/Random',folder,'/',CLIPtype,'/',TrasncLoc,'/',RunMode,'/downstream/scr')
  outdir_UP=paste0('/data/RBL_NCI/Wolin/Phil/mESC_clip/structure/Random',folder,'/',CLIPtype,'/',TrasncLoc,'/',RunMode,'/upstream')
  outdir_DN=paste0('/data/RBL_NCI/Wolin/Phil/mESC_clip/structure/Random',folder,'/',CLIPtype,'/',TrasncLoc,'/',RunMode,'/downstream')
  qgrs=paste0('/data/RBL_NCI/Phil/Tools/qgrs-cpp/qgrs')
  fastaRegexFinder='python /data/RBL_NCI/Phil/Tools/fastaRegexFinder/fastaRegexFinder.py'
  RNAfold='RNAfold'

  
  dir.create(file.path(scratch), showWarnings = FALSE,recursive = TRUE)
  dir.create(file.path(scratch_UP), showWarnings = FALSE,recursive = TRUE)
  dir.create(file.path(scratch_DN), showWarnings = FALSE,recursive = TRUE)
  dir.create(file.path(outdir_UP), showWarnings = FALSE,recursive = TRUE)
  dir.create(file.path(outdir_DN), showWarnings = FALSE,recursive = TRUE)
  dir.create(file.path(paste0(outdir_UP,'/Table')), showWarnings = FALSE,recursive = TRUE)
  dir.create(file.path(paste0(outdir_DN,'/Table')), showWarnings = FALSE,recursive = TRUE)
  
  
  setwd(paste0('/data/RBL_NCI/Wolin/Phil/mESC_clip/structure/Random',folder,'/',CLIPtype,'/',TrasncLoc,'/',RunMode))
  
  # 
  # fasta='../../../../Resources/ref/mm10/Gencode_VM23/fromGencode/GRCm38.primary_assembly.genome.fa'
  # intronbed="/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/Gencode_VM23/fromUCSC/KnownGene/KnownGene_GRCm38_exons.bed"
  # scratch='./structure/Random/Exons/scr'
  # scratch_UP="./structure/Random/Exons/upstream/scr"
  # scratch_DN="./structure/Random/Exons/downstream/scr"
  # outdir_UP='./structure/Random/Exons/upstream'
  # outdir_DN='./structure/Random/Exons/downstream'
  # qgrs='/Users/homanpj/Documents/Tools/qgrs-cpp/qgrs'
  # fastaRegexFinder= '/Users/homanpj/Documents/Tools/fastaRegexFinder/fastaRegexFinder.py'
  # RNAfold=  '/Users/homanpj/Documents/Tools/ViennaRNA-2.4.14/bin/RNAfold'
  # 
  # window=100
  # nrepeat=10
  # meanBPP=.2
  # maxBPP=.9
  # Nlocations=138
  # UseCores=8
    
  
  print(as.numeric(window))
  print(as.numeric(nrepeat))
  print(as.numeric(meanBPP))
  print(as.numeric(maxBPP))
  print(as.numeric(Nlocations))
  print(as.numeric(UseCores))
  print(CLIPtype)
  print(RunMode)
  print(TrasncLoc)
        
  #################################################
  # Function to process RNAfold output
  #################################################
  
  RNAfold_output=function(fls){
    chunkRout3=as.data.frame(matrix(nrow=window,ncol=4))
    chunkRname3=as.data.frame(matrix(nrow=1,ncol=3))
    chunkRout3[,1]=seq(1,window,1)
    chunkRout3_max=as.data.frame(matrix(nrow=window,ncol=4))
    
    # for (x in 1:length(fls)) {
    # if (x%in%seq(from=0,to=1000,by=25)){print(x/length(fls))}
    # fls=as.list(fls)
    # fls=fls[117]
    fls=unlist(fls)
    
    chunkR=read.delim(paste0(fls[1]),comment.char = "%",header = F,sep = "\t")
    chunkR=chunkR[grep("ubox",chunkR$V1),]
    chunkR=as.data.frame(chunkR[-1],stringsAsFactors=F)
    chunkR$`chunkR[-1]`=as.character(chunkR$`chunkR[-1]`);colnames(chunkR)="V1"
    chunkR=separate(chunkR,'V1',sep = " ",into = c('b1','b2','pvalue','ubox'))
    chunkR$pvalue=as.numeric(chunkR$pvalue)
    chunkR$b1=as.numeric(chunkR$b1)
    chunkR$b2=as.numeric(chunkR$b2)
    
    
    ### count each base
    for (y in 1:(window-9)) {
      ch=chunkR[chunkR$b1%in%seq(y,y+9,1)|chunkR$b2%in%seq(y,y+9,1),]
      ch2=ch
      #   
      if (nrow(ch2)>0) {
        ch2$br1=rowMins(as.matrix(ch[,c('b1','b2')]))
        ch2$br2=rowMaxs(as.matrix(ch[,c('b1','b2')]))
        ch$b1=ch2$br1
        ch$b2=ch2$br2
        ch=ch[duplicated(ch)==F,]
        
        chunkRout3[y,1]=y
        chunkRout3[y,2]=paste0(y,'-',y+9)
        chunkRout3[y,3]=mean(ch$pvalue) ##average all basepairinng accross chunck
        #     
        chunkRout3_max[y,1]=y
        chunkRout3_max[y,2]=paste0(y,'-',y+9)
        chunkRout3[y,4]=max(ch$pvalue) ##max all basepairinng accross chunck
      }
      if (nrow(ch2)==0) {
        chunkRout3[y,1]=y
        chunkRout3[y,2]=paste0(y,'-',y+9)
        chunkRout3[y,3]=0
        chunkRout3_max[y,1]=y
        chunkRout3_max[y,2]=paste0(y,'-',y+9)
        chunkRout3[y,4]=0
      }
      chunkRname3[1,3]=gsub("_dp.ps","",fls[1])
      
    }
    
    chunkRname3=t(chunkRname3)
    colnames(chunkRout3)[3:ncol(chunkRout3)]=chunkRname3[3:nrow(chunkRname3),1]
    
    chunkRout3=chunkRout3[is.na(chunkRout3[,3])==F,]
    chunkRout3[is.nan(as.matrix(chunkRout3))==T]=0
    
    return(chunkRout3)
  } 
  ######################################################################################################################################################################################################################################################################################################################
  
  # canonical=fread("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/Gencode_VM23/fromUCSC/KnownCanonical/KnownCanonical_GencodeM23_GRCm38.txt", header=T, sep="\t",stringsAsFactors = F,data.table=F)
  canonical=fread("/data/RBL_NCI/Phil/Reference/gtf/mouse/mm10/Gencode_VM23/KnownCanonical/KnownCanonical_GRCm38.txt", header=T, sep="\t",stringsAsFactors = F,data.table=F)
  canonical$transcript=removeVersion(canonical$transcript)
  canonical$protein=removeVersion(canonical$protein)
  
  # mm10all=fread("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/Gencode_VM23/fromGencode/gencode.vM23.basic.annotation.gtf.txt", header=T, sep="\t",stringsAsFactors = F,data.table=F)
  mm10all=fread("/data/RBL_NCI/Phil/Reference/gtf/mouse/mm10/Gencode_VM23/gencode.vM23.basic.annotation.gtf.txt", header=T, sep="\t",stringsAsFactors = F,data.table=F)
    mm10all$transcript_id=removeVersion(mm10all$transcript_id)
    mm10all$gene_id=removeVersion(mm10all$gene_id)
    mm10all$start=as.numeric(as.character(mm10all$start))
    mm10all$end=as.numeric(as.character(mm10all$end))
  
  mm10_transc=mm10all[mm10all$feature%in%c("transcript"),]
  mm10_transc=mm10_transc[mm10_transc$transcript_id%in%canonical$transcript,]
  mm10_transc=mm10_transc[mm10_transc$transcript_type%in%'protein_coding',]
  
  mm10_exon=mm10all[mm10all$feature%in%c("exon"),]
  exons=mm10_exon[mm10_exon$transcript_id%in%mm10_transc$transcript_id,]
  # exons$start=as.numeric(as.character(exons$start))
  # exons$end=as.numeric(as.character(exons$end))
  # exons=as.matrix(exons)
  # exons=exons[(exons$end-exons$start)> 5,]

  mm10_UTR=mm10all[mm10all$feature%in%c("UTR"),]
  mm10_UTR=mm10_UTR[mm10_UTR$transcript_id%in%exons$transcript_id,]
  mm10_UTR=mm10_UTR[mm10_UTR$exon_id%in%exons$exon_id,]
  mm10_UTR=mm10_UTR[(mm10_UTR$end-mm10_UTR$start)> 5,]
  
  mm10_3UTR=mm10_UTR[mm10_UTR$exon_number>1,]
  mm10_5UTR=mm10_UTR[mm10_UTR$exon_number==1,]
  
  # introns=fread("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/Gencode_VM23/fromUCSC/KnownGene/KnownGene_GRCm38_introns.bed", header=F, sep="\t",stringsAsFactors = F,data.table=F)
  introns=fread(intronbed, header=F, sep="\t",stringsAsFactors = F,data.table=F)
  introns=introns[grepl("_",introns$V1)==F,]
  colnames(introns)=c('seqname','start','end','attribute','V5','strand')
  introns=separate(introns,attribute,into=c('transcript_id','feature','exon_number','level','seqname2','intronnumber','dir'),remove = T,sep = "_")
  introns$transcript_id=removeVersion(introns$transcript_id)

    introns$start=introns$start+1
  introns$exon_number=as.numeric(introns$exon_number)+1
  introns=introns[introns$transcript_id%in%exons$transcript_id,]
  
  
  
  outchunck_rowaverage=matrix(nrow=nrepeat,ncol = (window-9) )
  outchunck_rowaverage_max=matrix(nrow=nrepeat,ncol = (window-9) )
  # outchunck_rowaverage=matrix(nrow=nrepeat,ncol = (window) )
  outchunck_frac=matrix(nrow=nrepeat,ncol = 2 )
  outchunck_frac_max=matrix(nrow=nrepeat,ncol = 2 )
  outchunck_rowaverage_dn=matrix(nrow=nrepeat,ncol = (window-9) )
  outchunck_rowaverage_maxdn=matrix(nrow=nrepeat,ncol = (window-9) )
  # outchunck_rowaverage_dn=matrix(nrow=nrepeat,ncol = (window) )
  outchunck_frac_dn=matrix(nrow=nrepeat,ncol = 2 )
  outchunck_frac_maxdn=matrix(nrow=nrepeat,ncol = 2 )
  
  
  n_UP=matrix(nrow=Nlocations,ncol=nrepeat)
  ID_UP=matrix(nrow=Nlocations,ncol=nrepeat)
  seq_UP=matrix(nrow=Nlocations,ncol=nrepeat)
  SingleExon_UP=matrix(nrow=Nlocations,ncol=nrepeat)
  Exon_UP=matrix(nrow=Nlocations,ncol=nrepeat)
  transcript_UP=matrix(nrow=Nlocations,ncol=nrepeat)
  LastExon_UP=matrix(nrow=Nlocations,ncol=nrepeat)
  Gscore_UP=matrix(nrow=Nlocations,ncol=nrepeat)
  GscoreDist_UP=matrix(nrow=Nlocations,ncol=nrepeat)
  GscoreMaxDist_UP=matrix(nrow=Nlocations,ncol=nrepeat)
  GscoreMax_UP=matrix(nrow=Nlocations,ncol=nrepeat)
  GscoreMax2_UP=matrix(nrow=Nlocations,ncol=nrepeat)
  FinderSeq_UP=matrix(nrow=Nlocations,ncol=nrepeat)
  qgrsSeq_UP=matrix(nrow=Nlocations,ncol=nrepeat)
  maxBPPmean_UP=matrix(nrow=Nlocations,ncol=nrepeat)
  maxBPPmax_UP=matrix(nrow=Nlocations,ncol=nrepeat)
  MFE_UP=matrix(nrow=Nlocations,ncol=nrepeat)
  
  n_DN=matrix(nrow=Nlocations,ncol=nrepeat)
  ID_DN=matrix(nrow=Nlocations,ncol=nrepeat)
  seq_DN=matrix(nrow=Nlocations,ncol=nrepeat)
  SingleExon_DN=matrix(nrow=Nlocations,ncol=nrepeat)
  Exon_DN=matrix(nrow=Nlocations,ncol=nrepeat)
  transcript_DN=matrix(nrow=Nlocations,ncol=nrepeat)
  LastExon_DN=matrix(nrow=Nlocations,ncol=nrepeat)
  Gscore_DN=matrix(nrow=Nlocations,ncol=nrepeat)
  GscoreDist_DN=matrix(nrow=Nlocations,ncol=nrepeat)
  GscoreMaxDist_DN=matrix(nrow=Nlocations,ncol=nrepeat)
  GscoreMax_DN=matrix(nrow=Nlocations,ncol=nrepeat)
  GscoreMax2_DN=matrix(nrow=Nlocations,ncol=nrepeat)
  FinderSeq_DN=matrix(nrow=Nlocations,ncol=nrepeat)
  qgrsSeq_DN=matrix(nrow=Nlocations,ncol=nrepeat)
  maxBPPmean_DN=matrix(nrow=Nlocations,ncol=nrepeat)
  maxBPPmax_DN=matrix(nrow=Nlocations,ncol=nrepeat)
  MFE_DN=matrix(nrow=Nlocations,ncol=nrepeat)
  
  checkoutR_summary=matrix(nrow=nrepeat,ncol=12)
  colnames(checkoutR_summary)=c('SingleExon','numberExCrossed','lastexon','lastexonCrossexExon','firstexon','firstexonCrossexExon','total Locations','Gquad','Avg_GquadScore','MeanBPP','MaxBPP','UTR')
  checkoutR_dn_summary=matrix(nrow=nrepeat,ncol=12)
  colnames(checkoutR_dn_summary)=c('SingleExon','numberExCrossed','lastexon','lastexonCrossexExon','firstexon','firstexonCrossexExon',"total Locations",'Gquad','Avg_GquadScore','MeanBPP','MaxBPP','UTR')
  
  # CLIPtype='Exons' #Introns or #Exons
  # RunMode='Transcriptomic' #Genomic or Transcriptomic
  # TrasncLoc='3UTR' # all, 3UTR, 5UTR, exons_CDS, introns
  # 
  if (CLIPtype=='Exons') {mm10_USE=exons }
  if (CLIPtype=='Introns') {mm10_USE=introns }

  system.time({
    for (r in 1:nrepeat) {
      print(r)
      
      checkoutR=as.data.frame(matrix(nrow = Nlocations,ncol=16));colnames(checkoutR)=c('n','ID','crosslink','SingleExon','numberExCrossed','crosslinkexon','s1exon','e1exon','lastexon','TranscChoise','numberTranscript','trancript','exon#','sequence','strand','DeltaG')
      checkoutR_dn=as.data.frame(matrix(nrow = Nlocations,ncol=16));colnames(checkoutR_dn)=c('n','ID','crosslink','SingleExon','numberExCrossed','crosslinkexon','s1dnexon','e1dnexon','lastexon','TranscChoise','numberTranscript','trancript','exon#','sequence','strand','DeltaG')
      faout=matrix(nrow = (Nlocations*2),ncol=1)
      faout_dn=matrix(nrow = (Nlocations*2),ncol=1)
      
      ##### select exons
      if (TrasncLoc!='3UTR') {
        mm10_USEselect=mm10_USE[(mm10_USE$end-mm10_USE$start)> 5,]
        rand=sample(1:nrow(mm10_USEselect), Nlocations, replace=F)
        exons_rand=mm10_USEselect[rand,] }
      
      if (TrasncLoc=='3UTR') {
        rand=sample(1:nrow(mm10_3UTR), Nlocations, replace=F)
        UTR_rand=mm10_3UTR[rand,] }
      
      for (x in 1:Nlocations) {
            ##### select location in exon
        if (TrasncLoc!='3UTR') {texn=exons_rand[x,]}
        if (TrasncLoc=='3UTR') {
          tutr=UTR_rand[x,]
          texn=exons[exons$exon_id%in%tutr$exon_id,]
          if (nrow(texn)==0) {noexon}
        }
        
        if (length(unique(texn$transcript_id))>1) {xxxxxxxxcxcxcx}
        
        strand=texn$strand
        chr=texn$seqname
        
        if (TrasncLoc!='3UTR') {texn$CL=sample(c((texn$start+1):(texn$end-1)), 1, replace=F)}
        if (TrasncLoc=='3UTR') {texn$CL=sample(c((tutr$start+1):(tutr$end-1)), 1, replace=F)}#; print(x); print(paste0(texn$seqname,':',texn$CL,'_',texn$strand))}
        
        if (strand=="+") {
          texn$s1=texn$CL-window#(selects base +1 of start site but ends at correct 5' site)
          texn$e1=texn$CL
          
          texn$s1dn=texn$CL-1
          texn$e1dn=texn$CL+window-1
        }
        if (strand=="-") {
          texn$s1=texn$CL-1#(selects base +1 of start site but ends at correct 5' site)
          texn$e1=texn$CL+window-1
          
          texn$s1dn=texn$CL-window
          texn$e1dn=texn$CL
        }
        
        s1=texn$s1
        e1=texn$e1
        
        s1dn=texn$s1dn
        e1dn=texn$e1dn
        
        checkoutR[x,'trancript']=unique(texn$transcript_id)
        checkoutR[x,'exon#']=(texn$exon_number)
        checkoutR[x,'strand']=strand
        checkoutR[x,'ID']=paste0(chr,":",texn$CL,"_",texn$strand)
        checkoutR[x,'crosslink']=texn$CL
        
        checkoutR_dn[x,'trancript']=unique(texn$transcript_id)
        checkoutR_dn[x,'exon#']=(texn$exon_number)
        checkoutR_dn[x,'strand']=strand
        checkoutR_dn[x,'ID']=paste0(chr,":",texn$CL,"_",texn$strand)
        checkoutR_dn[x,'crosslink']=texn$CL
        
        ### UPSTREAM
        #############################################################################################
        checkoutR[x,'SingleExon']=(texn$start<=texn$s1 & texn$end>=texn$e1)|
          (strand=="+"&(texn$start>=texn$s1&texn$end>=texn$e1)&texn$exon_number==1)| 
          (strand=="-"&(texn$start<=texn$s1&texn$end<=texn$e1)&texn$exon_number==1)
        if ((strand=="+"&(texn$start>=texn$s1&texn$end>=texn$e1)&texn$exon_number==1)|
            (strand=="-"&(texn$start<=texn$s1&texn$end<=texn$e1)&texn$exon_number==1)) {checkoutR[x,'numberTranscript']=paste0(checkoutR[x,'numberTranscript'],"_out of bounds exon upstream only")}
        
        
        if (RunMode=='Genomic') {checkoutR[x,'SingleExon']=T}
        
        
        
        #### window region falls in exon
        checkoutR[x,'crosslinkexon']=texn$start<=texn$CL&texn$end>=texn$CL #crosslink
        checkoutR[x,'s1exon']=(strand=="+"&texn$start<=s1&texn$end>=s1)|(strand=="-"&texn$start<=e1&texn$end>=e1) # 5` end in exon
        checkoutR[x,'e1exon']=(strand=="+"&texn$start<=e1&texn$end>=e1)|(strand=="-"&texn$start<=s1&texn$end>=s1) # 3` end in exon
        checkoutR[x,'lastexon']=checkoutR[x,'exon#']/max(mm10_USE[mm10_USE$transcript_id%in%checkoutR[x,'trancript'],'exon_number'])
        checkoutR[x,'crosslinkexon']=texn$start<=texn$CL&texn$end>=texn$CL #crosslink
        
        
      ###########################################################
        #run if not in single exon and not exon 1 (if crosslink is in exon then upstream has to be upstream out of gene for exon1)
        if ((checkoutR[x,'SingleExon']==F)&(checkoutR[x,'exon#']>1)) {
          bedpos=matrix(ncol=6,nrow=2)
          
          locx=mm10_USE[mm10_USE$transcript_id%in%texn$transcript_id,]
          # locx=locx[order(locx$exon_number,decreasing = F),]
          if (strand=="+") {
            if (checkoutR[x,'exon#',drop=F]==1) {next}
            
            # closest to peak
            bedpos[2,2]=locx[(locx$exon_number %in%texn$exon_number),'start']-1
            bedpos[2,3]=texn$CL
            
            # Upstream exon
            bedpos[1,2]=locx[(locx$exon_number %in%(texn$exon_number-1)),'end']-(window-abs(texn$CL-(locx[(locx$exon_number %in%texn$exon_number),'start'])+1))
            bedpos[1,3]=locx[(locx$exon_number %in%(texn$exon_number-1)),'end']
            
            checkoutR[x,'numberExCrossed']=nrow(bedpos)
            
            # exon number of bedpos[1,]
            exnum=(texn$exon_number)-(nrow(bedpos)-1)
            while ((exnum>1)&(locx[(locx$exon_number %in%(exnum)),'start']>bedpos[1,2])) {#xnxnxnx

              bedpos[1,2]=locx[(locx$exon_number %in%(exnum)),'start']-1
              bedpos[1,3]=locx[(locx$exon_number %in%(exnum)),'end']
              
              bedpos2=matrix(nrow=1,ncol = 6)
              bedpos2[1,2]=locx[(locx$exon_number %in%(exnum-1)),'end']-((window-1)-sum(rowDiffs(bedpos[,c(2,3)]))+1)
              bedpos2[1,3]=locx[(locx$exon_number %in%((exnum-1))),'end']
              
              bedpos=rbind(bedpos2,bedpos)
              checkoutR[x,'numberExCrossed']=nrow(bedpos)
              exnum=(texn$exon_number)-(nrow(bedpos)-1)-1
              if (exnum==0) {break}
            }
          }
          
          if (strand=="-") {
            if (checkoutR[x,'exon#']==1) {next}
            
            # closest to peak
            bedpos[2,2]=texn$CL-1
            bedpos[2,3]=locx[(locx$exon_number %in%texn$exon_number),'end']
            
            # Upstream exon
            bedpos[1,2]=locx[(locx$exon_number %in%(texn$exon_number-1)),'start']-1
            bedpos[1,3]=locx[(locx$exon_number %in%(texn$exon_number-1)),'start']+((window-1)-abs(texn$CL-locx[(locx$exon_number %in%(texn$exon_number)),'end']))-1
            
            checkoutR[x,'numberExCrossed']=nrow(bedpos)
            
            
                  exnum=(texn$exon_number)-(nrow(bedpos)-1)
            while ((exnum>1)&(locx[(locx$exon_number %in%(exnum)),'end']<bedpos[1,3])) {#xnegxneg
              
              bedpos[1,2]=locx[(locx$exon_number %in%(exnum)),'start']-1
              bedpos[1,3]=locx[(locx$exon_number %in%(exnum)),'end']
              bedpos2=matrix(nrow=1,ncol = 6)
              bedpos2[1,2]=locx[(locx$exon_number %in%((exnum-1))),'start']-1
              bedpos2[1,3]=locx[(locx$exon_number %in%((exnum-1))),'start']+((window)-sum(rowDiffs(bedpos[,c(2,3)]))-1)

              bedpos=rbind(bedpos2,bedpos)
              checkoutR[x,'numberExCrossed']=nrow(bedpos)
              exnum=(texn$exon_number)-(nrow(bedpos)-1)-1
              if (exnum==0) {break}
            }
          }
          
          bedposx=as.data.frame(bedpos)
          bedposx[,2]=as.numeric(as.character(bedposx[,2]))
          bedposx[,3]=as.numeric(as.character(bedposx[,3]))
          sum(rowDiffs(as.matrix(bedposx[,c(2,3)])) )       
          (rowDiffs(as.matrix(bedposx[,c(2,3)])) )       
          
          bedpos[,1]=chr
          bedpos[,4]=".";bedpos[,5]="."
          bedpos[,6]=strand
          
          write.table(bedpos, file=paste0(scratch_UP,'/',texn$seqname,'_',texn$CL,'_chunk.bed'), quote=F, sep="\t", row.names=F, col.names=F)
          # system(paste0('bedtools getfasta -s -fi ',fasta,' -bed structure/chunkwindow/exon/',texn$ID,'_chunk.bed > structure/chunkwindow/exon/',texn$ID,'_chunk.fa'))
          fa=system(paste0('bedtools getfasta -s -fi ',fasta,' -bed ',scratch_UP,'/',texn$seqname,'_',texn$CL,'_chunk.bed'),intern = T)
          
          
          if (nchar(fa[2])==window) {fa2=fa;checkoutR[x,'SingleExon']='Change-TRUE'}
          
          fa2=fa
          for (f in 2:nrow(bedpos)) {
          fa2[2]=paste(fa2[2],fa2[f*2],collapse = "")
          fa2[2]=gsub(" ","",fa2[2])
          
          if ((nchar(fa2[2])==window)&f<nrow(bedpos)){
            # checkoutR=TRUE
            checkoutR[x,'TranscChoise']='Change Multiple features'
            fa2=fa2[1:2]
            checkoutR[x,'numberExCrossed']=checkoutR[x,'numberExCrossed']-1
            break}
          }
          
          
          fa2=fa2[1:2]
          if (nchar(fa[2])==window){
            checkoutR[x,'SingleExon']=TRUE
            checkoutR[x,'TranscChoise']='Change'
            checkoutR[x,'numberExCrossed']=checkoutR[x,'numberExCrossed']-1
            fa2=fa[1:2]}
          
          if (nchar(fa2[2])!=window){print(nchar(fa2[2]));shrtstrand}
          fa=fa2
          write.table(fa2, file=paste0(scratch_UP,'/',texn$seqname,'_',texn$CL,'_chunk.fa'), quote=F, sep="\t", row.names=F, col.names=F)
          checkoutR[x,'sequence']=gsub('T','U',fa[2])
          checkoutR[x,'n']=fa[1]
          # faout=c(faout,fa[1:2])
                faout[(2*x-1):(2*x),1]=(as.matrix(fa[1:2]))
          remove(bedpos,fa2)
        }
        
        if ((checkoutR[x,'SingleExon']==T)|(checkoutR[x,'exon#']==1)) {
          bedpos=matrix(ncol=6,nrow=1)
          bedpos[1,1]=chr
          bedpos[1,2]=s1
          bedpos[1,3]=e1
          bedpos[1,4]=".";bedpos[,5]="."
          bedpos[1,6]=strand
          
          write.table(bedpos, file=paste0(scratch_UP,'/',texn$seqname,'_',texn$CL,'_chunk.bed'), quote=F, sep="\t", row.names=F, col.names=F)
          # system(paste0('bedtools getfasta -s -fi ',fasta,' -bed structure/chunkwindow/exon/',texn$ID,'_chunk.bed > structure/chunkwindow/exon/',texn$ID,'_chunk.fa'))
          fa=system(paste0('bedtools getfasta -s -fi ',fasta,' -bed ',scratch_UP,'/',texn$seqname,'_',texn$CL,'_chunk.bed'),intern = T)
          write.table(fa, file=paste0(scratch_UP,'/',texn$seqname,'_',texn$CL,'_chunk.fa'), quote=F, sep="\t", row.names=F, col.names=F)
          checkoutR[x,'sequence']=gsub('T','U',fa[2])
          checkoutR[x,'n']=fa[1]
          checkoutR[x,'numberExCrossed']=0
          
          # faout=c(faout,fa[1:2])
                faout[(2*x-1):(2*x),1]=(as.matrix(fa[1:2]))
          
          remove(bedpos)
        }
      
        if(file.exists(paste0(scratch_UP,'/',texn$seqname,'_',texn$CL,'_chunk.fa'))){file.remove(paste0(scratch_UP,'/',texn$seqname,'_',texn$CL,'_chunk.fa'))}
        if(file.exists(paste0(scratch_UP,'/',texn$seqname,'_',texn$CL,'_chunk.bed'))){file.remove(paste0(scratch_UP,'/',texn$seqname,'_',texn$CL,'_chunk.bed'))}
  
        ## DOWNSTREAM
      #################################################################################################################################################################################
      #################################################################################################################################################################################
        checkoutR_dn[x,'lastexon']=checkoutR_dn[x,'exon#']/max(mm10_USE[mm10_USE$transcript_id%in%checkoutR_dn[x,'trancript'],'exon_number'])
        checkoutR_dn[x,'crosslinkexon']=texn$start<=texn$CL&texn$end>=texn$CL #crosslink
        
        
        checkoutR_dn[x,'SingleExon']=(texn$start<=texn$s1dn & texn$end>=texn$e1dn)|
        (strand=="+"&(texn$start<=texn$s1dn&texn$end<=texn$e1dn)&checkoutR_dn[x,'lastexon']==1)| 
        (strand=="-"&(texn$start>=texn$s1dn&texn$end>=texn$e1dn)&checkoutR_dn[x,'lastexon']==1)
      
        
        if (RunMode=='Genomic') {checkoutR_dn[x,'SingleExon']=T}
        
        
        if ( ###start is out of bounds but window end is in
          checkoutR_dn[x,'SingleExon']==F&
          ((strand=="+"&texn$start>=texn$s1dn&texn$end>=texn$e1dn)|
           (strand=="-"&texn$end<=texn$e1dn&texn$start<=texn$s1dn))
        ) {
          checkoutR_dn[x,'SingleExon']=T
          checkoutR_dn[x,'numberTranscript']=paste0(checkoutR_dn[x,'numberTranscript'],"5` end is out of bounds but single exon")}
        #   
      
      #### window region faDNlls in exon
      checkoutR_dn[x,'s1dnexon']=(strand=="+"&texn$start<=s1dn&texn$end>=s1dn)|(strand=="-"&texn$start<=e1dn&texn$end>=e1dn) # 5` end in exon
      checkoutR_dn[x,'e1dnexon']=(strand=="+"&texn$start<=e1dn&texn$end>=e1dn)|(strand=="-"&texn$start<=s1dn&texn$end>=s1dn) # 3` end in exon
      
      ###########################################################
      #run if not in single exon and not exon 1 (if crosslink is in exon then upstream has to be upstream out of gene for exon1)
      if ((checkoutR_dn[x,'SingleExon']==F)&(checkoutR_dn[x,'lastexon']<1)) {
        bedposDN=matrix(ncol=6,nrow=2)
        
        locx=mm10_USE[mm10_USE$transcript_id%in%texn$transcript_id,]
        if (strand=="+") {
          if (checkoutR_dn[x,'lastexon',drop=F]==1) {next}
          
          # closest to peak
          bedposDN[2,2]=texn$CL-1
          bedposDN[2,3]=locx[(locx$exon_number %in%texn$exon_number),'end']
          
          
          # Upstream exon
          bedposDN[1,2]=locx[(locx$exon_number %in%(texn$exon_number+1)),'start']-1
          bedposDN[1,3]=locx[(locx$exon_number %in%(texn$exon_number+1)),'start']+(window-1-abs(texn$CL-(locx[(locx$exon_number %in%texn$exon_number),'end'])))-1
          
          checkoutR_dn[x,'numberExCrossed']=nrow(bedposDN)
          
          # exon number of bedposDN[1,]
          exnum=(texn$exon_number)+(nrow(bedposDN)-1)
          while ((exnum<max(locx$exon_number))&(locx[(locx$exon_number %in%(exnum)),'end']<bedposDN[1,3])) {#xnxnxnx
            
            bedposDN[1,2]=locx[(locx$exon_number %in%(exnum)),'start']-1
            bedposDN[1,3]=locx[(locx$exon_number %in%(exnum)),'end']
            
            bedposDN2=matrix(nrow=1,ncol = 6)
            bedposDN2[1,2]=locx[(locx$exon_number %in%((exnum+1))),'start']-1
            bedposDN2[1,3]=locx[(locx$exon_number %in%((exnum+1))),'start']+((window-1)-sum(rowDiffs(bedposDN[,c(2,3)])))
            
            bedposDN=rbind(bedposDN2,bedposDN)
            checkoutR_dn[x,'numberExCrossed']=nrow(bedposDN)
            exnum=exnum+1
            if (exnum==max(locx$exon_number)) {break}
          }
        }
        
        if (strand=="-") {
          if (checkoutR_dn[x,'lastexon',drop=F]==1) {next}
          
          # closest to peak
          bedposDN[2,2]=locx[(locx$exon_number %in%texn$exon_number),'start']-1
          bedposDN[2,3]=texn$CL
          
          # downstream exon
          bedposDN[1,2]=locx[(locx$exon_number %in%(texn$exon_number+1)),'end']-((window-1)-abs(texn$CL-locx[(locx$exon_number %in%(texn$exon_number)),'start']))
          bedposDN[1,3]=locx[(locx$exon_number %in%(texn$exon_number+1)),'end']
          
          checkoutR_dn[x,'numberExCrossed']=nrow(bedposDN)
          
          
          exnum=(texn$exon_number)+(nrow(bedposDN)-1)
          while ((exnum<max(locx$exon_number))&(locx[(locx$exon_number %in%(exnum)),'start']>bedposDN[1,2])) {#xnegxneg
            
            bedposDN[1,2]=locx[(locx$exon_number %in%(exnum)),'start']-1
            bedposDN[1,3]=locx[(locx$exon_number %in%(exnum)),'end']
            bedposDN2=matrix(nrow=1,ncol = 6)
            bedposDN2[1,2]=locx[(locx$exon_number %in%((exnum+1))),'end']-((window)-sum(rowDiffs(bedposDN[,c(2,3)])))
            bedposDN2[1,3]=locx[(locx$exon_number %in%((exnum+1))),'end']
            
            bedposDN=rbind(bedposDN2,bedposDN)
            checkoutR_dn[x,'numberExCrossed']=nrow(bedposDN)
            exnum=exnum+1
            if (exnum==max(locx$exon_number)) {break}
          }
        }
        
        bedposDNx=as.data.frame(bedposDN)
        bedposDNx[,2]=as.numeric(as.character(bedposDNx[,2]))
        bedposDNx[,3]=as.numeric(as.character(bedposDNx[,3]))
        sum(rowDiffs(as.matrix(bedposDNx[,c(2,3)])) )       
        (rowDiffs(as.matrix(bedposDNx[,c(2,3)])) )       
        
        bedposDN[,1]=chr
        bedposDN[,4]=".";bedposDN[,5]="."
        bedposDN[,6]=strand
        
        write.table(bedposDN, file=paste0(scratch_DN,'/',texn$seqname,'_',texn$CL,'_chunkDN.bed'), quote=F, sep="\t", row.names=F, col.names=F)
        # system(paste0('bedtools getfasta -s -fi ',fasta,' -bed structure/chunkwindow/exon/',texn$ID,'_chunkDN.bed > structure/chunkwindow/exon/',texn$ID,'_chunk.faDN'))
        faDN=system(paste0('bedtools getfasta -s -fi ',fasta,' -bed ',scratch_DN,'/',texn$seqname,'_',texn$CL,'_chunkDN.bed'),intern = T)
        if (nchar(faDN[2])==window) {
          faDN2=faDN[1:2];faDN=faDN[1:2];
          checkoutR_dn[x,'SingleExon']='Change-TRUE'
          checkoutR_dn[x,'numberExCrossed']=checkoutR_dn[x,'numberExCrossed']-1
        }
        
        if (length(faDN)>2) {
        faDN2=faDN[c(length(faDN):1)]
        for (f in seq(3,length(faDN2),2)) {
          faDN2[1]=paste(faDN2[1],faDN2[f],collapse = "")
          faDN2[1]=gsub(" ","",faDN2[1])
          
          if ((nchar(faDN2[1])==window)&f<max(seq(3,length(faDN2),2)) ){
            # checkoutR_dn=TRUE
            checkoutR_dn[x,'TranscChoise']='Change Multiple freatures'
            faDN2=faDN2[1:2]
            checkoutR_dn[x,'numberExCrossed']=checkoutR_dn[x,'numberExCrossed']-1
            break}
        }
        faDN2=faDN2[2:1]
        }
        
        if (nchar(faDN2[2])!=window){print(nchar(faDN2[2]));shrtstrand}
         
        write.table(faDN2, file=paste0(scratch_DN,'/',texn$seqname,'_',texn$CL,'_chunkDN.fa'), quote=F, sep="\t", row.names=F, col.names=F)
        checkoutR_dn[x,'sequence']=gsub('T','U',faDN2[2])
        checkoutR_dn[x,'n']=faDN2[1]
        
        # faout_dn=c(faout_dn,faDN2[1:2])
            faout_dn[(2*x-1):(2*x),1]=(as.matrix(faDN2[1:2]))
        
        remove(bedposDN,faDN2)
      }
 
      if ((checkoutR_dn[x,'SingleExon']==T)|(checkoutR_dn[x,'lastexon']==1)) {
        bedposDN=matrix(ncol=6,nrow=1)
        bedposDN[1,1]=chr
        bedposDN[1,2]=s1dn
        bedposDN[1,3]=e1dn
        bedposDN[1,4]=".";bedposDN[,5]="."
        bedposDN[1,6]=strand
        
        write.table(bedposDN, file=paste0(scratch_DN,'/',texn$seqname,'_',texn$CL,'_chunkDN.bed'), quote=F, sep="\t", row.names=F, col.names=F)
        # system(paste0('bedtools getfasta -s -fi ',fasta,' -bed structure/chunkwindow/exon/',texn$ID,'_chunkDN.bed > structure/chunkwindow/exon/',texn$ID,'_chunk.faDN'))
        faDN=system(paste0('bedtools getfasta -s -fi ',fasta,' -bed ',scratch_DN,'/',texn$seqname,'_',texn$CL,'_chunkDN.bed'),intern = T)
        write.table(faDN, file=paste0(scratch_DN,'/',texn$seqname,'_',texn$CL,'_chunkDN.fa'), quote=F, sep="\t", row.names=F, col.names=F)
        checkoutR_dn[x,'sequence']=gsub('T','U',faDN[2])
        checkoutR_dn[x,'n']=faDN[1]
        checkoutR_dn[x,'numberExCrossed']=0
        # faout_dn=c(faout_dn,faDN[1:2])
            faout_dn[(2*x-1):(2*x),1]=(as.matrix(faDN[1:2]))
        
        remove(bedposDN)
      }
      
      if(file.exists(paste0(scratch_DN,'/',texn$seqname,'_',texn$CL,'_chunkDN.fa'))){file.remove(paste0(scratch_DN,'/',texn$seqname,'_',texn$CL,'_chunkDN.fa'))}
      if(file.exists(paste0(scratch_DN,'/',texn$seqname,'_',texn$CL,'_chunkDN.bed'))){file.remove(paste0(scratch_DN,'/',texn$seqname,'_',texn$CL,'_chunkDN.bed'))}
     
    }
    
      checkoutR$length=apply(checkoutR[,'sequence',drop=F], 1, nchar)
      checkoutR$n=gsub(">","",checkoutR$n)
      
    checkoutR_dn$length=apply(checkoutR_dn[,'sequence',drop=F], 1, nchar)
    checkoutR_dn$n=gsub(">","",checkoutR_dn$n)
    
  
  ##########################################################################################  
  ##### Create Fasta file for all locations
  ##########################################################################################  
    
    faout=faout[is.na(faout)==F] 
    faout=faout[duplicated(faout)==F]
      faout=gsub('T','U',faout)
      faout=gsub('\\(\\+)',"",faout)
      faout=gsub('\\(\\-)',"",faout)
      
    faout_dn=faout_dn[is.na(faout_dn)==F]
    faout_dn=faout_dn[duplicated(faout_dn)==F]
      faout_dn=gsub('T','U',faout_dn)
      faout_dn=gsub('\\(\\+)',"",faout_dn)
      faout_dn=gsub('\\(\\-)',"",faout_dn)
      
    write.table(faout, file=paste0(outdir_UP, '/AllLocation_Random.fa'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(faout_dn, file=paste0(outdir_DN, '/AllLocation_DN_Random.fa'), quote=F, sep="\t", row.names=F, col.names=F)
    
    checkoutR$n=gsub('\\(\\+)',"",checkoutR$n)
    checkoutR$n=gsub('\\(\\-)',"",checkoutR$n) 
    
    checkoutR_dn$n=gsub('\\(\\+)',"",checkoutR_dn$n)
    checkoutR_dn$n=gsub('\\(\\-)',"",checkoutR_dn$n)
    
    
if (getseq==T) {next}
  

  ##########################################################################################  
  ##### Count Nucleotide Content
  ##########################################################################################  
   
  # upstream   
    checkoutR=checkoutR[is.na(checkoutR$sequence)==F,]
    checkoutR$sequence=gsub('T','U',checkoutR$sequence)
    for (ch in c(10,20,30,40,50)) { chadj=ch-1
    
    chucnkseqA=as.data.frame(matrix(nrow=nrow(checkoutR),ncol=(nchar(checkoutR[1,'sequence'])-chadj)))
    rownames(chucnkseqA)=checkoutR$ID;colnames(chucnkseqA)=seq(-(nchar(checkoutR[1,'sequence'])-chadj),-1)
    chucnkseqU=as.data.frame(matrix(nrow=nrow(checkoutR),ncol=(nchar(checkoutR[1,'sequence'])-chadj)))
    rownames(chucnkseqU)=checkoutR$ID;colnames(chucnkseqU)=seq(-(nchar(checkoutR[1,'sequence'])-chadj),-1)
    chucnkseqG=as.data.frame(matrix(nrow=nrow(checkoutR),ncol=(nchar(checkoutR[1,'sequence'])-chadj)))
    rownames(chucnkseqG)=checkoutR$ID;colnames(chucnkseqG)=seq(-(nchar(checkoutR[1,'sequence'])-chadj),-1)
    chucnkseqC=as.data.frame(matrix(nrow=nrow(checkoutR),ncol=(nchar(checkoutR[1,'sequence'])-chadj)))
    rownames(chucnkseqC)=checkoutR$ID;colnames(chucnkseqC)=seq(-(nchar(checkoutR[1,'sequence'])-chadj),-1)
    
    for (l in 1:(nchar(checkoutR[1,'sequence'])-chadj)) {
      ss= substr(checkoutR[,'sequence'], start = l, stop = (l+chadj))
      chucnkseqA[,l]=str_count(ss, "A") /ch
      chucnkseqU[,l]=str_count(ss, "U") /ch
      chucnkseqG[,l]=str_count(ss, "G") /ch
      chucnkseqC[,l]=str_count(ss, "C") /ch
    }
    
    if (r==1) {
    chucnkseqAout=as.data.frame(matrix(nrow=nrepeat,ncol=(nchar(checkoutR[1,'sequence'])-chadj)))
      colnames(chucnkseqAout)=seq(-(nchar(checkoutR[1,'sequence'])-chadj),-1)
    chucnkseqUout=as.data.frame(matrix(nrow=nrepeat,ncol=(nchar(checkoutR[1,'sequence'])-chadj)))
      colnames(chucnkseqUout)=seq(-(nchar(checkoutR[1,'sequence'])-chadj),-1)
    chucnkseqGout=as.data.frame(matrix(nrow=nrepeat,ncol=(nchar(checkoutR[1,'sequence'])-chadj)))
      colnames(chucnkseqGout)=seq(-(nchar(checkoutR[1,'sequence'])-chadj),-1)
    chucnkseqCout=as.data.frame(matrix(nrow=nrepeat,ncol=(nchar(checkoutR[1,'sequence'])-chadj)))
      colnames(chucnkseqCout)=seq(-(nchar(checkoutR[1,'sequence'])-chadj),-1)
    }  
   
    if (r>1) {  
    chucnkseqAout=fread(file=paste0(outdir_UP,'/chucnkseqA_window',ch,'.txt'),header = T,sep = "\t",stringsAsFactors = F,data.table=F)
    chucnkseqUout=fread(file=paste0(outdir_UP,'/chucnkseqU_window',ch,'.txt'),header = T,sep = "\t",stringsAsFactors = F,data.table=F)
    chucnkseqGout=fread(file=paste0(outdir_UP,'/chucnkseqG_window',ch,'.txt'),header = T,sep = "\t",stringsAsFactors = F,data.table=F)  
    chucnkseqCout=fread(file=paste0(outdir_UP,'/chucnkseqC_window',ch,'.txt'),header = T,sep = "\t",stringsAsFactors = F,data.table=F)
    }
    chucnkseqAoutadd=rbind(chucnkseqAout,colMeans(chucnkseqA))
    chucnkseqUoutadd=rbind(chucnkseqUout,colMeans(chucnkseqU))
    chucnkseqGoutadd=rbind(chucnkseqGout,colMeans(chucnkseqG))
    chucnkseqCoutadd=rbind(chucnkseqCout,colMeans(chucnkseqC))
      
    chucnkseqAoutadd=chucnkseqAoutadd[is.na(chucnkseqAoutadd[,1])==F,]
    chucnkseqUoutadd=chucnkseqUoutadd[is.na(chucnkseqUoutadd[,1])==F,]
    chucnkseqGoutadd=chucnkseqGoutadd[is.na(chucnkseqGoutadd[,1])==F,]
    chucnkseqCoutadd=chucnkseqCoutadd[is.na(chucnkseqCoutadd[,1])==F,]
    
      write.table(chucnkseqAoutadd, file=paste0(outdir_UP,'/chucnkseqA_window',ch,'.txt'), quote=F, sep="\t", row.names=F, col.names=T)
      write.table(chucnkseqUoutadd, file=paste0(outdir_UP,'/chucnkseqU_window',ch,'.txt'), quote=F, sep="\t", row.names=F, col.names=T)
      write.table(chucnkseqGoutadd, file=paste0(outdir_UP,'/chucnkseqG_window',ch,'.txt'), quote=F, sep="\t", row.names=F, col.names=T)
      write.table(chucnkseqCoutadd, file=paste0(outdir_UP,'/chucnkseqC_window',ch,'.txt'), quote=F, sep="\t", row.names=F, col.names=T)
    
      remove(chucnkseqAoutadd,chucnkseqUoutadd,chucnkseqGoutadd,chucnkseqCoutadd,chucnkseqAout,chucnkseqUout,chucnkseqGout,chucnkseqCout)  
    }

    
    # downstream   
    checkoutR_dn=checkoutR_dn[is.na(checkoutR_dn$sequence)==F,]
    checkoutR_dn$sequence=gsub('T','U',checkoutR_dn$sequence)
    for (ch in c(10,20,30,40,50)) { chadj=ch-1
    
    chucnkseqA=as.data.frame(matrix(nrow=nrow(checkoutR_dn),ncol=(nchar(checkoutR_dn[1,'sequence'])-chadj)))
    rownames(chucnkseqA)=checkoutR_dn$ID;colnames(chucnkseqA)=seq(1,(nchar(checkoutR_dn[1,'sequence'])-chadj))
    chucnkseqU=as.data.frame(matrix(nrow=nrow(checkoutR_dn),ncol=(nchar(checkoutR_dn[1,'sequence'])-chadj)))
    rownames(chucnkseqU)=checkoutR_dn$ID;colnames(chucnkseqU)=seq(1,(nchar(checkoutR_dn[1,'sequence'])-chadj))
    chucnkseqG=as.data.frame(matrix(nrow=nrow(checkoutR_dn),ncol=(nchar(checkoutR_dn[1,'sequence'])-chadj)))
    rownames(chucnkseqG)=checkoutR_dn$ID;colnames(chucnkseqG)=seq(1,(nchar(checkoutR_dn[1,'sequence'])-chadj))
    chucnkseqC=as.data.frame(matrix(nrow=nrow(checkoutR_dn),ncol=(nchar(checkoutR_dn[1,'sequence'])-chadj)))
    rownames(chucnkseqC)=checkoutR_dn$ID;colnames(chucnkseqC)=seq(1,(nchar(checkoutR_dn[1,'sequence'])-chadj))
    
    for (l in 1:(nchar(checkoutR_dn[1,'sequence'])-chadj)) {
      ss= substr(checkoutR_dn[,'sequence'], start = l, stop = (l+chadj))
      chucnkseqA[,l]=str_count(ss, "A") /ch
      chucnkseqU[,l]=str_count(ss, "U") /ch
      chucnkseqG[,l]=str_count(ss, "G") /ch
      chucnkseqC[,l]=str_count(ss, "C") /ch
    }
    
    if (r==1) {
      chucnkseqAout=as.data.frame(matrix(nrow=nrepeat,ncol=(nchar(checkoutR_dn[1,'sequence'])-chadj)))
      colnames(chucnkseqAout)=seq(1,(nchar(checkoutR_dn[1,'sequence'])-chadj))
      chucnkseqUout=as.data.frame(matrix(nrow=nrepeat,ncol=(nchar(checkoutR_dn[1,'sequence'])-chadj)))
      colnames(chucnkseqUout)=seq(1,(nchar(checkoutR_dn[1,'sequence'])-chadj))
      chucnkseqGout=as.data.frame(matrix(nrow=nrepeat,ncol=(nchar(checkoutR_dn[1,'sequence'])-chadj)))
      colnames(chucnkseqGout)=seq(1,(nchar(checkoutR_dn[1,'sequence'])-chadj))
      chucnkseqCout=as.data.frame(matrix(nrow=nrepeat,ncol=(nchar(checkoutR_dn[1,'sequence'])-chadj)))
      colnames(chucnkseqCout)=seq(1,(nchar(checkoutR_dn[1,'sequence'])-chadj))
    }  
    
    if (r>1) {  
      chucnkseqAout=fread(file=paste0(outdir_DN,'/chucnkseqA_window',ch,'.txt'),header = T,sep = "\t",stringsAsFactors = F,data.table=F)
      chucnkseqUout=fread(file=paste0(outdir_DN,'/chucnkseqU_window',ch,'.txt'),header = T,sep = "\t",stringsAsFactors = F,data.table=F)
      chucnkseqGout=fread(file=paste0(outdir_DN,'/chucnkseqG_window',ch,'.txt'),header = T,sep = "\t",stringsAsFactors = F,data.table=F)  
      chucnkseqCout=fread(file=paste0(outdir_DN,'/chucnkseqC_window',ch,'.txt'),header = T,sep = "\t",stringsAsFactors = F,data.table=F)
    }
    chucnkseqAoutadd=rbind(chucnkseqAout,colMeans(chucnkseqA))
    chucnkseqUoutadd=rbind(chucnkseqUout,colMeans(chucnkseqU))
    chucnkseqGoutadd=rbind(chucnkseqGout,colMeans(chucnkseqG))
    chucnkseqCoutadd=rbind(chucnkseqCout,colMeans(chucnkseqC))
    
    chucnkseqAoutadd=chucnkseqAoutadd[is.na(chucnkseqAoutadd[,1])==F,]
    chucnkseqUoutadd=chucnkseqUoutadd[is.na(chucnkseqUoutadd[,1])==F,]
    chucnkseqGoutadd=chucnkseqGoutadd[is.na(chucnkseqGoutadd[,1])==F,]
    chucnkseqCoutadd=chucnkseqCoutadd[is.na(chucnkseqCoutadd[,1])==F,]
    
    write.table(chucnkseqAoutadd, file=paste0(outdir_DN,'/chucnkseqA_window',ch,'.txt'), quote=F, sep="\t", row.names=F, col.names=T)
    write.table(chucnkseqUoutadd, file=paste0(outdir_DN,'/chucnkseqU_window',ch,'.txt'), quote=F, sep="\t", row.names=F, col.names=T)
    write.table(chucnkseqGoutadd, file=paste0(outdir_DN,'/chucnkseqG_window',ch,'.txt'), quote=F, sep="\t", row.names=F, col.names=T)
    write.table(chucnkseqCoutadd, file=paste0(outdir_DN,'/chucnkseqC_window',ch,'.txt'), quote=F, sep="\t", row.names=F, col.names=T)
    
    }
    
    
    
    
    ##########################################################################################  
    ##### Gquad idnettification
    ##########################################################################################  
    
    ###############################################################
    ## Upstream ###
    
    FAupstream=fread(paste0(outdir_UP,'/AllLocation_Random.fa') ,header = F,sep = "\t",stringsAsFactors = F,data.table=F)
    Gquadout=as.data.frame(matrix(nrow=(nrow(FAupstream)/2),ncol=6))
    colnames(Gquadout)=c('ID','qgrsSEQ','qgrsScore','qgrsStart','qgrsX','qgrs_maxScore')
    
    
    for (f in seq(1,nrow(FAupstream),2)) {
      fax=FAupstream[c(f,f+1),]
      nm=fax[1]
      nm=gsub(">","",nm);
      # nm=gsub("\\+","",nm);
      # nm=gsub("-","",nm);
      # nm=gsub("\\()","",nm);
      write.table(fax, file=paste0(scratch_UP,'/',nm,'.fa'), quote=F, sep="\t", row.names=F, col.names=F)
      ###############
      ## qgrs
      ###############
      (system(paste0(qgrs,' -s -csv -i ',scratch_UP,'/',nm,'.fa -o ',scratch_UP,'/',nm,'.txt'),intern = F,ignore.stderr = T));
    
        Gquadout[f,'ID']=fax[1]
      Gquadout$qgrs_maxScore=NA
      quad_qgrs=read.delim(paste0(scratch_UP,'/',nm,'.txt'),comment.char = "%",header = T,sep = ",") 
      if (nrow(quad_qgrs)>0){
        Gquadout[f,'qgrsSEQ']=paste(quad_qgrs$SEQ,collapse = ", ")
        Gquadout[f,'qgrsScore']=paste(quad_qgrs$GS,collapse = ", ")
        Gquadout[f,'qgrs_maxScore']=unique(max(quad_qgrs$GS))
        Gquadout[f,'qgrsStart']=paste(quad_qgrs$T1,collapse = ", ")
        Gquadout[f,'qgrsX']=paste(quad_qgrs$X,collapse = ", ")}
      
      if(file.exists(paste0(scratch_UP,'/',nm,'.txt'))){file.remove(paste0(scratch_UP,'/',nm,'.txt'))}
      if(file.exists(paste0(scratch_UP,'/',nm,'.fa'))){file.remove(paste0(scratch_UP,'/',nm,'.fa'))}
    }
    Gquadout=Gquadout[is.na(Gquadout$ID)==F,]
    Gquadout$ID=gsub(">","",Gquadout$ID);
    
    
    ############### 
    ### fastaRegexFinder
    ###############
    system(paste0(fastaRegexFinder,' -f ',outdir_UP,'/AllLocation_Random.fa -r "([gG]{2,4}\\w{1,7}){3,}[gG]{2,4}" --noreverse > ',outdir_UP,'/out.Random.fold.txt'),intern = T,ignore.stderr = T)
    # system(paste0(fastaRegexFinder,' -f ',outdir_UP,'/AllLocation_Random.fa -r "([gG]{2,4}\\w{1,7}){3,}[gG]{2,4}" --noreverse '),intern = T,ignore.stderr = T)
    
    info = file.info(paste0(outdir_UP,'/out.Random.fold.txt'))
    if (info$size>0){
      quad=read.delim(paste0(outdir_UP,'/out.Random.fold.txt'),comment.char = "%",header = F,sep = "\t",stringsAsFactors = F)
    
    quadx2=quad
    dup=quadx2[duplicated(quadx2$V1),]
    
    for (d in 1:nrow(dup)) {
      dd=dup[d,]
      quadx=quadx2[quad$V1%in%dd$V1,]
      quadx2[quadx2$V1%in%dd$V1,'V2']=paste(quadx$V2,collapse = ",")
      quadx2[quadx2$V1%in%dd$V1,'V3']=paste(quadx$V3,collapse = ",")
      quadx2[quadx2$V1%in%dd$V1,'V5']=paste(quadx$V5,collapse = ",")
      quadx2[quadx2$V1%in%dd$V1,'V7']=paste(quadx$V7,collapse = ",")
    }
    quadx2=quadx2[duplicated(quadx2$V1)==F,]
      # if (nrow(quad)!=nrow(quadx2)){badquad}
    quad=quadx2
    Gquadout=merge(quad[,c('V1','V2','V3','V5','V6','V7')],Gquadout,by.x='V1',by.y='ID',all=T)
    colnames(Gquadout)[colnames(Gquadout)%in%'V1']='ID'
    
    if(file.exists(paste0(outdir_UP,'/out.Random.fold.txt'))){file.remove(paste0(outdir_UP,'/out.Random.fold.txt'))}
    
    }
    
    checkoutR=merge(checkoutR,Gquadout,by.x="n",by.y="ID",all.x=T)
    
    
          checkoutR$qgrsSEQ= gsub(" ","",checkoutR$qgrsSEQ)
          checkoutR=cbind(checkoutR,str_locate(string=checkoutR[,'sequence'],pattern=checkoutR[,'qgrsSEQ']))
    
    ## upstream - start of sequence is farthest away from crosslink start-------------------------CL
          checkoutR$qgrsdist=window-checkoutR$end
          checkoutR$qgrsdist_maxScore=window-checkoutR$end
          checkoutR$qgrs_maxScore2=checkoutR$qgrsScore
            
          com=checkoutR[grep(",",checkoutR$qgrsSEQ),]
          if (nrow(com)>0) {
          for (cm in 1:nrow(com)) {
            com1=com[cm,]
            sq=str_split(com1$qgrsSEQ,",")[[1]]
            scr=str_split(com1$qgrsScore,",")[[1]]
            maxn=which(as.numeric(scr)%in%max(as.numeric(scr)))
            
            loc=str_locate(string=as.character(com1$sequence),pattern=sq)
            checkoutR[checkoutR$ID%in%com1$ID,'start'] = paste(loc[,'start'],collapse = ',')
            checkoutR[checkoutR$ID%in%com1$ID,'end'] = paste(loc[,'end'],collapse = ',')
            
            checkoutR[checkoutR$ID%in%com1$ID,'qgrsdist'] = paste(window-loc[,'end'],collapse = ',')
            checkoutR[checkoutR$ID%in%com1$ID,'qgrsdist_maxScore']=min(window-loc[maxn,'end'])
            checkoutR[checkoutR$ID%in%com1$ID,'qgrs_maxScore2']=unique(max(as.numeric(scr)))
            
                      }}
    
    ########################################################################
    ## downstream ###
    
    FAdownstream=fread(paste0(outdir_DN,'/AllLocation_DN_Random.fa'),header = F,sep = "\t",stringsAsFactors = F,data.table=F)
    Gquadout_dn=as.data.frame(matrix(nrow=(nrow(FAdownstream)/2),ncol=6))
    colnames(Gquadout_dn)=c('ID','qgrsSEQ','qgrsScore','qgrsStart','qgrsX','qgrs_maxScore')
    
    
    for (f in seq(1,nrow(FAdownstream),2)) {
      fax=FAdownstream[c(f,f+1),]
      nm=fax[1]
      nm=gsub(">","",nm);
      nm=gsub("\\+","",nm);
      nm=gsub("-","",nm);
      nm=gsub("\\()","",nm);
      write.table(fax, file=paste0(scratch_DN,'/',nm,'.fa'), quote=F, sep="\t", row.names=F, col.names=F)
      (system(paste0(qgrs,' -s -csv -i ',scratch_DN,'/',nm,'.fa -o ',scratch_DN,'/',nm,'.txt'),intern = F,ignore.stderr = T));
      
      Gquadout_dn[f,'ID']=fax[1]
      quad_qgrs_dn=read.delim(paste0(scratch_DN,'/',nm,'.txt'),comment.char = "%",header = T,sep = ",") 
      if (nrow(quad_qgrs_dn)>0){
        Gquadout_dn[f,'qgrsSEQ']=paste(quad_qgrs_dn$SEQ,collapse = ", ")
        Gquadout_dn[f,'qgrsScore']=paste(quad_qgrs_dn$GS,collapse = ", ")
        Gquadout_dn[f,'qgrs_maxScore']=unique(max(quad_qgrs_dn$GS))
        Gquadout_dn[f,'qgrsStart']=paste(quad_qgrs_dn$T1,collapse = ", ")
        Gquadout_dn[f,'qgrsX']=paste(quad_qgrs_dn$X,collapse = ", ")}
      
      if(file.exists(paste0(scratch_DN,'/',nm,'.txt'))){file.remove(paste0(scratch_DN,'/',nm,'.txt'))}
      if(file.exists(paste0(scratch_DN,'/',nm,'.fa'))){file.remove(paste0(scratch_DN,'/',nm,'.fa'))}
    }
    Gquadout_dn=Gquadout_dn[is.na(Gquadout_dn$ID)==F,]
    Gquadout_dn$ID=gsub(">","",Gquadout_dn$ID);
    
    
    system(paste0(fastaRegexFinder,' -f ',outdir_DN,'/AllLocation_DN_Random.fa -r "([gG]{2,4}\\w{1,7}){3,}[gG]{2,4}" --noreverse > ',outdir_DN,'/out.Random.fold.txt'),intern = T,ignore.stderr = T)
    # system(paste0(fastaRegexFinder,' -f ',outdir_DN,'/AllLocation_DN_Random.fa -r "([gG]{2,4}\\w{1,7}){3,}[gG]{2,4}" --noreverse '),intern = T,ignore.stderr = T)
    
    info = file.info(paste0(outdir_DN,'/out.Random.fold.txt'))
    if (info$size>0){quad_dn=read.delim(paste0(outdir_DN,'/out.Random.fold.txt'),comment.char = "%",header = F,sep = "\t",stringsAsFactors = F)
    
    dup=quad_dn[duplicated(quad_dn$V1),]
    quad_dnx2=quad_dn
     
    for (d in 1:nrow(dup)) {
      dd=dup[d,]
      quad_dnx=quad_dn[quad_dn$V1%in%dd$V1,]
      quad_dnx2[quad_dnx2$V1%in%dd$V1,'V2']=paste(quad_dnx$V2,collapse = ",")
      quad_dnx2[quad_dnx2$V1%in%dd$V1,'V3']=paste(quad_dnx$V3,collapse = ",")
      quad_dnx2[quad_dnx2$V1%in%dd$V1,'V5']=paste(quad_dnx$V5,collapse = ",")
      quad_dnx2[quad_dnx2$V1%in%dd$V1,'V7']=paste(quad_dnx$V7,collapse = ",")
    }
    quad_dn=quad_dnx2[duplicated(quad_dnx2$V1)==F,]
  
      Gquadout_dn=merge(quad_dn[,c('V1','V2','V3','V5','V6','V7')],Gquadout_dn,by.x='V1',by.y='ID',all=T)
      colnames(Gquadout_dn)[colnames(Gquadout_dn)%in%'V1']='ID'
      
    if(file.exists(paste0(outdir_DN,'/out.Random.fold.txt'))){file.remove(paste0(outdir_DN,'/out.Random.fold.txt'))}
    
    }
    
    checkoutR_dn=merge(checkoutR_dn,Gquadout_dn,by.x="n",by.y="ID",all.x=T)
    
    checkoutR_dn$qgrsSEQ= gsub(" ","",checkoutR_dn$qgrsSEQ)
    checkoutR_dn=cbind(checkoutR_dn,str_locate(string=checkoutR_dn[,'sequence'],pattern=checkoutR_dn[,'qgrsSEQ']))
    
    ## upstream - start of sequence is farthest away from crosslink CL-------------------------end
    checkoutR_dn$qgrsdist=checkoutR_dn$start
    checkoutR_dn$qgrsdist_maxScore=checkoutR_dn$start
    checkoutR_dn$qgrs_maxScore2=checkoutR_dn$qgrsScore
    
    com=checkoutR_dn[grep(",",checkoutR_dn$qgrsSEQ),]
    if (nrow(com)>0) {
    for (cm in 1:nrow(com)) {
      com1=com[cm,]
      sq=str_split(com1$qgrsSEQ,",")[[1]]
      scr=str_split(com1$qgrsScore,",")[[1]]
      maxn=which(as.numeric(scr)%in%max(as.numeric(scr)))
      
      loc=str_locate(string=as.character(com1$sequence),pattern=sq)
      checkoutR_dn[checkoutR_dn$ID%in%com1$ID,'start'] = paste(loc[,'start'],collapse = ',')
      checkoutR_dn[checkoutR_dn$ID%in%com1$ID,'end'] = paste(loc[,'end'],collapse = ',')
      
      checkoutR_dn[checkoutR_dn$ID%in%com1$ID,'qgrsdist'] = paste(loc[,'start'],collapse = ',')
      checkoutR_dn[checkoutR_dn$ID%in%com1$ID,'qgrsdist_maxScore']=min(loc[maxn,'start'])
      checkoutR_dn[checkoutR_dn$ID%in%com1$ID,'qgrs_maxScore2']=unique(max(as.numeric(scr)))
    }}
    
    # } #eliminate to calculate structure : repeat loop
    # }) #eliminate to calculate structure : RUNtime


    
    ##########################################################################################  
    ##### UTR overlap
    ##########################################################################################  
    
    TranscUTR=mm10all[mm10all$feature%in%'UTR',]
    TranscUTR=TranscUTR[TranscUTR$transcript_id%in%exons$transcript_id,]

    TranscUTR_GR=GRanges(seqnames=TranscUTR$seqname,
                         ranges=IRanges(start=TranscUTR$start,end=TranscUTR$end),
                         strand=TranscUTR$strand,
                         transcript_id=TranscUTR$transcript_id,
                         gene_id=TranscUTR$gene_id
    )


    checkoutR=separate(checkoutR,ID,into = c('chr','XLstrnd'),sep=':',remove = F)

    CrossLink_GR=GRanges(seqnames=checkoutR$chr,
                         ranges=IRanges(start=checkoutR$crosslink, end=checkoutR$crosslink),
                         strand=checkoutR$strand,
                         ID=checkoutR$ID,
                         trancript=checkoutR$trancript)


    gr = CrossLink_GR
    gr2 = TranscUTR_GR
    xo=as.data.frame(GenomicRanges::findOverlaps(gr,gr2,type = "any",ignore.strand=T))

   UTR_CrossLink_GR=gr2[xo$subjectHits]
    # qh=as.data.frame(gr[xo$queryHits],row.names = NULL)
    # colnames(qh)=paste0(colnames(qh),"_mm10")
    # sh=as.data.frame(gr2[xo$subjectHits],row.names = NULL)
    Crosslink_UTR_GR=gr[xo$queryHits]
    Crosslink_UTR_DF=as.data.frame(Crosslink_UTR_GR)
    
    checkoutR[which(checkoutR$ID%in%Crosslink_UTR_DF$ID),'UTR']="UTR"
    
    # if (nrow(as.data.frame(Crosslink_UTR_GR))>0) {
    # Crosslink_UTR_DF$UTR='UTR'
    # checkoutR=merge(checkoutR, Crosslink_UTR_DF[,c('ID','UTR')],by="ID",all.x=T)
    # }
    # 
    # if (nrow(as.data.frame(Crosslink_UTR_GR))==0){checkoutR$UTR=NA}
    # 
      ####################################################################
      #### Upstream Base Pairing Probability
      ####################################################################

    system(paste0('rm -f ',scratch_UP,'/*.ps'),intern = F)
    system(paste0('rm -f ',scratch_UP,'/out.up.fold'),intern = F)
    
      # system(paste0('RNAfold -p -i ',outdir_UP,'/AllLocation_Random.fa > out.up.fold'),intern = F)
      system(paste0(RNAfold,' -p -i ',outdir_UP,'/AllLocation_Random.fa > out.up.fold'),intern = F)
      system(paste0('mv *.ps ',scratch_UP),intern = F)
      system(paste0('mv out.up.fold ',scratch_UP),intern = F)
      
      enrg=read.delim(paste0(scratch_UP,'/out.up.fold'),comment.char = "%",header = F,sep = "\t",stringsAsFactors = F)
    
      for (e in 1:nrow(checkoutR)) {
        s=as.character(checkoutR[e,'n'])
        s=gsub("\\(-)","",s)
        s=gsub("\\(+)","",s)
        loc=grep(s,enrg$V1)
          if (isEmpty(loc)) {next ;xxx  }
      enrgx=enrg[c(loc:(loc+5)),]
      enrgx=as.data.frame(enrgx[4])
      enrgx=separate(enrgx,col=1,into = c('seq','DG'),window+2)
      enrgx$DG=gsub("]","",enrgx$DG)
      checkoutR[e,'DeltaG']=as.numeric(enrgx$DG)
      }
      
      
      fls=list.files(paste0(scratch_UP,'/'), pattern="dp.ps", all.files=FALSE,full.names=FALSE)
      fls=paste0(scratch_UP,'/',fls)
      fls=as.list(fls)

      # xxx=gsub("\\(","\\\\(",fls)
      # xxx=str_replace_all(fls,"\\(","/\\(")
      # str_replace_all(xxx,"/","\\\\")
      

      out=mclapply(fls,RNAfold_output,mc.cores=UseCores)
      out2=out[[1]][,1:2,drop=F]
      out2_max=out[[1]][,1:2,drop=F]
      for (o in 1:length(out)) {
        if (nrow(out[[o]][,3,drop=F])==0) {
          out2$Xx=0;colnames(out2)[colnames(out2)%in%"Xx"]=colnames(out[[o]][,3,drop=F])
          out2_max$Xx=0;colnames(out2_max)[colnames(out2_max)%in%"Xx"]=colnames(out[[o]][,4,drop=F])
        }
        if (nrow(out[[o]][,3,drop=F])>0) {
          out2=cbind(out2,out[[o]][,3,drop=F])
          out2_max=cbind(out2_max,out[[o]][,4,drop=F])
        }
      }

      system(paste0('rm -f ',scratch_UP,'/*.ps'))
      system(paste0('rm -f ',scratch_UP,'/*.fa'))
      system(paste0('rm -f ',scratch_UP,'/*.bed'))
      system(paste0('rm -f ',scratch_UP,'/out.dn.fold'))

      # #########################################
      ph=as.matrix(out2[,-c(1:2)])
      ph=ph[is.na(ph[,3])==F,]
      ph[is.nan(as.matrix(ph))==T]=0

      ph=t(ph)
      ph=as.data.frame(ph)

      #### sort by row diff
      ph=ph[order(rowMaxs(as.matrix(ph)),decreasing = T),]
      colnames(ph)=c(-window:(ncol(ph)-(window+1)))
      ph= cbind(rowMaxs(as.matrix(ph)),ph);colnames(ph)[1]='diff'
      
      
      outchunck_rowaverage[r,]=colMeans(ph[ph$diff>=meanBPP,-1,drop=F])
      outchunck_frac[r,1]=nrow(ph[ph$diff>=meanBPP,-1,drop=F])
      outchunck_frac[r,2]=nrow(ph)

      
      ph$n=rownames(ph)
      ph$n=gsub(paste0(scratch_UP,'/'),"",ph$n)
      ph$n=gsub("_",":",ph$n)
      
      checkoutR=merge(checkoutR,ph[,c("n","diff"),drop=F],by='n',all.x=T)
      colnames(checkoutR)[colnames(checkoutR)%in%'diff']='RowMax_MeanChuckBPP'
      
      
      #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
      #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

      ph_max=as.matrix(out2_max[,-c(1:2)])
      ph_max=ph_max[is.na(ph_max[,3])==F,]
      ph_max[is.nan(as.matrix(ph_max))==T]=0

      ph_max=t(ph_max)
      ph_max=as.data.frame(ph_max)
      
      
      #### sort by row diff
      ph_max=ph_max[order(rowMaxs(as.matrix(ph_max)),decreasing = T),]
      colnames(ph_max)=c(-window:(ncol(ph_max)-(window+1)))
      ph_max= cbind(rowMaxs(as.matrix(ph_max)),ph_max);colnames(ph_max)[1]='diff'

      
      
      outchunck_rowaverage_max[r,]=colMeans(ph_max[ph_max$diff>=maxBPP,-1,drop=F])
      outchunck_frac_max[r,1]=nrow(ph_max[ph_max$diff>=maxBPP,-1,drop=F])
      outchunck_frac_max[r,2]=nrow(ph_max)

      
      ph_max$n=rownames(ph_max)
      ph_max$n=gsub(paste0(scratch_UP,'/'),"",ph_max$n)
      ph_max$n=gsub("_",":",ph_max$n)
      
      checkoutR=merge(checkoutR,ph_max[,c("n","diff"),drop=F],by='n',all.x=T)
      colnames(checkoutR)[colnames(checkoutR)%in%'diff']='RowMax_MaxChuckBPP'
      
      
      
      ####################################################################
      #### downstream Base Pairing Probability
      ####################################################################
      ### High base pair probability chr8:18628170-18628210
      ### low basepair probibility chr9:72581297-72581336

      system(paste0('rm -f ',scratch_DN,'/*.ps'))
      system(paste0('rm -f ',scratch_DN,'/out.dn.fold'))
      
      system(paste0(RNAfold,' -p -i ',outdir_DN,'/AllLocation_DN_Random.fa > out.dn.fold'),intern=F)
      system(paste0('mv *.ps ',scratch_DN))
      system(paste0('mv out.dn.fold ',scratch_DN),intern = F)
      
      energDN=read.delim(paste0(scratch_DN,'/out.dn.fold'),comment.char = "%",header = F,sep = "\t",stringsAsFactors = F)
       
      for (e in 1:nrow(checkoutR_dn)) {
        s=as.character(checkoutR_dn[e,'n'])
        s=gsub("\\(-)","",s)
        s=gsub("\\(+)","",s)
        loc=grep(s,energDN$V1)
        if (isEmpty(loc)) {next ;xxx  }
        
        energDNx=energDN[c(loc:(loc+5)),]
        energDNx=as.data.frame(energDNx[4])
        energDNx=separate(energDNx,col=1,into = c('seq','DG'),window+2)
        energDNx$DG=gsub("]","",energDNx$DG)
        checkoutR_dn[e,'DeltaG']=as.numeric(energDNx$DG)
      }
    
      
      
      flsdn=list.files(paste0(scratch_DN,'/'), pattern="dp.ps", all.files=FALSE, full.names=FALSE)
      flsdn=paste0(scratch_DN,'/',flsdn)
      
      flsdn= as.list(flsdn)

      outdn=mclapply((flsdn),RNAfold_output,mc.cores=UseCores)
      outdn2=outdn[[1]][,1:2,drop=F]
      outdn2_max=outdn[[1]][,1:2,drop=F]
      for (o in 1:length(outdn)) {
        if (nrow(outdn[[o]][,3,drop=F])==0) {
          outdn2$Xx=0;colnames(outdn2)[colnames(outdn2)%in%"Xx"]=colnames(outdn[[o]][,3,drop=F])
          outdn2_max$Xx=0;colnames(outdn2_max)[colnames(outdn2_max)%in%"Xx"]=colnames(outdn[[o]][,4,drop=F])
        }
        if (nrow(outdn[[o]][,3,drop=F])>0) {
          outdn2=cbind(outdn2,outdn[[o]][,3,drop=F])
          outdn2_max=cbind(outdn2_max,outdn[[o]][,4,drop=F])
        }
      }

      system(paste0('rm -f ',scratch_DN,'/*.ps'))
      system(paste0('rm -f ',scratch_DN,'/*.fa'))
      system(paste0('rm -f ',scratch_DN,'/*dn.bed'))
      system(paste0('rm -f ',scratch_DN,'/out.dn.fold'))
      
      # #########################################
      phdn=as.matrix(outdn2[,-c(1:2)])
      phdn=phdn[is.na(phdn[,3])==F,]
      phdn[is.nan(as.matrix(phdn))==T]=0
      phdn=t(phdn)
      phdn=as.data.frame(phdn)

      #### sort by row diff
      phdn=phdn[order(rowMaxs(as.matrix(phdn)),decreasing = T),]
      colnames(phdn)=c(-window:(ncol(phdn)-(window+1)))
      phdn= cbind(rowMaxs(as.matrix(phdn)),phdn);colnames(phdn)[1]='diff'

      
      outchunck_rowaverage_dn[r,]=colMeans(phdn[phdn$diff>=meanBPP,-1,drop=F])
      outchunck_frac_dn[r,1]=nrow(phdn[phdn$diff>=meanBPP,-1,drop=F])
      outchunck_frac_dn[r,2]=nrow(phdn)

      
      phdn$n=rownames(phdn)
      phdn$n=gsub(paste0(scratch_DN,'/'),"",phdn$n)
      phdn$n=gsub("_",":",phdn$n)
      
      checkoutR_dn=merge(checkoutR_dn,phdn[,c("n","diff"),drop=F],by='n',all.x=T)
      colnames(checkoutR_dn)[colnames(checkoutR_dn)%in%'diff']='RowMax_MeanChuckBPP'
      
      # #########################################
      # #########################################
      phdn_max=as.matrix(outdn2_max[,-c(1:2)])
      phdn_max=phdn_max[is.na(phdn_max[,3])==F,]
      phdn_max[is.nan(as.matrix(phdn_max))==T]=0
      phdn_max=t(phdn_max)
      phdn_max=as.data.frame(phdn_max)

      #### sort by row diff
      phdn_max=phdn_max[order(rowMaxs(as.matrix(phdn_max)),decreasing = T),]
      colnames(phdn_max)=c(-window:(ncol(phdn_max)-(window+1)))
      phdn_max= cbind(rowMaxs(as.matrix(phdn_max)),phdn_max);colnames(phdn_max)[1]='diff'
      
      outchunck_rowaverage_maxdn[r,]=colMeans(phdn_max[phdn_max$diff>=maxBPP,-1,drop=F])
      outchunck_frac_maxdn[r,1]=nrow(phdn_max[phdn_max$diff>=maxBPP,-1,drop=F])
      outchunck_frac_maxdn[r,2]=nrow(phdn_max)

      phdn_max$n=rownames(phdn_max)
            phdn_max$n=gsub(paste0(scratch_DN,'/'),"",phdn_max$n)
            phdn_max$n=gsub("_",":",phdn_max$n)
      
      checkoutR_dn=merge(checkoutR_dn,phdn_max[,c("n","diff"),drop=F],by='n',all.x=T)
      colnames(checkoutR_dn)[colnames(checkoutR_dn)%in%'diff']='RowMax_MaxChuckBPP'
      

  

  rownames(outchunck_rowaverage)=c(1:nrepeat)
  # colnames(outchunck_rowaverage)=c(-window:-1)
  colnames(outchunck_rowaverage)=c(-window:-10)

  write.table(outchunck_rowaverage, file=paste0(outdir_UP,'/Random_outchunck_rowaverage_',window,'chuck4.txt'), quote=F, sep="\t", row.names=T, col.names=NA)
  write.table(outchunck_frac, file=paste0(outdir_UP,'/Random_outchunck_frac_',window,'chunk4.txt'), quote=F, sep="\t", row.names=F, col.names=F)

  rownames(outchunck_rowaverage_max)=c(1:nrepeat)
  # colnames(outchunck_rowaverage)=c(-window:-1)
  colnames(outchunck_rowaverage_max)=c(-window:-10)

  write.table(outchunck_rowaverage_max, file=paste0(outdir_UP,'/Random_outchunck_rowaverage_max_',window,'chuck4.txt'), quote=F, sep="\t", row.names=T, col.names=NA)
  write.table(outchunck_frac_max, file=paste0(outdir_UP,'/Random_outchunck_frac_max_',window,'chunk4.txt'), quote=F, sep="\t", row.names=F, col.names=F)



  rownames(outchunck_rowaverage_dn)=c(1:nrepeat)
  # colnames(outchunck_rowaverage_dn)=c(1:window)
  colnames(outchunck_rowaverage_dn)=c(10:window)

  write.table(outchunck_rowaverage_dn, file=paste0(outdir_DN,'/Random_dn_outchunck_rowaverage_',window,'chucnk4.txt'), quote=F, sep="\t", row.names=T, col.names=NA)
  write.table(outchunck_frac_dn, file=paste0(outdir_DN,'/Random_dn_outchunck_frac_',window,'chunk4.txt'), quote=F, sep="\t", row.names=F, col.names=F)


  rownames(outchunck_rowaverage_maxdn)=c(1:nrepeat)
  # colnames(outchunck_rowaverage_dn)=c(1:window)
  colnames(outchunck_rowaverage_maxdn)=c(10:window)

  write.table(outchunck_rowaverage_maxdn, file=paste0(outdir_DN,'/Random_dn_outchunck_rowaverage_max_',window,'chucnk4.txt'), quote=F, sep="\t", row.names=T, col.names=NA)
  write.table(outchunck_frac_maxdn, file=paste0(outdir_DN,'/Random_dn_outchunck_frac_max_',window,'chunk4.txt'), quote=F, sep="\t", row.names=F, col.names=F)

  
  
  ##########################################################################################  
  ##### Sumarize Exon sequences
  ##########################################################################################  
  
  
  checkoutRx=checkoutR
  checkoutRx[,"SingleExon"]=gsub("Change-TRUE","TRUE",checkoutRx[,"SingleExon"])
  checkoutRx[,"SingleExon"]=as.character(checkoutRx[,"SingleExon"])
  
  checkoutR_summary[r,"total Locations"]=nrow(checkoutRx)
  checkoutR_summary[r,"SingleExon"]=sum(checkoutRx[,"SingleExon"]=='TRUE')
  checkoutR_summary[r,"numberExCrossed"]=sum(checkoutRx[,'numberExCrossed']>1)
  checkoutR_summary[r,"lastexon"]=sum(checkoutRx[,'lastexon']==1)
  checkoutR_summary[r,"lastexonCrossexExon"]=sum(checkoutRx[,'lastexon']==1&checkoutRx[,'SingleExon']=='FALSE')
  checkoutR_summary[r,"firstexon"]=sum(checkoutRx[,'exon#']==1)
  checkoutR_summary[r,"firstexonCrossexExon"]=sum(checkoutRx[,'exon#']==1&checkoutRx[,'SingleExon']=='FALSE')
  checkoutR_summary[r,"Gquad"]=nrow(checkoutRx[is.na(checkoutRx[,'qgrsScore'])==F,])
  checkoutR_summary[r,"Avg_GquadScore"]=mean(as.numeric(unlist(str_split(checkoutRx$qgrsScore,", "))),na.rm=T)
  checkoutR_summary[r,"MeanBPP"]=nrow(checkoutRx[(checkoutRx[,'RowMax_MeanChuckBPP'])>=meanBPP,])
  checkoutR_summary[r,"MaxBPP"]=nrow(checkoutRx[(checkoutRx[,'RowMax_MaxChuckBPP'])>=maxBPP,])
  checkoutR_summary[r,"UTR"]=nrow(checkoutRx[checkoutRx$UTR%in%"UTR",])
  
  write.table(checkoutR_summary, file=paste0(outdir_UP,'/checkoutR_summary.txt'), quote=F, sep="\t", row.names=F, col.names=T)
  
  
   checkoutR_dnx=checkoutR_dn
  checkoutR_dnx[,"SingleExon"]=gsub("Change-TRUE","TRUE",checkoutR_dnx[,"SingleExon"])
  checkoutR_dnx[,"SingleExon"]=as.character(checkoutR_dnx[,"SingleExon"])
  
  checkoutR_dn_summary[r,"total Locations"]=nrow(checkoutR_dnx)
  checkoutR_dn_summary[r,"SingleExon"]=sum(checkoutR_dnx[,"SingleExon"]=='TRUE')
  checkoutR_dn_summary[r,"numberExCrossed"]=sum(checkoutR_dnx[,'numberExCrossed']>1)
  checkoutR_dn_summary[r,"lastexon"]=sum(checkoutR_dnx[,'lastexon']==1)
  checkoutR_dn_summary[r,"lastexonCrossexExon"]=sum(checkoutR_dnx[,'lastexon']==1&checkoutR_dnx[,'SingleExon']=='FALSE')
  checkoutR_dn_summary[r,"firstexon"]=sum(checkoutR_dnx[,'exon#']==1)
  checkoutR_dn_summary[r,"firstexonCrossexExon"]=sum(checkoutR_dnx[,'exon#']==1&checkoutR_dnx[,'SingleExon']=='FALSE')
  checkoutR_dn_summary[r,"Gquad"]=nrow(checkoutR_dnx[is.na(checkoutR_dnx[,'qgrsScore'])==F,])
  checkoutR_dn_summary[r,"Avg_GquadScore"]=mean(as.numeric(unlist(str_split(checkoutR_dnx$qgrsScore,", "))),na.rm=T)
  checkoutR_dn_summary[r,"MeanBPP"]=nrow(checkoutR_dnx[(checkoutR_dnx[,'RowMax_MeanChuckBPP'])>=meanBPP,])
  checkoutR_dn_summary[r,"MaxBPP"]=nrow(checkoutR_dnx[(checkoutR_dnx[,'RowMax_MaxChuckBPP'])>=maxBPP,])
  
  
  write.table(checkoutR_dn_summary, file=paste0(outdir_DN,'/checkoutR_dn_summary.txt'), quote=F, sep="\t", row.names=F, col.names=T)
  
  
  
  
  
  
  
  ########################################################
  ### save columns from checkout_table
  ############################################################
  
  n_UP[,r]=checkoutR[,'n']
  ID_UP[,r]=checkoutR[,'ID']
  seq_UP[,r]=checkoutR[,'sequence']
  SingleExon_UP[,r]=checkoutR[,'SingleExon']
  Exon_UP[,r]=checkoutR[,'exon#']
  transcript_UP[,r]=checkoutR[,'trancript']
  LastExon_UP[,r]=checkoutR[,'lastexon']
  Gscore_UP[,r]=checkoutR[,'qgrsScore']
  GscoreDist_UP[,r]=checkoutR[,'qgrsdist']
  GscoreMaxDist_UP[,r]=checkoutR[,'qgrsdist_maxScore']
  GscoreMax_UP[,r]=checkoutR[,'qgrs_maxScore']
  GscoreMax2_UP[,r]=checkoutR[,'qgrs_maxScore2']
  FinderSeq_UP[,r]=checkoutR[,'V7']
  qgrsSeq_UP[,r]=checkoutR[,'qgrsSEQ']
    maxBPPmean_UP[,r]=checkoutR[,'RowMax_MeanChuckBPP']
    maxBPPmax_UP[,r]=checkoutR[,'RowMax_MaxChuckBPP']
    MFE_UP[,r]=checkoutR[,'DeltaG']
    
    
    write.table(n_UP, file=paste0(outdir_UP,'/Table/Col_n_UP.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(ID_UP, file=paste0(outdir_UP,'/Table/Col_ID_UP.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(seq_UP, file=paste0(outdir_UP,'/Table/Col_seq_UP.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(SingleExon_UP, file=paste0(outdir_UP,'/Table/Col_SingleExon_UP.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(Exon_UP, file=paste0(outdir_UP,'/Table/Col_Exon_UP.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(transcript_UP, file=paste0(outdir_UP,'/Table/Col_transcript_UP.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(LastExon_UP, file=paste0(outdir_UP,'/Table/Col_LastExon_UP.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(Gscore_UP, file=paste0(outdir_UP,'/Table/Col_Gscore_UP.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(GscoreDist_UP, file=paste0(outdir_UP,'/Table/Col_GscoreDist_UP.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(GscoreMaxDist_UP, file=paste0(outdir_UP,'/Table/Col_GscoreMaxDist_UP.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(GscoreMax_UP, file=paste0(outdir_UP,'/Table/Col_GscoreMax_UP.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(GscoreMax2_UP, file=paste0(outdir_UP,'/Table/Col_GscoreMax2_UP.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(FinderSeq_UP, file=paste0(outdir_UP,'/Table/Col_FinderSeq_UP.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(qgrsSeq_UP, file=paste0(outdir_UP,'/Table/Col_qgrsSeq_UP.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(maxBPPmean_UP, file=paste0(outdir_UP,'/Table/Col_maxBPPmean_UP.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(maxBPPmax_UP, file=paste0(outdir_UP,'/Table/Col_maxBPPmax_UP.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(MFE_UP, file=paste0(outdir_UP,'/Table/Col_MFE_UP.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    
    write.table(checkoutR, file=paste0(outdir_UP,'/checkoutR.txt'), quote=F, sep="\t", row.names=F, col.names=T)
    
    
    

    
    n_DN[,r]=checkoutR_dn[,'n']
    ID_DN[,r]=checkoutR_dn[,'ID']
    seq_DN[,r]=checkoutR_dn[,'sequence']
    SingleExon_DN[,r]=checkoutR_dn[,'SingleExon']
    Exon_DN[,r]=checkoutR_dn[,'exon#']
    transcript_DN[,r]=checkoutR_dn[,'trancript']
    LastExon_DN[,r]=checkoutR_dn[,'lastexon']
    Gscore_DN[,r]=checkoutR_dn[,'qgrsScore']
    GscoreDist_DN[,r]=checkoutR_dn[,'qgrsdist']
    GscoreMaxDist_DN[,r]=checkoutR_dn[,'qgrsdist_maxScore']
    GscoreMax_DN[,r]=checkoutR_dn[,'qgrs_maxScore']
    GscoreMax2_DN[,r]=checkoutR_dn[,'qgrs_maxScore2']
    FinderSeq_DN[,r]=checkoutR_dn[,'V7']
    qgrsSeq_DN[,r]=checkoutR_dn[,'qgrsSEQ']
    maxBPPmean_DN[,r]=checkoutR_dn[,'RowMax_MeanChuckBPP']
    maxBPPmax_DN[,r]=checkoutR_dn[,'RowMax_MaxChuckBPP']
    MFE_DN[,r]=checkoutR_dn[,'DeltaG']
    
    write.table(n_DN, file=paste0(outdir_DN,'/Table/Col_n_DN.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(ID_DN, file=paste0(outdir_DN,'/Table/Col_ID_DN.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(seq_DN, file=paste0(outdir_DN,'/Table/Col_seq_DN.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(SingleExon_DN, file=paste0(outdir_DN,'/Table/Col_SingleExon_DN.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(Exon_DN, file=paste0(outdir_DN,'/Table/Col_Exon_DN.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(transcript_DN, file=paste0(outdir_DN,'/Table/Col_transcript_DN.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(LastExon_DN, file=paste0(outdir_DN,'/Table/Col_LastExon_DN.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(Gscore_DN, file=paste0(outdir_DN,'/Table/Col_Gscore_DN.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(GscoreDist_DN, file=paste0(outdir_DN,'/Table/Col_GscoreDist_DN.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(GscoreMaxDist_DN, file=paste0(outdir_UP,'/Table/Col_GscoreMaxDist_DN.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(GscoreMax_DN, file=paste0(outdir_UP,'/Table/Col_GscoreMax_DN.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(GscoreMax2_DN, file=paste0(outdir_UP,'/Table/Col_GscoreMax2_DN.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(FinderSeq_DN, file=paste0(outdir_DN,'/Table/Col_FinderSeq_DN.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(qgrsSeq_DN, file=paste0(outdir_DN,'/Table/Col_qgrsSeq_DN.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(maxBPPmean_DN, file=paste0(outdir_DN,'/Table/Col_maxBPPmean_DN.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(maxBPPmax_DN, file=paste0(outdir_DN,'/Table/Col_maxBPPmax_DN.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    write.table(MFE_DN, file=paste0(outdir_DN,'/Table/Col_MFE_DN.txt'), quote=F, sep="\t", row.names=F, col.names=F)
    
    }
  
  })
  

}
RandomIntron_chuck()