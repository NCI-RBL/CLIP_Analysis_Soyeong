CLIPannotation=function(peaks,WriteClassTable){
  

  test=F
# test=T  

library(VariantAnnotation,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
library(GenomicRanges,quietly = T,verbose = F,warn.conflicts = F,logical.return = F)
#library(WhopGenome)
library(trackViewer)
library("pheatmap")

library(vcfR)
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
library(knitr)
library(stringi)
library(GeneStructureTools)
library(biomaRt)
library(plotly)
library(tidyr)
library(GenomicRanges)
library(RColorBrewer)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )


if (test==T) {
peaks=fread(paste0("./NovoAlign_umi/All/recomb_NHtag/peaks_Dist50nt_postUniq_postMM_2primary/Ro_Clip_iCountcutadpt_all.unique.NH.mm.ddup.s.unique.bam.peaks.bedtools.bed"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
peaks=peaks[peaks$V4>=20,]
peaks=peaks[,c("V1","V2","V3",'V6')]
peaks=Peaksdata2[,c('chr','start','end','strand')]
WriteClassTable=F
}


##### Begin script ###
outdir="./annotation/"


colnames(peaks)=c('chr','start','end','strand')

peaks$ID=paste0(peaks$chr,":",peaks$start,"-",peaks$end)
peaks$ID2=paste0(peaks$chr,":",peaks$start,"-",peaks$end,'_',peaks$strand)

peaks_oppo=peaks
peaks_oppo$strand=gsub("\\+","pos",peaks_oppo$strand) 
peaks_oppo$strand=gsub("\\-","+",peaks_oppo$strand)
peaks_oppo$strand=gsub("pos","-",peaks_oppo$strand)





###################################################################################################################################################################################################
###################################################################################################################################################################################################

# 1. Read in Annotation files

###################################################################################################################################################################################################
###################################################################################################################################################################################################


#################################################
#### FROM GENCODE
#################################################
  library(Rgb)
  
  # for biotypes https://www.gencodegenes.org/pages/biotypes.html
  #################################
  ### annotations will most closely match refseq 
  # mm10all=Rgb::read.gtf("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/Gencode_VM23/fromGencode/gencode.vM23.basic.annotation.gtf",attr='split')
  # write.table(mm10all,file=("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/Gencode_VM23/fromGencode/gencode.vM23.basic.annotation.gtf.txt"), sep = "\t", row.names = FALSE, col.names = T, append = F, quote= FALSE)
  mm10all=fread("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/Gencode_VM23/fromGencode/gencode.vM23.basic.annotation.gtf.txt", header=T, sep="\t",stringsAsFactors = F,data.table=F)
  
  
  #################################
  ### extra annooatations   
  # mm10all2=Rgb::read.gtf("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/Gencode_VM23/fromGencode/gencode.vM23.annotation.gtf",attr='split')
  # write.table(mm10all2,file=("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/Gencode_VM23/fromGencode/gencode.vM23.annotation.gtf.txt"), sep = "\t", row.names = FALSE, col.names = T, append = F, quote= FALSE)
  # mm10all2=fread("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/Gencode_VM23/fromGencode/gencode.vM23.annotation.gtf.txt", header=T, sep="\t",stringsAsFactors = F,data.table=F)
  
  
  colnames(mm10all)[colnames(mm10all)%in%'seqname']='chr'
  colnames(mm10all)[colnames(mm10all)%in%'gene_id']='ensembl_gene_id'
  colnames(mm10all)[colnames(mm10all)%in%'gene_name']='external_gene_name'
  mm10all$transcript_id=removeVersion(mm10all$transcript_id)
  mm10all$ensembl_gene_id=removeVersion(mm10all$ensembl_gene_id)
  mm10all$chr=as.character(mm10all$chr)


  
  ##remove TEC genes
  mm10all=mm10all[!(mm10all$transcript_type%in%'TEC'),]
  mm10all=mm10all[!(mm10all$gene_type%in%'TEC'),]  
  
  
  ### combine biotypesinto more general catagories
  mm10all$gene_biotype_ALL=mm10all$gene_type
  
  
  ## combine all Pseudogenes
  # unique(mm10all[grep('pseudogene',mm10all$gene_type),'gene_type'])
  mm10all[grep('pseudogene',mm10all$gene_type),'gene_type']='pseudogene'
  
  
  ## combine all ncRNA
  mm10all$gene_type=gsub('miRNA','ncRNA',mm10all$gene_type)
  mm10all$gene_type=gsub('miscRNA','ncRNA',mm10all$gene_type)
  mm10all$gene_type=gsub('misc_RNA','ncRNA',mm10all$gene_type)
  mm10all$gene_type=gsub('piRNA','ncRNA',mm10all$gene_type)
  mm10all$gene_type=gsub('rRNA','ncRNA',mm10all$gene_type)
  mm10all$gene_type=gsub('siRNA','ncRNA',mm10all$gene_type)
  mm10all$gene_type=gsub('snRNA','ncRNA',mm10all$gene_type)
  mm10all$gene_type=gsub('snoRNA','ncRNA',mm10all$gene_type)
  mm10all$gene_type=gsub('tRNA','ncRNA',mm10all$gene_type)
  mm10all$gene_type=gsub('ribozyme','ncRNA',mm10all$gene_type)
  mm10all$gene_type=gsub('lncRNA','ncRNA',mm10all$gene_type)


  ##############################
  mm10=mm10all[mm10all$feature%in%'transcript',]
  mm10_exon=mm10all[mm10all$feature%in%c("exon"),]
  mm10_FTR=mm10all[mm10all$feature%in%c("3UTR","5UTR",'CDS','UTR'),]
  
  ##############################
  canonical=fread("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/Gencode_VM23/fromUCSC/KnownCanonical/KnownCanonical_GencodeM23_GRCm38.txt", header=T, sep="\t",stringsAsFactors = F,data.table=F)
  canonical$transcript=removeVersion(canonical$transcript)
  canonical$protein=removeVersion(canonical$protein)
  
  
  ########
  introns=fread("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/Gencode_VM23/fromUCSC/KnownGene/KnownGene_GRCm38_introns.bed", header=F, sep="\t",stringsAsFactors = F,data.table=F)
  introns=introns[grepl("_",introns$V1)==F,]
  colnames(introns)=c('chr','start','end','attribute','V5','strand')
  introns=separate(introns,attribute,into=c('transcript_id','feature','exon_number','level','chr2','intronnumber','dir'),remove = T,sep = "_")
  introns$transcript_id=removeVersion(introns$transcript_id)
  # introns_canon=introns[introns$transcript_id%in%canonical$transcript_id,]
  introns$start=introns$start+1
  introns$exon_number=as.numeric(introns$exon_number)+1
  

  mm10_ExnItrn=rbind(mm10_exon[,c('chr','feature','start','end','strand','transcript_id','exon_number')],
                     introns[,c('chr','feature','start','end','strand','transcript_id','exon_number')])
  mm10_ExnItrn$ID=paste0(mm10_ExnItrn$chr,':',mm10_ExnItrn$start,'-',mm10_ExnItrn$end)
  
# mm10 exon_number matches strand direction Intron does not 
# View(mm10_exon[,c('chr','start','end','strand',"feature","transcript_id","exon_number","level")])

calcIntron=0
if (calcIntron==1){
  #####################################################
  library(GenomicFeatures)
  library(rtracklayer)
  
  gtf <- makeTxDbFromGFF("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/Gencode_VM23/mm10.KnownGene.gtf") #change me!
  exons <- exonsBy(gtf, by="gene")
  #make introns
  exons <- reduce(exons)
  exons <- exons[sapply(exons, length) > 1]
  
  introns <- lapply(exons, function(x) {
    #Make a "gene" GRange object
    gr = GRanges(seqnames=seqnames(x)[1], ranges=IRanges(start=min(start(x)),
                                                         end=max(end(x))),
                 strand=strand(x)[1])
    db = disjoin(c(x, gr))
    ints = db[countOverlaps(db, x) == 0]
    #Add an ID
    if(as.character(strand(ints)[1]) == "-") {
      ints$exon_id = c(length(ints):1)
    } else {
      ints$exon_id = c(1:length(ints))
    }
    ints
  })
  introns <- GRangesList(introns)
  as.data.frame(introns)
  ######################################
}
#### for possible more effecient method
####https://genomicsclass.github.io/book/pages/bioc1_igranges.html


mm10_ExnItrn=merge(mm10_ExnItrn,mm10[,c('transcript_id','gene_type')],by='transcript_id',all.x=T)


########################################################################################
########################################################################################

calcbiomart=0
if (calcbiomart==1) {
  
  library(biomaRt)
  #############################################
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  dataset="mmusculus_gene_ensembl")
  # View(listAttributes(mart))
  # View(listFilters(mart))
  # View(listMarts())
  # View(listDatasets(mart))
  
  genename=getBM(attributes = c('ensembl_gene_id', 'external_gene_name','ensembl_transcript_id'),
                 mart = mart,uniqueRows = TRUE)
  genebiotype=getBM(attributes = c('ensembl_gene_id', 'external_gene_name','ensembl_transcript_id','gene_biotype','transcript_biotype'),
                    mart = mart,uniqueRows = TRUE)
  write.table(genename,file=("genename_mart.txt"), sep = "\t", row.names = FALSE, col.names = T, append = F, quote= FALSE)
  write.table(genebiotype,file=("genebiotype_mart.txt"), sep = "\t", row.names = FALSE, col.names = T, append = F, quote= FALSE)
  
}  


if (calcbiomart==0) {
  
  genename=fread("genename_mart.txt", header=T, sep="\t",stringsAsFactors = F,data.table=F)
  genename$ensembl_transcript_id=removeVersion(genename$ensembl_transcript_id)
  genebiotype=fread("genebiotype_mart.txt", header=T, sep="\t",stringsAsFactors = F,data.table=F)
  genebiotype$ensembl_transcript_id=removeVersion(genebiotype$ensembl_transcript_id)
}

###################################
## GENCODE
## transcript_type

# unique(mm10$transcript_type)
# snRNA
# snoRNA
# scRNA
# scaRNA 
# miRNA
# rRNA  
# lncRNA
# 
# ribozyme 
# sRNA                              
# misc_RNA                          
#                   
# protein_coding  
# TEC
# Mt_tRNA                           
# Mt_rRNA
# nonsense_mediated_decay
# retained_intron
# c("processed_pseudogene","transcribed_unitary_pseudogene","polymorphic_pseudogene","translated_unprocessed_pseudogene","pseudogene","unitary_pseudogene","transcribed_processed_pseudogene","transcribed_unprocessed_pseudogene","unprocessed_pseudogene")
# c("TR_V_gene","TR_V_pseudogene","TR_D_gene","TR_J_gene","TR_C_gene","TR_J_pseudogene")
# c("IG_D_pseudogene","IG_C_pseudogene","IG_D_gene","IG_LV_gene","IG_V_gene","IG_V_pseudogene","IG_J_gene","IG_C_gene","IG_pseudogene")


snRNA_mm10=mm10[mm10$transcript_type%in%'snRNA',]
snoRNA_mm10=mm10[mm10$transcript_type%in%'snoRNA',]
scRNA_mm10=mm10[mm10$transcript_type%in%'scRNA',]
scaRNA_mm10=mm10[mm10$transcript_type%in%'scaRNA',]
miRNA_mm10=mm10[mm10$transcript_type%in%'miRNA',]
rRNA_mm10=mm10[mm10$transcript_type%in%'rRNA',]
lncRNA_mm10=mm10[mm10$transcript_type%in%'lncRNA',]

lincRNA_mm10=lncRNA_mm10[lncRNA_mm10$transcript_id%in%unique(introns$transcript_id),]
lincRNA_mm10$transcript_type='lincRNA'
# lncRNA_mm10=lncRNA_mm10[!lncRNA_mm10$transcript_id%in%lincRNA_mm10$transcript_id,]
# mm10[mm10$transcript_id%in%lincRNA_mm10$transcript_id,'transcript_type']='lincRNA'


sRNA_mm10=mm10[mm10$transcript_type%in%'sRNA',]
misc_RNA_mm10=mm10[mm10$transcript_type%in%'misc_RNA',]
ribozyme_RNA_mm10=mm10[mm10$transcript_type%in%'ribozyme',]



# lncRNA_mm10_G=mm10[mm10$gene_biotype_ALL%in%'lncRNA',]


#######################################################
#### get  repeat regions and extar small RNA locations
newRmsk=1
if (newRmsk==1) {
  
  rmsk_GRCm38=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/repeatmasker/rmsk_GRCm38.txt"), header=T, sep="\t",stringsAsFactors = F,data.table=F)
  # unique(rmsk_GRCm38$repClass)
  rmsk_GRCm38all=rmsk_GRCm38
  
  # r1=rmsk_GRCm38[rmsk_GRCm38$repFamily%in%c('scRNA','snRNA','srpRNA','tRNA','rRNA','RNA'),]
  # r2=rmsk_GRCm38[rmsk_GRCm38$repClass%in%c('LINE','SINE','LTR','DNA','Satellite','Simple_repeat','Low_complexity','Other','Unknown'),]
  # rmsk_GRCm38=rbind(r1,r2)
  # 
  # write.table(rmsk_GRCm38,file=("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/repeatmasker/rmsk_GRCm38_subClip.txt"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
}


#################################################
###### RNA SUBTYPES FROM GENCODE AND REPEATMASKER
#################################################

## REPEATMASKER  ## repClass
# unique(rmsk_GRCm38all$repClass)
aa=c("LTR", "LINE", "SINE", "Simple_repeat", "DNA", "Satellite", "Low_complexity", "Other", "Unknown", "RC")
bb=c("SINE?", "LTR?", "DNA?", "RC?", "LINE?")  
# scRNA
# snRNA
# srpRNA
# tRNA
# rRNA
# RNA -- all RNA are 7SK (all 7SK are RNA)

YRNA_rmsk=rmsk_GRCm38all[(rmsk_GRCm38all$repFamily%in%'scRNA')&(grepl("HY",rmsk_GRCm38all$repName)==T),]; 
YRNA_rmsk=YRNA_rmsk[,c('genoName','genoStart','genoEnd','repName','swScore','strand')]
colnames(YRNA_rmsk)=c('chr','start','end','name','swScore','strand')
YRNA_rmsk$type='yRNA'
scRNA_rmsk=rmsk_GRCm38all[(rmsk_GRCm38all$repFamily%in%'scRNA')&(grepl("HY",rmsk_GRCm38all$repName)==F),]; 
scRNA_rmsk=scRNA_rmsk[,c('genoName','genoStart','genoEnd','repName','swScore','strand')]
colnames(scRNA_rmsk)=c('chr','start','end','name','swScore','strand')
scRNA_rmsk$type='scRNA'
snRNA_rmsk=rmsk_GRCm38all[rmsk_GRCm38all$repFamily%in%'snRNA',]; 
snRNA_rmsk=snRNA_rmsk[,c('genoName','genoStart','genoEnd','repName','swScore','strand')]
colnames(snRNA_rmsk)=c('chr','start','end','name','swScore','strand')
snRNA_rmsk$type='snRNA'
srpRNA_rmsk=rmsk_GRCm38all[rmsk_GRCm38all$repFamily%in%'srpRNA',]; 
srpRNA_rmsk=srpRNA_rmsk[,c('genoName','genoStart','genoEnd','repName','swScore','strand')]
colnames(srpRNA_rmsk)=c('chr','start','end','name','swScore','strand')
srpRNA_rmsk$type='srpRNA'
tRNA_rmsk=rmsk_GRCm38all[rmsk_GRCm38all$repFamily%in%'tRNA',]; 
tRNA_rmsk=tRNA_rmsk[,c('genoName','genoStart','genoEnd','repName','swScore','strand')]
colnames(tRNA_rmsk)=c('chr','start','end','name','swScore','strand')
tRNA_rmsk$type='tRNA' 
rRNA_rmsk=rmsk_GRCm38all[rmsk_GRCm38all$repFamily%in%'rRNA',]; 
rRNA_rmsk=rRNA_rmsk[,c('genoName','genoStart','genoEnd','repName','swScore','strand')]
colnames(rRNA_rmsk)=c('chr','start','end','name','swScore','strand')
rRNA_rmsk$type='rRNA'
SKRNA_rmsk=rmsk_GRCm38all[rmsk_GRCm38all$repFamily%in%'RNA',]; 
SKRNA_rmsk=SKRNA_rmsk[,c('genoName','genoStart','genoEnd','repName','swScore','strand')]
colnames(SKRNA_rmsk)=c('chr','start','end','name','swScore','strand')
SKRNA_rmsk$type='7SKRNA'


### old version
# RNA_rmsk_GRCm38=rmsk_GRCm38[rmsk_GRCm38$repClass%in%'RNA',]
# rRNA_rmsk_GRCm38=rmsk_GRCm38[rmsk_GRCm38$repClass%in%'rRNA',]
# tRNA_rmsk_GRCm38=rmsk_GRCm38[rmsk_GRCm38$repClass%in%'tRNA',]
# scRNA_rmsk_GRCm38=rmsk_GRCm38[rmsk_GRCm38$repClass%in%'scRNA',]
# snRNA_rmsk_GRCm38=rmsk_GRCm38[rmsk_GRCm38$repClass%in%'snRNA',]
# rmsk_GRCm38=rmsk_GRCm38[!rmsk_GRCm38$repClass%in%c('RNA','rRNA','tRNA','scRNA','snRNA'),]

###################################
## REPEATMASKER
## repFamily

# unique(rmsk_GRCm38all$repFamily)
# scRNA
# snRNA
# srpRNA
# tRNA
# rRNA
# RNA           
# MIR -- All MIR are SINES
miRNA_rmsk=rmsk_GRCm38all[rmsk_GRCm38all$repFamily%in%'MIR',]

# 
# Simple_repeat
# Alu
# Satellite     
# Low_complexity 
# Other
# LTR
# Unknown
# DNA

a=c("ID", "ERVK", "L1", "B4", "B2", "L2", "CR1", "Deu", "ERV1", "ERVL-MaLR", "ERVL", "hAT-Charlie", "TcMar-Tigger", "RTE-X", "hAT-Tip100", "Gypsy", "hAT-Blackjack", "RTE-BovB", "TcMar-Mariner", "TcMar-Tc2", "Y-chromosome", "Helitron", "PiggyBac", "MULE-MuDR", "hAT", "Dong-R4", "TcMar", "MuDR", "TcMar-Pogo", "centr")
b=c("ERVK?", "SINE?", "Gypsy?", "LTR?", "ERVL?", "TcMar?", "ERV1?", "DNA?", "Helitron?", "hAT?", "hAT-Tip100?", "PiggyBac?", "L1?", "Penelope?")

rmsk_GRCm38_LISI=rmsk_GRCm38[rmsk_GRCm38$repClass%in%c('LINE','SINE'),]
rmsk_GRCm38_LTR=rmsk_GRCm38[rmsk_GRCm38$repClass%in%'LTR',]
rmsk_GRCm38_DNA=rmsk_GRCm38[rmsk_GRCm38$repClass%in%'DNA',]
rmsk_GRCm38_sat=rmsk_GRCm38[rmsk_GRCm38$repClass%in%'Satellite',]

rmsk_GRCm38_SR=rmsk_GRCm38[rmsk_GRCm38$repClass%in%'Simple_repeat',]
rmsk_GRCm38_LC=rmsk_GRCm38[rmsk_GRCm38$repClass%in%'Low_complexity',]
rmsk_GRCm38_Other=rmsk_GRCm38[rmsk_GRCm38$repClass%in%'Other',]
rmsk_GRCm38_unknown=rmsk_GRCm38[rmsk_GRCm38$repClass%in%'Unknown',]
rmsk_GRCm38_LowComplx=rmsk_GRCm38[rmsk_GRCm38$repClass%in%'Low_complexity',]
# unique(rmsk_GRCm38$repClass)
# rmsk_GRCm38_repeat=rmsk_GRCm38[rmsk_GRCm38$repClass%in%c('Simple_repeat','Low_complexity','Other','Unknown'),]


#################################################
#### SOYEONG
#################################################

##### yRNA 
yRNA_sy=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mm10_YRNA.bed"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
colnames(yRNA_sy)=c('chr','start','end','name','swScore','strand')
yRNA_sy$type='yRNA'


##### sncRNA (mrp RNA, Rpph1, vault RNA)	mm10_snc.gtf
srpRNA_sy=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mm10_srpRNA.bed"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
colnames(srpRNA_sy)=c('chr','start','end','name','swScore','strand')
srpRNA_sy$type='srpRNA'

tRNA_sy=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mm10_tRNA.bed"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
colnames(tRNA_sy)=c('chr','start','end','name','swScore','strand')
tRNA_sy$type='tRNA'


##### Soyeong Files Bed - from Repeatmasker

##### scRNA (4.5S and BC1 RNA) (Gencode)
scRNA_sy=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mm10_scRNA.bed"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
# colnames(scRNA_sy)=c('chr','start','end','name','swScore','strand')
# scRNA_sy$type='scRNA'

##### 7SK RNA	mm10_7SKRNA.bed
SKRNA_sy=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mm10_7SKRNA.bed"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
# colnames(SKRNA_sy)=c('chr','start','end','name','swScore','strand')
# SKRNA_sy$type='7SKRNA'


###### rRNA pseudogenes and 5S	mm10_rRNA.bed and mouse_gencode_rRNA.gtf.1
rRNA_sy=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mm10_rRNA.bed"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
# colnames(rRNA_sy)=c('chr','start','end','name','swScore','strand')
# rRNA_sy$type='rRNA'

###### rRNA BK00964.3 'https://www.ncbi.nlm.nih.gov/nuccore/NR_046233.1'
rRNA_BK00964=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/rCDNA_nacentproject/BK000964.3_TPA_rRNA_repeats2.bed.txt"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
rRNA_BK00964=rRNA_BK00964[,c('V1','V2','V3','V4','V5','V6')]
colnames(rRNA_BK00964)=c('chr','start','end','name','swScore','strand')
rRNA_BK00964$strand="*"
rRNA_BK00964$type='rRNA'

##### Introns		mm10_intron.bed
introns_SY=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mm10_intron.bed"), header=F, sep="\t",stringsAsFactors = F,data.table=F)

######## 
# Anno_bed_comb=rbind(yRNA_SY,tRNA_SY,srpRNA_SY,SKRNA_SY,scRNA_SY,rRNA_SY)
###############################################

##### Soyeong Files GTF - from Gencode

##### sncRNA (mrp RNA, Rpph1, vault RNA)	mm10_snc.gtf
sncRNA_sy=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mm10_snc.gtf"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
sncRNA_sy=separate(sncRNA_sy,col = "V9",into = c("gene","trans",'x'),sep = ";");sncRNA_sy$gene=gsub("gene_id ","",sncRNA_sy$gene);sncRNA_sy$trans=gsub("transcript_id ","",sncRNA_sy$trans)
sncRNA_sy=as.data.frame(sncRNA_sy)
sncRNA_sy=sncRNA_sy[,((colnames(sncRNA_sy)%in%'x')==F)]

colnames(sncRNA_sy)=c('chr','source','gloc','start','end','V6','strand','V8','Gene_Ensemble','Transcript_id')
sncRNA_sy$type='sncRNA'
colnames(sncRNA_sy)[colnames(sncRNA_sy)%in%'Gene_Ensemble']='name'


###### snRNA (Gencode)
# snRNA_sy=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mouse_gencode_snRNA.gtf.1"), header=F, sep="\t",stringsAsFactors = F,data.table=F )
snRNA_sy=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mouse_gencode_snRNA.gtf.1"), header=F, sep="\t",stringsAsFactors = F,data.table=F,quote="" )
colnames(snRNA_sy)=c('chr','source','gloc','start','end','V6','strand','V8','V9','Gene_Ensemble','V11','V12','V13','Transcript_id','V15','version')
snRNA_sy=snRNA_sy[snRNA_sy$gloc%in%'gene',]


###### snoRNA (Gencode)
# snoRNA_sy=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mouse_gencode_sno.gtf.1"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
snoRNA_sy=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mouse_gencode_sno.gtf.1"), header=F, sep="\t",stringsAsFactors = F,data.table=F,quote="")
colnames(snoRNA_sy)=c('chr','source','gloc','start','end','V6','strand','V8','V9','Gene_Ensemble','V11','V12','V13','Transcript_id','V15','version')
snoRNA_sy=snoRNA_sy[snoRNA_sy$gloc%in%'gene',]


###### miRNA		mouse_gencode_miRNA.gtf.1 (Gencode)
miRNA_sy=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mouse_gencode_miRNA.gtf.1"), header=F, sep="\t",stringsAsFactors = F,data.table=F,quote="")
colnames(miRNA_sy)=c('chr','source','gloc','start','end','V6','strand','V8','V9','Gene_Ensemble','V11','V12','V13','Transcript_id','V15','version')
miRNA_sy=miRNA_sy[miRNA_sy$gloc%in%'gene',]


###### rRNA pseudogenes and 5S	mm10_rRNA.bed and mouse_gencode_rRNA.gtf.1 (Gencode)
rRNA_gtf=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mouse_gencode_rRNA.gtf.1"), header=F, sep="\t",stringsAsFactors = F,data.table=F,quote="")
colnames(rRNA_gtf)=c('chr','source','gloc','start','end','V6','strand','V8','V9','Gene_Ensemble','V11','V12','V13','Transcript_id','V15','version')
rRNA_gtf=rRNA_gtf[rRNA_gtf$gloc%in%'gene',]
## an examleof rRNA and rRNA_gtf chrY:991632-991664



###### linc RNA	mouse_gencode_linc.gtf.1 (Repeatmasker)
lincRNA=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mouse_gencode_linc.gtf.1"), header=F, sep="\t",stringsAsFactors = F,data.table=F,quote="")
colnames(lincRNA)=c('chr','source','gloc','start','end','V6','strand','V8','V9','Gene_Ensemble','V11','V12','V13','Transcript_id','V15','version')
lincRNA=lincRNA[lincRNA$gloc%in%'gene',]


###### mRNA	mm10_refFlat.gtf
# mRNA=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mm10_refFlat.gtf"), header=F, sep="\t",stringsAsFactors = F,data.table=F,quote="")
# mRNA=separate(mRNA,col = "V9",into = c("gene","trans","x"),sep = ";");mRNA$gene=gsub("gene_id ","",mRNA$gene);mRNA$trans=gsub("transcript_id ","",mRNA$trans)
# mRNA$class=substring(separate(mRNA,col="gene",into = c("num","class"),sep = "-")$class,1,3)
# mRNA=separate(data=mRNA,col=trans,into=c('x','trans'),remove=T,sep = "\\.")


SY_LTR=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mm10_LTR.bed"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
SY_DNA=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mm10_DNA.bed"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
SY_sat=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mm10_sat.bed"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
SY_SR=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mm10_simple.bed"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
SY_LC=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mm10_LC.bed"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
SY_other=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mm10_other.bed"), header=F, sep="\t",stringsAsFactors = F,data.table=F)
SY_unknown=fread(paste0("/Users/homanpj/OneDrive - National Institutes of Health/Resources/ref/mm10/from_Soyeong/mm10_annotation/mm10_unknown.bed"), header=F, sep="\t",stringsAsFactors = F,data.table=F)

# rmsk_GRCm38_repeat=rmsk_GRCm38[rmsk_GRCm38$repClass%in%c('Simple_repeat','Low_complexity','Other','Unknown'),]
# repeats_GRCm38=fread(paste0("./rbl.wolin.parclip/inst/extdata/rmsk_GRCm38_LSlc.txt"), header=T, sep="\t",stringsAsFactors = F,data.table=F)




###########################################################################################################################################################################################
###########################################################################################################################################################################################

#### Create CLASSIFICATION table 

###########################################################################################################################################################################################
###########################################################################################################################################################################################


rnames=c('45S','chrM','yRNA','snRNA','snoRNA','srpRNA','tRNA','7SK RNA','scRNA','sncRNA','miRNA','rRNA_gencode','rRNA_rmsk','rRNA_DNA','lncRNA','lincRNA','mRNA','LINE','SINE','LTR','DNA','Satalites','Simple Repeats','Low Complexity','Other repeats','Unknown repeats','Introns')
cnames=c('SY_count','Gencode','Repeatmasker','Comb','contents','source','Description','notes','Annotation_file')
Classification=as.data.frame(matrix(nrow=length(rnames),ncol = length(cnames)))
rownames(Classification)=rnames
colnames(Classification)=cnames

Classification['yRNA','SY_count']=nrow(yRNA_sy)
Classification['snRNA','SY_count']=nrow(snRNA_sy)
Classification['snoRNA','SY_count']=nrow(snoRNA_sy)
Classification['srpRNA','SY_count']=nrow(srpRNA_sy)
Classification['tRNA','SY_count']=nrow(tRNA_sy)
Classification['7SK RNA','SY_count']=nrow(SKRNA_sy)
Classification['scRNA','SY_count']=nrow(scRNA_sy)
Classification['sncRNA','SY_count']=nrow(sncRNA_sy)
Classification['miRNA','SY_count']=nrow(miRNA_sy)
Classification['rRNA','SY_count']=nrow(rRNA_sy)+nrow(rRNA_gtf)
# Classification['lincRNA','SY_count']=nrow(lincRNA_mm10)
Classification['LTR','SY_count']=nrow(SY_LTR)
Classification['DNA','SY_count']=nrow(SY_DNA)
Classification['Satalites','SY_count']=nrow(SY_sat)
Classification['Simple Repeats','SY_count']=nrow(SY_SR)
Classification['Low Complexity','SY_count']=nrow(SY_LC)
Classification['Other repeats','SY_count']=nrow(SY_other)
Classification['Unknown repeats','SY_count']=nrow(SY_unknown)

Classification['snRNA','Gencode']=nrow(snRNA_mm10)
Classification['snoRNA','Gencode']=nrow(snoRNA_mm10)
Classification['scRNA','Gencode']=nrow(scRNA_mm10)
Classification['miRNA','Gencode']=nrow(miRNA_mm10)
Classification['rRNA_gencode','Gencode']=nrow(rRNA_mm10)
Classification['lncRNA','Gencode']=nrow(lncRNA_mm10)


Classification['yRNA','Repeatmasker']=nrow(YRNA_rmsk)
Classification['snRNA','Repeatmasker']=nrow(snRNA_rmsk)
Classification['srpRNA','Repeatmasker']=nrow(srpRNA_rmsk)
Classification['tRNA','Repeatmasker']=nrow(tRNA_rmsk)
Classification['scRNA','Repeatmasker']=nrow(scRNA_rmsk)
Classification['miRNA','Repeatmasker']=nrow(miRNA_rmsk)
Classification['rRNA_rmsk','Repeatmasker']=nrow(rRNA_rmsk)
Classification['7SK RNA','Repeatmasker']=nrow(SKRNA_rmsk)

Classification['LTR','Repeatmasker']=nrow(rmsk_GRCm38_LTR)
Classification['DNA','Repeatmasker']=nrow(rmsk_GRCm38_DNA)
Classification['Satalites','Repeatmasker']=nrow(rmsk_GRCm38_sat)
Classification['Simple Repeats','Repeatmasker']=nrow(rmsk_GRCm38_SR)
Classification['Low Complexity','Repeatmasker']=nrow(rmsk_GRCm38_LC)
Classification['Other repeats','Repeatmasker']=nrow(rmsk_GRCm38_Other)
Classification['Unknown repeats','Repeatmasker']=nrow(rmsk_GRCm38_unknown)

Classification$Comb=rowSums(Classification[,c('Gencode','Repeatmasker')])

Classification['yRNA','notes']='use Repeatmasker - subset of scRNA'
Classification['snRNA','notes']='Use Gencode'
Classification['snoRNA','notes']='Use Gencode'
Classification['srpRNA','notes']='Use Repeatmasker - can be (7SL, 6S, or 4.5S RNA) 4.5S is coveredin scRNA'
Classification['tRNA','notes']='use Soyeong'
Classification['7SK RNA','notes']='Use Repeatmakser'
Classification['scRNA','notes']='use Repeatmakser - removed yRNA'
Classification['sncRNA','notes']='use Soyeong'
Classification['miRNA','notes']='use Gencode'
Classification['rRNA','notes']='Use Gencode + Repeatmasker'
Classification['lincRNA','notes']='older versions of Gencode use linc'
Classification['lncRNA','notes']='Gencode'

Classification['yRNA','source']='Repeatmasker'
Classification['snRNA','source']='Gencode_VM23'
Classification['snoRNA','source']='Gencode_VM23'
Classification['srpRNA','source']='Repeatmasker'
Classification['tRNA','source']='GtRNAdb'
Classification['7SK RNA','source']='Repeatmakser'
Classification['scRNA','source']='Repeatmakser'
Classification['sncRNA','source']=''
Classification['miRNA','source']='Gencode_VM23'
Classification['rRNA_gencode','source']='Gencode_VM23'
Classification['rRNA_rmsk','source']='Repeatmasker'
Classification['rRNA_DNA','source']='BK000964.3'
Classification['lncRNA','source']='Gencode_VM23'

Classification['yRNA','Description']=''
Classification['snRNA','Description']='small nuclear RNA : Small RNA molecules that are found in the cell nucleus and are involved in the processing of pre messenger RNAs'
Classification['snoRNA','Description']='Small nucleolar RNAs : Small RNA molecules that are found in the cell nucleolus and are involved in the post-transcriptional modification of other RNAs'
Classification['srpRNA','Description']='signal recognition particle RNA' 
Classification['tRNA','Description']='transfer RNA, which acts as an adaptor molecule for translation of mRNA.'
Classification['7SK RNA','Description']='subset of small nuclear RNA and part of the small nuclear ribonucleoprotein complex (snRNP)'
Classification['scRNA','Description']='Small cytoplasmic RNA' 
Classification['sncRNA','Description']='small non-coading RNA' 
Classification['miRNA','Description']='Micro RNA : A small RNA (~22bp) that silences the expression of target mRNA'
Classification['rRNA_gencode','Description']='Ribosomal RNA'
Classification['rRNA_rmsk','Description']='Ribosomal RNA'
Classification['rRNA_DNA','Description']='Ribosomal DNA'
Classification['lncRNA','Description']='Generic long non-coding RNA biotype'
Classification['lincRNA','Description']='long non-coding RNA biotype with Intronic + Exonic Regions'

Classification['yRNA','contents']=paste0(unique(YRNA_rmsk$name),collapse = ', ')
Classification['snRNA','contents']='U1,U2,U5,U6,U7,U11,U12 and various predicted genes' ;#paste0(unique(snRNA_mm10$transcript_type),collapse = ', ')
Classification['snoRNA','contents']=paste0('Various ',paste0(unique(snoRNA_mm10$transcript_type),collapse = ', '))
Classification['srpRNA','contents']=paste0(unique(srpRNA_rmsk$name),collapse = ', ')
Classification['tRNA','contents']=paste0(unique(tRNA_sy$type),collapse = ', ')
Classification['7SK RNA','contents']=paste0(unique(SKRNA_rmsk$name),collapse = ', ')
Classification['scRNA','contents']=paste0(unique(scRNA_rmsk$name),collapse = ', ')
Classification['sncRNA','contents']=paste0('3 annotations: ',paste0(unique(sncRNA_sy$name),collapse = ', '))
Classification['miRNA','contents']=paste0(unique(miRNA_mm10$transcript_type),collapse = ', ')
Classification['rRNA_gencode','contents']='5S, 5.8s, predicted gene';
Classification['rRNA_rmsk','contents']=paste0(unique(rRNA_rmsk[order(rRNA_rmsk$name),'name']),collapse = ', ')
Classification['rRNA_DNA','contents']=paste0(unique(rRNA_BK00964$name),collapse = ', ')
Classification['lncRNA','contents']=paste0('Various ',paste0(unique(lncRNA_mm10$transcript_type),collapse = ', '))
Classification['lincRNA','contents']=paste0('Various ',paste0(unique(lincRNA_mm10$transcript_type),collapse = ', '))

Classification['yRNA','Annotation_file']='yRNA.bed'
Classification['snRNA','Annotation_file']='snRNA.bed'
Classification['snoRNA','Annotation_file']='snoRNA.bed'
Classification['srpRNA','Annotation_file']='srpRNA.bed'
Classification['tRNA','Annotation_file']='tRNA.bed'
Classification['7SK RNA','Annotation_file']='SKRNA.bed'
Classification['scRNA','Annotation_file']='scRNA.bed'
Classification['sncRNA','Annotation_file']='sncRNA.bed'
Classification['miRNA','Annotation_file']='miRNA.bed'
Classification['rRNA_gencode','Annotation_file']='rRNA_gencode.bed'
Classification['rRNA_rmsk','Annotation_file']='rRNA_rmsk.bed'
Classification['rRNA_DNA','Annotation_file']='rRNA_BK00964.bed'
Classification['lncRNA','Annotation_file']='lncRNA.bed'
Classification['lncRNA','Annotation_file']='lincRNA.bed'



# rRNA_gencode types unique(rRNA_mm10[order(rRNA_mm10$external_gene_name),'external_gene_name'])
# lncRMA_gencode types unique(lncRNA_mm10[order(lncRNA_mm10$external_gene_name),'external_gene_name'])


### add annotations to Gencode for:
# yRNA_sy
# srpRNA_rmsk
# tRNA_sy
# SKRNA_rmsk
# scRNA_rmsk
# sncRNA_sy
# rRNA_rmsk

if (WriteClassTable==T) {
sycol=c('chr','start','end','name','swScore','strand')
mm10col=c('chr','start','end','transcript_name','score','strand')
rmskcol=c('chr','start','end','name','swScore','strand')
write.table(YRNA_rmsk[,rmskcol],file=("/Users/homanpj/OneDrive - National Institutes of Health/RBL/Wolin/R0_KO/mESC_clip/annotation/yRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
write.table(snRNA_mm10[,mm10col],file=("/Users/homanpj/OneDrive - National Institutes of Health/RBL/Wolin/R0_KO/mESC_clip/annotation/snRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
write.table(snoRNA_mm10[,mm10col],file=("/Users/homanpj/OneDrive - National Institutes of Health/RBL/Wolin/R0_KO/mESC_clip/annotation/snoRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
write.table(srpRNA_rmsk[,rmskcol],file=("/Users/homanpj/OneDrive - National Institutes of Health/RBL/Wolin/R0_KO/mESC_clip/annotation/srpRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
write.table(tRNA_sy[,rmskcol],file=("/Users/homanpj/OneDrive - National Institutes of Health/RBL/Wolin/R0_KO/mESC_clip/annotation/tRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
write.table(SKRNA_rmsk[,rmskcol],file=("/Users/homanpj/OneDrive - National Institutes of Health/RBL/Wolin/R0_KO/mESC_clip/annotation/SKRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
write.table(scRNA_rmsk[,rmskcol],file=("/Users/homanpj/OneDrive - National Institutes of Health/RBL/Wolin/R0_KO/mESC_clip/annotation/scRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
write.table(sncRNA_sy[,c('chr','start','end','name','V6','strand')],file=("/Users/homanpj/OneDrive - National Institutes of Health/RBL/Wolin/R0_KO/mESC_clip/annotation/sncRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
write.table(miRNA_mm10[,mm10col],file=("/Users/homanpj/OneDrive - National Institutes of Health/RBL/Wolin/R0_KO/mESC_clip/annotation/mirNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
write.table(rRNA_rmsk[,rmskcol],file=("/Users/homanpj/OneDrive - National Institutes of Health/RBL/Wolin/R0_KO/mESC_clip/annotation/rRNA_rmsk.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
write.table(rRNA_mm10[,mm10col],file=("/Users/homanpj/OneDrive - National Institutes of Health/RBL/Wolin/R0_KO/mESC_clip/annotation/rRNA_gencode.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
write.table(rRNA_BK00964[,rmskcol],file=("/Users/homanpj/OneDrive - National Institutes of Health/RBL/Wolin/R0_KO/mESC_clip/annotation/rRNA_BK00964.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)
write.table(lncRNA_mm10[,mm10col],file=("/Users/homanpj/OneDrive - National Institutes of Health/RBL/Wolin/R0_KO/mESC_clip/annotation/lncRNA.bed"), sep = "\t", row.names = F, col.names = F, append = F, quote= FALSE)


rnames_ncRNA=c('yRNA','snRNA','snoRNA','srpRNA','tRNA','7SK RNA','scRNA','sncRNA','miRNA','rRNA_gencode','rRNA_rmsk','rRNA_DNA','lncRNA','lincRNA')
cnames_ncRNA=c('source','contents','Description','Annotation_file')
Classification_ncRNA=Classification[rnames_ncRNA,cnames_ncRNA]

  write.table(Classification,file=paste0(outdir,"Annotations.txt"), sep = "\t", row.names = T, col.names = T, append = F, quote= FALSE)
  write.table(Classification_ncRNA,file=paste0(outdir,"ncRNA_Annotations.txt"), sep = "\t", row.names = T, col.names = T, append = F, quote= FALSE)
  
}


###########################################################################################################################################################################################
###########################################################################################################################################################################################
## Identify CLIP peak Gene Location

### Identify CLIP peak Host gene   

# Using GTF file from GENCODE v24  
# Peaks were annotated with overlapping gene  
###########################################################################################################################################################################################
###########################################################################################################################################################################################

##################################################################################################
#### FUNCTIONS
##################################################################################################


mm10_anno=function(Annotable,peaksTable){
  
# Annotable=mm10
# peaksTable=peaks
ColumnName=c('ensembl_gene_id','external_gene_name','gene_type','gene_biotype_ALL')

s=peaksTable
s.GR <- GRanges(seqnames = as.character(s$chr), ranges=IRanges(start = as.numeric(s$start), end = as.numeric(s$end)),strand = s$strand,ID=s$ID )


anno.GR <- GRanges(seqnames = as.character(Annotable$chr), ranges=IRanges(start = as.numeric(Annotable$start), end = as.numeric(Annotable$end)),strand = as.character(Annotable$strand),
                   ensembl_gene_id=Annotable$ensembl_gene_id,
                   external_gene_name=Annotable$external_gene_name,
                   gene_type=Annotable$gene_type,
                   gene_biotype_ALL=Annotable$gene_biotype_ALL
)

q =s.GR
s=anno.GR
xo=as.data.frame(GenomicRanges::findOverlaps(q,s,type = "any",ignore.strand=F))

qh=as.data.frame(q[xo$queryHits],row.names = NULL)
sh=as.data.frame(s[xo$subjectHits],row.names = NULL);colnames(sh)=paste0(colnames(sh),"_anno")
rmskinfo=cbind(qh,sh)
rmskinfo=rmskinfo[,c('ID','ensembl_gene_id_anno','external_gene_name_anno','gene_type_anno','gene_biotype_ALL_anno')];
colnames(rmskinfo)[colnames(rmskinfo)%in%c('ensembl_gene_id_anno','external_gene_name_anno','gene_type_anno','gene_biotype_ALL_anno')]= c(ColumnName)      

peaksTable=merge(peaksTable[,!colnames(peaksTable)%in%c(ColumnName)],rmskinfo,by='ID',all.x=T)

dup=unique(peaksTable[duplicated(peaksTable$ID),'ID'])
peaksTable_single=peaksTable[!(peaksTable$ID%in%dup),]
peaksTable_double=peaksTable[(peaksTable$ID%in%dup),]

u=unique(peaksTable_double$ID)
peaksTable_colapsed=as.data.frame(matrix(nrow=length(u),ncol=ncol(peaksTable_double)));
colnames(peaksTable_colapsed)=colnames(peaksTable_double)

for(x in 1:length(u)){
  p=u[x]
  pam=peaksTable_double[peaksTable_double$ID%in%p,]
  peaksTable_colapsed[x,colnames(peaksTable_double)%in%c(ColumnName)==F]=(pam[1,colnames(peaksTable_double)%in%c(ColumnName)==F,drop=T])
  
  peaksTable_colapsed[x,ColumnName[1]]=paste(sort(unique((pam[,ColumnName[1]]))), collapse =",")
  peaksTable_colapsed[x,ColumnName[2]]=paste(sort(unique((pam[,ColumnName[2]]))), collapse =",")
  peaksTable_colapsed[x,ColumnName[3]]=paste(sort(unique((pam[,ColumnName[3]]))), collapse =",")
  peaksTable_colapsed[x,ColumnName[4]]=paste(sort(unique((pam[,ColumnName[4]]))), collapse =",")
  
}

peaksTable_colapsed=rbind(peaksTable_colapsed,peaksTable_single)
peaksTable_colapsed=peaksTable_colapsed[!is.na(peaksTable_colapsed$ID),]


# peaksTable_colapsed[(((peaksTable_colapsed$gene_type>0)==F)&(is.na(peaksTable_colapsed$gene_type)==F)),]

peaksTable_colapsed[(((peaksTable_colapsed$gene_type>0)==F)&(is.na(peaksTable_colapsed$gene_type)==F)),c('ensembl_gene_id','external_gene_name','gene_type','gene_biotype_ALL')]=NA
peaksTable=peaksTable_colapsed
return(peaksTable) 
}


##################################################################################################

bam_anno=function(ColumnName,Annotable,peaksTable){
  # peakstest=PeaksdataOut
  # ColumnName='test'
  # Annotable=tRNA_sy
  # peaksTable=peakstest
  
  s=peaksTable
  s=separate(s,ID,into=c('chr','start'),sep=":",remove=F)
  s=separate(s,start,into=c('start','end'),sep="-",remove=F)
  
  s.GR <- GRanges(seqnames = as.character(s$chr), ranges=IRanges(start = as.numeric(s$start), end = as.numeric(s$end)),strand = s$strand,ID=s$ID )
  
  anno.GR <- GRanges(seqnames = as.character(Annotable$chr), ranges=IRanges(start = as.numeric(Annotable$start), end = as.numeric(Annotable$end)),strand = Annotable$strand,name=Annotable$name,type=Annotable$type)
  
  q =s.GR
  s=anno.GR
  xo=as.data.frame(GenomicRanges::findOverlaps(q,s,type = "any",ignore.strand=F))
  
  qh=as.data.frame(q[xo$queryHits],row.names = NULL)
  sh=as.data.frame(s[xo$subjectHits],row.names = NULL);colnames(sh)=paste0(colnames(sh),"_anno")
  rmskinfo=cbind(qh,sh)
  rmskinfo=rmskinfo[,c('ID','type_anno','name_anno')];
  colnames(rmskinfo)[colnames(rmskinfo)%in%c('type_anno','name_anno')]= c(paste0('type_',ColumnName),paste0('name_',ColumnName))      
  
  peaksTable=merge(peaksTable[,!colnames(peaksTable)%in%c(paste0('type_',ColumnName),paste0('name_',ColumnName))],rmskinfo,by='ID',all.x=T)
  
  dup=unique(peaksTable[duplicated(peaksTable$ID),'ID'])
  peaksTable_single=peaksTable[!(peaksTable$ID%in%dup),]
  peaksTable_double=peaksTable[(peaksTable$ID%in%dup),]
  
  u=unique(peaksTable_double$ID)
  peaksTable_colapsed=as.data.frame(matrix(nrow=length(u),ncol=ncol(peaksTable_double)));
  colnames(peaksTable_colapsed)=colnames(peaksTable_double)
  
  for(x in 1:length(u)){
    p=u[x]
    pam=peaksTable_double[peaksTable_double$ID%in%p,]
    peaksTable_colapsed[x,colnames(peaksTable_double)%in%c(ColumnName)==F]=(pam[1,colnames(peaksTable_double)%in%c(paste0('type_',ColumnName),paste0('name_',ColumnName))==F,drop=T])
    
    peaksTable_colapsed[x,paste0('name_',ColumnName)]=paste(sort(unique((pam[,paste0('name_',ColumnName)]))), collapse =",")
    peaksTable_colapsed[x,paste0('type_',ColumnName)]=paste(sort(unique((pam[,paste0('type_',ColumnName)]))), collapse =",")
  }
  
  peaksTable_colapsed=rbind(peaksTable_colapsed,peaksTable_single)
  peaksTable_colapsed=peaksTable_colapsed[!is.na(peaksTable_colapsed$ID),]
  
  
  peaksTable_colapsed[is.na(peaksTable_colapsed[,paste0('name_',ColumnName)]),paste0('name_',ColumnName)]=NA
  peaksTable_colapsed[peaksTable_colapsed[,paste0('name_',ColumnName)]%in%"",paste0('name_',ColumnName)]=NA
  
  
  peaksTable_colapsed[is.na(peaksTable_colapsed[,paste0('type_',ColumnName)]),paste0('type_',ColumnName)]=NA
  peaksTable_colapsed[peaksTable_colapsed[,paste0('type_',ColumnName)]%in%"",paste0('type_',ColumnName)]=NA
  
  peaksTable=peaksTable_colapsed
  
  return(peaksTable) 
}


##################################################################################################

rpmsk_anno=function(ColumnName,Annotable,peaksTable){
  
  # ColumnName='test'
  # AnnoTable=rmsk_GRCm38_LISI
  # peaksTable=PeaksdataOut
  
  s=peaksTable
  s=separate(s,ID,into=c('chr','start'),sep=":",remove=F)
  s=separate(s,start,into=c('start','end'),sep="-",remove=F)
  
  s.GR <- GRanges(seqnames = as.character(s$chr), ranges=IRanges(start = as.numeric(s$start), end = as.numeric(s$end)),strand = s$strand,ID=s$ID )
  
  anno.GR <- GRanges(seqnames = as.character(Annotable$genoName), ranges=IRanges(start = as.numeric(Annotable$genoStart), end = as.numeric(Annotable$genoEnd)),strand = Annotable$strand,repClass=Annotable$repClass,repName=Annotable$repName )
  
  q =s.GR
  s=anno.GR
  xo=as.data.frame(GenomicRanges::findOverlaps(q,s,type = "any",ignore.strand=F))
  
  qh=as.data.frame(q[xo$queryHits],row.names = NULL)
  sh=as.data.frame(s[xo$subjectHits],row.names = NULL);colnames(sh)=paste0(colnames(sh),"_repeat")
  rmskinfo=cbind(qh,sh)
  rmskinfo=rmskinfo[,c('ID','repClass_repeat')];colnames(rmskinfo)[colnames(rmskinfo)%in%'repClass_repeat']= ColumnName      
  
  peaksTable=merge(peaksTable[,!colnames(peaksTable)%in%ColumnName],rmskinfo,by='ID',all.x=T)
  
  dup=unique(peaksTable[duplicated(peaksTable$ID),'ID'])
  peaksTable_single=peaksTable[!(peaksTable$ID%in%dup),]
  peaksTable_double=peaksTable[(peaksTable$ID%in%dup),]
  
  u=unique(peaksTable_double$ID)
  peaksTable_colapsed=as.data.frame(matrix(nrow=length(u),ncol=ncol(peaksTable_double)));
  colnames(peaksTable_colapsed)=colnames(peaksTable_double)
  for(x in 1:length(u)){
    p=u[x]
    # pd=as.numeric(peaks2[x,'V4'])
    
    pam=peaksTable_double[peaksTable_double$ID%in%p,]
    peaksTable_colapsed[x,colnames(peaksTable_double)%in%c(ColumnName)==F]=(pam[1,colnames(peaksTable_double)%in%c(ColumnName)==F,drop=T])
    
    peaksTable_colapsed[x,ColumnName]=paste(sort(unique((pam[,ColumnName]))), collapse =" | ")
  }
  
  peaksTable_colapsed=rbind(peaksTable_colapsed,peaksTable_single)
  # peaksTable_colapsed=peaksTable_colapsed[duplicated(peaksTable_colapsed[,colnames(peaksTable_colapsed)%in%c(ColumnName)==F])==F,]
  
  peaksTable_colapsed=peaksTable_colapsed[!is.na(peaksTable_colapsed$ID),]
  
  peaksTable_colapsed[is.na(peaksTable_colapsed[,ColumnName]),ColumnName]=NA
  peaksTable_colapsed[peaksTable_colapsed[,ColumnName]%in%"",ColumnName]=NA
  
  peaksTable=peaksTable_colapsed
  
  return(peaksTable) 
}


for (xopp in 1:2) {
  
  if (xopp==1) {peaks=peaks;nmeprfix='Same_'}
  if (xopp==2) {peaks=peaks_oppo;nmeprfix='Oppo_'}
  

##################################################################################################
# ANNOTATE
##################################################################################################

PeaksdataOut=mm10_anno(mm10,peaks)

###### idnividual columns
### add annotations for:
# yRNA_sy or yRNA_rmsk
# srpRNA_rmsk
# tRNA_sy
# SKRNA_rmsk
# scRNA_rmsk
# sncRNA_sy
# rRNA_rmsk
# PeaksdataOut=bam_anno('yRNA',yRNA_sy,PeaksdataOut)
# PeaksdataOut=bam_anno('srpRNA',srpRNA_rmsk,PeaksdataOut)
# PeaksdataOut=bam_anno('tRNA',tRNA_sy,PeaksdataOut)
# PeaksdataOut=bam_anno('SKRNA',SKRNA_rmsk,PeaksdataOut)
# PeaksdataOut=bam_anno('scRNA',scRNA_rmsk,PeaksdataOut)
# PeaksdataOut=gtf_anno('sncRNA',sncRNA_sy,PeaksdataOut)
# PeaksdataOut=bam_anno('rRNA',rRNA_rmsk,PeaksdataOut)

### Add Column that commbines all additional annotations
annocol=c('chr','start','end','strand','type','name')
Anno_RNA_comb=rbind(YRNA_rmsk[,annocol],srpRNA_rmsk[,annocol],tRNA_sy[,annocol],SKRNA_rmsk[,annocol],scRNA_rmsk[,annocol],sncRNA_sy[,annocol],rRNA_rmsk[,annocol],rRNA_BK00964[,annocol])

PeaksdataOut=bam_anno('RNA_anno',Anno_RNA_comb,PeaksdataOut)


##################################################################################################
#### CLEAN UP RNA TYPE
##################################################################################################

p=PeaksdataOut

##########################################################################################
## change lincRNA,rRNA to rRNA only
## do not change the name if there is actually a lincRNA name separate from the rRNA grepl(',',p$name_RNA_anno)==F)
p[(grepl(',',p$name_RNA_anno)==F)&(p$type_RNA_anno%in%'lincRNA,rRNA'),'type_RNA_anno']='rRNA'

##########################################################################################
## change all psueudogene classes to pseudogene
p[grep('pseudogene',p$gene_biotype_ALL),'gene_biotype_ALL']='pseudogene'

##########################################################################################
## change all sc RNA to yRNA names
p$type_RNA_anno=gsub('scRNA,yRNA','yRNA',p$type_RNA_anno)

PeaksdataOut=p


#########################
## Combine peak RNAtype info
########################
p=PeaksdataOut

p$type_simple_comb=NA
p$type_comb=NA
p$gene_name_comb=NA

#######################################################
dup=p[(is.na(p$gene_biotype_ALL)==F)|(is.na(p$type_RNA_anno)==F),'ID']
peaksTable_single=p[!(p$ID%in%dup),]
peaksTable_double=p[(p$ID%in%dup),]

u=unique(peaksTable_double$ID)
peaksTable_colapsed=as.data.frame(matrix(nrow=length(u),ncol=ncol(peaksTable_double)));
colnames(peaksTable_colapsed)=colnames(peaksTable_double)


for(x in 1:length(u)){
  p=u[x]
  pam=peaksTable_double[peaksTable_double$ID%in%p,]
  
  
  ###############################################
  ###############################################
  ### Comb Annotation
  
  pam_1=pam$gene_biotype_ALL ### Biotype from Gencode 
  pam_1=as.data.frame(strsplit(pam_1,','));colnames(pam_1)='a';pam_1$a=as.character(pam_1$a)
  pam_2=pam$type_RNA_anno ### Bitype from Annotation
  pam_2=as.data.frame(strsplit(pam_2,','));colnames(pam_2)='a';pam_2$a=as.character(pam_2$a)
  
  ##################################    
  
  ### remove misc_RNA catagory : these are covered by additional annotations (mostly yRNA)
  if (grepl('misc_RNA',pam_1)) {pam_1[pam_1$a%in%'misc_RNA','a']=NA}
  
  ### lincRNA are from an older Gencode version (VM18 or older) so don't use    
  if (grepl('lncRNA',pam_1)&grepl('lincRNA',pam_2)) {pam_1[pam_2$a%in%'lincRNA',]='lncRNA'}
  
  ### change ribozyme (gencode) to RNA type from additional anno
  if (grepl('ribozyme',pam_1)) {pam_1[pam_1$a%in%'ribozyme',]=unique(pam_2$a)}
  
  ##################################    
  ## combine all rna types subtypes into ncRNA
  pam_c=rbind(pam_1,pam_2)
  pam_c=pam_c[!is.na(pam_c$a),,drop=F]
  
  #### Annotation from Gencode : unique(sort(mm10$gene_biotype_ALL))
  # "miRNA","misc_RNA","snRNA","snoRNA","ribozyme","lncRNA","scRNA","sRNA","scaRNA","rRNA"
  pam_c2=pam_c
  pam_c2$a=gsub('miRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('miscRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('misc_RNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('piRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('rRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('siRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('snRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('snoRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('ribozyme','ncRNA',pam_c2$a)
  pam_c2$a=gsub('lncRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('lincRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('scRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('sRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('scaRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('rRNA','ncRNA',pam_c2$a)
  
  #### extra annotation unique(Anno_RNA_comb$type )
  pam_c2$a=gsub('yRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('srpRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('tRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('7SKRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('scRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('sncRNA','ncRNA',pam_c2$a)
  pam_c2$a=gsub('rRNA','ncRNA',pam_c2$a)
  
  if (length(pam_c$a)>1) {pam_c=pam_c[order(pam_c$a),,drop=F]}
  if (length(pam_c2$a)>1) {pam_c2=pam_c2[order(pam_c2$a),,drop=F]}
  
  pam_c=paste(unique(pam_c$a),collapse = ',')
  pam_c2=paste(unique(pam_c2$a),collapse = ',')
  
  
  ########################################################################################################################################
  ########################################################################################################################################
  ### Comb GeneName 
  
  pam_G1=pam$external_gene_name ### gene name from gencode
  pam_G1=as.data.frame(strsplit(pam_G1,','));colnames(pam_G1)='a';pam_G1$a=as.character(pam_G1$a)
  pam_G2=pam$name_RNA_anno ### gene name from annotation
  pam_G2=as.data.frame(strsplit(pam_G2,','));colnames(pam_G2)='a';pam_G2$a=as.character(pam_G2$a)
  
  ### if using Linc any lnc,linc get changed to linc
  if ((pam_G1%in%'lncRNA')&(pam_2%in%'lincRNA')) {pam_1$a='lincRNA'}
  
  ### if gene name is gm (gencode) change to HY1 ( soyeong db)
  # if (grepl('Gm',pam_G1)&grepl('HY1',pam_G2)) {pam_G2[pam_G2$a%in%'HY1',]=NA}
  
  ### if gene name is sk (gencode) change to sk ( repeatmakser db)
  if (grepl('sk',pam_G1)&grepl('7SK',pam_G2)) {pam_G2[pam_G2$a%in%'7SK',]=NA}
  
  pam_G3=rbind(pam_G1,pam_G2)
  pam_G3=pam_G3[!is.na(pam_G3$a),,drop=F]
  pam_G3=paste(unique(pam_G3$a),collapse = ',')
  
  
  ########################################################################################################################################
  ########################################################################################################################################
  ### Create Table 
  
  peaksTable_colapsed[x,colnames(peaksTable_double)]=(pam[1,colnames(peaksTable_double),drop=T])
  peaksTable_colapsed[x,'type_simple_comb']=pam_c2
  peaksTable_colapsed[x,'type_comb']=pam_c
  peaksTable_colapsed[x,'gene_name_comb']=pam_G3
  
  remove('pam_1','pam_2','pam_c','pam_c2','pam_G1','pam_G2','pam_G3','pam')
}
p=peaksTable_colapsed
rnatype=p[,c('ID','ensembl_gene_id','external_gene_name','gene_type',"gene_biotype_ALL",'type_RNA_anno','name_RNA_anno','type_simple_comb','type_comb','gene_name_comb')]


peaksTable_colapsed=rbind(peaksTable_colapsed,peaksTable_single)
peaksTable_colapsed=peaksTable_colapsed[!is.na(peaksTable_colapsed$ID),]

colnames(peaksTable_colapsed)[colnames(peaksTable_colapsed)%in%c("ensembl_gene_id","external_gene_name","gene_type","gene_biotype_ALL","type_RNA_anno","name_RNA_anno","type_simple_comb","type_comb","gene_name_comb")]=
  paste0(nmeprfix,c("ensembl_gene_id","external_gene_name","gene_type","gene_biotype_ALL","type_RNA_anno","name_RNA_anno","type_simple_comb","type_comb","gene_name_comb"))
PeaksdataOut=peaksTable_colapsed



##########################################################################################
## change lincRNA,RNA to RNA only
## do not change the name if there is actually a lincRNA name separate from the rRNA grepl(',',p$name_RNA_anno)==F)
##########################################################################################



#############################################################################################################
#############################################################################################################

### Identify if CLIP peak overlaps with Intron or Exonic region   
  
#   Using GTF file from GENCODE v24
# Peaks were annotated by whether they overlap with Host gene intron/exon region
# Intron coordinates were calculated from GTF file.
# 
# A second column was added to idenify if the peak also overlapped with the 5'UTR 3'UTR or CDS (Column: Featrue 2)

#############################################################################################################
#############################################################################################################

AnnoTable=mm10_ExnItrn[grep('protein_coding',mm10_ExnItrn$gene_type),]
peaksTable=PeaksdataOut
ColumnName=paste0(nmeprfix,c('ensembl_gene_id','external_gene_name'))

Annotable=AnnoTable
s=peaksTable
s=separate(s,ID,into=c('chr','start'),sep=":",remove=F)
s=separate(s,start,into=c('start','end'),sep="-",remove=F)
s.GR <- GRanges(seqnames = as.character(s$chr), ranges=IRanges(start = as.numeric(s$start), end = as.numeric(s$end)),strand = s$strand,ID=s$ID )

anno.GR <- GRanges(seqnames = as.character(Annotable$chr), ranges=IRanges(start = as.numeric(Annotable$start), end = as.numeric(Annotable$end)),strand = Annotable$strand,transcript_id=Annotable$transcript_id,feature=Annotable$feature,exon_number=Annotable$exon_number)

q =s.GR
s=anno.GR
xo=as.data.frame(GenomicRanges::findOverlaps(q,s,type = "any",ignore.strand=F))

qh=as.data.frame(q[xo$queryHits],row.names = NULL)
sh=as.data.frame(s[xo$subjectHits],row.names = NULL);colnames(sh)=paste0(colnames(sh),"_anno")


exoninof=cbind(qh,sh)


peaksTable[,paste0(nmeprfix,'feature')]=NA
peaksTable[,paste0(nmeprfix,'exon_number')]=NA
peaksTable[,paste0(nmeprfix,'intron_number')]=NA
peaksTable[,paste0(nmeprfix,'intron_5pStart')]=NA
peaksTable[,paste0(nmeprfix,'intron_length')]=NA
peaksTable[,paste0(nmeprfix,'exon_length')]=NA

if (xopp==1) {
###########
#### Calc distance of 5' peak to 5' intron/exondist
exoninof_pos=exoninof[exoninof$strand%in%'+',]
exoninof_neg=exoninof[exoninof$strand%in%'-',]

exoninof_pos[,paste0(nmeprfix,'feature_Distance')]=((exoninof_pos$start-exoninof_pos$start_anno)/abs(exoninof_pos$end_anno-exoninof_pos$start_anno))*100
exoninof_pos[,paste0(nmeprfix,'feature_5pStart')]=exoninof_pos$start_anno
exoninof_pos[,paste0(nmeprfix,'feature_length')]=abs(exoninof_pos$end_anno-exoninof_pos$start_anno)


exoninof_neg[,paste0(nmeprfix,'feature_Distance')]=((exoninof_neg$end_anno-exoninof_neg$end)/abs(exoninof_neg$end_anno-exoninof_neg$start_anno))*100
exoninof_neg[,paste0(nmeprfix,'feature_5pStart')]=exoninof_neg$end_anno
exoninof_neg[,paste0(nmeprfix,'feature_length')]=abs(exoninof_neg$end_anno-exoninof_neg$start_anno)


exoninof=rbind(exoninof_pos,exoninof_neg);
# exoninof[,paste0(nmeprfix,'feature_Distance')]=format(exoninof[,paste0(nmeprfix,'feature_Distance')], scientific=F)

peaksTable[,paste0(nmeprfix,'Exn_start_dist')]=NA
peaksTable[,paste0(nmeprfix,'Intron_start_dist')]=NA
peaksTable[,paste0(nmeprfix,'Intron_5pStart')]=NA
peaksTable[,paste0(nmeprfix,'Exn_5pStart')]=NA

}
# xxx1
########### 
# ID="chr19:5490528-5490574"
# g=exoninof[exoninof$ID%in%ID,]
for (x in 1:nrow(peaksTable)) {
  l=peaksTable[x,c('ID','ID2')]
  g=exoninof[exoninof$ID%in%l$ID,]
  
  if (nrow(g)>0) {
    g_e=g[g$feature_anno%in%'exon',]
    g_i=g[g$feature_anno%in%'intron',]
    
    
    gname=paste(unique(g$feature_anno),collapse = ",")
    gname=gsub('NA,',"",gname);gname=gsub(',NA',"",gname)
    peaksTable[x,paste0(nmeprfix,'feature')]=gsub('NA,',"",gname)
    
    if (nrow(g_e)>0) {
      gname=paste(unique(g_e$exon_number_anno),collapse = ",")
      gname=gsub('NA,',"",gname);gname=gsub(',NA',"",gname)
      peaksTable[x,paste0(nmeprfix,'exon_number')]=gsub('NA,',"",gname)
      
      if (xopp==1) {
        g_e[g_e[,paste0(nmeprfix,'feature_Distance')]<0,paste0(nmeprfix,'feature_Distance')]=0
        peaksTable[x,paste0(nmeprfix,'Exn_start_dist')]=mean(as.numeric(g_e[,paste0(nmeprfix,'feature_Distance')]))
        # peaksTable[x,paste0(nmeprfix,'Exn_5pStart')]=paste(unique(as.numeric(g_e[,paste0(nmeprfix,'feature_5pStart')])),collapse = ", ")
        # peaksTable[x,paste0(nmeprfix,'Exn_5pStart')]=g_e[g_e$transcript_id_anno%in%canonical$transcript,paste0(nmeprfix,'feature_5pStart')]
      }
    }
    
    if (nrow(g_i)>0) {
      gname=paste(unique(g_i$exon_number_anno),collapse = ",")
      gname=gsub('NA,',"",gname);gname=gsub(',NA',"",gname)
      peaksTable[x,paste0(nmeprfix,'intron_number')]=gsub('NA,',"",gname)
      
      if (xopp==1) {
        g_i[g_i[,paste0(nmeprfix,'feature_Distance')]<0,paste0(nmeprfix,'feature_Distance')]=0
        peaksTable[x,paste0(nmeprfix,'Intron_start_dist')]=mean(as.numeric(g_i[,paste0(nmeprfix,'feature_Distance')]))
        # peaksTable[x,paste0(nmeprfix,'Intron_5pStart')]=paste(unique(as.numeric(g_i[,paste0(nmeprfix,'feature_5pStart')])),collapse = ", ")
        # peaksTable[x,paste0(nmeprfix,'Intron_5pStart')]=g_i[g_i$transcript_id_anno%in%canonical$transcript,paste0(nmeprfix,'feature_5pStart')]
        
      }
    }
  }
    if (nrow(g)==0) {
      if (xopp==1) {peaksTable[x,paste0(nmeprfix,c("feature",'exon_number','intron_number','Exn_start_dist','Intron_start_dist'))]=NA}
      if (xopp==2) {peaksTable[x,paste0(nmeprfix,c("feature",'exon_number','intron_number'))]=NA}
  } 
}

PeaksdataOut=peaksTable
remove(qh,sh,xo,g,l,q,s)     


#############################################################################################################
#############################################################################################################
### IDENTIFY PEAKS IN REPEAT REGIONS   

# Annotate all repeat regions/Classes identified in Repeatmasker Annotation file (UCSC Table browser)  
# Data was not filtered based on any of the identified Repeats.  
#  1) LINE/SINE   
#  2) LTR   
#  3) DNA   
#  4) Satalites   
#  5) Simple Repeats   
#  6) Low Complexity   
#  7) Other   
#  8) Unknown  
#############################################################################################################
#############################################################################################################

#suppressWarnings()


PeaksdataOut=PeaksdataOut[!is.na(PeaksdataOut$ID),]

PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_LINE_SINE'),rmsk_GRCm38_LISI,PeaksdataOut)
PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_LTR'),rmsk_GRCm38_LTR,PeaksdataOut)
PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_DNA'),rmsk_GRCm38_DNA,PeaksdataOut)
PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_Satalites'),rmsk_GRCm38_sat,PeaksdataOut)
PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_Simple_Repeats'),rmsk_GRCm38_SR,PeaksdataOut)
PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_Low_Complexity'),rmsk_GRCm38_LowComplx,PeaksdataOut)
PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_Other'),rmsk_GRCm38_Other,PeaksdataOut)
PeaksdataOut=rpmsk_anno(paste0(nmeprfix,'Repeat_Unknown'),rmsk_GRCm38_unknown,PeaksdataOut)


p=PeaksdataOut
p[,paste0(nmeprfix,'repeat_comb')]=NA
repcol=paste0(nmeprfix,c('Repeat_LINE_SINE','Repeat_LTR','Repeat_DNA','Repeat_Satalites','Repeat_Simple_Repeats','Repeat_Low_Complexity','Repeat_Other','Repeat_Unknown'))

dup=p[rowSums((p[,repcol]>0),na.rm = T)>0,'ID']

peaksTable_single=p[!(p$ID%in%dup),]
peaksTable_double=p[(p$ID%in%dup),]

u=unique(peaksTable_double$ID)
peaksTable_colapsed=as.data.frame(matrix(nrow=length(u),ncol=ncol(peaksTable_double)));
colnames(peaksTable_colapsed)=colnames(peaksTable_double)


for(x in 1:length(u)){
  p=u[x]
  pam=peaksTable_double[peaksTable_double$ID%in%p,]
  pam_c=((pam[1,repcol,drop=F]))
  
  pmat=as.data.frame(matrix(nrow=length(pam_c),ncol=1));colnames(pmat)='a'
  pmat$a=t(pam_c)
  pmat=pmat[is.na(pmat)==F,]
  pmat=paste(unique(pmat),collapse = ',')
  
  ######################     
  
  peaksTable_colapsed[x,colnames(peaksTable_double)]=(pam[1,colnames(peaksTable_double),drop=T])
  peaksTable_colapsed[x,paste0(nmeprfix,'repeat_comb')]=pmat
  
  remove('p','pam','pam_c')
}


peaksTable_colapsed=rbind(peaksTable_colapsed,peaksTable_single)
peaksTable_colapsed=peaksTable_colapsed[!is.na(peaksTable_colapsed$ID),]

p=peaksTable_colapsed
# rnatype=p[,c('ID','ensembl_gene_id','external_gene_name','gene_type',"gene_biotype_ALL",'type_RNA_anno','name_RNA_anno','type_simple_comb','type_comb','gene_name_comb')]

PeaksdataOut=peaksTable_colapsed
# xxx







#############################################################################################################
#############################################################################################################
### Asigning Clip peak attributes   

# Not all Peaks overlap with a single feature so peak assignments were assigned by priority:  
# 
# ncRNA > Protein coding : Exonic > repeats > Pseudogene > Protein Coding : Intronic  
# 
# All annotations from RNA type, Repeat regions, and Intronic/exonic regions are annoted in the Table.   

#############################################################################################################
#############################################################################################################

# write.table(peaksTable,file=("Clip_Peaks_2.txt"), sep = "\t", row.names = FALSE, col.names = T, append = F, quote= FALSE)
# ncRNA > exonic > repeats > pseudogene > antisense > intronic > lncRNA > noFeature 


PeaksdataOut[,paste0(nmeprfix,'Comb_biotype_exon')]=NA

# 1. ncRNA        
comp=( grepl('ncRNA',PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')])& ((PeaksdataOut[,paste0(nmeprfix,'type_comb')]%in%'lncRNA')==F) )
# sum(as.numeric(comp))
PeaksdataOut[comp,paste0(nmeprfix,'Comb_biotype_exon')]='ncRNA'

# 2. protein coding - Exonic
comp=( is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_biotype_exon')])& grepl('protein_coding',PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')])& (PeaksdataOut[,paste0(nmeprfix,'feature')]%in%'exon') )
# sum(as.numeric(comp))
PeaksdataOut[comp,paste0(nmeprfix,'Comb_biotype_exon')]=paste0('protein_coding: ',PeaksdataOut[comp,paste0(nmeprfix,'feature')])

# 3. repeats
comp=( is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_biotype_exon')])& (is.na(PeaksdataOut[,paste0(nmeprfix,'repeat_comb')])==F) )
# sum(as.numeric(comp))
PeaksdataOut[comp,paste0(nmeprfix,'Comb_biotype_exon')]="Repeat Element"

# 4. Pseudogene
comp=( is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_biotype_exon')])& grepl('pseudogene',PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]) )
# sum(as.numeric(comp))
PeaksdataOut[comp,paste0(nmeprfix,'Comb_biotype_exon')]='pseudogene'

# 5. intron
comp=( is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_biotype_exon')])& grepl('protein_coding',PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')]) )
# sum(as.numeric(comp))
PeaksdataOut[comp,paste0(nmeprfix,'Comb_biotype_exon')]=paste0('protein_coding: Intron')

# 6. lncRNA
comp=( is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_biotype_exon')])& grepl('ncRNA',PeaksdataOut[,paste0(nmeprfix,'type_simple_comb')])& ((PeaksdataOut[,paste0(nmeprfix,'type_comb')]%in%'lncRNA')==T) )
# sum(as.numeric(comp))
PeaksdataOut[comp,paste0(nmeprfix,'Comb_biotype_exon')]='lncRNA'

PeaksdataOut[is.na(PeaksdataOut[,paste0(nmeprfix,'Comb_biotype_exon')]),paste0(nmeprfix,'Comb_biotype_exon')]='no Feature'


PeaksdataOut[,paste0(nmeprfix,'Comb_biotype_exon')]=factor(PeaksdataOut[,paste0(nmeprfix,'Comb_biotype_exon')], levels = c("ncRNA", "protein_coding: exon", "Repeat Element","pseudogene","Antisense Feature","protein_coding: Intron","lncRNA","no Feature"))


##############################
#### RNA subtypes
##############################

p=PeaksdataOut
p[,paste0(nmeprfix,'Comb_biotype_ncRNA')]=NA
p1=p[p[,paste0(nmeprfix,'Comb_biotype_exon')]%in%'ncRNA',]
p2=p[!p[,paste0(nmeprfix,'Comb_biotype_exon')]%in%'ncRNA',]

p1[,paste0(nmeprfix,'Comb_biotype_ncRNA')]=p1[,paste0(nmeprfix,'type_comb')]
### Protein coding + ncRNA annotations -> ncRNA subtype
p1[,paste0(nmeprfix,'Comb_biotype_ncRNA')]=gsub('protein_coding,',"",p1[,paste0(nmeprfix,'Comb_biotype_ncRNA')])
p1[,paste0(nmeprfix,'Comb_biotype_ncRNA')]=gsub(',protein_coding',"",p1[,paste0(nmeprfix,'Comb_biotype_ncRNA')])

### any double annotations with lncRNA become second annotation only
p1[,paste0(nmeprfix,'Comb_biotype_ncRNA')]=gsub('lncRNA,',"",p1[,paste0(nmeprfix,'Comb_biotype_ncRNA')])

### any double annotations with yRNA become yRNA only
p1[grep('yRNA',p1[,paste0(nmeprfix,'Comb_biotype_ncRNA')]),paste0(nmeprfix,'Comb_biotype_ncRNA')]='yRNA'

### tRNA takes priority over miRNA
p1[,paste0(nmeprfix,'Comb_biotype_ncRNA')]=gsub('miRNA,tRNA',"miRNA",p1[,paste0(nmeprfix,'Comb_biotype_ncRNA')])

### miRNA takes priority over rRNA
p1[,paste0(nmeprfix,'Comb_biotype_ncRNA')]=gsub('miRNA,rRNA',"miRNA",p1[,paste0(nmeprfix,'Comb_biotype_ncRNA')])

PeaksdataOut=rbind(p1,p2)


if (xopp==1) {PeaksdataOut_same=PeaksdataOut}
if (xopp==2) {PeaksdataOut=merge(PeaksdataOut_same,PeaksdataOut[,colnames(PeaksdataOut)[!colnames(PeaksdataOut)%in%c("chr","start","end","strand","ID2" )]],by='ID')}

}### for pos or neg anno



PeaksdataOut$Comb_biotype_exon_Oppo=NA

# 1. ncRNA        
comp=(grepl('ncRNA',PeaksdataOut[,paste0('Same_','type_simple_comb')])& (PeaksdataOut[,paste0('Same_','type_comb')]%in%'lncRNA')==F)
# sum(as.numeric(comp))
PeaksdataOut[comp,'Comb_biotype_exon_Oppo']='ncRNA'

# 2. protein coding - Exonic
comp=(is.na(PeaksdataOut[,'Comb_biotype_exon_Oppo'])&grepl('protein_coding',PeaksdataOut[,paste0('Same_','type_simple_comb')])&(PeaksdataOut[,paste0('Same_','feature')]%in%'exon'))
# sum(as.numeric(comp))
PeaksdataOut[comp,'Comb_biotype_exon_Oppo']=paste0('protein_coding: ',PeaksdataOut[comp,paste0('Same_','feature')])

# 3. repeats
comp=(is.na(PeaksdataOut[,'Comb_biotype_exon_Oppo'])&(is.na(PeaksdataOut[,paste0('Same_','repeat_comb')])==F) )
# sum(as.numeric(comp))
PeaksdataOut[comp,'Comb_biotype_exon_Oppo']="Repeat Element"

# 4. Pseudogene
comp=(is.na(PeaksdataOut[,'Comb_biotype_exon_Oppo'])& grepl('pseudogene',PeaksdataOut[,paste0('Same_','type_simple_comb')]))
# sum(as.numeric(comp))
PeaksdataOut[comp,'Comb_biotype_exon_Oppo']='pseudogene'

# 5. Antisense
comp=( is.na(PeaksdataOut[,'Comb_biotype_exon_Oppo'])& ((PeaksdataOut[,paste0('Oppo_','Comb_biotype_exon')]%in%'no Feature')==F) )
# sum(as.numeric(comp))
PeaksdataOut[comp,'Comb_biotype_exon_Oppo']='Antisense Feature'

# 6. intron
comp=(is.na(PeaksdataOut[,'Comb_biotype_exon_Oppo'])& grepl('protein_coding',PeaksdataOut[,paste0('Same_','type_simple_comb')]))
# sum(as.numeric(comp))
PeaksdataOut[comp,'Comb_biotype_exon_Oppo']=paste0('protein_coding: Intron')

# 7. lncRNA
comp=(is.na(PeaksdataOut[,'Comb_biotype_exon_Oppo'])& grepl('ncRNA',PeaksdataOut[,paste0('Same_','type_simple_comb')])& ((PeaksdataOut[,paste0('Same_','type_comb')]%in%'lncRNA')==T) )
# sum(as.numeric(comp))
PeaksdataOut[comp,'Comb_biotype_exon_Oppo']='lncRNA'


PeaksdataOut[is.na(PeaksdataOut[,'Comb_biotype_exon_Oppo']),'Comb_biotype_exon_Oppo']='no Feature'

PeaksdataOut[,'Comb_biotype_exon_Oppo']=factor(PeaksdataOut[,'Comb_biotype_exon_Oppo'], levels = c("ncRNA", "protein_coding: exon", "Repeat Element","pseudogene","Antisense Feature","protein_coding: Intron","lncRNA","no Feature"))


### only distances for final classification
# PeaksdataOut[(PeaksdataOut$Comb_biotype_exon_Oppo%in%c("protein_coding: exon")==F),c("Same_Exn_start_dist")]=NA
# PeaksdataOut[(PeaksdataOut$Comb_biotype_exon_Oppo%in%c("protein_coding: Intron")==F),c("Same_Intron_start_dist")]=NA
### distances to any protein coading annotation
PeaksdataOut[grepl('protein_coding',PeaksdataOut$Same_type_comb)==F,c("Same_Exn_start_dist")]=NA
PeaksdataOut[grepl('protein_coding',PeaksdataOut$Same_type_comb)==F,c("Same_Intron_start_dist")]=NA

# View(PeaksdataOut[(PeaksdataOut$Comb_biotype_exon_Oppo%in%c("protein_coding: exon","protein_coding: Intron")==F),])

return(PeaksdataOut)
}

