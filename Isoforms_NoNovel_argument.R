# remove variables
rm(list=ls(all=TRUE))
library(cummeRbund)
setwd("~/Documents/NCSU/RNAseq_CMM_BL_Unique/Analysis_ALL/05_BL_Uni_cuffdiff_ALL")
getwd()
cuff <- readCufflinks(rebuild=F)
cuff

TSS.filename <- "TSS_TAIR_only.data"
ISO.filename <- "Isoform_TAIR_only.data"


##############################
# some basic quality control #
##############################

#d<-dispersionPlot(genes(cuff))
#d
# CummeRbund has many QC plots, see manual
# here just some examples that focus on fpkm
# distribution and overall sample consistency

# fpkm distribution of genes
# (looks ok for me)
#densRep <- csDensity(genes(cuff),replicates=F);
#densRep

# fpkm box plot
#brep = csBoxplot(genes(cuff), replicates=T);
#brep

# fpkm distribution of isoforms
#idensRep <- csDensity(isoforms(cuff),replicates=T)
#idensRep

# fpkm of TSS
#TSSRep <- csDensity(TSS(cuff),replicates=T)
#TSSRep

# scatter plots of fpkm and count values
# normalization seems ok, but higher variability.
#s1<-csScatterMatrix(genes(cuff), replicates=T)
#s1

#s2<-csScatterMatrix(genes(cuff), replicates=T, useCounts=TRUE)
#s2

# compare samples & replicates

# dendrogram of replicates
# is consistent with experiment design
# rep 1 is suspicious
#dend1.rep <- csDendro(genes(cuff),replicates=T)

# compare replicates with Jensen-Shannon Distance
#RepDistHeat<-csDistHeat(genes(cuff),replicates=T)
#RepDistHeat

######################
# expressed features #
######################

##########################
# analysis on gene level #
##########################

# define expressed features
g.fpkm <- fpkm(genes(cuff))

# compile gene information
# cuff id,TAIR id, locus location
gene.info <- annotation(genes(cuff))[,c(1,4,5)]
#
# fpkm values
sample.names<-samples(cuff)$sample_name
for (i in sample.names) {
  fpkm.data <- g.fpkm[g.fpkm$sample_name==i,c(1,3)]
  names(fpkm.data) <- c("gene_id", paste("gene.fpkm",i,sep="."))
  gene.info <- merge(gene.info, fpkm.data, by="gene_id")
}
#
# expressed genes
# I have chosen conf_lo>0 & quant_status=="OK"
# other defs are possible, other defs are possible.
for (i in sample.names) {
  g.expressed <- g.fpkm[g.fpkm$sample_name==i,]
  g.expressed$fpkm <- 0
  g.expressed[g.expressed$conf_lo>0 & g.expressed$quant_status=="OK",3]<-1
  g.expressed <- g.expressed[,c(1,3)]
  names(g.expressed) <- c("gene_id", paste("gene.expressed",i,sep="."))
  gene.info <- merge(gene.info, g.expressed, by="gene_id")
}

# differentially expressed genes
# I have chosen 0.01 as
# significance level
mylevel <- 'genes';
myalpha <- 0.01;

# differential only in YFP_POS/NEG comparison
G.id.1.sig <- getSig(cuff,alpha=myalpha,level=mylevel,x="ALL_SORT",y="NO_SORT", useCuffMTC=F)
# x="ALL_SORT",y="NO_SORT": 4418
G.id.2.sig <- getSig(cuff,alpha=myalpha,level=mylevel,x="YFP_POS",y="YFP_NEG", useCuffMTC=F)
# x="YFP_POS",y="YFP_NEG" : 7195

# add to gene info
gene.info$gene.sig.YFP <-0
gene.info[gene.info$gene_id %in% G.id.2.sig, dim(gene.info)[[2]]]<-1
gene.info$gene.sig.SORT <-0
gene.info[gene.info$gene_id %in% G.id.1.sig, dim(gene.info)[[2]]]<-1

# TODO: would be better ot build new cuffset object
# that consists of only expressed genes

###############################
# differential promoter level #
###############################

# used higher significance level!
myalpha <- 0.1
mylevel <- 'promoters'

P.id.1.sig <- getSig(cuff,alpha=myalpha,level=mylevel,x="ALL_SORT",y="NO_SORT", useCuffMTC=F)
P.id.2.sig <- getSig(cuff,alpha=myalpha,level=mylevel,x="YFP_NEG",y="YFP_POS", useCuffMTC=F)
P.id.3.sig <- getSig(cuff,alpha=myalpha,level=mylevel,x="NO_SORT",y="YFP_POS", useCuffMTC=F)
# x="ALL_SORT",y="NO_SORT": 231
# x="YFP_POS",y="YFP_NEG" : 30

# differential only in YFP_POS/NEG comparison
P.sig <- setdiff(P.id.2.sig,P.id.1.sig)
# 24

# add to gene info
gene.info$gene.promoter.sig.YFP <-0
gene.info[gene.info$gene_id %in% P.id.2.sig, dim(gene.info)[[2]]]<-1
gene.info$gene.promoter.sig.YFP.NO_SORT <-0
gene.info[gene.info$gene_id %in% P.id.3.sig, dim(gene.info)[[2]]]<-1
gene.info$gene.promoter.sig.SORT <-0
gene.info[gene.info$gene_id %in% P.id.1.sig, dim(gene.info)[[2]]]<-1

# visual inspection of significant genes
l <- length(P.sig)
#old.par <- par(ask=TRUE)
for (i in 1:l) {
  gene <- getGenes(cuff, P.sig[i])
  p <- expressionPlot(TSS(gene), replicates=T) + theme_minimal()
  print(p) 
}
#par(old.par)

#this is a list so you cant merge 
#motif = isoform.gene.info[isoform.gene.info[,1] %in% P.id.2.sig,]



#############################
# analysis on isoform level #
#############################

# compile isoform information
# isoform id,gene id, TAIR gene name, TSS, class code, TAIR transcript, length
isoform.info<-annotation(isoforms(cuff))[,c(1,2,4,5,6,7,9)]

#
# fpkm values
i.fpkm = fpkm(isoforms(cuff))

sample.names<-samples(cuff)$sample_name
for (i in sample.names) {
  fpkm.data <- i.fpkm[i.fpkm$sample_name==i,c(1,3)]
  names(fpkm.data) <- c("isoform_id", paste("isoform.fpkm",i,sep="."))
  isoform.info <- merge(isoform.info, fpkm.data, by="isoform_id")
}
#
# expressed isoforms
for (i in sample.names) {
  i.expressed <- i.fpkm[i.fpkm$sample_name==i,]
  i.expressed$fpkm <- 0
  i.expressed[i.expressed$conf_lo>0 & i.expressed$quant_status=="OK",3]<-1
  i.expressed <- i.expressed[,c(1,3)]
  names(i.expressed) <- c("isoform_id", paste("isoform.expressed",i,sep="."))
  isoform.info <- merge(isoform.info, i.expressed, by="isoform_id")
}

# differentially expressed isoforms
# I have chosen 0.01 as
# significance level
mylevel <- 'isoforms';
myalpha <- 0.01;

# differential in YFP_POS/NEG comparison
I.id.1.sig <- getSig(cuff,alpha=myalpha,level=mylevel,x="ALL_SORT",y="NO_SORT", useCuffMTC=F)
# x="ALL_SORT",y="NO_SORT": 3938
I.id.2.sig <- getSig(cuff,alpha=myalpha,level=mylevel,x="YFP_POS",y="YFP_NEG", useCuffMTC=F)
# x="YFP_POS",y="YFP_NEG" : 6411
# differential in YFP_POS/NO_SORT comparison
I.id.3.sig <- getSig(cuff,alpha=myalpha,level=mylevel,x="YFP_POS",y="NO_SORT", useCuffMTC=F)

# add to isoform info
isoform.info$isoform.sig.YFP <-0
isoform.info[isoform.info$isoform_id %in% I.id.2.sig, dim(isoform.info)[[2]]]<-1
isoform.info$isoform.sig.YFP.NO_SORT <-0
isoform.info[isoform.info$isoform_id %in% I.id.3.sig, dim(isoform.info)[[2]]]<-1
isoform.info$isoform.sig.SORT <-0
isoform.info[isoform.info$isoform_id %in% I.id.1.sig, dim(isoform.info)[[2]]]<-1

# merge isoform and gene data
isoform.gene.info <- merge(isoform.info, gene.info[,-2], by="gene_id")
# add relative isoform fraction
isoform.gene.info$isoform.rel.fpkm.YFP_NEG <- isoform.gene.info$isoform.fpkm.YFP_NEG/isoform.gene.info$gene.fpkm.YFP_NEG
isoform.gene.info$isoform.rel.fpkm.YFP_POS <- isoform.gene.info$isoform.fpkm.YFP_POS/isoform.gene.info$gene.fpkm.YFP_POS
isoform.gene.info$isoform.rel.fpkm.NO_SORT <- isoform.gene.info$isoform.fpkm.NO_SORT/isoform.gene.info$gene.fpkm.NO_SORT
isoform.gene.info$isoform.rel.fpkm.ALL_SORT <- isoform.gene.info$isoform.fpkm.ALL_SORT/isoform.gene.info$gene.fpkm.ALL_SORT
# add fpkm fraction YFP+/YFP-, YFP+/No_SORT, YFP+/All_SORT
isoform.gene.info$isoform.frac.fpkm.POS.NEG <- isoform.gene.info$isoform.fpkm.YFP_NEG/isoform.gene.info$isoform.fpkm.YFP_POS
isoform.gene.info$isoform.frac.fpkm.POS.NO <- isoform.gene.info$isoform.fpkm.NO_SORT/isoform.gene.info$isoform.fpkm.YFP_POS
isoform.gene.info$isoform.frac.fpkm.POS.ALL <- isoform.gene.info$isoform.fpkm.ALL_SORT/isoform.gene.info$isoform.fpkm.YFP_POS
# replace NaN by 0
isoform.gene.info[is.na(isoform.gene.info)] <- 0
# add diferences of rel.fpkms
isoform.gene.info$isoform.rel.dif.YFP <- isoform.gene.info$isoform.rel.fpkm.YFP_POS-isoform.gene.info$isoform.rel.fpkm.YFP_NEG
isoform.gene.info$isoform.rel.dif.SORT <- isoform.gene.info$isoform.rel.fpkm.ALL_SORT-isoform.gene.info$isoform.rel.fpkm.NO_SORT

# compute number of isoforms per gene
library(dplyr)
attach(isoform.gene.info)

transcripts_per_gene <- group_by(isoform.gene.info,gene_id)
n.trans <-summarise(transcripts_per_gene,
          n_transcripts = n_distinct(isoform_id)
          )

r.trans <- mutate(transcripts_per_gene, 
                  t_rank.YFP_NEG = rank(desc(isoform.rel.fpkm.YFP_NEG),ties.method = "random"),
                  t_rank.YFP_POS = rank(desc(isoform.rel.fpkm.YFP_POS),ties.method = "random"),
                  t_rank.NO_SORT = rank(desc(isoform.rel.fpkm.NO_SORT),ties.method = "random"),
                  t_rank.ALL_SORT = rank(desc(isoform.rel.fpkm.ALL_SORT),ties.method = "random")
                  )
isoform.gene.info2 <- merge(r.trans, n.trans, by="gene_id")

detach(isoform.gene.info)

# TODO: would be better ot build new cuffset object
# that consists of only expressed genes

###################
# analysis of TSS #
###################

# compile TSS information
# TSS id,gene id
TSS.info <- unique(isoform.info[,c(4,2)])

#
# fpkm values
TSS.fpkm = fpkm(TSS(cuff))

sample.names<-samples(cuff)$sample_name
for (i in sample.names) {
  fpkm.data <- TSS.fpkm[TSS.fpkm$sample_name==i,c(1,3)]
  names(fpkm.data) <- c("TSS_group_id", paste("TSS.fpkm",i,sep="."))
  TSS.info <- merge(TSS.info, fpkm.data, by="TSS_group_id")
}
#
# expressed TSSs
for (i in sample.names) {
  TSS.expressed <- TSS.fpkm[TSS.fpkm$sample_name==i,]
  TSS.expressed$fpkm <- 0
  TSS.expressed[TSS.expressed$conf_lo>0 & TSS.expressed$quant_status=="OK",3]<-1
  TSS.expressed <- TSS.expressed[,c(1,3)]
  names(TSS.expressed) <- c("TSS_group_id", paste("TSS.expressed",i,sep="."))
  TSS.info <- merge(TSS.info, TSS.expressed, by="TSS_group_id")
}

# differentially expressed isoforms
# I have chosen 0.01 as
# significance level
mylevel <- 'TSS';
myalpha <- 0.01;

# differential only in YFP_POS/NEG comparison
TSS.id.1.sig <- getSig(cuff,alpha=myalpha,level=mylevel,x="ALL_SORT",y="NO_SORT", useCuffMTC=F)
# x="ALL_SORT",y="NO_SORT": 3938
TSS.id.2.sig <- getSig(cuff,alpha=myalpha,level=mylevel,x="YFP_POS",y="YFP_NEG", useCuffMTC=F)
# x="YFP_POS",y="YFP_NEG" : 6411

# add to isoform info
TSS.info$TSS.sig.YFP <-0
TSS.info[TSS.info$TSS_group_id %in% TSS.id.2.sig, dim(TSS.info)[[2]]]<-1
TSS.info$TSS.sig.SORT <-0
TSS.info[TSS.info$TSS_group_id %in% TSS.id.1.sig, dim(TSS.info)[[2]]]<-1

##################
# splicing level #
##################

# differentially splicing
myalpha <- 0.05
mylevel <- 'splicing'

S.id.1.sig <- getSig(cuff,alpha=myalpha,level=mylevel,x="ALL_SORT",y="NO_SORT", useCuffMTC=F)
S.id.2.sig <- getSig(cuff,alpha=myalpha,level=mylevel,x="YFP_NEG",y="YFP_POS", useCuffMTC=F)
# x="ALL_SORT",y="NO_SORT": 183
# x="YFP_POS",y="YFP_NEG" : 23

# differential only in YFP_POS/NEG comparison
S.sig <- setdiff(S.id.2.sig,S.id.1.sig)
# 21

# add to isoform info
TSS.info$TSS.splicing.sig.YFP <-0
TSS.info[TSS.info$TSS_group_id %in% S.id.2.sig, 13]<-1
TSS.info$TSS.splicing.sig.SORT <-0
TSS.info[TSS.info$TSS_group_id %in% S.id.1.sig, 14]<-1

# visual inspection of AS genes
l <- length(S.sig)
#old.par <- par(ask=TRUE)
for (i in 1:l) {
  gene <- getGenes(cuff, S.sig[i])
  p <- expressionPlot(isoforms(gene), replicates=F, showStatus=FALSE,drawSummary=FALSE)
  print(p)
}
#par(old.par)

# merge TSS and gene data
TSS.gene.info <- merge(TSS.info, gene.info, by="gene_id")
# add relative isoform fraction
TSS.gene.info$TSS.rel.fpkm.YFP_NEG <- TSS.gene.info$TSS.fpkm.YFP_NEG/TSS.gene.info$gene.fpkm.YFP_NEG
TSS.gene.info$TSS.rel.fpkm.YFP_POS <- TSS.gene.info$TSS.fpkm.YFP_POS/TSS.gene.info$gene.fpkm.YFP_POS
TSS.gene.info$TSS.rel.fpkm.NO_SORT <- TSS.gene.info$TSS.fpkm.NO_SORT/TSS.gene.info$gene.fpkm.NO_SORT
TSS.gene.info$TSS.rel.fpkm.ALL_SORT <- TSS.gene.info$TSS.fpkm.ALL_SORT/TSS.gene.info$gene.fpkm.ALL_SORT
# replace NaN by 0
TSS.gene.info[is.na(TSS.gene.info)] <- 0
# add diferences of rel.fpkms
TSS.gene.info$TSS.rel.dif.YFP <- TSS.gene.info$TSS.rel.fpkm.YFP_POS-TSS.gene.info$TSS.rel.fpkm.YFP_NEG
TSS.gene.info$TSS.rel.dif.SORT <- TSS.gene.info$TSS.rel.fpkm.ALL_SORT-TSS.gene.info$TSS.rel.fpkm.NO_SORT

attach(TSS.gene.info)

TSS_per_gene <- group_by(TSS.gene.info,gene_id)
n.TSS_groups <-summarise(TSS_per_gene,
                         n_TSS = n_distinct(TSS_group_id)
)

r.TSS <- mutate(TSS_per_gene, 
                TSS_rank.YFP_NEG = rank(desc(TSS.rel.fpkm.YFP_NEG),ties.method = "random"),
                TSS_rank.YFP_POS = rank(desc(TSS.rel.fpkm.YFP_POS),ties.method = "random"),
                TSS_rank.NO_SORT = rank(desc(TSS.rel.fpkm.NO_SORT),ties.method = "random"),
                TSS_rank.ALL_SORT = rank(desc(TSS.rel.fpkm.ALL_SORT),ties.method = "random")
)
TSS.gene.info2 <- merge(r.TSS, n.TSS_groups, by="gene_id")

#write.table(isoform.gene.info2, "ISO.filename_v2", sep=",", row.names=FALSE)

####################################
# Graph for paper these 3 into 1 graph   #
####################################
#  XLOC_012987 = AT2G33815 = cis-natural antisense       
#  XLOC_006584 = AT1G48742 = MICRORNA157D,
#  XLOC_010258 = At2g33810 = SPL3 

myGeneId1<-"XLOC_012987" 
myGeneId1<-getGene(cuff,myGeneId1) 
myGeneId1 
gl.rep<-expressionPlot(myGeneId1,replicates=TRUE) 
gl.rep 

myGeneId2<-"XLOC_006584" 
myGeneId2<-getGene(cuff,myGeneId2) 
myGeneId2 
gl.rep<-expressionPlot(myGeneId2,replicates=TRUE) 
gl.rep 

myGeneId3<-"XLOC_010258" 
myGeneId3<-getGene(cuff,myGeneId3) 
myGeneId3 
gl.rep<-expressionPlot(myGeneId3,replicates=TRUE) 
gl.rep 


### graph with the master Qiwen

fpkm.graph<- read.table("isoforms.fpkm_tracking", header=TRUE)

myGeneId1<- fpkm.graph[fpkm.graph$gene_id=="XLOC_012987",]
myGeneId2 <- fpkm.graph[fpkm.graph$gene_id=="XLOC_010258",]

id_fpkm = c(myGeneId1$YFP_NEG_FPKM, myGeneId1$YFP_POS_FPKM, myGeneId1$NO_SORT_FPKM, myGeneId1$ALL_SORT_FPKM,myGeneId2$YFP_NEG_FPKM, myGeneId2$YFP_POS_FPKM, myGeneId2$NO_SORT_FPKM, myGeneId2$ALL_SORT_FPKM)
id_celltype = factor(c("YFP_NEG","YFP_POS","NO_SORT","ALL_SORT","YFP_NEG","YFP_POS","NO_SORT","ALL_SORT"), levels=c("YFP_NEG","YFP_POS","NO_SORT","ALL_SORT"))

dat1 = data.frame(celltype=id_celltype,FPKM=id_fpkm)
dat1$genes = c("Antisense","Antisense","Antisense","Antisense","SPL3","SPL3","SPL3","SPL3")



ggplot(data=dat1, aes(x=celltype, y=FPKM, group=genes, colour=genes)) +
        geom_line() +
        geom_point()+geom_errorbar(aes(ymin=c(myGeneId1$YFP_NEG_conf_lo, myGeneId1$YFP_POS_conf_lo, myGeneId1$NO_SORT_conf_lo, myGeneId1$ALL_SORT_conf_lo,myGeneId2$YFP_NEG_conf_lo, myGeneId2$YFP_POS_conf_lo, myGeneId2$NO_SORT_conf_lo, myGeneId2$ALL_SORT_conf_lo), ymax=c(myGeneId1$YFP_NEG_conf_hi, myGeneId1$YFP_POS_conf_hi, myGeneId1$NO_SORT_conf_hi, myGeneId1$ALL_SORT_conf_hi,myGeneId2$YFP_NEG_conf_hi, myGeneId2$YFP_POS_conf_hi, myGeneId2$NO_SORT_conf_hi, myGeneId2$ALL_SORT_conf_hi)),
                                  size=.3,    # Thinner lines
                                  width=.2,
                                  position=position_dodge(0))

### with 3 genes

myGeneId3 <- fpkm.graph[fpkm.graph$gene_id=="XLOC_006584",]

id_fpkm = c(myGeneId1$YFP_NEG_FPKM, myGeneId1$YFP_POS_FPKM, myGeneId1$NO_SORT_FPKM, myGeneId1$ALL_SORT_FPKM,
            myGeneId2$YFP_NEG_FPKM, myGeneId2$YFP_POS_FPKM, myGeneId2$NO_SORT_FPKM, myGeneId2$ALL_SORT_FPKM,
            myGeneId3$YFP_NEG_FPKM, myGeneId3$YFP_POS_FPKM, myGeneId3$NO_SORT_FPKM, myGeneId3$ALL_SORT_FPKM)
id_celltype = factor(c("YFP_NEG","YFP_POS","NO_SORT","ALL_SORT","YFP_NEG","YFP_POS","NO_SORT","ALL_SORT","YFP_NEG","YFP_POS","NO_SORT","ALL_SORT"), levels=c("YFP_NEG","YFP_POS","NO_SORT","ALL_SORT"))

dat2 = data.frame(celltype=id_celltype,FPKM=id_fpkm)
dat2$genes = c("Antisense","Antisense","Antisense","Antisense","SPL3","SPL3","SPL3","SPL3","miRNA","miRNA","miRNA","miRNA")



ggplot(data=dat2, aes(x=celltype, y=FPKM, group=genes, colour=genes)) +
        geom_bar(stat="identity", position=position_dodge()) +
        scale_fill_brewer(palette="Paired")  +
        geom_errorbar(aes(ymin=c(myGeneId1$YFP_NEG_conf_lo, myGeneId1$YFP_POS_conf_lo, myGeneId1$NO_SORT_conf_lo, myGeneId1$ALL_SORT_conf_lo,
                                              myGeneId2$YFP_NEG_conf_lo, myGeneId2$YFP_POS_conf_lo, myGeneId2$NO_SORT_conf_lo, myGeneId2$ALL_SORT_conf_lo,
                                              myGeneId3$YFP_NEG_conf_lo, myGeneId3$YFP_POS_conf_lo, myGeneId3$NO_SORT_conf_lo, myGeneId3$ALL_SORT_conf_lo), 
                                       ymax=c(myGeneId1$YFP_NEG_conf_hi, myGeneId1$YFP_POS_conf_hi, myGeneId1$NO_SORT_conf_hi, myGeneId1$ALL_SORT_conf_hi,
                                              myGeneId2$YFP_NEG_conf_hi, myGeneId2$YFP_POS_conf_hi, myGeneId2$NO_SORT_conf_hi, myGeneId2$ALL_SORT_conf_hi,
                                              myGeneId3$YFP_NEG_conf_hi, myGeneId3$YFP_POS_conf_hi, myGeneId3$NO_SORT_conf_hi, myGeneId3$ALL_SORT_conf_hi)),
                                   size=.5,    # Thinner lines
                                   width=.2,
                                   position=position_dodge(1))


geom_bar((data=dat2, aes(x=celltype, y=FPKM, group=genes, colour=genes))

         
         ggplot(dat2, aes(x = factor(genes))) + geom_bar(stat = "bin")
   
         ggplot(data=dat2, aes(x=celltype, y=fpkm, group=genes, fill=genes)) +
                 geom_bar(stat="identity", position=position_dodge()) +
                 scale_fill_brewer(palette="Paired") +
                 geom_errorbar(aes(ymin= min,ymax= max),size=.5,
                               width=.2,
                               position=position_dodge(0.9)) + theme_minimal()
### graph trans-acting siRNA (tasi-RNA) AT5G49615 and AT3G17185 
        # XLOC_028213 = AT5G49615 = TAS3B
        # XLOC_014682 = AT3G17185 = TAS3
        # XLOC_028933 = AT5G62000 = ARF2 
        # XLOC_012994 = AT2G33860 = ARF3 
        # XLOC_032523 = AT5G60450 = ARF4

myGeneId3<-"XLOC_032523" 
myGeneId3<-getGene(cuff,myGeneId3) 
myGeneId3 
gl.rep<-expressionPlot(myGeneId3,replicates=TRUE) 
gl.rep 


