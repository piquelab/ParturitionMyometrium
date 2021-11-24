
#############################################################
###  Calculating differentially expressed genes
### 
#############################################################

library(tidyverse)
library(DESeq2)
library(qqman)
library(Seurat)
library(reshape)
library("BiocParallel")
register(MulticoreParam(12))


outFolder <- paste0("7_outputs_DESeq_ConditionsByCluster_bath_library/")
system(paste0("mkdir -p ",outFolder))

## Load cells
sc <- read_rds("6_harmony_cellClass_plots_res0.8_final/SeuratObject.rds")

md <- sc@meta.data

## Load gene annotations. 
anno <- read_rds("3_MergeDemux_Output/anno.rds")

stopifnot(identical(rownames(sc),anno$kbid))

gene_symbol <- anno$gene_name
names(gene_symbol) <- anno$kbid

##afile <- paste0(basefolder,"/../ED2/80K-all.rds")
##all <- readr::read_rds(afile)

##cSize <- md %>% dplyr::count(Location,cluster_name,SNG.BEST.GUESS,Condition)


# seurat_clusters : cluster_name
cSize <- md %>% dplyr::count(seurat_clusters,Pregnancy_ID,Origin,Condition)

#filtering >
sSize <- cSize %>% filter(n>20) %>% dplyr::count(seurat_clusters,Origin,Condition) %>% filter(n>2)

kSize <- sSize %>% dplyr::count(seurat_clusters,Origin) %>% filter(n>=2) %>%
    mutate(cname=paste(seurat_clusters,Origin,sep="_"))

kSize


md<-md %>% mutate(cname=paste(seurat_clusters,Origin,sep="_"))
md<-md %>% mutate(GroupPregnancy_ID=paste(Group,Pregnancy_ID,sep="_"))

resList <- lapply(1:nrow(kSize),function(ii){
    ##ii <- 1
    cname1<-cname <- kSize$cname[ii]
    cat("#",cname,"\n")
    ##
    md_i <- md %>%
        filter(
            seurat_clusters==kSize$seurat_clusters[ii],
            Origin==kSize$Origin[ii]) 
    
    ## samples (control/case)
    bti <- md_i %>% transmute(bti=paste(Group,Pregnancy_ID,Library,sep="_")) %>% unlist %>% factor
    
    ## data selected
    all_i <- sc@assays$RNA@data[,rownames(md_i)]
    ##
    
    ##
    X <- model.matrix( ~ 0 + bti)
    qr.X <- qr(X)
    qr.X$rank
    dim(X)
    YtX <- all_i %*% X
    YtX <- as.matrix(YtX)
    dim(YtX)
    ##
    ## Running DEseq2 on the aggregate
    ## ========================================================
    ##
    bti2 <- gsub("bti","",colnames(YtX))
    colnames(YtX) <- bti2
    ##
    cmat<-YtX
    
     anno_i <- tibble(kbid=rownames(cmat),rs=rowSums(cmat),nz=rowSums(cmat>0)) %>%
        inner_join(anno %>%  dplyr::select(kbid,Chr,TSS,Strand,gene_name)) %>%
        filter(Chr %in% paste0("chr", 1:22), rs > 10, nz > 3) ## keep only autosomal
    ##
    ##
    table(anno_i$Chr)
    ##
    genesel <- (unique(anno_i$kbid))
    cmat <- cmat[genesel,]
    dim(cmat)
    ##
    ##Create sample table
    cn<-colnames(cmat)
    x<-strsplit(cn,"_")
    ##
    cvt <- data.frame(matrix(unlist(x), nrow=length(x), byrow=T),stringsAsFactors=FALSE)
    colnames(cvt)<-c("Group","Indiv","Library")
    ##
    cvt$Group = relevel(factor(cvt$Group),"Control")
    
    ##
    ##
    ## Running DEseq2 on the aggregate
    ## ========================================================
    ##
    dds <- DESeqDataSetFromMatrix(round(cmat),cvt, ~ Library+Group)    
    dds <- DESeq(dds,parallel=TRUE)
    ##
    ## save DEseq object.    
    fname <- paste0(outFolder,cname,"_dds_",Sys.Date(),".rds")
    fname
    write_rds(dds,fname)
    ##
    ## Parse the results 
    ## ========================================================
    ##
    res <- results(dds)
    myres <- as.data.frame(res) %>%
        rownames_to_column("kbid") %>%
        left_join(anno_i)
    myres$cname <- cname
    ##
    write_tsv(myres,paste0(outFolder,cname,".",Sys.Date(),".txt"))
    ##
    nres<-nrow(myres %>% filter(padj<.1,abs(log2FoldChange)>1))
    cat("# Total sig for",cname,": ",nres,"\n")
    myres
})

res <- do.call(rbind,resList)

sum(res$padj<0.1,na.rm=TRUE)



Sys.Date()

system(paste0("mkdir -p ",outFolder))

res %>% filter(padj<0.1,abs(log2FoldChange)>0) %>% dplyr::count(cname) %>%
    write_tsv(paste0(outFolder,"Summary.FDR.",Sys.Date(),".tsv"))

res %>% write_tsv(paste0(outFolder,"ALL.combined.",Sys.Date(),".tsv"))

res %>% filter(padj<0.1,abs(log2FoldChange)>0.0) %>%
    write_tsv(paste0(outFolder,"SIG.combined.2021-02-17.tsv"))

res %>% filter(padj<0.1,abs(log2FoldChange)>1.0) %>% select(kbid,gene_name,cname,log2FoldChange,padj) 

res<- read_tsv(paste0(outFolder,"ALL.combined.2021-02-17.tsv"))
res %>% filter(padj<0.1,abs(log2FoldChange)>0.5) %>% dplyr::count(cname) %>%
    write_tsv(paste0(outFolder,"Summary.FDR.FC.0.5_",Sys.Date(),".tsv"))

res %>% filter(padj<0.1,abs(log2FoldChange)>1) %>% dplyr::count(cname) %>%
    write_tsv(paste0(outFolder,"Summary.FDR.FC.1.0_",Sys.Date(),".tsv"))




#### barplot DE genes ###

outFolder <- paste0("7_outputs_DESeq_ConditionsByCluster/")
system(paste0("mkdir -p ",outFolder))

res <-read_tsv(paste0(outFolder,"ALL.combined.2021-02-17.tsv"))
res_up<-res %>% filter(padj<0.1,log2FoldChange>0) %>% dplyr::count(cname)
colnames(res_up)[2]<-"up-regulated"
res_down<-res %>% filter(padj<0.1,log2FoldChange<0) %>% dplyr::count(cname)
colnames(res_down)[2]<-"down-regulated"
res_down[,2]<--1*res_down[,2]
res_total<-res_down  %>% left_join(res_up)

clust2Name<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Name<-paste0(c(0:23),"_",clust2Name)
names(clust2Name)<-paste0(c(0:23),"_M")
res_total$sum<-res_total$`down-regulated`+res_total$`up-regulated`
rw<-as.character(clust2Name[res_total$cname])
rownames(res_total)<-rw
res_total[is.na(res_total)]<-0
rw<-clust2Name[res_total$cname]
res_total<-res_total[,c(-1,-4)]
res_total<-as.matrix(res_total)


library(reshape2)
library(ggplot2)
dat2 <- melt(res_total)
dat2$Var1<-as.character(rw)
colnames(dat2)<-c("cluster","DE","DEnumber")
#dat2<-dat2[order(dat2$cluster),]

names(clust2Name)<-c(0:23)

dat2$clusternumber<-as.numeric(sapply(dat2$cluster,function(x){return(unlist(strsplit(x,"_"))[1])}))
#dat2$cluster <- factor(dat2$cluster,levels=unique(dat2$cluster))
dat2<-dat2[order(dat2$clusternumber,decreasing = FALSE),]

fname=paste0(outFolder,"barplot-DE-up-down_v2.pdf");
pdf(fname,width=10,height=7)
ggplot(dat2, aes(x=reorder(cluster,-clusternumber), y=DEnumber,fill=DE))+#,color=mycolor)+#,fill=DE) +
    geom_bar(stat='identity') +
    theme_bw()+
    scale_fill_manual("legend",values=c("down-regulated"="#333399","up-regulated"="#A50021") )+
    coord_flip()+
    theme(legend.title=element_blank())+
    ylab("Number of differentially expressed genes")+
    xlab("")
dev.off() 

