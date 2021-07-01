
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


outFolder <- paste0("7_outputs_DESeq_ConditionsByCluster/")
system(paste0("mkdir -p ",outFolder))

## Load cells
sc <- read_rds("6_harmony_cellClass_plots_res0.8_final/SeuratObject.rds")

md <- sc@meta.data

# filter sample "HPL20874"


# md <- filter(md,!Pregnancy_ID %in% filter_sample)
# #md <- filter(md,!Pregnancy_ID %in% c("HPL20874","HPL20875"))
# 
# # 
# length(md$Pregnancy_ID=="HPL20874")
# 
# length(md$Pregnancy_ID %in% c("HPL20874","HPL20875"))

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

resList <- lapply(1:nrow(kSize),function(ii){
    ##ii <- 1
    cname <- kSize$cname[ii]
    cat("#",cname,"\n")
    ##
    md_i <- md %>%
        filter(
            #Location==kSize$Location[ii],
            seurat_clusters==kSize$seurat_clusters[ii],
            Origin==kSize$Origin[ii]) 
    
    ## samples (control/case)
    bti <- md_i %>% transmute(bti=paste(Group,Pregnancy_ID,sep="_")) %>% unlist %>% factor
    
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
    ##
    ##
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
    colnames(cvt)<-c("Group","Indiv")
    ##
    cvt$Group = relevel(factor(cvt$Group),"Control")
    ##
    ##
    ## Running DEseq2 on the aggregate
    ## ========================================================
    ##
    dds <- DESeqDataSetFromMatrix(round(cmat),cvt, ~ Group)    
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

#outFolder <- paste0("./7_outputs_DESeq_ConditionsByCluster_",Sys.Date(),"/")

Sys.Date()

#outFolder <- paste0("./7_outputs_DESeq_ConditionsByCluster_filtered_HPL20874_HPL20874_2020-11-27/")
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



res<-read_tsv(paste0(outFolder,"SIG.combined.2021-02-17.tsv"))
res<-res %>% filter(padj<0.1,abs(log2FoldChange)>0.5) %>% select(kbid,gene_name,cname,log2FoldChange,padj) 
res_oxy<-res[which(res$gene_name %in% c("OXT.L","OXT.S","AVP","LZTS3","OXTR","OXTRL","OXTR.L","OXTR.S","RAD18","VT1")),]    

oxy_genes<-c("OXT.L","OXT.S","AVP","LZTS3","OXTR","OXTRL","OXTR.L","OXTR.S","RAD18","VT1")
oxy_genes[which(oxy_genes%in% anno_i$gene_name)]
