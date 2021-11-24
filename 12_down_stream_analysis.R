library(UpSetR)
library(tidyverse)
library(ComplexHeatmap)
library(mashr)

###################################################
### Multivariate Adaptive Shrinkage

###################################################


library(rmeta)
library(ashr)
library(mashr)
library(pheatmap)

#outFolder <- paste0("12_downstream_analysis/")
outFolder <- paste0("12_downstream_analysis_batch_corrected/")
system(paste0("mkdir -p ",outFolder))


#res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/SIG.combined.2021-02-17.tsv")
res <- read_tsv("7_outputs_DESeq_ConditionsByCluster_bath_library/SIG.combined.2021-10-18.tsv")
res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))
clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(c(0:23),"_",clust2Names)
names(clust2Names)<-c(0:23)

cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D", "#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","#2C3E18","#745B48", "#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-clust2Names

res$cluster_colors<-cluster.Colors[res$Cell_type]

res$Cell_type<-clust2Names[res$Cell_type]



################################################ 
### upset plot
################################################
res<-res %>% filter(abs(log2FoldChange)>0 & padj< 0.1) #%>% select(gene_name,Cell_type)

res <- res%>% mutate (DE=ifelse(abs(log2FoldChange)>0 & padj< 0.1 ,1,0))

res_df<-matrix(0,nrow=length(unique(res$gene_name)),ncol=length(unique(res$Cell_type)))
colnames(res_df)<-unique(res$Cell_type)
rownames(res_df)<-unique(res$gene_name)
cell_types<-unique(res$Cell_type)
for (i in 1:length(cell_types))
{
  res_celltype<-res %>% filter(Cell_type %in% cell_types[i]) %>% dplyr::select(gene_name)
  gn<-res_celltype$gene_name
  res_df[gn,cell_types[i]]<-1
}
res_df<-as.data.frame(res_df)
mt = make_comb_mat(res_df)#, top_n_sets = 25)
m2 <- mt[,comb_size(mt)>20]
fig0 <- UpSet(m2, set_order=colnames(res_df),comb_order=order(comb_size(m2)))#sets.bar.color = cluster.Colors[colnames(res_df)]
pdf(file=paste0(outFolder,"upsetplot_comb_size25.pdf"))#,width=45,height=20)# or other device
fig0
dev.off()


m2 <- mt[,comb_size(mt)>10]
fig0 <- UpSet(m2, set_order=colnames(res_df),comb_order=order(comb_size(m2)))#sets.bar.color = cluster.Colors[colnames(res_df)]
pdf(file=paste0(outFolder,"upsetplot_comb_size10.pdf"))#,width=45,height=20)# or other device
fig0
dev.off()


##############

set.seed(1)



########################################################################
#load data 
########################################################################
#res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/SIG.combined.2021-02-17.tsv")
#res <- read_tsv("7_outputs_DESeq_ConditionsByCluster/ALL.combined.2021-02-17.tsv")
res <- read_tsv("7_outputs_DESeq_ConditionsByCluster_bath_library/ALL.combined.2021-10-18.tsv")

res <- res %>% separate(cname,c("Cell_type","Origin"),sep="_",remove=FALSE)
res <- res %>% filter(!is.na(pvalue))
clust2Names<-c("Stromal-1","Macrophage-2","Macrophage-1","Endothelial-1","Monocyte","CD4_T-cell","Decidual","CD8_T-cell","LED","Stromal-2","ILC","NK-cell","Smooth muscle cells-1","Myofibroblast","Macrophage-3","Endothelial-2","DC","Smooth muscle cells-2","EVT","Plasmablast","Smooth muscle cells-3","Macrophage-4","B-cell","Unciliated Epithelial")
clust2Names<-paste0(c(0:23),"_",clust2Names)
names(clust2Names)<-c(0:23)

cluster.Colors<-c("#DF7D99","#838EDF","#4E65A6","#FFC000","#2BA3D3","#9ABF5C","#D14357","#329B2D","#D5438E","#ED4315","#76956C","#7BC791","#CA8588","#F88091","#72C6C8","#E4652C","#9B91B9","#A37584","2C3E18","#745B48","#AA5485","#4E747A","#C59A89","#C9C76F")   
names(cluster.Colors)<-clust2Names
res$cluster_colors<-cluster.Colors[res$Cell_type]
res$Cell_type<-clust2Names[res$Cell_type]





########################################################################
#### mash analysis 
########################################################################

hist(table(res$kbid))
max(table(res$kbid))
sum(table(res$kbid)==16)
selgenes <- names(which(table(res$kbid)==16))

Bhat<-res %>% filter(kbid %in% selgenes) %>% dplyr::select(kbid,Cell_type,log2FoldChange) %>%
  pivot_wider(names_from = Cell_type, values_from = log2FoldChange)
rw<-Bhat$kbid
Bhat<-as.matrix(Bhat)
rownames(Bhat)<-Bhat[,"kbid" ]
Bhat<-Bhat[,-1]
Bhat<-apply(Bhat,c(1,2),as.numeric)

Shat<-res %>% filter(kbid %in% selgenes) %>% dplyr::select(kbid,Cell_type,lfcSE) %>%
  pivot_wider(names_from = Cell_type, values_from = lfcSE)

rw<-Shat$kbid
Shat<-as.matrix(Shat)
rownames(Shat)<-Shat[,"kbid" ]
Shat<-Shat[,-1]
Shat<-apply(Shat,c(1,2),as.numeric)

########################################################################
# data
########################################################################
data = mash_set_data(Bhat, Shat)


########################################################################
#covariance matrices
########################################################################
U.c = cov_canonical(data)
print(names(U.c))


########################################################################
#fit the model
########################################################################
m.c = mash(data, U.c)

write_rds(m.c, paste0(outFolder,"m.c.rds"))
########################################################################
#Extract Posterior Summaries
########################################################################

#local false sign rate

lfsr_m.c<-get_lfsr(m.c)

head(get_lfsr(m.c))

write_rds(lfsr_m.c, paste0(outFolder,"lfsr_m.c.rds"))


# comparison based on negative and positive effects 
######################
PosteriorMean<-m.c$result$PosteriorMean
PosteriorMean_positive<-rownames(PosteriorMean)[which(PosteriorMean>0)]
PosteriorMean_negative<-rownames(PosteriorMean)[which(PosteriorMean<0)]

outFolder2 <- paste0("12_mash_analysis/")
system(paste0("mkdir -p ",outFolder2))
write_rds(m.c, paste0(outFolder2,"m.c.rds"))
write_rds(PosteriorMean_positive, paste0(outFolder2,"PosteriorMean_positive.rds"))
write_rds(PosteriorMean_negative, paste0(outFolder2,"PosteriorMean_negative.rds"))




##################################################################################
# wider matrix generation

selgenes <- names(which(table(res$kbid)==16))
padj<-res %>% filter(!is.na(pvalue) & kbid %in% selgenes) %>% select(kbid,Cell_type,padj) %>%
  pivot_wider(names_from = Cell_type, values_from = padj)
rw<-padj$kbid
padj<-as.matrix(padj)
rownames(padj)<-padj[,"kbid" ]
padj<-padj[,-1]
padj<-apply(padj,c(1,2),as.numeric)
padj[padj>=0.1]<--1
padj[padj<=0.1 & padj>=0]<-1
padj[padj==-1]<-0
padj<-na.omit(padj)
DEgene_cluster<-padj
hist(rowSums(DEgene_cluster))


#############################################
## plots from result of mash analysis
#############################################

## heatmap 

pairwise_sharing_data<-get_pairwise_sharing(m.c)

fname=paste0(outFolder,"heatmap_pairwise_sharing.pdf");
pdf(fname,width=7,height=7)
pheatmap(pairwise_sharing_data,cluster_rows=TRUE,cluster_cols=TRUE,scale="none")
dev.off()

#############################################
### upset plot 
#############################################

lfsr<- get_lfsr(m.c)
lfsr[which(lfsr<0.1)]<--1
lfsr[which(lfsr>=0.1)]<-0
lfsr[which(lfsr==-1)]<-1

lfsr<-as.data.frame(lfsr)
mt = make_comb_mat(lfsr,complement_size = 0)#, top_n_sets = 25)

m2 <- mt[,comb_size(mt)>10]
#m2 <- mt[comb_degree(mt)>0]
#fig0 <- UpSet(m2, set_order=colnames(lfsr),comb_order=order(comb_size(m2)))#sets.bar.color = cluster.Colors[colnames(res_df)]
fig0 <- UpSet(m2,set_order=colnames(lfsr),comb_order=order(comb_size(m2)))
pdf(file=paste0(outFolder,"upsetplot_lfsr_comb_size25.pdf"))
fig0
dev.off()



##########################################################################################
#meta plot 
##########################################################################################



#outFolder<-"12_downstream_analysis/"
outFolder<-"12_downstream_analysis_batch_corrected/"
# load lfsr_m.c and m.c from mashr result
lfsr_m.c<-read_rds(paste0(outFolder,"lfsr_m.c.rds"))
m.c<-read_rds(paste0(outFolder,"m.c.rds"))


celltype<-"4_Monocyte"  
otherDEGs<-res$gene_name [which(res$Cell_type !=celltype)]
res_specific<-res %>% filter (Cell_type==celltype & !gene_name %in% otherDEGs)
res_specific<-res_specific %>% arrange(-log2FoldChange)
sample_genes<-res$kbid[which(res$gene_name %in%  res_specific$gene_name[1:20] )]




# meta plot
plot_func<-function(m.c,mytop,k=1,outFolder)
{
  for ( k in 1:length(which(mytop)))
  {
    par(mar=c(2, 1 ,4 ,3))
    i<-which(mytop)[k]
    plot.title<-res$gene_name[which(res$kbid==names(i))[1]]
    print(plot.title)
    system(paste0("mkdir -p ",outFolder))
    fname=paste0(outFolder,plot.title,".pdf");
    pdf(fname,width=10,height=7)
    
    print(plot.title)
    
    metaplot(get_pm(m.c)[i,],get_psd(m.c)[i,],colors = meta.colors(box = as.character(cluster.Colors))
             ,xlim = c(-1,1),xlab = "",ylab = "")
    #legend("topleft",legend=as.character(unique(res$Cell_type)), fill=cluster.Colors[unique(res$Cell_type)],bty = "n",cex=0.8) #
    title(plot.title)
    dev.off()
  }
  
  
}

mytop<-rep(FALSE,length(rownames(lfsr_m.c)))
names(mytop)<-rownames(lfsr_m.c)
mytop[which(names(mytop) %in% unique(sample_genes))]<-TRUE


system(paste0("mkdir -p ",paste0(outFolder,celltype,"/")))
# call to generate metaplots for all genes TRUE in mytop:
plot_func(m.c,mytop,outFolder=paste0(outFolder,celltype,"/"))




##################################################################################

#############################################
# load ref data 
#############################################
load_ref_data<-function(outFolder,fl)
{
  
  ref_data <- read.csv(paste0(fl,".csv"),stringsAsFactors = FALSE)
  ref_data<-ref_data %>% select(SYMBOL,logFC,P.Value,adj.P.Val,ENTREZ,t)
  colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","ENTREZID","Rt")
  #colnames(ref_data) <- c("R.gene_name","R.Log2FC","Rpadj","ENTREZID","Rt")
  # colnames(ref_data)<-c("R.gene_name","R.Log2FC","Rpvalue","Rpadj","Rt")
  ref_data$ENTREZID<-as.character(ref_data$ENTREZID)
  return(ref_data)
}    


#experiment<-"myometrium_term_TL-TNL_ALLList"
experiment<-"TLvsTNL_blood_ENTREZ"

ref_data<-load_ref_data(fl=experiment) 
ref_data2<-load_ref_data(fl="myometrium_term_TL-TNL_ALLList") 
colnames(ref_data)[which(colnames(ref_data)=="R.gene_name")]<-"gene_name"
ref_data2 <- ref_data %>% filter(Rpadj<0.1 &  !is.na(R.Log2FC) & !is.na(Rpadj))

res4 <- res %>% filter(padj<0.1 &  !is.na(log2FoldChange) & !is.na(padj))

resJoin<-res4 %>% inner_join(ref_data2)


resJoin<-res4 %>% inner_join(ref_data)
resJoin<-resJoin[!duplicated(resJoin[ , c("ENTREZID")]),]
resJoin_stratified<-resJoin%>% mutate(Rpadj.stratified = p.adjust(Rpvalue,"fdr"))

resJoin_stratified <- resJoin_stratified %>% filter(Rpadj.stratified<0.1 & !is.na(padj))

resJoin_stratified_ref<-resJoin_stratified %>% dplyr::select(gene_name,R.Log2FC,Rpadj,Rpadj.stratified,ENTREZID)

resJoin_stratified_full_list<-res4 %>% inner_join(resJoin_stratified_ref)

resJoin_stratified_full_list<-resJoin_stratified_full_list %>% dplyr::select(kbid,ENTREZID,gene_name, Cell_type, log2FoldChange,lfcSE,pvalue,padj,R.Log2FC,Rpadj,Rpadj.stratified)


write.csv( resJoin_stratified_full_list, file="12_downstream_analysis/TLvsTNL_blood_ENTREZ/resJoin_stratified_full_list_v2.csv")

resJoin<-resJoin_stratified


gene_select<-unique(resJoin_stratified$gene_name)
mytop<-rep(FALSE,length(rownames(lfsr_m.c)))
names(mytop)<-rownames(lfsr_m.c)
sample_genes<-res$kbid[which(res$gene_name %in% unique(gene_select)) ]
mytop[which(names(mytop) %in%sample_genes)]<-TRUE


gene_select<-unique(resJoin_stratified$gene_name)

system(paste0("mkdir -p ",paste0(outFolder,"smc-1/")))

plot_func(m.c,mytop,outFolder=paste0("mkdir -p ",paste0(outFolder,"smc-1/")))

system(paste0("mkdir -p ",paste0(outFolder,experiment,"/stratified_v2/")))
plot_func(m.c,mytop,outFolder=paste0(outFolder,experiment,"/stratified_v2/"))

