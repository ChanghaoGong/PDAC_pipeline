

###
 # @Description: 
 # @Author: gongchanghao
 # @Date: 2023-09-01 16:09:51
 # @LastEditTime: 2023-09-20 20:06:01
 # @LastEditors: gongchanghao
 # @E-mail: gongchanghao@genomics.cn
 # @FilePath: /gongchanghao/script/spatial_pipeline/bin/cellbin_pipeline_v2/7.ADM/7.3_ADM_Cell_monocle2.R
 ###


###
 # @Description: 
 # @Author: gongchanghao
 # @Date: 2023-05-25 11:10:10
 # @LastEditTime: 2023-08-24 18:15:48
 # @LastEditors: gongchanghao
 # @E-mail: gongchanghao@genomics.cn
 # @FilePath: /gongchanghao/script/spatial_pipeline/bin/cellbin_pipeline_v2/5.CAF/5.7_CAF_scRNA_monocle2.R
 ###

library(Seurat)
library(ggplot2)
library(monocle)
library(patchwork)
library(dplyr)
library(magrittr)
library(ggsci)
args <- commandArgs(T)
sample <- args[1]
sample <- 'SS200000495BR_C5'
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/script/wangdi_code/SelfFunctions/colors.R")
wd <- sprintf("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/PDAC/spatial/cellbin_v3/7.3_ADM_region_monocle2/%s",sample)
dir.create(wd,recursive = T)
setwd(wd)
seurat_obj_all <- readRDS('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/PDAC/spatial/cellbin_v3/2.1_merge_big_cell_community/all_merged_adata.rds')
is_duct <- grep("^Duct_|Acinar_", seurat_obj_all$celltype)
seurat_obj <- seurat_obj_all[, (is_duct) ]
seurat_obj <- seurat_obj[,seurat_obj$batch==sample]
# seurat_obj <- readRDS(opt$infile)
Idents(seurat_obj) <- seurat_obj$community_type
removeBiasGenes <- function(mat){

  RPgenes <- rownames(mat)[intersect(grep("^RP", rownames(mat)), grep("-", rownames(mat)))]
  RPgenes2 <- rownames(mat)[grep("^RP[SL]", rownames(mat))]
  MTgenes <- rownames(mat)[grep("^MT-", rownames(mat))]
  CTCgenes <- rownames(mat)[intersect(grep("^CTC", rownames(mat)), grep("-", rownames(mat)))]
  MIRgenes <- rownames(mat)[grep("^MIR", rownames(mat))]
  ACgenes <- rownames(mat)[intersect(grep("^AC[0-9]", rownames(mat)), grep(".", rownames(mat)))]
  CTgenes <- rownames(mat)[intersect(grep("^CT", rownames(mat)), grep("-", rownames(mat)))]
  LINCgenes <- rownames(mat)[grep("^LINC[0-9]", rownames(mat))]
  ALgenes <- rownames(mat)[intersect(grep("^AL", rownames(mat)), grep(".", rownames(mat)))]

  rmgenes <- c(RPgenes, RPgenes2, MTgenes, CTCgenes, MIRgenes, ACgenes, CTgenes, LINCgenes, ALgenes)

  # datacount <- mat[!rownames(mat)%in%rmgenes,]
  # datacount <- datacount[rowSums(datacount > 0) > 1,]
  return(rmgenes)
}
rmgenes <- removeBiasGenes(seurat_obj)
filetered_genes <- rownames(seurat_obj)[!rownames(seurat_obj)%in%rmgenes]
seurat_obj <- subset(seurat_obj, features = filetered_genes)
seurat_obj <- NormalizeData(seurat_obj, assay = "Spatial", verbose = FALSE)
# 计算每个细胞类型的细胞数
celltype_counts <- table(seurat_obj$community_type)
# 获取细胞数大于等于20的细胞类型
celltype_to_keep <- names(celltype_counts[celltype_counts >= 50])
# 过滤细胞数大于等于20的细胞类型
seurat_obj <- subset(seurat_obj, subset = community_type %in% celltype_to_keep)
# DEGs
DEG_wilcox <- FindAllMarkers(seurat_obj, min.pct = 0.1, logfc.threshold = 0.2, test.use = "wilcox")
DEG_wilcox_up <- subset(DEG_wilcox,avg_log2FC > 0.263 & p_val_adj < 0.05)
DEG_wilcox_up <- unique(DEG_wilcox_up$gene)

data <- GetAssayData(seurat_obj, assay = 'Spatial', slot = 'counts')
pd <- new('AnnotatedDataFrame', data = seurat_obj@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new("AnnotatedDataFrame", data = fData)

#创建monocle对象cds
cd <- newCellDataSet(data,
                     phenoData = pd,
                     featureData = fd,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())

#Estimate size factors and dispersions
cd <- estimateSizeFactors(cd)
cd <- estimateDispersions(cd)

#Select genes that differ between clusters/stages
diff_test_res <- differentialGeneTest(cd, fullModelFormulaStr = "~community_type", cores = 20)

# cd <- detectGenes(cd, min_expr = 0.1)
# fData(cd)$use_for_ordering <-
#     fData(cd)$num_cells_expressed > 0.05 * ncol(cd)
# cd <- reduceDimension(cd,
#                               max_components = 2,
#                               norm_method = 'log',
#                               num_dim = 3,
#                               reduction_method = 'tSNE',
#                               verbose = T)
# cd <- clusterCells(cd, verbose = F)
# pdf("plot_cell_clusters.pdf")
# plot_cell_clusters(cd, color_by = 'as.factor(Cluster)')
# plot_cell_clusters(cd, color_by = 'as.factor(celltype)')
# dev.off()
# cd_expressed_genes <-  row.names(subset(fData(cd),
# num_cells_expressed >= 10))
# pdf("plot_rho_delta.pdf")          
# plot_rho_delta(cd, rho_threshold = 2, delta_threshold = 3 )
# dev.off()
# cd <- clusterCells(cd,
#                  rho_threshold = 2,
#                  delta_threshold = 3,
#                  skip_rho_sigma = T,
#                  verbose = F)

# pdf("plot_cell_clusters.pdf")
# plot_cell_clusters(cd, color_by = 'as.factor(Cluster)')
# plot_cell_clusters(cd, color_by = 'as.factor(celltype)')
# dev.off()
# clustering_DEG_genes <-
#     differentialGeneTest(cd[cd_expressed_genes,],
#           fullModelFormulaStr = '~Cluster',
#           cores = 20)

# ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
length(ordering_genes)
# ordering_genes <- row.names(subset(diff_test_res, qval < 0.05))
# if(length(ordering_genes) > 1000){
#     ordering_genes <- row.names(subset(diff_test_res, qval < 0.05))[1:1000]
# }else{
#     ordering_genes <- row.names(subset(diff_test_res, qval < 0.05))
# }
ordering_genes <- intersect(ordering_genes, DEG_wilcox_up)
# ordering_genes <- VariableFeatures(seurat_obj)
cd <- setOrderingFilter(cd, ordering_genes)
pdf("plot_ordergenes_cluster.pdf")
plot_ordering_genes(cd)
dev.off()

#Order cells by progress
cd <- reduceDimension(cd, max_components = 2, method = 'DDRTree')
cd <- orderCells(cd, reverse = FALSE)
#调整细胞排序方向
# cd <- orderCells(cd, reverse = TRUE)
#保存结果
saveRDS(cd, 'ADM_monocle2.rds')

#排序好的细胞可以进行可视化，可标注细胞的各注释信息

rename_vector <- c(
  "Duct1_NTRK2_NMF_2" = "Duct_NTRK2",
  "Acinar_NMF_8" = "Acinar_NMF_8",
  "Duct_MMP7_NMF_10" = "Duct_MMP7",
  "Acinar_REG3G_NMF_14" = "Acinar_REG3G",
  "Duct2_CEACAM6_NMF_7" = "Duct_like_CEACAM6",
  "Acinar_NMF_15" = "Acinar_NMF_15",
  "Acinar_RBPJL_NMF_11" = "Acinar_RBPJL",
  "Endocrine_NMF_12" = "Endocrine",
  "Duct_FTH1_NMF_6" = "Duct_FTH1",
  "Duct_FOSB_NMF_5" = "Duct_FOSB"
)
region_map <- c(
    "Tumor_Cell_region_1"= 'TCR1',
    'Tumor_Cell_region_2'= 'TCR2',
    'Tumor_Cell_region_3'= 'TCR3',
    'Tumor_Cell_region_4'= 'TCR4',
    'Tumor_around_region_1'= 'TAR1',
    'Tumor_around_region_2'= 'TAR2',
    'Tumor_around_region_3'= 'TAR3',
    'Tumor_around_region_4'= 'TAR4',
    'Tumor_around_region_5'= 'TAR5',
    'Mesenchymal_Cell_region_1'= 'MCR1',
    'Mesenchymal_Cell_region_2'= 'MCR2',
    'Mesenchymal_Cell_region_3'= 'MCR3',
    'ADM_region_1'= 'ADR1',
    'ADM_region_2'= 'ADR2',
    'ADM_region_3'= 'ADR3',
    'ADM_region_4'= 'ADR4',
    'Acinar_Cell_region'='ACR',
    'Immune_Cell_region'='ICR'
)
cd$cellsubtype <- rename_vector[cd$cellsubtype]
cd$region <- region_map[cd$community_type]

pdf("cell_track.pdf",width=7,height=5)
plot_cell_trajectory(cd, color_by = "celltype")+
  theme(legend.position="right", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12, angle = 90)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12, colour = "black")) +
  theme(axis.text.x = element_text(size = 12, colour = "black")) +
  theme(legend.position = "right", 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+
  scale_color_manual(values=col16)
plot_cell_trajectory(cd, color_by = "cellsubtype")+
  theme(legend.position="right", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12, angle = 90)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12, colour = "black")) +
  theme(axis.text.x = element_text(size = 12, colour = "black")) +
  theme(legend.position = "right", 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+
  scale_color_manual(values=col16)
plot_cell_trajectory(cd, color_by = "community_type") +
  theme(legend.position="right", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12, angle = 90)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12, colour = "black")) +
  theme(axis.text.x = element_text(size = 12, colour = "black")) +
  theme(legend.position = "right", 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+
  scale_color_manual(values=col16)
plot_cell_trajectory(cd, color_by = "community_type") + facet_wrap(~community_type, nrow = 3) + scale_color_manual(values=col16)
plot_cell_trajectory(cd, color_by = "State") +
  theme(legend.position="right", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12, angle = 90)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12, colour = "black")) +
  theme(axis.text.x = element_text(size = 12, colour = "black")) +
  theme(legend.position = "right", 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+
  scale_color_npg()
plot_cell_trajectory(cd, color_by = "State") + facet_wrap(~State, nrow = 1) + scale_color_npg()
plot_cell_trajectory(cd, cell_size = 4, color_by = "Pseudotime")
plot_complex_cell_trajectory(cd, color_by = 'State', show_branch_points = T, cell_size = 2, cell_link_size = 0.5) + scale_color_npg()
plot_complex_cell_trajectory(cd, color_by = 'community_type', show_branch_points = T, cell_size = 2, cell_link_size = 0.5) + scale_color_manual(values=col16)
dev.off()

cols <- col16[c(1,2,5,6,7,8,9)]
pdf("cell_track_subtype_region.pdf",width=6,height=4)
plot_cell_trajectory(cd, color_by = "cellsubtype") +
  theme(legend.position="right", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12, angle = 90)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12, colour = "black")) +
  theme(axis.text.x = element_text(size = 12, colour = "black")) +
  theme(legend.position = "right", 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+
  scale_color_manual(values=cols)
plot_cell_trajectory(cd, color_by = "region") +
  theme(legend.position="right", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12, angle = 90)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12, colour = "black")) +
  theme(axis.text.x = element_text(size = 12, colour = "black")) +
  theme(legend.position = "right", 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+
  scale_color_manual(values=getDefaultColors(n=6,type=6))
plot_cell_trajectory(cd, color_by = "State") +
  theme(legend.position="right", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12, angle = 90)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12, colour = "black")) +
  theme(axis.text.x = element_text(size = 12, colour = "black")) +
  theme(legend.position = "right", 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+
  scale_color_npg()
plot_cell_trajectory(cd, cell_size = 4, color_by = "Pseudotime")
dev.off()

pdf("cell_track_subtype_fct.pdf",width=7,height=7)
plot_cell_trajectory(cd, color_by = "cellsubtype") + 
    facet_wrap(~cellsubtype, nrow = 3) +   
    theme(legend.position = "none") + 
    scale_color_manual(values=cols)
plot_cell_trajectory(cd, color_by = "region") + 
    facet_wrap(~region, nrow = 3) +   
    theme(legend.position = "none") + 
    scale_color_manual(values=getDefaultColors(n=6,type=6))
dev.off()

#鉴定随时间变化的基因Finding genes that change as a function of pseudotime
diff_test_res <- differentialGeneTest(cd,fullModelFormulaStr ="~sm.ns(Pseudotime)", cores = 40)
diff_test_res <- diff_test_res[order(diff_test_res$qval),]
saveRDS(diff_test_res, 'diff_gene_over_time.rds')
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01))
pdf("pseudotime_diff_gene_heatmap_top100.pdf",height=9)
p <- plot_pseudotime_heatmap2(cd[sig_gene_names[1:100],],
                num_clusters = 5,
                cores = 20,
                show_rownames = T)
print(p)
dev.off()

   


# 绘制pseudotime热图+富集分析合图
n_cluster = 4 # number of clusters
n_top_genes = 5 # number of top genes to show in each cluster
n_top_pathways = 5 # number of top pathways to show in each cluster
library(ClusterGVis)
# plot_pseudotime_heatmap2(cd[sig_gene_names,],
#                          num_clusters = n_cluster,
#                          cores = 20,
#                          show_rownames = T,
#                          return_heatmap = T)
p <- plot_pseudotime_heatmap2(cd[sig_gene_names,],
                         num_clusters = n_cluster,
                         cores = 20,
                         show_rownames = T,
                         return_heatmap = F)

library(org.Hs.eg.db)
library(msigdbr)
library(clusterProfiler)
m_df <- msigdbr(species = "Homo sapiens", category = "H")%>% 
      dplyr::select(gs_name, gene_symbol)
m_df$gs_name <- gsub("HALLMARK_","",m_df$gs_name)
m_df_subset <- m_df[m_df$gene_symbol %in% rownames(seurat_obj),]
enrich.data <- p$wide.res
gc_hyper <- list()
for(i in unique(enrich.data$cluster)){
    gc_hyper[[i]] <- rownames(enrich.data)[enrich.data$cluster==i]
}
ehm <- compareCluster(gc_hyper, fun = "enricher", TERM2GENE = m_df, pvalueCutoff = 0.05)
write.table(ehm,'enrich_hm.txt',sep='\t',col.names=NA)     


# do GO and HALLMARK enrichment
library(org.Hs.eg.db)
library(msigdbr)
m_df <- msigdbr(species = "Homo sapiens", category = "H")%>% 
      dplyr::select(gs_name, gene_symbol)
m_df$gs_name <- gsub("HALLMARK_","",m_df$gs_name)
m_df_subset <- m_df[m_df$gene_symbol %in% rownames(seurat_obj),]
enrich_HM <- enrichCluster(object = p,
                        type = "ownSet",
                        id.trans = FALSE,
                        TERM2GENE = m_df_subset,
                        TERM2NAME = NA,
                        pvalueCutoff = 1,
                        topn = n_top_pathways,
                        seed = 42)
enrich_GO <- enrichCluster(object = p,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = n_top_pathways,
                        seed = 42)

# get top n DEGs in each cell state              
seurat_obj$mn2_state <- cd$State
Idents(seurat_obj) <- seurat_obj$mn2_state
seurat_obj_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,group.by='mn2_state',test.use='wilcox')
filtered_markers <- seurat_obj_markers[seurat_obj_markers$p_val_adj < 0.05, ]
filtered_markers %>%
    group_by(cluster) %>%
    top_n(n = n_top_genes, wt = avg_log2FC) -> topn

# get a color vector according to enrich_HM$group
colors = jjAnno::useMyCol("calm",n = length(unique(enrich_HM$group)))
go.col = c()
i = 1
for(n in table(enrich_HM$group)){
    go.col = c(go.col,rep(colors[i],n))
    i = i + 1
}
pdf("pseudotime_heatmap_HM.pdf",width=12,height=7)                       
visCluster(object = p,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = unique(topn$gene),
           markGenes.side = "left",
           annoTerm.data = enrich_HM,
           go.col = go.col,
           add.bar = T,
           line.side = "left")
dev.off()
colors = jjAnno::useMyCol("calm",n = length(unique(enrich_GO$group)))
go.col = c()
i = 1
for(n in table(enrich_GO$group)){
    go.col = c(go.col,rep(colors[i],n))
    i = i + 1
}
pdf("pseudotime_heatmap_GO_BP.pdf",width=12,height=7)                       
visCluster(object = p,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = unique(topn$gene),
           markGenes.side = "left",
           annoTerm.data = enrich_GO,
           go.col = go.col,
           add.bar = T,
           line.side = "left")
dev.off()


cd <- readRDS('ADM_monocle2.rds')
table(cd$State)
seurat_obj$mn2_state <- cd$State
pdf(sprintf("mn2_state_spatial1_%s.pdf",sample),width=7,height=5)
print(SpatialDimPlot(seurat_obj, label = FALSE, label.size = 3, stroke = 0,pt.size.factor = 150, group.by = 'mn2_state')+scale_fill_npg())
dev.off()

Idents(seurat_obj) <- seurat_obj$mn2_state
seurat_obj_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,group.by='mn2_state',test.use='wilcox')
filtered_markers <- seurat_obj_markers[seurat_obj_markers$p_val_adj < 0.05, ]
filtered_markers %>%
    group_by(cluster) %>%
    top_n(n = 200, wt = avg_log2FC) -> top200
# get a df for each cluster in top200 as a colomn
dftop200 <- top200[,c('cluster','gene')] %>%
  pivot_wider(names_from = cluster, values_from = gene,names_repair='minimal')
dftop200 <- as.data.frame(dftop200)
dftop200 <- lapply(dftop200, function(col) col[[1]])
output_file <- "DEGs_mn2_state.txt"
# 打开文件以写入数据
file_conn <- file(output_file, "w")
# 遍历每个基因列表，写入文件
for (i in 1:length(dftop200)) {
  gene_list <- dftop200[[i]]
  gene_list_name <- names(dftop200)[i]  
  # 将基因列表的名称和基因列表内容连接成一个字符串
  gene_list_entry <- paste(gene_list_name, paste(gene_list, collapse = ","), sep = "\t")
  
  # 写入文件
  cat(gene_list_entry,sep='',file = file_conn, "\n")
}

# 关闭文件连接
close(file_conn)



library(clusterProfiler)
library(msigdbr)
m_df <- msigdbr(species = "Homo sapiens", category = "H")%>% 
      dplyr::select(gs_name, gene_symbol)
m_df$gs_name <- gsub("HALLMARK_","",m_df$gs_name)

m_df_subset <- m_df[m_df$gene_symbol %in% rownames(seurat_obj),]

ehm <- compareCluster(dftop200, fun = "enricher", TERM2GENE = m_df_subset, pvalueCutoff = 0.05)
write.table(as.data.frame(ehm), "compare_hyper_HALLMARK.txt", sep = "\t", row.names = FALSE)
p <- dotplot(ehm,showCategory = 20)
cluster_factor <- p$data$Cluster
cleaned_levels <- gsub("\n.*", "", levels(cluster_factor))
cleaned_factor <- factor(gsub("\n.*", "", cluster_factor),levels=cleaned_levels)
p$data$Cluster <- cleaned_factor
p <- p + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
pdf("compare_hyper_HALLMARK_dotplot.pdf", width = 7, height = 10)
print(p)
dev.off()
library(org.Hs.eg.db)
for (term in c("CC", "BP", "MF")) {
# for (term in c("BP")) {
    ego <- compareCluster(dftop200,
        fun = 'enrichGO',
        universe = rownames(seurat_obj),
        keyType = "SYMBOL",
        OrgDb = org.Hs.eg.db,
        ont = term,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
    )
    ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
    write.table(as.data.frame(ego2), sprintf("compare_hyper_GO_%s.txt",term), sep = "\t", row.names = FALSE)
    p <- dotplot(ego2,showCategory = 20)
    cluster_factor <- p$data$Cluster
    cleaned_levels <- gsub("\n.*", "", levels(cluster_factor))
    cleaned_factor <- factor(gsub("\n.*", "", cluster_factor),levels=cleaned_levels)
    p$data$Cluster <- cleaned_factor
    p <- p + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
    pdf(sprintf("compare_hyper_GO_%s_dotplot.pdf",term), width = 7, height = length(unique(p$data$Description))*0.6)
    print(p)
    dev.off()
}

