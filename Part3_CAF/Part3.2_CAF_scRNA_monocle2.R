

###
 # @Description: 
 # @Author: gongchanghao
 # @Date: 2023-05-25 11:10:10
 # @LastEditTime: 2023-10-18 11:31:59
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
# sample <- args[1]
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/script/wangdi_code/SelfFunctions/colors.R")
wd <- sprintf("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/PDAC/spatial/cellbin_v3/5.7_CAF_scRNA_monocle2/")

dir.create(wd,recursive = T)
setwd(wd)
seurat_obj_all <- readRDS('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/yixianai/SC_reanalysis/result/3_cluster/stromal/seurat_obj.rds')
is_stromal <- grepl("^Fibro_|PVL_", seurat_obj_all$cellsubtype)
is_tumor <- grepl("^Fibro_CLDN18", seurat_obj_all$cellsubtype)
is_csCAF <- grepl("^Fibro_CFD|Fibro_APOD|Fibro_USP53", seurat_obj$cellsubtype)
is_myCAF <- grepl("^Fibro_COL11A1|Fibro_F2RL2|Fibro_MME", seurat_obj$cellsubtype)
seurat_obj <- seurat_obj_all[, (is_stromal) & (!is_tumor)]
seurat_obj <- seurat_obj[, is_myCAF]
# seurat_obj <- readRDS(opt$infile)
Idents(seurat_obj) <- seurat_obj$cellsubtype
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

# DEGs
DEG_wilcox <- FindAllMarkers(seurat_obj, min.pct = 0.1, logfc.threshold = 0.2, test.use = "wilcox")
DEG_wilcox_up <- subset(DEG_wilcox,avg_log2FC > 0.2 & p_val_adj < 0.05)
DEG_wilcox_up <- unique(DEG_wilcox_up$gene)

data <- GetAssayData(seurat_obj, assay = 'RNA', slot = 'counts')
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
diff_test_res <- differentialGeneTest(cd, fullModelFormulaStr = "~cellsubtype", cores = 20)

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
# plot_cell_clusters(cd, color_by = 'as.factor(cellsubtype)')
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
# plot_cell_clusters(cd, color_by = 'as.factor(cellsubtype)')
# dev.off()
# clustering_DEG_genes <-
#     differentialGeneTest(cd[cd_expressed_genes,],
#           fullModelFormulaStr = '~Cluster',
#           cores = 20)

# ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

ordering_genes <- row.names(subset(diff_test_res, qval < 0.05))
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
saveRDS(cd, 'CAF_monocle2.rds')

#排序好的细胞可以进行可视化，可标注细胞的各注释信息

pdf("cell_track.pdf",width=7,height=5)
plot_cell_trajectory(cd, color_by = "sample")+
  theme(legend.position="right", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12, angle = 90)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12, colour = "black")) +
  theme(axis.text.x = element_text(size = 12, colour = "black")) +
  theme(legend.position = "right", 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))
plot_cell_trajectory(cd, color_by = "cellsubtype") +
  theme(legend.position="right", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12, angle = 90)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12, colour = "black")) +
  theme(axis.text.x = element_text(size = 12, colour = "black")) +
  theme(legend.position = "right", 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+
  scale_color_manual(values=col16)
plot_cell_trajectory(cd, color_by = "cellsubtype") + facet_wrap(~cellsubtype, nrow = 3) + scale_color_manual(values=col16) + theme(legend.position="none")
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
plot_complex_cell_trajectory(cd, color_by = 'cellsubtype', show_branch_points = T, cell_size = 2, cell_link_size = 0.5) + scale_color_manual(values=col16)
dev.off()

genes <- c('C7','CXCL12','LAMB1','COL15A1','')
cds_subset <- cd[genes,]
pdf('plot_genes_in_pseudotime.pdf',width=7,height=5)
plot_genes_in_pseudotime(cds_subset, color_by = "cellsubtype")
plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
plot_genes_in_pseudotime(cds_subset, color_by = "State")
dev.off()

#鉴定随时间变化的基因Finding genes that change as a function of pseudotime
diff_test_res <- differentialGeneTest(cd,fullModelFormulaStr ="~sm.ns(Pseudotime)", cores = 40)
diff_test_res <- diff_test_res[order(diff_test_res$qval),]
saveRDS(diff_test_res, 'diff_gene_over_time.rds')
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01))
pdf("pseudotime_diff_gene_heatmap_top100.pdf",height=9)
p <- plot_pseudotime_heatmap(cd[sig_gene_names[1:100],],
                num_clusters = 4,
                cores = 20,
                show_rownames = T)
print(p)
dev.off()
n_cluster = 3
pdf("pseudotime_diff_gene_heatmap.pdf")
p <- plot_pseudotime_heatmap(cd[sig_gene_names,],
                # num_clusters = n_cluster,
                cores = 20,
                return_heatmap = T,
                show_rownames = F)
print(p)
dev.off()

k = length(unique(cd$State))
diff_test_res_state <- differentialGeneTest(cd,fullModelFormulaStr ="~State", cores = 40)
sig_gene_names <- row.names(subset(diff_test_res_state, qval < 1e-5))
pdf("state_diff_gene_heatmap_top100.pdf",height=9)
p <- plot_pseudotime_heatmap(cd[sig_gene_names[1:100],],
                num_clusters = k,
                cores = 20,
                show_rownames = T)
print(p)
dev.off()

pdf("state_diff_gene_heatmap.pdf")
p <- plot_pseudotime_heatmap(cd[sig_gene_names,],
                num_clusters = k,
                cores = 20,
                return_heatmap = T,
                show_rownames = F)
print(p)
dev.off()


 
cd <- readRDS('CAF_monocle2.rds')
table(cd$State)
seurat_obj$mn2_state <- cd$State
# remove cluster 11
seurat_obj <- subset(seurat_obj, idents = '11', invert = T)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), verbose = F)
library(harmony)
seurat_obj <- seurat_obj %>% 
    RunHarmony("Patients", plot_convergence = FALSE,reduction = "pca")
seurat_obj <- RunUMAP(seurat_obj, reduction = 'harmony', dims = 1:30)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, verbose = FALSE, algorithm=4, resolution = 0.4)
seurat_obj <- FindClusters(seurat_obj, verbose = FALSE, algorithm=4, resolution = 0.8)
seurat_obj <- FindClusters(seurat_obj, verbose = FALSE, algorithm=4, resolution = 1)
seurat_obj <- FindClusters(seurat_obj, verbose = FALSE, algorithm=4, resolution = 1.2)
pdf(file = "umap_cluster.pdf", width = 7, height = 6)
print(DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.4", label = TRUE))
print(DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.8", label = TRUE))
print(DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1", label = TRUE))
print(DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.2", label = TRUE))
print(DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "cellsubtype", label = TRUE))
print(DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "sample", label = TRUE))
dev.off()

# DEGs
Idents(seurat_obj) <- seurat_obj$RNA_snn_res.1
seurat_obj_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat_obj_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
pdf(file = "seurat_obj_markers_heatmap.pdf", width = 10, height = 24)
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()
dev.off()
hei <- ceiling(length(unique(top10$gene))/4) * 3
pdf(file = "DEG_featureplots.pdf", width = 14.5, height = hei)
print(FeaturePlot(seurat_obj, features = unique(top10$gene), reduction = "umap", ncol = 4,raster=TRUE))
dev.off()
pdf(file = "DEG_vlnplots.pdf", width = 10, height = 10)
print(VlnPlot(seurat_obj, features = unique(top10$gene),
    ncol = 4, pt.size = 0))
dev.off()

# CAF markers
genes <- c("CFD", "C3","C7","APOD","PLA2G2A","CRABP2","LDHB", "COL11A1","COL10A1","POSTN","RGS5","COL1A1","COL1A2","WNT5A","F2RL2",'PAG1','NR2F2','NR2F1')
pdf(file = "strommarker_vlnplots.pdf", width = 10, height = 10)
print(VlnPlot(seurat_obj, features = genes,
    ncol = 4, pt.size = 0))
dev.off()
hei <- ceiling(length(genes)/4) * 3
pdf(file = "strommarker_featureplot.pdf", width = 14.5, height = hei)
print(FeaturePlot(seurat_obj, features = genes, reduction = "umap", ncol = 4,raster=TRUE))
dev.off()

# rename label
renamelabelfl <- "/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/PDAC/spatial/cellbin_v3/5.7_CAF_scRNA_monocle2/rename.txt"
renamelabel <- read.delim(renamelabelfl,header=F)
labels_cluster <- renamelabel[,2]
names(labels_cluster) <- levels(seurat_obj)
Idents(seurat_obj) <- seurat_obj$RNA_snn_res.1
table(Idents(seurat_obj))
seurat_obj <- RenameIdents(seurat_obj, labels_cluster)
seurat_obj$cellsubtype <- Idents(seurat_obj)
table(Idents(seurat_obj))
# DEGs_rename
seurat_obj_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat_obj_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
pdf(file = "DEG_heatmap.pdf", width = 10, height = 24)
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()
dev.off()
hei <- ceiling(length(unique(top10$gene))/4) * 3
pdf(file = "DEG_featureplots.pdf", width = 14.5, height = hei)
print(FeaturePlot(seurat_obj, features = unique(top10$gene), reduction = "umap", ncol = 4,raster=TRUE))
dev.off()
pdf(file = "DEG_vlnplots.pdf", width = 10, height = 30)
print(VlnPlot(seurat_obj, features = unique(top10$gene),
     pt.size = 0, flip = TRUE, stack = TRUE,log=FALSE)+ NoLegend())
dev.off()
pdf(file = "umap_cluster.pdf", width = 7, height = 6)
print(DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "cellsubtype", label = TRUE,repel=TRUE))
print(DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "sample", label = TRUE,repel=TRUE))
dev.off()
saveRDS(seurat_obj,'CAF_recluster.rds')

source("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/script/wangdi_code/SelfFunctions/colors.R")
pdf("cell_state_umap.pdf")
Idents(seurat_obj) <- seurat_obj$mn2_state
write.table(seurat_obj@meta.data,'metadata.txt',sep='\t',col.names=NA)
colors <- getDefaultColors(length(names(table(Idents(seurat_obj)))))
names(colors) <- unique(Idents(seurat_obj))
DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = "mn2_state" )
dev.off()
library(SeuratDisk)
i <- sapply(seurat_obj@meta.data, is.factor)
seurat_obj@meta.data[i] <- lapply(seurat_obj@meta.data[i], as.character)
seurat_obj2 <- DietSeurat(seurat_obj,
                            counts = TRUE,
                            data = TRUE,
                            scale.data = FALSE,
                            features = NULL,
                            assays = NULL,
                            dimreducs = c('pca','umap','nmf','harmony'),
                            graphs = NULL,
                            misc = TRUE)
SaveH5Seurat(seurat_obj2, filename = "CAF_recluster.h5Seurat",overwrite=TRUE)
Convert("CAF_recluster.h5Seurat", dest = "h5ad",overwrite=TRUE)


filtered_markers <- seurat_obj_markers[seurat_obj_markers$p_val_adj < 0.05, ]
filtered_markers %>%
    group_by(cluster) %>%
    top_n(n = 200, wt = avg_log2FC) -> top200
# get a df for each cluster in top200 as a colomn
dftop200 <- top200[,c('cluster','gene')] %>%
  pivot_wider(names_from = cluster, values_from = gene,names_repair='minimal')
dftop200 <- as.data.frame(dftop200)
dftop200 <- lapply(dftop200, function(col) col[[1]])
output_file <- "DEGs.txt"
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




# markers <- FindMarkers(seurat_obj, ident.1 = "7", verbose = FALSE, only.pos = TRUE, min.pct = 0.1, group.by='mn2_state',logfc.threshold=0.25)

markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,group.by='mn2_state',test.use='LR')
filtered_markers <- markers[markers$p_val_adj < 0.05, ]
filtered_markers %>%
    group_by(cluster) %>%
    top_n(n = 200, wt = avg_log2FC) -> top200


# get a df for each cluster in top200 as a colomn
dftop200 <- top200[,c('cluster','gene')] %>%
  pivot_wider(names_from = cluster, values_from = gene,names_repair='minimal')
dftop200 <- as.data.frame(dftop200)
dftop200 <- lapply(dftop200, function(col) col[[1]])
output_file <- "DEGs.txt"
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

library(ClusterGVis)
library(NbClust)
# 绘制pseudotime热图+富集分析合图
n_cluster = 4 # number of clusters
n_top_genes = 5 # number of top genes to show in each cluster
n_top_pathways = 5 # number of top pathways to show in each cluster
BEAM_res <- BEAM(cd, branch_point = 1, cores = 20)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
dat <- cd[row.names(subset(BEAM_res,qval < 1e-10)),]
y <- NbClust(cor_mat_module, distance = "euclidean", min.nc = 2, max.nc = 40 , method = "kmeans", index = "kl", alphaBeale = 0.1)

# get top n DEGs in each cell state              
seurat_obj$mn2_state <- cd$State
Idents(seurat_obj) <- seurat_obj$cellsubtype
seurat_obj_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,group.by='cellsubtype',test.use='wilcox')
filtered_markers <- seurat_obj_markers[seurat_obj_markers$p_val_adj < 0.05, ]
filtered_markers %>%
    group_by(cluster) %>%
    top_n(n = n_top_genes, wt = avg_log2FC) -> topn

df <- plot_genes_branched_heatmap2(cd[row.names(subset(BEAM_res,qval < 1e-20)),],
                                   branch_point = 1,
                                   num_clusters = 6,
                                   cores = 20,
                                   use_gene_short_name = T,
                                   show_rownames = T)
dat <- df$wide.res[,c(3:ncol(df$wide.res)-1)]                           
y <- NbClust(dat, distance = "euclidean", min.nc = 3, max.nc = 10 , method = "kmeans", index = "kl", alphaBeale = 0.1)
num_clusters <- y$Best.nc[[1]]
num_clusters <- 5
df <- plot_genes_branched_heatmap2(cd[row.names(subset(BEAM_res,qval < 1e-10)),],
                                   branch_point = 1,
                                   num_clusters = num_clusters,
                                   cores = 20,
                                   use_gene_short_name = T,
                                   show_rownames = T)
library(msigdbr)
library(org.Hs.eg.db)
library(clusterProfiler)
m_df <- msigdbr(species = "Homo sapiens", category = "H")%>% 
      dplyr::select(gs_name, gene_symbol)
m_df$gs_name <- gsub("HALLMARK_","",m_df$gs_name)
m_df_subset <- m_df[m_df$gene_symbol %in% rownames(seurat_obj),]
path_selected <- read.gmt('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/database/signature/signature_BP.gmt')
enrich_HM <- enrichCluster(object = df,
                        type = "ownSet",
                        id.trans = FALSE,
                        TERM2GENE = m_df_subset,
                        TERM2NAME = NA,
                        pvalueCutoff = 1,
                        topn = n_top_pathways,
                        seed = 42)
enrich_paths <- enrichCluster(object = df,
                        type = "ownSet",
                        id.trans = FALSE,
                        TERM2GENE = path_selected,
                        TERM2NAME = NA,
                        pvalueCutoff = 1,
                        topn = n_top_pathways,
                        seed = 42)
enrich_GO <- enrichCluster(object = df,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = n_top_pathways,
                        seed = 42)
# get a color vector according to enrich_HM$group
# colors = jjAnno::useMyCol("stallion2",n = length(unique(enrich_HM$group)))
col16 <- c('#a2cffe', '#87a922', '#ffa62b', '#f8481c',
           '#a6814c', '#a484ac', '#fc86aa', '#952e8f', '#02ccfe',
           '#2000b1', '#009337', '#ad0afd', '#3c9992', '#d8dcd6',
           '#cb6843')
colors = col16[1:length(unique(enrich_HM$group))]
go.col = c()
i = 1
for(n in table(enrich_HM$group)){
    go.col = c(go.col,rep(colors[i],n))
    i = i + 1
}
pdf("branch_heatmap_HM.pdf",width=14,height=7)                       
visCluster(object = df,
           plot.type = "both",
           pseudotime_col = c("powderblue","grey","lightcoral"),
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = c(unique(topn$gene),'USP53','MME','WNT5A'),
           markGenes.side = "left",
           annoTerm.data = enrich_HM,
           ctAnno.col = col16,
           go.col = go.col,
           cluster.order = c(2,5,1,6,3,4),
           add.bar = T,
           line.side = "left")
dev.off()

# selected pathways
colors = col24[1:length(unique(enrich_paths$group))]
go.col = c()
i = 1
for(n in table(factor(enrich_paths$group,levels = unique(enrich_paths$group)))){
    go.col = c(go.col,rep(colors[i],n))
    i = i + 1
}
pdf("branch_heatmap_pathways.pdf",width=14,height=7)                       
visCluster(object = df,
           plot.type = "both",
           pseudotime_col = c("powderblue","grey","lightcoral"),
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = unique(topn$gene),
           markGenes.side = "left",
           annoTerm.data = enrich_paths,
           ctAnno.col = col24,
           go.col = go.col,
           cluster.order = c(2,5,1,6,3,4),
           add.bar = T,
           line.side = "left")
dev.off()


colors = jjAnno::useMyCol("stallion2",n = length(unique(enrich_GO$group)))
go.col = c()
i = 1
for(n in table(enrich_GO$group)){
    go.col = c(go.col,rep(colors[i],n))
    i = i + 1
}
pdf("branch_heatmap_GO_BP.pdf",width=14,height=7)                       
visCluster(object = df,
           plot.type = "both",
           column_names_rot = 45,
           pseudotime_col = c("powderblue","grey","lightcoral"),
           show_row_dend = F,
           markGenes = unique(topn$gene),
           markGenes.side = "left",
           annoTerm.data = enrich_GO,
           go.col = go.col,
           add.bar = T,
           line.side = "left")
dev.off()
