
library(CellChat)
library(Seurat)
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)

source("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/script/wangdi_code/SelfFunctions/colors.R")
wd <- sprintf("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/PDAC/spatial/cellbin_v3/3.11_region_cellchat_compare/all")
dir.create(wd,recursive = TRUE)
setwd(wd)
regionlist <- c('Tumor_around_region_1','Tumor_around_region_2','Tumor_around_region_3','Tumor_around_region_4','Tumor_around_region_5','Tumor_Cell_region_1','Tumor_Cell_region_2','Tumor_Cell_region_3','Tumor_Cell_region_4')
regionshort <- c('TAR1','TAR2','TAR3','TAR4','TAR5','TCR1','TCR2','TCR3','TCR4')
object.list <- list()
i <- 1
for(region in regionlist){
    cellchat <- readRDS(paste0('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/PDAC/spatial/cellbin_v3/3.10_region_cellchat/',region,'/cellchat.rds'))
    cellchat <- netAnalysis_computeCentrality(cellchat)
    object.list[[regionshort[i]]]<- cellchat
    i <- i+1
}

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat1 <- readRDS("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/PDAC/spatial/cellbin_v3/3.10_region_cellchat/Tumor_Cell_region_1/cellchat.rds")
cellchat2 <- readRDS("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/PDAC/spatial/cellbin_v3/3.10_region_cellchat/Tumor_Cell_region_2/cellchat.rds")
interaction_type <- intersect(unique(cellchat1@idents),unique(cellchat2@idents))
source('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/PDAC/spatial/cellbin_v3/3.11_region_cellchat_compare/CellChat_class.R')
cellchat1 <- subsetCellChat(cellchat1,idents.use = interaction_type)
cellchat2 <- subsetCellChat(cellchat2,idents.use = interaction_type)
object.list <- list(TCR1 = cellchat1, TCR2 = cellchat2)
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
}
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# 比较交互总数和交互强度
df <- as.data.frame(sapply(cellchat@net, function(x) sum(x$count)))
df$group = rownames(df)
colnames(df)[1] <- 'count'
gg <- ggplot(df, aes(x=group, y=count)) +
        geom_bar(stat="identity", width=0.6) + scale_color_npg()
colors <- col16
library(ggsci)
gg1 <- compareInteractions(cellchat, show.legend = F,angle.x=45,x.lab.rot = TRUE)
gg2 <- compareInteractions(cellchat, show.legend = F,angle.x=45,x.lab.rot = TRUE, measure = "weight")
pdf('compareInteractions.pdf',width=6,height=4)
print(gg1+gg2)
dev.off()

# 比较不同细胞群之间的相互作用数量和相互作用强度
par(mfrow = c(1,2), xpd=TRUE)
pdf("netVisual_diffInteraction.pdf")
print(netVisual_diffInteraction(cellchat, weight.scale = T))
print(netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight"))
dev.off()

# 热图
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
pdf("netVisual_heatmap.pdf",height=4,width=6)
print(gg1 + gg2)
dev.off()

# 计算每个细胞组的最大细胞数以及所有数据集的最大交互数（或交互权重）
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("netVisual_circle_max_weight.pdf")
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

# 不同细胞类型之间相互作用的数量或相互作用强度不同
# group.cellType1 <- c(rep("Mesenchymal", 4), rep("Neoplastic", 6))
# group.cellType2 <- c(rep("Mesenchymal", 6), rep("Neoplastic", 4))
# group.cellType1 <- factor(group.cellType1, levels = c("Mesenchymal", "Neoplastic"))
# group.cellType2 <- factor(group.cellType2, levels = c("Mesenchymal", "Neoplastic"))
# object.list[[1]] <- mergeInteractions(object.list[[1]], group.cellType1)
# object.list[[2]] <- mergeInteractions(object.list[[1]], group.cellType2)
# # object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
# cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.

# 比较二维空间中的主要源和目标
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf("source_targets_signalingRole_scatter.pdf",width=12)
patchwork::wrap_plots(plots = gg)
dev.off()

gg0 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Fibroblast_Cell")
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Tumor_Cell")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Macrophage_Cell")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Duct_like_Cell")
pdf("source_signalingChanges_scatter.pdf",width=15,height=9)
patchwork::wrap_plots(ncol=2,nrow=2,plots = list(gg0,gg1,gg2,gg3))
dev.off()

# 根据功能相似性识别信号基团
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
pdf("embeddingPairwise_functional.pdf")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
dev.off()

# 根据结构相似性识别信号基团
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
pdf("embeddingPairwise_structural.pdf")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
dev.off()

# 比较与每个细胞群相关的传出（或传入）信号
library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
pathway.union <-c("COLLAGEN",
                    "FN1",
                    "SPP1",
                    "LAMININ",
                    "THBS",
                    "MK",
                    "PERIOSTIN",
                    "ANGPTL",
                    "HSPG",
                    "PARs",
                    "VEGF",
                    "AGRN",
                    "GALECTIN",
                    "CDH",
                    "EPHA",
                    "NEGR",
                    "APP",
                    "GRN",
                    "COMPLEMENT",
                    "SEMA4",
                    "NOTCH",
                    "ADGRE5",
                    "ITGB2",
                    "TGFb",
                    "CEACAM",
                    "CD45",
                    "SEMA3",
                    "CDH5",
                    "EGF",
                    "ESAM",
                    "CXCL",
                    "PTPRM",
                    "CCL")
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 3, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 3, height = 6)
pdf("outgoing_signalingRole_heatmap.pdf",width=6,height=6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 4, height = 12, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 4, height = 12, color.heatmap = "GnBu")
pdf("incoming_signalingRole_heatmap.pdf",width=6,height=8)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
library(ggplot2)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 4, height = 8, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 4, height = 8, color.heatmap = "OrRd")
pdf("overall_signalingRole_heatmap.pdf",width=8,height=8)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

# 通过比较通信概率来识别功能失调的信号

gg1 <- netVisual_bubble(cellchat, sources.use = c('Tumor_Cell'), targets.use = c('Fibroblast_Cell'),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in TCR2", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c('Tumor_Cell'), targets.use = c('Fibroblast_Cell'),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in TCR2", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("Increased_Decrease_signaling_Tumor_CAF.pdf",width=7,height=10)
print(gg1 + gg2)
dev.off()

gg1 <- netVisual_bubble(cellchat, sources.use = c('Fibroblast_Cell'), targets.use = c('Tumor_Cell'),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in TCR2", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c('Fibroblast_Cell'), targets.use = c('Tumor_Cell'),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in TCR2", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("Increased_Decrease_signaling_CAF_Tumor.pdf",width=7,height=18)
print(gg1 + gg2)
dev.off()



# 使用差异表达分析识别功能失调的信号传导
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "TCR2"
neg.dataset = "TCR1"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.263, thresh.p =1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = pos.dataset,ligand.logFC = 0.263, receptor.logFC = 0.1)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = neg.dataset,ligand.logFC = -0.263, receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c('Fibroblast_Cell'), targets.use = c('Tumor_Cell'), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", pos.dataset))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('Fibroblast_Cell'), targets.use = c('Tumor_Cell'), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", pos.dataset))
#> Comparing communications on a merged object
pdf("up_down_signaling_CAF_Tumor.pdf",width=7,height=12)
print(gg1 + gg2)
dev.off()


par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c('Fibroblast_Cell'), targets.use = c('Tumor_Cell'),slot.name = "netP", title.name = paste0("Signaling pathways sending from fibroblast - ", names(object.list)[i]), legend.pos.x = 10)
}

gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c('Tumor_Cell'), targets.use = c('Fibroblast_Cell'), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", pos.dataset))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('Tumor_Cell'), targets.use = c('Fibroblast_Cell'), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", pos.dataset))
#> Comparing communications on a merged object
pdf("up_down_signaling_Tumor_CAF.pdf",width=7,height=9)
print(gg1 + gg2)
dev.off()

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
pdf("up_down_signaling_CAF_Tumor_chord.pdf",width=12,height=12)
netVisual_chord_gene(object.list[[2]], sources.use = c('Fibroblast_Cell'), targets.use = c('Tumor_Cell'), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", pos.dataset))
netVisual_chord_gene(object.list[[1]], sources.use = c('Fibroblast_Cell'), targets.use = c('Tumor_Cell'), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", pos.dataset))
dev.off()

par(mfrow = c(1,2), xpd=TRUE)
pdf("up_down_signaling_Tumor_CAF_chord.pdf",width=12,height=12)
netVisual_chord_gene(object.list[[2]], sources.use = c('Tumor_Cell'), targets.use = c('Fibroblast_Cell'), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", pos.dataset))
netVisual_chord_gene(object.list[[1]], sources.use = c('Tumor_Cell'), targets.use = c('Fibroblast_Cell'), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", pos.dataset))
dev.off()

# 比较不同数据集之间的信号基因表达分布
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("TCR1", "TCR2")) # set factor level
pdf('ncWNT_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "ncWNT", split.by = "datasets", colors.ggplot = T)
dev.off()
pdf('COLLAGEN_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "COLLAGEN", split.by = "datasets", colors.ggplot = T)
dev.off()
pdf('SPP1_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "SPP1", split.by = "datasets", colors.ggplot = T)
dev.off()
pdf('LAMININ_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "LAMININ", split.by = "datasets", colors.ggplot = T)
dev.off()
pdf('FN1_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "FN1", split.by = "datasets", colors.ggplot = T)
dev.off()
pdf('COMPLEMENT_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "COMPLEMENT", split.by = "datasets", colors.ggplot = T)
dev.off()
pdf('CXCL_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)
dev.off()
