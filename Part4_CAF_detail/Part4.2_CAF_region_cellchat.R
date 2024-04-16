

###
 # @Description: 
 # @Author: gongchanghao
 # @Date: 2023-08-23 15:11:58
 # @LastEditTime: 2023-09-15 17:48:15
 # @LastEditors: gongchanghao
 # @E-mail: gongchanghao@genomics.cn
 # @FilePath: /gongchanghao/script/spatial_pipeline/bin/cellbin_pipeline_v2/5.CAF/5.16_CAF_region_cellchat.R
 ###


library(CellChat)
library(Seurat)
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)

source("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/script/wangdi_code/SelfFunctions/colors.R")

regionstr <- args[1]
regions <- strsplit(regionstr,split = '-')[[1]]
# regions <- c('Tumor_around_region_1','Tumor_Cell_region_2')
wd <- sprintf("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/PDAC/spatial/cellbin_v3/5.16_CAF_region_cellchat/%s/",regionstr)
dir.create(wd)
setwd(wd)
seurat_obj_all <- readRDS('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/PDAC/spatial/cellbin_v3/2.1_merge_big_cell_CAF_new/all_merged_adata.rds')
seurat_obj <- seurat_obj_all[,seurat_obj_all$community_type %in% regions]
seurat_obj <- NormalizeData(seurat_obj, assay = "Spatial", verbose = FALSE)
seurat_obj$cluster <- seurat_obj$cellsubtype
seurat_obj$cluster[grep('^Fibro_|^Stellate_',seurat_obj$cellsubtype,invert=TRUE)] <- seurat_obj$celltype[grep('Fibro_|Stellate_',seurat_obj$cellsubtype,invert=TRUE)]

# 计算每个细胞类型的细胞数
celltype_counts <- table(seurat_obj$cluster)
# 获取细胞数大于等于20的细胞类型
celltype_to_keep <- names(celltype_counts[celltype_counts >= 10])
# 过滤细胞数大于等于20的细胞类型
seurat_obj <- subset(seurat_obj, subset = cluster %in% celltype_to_keep)
saveRDS(seurat_obj,'seurat_obj_for_communication.rds')
options(future.globals.maxSize= 891289600000)
cellchat <- createCellChat(object = seurat_obj, group.by = "cluster", assay = "Spatial")
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- readRDS('/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/database/cellchat/human/CellChatDB.human_nichenet.rds')
# CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
# interaction_input <- CellChatDB$interaction
# complex_input <- CellChatDB$complex
# cofactor_input <- CellChatDB$cofactor
# geneInfo <- CellChatDB$geneInfo
# write.csv(interaction_input, file = "interaction_input_CellChatDB.csv")
# write.csv(complex_input, file = "complex_input_CellChatDB.csv")
# write.csv(cofactor_input, file = "cofactor_input_CellChatDB.csv")
# write.csv(geneInfo, file = "geneInfo_input_CellChatDB.csv")
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
# future::plan("multicore", workers = 30) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
# cellchat <- computeCommunProb(cellchat,
#     type = "truncatedMean", trim = 0.1,
#     distance.use = TRUE, interaction.length = 200, scale.distance = 0.1, nboot = 20
# )
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.table(df.net, "CellCommunication.xls", sep = "\t", quote = FALSE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, "cellchat.rds")
# visualize the aggregated cell-cell communication network.
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd = TRUE) # nolint
pdf("circlePlot.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")
dev.off()

# Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. Here we also control the parameter edge.weight.max so that we can compare edge weights between differet network.
mat <- cellchat@net$weight
par(mfrow = c(3, 4), xpd = TRUE)
pdf("circlePlot_subset.pdf")
for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    try(netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i]))
}
dev.off()

pathways <- cellchat@netP$pathways
# dir.create("pathway_netVisual_aggregate", recursive = TRUE)
# for (pathways.show in pathways) {
#     # Hierarchy plot
#     # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells
#     vertex.receiver <- seq(1, 4) # a numeric vector.
#     pdf(paste0("pathway_netVisual_aggregate/", pathways.show, ".pdf"))
#     # print(netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.4, vertex.label.cex = 2, point.size = 1))
#     netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver)
#     # Circle plot
#     netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
#     netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#     print(netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds"))
#     dev.off()
# }
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pathways <- c("COLLAGEN",
                "FN1",
                "LAMININ",
                "SPP1",
                "THBS",
                "MK",
                "PERIOSTIN",
                "JAM",
                "PARs",
                "MPZ",
                "PDGF",
                "VEGF",
                "TENASCIN",
                "COMPLEMENT",
                "NEGR",
                "ANGPTL",
                "SEMA4",
                "NOTCH",
                "EPHA",
                "BMP",
                "ITGB2",
                "ADGRE5",
                "SEMA3",
                "TGFb",
                "ncWNT",
                "EPHB",
                "AGRN",
                "CXCL",
                "CD45",
                "EGF",
                "ICAM",
                "CEACAM",
                "CCL",
                "PTPRM",
                "SEMA5")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat,height=9,width=4, pattern = "outgoing",signaling = pathways)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat,height=9,width=4,pattern = "incoming",signaling = pathways)
pdf("signalingRole_heatmap2.pdf",width=6,height=7)
print(ht1 + ht2)
dev.off()
library(NMF)
library(ggalluvial)
pdf("selectK_outgoing.pdf")
print(selectK(cellchat, pattern = "outgoing"))
dev.off()
pdf("selectK_incoming.pdf")
print(selectK(cellchat, pattern = "incoming"))
dev.off()
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
pdf("netVisual_embedding_functional.pdf",height=7,width=7)
print(netVisual_embedding(cellchat, type = "functional", label.size = 3.5))
print(netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2))
dev.off()
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
pdf("netVisual_embedding_structural.pdf",height=7,width=7)
print(netVisual_embedding(cellchat, type = "structural", label.size = 3.5))
print(netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2))
dev.off()
saveRDS(cellchat, "cellchat.rds")

cellchat <- readRDS('cellchat.rds')

# do CellChat
# coords <- seurat_obj@meta.data[, c("x", "y")]
# cellchat <- createCellChat(object = seurat_obj, group.by = "community_type", assay = "Spatial", coordinates = coords, scale.factors = list(spot.diameter = 2.2 spot = 1), datatype = "spatial")
# cellchat <- createCellChat(object = seurat_obj, group.by = "community_type", assay = "Spatial")
# groupSize <- as.numeric(table(cellchat@idents))
# CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
# interaction_input <- CellChatDB$interaction
# complex_input <- CellChatDB$complex
# cofactor_input <- CellChatDB$cofactor
# geneInfo <- CellChatDB$geneInfo
# write.csv(interaction_input, file = "interaction_input_CellChatDB.csv")
# write.csv(complex_input, file = "complex_input_CellChatDB.csv")
# write.csv(cofactor_input, file = "cofactor_input_CellChatDB.csv")
# write.csv(geneInfo, file = "geneInfo_input_CellChatDB.csv")


# Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. Here we also control the parameter edge.weight.max so that we can compare edge weights between differet network.
# mat <- cellchat@net$weight
# par(mfrow = c(3, 4), xpd = TRUE)
# pdf(paste0(od, "/circlePlot_subset.pdf"))
# for (i in 1:nrow(mat)) {
#     mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#     mat2[i, ] <- mat[i, ]
#     try(netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i]))
# }
# dev.off()

# pathways <- cellchat@netP$pathways
# dir.create(paste0(od, "/pathway_netVisual_aggregate"), recursive = TRUE)
# for (pathways.show in pathways) {
#     # Hierarchy plot
#     # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells
#     vertex.receiver <- seq(1, 4) # a numeric vector.
#     pdf(paste0(od, "/pathway_netVisual_aggregate/", pathways.show, ".pdf"))
#     # print(netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.4, vertex.label.cex = 2, point.size = 1))
#     netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver)
#     # Circle plot
#     netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
#     netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#     print(netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds"))
#     dev.off()
#     pdf("Neoplastic_Cell_Mesenchymal_Cell_gene.pdf", height = 12, width = 12)
#     netVisual_chord_gene(cellchat, sources.use = c("Neoplastic_Cell"), targets.use = c("Mesenchymal_Cell"), lab.cex = 0.5, legend.pos.y = 30)
#     dev.off()
#     pdf("Mesenchymal_Cell_Neoplastic_Cell_gene.pdf", height = 12, width = 12)
#     netVisual_chord_gene(cellchat, sources.use = c("Mesenchymal_Cell"), targets.use = c("Neoplastic_Cell"), lab.cex = 0.5, legend.pos.y = 30)
#     dev.off()
#     pdf("Immune_Cell_Mesenchymal_Cell_gene.pdf", height = 12, width = 12)
#     netVisual_chord_gene(cellchat, sources.use = c("Immune_Cell"), targets.use = c("Mesenchymal_Cell"), lab.cex = 0.5, legend.pos.y = 30)
#     dev.off()
#     pdf("Immune_Cell_Neoplastic_Cell_gene.pdf", height = 12, width = 12)
#     netVisual_chord_gene(cellchat, sources.use = c("Immune_Cell"), targets.use = c("Neoplastic_Cell"), lab.cex = 0.5, legend.pos.y = 30)
#     dev.off()

# }
celltypes <- unique(cellchat@idents)
CAF_celltypes <- celltypes[grep('^Fibro_|Stellate_',celltypes)]
Tumor_celltypes <- celltypes[grep('^Tumor_',celltypes)]
# Macro_celltypes <- celltypes[grep('^Macrophage_',celltypes)]
Immune_celltypes <- celltypes[celltypes %in% c('T_Cell','B_Cell','T/NK_Cell','Macrophage_Cell','Mast_Cell','Plasma_Cell')]
source_celltypes <- 'Fibro_USP53'
source_celltypes_str <- paste(source_celltypes,collapse = '_')
source_celltypes_str <- gsub('/','_',source_celltypes_str)
target_celltypes <- c('Stellate_MYH11','Macrophage_Cell','Fibro_USP53','Duct_like_Cell','T/NK_Cell')
target_celltypes_str <- paste(target_celltypes,collapse = '_')
target_celltypes_str <- gsub('/','_',target_celltypes_str)
pdf(sprintf("%s_%s_gene_chord.pdf",source_celltypes_str,target_celltypes_str), height = 20, width = 20)
netVisual_chord_gene(cellchat, sources.use = source_celltypes, targets.use = target_celltypes, lab.cex = 0.8, legend.pos.y = 30,thres=1e-5)
dev.off()
pdf(sprintf("%s_%s_gene_bubble.pdf",source_celltypes_str,target_celltypes_str), height = 40, width = 5)
netVisual_bubble(cellchat, sources.use = source_celltypes, targets.use = target_celltypes, remove.isolate = FALSE)
dev.off()
pdf("Fibro_F2RL2_other_Cell_gene_chord.pdf", height = 40, width = 40)
netVisual_chord_gene(cellchat, sources.use = CAF_celltypes, targets.use = Tumor_celltypes, lab.cex = 0.8, legend.pos.y = 30,thres=0.01)
dev.off()
pdf("Tumor_CAF_Cell_gene_chord.pdf", height = 12, width = 12)
netVisual_chord_gene(cellchat, sources.use = Tumor_celltypes, targets.use = CAF_celltypes, lab.cex = 0.8, legend.pos.y = 30)
dev.off()
pdf("CAF_Immune_Cell_gene_chord.pdf", height = 12, width = 12)
netVisual_chord_gene(cellchat, sources.use = CAF_celltypes, targets.use = Immune_celltypes, lab.cex = 0.8, legend.pos.y = 30)
dev.off()
pdf("Immune_CAF_Cell_gene_chord.pdf", height = 12, width = 12)
netVisual_chord_gene(cellchat, sources.use = Immune_celltypes, targets.use = CAF_celltypes, lab.cex = 0.8, legend.pos.y = 30)
dev.off()

target_CAFs <- c('Fibro_F2RL2','Fibro_MME')
df.net <- subsetCommunication(cellchat,sources.use = target_CAFs, targets.use = Tumor_celltypes)
pairLR.use <- df.net[,'interaction_name',drop=F]
pdf('CAF_Tumor_Cell_bubble.pdf',height=40,width=7)
netVisual_bubble(cellchat, sources.use = CAF_celltypes, targets.use = Tumor_celltypes, remove.isolate = FALSE,pairLR.use=pairLR.use)
dev.off()

df.net <- subsetCommunication(cellchat,targets.use = target_CAFs, sources.use = Tumor_celltypes)
pairLR.use <- df.net[,'interaction_name',drop=F]

pdf('Tumor_CAF_Cell_bubble.pdf',height=40,width=7)
netVisual_bubble(cellchat, sources.use = Tumor_celltypes, targets.use = CAF_celltypes,remove.isolate = FALSE,pairLR.use=pairLR.use)
dev.off()


pdf('CAF_Immune_Cell_bubble.pdf',height=40,width=7)
netVisual_bubble(cellchat, sources.use = CAF_celltypes, targets.use = Immune_celltypes, remove.isolate = FALSE)
dev.off()
pdf('Immune_CAF_Cell_bubble.pdf',height=40,width=7)
netVisual_bubble(cellchat, sources.use = Immune_celltypes, targets.use = CAF_celltypes, remove.isolate = FALSE)
dev.off()



# pdf("Macro_ADM_Cell_gene_chord.pdf", height = 7, width = 7)
# netVisual_chord_gene(cellchat, sources.use = Macro_celltypes, targets.use = Tumor_celltypes, lab.cex = 0.8, legend.pos.y = 30)
# dev.off()
# pdf("Immune_Cell_Mesenchymal_Cell_gene.pdf", height = 7, width = 7)
# netVisual_chord_gene(cellchat, sources.use = c("Immune_Cell"), targets.use = c("Mesenchymal_Cell"), lab.cex = 0.5, legend.pos.y = 30)
# dev.off()
# pdf("Immune_Cell_Neoplastic_Cell_gene.pdf", height = 7, width = 7)
# netVisual_chord_gene(cellchat, sources.use = c("Immune_Cell"), targets.use = c("Neoplastic_Cell"), lab.cex = 0.5, legend.pos.y = 30)
# dev.off()

# pdf('Neoplastic_Cell_Mesenchymal_Cell_gene_bubble.pdf',height=18,width=7)
# netVisual_bubble(cellchat, sources.use = c(7,9,11), targets.use = c(1,2,3,5), remove.isolate = FALSE)
# dev.off()
# pdf('Mesenchymal_Cell_Neoplastic_Cell_gene_bubble.pdf',height=25,width=7)
# netVisual_bubble(cellchat, sources.use = c(1,2,3,5), targets.use = c(7,9,11), remove.isolate = FALSE)
# dev.off()




# pathways.show <- "COLLAGEN"
# pdf("COLLAGEN_contribution.pdf")
# print(netAnalysis_contribution(cellchat, signaling = pathways.show))
# dev.off()
# pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
# LR.show <- pairLR.CXCL[43,] # show one ligand-receptor pair
# vertex.receiver = seq(1,4) # a numeric vector
# pdf("COLLAGEN_SDC4.pdf",width=12,height=12)
# print(netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord"))
# dev.off()


# pdf("COLLAGEN_GeneExpression.pdf")
# print(plotGeneExpression(cellchat, signaling = "COLLAGEN"))
# dev.off()
# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# pdf("netAnalysis_signalingRole_network.pdf")
# print(netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10))
# dev.off()

library(NMF)
library(ggalluvial)
# pdf("selectK_outgoing.pdf")
# print(selectK(cellchat, pattern = "outgoing"))
# dev.off()
nPatterns = 7
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
pdf("netAnalysis_outgoing.pdf",width=12,height=7)
print(netAnalysis_river(cellchat, pattern = "outgoing"))
print(netAnalysis_dot(cellchat, pattern = "outgoing"))
dev.off()
# pdf("selectK_incoming.pdf")
# print(selectK(cellchat, pattern = "incoming"))
# dev.off()
nPatterns = 8
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
pdf("netAnalysis_incoming.pdf",width=12,height=7)
print(netAnalysis_river(cellchat, pattern = "incoming"))
print(netAnalysis_dot(cellchat, pattern = "incoming"))
dev.off()

pdf('LAMININ_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "LAMININ")
dev.off()
pdf('WNT_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "WNT")
dev.off()
pdf('ncWNT_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "ncWNT")
dev.off()
pdf('EGF_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "EGF")
dev.off()
pdf('ITGB2_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "ITGB2")
dev.off()
pdf('COLLAGEN_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "COLLAGEN")
dev.off()
pdf('SPP1_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "SPP1")
dev.off()
pdf('COMPLEMENT_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "COMPLEMENT")
dev.off()
pdf('APP_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "APP")
dev.off()
pdf('GALECTIN_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "GALECTIN")
dev.off()
pdf('GRN_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "GRN")
dev.off()
pdf('AGRN_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "AGRN")
dev.off()
pdf('NEGR_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "NEGR")
dev.off()
pdf('CSF_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "CSF")
dev.off()
pdf('THBS_gene_expr_diff.pdf')
plotGeneExpression(cellchat, signaling = "THBS")
dev.off()
pdf('LGALS3_ANXA2_gene_expr_diff.pdf')
plotGeneExpression(cellchat, features = c("LGALS3",'ANXA2'))
dev.off()
pdf('TIMP1_TIMP3_CD63_MMP2_CD44_gene_expr_diff.pdf')
plotGeneExpression(cellchat, features = c("TIMP1",'TIMP3','CD63','MMP2','CD44'))
dev.off()




