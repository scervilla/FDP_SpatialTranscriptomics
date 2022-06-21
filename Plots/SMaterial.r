################################################################################

#                             Supplementary figures

################################################################################
library(Seurat)
library(ggplot2)
library(ggthemes)
library(STutility)
library(pheatmap)
library(gprofiler2)
# SFig 1
source("Desktop/IJC/datasets/IGTP/figuresPaper/scripts/utilities/plotfunct.r")
## prefiltering

### OV4A

infoTable <- data.frame(samples="Desktop/IJC/datasets/IGTP/4A/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5")
infoTable$spotfiles="Desktop/IJC/datasets/IGTP/4A/spatial/tissue_positions_list.csv"
infoTable$imgs="Desktop/IJC/datasets/IGTP/4A/spatial/tissue_hires_image.png"
infoTable$json="Desktop/IJC/datasets/IGTP/4A/spatial/scalefactors_json.json"

OV4A <- InputFromTable(infotable = infoTable, 
                       platform =  "Visium",
                       minUMICountsPerGene = 5)
annotLookup <- read.table("Desktop/IJC/datasets/annotLook.txt", header = T)
protein_genes <- annotLookup$gene_name[annotLookup$gene_type %in% c("protein_coding", "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene", "IG_LV_gene", "IG_V_gene", "IG_J_gene", "IG_C_gene" , "IG_D_gene")]
selected_genes <- rownames(OV4A)
selected_genes <- protein_genes[protein_genes %in% selected_genes]
selected_genes <- selected_genes[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", selected_genes)]
selected_genes <- selected_genes[!grepl("^MT-", selected_genes)]

OV4A <- SubsetSTData(OV4A, features = selected_genes)

### OV5B

infoTable <- data.frame(samples="Desktop/IJC/datasets/IGTP/5A/filtered_feature_bc_matrix.h5")
infoTable$spotfiles="Desktop/IJC/datasets/IGTP/5A/spatial//tissue_positions_list.csv"
infoTable$imgs="Desktop/IJC/datasets/IGTP/5A/spatial/tissue_hires_image.png"
infoTable$json="Desktop/IJC/datasets/IGTP/5A/spatial/scalefactors_json.json"

OV5A <- InputFromTable(infotable = infoTable, 
                       platform =  "Visium",
                       minUMICountsPerGene = 5)


protein_genes <- annotLookup$gene_name[annotLookup$gene_type %in% c("protein_coding", "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene", "IG_LV_gene", "IG_V_gene", "IG_J_gene", "IG_C_gene" , "IG_D_gene")]
selected_genes <- rownames(OV5A)
selected_genes <- protein_genes[protein_genes %in% selected_genes]
selected_genes <- selected_genes[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", selected_genes)]
selected_genes <- selected_genes[!grepl("^MT-", selected_genes)]

OV5A <- SubsetSTData(OV5A, features = selected_genes)

### OV5B

infoTable <- data.frame(samples="Desktop/IJC/datasets/IGTP/5B/filtered_feature_bc_matrix.h5")
infoTable$spotfiles="Desktop/IJC/datasets/IGTP/5B/images/tissue_positions_list.csv"
infoTable$imgs="Desktop/IJC/datasets/IGTP/5B/images/tissue_hires_image.png"
infoTable$json="Desktop/IJC/datasets/IGTP/5B/images/scalefactors_json.json"

OV5B <- InputFromTable(infotable = infoTable, 
                       platform =  "Visium",
                       minUMICountsPerGene = 5)


protein_genes <- annotLookup$gene_name[annotLookup$gene_type %in% c("protein_coding", "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene", "IG_LV_gene", "IG_V_gene", "IG_J_gene", "IG_C_gene" , "IG_D_gene")]
selected_genes <- rownames(OV5B)
selected_genes <- protein_genes[protein_genes %in% selected_genes]
selected_genes <- selected_genes[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", selected_genes)]
selected_genes <- selected_genes[!grepl("^MT-", selected_genes)]

OV5B <- SubsetSTData(OV5B, features = selected_genes)

OV4A$ident <- "P1"
OV5A$ident <- "P2 (A)"
OV5B$ident <- "P2 (B)"

OV <- merge(OV4A, OV5A)
OV <- merge(OV, OV5B)
Idents(OV) <- OV$ident
p <- VlnPlot(OV, features = c("nCount_RNA", "nFeature_RNA"), 
             group.by = "ident", cols = colorblind_pal()(4)[2:4]) & labs(x="") 
pdf("Desktop/IJC/TFG/SuppFigures/SFig1_prefiltered.pdf", width = 9, height = 5)
plot(p)
dev.off()

rm(OV4A, OV5A, OV5B, OV)
gc()

### Post filtering

OV4A <- readRDS("Desktop/IJC/TFG/RDS/OV4A.rds")
OV5A <- readRDS("Desktop/IJC/TFG/RDS/OV5A.rds")
OV5B <- readRDS("Desktop/IJC/TFG/RDS/OV5B.rds")

OV4A$ident <- "P1"
OV5A$ident <- "P2 (A)"
OV5B$ident <- "P2 (B)"

OV <- merge(OV4A, OV5A)
OV <- merge(OV, OV5B)
Idents(OV) <- OV$ident
p <- VlnPlot(OV, features = c("nCount_RNA", "nFeature_RNA"), 
             group.by = "ident", cols = colorblind_pal()(4)[2:4]) & labs(x="") 
pdf("Desktop/IJC/TFG/SuppFigures/SFig1_postfiltered.pdf", width = 9, height = 5)
plot(p)
dev.off()

rm(OV4A, OV5A, OV5B, OV)
gc()


# SFig 2
samples.combined <- readRDS("Desktop/IJC/TFG/")


#SFig 3
OV4A <- readRDS("Desktop/IJC/TFG/RDS/OV4A.rds")
OV5A <- readRDS("Desktop/IJC/TFG/RDS/OV5A.rds")
OV5B <- readRDS("Desktop/IJC/TFG/RDS/OV5B.rds")

unique(rownames(OV4A[["Deconvolution_Major"]]))
celltypes <- unique(rownames(OV4A[["Deconvolution_Major"]]))
names_major <- c("CAFs", "Endothelial cells", "Epithelial", "TAMs", "M1 Macrophages", "T cells",
                 "Proliferative Epithelial", "Myofibroblast", "Dendritic cells", "Plasma B cells",
                 "Natural Killer", "Plasmacytoid Dendritic cells", "B cells", "Mast cells")
## OV4A
DefaultAssay(OV4A) <- "Deconvolution_Major"
for (n in 1:length(celltypes)) {
  p <- SpatialFeaturePlot(OV4A, features = celltypes[n]) + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
    theme(legend.position="bottom", plot.title = element_text(hjust = 0.5, size = 25),
          legend.key.height = unit(0.4, 'cm'), legend.key.width  = unit(2.75, 'cm'),
          legend.direction="horizontal", legend.margin=margin(t=-15, ),
          legend.text = element_text(size=15)) + ggtitle(names_major[n])
  pdf(paste0("Desktop/IJC/TFG/SuppFigures/SFig3_OV4A_Major_",celltypes[n], ".pdf"), width = 10, height = 7)
  plot(p)
  dev.off()  
}

## OV5A
DefaultAssay(OV5A) <- "Deconvolution_Major"
for (n in 1:length(celltypes)) {
  p <- SpatialFeaturePlot(OV5A, features = celltypes[n]) + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
    theme(legend.position="bottom", plot.title = element_text(hjust = 0.5, size = 25),
          legend.key.height = unit(0.4, 'cm'), legend.key.width  = unit(2.75, 'cm'),
          legend.direction="horizontal", legend.margin=margin(t=-15),
          legend.text = element_text(size=15)) + ggtitle(names_major[n])
  pdf(paste0("Desktop/IJC/TFG/SuppFigures/SFig3_OV5A_Major_",celltypes[n], ".pdf"), width = 10, height = 7)
  plot(p)
  dev.off()  
}

## OV5B
DefaultAssay(OV5B) <- "Deconvolution_Major"
for (n in 1:length(celltypes)) {
  p <- SpatialFeaturePlot(OV5B, features = celltypes[n]) + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
    theme(legend.position="bottom", plot.title = element_text(hjust = 0.5, size = 25),
          legend.key.height = unit(0.4, 'cm'), legend.key.width  = unit(2.75, 'cm'),
          legend.direction="horizontal", legend.margin=margin(t=-15),
          legend.text = element_text(size=15)) + ggtitle(names_major[n])
  pdf(paste0("Desktop/IJC/TFG/SuppFigures/SFig3_OV5B_Major_",celltypes[n], ".pdf"), width = 10, height = 7)
  plot(p)
  dev.off()  
}


# SFig6
OV4A <- readRDS("Desktop/IJC/TFG/RDS/OV4A.rds")
OV5A <- readRDS("Desktop/IJC/TFG/RDS/OV5A.rds")
OV5B <- readRDS("Desktop/IJC/TFG/RDS/OV5B.rds")
DefaultAssay(OV4A) <- "Deconvolution_Major"
DefaultAssay(OV5A) <- "Deconvolution_Major"
DefaultAssay(OV5B) <- "Deconvolution_Major"
d <- data.frame(matrix(ncol = nrow(OV4A), nrow = 0))
colnames(d) <- rownames(OV4A)

tmp <- as.data.frame(t(as.data.frame(OV4A@assays$Deconvolution_Major@data)))
tmp$sample <- "OV4A"
d <- rbind(d, tmp)

tmp <- as.data.frame(t(as.data.frame(OV5A@assays$Deconvolution_Major@data)))
tmp$sample <- "OV5A"
d <- rbind(d, tmp)

tmp <- as.data.frame(t(as.data.frame(OV5B@assays$Deconvolution_Major@data)))
tmp$sample <- "OV5B"
d <- rbind(d, tmp)

colnames(d) <- make.names(colnames(d))
celltypes2 <- make.names(celltypes)

for (n in 1:length(celltypes2)) {
  p <- ggplot(d) + geom_density(aes_string(x=celltypes2[n], col="sample"), alpha=0.3) + theme_classic() +
    labs(x="Cell type proportion", y="Density", title=names_major[n]) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20,),axis.title = element_text(size = 20),
          legend.position = "none") + scale_color_manual(values = colorblind_pal()(3))
  pdf(paste0("Desktop/IJC/TFG/SuppFigures/SFig_6_",names_major[n],".pdf"), width = 4, height = 6)
  plot(p)
  dev.off()
}

#SFig7
OV4A <- readRDS("Desktop/IJC/TFG/RDS/OV4A.rds")
OV5 <- readRDS("Desktop/IJC/TFG/RDS/5merged.rds")

DefaultAssay(OV4A) <- "Deconvolution_Major"
cor4A <- sapply(rownames(OV4A), function(celltype){
  sapply(1:15, function(factor){
    cor(OV4A@assays$Deconvolution_Major@data[celltype,], OV4A@reductions$NMF@cell.embeddings[,factor])
  })
})
rownames(cor4A) <- paste0("factor_", 1:15)
pheatmap(t(cor4A), angle_col = 315, color =   rev(paletteer::paletteer_c("grDevices::RdBu", 60) ),fontsize = 12)



DefaultAssay(OV5) <- "Deconvolution_Major"
cor4A <- sapply(rownames(OV5), function(celltype){
  sapply(1:15, function(factor){
    cor(OV5@assays$Deconvolution_Major@data[celltype,], OV5@reductions$NMF@cell.embeddings[,factor])
  })
})
rownames(cor4A) <- paste0("factor_", 1:15)
pheatmap(t(cor4A), angle_col = 315, color =   rev(paletteer::paletteer_c("grDevices::RdBu", 60) ),fontsize = 12)


#Tables

## Functional enrichement



select_top <- function(st, factor) {
  vals <- st@reductions$NMF@feature.loadings[, factor]
  logvals <- log(vals)
  logvals <- logvals[is.finite(logvals)]
  thr <- mean(logvals) + 1.645*sd(logvals)
  names(logvals[logvals > thr])
}


top.genes4 <- list()
for (factor in 1:15) {
  top.genes4[[factor]] <- select_top(OV4A, factor)
}


#compute gene ontology (GO:BP) for each factor
pathways4 <- lapply(seq_along(top.genes4), function(i) {
  gset <- top.genes4[[i]]
  df <- gost(query = gset, organism = "hsapiens", sources = "GO:BP", evcodes = TRUE, )$result
  #df2 <- gost(query = gset, organism = 'gp__aJd5_EBf2_fFc', evcodes = TRUE)$result     #hallmarks from DBsig
  #df <- rbind(df1, df2)
  if (is.null(df)) return(NULL)
  df$factor <- paste0("factor_", i)
  return(df)
})
pathways4 <- do.call(rbind, pathways4)
pathways4$patient <- "4"



top.genes5 <- list()
for (factor in 1:15) {
  top.genes5[[factor]] <- select_top(OV5, factor)
}


#compute gene ontology (GO:BP) for each factor
pathways5 <- lapply(seq_along(top.genes5), function(i) {
  gset <- top.genes5[[i]]
  df <- gost(query = gset, organism = "hsapiens", sources = "GO:BP", evcodes = TRUE, )$result
  #df2 <- gost(query = gset, organism = 'gp__aJd5_EBf2_fFc', evcodes = TRUE)$result     #hallmarks from DBsig
  #df <- rbind(df1, df2)
  if (is.null(df)) return(NULL)
  df$factor <- paste0("factor_", i)
  return(df)
})
pathways5 <- do.call(rbind, pathways5)
pathways5$patient <- "5"


pathways <- rbind(pathways4, pathways5)
write.table(apply(pathways,2,as.character), "Desktop/IJC/TFG/STbl_enrichment.txt")


## Factors

for (i in paste0("factor_",1:15)) {
  p <- SpatialFeaturePlot(OV4A, features = i) + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0(i))
  pdf(paste0("Desktop/IJC/TFG/SuppFigures/Factors/OV4A_", i, ".pdf"), width = 10, height = 7)
  plot(p)
  dev.off()
}

for (i in paste0("factor_",1:15)) {
  p <- SpatialFeaturePlot(OV5, features = i, images = "OC") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle("")
  pdf(paste0("Desktop/IJC/TFG/SuppFigures/Factors/OV5A_", i, ".pdf"), width = 10, height = 7)
  plot(p)
  dev.off()
}

for (i in paste0("factor_",1:15)) {
  p <- SpatialFeaturePlot(OV5, features = i, images = "OC.1") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle("")
  pdf(paste0("Desktop/IJC/TFG/SuppFigures/Factors/OV5B_", i, ".pdf"), width = 10, height = 7)
  plot(p)
  dev.off()
}






##
#SFig 3
OV4A <- readRDS("Desktop/IJC/TFG/RDS/OV4A.rds")
OV5A <- readRDS("Desktop/IJC/TFG/RDS/OV5A.rds")
OV5B <- readRDS("Desktop/IJC/TFG/RDS/OV5B.rds")

unique(rownames(OV4A[["Deconvolution_Sub"]]))
celltypes <- unique(rownames(OV4A[["Deconvolution_Sub"]]))
names_Sub <- c("iCAFs2", "Endothelial cells", "CAFs", "iCAFs1", "Epithelial", "TAMs", "Proliferative Macrophages",
               "M1 Macrophages", "Proliferative Epithelial", "Proliferative T cells", "Proliferative Fibroblast", "Myofibroblast",
               "Dendritic cells", "Plasma B cells", "CD4+ T cells", "CD8+ T cells", "Natural Killer", "T regulatory cells",
               "Plasmacytoid Dendritic cells", "B cells", "Mast cells")
## OV4A
DefaultAssay(OV4A) <- "Deconvolution_Sub"
for (n in 1:length(celltypes)) {
  p <- SpatialFeaturePlot(OV4A, features = celltypes[n]) + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
    theme(legend.position="bottom", plot.title = element_text(hjust = 0.5, size = 25),
          legend.key.height = unit(0.4, 'cm'), legend.key.width  = unit(2.75, 'cm'),
          legend.direction="horizontal", legend.margin=margin(t=-15, ),
          legend.text = element_text(size=15)) + ggtitle(names_Sub[n])
  pdf(paste0("Desktop/IJC/TFG/SuppFigures/SFig3_OV4A_Sub_",celltypes[n], ".pdf"), width = 10, height = 7)
  plot(p)
  dev.off()  
}

## OV5A
DefaultAssay(OV5A) <- "Deconvolution_Sub"
for (n in 1:length(celltypes)) {
  p <- SpatialFeaturePlot(OV5A, features = celltypes[n]) + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
    theme(legend.position="bottom", plot.title = element_text(hjust = 0.5, size = 25),
          legend.key.height = unit(0.4, 'cm'), legend.key.width  = unit(2.75, 'cm'),
          legend.direction="horizontal", legend.margin=margin(t=-15),
          legend.text = element_text(size=15)) + ggtitle(names_Sub[n])
  pdf(paste0("Desktop/IJC/TFG/SuppFigures/SFig3_OV5A_Sub_",celltypes[n], ".pdf"), width = 10, height = 7)
  plot(p)
  dev.off()  
}

## OV5B
DefaultAssay(OV5B) <- "Deconvolution_Sub"
for (n in 1:length(celltypes)) {
  p <- SpatialFeaturePlot(OV5B, features = celltypes[n]) + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
    theme(legend.position="bottom", plot.title = element_text(hjust = 0.5, size = 25),
          legend.key.height = unit(0.4, 'cm'), legend.key.width  = unit(2.75, 'cm'),
          legend.direction="horizontal", legend.margin=margin(t=-15),
          legend.text = element_text(size=15)) + ggtitle(names_Sub[n])
  pdf(paste0("Desktop/IJC/TFG/SuppFigures/SFig3_OV5B_Sub_",celltypes[n], ".pdf"), width = 10, height = 7)
  plot(p)
  dev.off()  
}

####
df <- as.data.frame(matrix(nrow = 0, ncol = 3))
for (factor in 1:15) {
  for (gene in top.genes4[[factor]]){
    df <- rbind(df, c(gene, paste0("factor_", factor), "Patient 4")) 
  }
}
for (factor in 1:15) {
  for (gene in top.genes5[[factor]]){
    df <- rbind(df, c(gene, paste0("factor_", factor), "Patient 5")) 
  }
}

colnames(df) <- c("Gene", "Factor", "Patient")
write.table(df, "Desktop/IJC/TFG/SuppFigures/STbl_topgenes.r", sep = "\t")


View(samples.combined@misc$enriched_pathways$cerebro_seurat_enrichr$seurat_clusters)
write.table(samples.combined@misc$enriched_pathways$cerebro_seurat_enrichr$seurat_clusters,
            "Desktop/IJC/TFG/SuppFigures/STbl_enrichmentSC.txt", sep = "\t")

