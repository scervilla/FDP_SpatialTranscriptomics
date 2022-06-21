library(Seurat)
library(ggplot2)
library(ggthemes)

source("plotfunct.r")


# 1B 
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
p <- SpatialFeaturePlot(OV4A, features = "nCount_RNA")  + scale_fill_gradientn("", colours = SpatialColors(100)) + 
  theme(legend.position = "right", legend.text = element_text(size = 18), plot.title = element_text(hjust = 0.5, size = 35)) + ggtitle(paste0("UMI Count"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/nCount_OV4A.pdf"), width = 10, height = 7)
plot(p)
dev.off()
p <- SpatialFeaturePlot(OV5A, features = "nCount_RNA")  + scale_fill_gradientn("", colours = SpatialColors(100)) + 
  theme(legend.position = "right", legend.text = element_text(size = 18), plot.title = element_text(hjust = 0.5, size = 35)) + ggtitle(paste0("UMI Count"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/nCount_OV5A.pdf"), width = 10, height = 7)
plot(p)
dev.off()
p <- SpatialFeaturePlot(OV5B, features = "nCount_RNA")  + scale_fill_gradientn("", colours = SpatialColors(100)) + 
  theme(legend.position = "right", legend.text = element_text(size = 18), plot.title = element_text(hjust = 0.5, size = 35)) + ggtitle(paste0("UMI Count"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/nCount_OV5B.pdf"), width = 10, height = 7)
plot(p)
dev.off()

#1C

#1C
## p53
DefaultAssay(OV4A) <- "progeny"
DefaultAssay(OV5) <- "progeny"
p <- SpatialFeaturePlot(OV4A, features = "p53") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("p53"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/p53_4A.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV5, features = "p53", images = "OC") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("p53"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/p53_5A.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV5, features = "p53", images = "OC.1") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("p53"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/p53_5B.pdf"), width = 10, height = 7)
plot(p)
dev.off()

#CLU
DefaultAssay(OV4A) <- "SCT"
DefaultAssay(OV5) <- "SCT"
p <- SpatialFeaturePlot(OV4A, features = "CLU") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("CLU"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/CLU_4A.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV5, features = "CLU", images = "OC") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("CLU"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/CLU_5A.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV5, features = "CLU", images = "OC.1") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("CLU"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/CLU_5B.pdf"), width = 10, height = 7)
plot(p)
dev.off()


#WFDC2
p <- SpatialFeaturePlot(OV4A, features = "WFDC2") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("WFDC2"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/WFDC2_4A.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV5, features = "WFDC2", images = "OC") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("WFDC2"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/WFDC2_5A.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV5, features = "WFDC2", images = "OC.1") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("WFDC2"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/WFDC2_5B.pdf"), width = 10, height = 7)
plot(p)
dev.off()


#COL1A1
p <- SpatialFeaturePlot(OV4A, features = "COL1A1") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("COL1A1"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/COL1A1_4A.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV5, features = "COL1A1", images = "OC") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("COL1A1"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/COL1A1_5A.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV5, features = "COL1A1", images = "OC.1") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("COL1A1"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/COL1A1_5B.pdf"), width = 10, height = 7)
plot(p)
dev.off()



#CXCL9
p <- SpatialFeaturePlot(OV4A, features = "CXCL9") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("CXCL9"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/CXCL9_4A.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV5, features = "CXCL9", images = "OC") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("CXCL9"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/CXCL9_5A.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV5, features = "CXCL9", images = "OC.1") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("CXCL9"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure1/CXCL9_5B.pdf"), width = 10, height = 7)
plot(p)
dev.off()

#2A
samples.combined <- RenameIdents(samples.combined, `Proliferative_Macrophages`="Proliferative Macrophages",
                                 `PlasmaB`= "Plasma B cells", `Mast_cells` = "Mast cells",
                                 `Fibroblast_Proliferative`= "Proliferative Fibroblast", `Epithelial_Proliferative`="Proliferaive Epithelial",
                                 `Bcells`="B cells", `CAFs_Immune1`="iCAFs1", `CAFs_Immune2`="iCAFs2",
                                 `ProliferativeTcells`="Proliferative T cells")
samples.combined$SubTypes <- as.character(Idents(samples.combined))
p <- DimPlot(samples.combined, reduction="umap", label=T, repel = T, cols = DiscretePalette(22)[-5]) 
pdf("Desktop/IJC/TFG/FinalFigures/Figure2/celltypes.pdf", width = 14, height = 5)
plot(p)
dev.off()

#2B

#2C

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



library(ggplot2)
library(ggthemes)
p <- ggplot(d) + geom_density(aes(x=TAMs, col=sample), alpha=0.3) + theme_classic() +
  labs(x="Cell type proportion", y="Density", title="TAMs") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20,),axis.title = element_text(size = 20),
        legend.position = "none") + scale_color_manual(values = colorblind_pal()(3))
pdf("Desktop/IJC/TFG/FinalFigures/Figure2/dist_TAMs.pdf", width = 4, height = 6)
plot(p)
dev.off()
p <- ggplot(d) + geom_density(aes(x=Epithelial, color=sample, col=sample), alpha=0.3) + theme_classic() +
  labs(x="Cell type proportion", y="Density", title="Epithelial") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20),axis.title = element_text(size = 20), legend.position = "none") + scale_color_manual(values = colorblind_pal()(3))
pdf("Desktop/IJC/TFG/FinalFigures/Figure2/dist_Epithelial.pdf", width = 4, height = 6)
plot(p)
dev.off()


p <- ggplot(d) + geom_density(aes(x=CAFs, color=sample, col=sample), alpha=0.3) + theme_classic() +
  labs(x="Cell type proportion", y="Density", title="CAFs") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20 ),axis.title = element_text(size = 20), legend.position = "none") + scale_color_manual(values = colorblind_pal()(3))
pdf("Desktop/IJC/TFG/FinalFigures/Figure2/dist_CAFs.pdf", width = 4, height = 6)
plot(p)
dev.off()

p <- ggplot(d) + geom_density(aes(x=Tcells, color=sample, col=sample), alpha=0.3) + theme_classic() +
  labs(x="Cell type proportion", y="Density", title="T cells") + ylim(0,20)+ 
  theme(plot.title = element_text(hjust = 0.5, size = 20),axis.title = element_text(size = 20), legend.position = "none") + scale_color_manual(values = colorblind_pal()(3))
pdf("Desktop/IJC/TFG/FinalFigures/Figure2/dist_Tcell.pdf", width = 4, height = 6)
plot(p)
dev.off()

#2E
corr4 <- lapply(rownames(OV4A[["Deconvolution_Major"]]), function(celltype){
  sapply(rownames(OV4A[["progeny"]]), function(pathway){
    cor(OV4A[["Deconvolution_Major"]]@data[celltype,], OV4A[["progeny"]]@data[pathway,])
  })
})
corr4 <- do.call(rbind, corr4)
rownames(corr4) <- rownames(OV4A[["Deconvolution_Major"]])
pheatmap(t(corr4), angle_col = 315, color =   rev(paletteer::paletteer_c("grDevices::RdBu", 60) ),fontsize = 12)


DefaultAssay(OV5) <- "SCT"
pathways <- progeny::progeny(as.matrix(GetAssayData(OV5)), scale=TRUE, 
                             organism="Human", top = 1000, perm = 1)
OV5[["progeny"]] <- CreateAssayObject(t(pathways))

corr5 <- lapply(rownames(OV5[["Deconvolution_Major"]]), function(celltype){
  sapply(rownames(OV5[["progeny"]]), function(pathway){
    cor(OV5[["Deconvolution_Major"]]@data[celltype,], OV5[["progeny"]]@data[pathway,])
  })
})
corr5 <- do.call(rbind, corr5)
rownames(corr5) <- rownames(OV5[["Deconvolution_Major"]])
pheatmap(t(corr5), angle_col = 315, color =   rev(paletteer::paletteer_c("grDevices::RdBu", 60) ),fontsize = 12)

#3A
p <- SpatialFeaturePlot(OV5, features = "factor_14", images = "OC") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("factor_14"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure3/facto14_5A.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV5, features = "factor_14", images = "OC.1") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("factor_14"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure3/facto14_5B.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV4A, features = "factor_12") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("factor_12"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure3/facto12_4A.pdf"), width = 10, height = 7)
plot(p)
dev.off()

#3B

df4 <- data.frame(genes=factor(genes4, levels = rev(genes4)), weight=OV4A@reductions$NMF@feature.loadings[genes4,12])
df5 <- data.frame(genes=factor(genes5, levels = genes5), weight=OV5@reductions$NMF@feature.loadings[genes5,14])


g <- intersect(genes4, genes5)
f <- c(rep(4,100))
f[which(genes4 %in% g)] <- 7

ggplot(df4, aes(y=genes, x=weight)) + geom_bar(stat = "identity", width = 0.5) +
  theme_classic() + theme(axis.ticks = element_blank(), axis.text.y = element_text(size = rev(f)),
                          plot.title = element_text(hjust = 0.5)) +
  labs(x = "Gene weight", y="") + ggtitle("factor_12 (P4)")



f <- c(rep(5,100))
f[which(genes5 %in% g)] <- 12

ggplot(df5, aes(y=genes, x=weight)) + geom_bar(stat = "identity", width = 0.5) +
  theme_classic() + theme(axis.ticks = element_blank(), axis.text.y = element_text(size = rev(f)),
                          plot.title = element_text(hjust = 0.5)) +
  labs(x = "Gene weight", y="") + ggtitle("factor_14 (P5)")

hjust <- rep(0, 100)
hjust[seq(1, 100, by=2)] <- 0.2
ggplot(df5, aes(x=genes, y=weight)) + geom_bar(stat = "identity", width = 0.5) +
  theme_classic() + theme(axis.ticks = element_blank(), axis.text.x = element_text(size = f, angle = -45, hjust = 0),
                          plot.title = element_text(hjust = 0.5),) +
  labs(y = "Gene weight", x="") + ggtitle("factor_14 (P5)") + scale_x_discrete(guide = guide_axis(n.dodge=2))

samples.combined$SubTypes <- Idents(samples.combined)
samples.combined$MajorTypes <- as.character(samples.combined$MajorTypes)

Idents(samples.combined) <- samples.combined$MajorTypes
samples.combined <- RenameIdents(samples.combined, `Proliferative_Macrophages`="Proliferative Macrophages",
                                 `PlasmaB`= "Plasma B cells", `Mast_cells` = "Mast cells",
                                 `Fibroblast_Proliferative`= "Proliferative Fibroblast", `Epithelial_Proliferative`="Proliferaive Epithelial",
                                 `Bcells`="B cells")
samples.combined$MajorTypes <- as.character(Idents(samples.combined))
p <- DotPlot(samples.combined, features = Reduce(intersect, l), cols = "RdBu", group.by = "MajorTypes") + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = -45, hjust = 0)) 
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figures3/DotPlot.pdf"), width = 10, height = 7)
plot(p)
dev.off()

#3C
p <- SpatialFeaturePlot(OV4A, features = "SLPI") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 30),      legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("SLPI"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure3/OV4A_SLPI.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV5, features = "SLPI", images = "OC") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 30),      legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("SLPI"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure3/OV5A_SLPI.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV5, features = "SLPI", images = "OC.1") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 30),      legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("SLPI"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure3/OV5B_SLPI.pdf"), width = 10, height = 7)
plot(p)
dev.off()
p <- SpatialFeaturePlot(OV4A, features = "MIF") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 30),       legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("MIF"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure3/MIF_4A.pdf"), width = 10, height = 7)
plot(p)
dev.off()
p <- SpatialFeaturePlot(OV5, features = "MIF", images="OC.1") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 30),       legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("MIF"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure3/MIF_5B.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV5, features = "MIF", images="OC") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 30), legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("MIF"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure3/MIF_5A.pdf"), width = 10, height = 7)
plot(p)
dev.off()
#3D
df <- as.data.frame(do.call(rbind, univ_results))
df$group <- factor(covariates, levels = rev(covariates))

df$Label <- NA
df$Label[df$p.value<0.001]<-'***'
df$Label[df$p.value<0.01 & is.na(df$Label)]<-'**'
df$Label[df$p.value<0.05 & is.na(df$Label)]<-'*'
df$Label[is.na(df$Label)] <-""

df$p.value <- format.pval(df$p.value)
df$p <- paste0(df$p.value, df$Label)

p <- ggplot(df, aes(y=group, x=HR)) +
  geom_point(shape=18, size=4, fill="white") +
  geom_errorbar(width=0.3, aes(xmin=CIlower, xmax=CIupper), colour="black") +
  annotate(geom="text", label=df$p, x = max(df$CIupper)+0.25, y=18:1, size=3) +
  annotate(geom="text", label="P-value", x = max(df$CIupper)+0.25, y=19, size =3.5) +
  geom_vline(xintercept=1, linetype="dashed", color ="gray") +
  labs(y="High expression group", x = "Hazard Ratio") +
  coord_cartesian(clip = 'off') +
  theme_tufte()

pdf("Desktop/HR.pdf", width = 6, height = 7)
plot(p)
dev.off()

#4A
p <- SpatialFeaturePlot(OV4A, features = "HLA-A—LILRB1") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("HLA-A—LILRB1"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure4/HLA-A—LILRB1.pdf"), width = 10, height = 7)
plot(p)
dev.off()

#4B
labels <- c("Distant TAMs", "Cancer - TAMs crosstalk", "Cancer")
palette <- c("#E69F00", "#56B4E9", "#009E73")
names(palette) <- labels
p <- SpatialDimPlot(OV4A, group.by = "interest2", cols = c("#E69F00", "#56B4E9", "#009E73"))  + 
  scale_fill_manual(labels = labels, values = palette) + 
  theme(legend.key = element_blank(), legend.direction = "horizontal", legend.position = "top", legend.title = element_blank(), legend.text=element_text(size=22)) + 
  guides(fill = guide_legend(override.aes = list(size = 6))) 
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure4/Clusters.pdf"), width = 10, height = 7)
plot(p)
dev.off()
#4C
Idents(OV4A) <- OV4A$interest2
markers <- FindAllMarkers(OV4A, test.use = "roc", only.pos = T)
top <- markers %>% group_by(cluster) %>% top_n(15, myAUC)
p <- DoHeatmap(OV4A[,!is.na(OV4A$interest2)], features = top$gene, group.colors = c("#E69F00", "#56B4E9", "#009E73")) + 
  scale_color_manual(values =  c("#E69F00", "#56B4E9", "#009E73"))
pdf("Desktop/IJC/TFG/FinalFigures/Figure4/LRheatmap.pdf", width = 7, height = 9)
plot(p)
dev.off()

#4D
p <- SpatialFeaturePlot(OV4A, features = "APOE—TREM2") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("APOE—TREM2"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure4/APOE—TREM2.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV4A, features = "MIF—CD74") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("MIF—CD74"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure4/MIF—CD74.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV4A, features = "EPCAM—CLDN1") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("EPCAM—CLDN1"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure4/EPCAM—CLDN1.pdf"), width = 10, height = 7)
plot(p)
dev.off()

#5A
p <- SpatialFeaturePlot(OV4A, features = "MYC") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("MYC"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure4/MYC.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV5, features = "MYC", images="OC.1") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("MYC"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure4/MYC_5B.pdf"), width = 10, height = 7)
plot(p)
#5B
p <- SpatialFeaturePlot(OV4A, features = "HIF1A") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("HIF1A"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure4/HIF1A_4A.pdf"), width = 10, height = 7)
plot(p)
dev.off()

p <- SpatialFeaturePlot(OV5, features = "HIF1A", images = "OC.1") + scale_fill_gradientn("", colours = viridis::viridis(100)) + 
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 25),         legend.key.height = unit(3, 'cm'), legend.key.width  = unit(0.4, 'cm'), legend.text = element_text(size=25),         legend.margin=margin(l=-10)) + ggtitle(paste0("HIF1A"))
pdf(paste0("Desktop/IJC/TFG/FinalFigures/Figure4/HIF1A_5B.pdf"), width = 10, height = 7)
plot(p)
dev.off()
#5C

