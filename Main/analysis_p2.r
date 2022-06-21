################################################################################

#Dataset: P1

################################################################################
library(Seurat)
library(STutility)

source("Desktop/IJC/TFG/scripts/utilitites/AddAssays.r")

## Sample A
#path of each file
infoTable <- data.frame(samples="Desktop/IJC/datasets/IGTP/5A/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5")
infoTable$spotfiles="Desktop/IJC/datasets/IGTP/5A/spatial/tissue_positions_list.csv"
infoTable$imgs="Desktop/IJC/datasets/IGTP/5A/spatial/tissue_hires_image.png"
infoTable$json="Desktop/IJC/datasets/IGTP/5A/spatial/scalefactors_json.json"

#create Seurat object through STutility, keep genes that have at least 5 counts in the whole tissue
STobject <- InputFromTable(infotable = infoTable, 
                           platform =  "Visium", minUMICountsPerGene = 5)


annotLookup <- read.table("Desktop/IJC/datasets/annotLook.txt", header = T, sep = "\t")

# remove mitochondrial, ribosomal and non-coding genes
protein_genes <- annotLookup$gene_name[annotLookup$gene_type %in% c("protein_coding", "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene", "IG_LV_gene", "IG_V_gene", "IG_J_gene", "IG_C_gene" , "IG_D_gene")]

selected_genes <- rownames(STobject)
#keep protein coding genes
selected_genes <- protein_genes[protein_genes %in% selected_genes]
#remove ribosomal genes
selected_genes <- selected_genes[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", selected_genes)]
#remove mitochondrial genes
selected_genes <- selected_genes[!grepl("^MT-", selected_genes)]


STobject <- SubsetSTData(STobject, features = selected_genes)


#Load the image in Staffli object to use ManualAnnotation
STobject <- LoadImages(STobject, verbose = TRUE) %>% MaskImages()
ImagePlot(STobject)

#Select spots that are out of the tissue and filter out
STobject<- ManualAnnotation(STobject)
STobject<- SubsetSTData(STobject, labels == "Default")

#Quality Control -> remove spots having a low number of Features
VlnPlot(STobject, features = c("nFeature_RNA", "nCount_RNA"))
STobject<- subset(STobject, nFeature_RNA > 0)
VlnPlot(STobject, features = c("nFeature_RNA", "nCount_RNA"))

#Load image in Seurat object to use Seurat plotting functions
img <- Seurat::Read10X_Image(image.dir = 'Desktop/IJC/datasets/IGTP/5A/images/')
Seurat::DefaultAssay(object = img) <- 'RNA'
rownames(img@coordinates) <- paste0(rownames(img@coordinates), "_1")
img <- img[colnames(x = OV5A)]
OV5A[['OC']] <- img

#Deconvolution
## CARD -> script: deconvolution.r


## Sample B
#path of each file
infoTable <- data.frame(samples="Desktop/IJC/datasets/IGTP/5B/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5")
infoTable$spotfiles="Desktop/IJC/datasets/IGTP/5B/spatial/tissue_positions_list.csv"
infoTable$imgs="Desktop/IJC/datasets/IGTP/5B/spatial/tissue_hires_image.png"
infoTable$json="Desktop/IJC/datasets/IGTP/5B/spatial/scalefactors_json.json"

#create Seurat object through STutility, keep genes that have at least 5 counts in the whole tissue
STobject <- InputFromTable(infotable = infoTable, 
                           platform =  "Visium", minUMICountsPerGene = 5)


annotLookup <- read.table("Desktop/IJC/datasets/annotLook.txt", header = T, sep = "\t")

# remove mitochondrial, ribosomal and non-coding genes
protein_genes <- annotLookup$gene_name[annotLookup$gene_type %in% c("protein_coding", "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene", "IG_LV_gene", "IG_V_gene", "IG_J_gene", "IG_C_gene" , "IG_D_gene")]

selected_genes <- rownames(STobject)
#keep protein coding genes
selected_genes <- protein_genes[protein_genes %in% selected_genes]
#remove ribosomal genes
selected_genes <- selected_genes[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", selected_genes)]
#remove mitochondrial genes
selected_genes <- selected_genes[!grepl("^MT-", selected_genes)]


STobject <- SubsetSTData(STobject, features = selected_genes)


#Load the image in Staffli object to use ManualAnnotation
STobject <- LoadImages(STobject, verbose = TRUE) %>% MaskImages()
ImagePlot(STobject)

#Select spots that are out of the tissue and filter out
STobject<- ManualAnnotation(STobject)
STobject<- SubsetSTData(STobject, labels == "Default")

#Quality Control -> remove spots having a low number of Features
VlnPlot(STobject, features = c("nFeature_RNA", "nCount_RNA"))
STobject<- subset(STobject, nFeature_RNA > 0)
VlnPlot(STobject, features = c("nFeature_RNA", "nCount_RNA"))

#Load image in Seurat object to use Seurat plotting functions
img <- Seurat::Read10X_Image(image.dir = 'Desktop/IJC/datasets/IGTP/5B/images/')
Seurat::DefaultAssay(object = img) <- 'RNA'
rownames(img@coordinates) <- paste0(rownames(img@coordinates), "_1")
img <- img[colnames(x = OV5B)]
OV5B[['OC']] <- img

#Deconvolution
## CARD -> script: deconvolution.r

#merge samples
OV5 <- merge(OV5A, OV5B)

#Normalization of the data returning all genes (SCT)
OV5 <- SCTransform(OV5, return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.1, assay = "RNA")


#Ligand Receptor activities
OV5 <- AddNiches(OV5)
#Progeny pathways
OV5 <- AddProgeny(OV5)
#Transcription factor activities
OV5 <- AddDorothea(OV5)

#NNMF
DefaultAssay(OV5) <- "SCT"
OV5 <- RunNMF(OV5, nfactors = 15)

#AssignClusters
OV5$interest2 <- NA

groups <- kmeans(OV5$factor_12, centers = 2)$cluster
if (mean(OV5@reductions$NMF@cell.embeddings[groups==1, 1]) >
    mean(OV5@reductions$NMF@cell.embeddings[groups==2, 1])) OV5$interest2[groups==1] <- "Cancer-TAMs crosstalk"
else {OV5$interest2[groups==2] <- "Cancer-TAMs crosstalk"}

