################################################################################

#Dataset: Ductal_FFPE

################################################################################



#path of each file
infoTable <- data.frame(samples="Desktop/IJC/datasets/IGTP/4A/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5")
infoTable$spotfiles="Desktop/IJC/datasets/IGTP/4A/spatial/tissue_positions_list.csv"
infoTable$imgs="Desktop/IJC/datasets/IGTP/4A/spatial/tissue_hires_image.png"
infoTable$json="Desktop/IJC/datasets/IGTP/4A/spatial/scalefactors_json.json"

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
img <- Seurat::Read10X_Image(image.dir = 'Desktop/IJC/datasets/IGTP/4A/images/')
Seurat::DefaultAssay(object = img) <- 'RNA'
rownames(img@coordinates) <- paste0(rownames(img@coordinates), "_1")
img <- img[colnames(x = OV4A)]
OV4A[['OC']] <- img



#Normalization of the data returning all genes (SCT)
STobject <- SCTransform(STobject,return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.1, assay = "RNA")


#Deconvolution
## CARD -> script: deconvolution.r



## Seurat (probabilistic)

#Ligand Receptor activities
OV4A <- AddNiches(OV4A)

#Progeny pathways
OV4A <- AddProgeny(OV4A)
#Transcription factor activities

OV4A <- AddDorothea(OV4A)