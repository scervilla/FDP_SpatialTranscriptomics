################################################################################
#                          Functions to add assays
################################################################################

#### pathway activities ####
AddProgeny <- function(STobject) {
  expr <- progeny::progeny(as.matrix(t(GetAssayData(STobject, slot = "SCT"))), top=1000)
  STobject[["progeny"]] <- CreateAssayObject(t(expr))
  return(STobject)
}

#### Ligand-receptor activity ####
AddNiches <- function(STobject) {
  STobject@meta.data$x <- STobject@images[[1]]@coordinates$row
  STobject@meta.data$y <- STobject@images[[1]]@coordinates$col
  DefaultAssay(STobject) <- "SCT"
  STobject <- SeuratWrappers::RunALRA(STobject)
  
  
  NICHES_output <- RunNICHES(object = STobject,
                             LR.database = "omnipath",
                             species = "human",
                             assay = "alra",
                             position.x = 'x',
                             position.y = 'y',
                             rad.set = 2, # Geometry dependent
                             min.cells.per.ident = 0,
                             min.cells.per.gene = NULL,
                             meta.data.to.map = c('bayes.cluster'),
                             CellToCell = F,CellToSystem = F,SystemToCell = F,
                             CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = T)
  niche <- NICHES_output[['NeighborhoodToCell']]
  Idents(niche) <- niche[['ReceivingType']]
  niches.data <- GetAssayData(object =  niche[['NeighborhoodToCell']], slot = 'data')
  colnames(niches.data) <- niche[['ReceivingCell']]$ReceivingCell
  STobject[["NeighborhoodToCell"]] <- CreateAssayObject(data = niches.data )
  DefaultAssay(STobject) <- "NeighborhoodToCell"
  STobject <- ScaleData(STobject)
  return(STobject)
}

#### Transcription factors ####
AddDorothea <- function(STobject) {
  #Run dorothea
  df2regulon = function(df) {
    regulon = df %>%
      split(.$tf) %>%
      map(function(dat) {
        tf = dat %>% distinct(tf) %>% pull()
        targets = setNames(dat$mor, dat$target)
        likelihood = dat$likelihood
        list(tfmode =targets, likelihood = likelihood)
      })
    return(regulon)
  }
  
  ## <<DOROTHEA>>  
  
  dorothea_regulon_human = read_csv("https://raw.githubusercontent.com/saezlab/ConservedFootprints/master/data/dorothea_benchmark/regulons/dorothea_regulon_human_v1.csv")
  
  # We obtain the regulons based on interactions with confidence level A, B and C
  DefaultAssay(STobject) <- "SCT"
  
  regulon = dorothea_regulon_human %>%
    dplyr::filter(confidence %in% c("A","B","C")) %>%
    
    df2regulon()
  
  
  tf_act_mat = viper(eset = as.matrix(STobject[["SCT"]]@data),
                     regulon = regulon, nes = TRUE,
                     method = "scale", minsize = 4,
                     eset.filter = FALSE,
                     verbose = FALSE, cores = 14)  
  
  
  
  STobject[['dorothea']] = CreateAssayObject(counts = tf_act_mat)
  DefaultAssay(STobject) <- "dorothea"
  rm(tf_act_mat)
  STobject <- ScaleData(STobject)
  return(STobject)
}

