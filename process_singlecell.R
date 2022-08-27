options(timeout = 1000)
## package 
library(Seurat)
library(tidyverse)
library(data.table)
library(SeuratDisk)
library(ggplot2)
library(patchwork)

set.seed(221)
raw_count_list <- readRDS("~/Desktop/rds_raw_count.rds")
#raw_count_list <- as.matrix(fread("~/Desktop/GSE189125_5prime_scRNAseq_seqbatchA_expression.txt"),rownames=1)
#raw_count_list <- as(raw_count_list, "dgCMatrix") 
metadata <- read_table("~/Desktop/GSE189125_5prime_scRNAseq_seqbatchA_metadata.txt")
rownames(metadata) <- metadata$barcode

state_seurat <- CreateSeuratObject(raw_count_list, project = "SeuratProject", assay = "RNA", meta.data = NULL)
state_seurat[["percent.mt"]] <- PercentageFeatureSet(state_seurat, pattern = "^MT-")

#data has been QCed
pdf(file= "~/Desktop/QC_metrics.pdf", width = 20, height = 11)
p <- VlnPlot(state_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt") , ncol = 3)
p
dev.off()

# Add metadata
state_seurat@meta.data$orig.ident <- rownames(state_seurat@meta.data)

metadata$orig.ident <- metadata$barcode
state_seurat@meta.data <- left_join(state_seurat@meta.data, metadata)
rownames(state_seurat@meta.data) <- state_seurat@meta.data$orig.ident

table(rownames(state_seurat@meta.data) == colnames(state_seurat@assays$RNA@counts)) 
# make sure metadata and assay fits well

# preprocess
tate_seurat <- state_seurat

  ##Variable Feature Selection
tate_seurat <-  FindVariableFeatures(tate_seurat,
                       selection.method = "vst",
                       nfeatures=500)
## scale
tate_seurat <- ScaleData(
  tate_seurat,
  features=rownames(tate_seurat),
  )
  ##Dimensionality Reduction
tate_seurat <-  RunPCA(
    tate_seurat,
    dims = 1:30
    ) 

bm <- tate_seurat

bm <- RunUMAP(
   bm,
   dims = 1:30,
    )
 
  ##Computing a cached neighbor index and Defining clusters
 bm <- FindNeighbors(
                object = bm,
                reduction = "pca",
                dims = 1:30,
                )
bm <- FindClusters(bm, resolution = 0.5)
DimPlot(bm)



##reference mapping
#reference <- LoadH5Seurat("~/Desktop/pbmc_multimodal.h5seurat")
#DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", 
        #label = TRUE, label.size = 3, repel = TRUE) + NoLegend()

#library(SeuratData)
#InstallData('pbmc3k')

#pbmc3k <- SCTransform(pbmc3k, verbose = FALSE)
#scTransform: sct normalization + scaling + find variable features

#anchors <- FindTransferAnchors(
  #reference = reference,
  #query = pbmc3k,
  #normalization.method = "SCT",
  #reference.reduction = "pca",
  #dims = 1:50
#)

#pbmc3k <- MapQuery(
  #anchorset = anchors,
  #query = pbmc3k,
  #reference = reference,
  #refdata = list(
    #celltype.l1 = "celltype.l1",
    #celltype.l2 = "celltype.l2",
    #predicted_ADT = "ADT"
  #),
  #reference.reduction = "pca", 
  #reduction.model = "wnn.umap"
#)
 

pdf(file = "~/Desktop/UMAP_Patient.pdf", width = 25, height = 18)
DimPlot(tate_seurat, group.by = "patient")
dev.off()
pdf(file = "~/Desktop/UMAP_ID.pdf", width = 25, height = 18)
DimPlot(tate_seurat, group.by = "ID", label = T, repel = T)
dev.off()

## add patient info into metadata
addInfo <- data.frame(
  "patient" = c('YUCARD',"YUENZO","YUFUB","YUFURL","YUGRUS","YUHERN","YUKEND","YUPIXEL",
                "YUROD","YUTAUR","YUTHEA","YUTORY","YUVARDO","YUALOE","YUHONEY","YUNANCY","YUCARD","YUENZO",
                "YUFUB","YUFURL","YUGRUS","YUHERN","YUKEND","YUPIXEL","YUROD","YUTAUR","YUTHEA","YUTORY","YUVARDO"),
  "grade" = c("3","1","0","3","2","2","0","3","2","3",'1',"3","3","3","4","3","3","1","0","3","2","2","0","3","2","3","1","3","3"),
  "Immunothrepy" = c("Combination","Anti-PD1","Combination","Combination","Combination",
                     "Combination","Anti-PD1","Combination","Combination","Combination",
                     "Combination","Combination","Combination","Combination","Combination","Combination","Combination",
                     "Anti-PD1","Combination","Combination","Combination","Combination","Anti-PD1","Combination","Combination",
                     "Combination","Combination","Combination","Combination")
)
addInfo

DimPlot(bm, group.by = "ID")

## remove the repeat patient
addInfo <- dplyr::distinct(addInfo, patient, .keep_all = T)

metadata <- bm@meta.data
metadata <- left_join(metadata,addInfo,by="patient")

bm@meta.data$grade <- metadata$grade
bm@meta.data$immunotherapy <- metadata$Immunothrepy
##Composition
#result
DimPlot(bm, group.by = "grade", label = F)
table(metadata$grade) %>% prop.table() * 100
table(metadata$grade, metadata$ID) %>% prop.table(margin = 2) * 100








