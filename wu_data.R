library(Seurat)

# data = readRDS("mtx_complete.rds")
# meta = read.table("metadata.csv", header = T, row.names=1, sep = ",")


# so = CreateSeuratObject(data, meta.data = meta)
# so <- NormalizeData(so, assay = "RNA")
# so <- FindVariableFeatures(so, nfeatures = 2000)
# so <- ScaleData(so)
# so <- RunPCA(so)
# so <- RunUMAP(so, dims = 1:30)
# so <- FindNeighbors(so, dims = 1:30)
# so <- FindClusters(so, resolution = c(0,seq(0.1, 1, by=0.1)))
   

    
# saveRDS(so, "Wu_so.rds")

# so_list = SplitObject(so, split.by = "orig.ident")
# features = SelectIntegrationFeatures(so_list)
    # so_list <- lapply(X = l, FUN = function(x) {
        # x <- ScaleData(x, features = features)
        # x <- RunPCA(x, features = features)

# anchors <- FindIntegrationAnchors(object.list = so_list, 
                # anchor.features = features, reduction = "rpca")
    # int_set = IntegrateData(anchorset = anchors)
    # DefaultAssay(int_set) <- "integrated"

    ## Run the standard workflow for visualization and clustering
    # int_set <- ScaleData(int_set, verbose = FALSE)
    # int_set <- RunPCA(int_set, npcs = 30, verbose = FALSE)
    # int_set <- RunUMAP(int_set, reduction = "pca", dims = 1:30)
    # int_set <- FindNeighbors(int_set, reduction = "pca", dims = 1:30)
    # int_set <- FindClusters(int_set, resolution = c( 0.5)
    # saveRDS(int_set, paste("int_group_", i, sep=""))
	
	
# integration using CCA (canonical correlation analysis)
so = readRDS("Wu_so.rds")
so_list = SplitObject(so, split.by = "orig.ident")

so_list <- lapply(X = so_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = so_list)
anchors <- FindIntegrationAnchors(object.list = so_list, anchor.features = features)
int_CCA_set <- IntegrateData(anchorset = anchors)
DefaultAssay(int_CCA_set) <- "integrated"

# Run the standard workflow for visualization and clustering
int_CCA_set <- ScaleData(int_CCA_set, verbose = FALSE)
int_CCA_set <- RunPCA(int_CCA_set, npcs = 30, verbose = FALSE)
int_CCA_set <- RunUMAP(int_CCA_set, reduction = "pca", dims = 1:30)
int_CCA_set <- FindNeighbors(int_CCA_set, reduction = "pca", dims = 1:30)
int_CCA_set <- FindClusters(int_CCA_set, resolution = 0.5)

saveRDS(int_CCA_set, "wu_int_CCA_set.rds")

pdf("dimplot_int_cca.pdf")
DimPlot(int_CCA_set, group.by="celltype_major")
dev.off()


