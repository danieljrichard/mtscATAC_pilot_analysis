####
##Processing MGATK data
##and Seurat objects
##July 17th 2024
####

##working through excellent tutorial on signac website: https://stuartlab.org/signac/articles/mito

library(Seurat)
library(ggplot2)
library(patchwork)
library(Signac)

library(Seurat)
library(ggplot2)
library(patchwork)
library(Signac)
library(EnsDb.Mmusculus.v79)
#annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
#seqlevelsStyle(annotations) <- 'UCSC'
#genome(annotations) <- "mm10"
#saveRDS(annotations, "mm10_seurat_annotations.rds")

annotations <- readRDS("../C1/mm10_seurat_annotations.rds")
counts <- Read10X_h5("filtered_peak_bc_matrix.h5")

grange.counts <- StringToGRanges(rownames(counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- counts[as.vector(grange.use), ]

chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments = 'fragments.tsv.gz',
  min.cells = 10,
  min.features = 200, annotation = annotations)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"#, meta.data = metadata
)

svz <- pbmc

##read in MGATK results
mgatk_dir <- "C1_test_MGATK_FINAL/final"

mito.data <- ReadMGATK(dir = mgatk_dir)

##cut down to cells passing QC 
mito <- CreateAssayObject(counts = mito.data$counts)
mito <- subset(mito, cells = Cells(svz))
svz[["mito"]] <- mito
svz <- AddMetaData(svz, metadata = mito.data$depth[Cells(mito), ], col.name = "mtDNA_depth")

##I end up using ArchR clusters instead.
svz <- RunTFIDF(svz)
svz <- FindTopFeatures(svz, min.cutoff = 10)
svz <- RunSVD(svz)
svz <- RunUMAP(svz, reduction = "lsi", dims = 2:50)
svz <- FindNeighbors(svz, reduction = "lsi", dims = 2:50)
svz <- FindClusters(svz, resolution = 0.7, algorithm = 3)

pdf("C1_dimplot_UMAP.pdf", width = 12, height = 12)
DimPlot(svz, label = TRUE) + NoLegend()
dev.off()

gene.activities <- GeneActivity(svz)

# add to the Seurat object as a new assay
svz[['RNA']] <- CreateAssayObject(counts = gene.activities)

svz <- NormalizeData(
  object = svz,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(svz$nCount_RNA))

DefaultAssay(svz) <- "peaks"

###############
##For simplicity's sake, I've omitted the code pertaining to label transfer
##from a reference Multiome dataset.

saveRDS(svz, "C1_celltype_LABELTRANSFER_withmito.rds")
