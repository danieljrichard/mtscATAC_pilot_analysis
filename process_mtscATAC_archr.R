######
##Processing mtscATAC data with ArchR
##July 2024

##need to have installed ArchR

library(ArchR)
##defaults to half threads

addArchRGenome("mm10")

addArchRThreads(threads = 16) 


frag_files <- system("find . -name '*frag*.gz'", intern = T)
sample_name <- "C1"

##mgatk output
het_matrix_file <- "C1_test1_RAW.cell_heteroplasmic_df.tsv.gz"

ArrowFiles <- createArrowFiles(inputFiles = frag_files,
    sampleNames = sample_name,
    filterTSS = 4, 
    filterFrags = 1000,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE) 

#ArrowFiles will contain a list of .arrow files, which are saved to the current work directory

batch_files <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = paste0(sample_name, "_ARCHR"),
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

saveArchRProject(ArchRProj = batch_files, outputDirectory = paste0("Save-", sample_name, "_ARCHR"), load = FALSE)

##mostly left to ArchR tutorial defaults.
batch_files <- addIterativeLSI(
    ArchRProj = batch_files,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

batch_files <- addClusters(input = batch_files, reducedDims = "IterativeLSI")
batch_files <- addUMAP(ArchRProj = batch_files, reducedDims = "IterativeLSI")

p2 <- plotEmbedding(ArchRProj = batch_files, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

library(ggplot2)
pdf(paste0("ARCHR_", sample_name, "_clustering.pdf"), width = 10, height = 10)
print(p2)
dev.off()

saveArchRProject(ArchRProj = batch_files, outputDirectory = paste0("Save-", sample_name, "_ARCHR"), load = FALSE)
###and later:
#batch_files <- loadArchRProject(path = paste0("Save-", sample_name, "_ARCHR"))

####
##Using a previous multiome dataset, I assigned labels using Seurat and FindAnchors on the ATAC data. Majority of cells had > 0.6 prediction for the assigned celltype.
metadata <- fread("C1_celltype_LABELTRANSFER_withmito_metadata.txt")

##align seurat data
seurat_cells <- metadata$CB
seurat_cells <- paste0(sample_name, "#", seurat_cells)
seurat_cells <- gsub("_[0-9]", "", seurat_cells)

cellname_overlap <- intersect(seurat_cells, batch_files$cellNames)
##some cells may be lost as the Seurat script I used applied different quality cutoffs.

##then we'll subset both objects

subset_arch <- batch_files[cellname_overlap,]
metadata$CELLNAME <- seurat_cells
subset_seurat_meta <- metadata[seurat_cells %in% cellname_overlap, ]

mapper <- subset_seurat_meta$TRANSFER
names(mapper) <- subset_seurat_meta$CELLNAME

map_arch_cell <- unlist(lapply(subset_arch$cellNames, function(x) mapper[[x]]))
names(map_arch_cell) <- subset_arch$cellNames

##add in transferred labels
subset_arch <- addCellColData(
    ArchRProj = subset_arch,
    data = map_arch_cell,
    name = "LABEL_TRANSFER",
    cells = names(map_arch_cell)
)

####
##see if ArchrR clusters are largely comprised of similar cell-types

labeltrans <- plotEmbedding(ArchRProj = subset_arch, colorBy = "cellColData", name = "LABEL_TRANSFER", embedding = "UMAP")
clust_plot <- plotEmbedding(ArchRProj = subset_arch, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

pdf(paste0("ARCHR_", sample_name, "_attempt_cluster_labeltransfer.pdf"), width = 12, height = 12)
print(ggAlignPlots(labeltrans, clust_plot, type = "h"))
dev.off()

##Now, adding the MGATK output to the ArchR object

het_matrix <- as.data.frame(fread(het_matrix_file))
het_matrix$V1 <- paste0(sample_name, "#", het_matrix$V1)
het_matrix$V1 <- gsub("_[0-9]", "", het_matrix$V1)
rownames(het_matrix) <- het_matrix$V1
het_matrix <- het_matrix[, 2:dim(het_matrix)[2]]

het_matrix_sub <- het_matrix[rownames(subset_arch@cellColData),]
het_matrix_sub[is.na(het_matrix_sub)] <- 0
##changing variable names, as downstream UMAP plotting doesn't like special characters or numbers at the start of variable names
colnames(het_matrix_sub) <- gsub(">", '', colnames(het_matrix_sub))
colnames(het_matrix_sub) <- paste0("M", colnames(het_matrix_sub))
super_col <- subset_arch@cellColData
super_col <- cbind(super_col, het_matrix_sub)

subset_arch@cellColData <- super_col

######
##Prioritizing variants that are more specific to a few clusters.
##Calculating the average heteroplasmy for a variant across a cluster.

arch_meta <- subset_arch@cellColData

##we'll go variant-by-variant

###We want to count the number of cells having the mutation detected
##as well as the average heteroplasmy of a mutation within a cluster.

variantwise_average_rows <- list()
variantwise_hits_rows <- list()
cluster_split <- split(arch_meta, arch_meta$Clusters)
for (var in colnames(het_matrix_sub)) {
    clusterwise_average <- unlist(lapply(cluster_split, function(x) mean(x[, var], na.rm = T)))
    clusterwise_hits <- unlist(lapply(cluster_split, function(x) length(which(x[, var] != 0))/dim(x)[1]))
    variantwise_average_rows[[var]] <- clusterwise_average
    variantwise_hits_rows[[var]] <- clusterwise_hits
    print(var)
}

var_avg_frame <- do.call("rbind", variantwise_average_rows)
rownames(var_avg_frame) <- colnames(het_matrix_sub)

#####
##Sort variants based on their skew
var_avg_frame_skew <- apply(var_avg_frame, 1, function(x) max(x) - mean(x, na.rm =T))
sorted_indices <- order(var_avg_frame_skew, decreasing = TRUE)
var_avg_frame_sorted <- var_avg_frame[sorted_indices,]

var_hits_frame <- do.call("rbind", variantwise_hits_rows)
rownames(var_hits_frame) <- colnames(het_matrix_sub)

hits_frame_skew <- apply(var_hits_frame, 1, function(x) max(x) - mean(x, na.rm =T))
sorted_indices2 <- order(hits_frame_skew, decreasing = TRUE)
var_hits_frame_sorted <- var_hits_frame[sorted_indices2,]

library(ComplexHeatmap)
library(circlize)

avg_frame_heat_SORTED <- Heatmap(var_avg_frame_sorted, name = "Sorted Average heteroplasmy",
    col = colorRamp2(c(0, max(var_avg_frame)), c("white", "red")), show_row_names = T, show_column_names = T, cluster_rows = F)

var_frame_heat_SORTED <- Heatmap(var_hits_frame_sorted, name = "Sorted Percent Cells hit",
    col = colorRamp2(c(0, max(var_hits_frame)), c("white", "red")), show_row_names = T, show_column_names = T, cluster_rows = F)

avg_frame_heat_scale_SORTED <- Heatmap(var_avg_frame_sorted/rowSums(var_avg_frame_sorted), name = "Scaled Average heteroplasmy",
    col = colorRamp2(c(0, 1), c("white", "red")), show_row_names = T, show_column_names = T, cluster_rows = FALSE)

var_frame_heat_scale_SORTED <- Heatmap(var_hits_frame_sorted / rowSums(var_hits_frame_sorted), name = "Sorted Scaled Percent Cells hit",
    col = colorRamp2(c(0, 1), c("white", "red")), show_row_names = T, show_column_names = T, cluster_rows = FALSE)

pdf(paste0("ARCHR_", sample_name, "_clusterwise_hetero_heatmaps.pdf"), width = 10, height = 10)
print(avg_frame_heat_SORTED)
print(var_frame_heat_SORTED)
print(avg_frame_heat_scale_SORTED)
print(var_frame_heat_scale_SORTED)

dev.off()

######
##Explicitly looking for variants which are at low average heteroplasmy in a cell cluster from a different germ layer (immune cells)
##

to_exclude <- "C1"

var_hits_frame_frame <- as.data.frame(var_hits_frame)

exclusive_hits <- var_hits_frame_frame[var_hits_frame_frame[, to_exclude] < 0.05,]
exclusive_hits <- exclusive_hits[order(exclusive_hits[,to_exclude]),]
##leaves 12 variants.

pdf(paste0("ARCHR_", sample_name, "_clusterwise_hetero_heatmaps_EXCLUSIVE.pdf"), width = 10, height = 10)
avg_frame_heat_EX <- Heatmap(var_avg_frame[rownames(exclusive_hits),], name = "Average heteroplasmy",
    col = colorRamp2(c(0, max(var_avg_frame[rownames(exclusive_hits),])), c("white", "red")), show_row_names = T, show_column_names = T)

var_frame_heat_EX <- Heatmap(var_hits_frame[rownames(exclusive_hits),], name = "Percent Cells hit",
    col = colorRamp2(c(0, max(var_hits_frame[rownames(exclusive_hits),])), c("white", "red")), show_row_names = T, show_column_names = T)
print(avg_frame_heat_EX)
print(var_frame_heat_EX)
dev.off()

##Visualize these variants 
variants_of_interest <- rownames(exclusive_hits)

pdf(paste0("ARCHR_", sample_name, "_clusterwise_heteroplasmy_UMAPS_EXCLUSIVE.pdf"), width = 10, height = 10)
raster_set <- list()
for (var in variants_of_interest) {
    curr_plot <- plotEmbedding(ArchRProj = subset_arch, colorBy = "cellColData", name = var, embedding = "UMAP")
    print(curr_plot)
        raster_set[[var]] <- curr_plot
print(var)
}
dev.off()

library(gridExtra)
pdf(paste0("ARCHR_", sample_name, "_clusterwise_heteroplasmy_UMAPS_EXCLUSIVE_COMPLETE.pdf"), width = 20, height = 20)
grid.arrange(grobs = raster_set, nrow  = 4)
dev.off()

####
##For each variant, show the breakdown of clusters and assigned celltypes

arch_meta <- subset_arch@cellColData

variant_out <- list()
variant_out_cluster <- list()
for (var in variants_of_interest) {
    curr_meta <- arch_meta[arch_meta[, var] != 0,]
    curr_table <- unlist(lapply(unique(arch_meta$LABEL_TRANSFER), function(x) length(which(curr_meta$LABEL_TRANSFER == x))))
    curr_table_cluster <- unlist(lapply(unique(arch_meta$Clusters), function(x) length(which(curr_meta$Clusters == x))))
    curr_frame <- data.frame(variant = var, cell_type = unique(arch_meta$LABEL_TRANSFER), count = curr_table)
    curr_frame_cluster <- data.frame(variant = var, cell_type = unique(arch_meta$Clusters), count = curr_table_cluster)
    variant_out[[var]] <- curr_frame
    variant_out_cluster[[var]] <- curr_frame_cluster
    print(var)
}

final_cell_summary <- do.call("rbind", variant_out)
final_cluster_summary <- do.call("rbind", variant_out_cluster)

library(ggplot2)

pdf(paste0("ARCHR_", sample_name, "_exclusive_variants_analysis_celltype_breakdowns.pdf"), width = 6, height = 6)
print(ggplot(final_cell_summary, aes(fill=cell_type, y=count, x=variant)) + 
    geom_bar(position="fill", stat="identity") + ggtitle("Variant assignments"))
print(ggplot(final_cell_summary, aes(fill=cell_type, y=count, x=variant)) + 
    geom_bar(position="stack", stat="identity") + ggtitle("Variant assignments"))
dev.off()

pdf(paste0("ARCHR_", sample_name, "_exclusive_variants_analysis_cluster_breakdowns.pdf"), width = 6, height = 6)
print(ggplot(final_cluster_summary, aes(fill=cell_type, y=count, x=variant)) + 
    geom_bar(position="fill", stat="identity") + ggtitle("Variant assignments"))
print(ggplot(final_cluster_summary, aes(fill=cell_type, y=count, x=variant)) + 
    geom_bar(position="stack", stat="identity") + ggtitle("Variant assignments"))
dev.off()

#########
##July 31st 2024
##Finally, can we define clonotypes using Seurat/Signac's FindClonotypes function

##in a separate script, loaded in mgatk data and added to Seurat object for this 10X run.
library(Signac)
svz <- readRDS("C1_celltype_LABELTRANSFER_withmito.rds")

##revert back to original variant formatting
variants_of_interest_convert <- gsub("M", "", variants_of_interest)
variants_of_interest_convert <- unlist(lapply(variants_of_interest_convert, function(x)
    paste0(substr(x, 1, nchar(x) - 1), ">", substr(x, nchar(x), nchar(x)))))

####
##we can really only define clonotypes for cells having at least one of these variants.
##otherwise, the clonotypes devolve into 1 group of 'no mutations' and 1 group with 'some mutations'.
arch_meta$CB <- gsub(paste0(sample_name, "#"), "", rownames(arch_meta))
arch_meta_subset <- arch_meta[, variants_of_interest]
rownames(arch_meta_subset) <- arch_meta$CB
per_cell <- apply(arch_meta_subset, 1, function(x) length(which(x != 0)))
arch_meta_cut <- arch_meta_subset[per_cell >= 1,]

##cut down seurat object to match.
svz_sub <- svz[, rownames(arch_meta_cut)]

##calculate allele frequencies for variants of interest.
svz_sub <- AlleleFreq(
  object = svz_sub,
  variants = variants_of_interest_convert,
  assay = "mito"
)

svz_sub[["alleles"]]

DefaultAssay(svz_sub) <- "alleles"

svz_sub <- FindClonotypes(svz_sub)

svz_sub$CLONOTYPE <- Idents(svz_sub)

##to visualize variant heatmaps, add in cluster information.
svz_sub$ARCH_CLUSTER <- unlist(lapply(colnames(svz_sub), function(x) arch_meta$Clusters[arch_meta$CB == x]))

table(Idents(svz_sub))

pdf(paste0("ARCHR_", sample_name, "_test_clonotype_image.pdf"), width = 12, height = 12)
Idents(svz_sub) <- svz_sub$CLONOTYPE
print(DoHeatmap(svz_sub, features = VariableFeatures(svz_sub), slot = "data", disp.max = 0.1) +
  scale_fill_viridis_c())

Idents(svz_sub) <- svz_sub$ARCH_CLUSTER
print(DoHeatmap(svz_sub, features = VariableFeatures(svz_sub), slot = "data", disp.max = 0.1) +
  scale_fill_viridis_c())
dev.off()

##what are the cluster and cell-type breakdowns of these clonotypes

cluster_breakdown <- list()
svz_meta <- svz_sub@meta.data

svz_meta$CLONOTYPE <- as.character(svz_meta$CLONOTYPE)
for (x in unique(svz_meta$CLONOTYPE)) {
    curr_cells <- svz_meta[svz_meta$CLONOTYPE == x,]
    curr_table <- unlist(lapply(unique(svz_meta$SIMPLE), function(x) length(which(curr_cells$SIMPLE == x))))
    curr_frame <- data.frame(CLONOTYPE = x, cell_type = unique(svz_meta$SIMPLE), count = curr_table)
    curr_frame$perc <- paste0(round((curr_frame$count / sum(curr_frame$count))*100), "%")
    cluster_breakdown[[x]] <- curr_frame
}

final_summary <- do.call("rbind", cluster_breakdown)
write.csv(final_summary, paste0("ARCHR_", sample_name, "_clonotype_analysis_celltype_breakdowns.csv"))

library(ggplot2)

pdf(paste0("ARCHR_", sample_name, "_clonotype_analysis_celltype_breakdowns.pdf"), width = 6, height = 6)
print(ggplot(final_summary, aes(fill=cell_type, y=count, x=CLONOTYPE)) + 
    geom_bar(position="fill", stat="identity") + ggtitle(paste0(sample_name, " Clonotype assignments")))
print(ggplot(final_summary, aes(fill=cell_type, y=count, x=CLONOTYPE)) + 
    geom_bar(position="stack", stat="identity") + ggtitle(paste0(sample_name, " Clonotype assignments")))
dev.off()

########
##And finally, we add this information BACK into the archr object

arch_meta_match <- arch_meta[arch_meta$CB %in% colnames(svz_sub),]

subset_arch_clono <- subset_arch[rownames(arch_meta_match),]
subset_arch_clono_meta <- subset_arch_clono@cellColData
subset_arch_clono_meta$CB <- gsub(paste0(sample_name, "#"), "", rownames(subset_arch_clono_meta))
subset_arch_clono_meta$CLONOTYPE <- unlist(lapply(subset_arch_clono_meta$CB, function(x) svz_sub@meta.data[x, "CLONOTYPE"]))

subset_arch_clono@cellColData <- subset_arch_clono_meta

##visualize clonotypes in UMAP space:
pdf(paste0("ARCHR_", sample_name, "_CLONOTYPE_UMAPS.pdf"), width = 10, height = 10)
raster_set <- list()
for (CLONOTYPE in unique(subset_arch_clono_meta$CLONOTYPE)) {
    curr_plot <- plotEmbedding(ArchRProj = subset_arch[rownames(subset_arch_clono_meta[subset_arch_clono_meta$CLONOTYPE == CLONOTYPE,]),], colorBy = "cellColData", name = "LABEL_TRANSFER", embedding = "UMAP") + ggtitle(paste0("CLONOTYPE - ", CLONOTYPE))
    print(curr_plot)
        raster_set[[CLONOTYPE]] <- curr_plot
print(var)
}
dev.off()

library(gridExtra)
pdf(paste0("ARCHR_", sample_name, "_CLONOTYPE_UMAPS_AGGREGATED.pdf"), width = 20, height = 20)
grid.arrange(grobs = raster_set, nrow  = 4)
dev.off()

##Plot a single variant of interest too...
example_plot <- plotEmbedding(ArchRProj = subset_arch, colorBy = "cellColData", name = "M16012GA", embedding = "UMAP") + ggtitle("M16012GA")
pdf(paste0("ARCHR_", sample_name, "_example_mitoUMAP.pdf"), width = 10, height = 10)
print(example_plot)
dev.off()
