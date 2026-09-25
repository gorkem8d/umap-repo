# region_set_umap.R
# Reproduce the ATAC-only UMAP + clustering of Zamboni et al. (Fig. 1d/e) on the
# lab microglia SHARE-seq data, using Alena's lifted-over (hg -> mm10) region set
# as the feature space, then check where the original clusters end up.
#
# Three embeddings, identical processing, following the repo code for Fig. 1
# (LlorensLab/zamboni_et_al, code/preprocessing/merge_sample_objects.R), which
# differs from the Methods text in places:
#   FindTopFeatures(min.cutoff = "q75") -> RunTFIDF -> RunSVD   (Methods says min.cutoff = 10)
#   FindNeighbors on LSI 2:30 with ONE graph name -> only the kNN graph is stored,
#   FindClusters(algorithm = 3 (SLM), resolution = 1.1) on that kNN graph
#   No batch correction/integration: samples are merged and embedded as is.
# ATAC-only (LSI) on purpose, not WNN: the RNA half of a WNN is unchanged by the region
# choice and would carry the original (WNN-defined) clusters by itself. The original WNN
# clusters are used only as labels.
#   all  - all consensus peaks            -> our version of Fig. 1d (baseline)
#   reg  - Alena's regions                -> the actual question
#   rand - random consensus peaks, same n -> control: how much does ANY subset recover?
#
# UNTESTED against the real object - check the config block first.

suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(qs)
  library(GenomicRanges); library(rtracklayer)
  library(ggplot2); library(patchwork); library(mclust); library(Matrix)
})
set.seed(1234)

## ---- config: check these against the object ----
obj_path   <- "SHARE_mouse_data/MGobject_processed.qs"
bed_path   <- "regions_mm10_clean.bed"   # output of prep_regions.py, not the raw liftover BED
atac_assay <- "ATAC"                   # names(obj@assays)
orig_col   <- "seurat_clusters"        # column holding the original WNN clusters (MG1..MG12)
cond_col   <- "condition_memory_main"  # condition column, for colouring
sample_col <- "orig.ident"             # sample / mouse column, to spot batch-driven structure
top_cutoff <- "q75"                    # repo value; the Methods text says 10
dims       <- 2:30
res        <- 1.1
algorithm  <- 3                        # SLM, as in the repo
n_test     <- NULL                     # e.g. 2000 for the dry run; NULL = all cells
out_dir    <- "results/region_umap"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- qread(obj_path)
if (!is.null(n_test)) {
  obj <- subset(obj, cells = sample(colnames(obj), min(n_test, ncol(obj))))
  out_dir <- paste0(out_dir, "_test", n_test)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
obj$orig_cluster <- obj[[orig_col, drop = TRUE]]   # save before FindClusters overwrites seurat_clusters
DefaultAssay(obj) <- atac_assay

## ---- 1. load and sanity-check the region set ----
regions <- unique(import(bed_path, format = "BED"))
regions <- keepStandardChromosomes(regions, pruning.mode = "coarse")
peaks   <- granges(obj[[atac_assay]])
stopifnot("chromosome names do not match (chr1 vs 1? mm10 vs GRCm39?)" =
            any(seqlevels(regions) %in% seqlevels(peaks)))
n_ovl <- sum(overlapsAny(regions, peaks))
message(sprintf("%d regions; %d (%.1f%%) overlap a called peak; median width %d bp",
                length(regions), n_ovl, 100 * n_ovl / length(regions),
                as.integer(median(width(regions)))))

## ---- 2. count fragments directly over the regions ----
# Lifted regions rarely match called peaks exactly, so count over them directly.
# If this fails on file paths, the fragment paths stored in the object need UpdatePath().
frags  <- Fragments(obj[[atac_assay]])
counts <- FeatureMatrix(fragments = frags, features = regions, cells = colnames(obj))
counts <- counts[rowSums(counts > 0) > 0, ]
obj[["REG"]] <- CreateChromatinAssay(counts = counts)

## random control: same number of features, drawn from the consensus peaks.
## Autosomes only, because the region set has almost no chrX (source table lacked sex chromosomes).
all_counts <- GetAssayData(obj[[atac_assay]], slot = "counts")
auto_idx   <- which(as.character(seqnames(peaks)) %in% paste0("chr", 1:19))
set.seed(1234)   # reset so the random draw does not depend on the dry-run subsampling
rand_idx   <- sample(auto_idx, min(nrow(counts), length(auto_idx)))
obj[["RAND"]] <- CreateChromatinAssay(counts = all_counts[rand_idx, ], ranges = peaks[rand_idx])

## per-cell depth for each feature set (DepthCor reads nCount_<assay>)
obj$nCount_REG  <- colSums(GetAssayData(obj[["REG"]],  slot = "counts"))
obj$nCount_RAND <- colSums(GetAssayData(obj[["RAND"]], slot = "counts"))
if (!paste0("nCount_", atac_assay) %in% colnames(obj@meta.data))
  obj[[paste0("nCount_", atac_assay)]] <- colSums(all_counts)

## ---- 3. the paper's ATAC processing ----
run_atac <- function(obj, assay, prefix) {
  DefaultAssay(obj) <- assay
  obj <- FindTopFeatures(obj, min.cutoff = top_cutoff, verbose = FALSE)
  obj <- RunTFIDF(obj, verbose = FALSE)
  message(prefix, ": ", length(VariableFeatures(obj)), " features used for LSI")
  lsi <- paste0(prefix, "_lsi")
  obj <- RunSVD(obj, reduction.name = lsi, reduction.key = paste0(toupper(prefix), "LSI_"),
                verbose = FALSE)
  # reduction given explicitly: in the repo's Fig. 1 block RunUMAP has no
  # reduction argument, so Seurat's default ("pca", i.e. RNA) is used there
  obj <- RunUMAP(obj, reduction = lsi, dims = dims, verbose = FALSE,
                 reduction.name = paste0(prefix, "_umap"),
                 reduction.key  = paste0(toupper(prefix), "UMAP_"))
  obj <- FindNeighbors(obj, reduction = lsi, dims = dims, verbose = FALSE,
                       graph.name = paste0(prefix, c("_nn", "_snn")))
  # repo clusters on the kNN graph (single graph.name), not the SNN graph
  obj <- FindClusters(obj, graph.name = paste0(prefix, "_nn"), resolution = res,
                      algorithm = algorithm, verbose = FALSE)
  obj[[paste0(prefix, "_clusters")]] <- Idents(obj)
  obj
}
for (p in list(c("ATAC_all", "all"), c("REG", "reg"), c("RAND", "rand"))) {
  a <- if (p[1] == "ATAC_all") atac_assay else p[1]
  obj <- run_atac(obj, a, p[2])
}

## LSI_1 should track depth (that is why the paper drops it) - check, per embedding
pdf(file.path(out_dir, "depth_correlation.pdf"), width = 10, height = 3.5)
print(DepthCor(obj, assay = atac_assay, reduction = "all_lsi")  + ggtitle("all peaks") |
      DepthCor(obj, assay = "REG",      reduction = "reg_lsi")  + ggtitle("region set") |
      DepthCor(obj, assay = "RAND",     reduction = "rand_lsi") + ggtitle("random"))
dev.off()

## ---- 4. how well are the original clusters preserved? ----
knn_purity <- function(obj, graph, labels) {        # fraction of kNN sharing the cell's original label
  s <- summary(as(obj@graphs[[graph]], "dgCMatrix"))
  s <- s[s$i != s$j, ]
  mean(tapply(labels[s$i] == labels[s$j], s$i, mean))
}
lab <- as.character(obj$orig_cluster)
purity_chance <- sum(prop.table(table(lab))^2)   # expected kNN purity if labels were random
depth_col <- c(all = paste0("nCount_", atac_assay), reg = "nCount_REG", rand = "nCount_RAND")
metrics <- do.call(rbind, lapply(c("all", "reg", "rand"), function(p) data.frame(
  embedding    = p,
  n_features   = nrow(obj[[switch(p, all = atac_assay, reg = "REG", rand = "RAND")]]),
  n_lsi_feat   = length(VariableFeatures(obj[[switch(p, all = atac_assay, reg = "REG", rand = "RAND")]])),
  n_clusters   = nlevels(obj[[paste0(p, "_clusters"), drop = TRUE]]),
  ARI_vs_orig  = adjustedRandIndex(lab, obj[[paste0(p, "_clusters"), drop = TRUE]]),
  knn_purity   = knn_purity(obj, paste0(p, "_nn"), lab),
  purity_chance = purity_chance,
  median_counts_per_cell = median(obj[[depth_col[[p]], drop = TRUE]]))))
print(metrics)
write.csv(metrics, file.path(out_dir, "metrics.csv"), row.names = FALSE)

## Fig. 1e analogue: per original cluster, proportion of cells in each new cluster
agree_plot <- function(p) {
  tab  <- table(orig = lab, new = obj[[paste0(p, "_clusters"), drop = TRUE]])
  prop <- as.data.frame(sweep(tab, 1, rowSums(tab), "/"))
  ggplot(prop, aes(new, orig, fill = Freq)) + geom_tile() +
    scale_fill_viridis_c(limits = c(0, 1), name = "proportion") +
    labs(title = p, x = paste(p, "clusters"), y = "original clusters") + theme_minimal()
}

## ---- 5. plots ----
pdf(file.path(out_dir, "umaps.pdf"), width = 16, height = 11)
for (col in intersect(c("orig_cluster", cond_col, sample_col), colnames(obj@meta.data))) {
  print(wrap_plots(lapply(c("all", "reg", "rand"), function(p)
    DimPlot(obj, reduction = paste0(p, "_umap"), group.by = col, label = col == "orig_cluster",
            shuffle = TRUE) + ggtitle(paste(p, "-", col))), ncol = 3))
}
print(wrap_plots(lapply(c("all", "reg", "rand"), function(p)
  DimPlot(obj, reduction = paste0(p, "_umap"), group.by = paste0(p, "_clusters"), label = TRUE) +
    ggtitle(paste(p, "- own clusters"))), ncol = 3))
print(wrap_plots(lapply(c("all", "reg", "rand"), agree_plot), ncol = 3))
dev.off()

## ---- 6. export for Python ----
emb <- do.call(cbind, lapply(c("all", "reg", "rand"), function(p)
  Embeddings(obj, paste0(p, "_umap"))))
meta <- obj@meta.data[, intersect(c("orig_cluster", cond_col, sample_col,
                                    "all_clusters", "reg_clusters", "rand_clusters"),
                                  colnames(obj@meta.data))]
write.csv(cbind(barcode = colnames(obj), meta, emb), file.path(out_dir, "umap_coords_clusters.csv"),
          row.names = FALSE)
