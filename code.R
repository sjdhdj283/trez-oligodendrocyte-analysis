# Code #### 
rm(list = ls())
# analysis_utils ####
options(Seurat.warn.umap.uwot = FALSE)

pdf_device_for_output <- function() {
  if (isTRUE(capabilities("cairo"))) {
    return(grDevices::cairo_pdf)
  }
  grDevices::pdf
}

safe_ggsave <- function(filename, plot, ...) {
  args <- list(...)
  if (grepl("\\.pdf$", filename, ignore.case = TRUE) && is.null(args$device)) {
    args$device <- pdf_device_for_output()
  }
  ok <- tryCatch({
    do.call(ggplot2::ggsave,
            c(list(filename = filename, plot = plot), args))
    TRUE
  }, error = function(e) {
    alt <- sub("\\.pdf$", paste0("_", format(Sys.time(), "%H%M%S"), ".pdf"), filename)
    message("Cannot write ", basename(filename), "; saving as ", basename(alt))
    do.call(ggplot2::ggsave,
            c(list(filename = alt, plot = plot), args))
    FALSE
  })
  invisible(ok)
}

safe_pdf <- function(filename, width, height, ...) {
  open_device <- function(path) {
    pdf_device_for_output()(path, width = width, height = height, ...)
  }
  tryCatch({
    open_device(filename)
    TRUE
  }, error = function(e) {
    alt <- sub("\\.pdf$", paste0("_", format(Sys.time(), "%H%M%S"), ".pdf"), filename)
    message("Cannot write ", basename(filename), "; saving as ", basename(alt))
    open_device(alt)
    FALSE
  })
}

ucell_bpparam <- function() {
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    return(BiocParallel::SerialParam())
  }
  NULL
}

sanitize_feature_names <- function(mat) {
  rn <- rownames(mat)
  if (!is.null(rn)) rownames(mat) <- gsub("_", "-", rn, fixed = TRUE)
  mat
}

drop_scale_data_layer <- function(object, assay = "RNA") {
  if (!assay %in% names(object@assays)) return(object)
  assay_obj <- object[[assay]]
  if (inherits(assay_obj, "Assay5")) {
    assay_obj$scale.data <- NULL
    object[[assay]] <- assay_obj
  } else if ("scale.data" %in% slotNames(assay_obj)) {
    assay_obj@scale.data <- matrix(numeric(0), nrow = 0, ncol = 0)
    object[[assay]] <- assay_obj
  }
  object
}

floor_zero_pvalues <- function(x) {
  x <- as.numeric(x)
  pos <- x[is.finite(x) & x > 0]
  if (length(pos) == 0) return(x)
  x[is.finite(x) & x <= 0] <- min(pos, na.rm = TRUE) * 0.1
  x
}

quiet_known_package_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      known <- grepl("Dropping illegal slot\\(s\\) norm.method", msg) ||
        grepl("Using `size` aesthetic", msg, fixed = TRUE) ||
        grepl("The `size` argument of `element_", msg, fixed = TRUE) ||
        grepl("element_rect\\(\\)", msg) ||
        grepl("One or more p-values is 0", msg, fixed = TRUE)
      if (known) invokeRestart("muffleWarning")
    }
  )
}

fetch_seurat_data <- function(object, vars, layer = "data", assay = NULL) {
  args <- list(object = object, vars = vars)
  if (!is.null(assay)) args$assay <- assay
  args$layer <- layer
  out <- tryCatch(do.call(Seurat::FetchData, args), error = function(e) NULL)
  if (!is.null(out)) return(out)
  args$layer <- NULL
  args$slot <- layer
  do.call(Seurat::FetchData, args)
}

safe_wilcox_p <- function(x, group, alternative = "two.sided") {
  keep <- is.finite(x) & !is.na(group)
  x <- x[keep]
  group <- as.factor(group[keep])
  if (length(levels(group)) != 2 || any(table(group) < 1)) return(NA_real_)
  tryCatch(
    stats::wilcox.test(x ~ group, alternative = alternative, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

safe_wilcox_two_sample <- function(x1, x2, alternative = "two.sided") {
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  if (length(x1) < 1 || length(x2) < 1) return(NA_real_)
  tryCatch(
    stats::wilcox.test(x1, x2, alternative = alternative, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

format_p_label <- function(p) {
  if (length(p) == 0 || !is.finite(p)) return("p=NA")
  if (p < 1e-4) return("p<0.0001")
  paste0("p=", sprintf("%.4f", p))
}

resolve_feature_names <- function(features, available_features) {
  features <- unique(as.character(features))
  available_features <- as.character(available_features)
  available_lower <- tolower(available_features)
  
  resolved <- rep(NA_character_, length(features))
  for (i in seq_along(features)) {
    exact_idx <- match(features[i], available_features)
    if (!is.na(exact_idx)) {
      resolved[i] <- available_features[exact_idx]
      next
    }
    case_idx <- match(tolower(features[i]), available_lower)
    if (!is.na(case_idx)) resolved[i] <- available_features[case_idx]
  }
  
  mapping <- data.frame(
    requested = features,
    resolved = resolved,
    stringsAsFactors = FALSE
  )
  out <- unique(stats::na.omit(resolved))
  attr(out, "feature_mapping") <- mapping
  attr(out, "missing_features") <- mapping$requested[is.na(mapping$resolved)]
  out
}

make_hc_tn_gene_violin <- function(object, features, group_col = "group",
                                   group_levels = c("HC", "TN"),
                                   assay = "RNA", layer = "data",
                                   cols = NULL, ncol = 3,
                                   title = NULL, subtitle = NULL,
                                   y_label = "RNA log-normalized expression") {
  if (!assay %in% names(object@assays)) {
    stop("Assay not found: ", assay, call. = FALSE)
  }
  features_resolved <- resolve_feature_names(features, rownames(object[[assay]]))
  if (length(features_resolved) == 0) {
    stop("None of the requested features were found in assay ", assay, call. = FALSE)
  }
  
  expr_df <- fetch_seurat_data(object, vars = features_resolved,
                               layer = layer, assay = assay)
  meta <- object@meta.data[rownames(expr_df), , drop = FALSE]
  if (!group_col %in% colnames(meta)) {
    stop("Metadata column not found: ", group_col, call. = FALSE)
  }
  
  group <- factor(as.character(meta[[group_col]]), levels = group_levels)
  keep <- !is.na(group)
  expr_mat <- as.matrix(expr_df[keep, features_resolved, drop = FALSE])
  storage.mode(expr_mat) <- "numeric"
  group <- group[keep]
  
  plot_df <- data.frame(
    gene = factor(rep(features_resolved, each = nrow(expr_mat)),
                  levels = features_resolved),
    expression = as.vector(expr_mat),
    group = rep(group, times = length(features_resolved)),
    stringsAsFactors = FALSE
  )
  plot_df$group <- factor(plot_df$group, levels = group_levels)
  plot_df$group_x <- as.numeric(plot_df$group)
  
  annot_df <- do.call(rbind, lapply(features_resolved, function(gene) {
    sub <- plot_df[plot_df$gene == gene & is.finite(plot_df$expression), ,
                   drop = FALSE]
    y_range <- if (nrow(sub) > 0) range(sub$expression, na.rm = TRUE) else c(0, 0)
    y_span <- diff(y_range)
    if (!is.finite(y_span) || y_span == 0) y_span <- 0.1
    p_val <- safe_wilcox_p(sub$expression, sub$group)
    data.frame(
      gene = factor(gene, levels = features_resolved),
      x = 1.08,
      xend = 1.92,
      xmid = 1.5,
      line_y = y_range[2] + 0.22 * y_span,
      label_y = y_range[2] + 0.32 * y_span,
      panel_top = y_range[2] + 0.50 * y_span,
      p_value = p_val,
      p_label = format_p_label(p_val),
      stringsAsFactors = FALSE
    )
  }))
  
  if (is.null(cols)) cols <- stats::setNames(c("#C97A74", "#5F7FA8"), group_levels)
  if (!is.null(names(cols))) {
    cols <- cols[group_levels]
    cols[is.na(cols)] <- "#B8B2AC"
  }
  
  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = group_x, y = expression, fill = group)) +
    ggplot2::geom_violin(
      ggplot2::aes(group = group),
      trim = FALSE, scale = "width", alpha = 0.78,
      color = "grey25", linewidth = 0.35, width = 0.72,
      na.rm = TRUE
    ) +
    ggplot2::geom_blank(
      data = annot_df,
      ggplot2::aes(x = xmid, y = panel_top),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      data = annot_df,
      ggplot2::aes(x = x, xend = xend, y = line_y, yend = line_y),
      inherit.aes = FALSE, linewidth = 0.45
    ) +
    ggplot2::geom_text(
      data = annot_df,
      ggplot2::aes(x = xmid, y = label_y, label = p_label),
      inherit.aes = FALSE, size = 3.1, fontface = "bold"
    ) +
    ggplot2::facet_wrap(stats::as.formula("~gene"), scales = "free_y", ncol = ncol) +
    ggplot2::scale_fill_manual(values = cols, breaks = group_levels,
                               drop = FALSE, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq_along(group_levels),
                                labels = group_levels,
                                expand = ggplot2::expansion(add = 0.35)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.03))) +
    ggplot2::labs(x = NULL, y = y_label, title = title, subtitle = subtitle) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "top",
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.major.x = ggplot2::element_blank()
    )
  
  attr(p, "features") <- features_resolved
  attr(p, "missing_features") <- attr(features_resolved, "missing_features")
  attr(p, "feature_mapping") <- attr(features_resolved, "feature_mapping")
  attr(p, "p_values") <- annot_df[, c("gene", "p_value", "p_label")]
  p
}

check_groups_present <- function(meta, group_col = "group", required = c("HC", "TN")) {
  missing_groups <- setdiff(required, unique(as.character(meta[[group_col]])))
  if (length(missing_groups) > 0) {
    stop("Missing required group(s): ", paste(missing_groups, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

scale_matrix_safe <- function(mat) {
  mat <- as.matrix(mat)
  mu <- colMeans(mat, na.rm = TRUE)
  sdv <- apply(mat, 2, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  sweep(sweep(mat, 2, mu, "-"), 2, sdv, "/")
}

assign_dominant_stage <- function(meta, stage_order, margin_cutoff = 0.05, use_z = TRUE) {
  stage_mat <- as.matrix(meta[, stage_order, drop = FALSE])
  score_mat <- if (use_z) scale_matrix_safe(stage_mat) else stage_mat
  dominant <- apply(score_mat, 1, function(x) stage_order[which.max(x)])
  margin <- apply(score_mat, 1, function(x) {
    sx <- sort(x, decreasing = TRUE)
    if (length(sx) < 2) return(NA_real_)
    sx[1] - sx[2]
  })
  list(
    dominant_stage = factor(dominant, levels = stage_order),
    stage_margin = margin,
    stage_ambiguous = margin < margin_cutoff,
    score_matrix = score_mat
  )
}

complete_stage_proportions <- function(meta, stage_col = "dominant_stage", stage_order) {
  sample_group <- dplyr::distinct(meta[, c("sample", "group")])
  grid <- merge(sample_group,
                data.frame(dominant_stage = factor(stage_order, levels = stage_order)),
                all = TRUE)
  counts <- meta %>%
    dplyr::count(sample, group, !!rlang::sym(stage_col), name = "n") %>%
    dplyr::rename(dominant_stage = !!rlang::sym(stage_col))
  grid %>%
    dplyr::left_join(counts, by = c("sample", "group", "dominant_stage")) %>%
    dplyr::mutate(n = dplyr::coalesce(n, 0L)) %>%
    dplyr::group_by(sample, group) %>%
    dplyr::mutate(pct = 100 * n / sum(n)) %>%
    dplyr::ungroup()
}

sig_star <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05,  "*", "ns"))))
}

select_volcano_labels <- function(data, gene_col = "gene", p_col, fc_col,
                                  candidate_genes = NULL, max_labels = 12,
                                  max_per_direction = ceiling(max_labels / 2)) {
  needed <- c(gene_col, p_col, fc_col)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing volcano label column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  
  df <- data.frame(
    gene = as.character(data[[gene_col]]),
    p = as.numeric(data[[p_col]]),
    fc = as.numeric(data[[fc_col]]),
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$gene) & nzchar(df$gene) &
             is.finite(df$p) & is.finite(df$fc), , drop = FALSE]
  if (!is.null(candidate_genes)) {
    candidate_genes <- unique(as.character(candidate_genes))
    candidate_df <- df[df$gene %in% candidate_genes, , drop = FALSE]
    if (nrow(candidate_df) > 0) df <- candidate_df
  }
  if (nrow(df) == 0) return(character())
  
  df$p_rank <- pmax(df$p, .Machine$double.xmin)
  df$abs_fc <- abs(df$fc)
  pick_top <- function(x, n) {
    if (nrow(x) == 0 || n <= 0) return(x[0, , drop = FALSE])
    x[order(x$p_rank, -x$abs_fc, x$gene), , drop = FALSE][
      seq_len(min(n, nrow(x))), , drop = FALSE]
  }
  
  picked <- rbind(
    pick_top(df[df$fc >= 0, , drop = FALSE], max_per_direction),
    pick_top(df[df$fc < 0, , drop = FALSE], max_per_direction)
  )
  if (nrow(picked) < max_labels) {
    rest <- df[!df$gene %in% picked$gene, , drop = FALSE]
    picked <- rbind(picked, pick_top(rest, max_labels - nrow(picked)))
  }
  picked <- picked[order(picked$p_rank, -picked$abs_fc, picked$gene), , drop = FALSE]
  unique(head(picked$gene, max_labels))
}

# plot_palette ####
COL_CONTROL <- "#5F7FA8"
COL_DISEASE <- "#C97A74"
COL_NEUTRAL_DARK  <- "#7E7974"
COL_NEUTRAL_MID   <- "#B8B2AC"
COL_NEUTRAL_LIGHT <- "#F2EEE8"
COL_NEUTRAL_PALE  <- "#E3DED8"
COL_GROUP <- c(
  HC  = COL_CONTROL,
  TN  = COL_DISEASE,
  Con = COL_CONTROL,
  CSD = COL_DISEASE
)
COL_SIGN <- c(
  "FALSE" = COL_DISEASE,
  "TRUE"  = COL_CONTROL
)
COL_REGULATION <- c(
  Down = COL_DISEASE,
  ns   = COL_NEUTRAL_MID,
  Up   = COL_CONTROL
)
COL_DIVERGING <- c(
  low  = COL_DISEASE,
  mid  = COL_NEUTRAL_LIGHT,
  high = COL_CONTROL
)
COL_STAGE <- c(
  OPC  = "#6B89AE",
  COP  = "#93A7BE",
  NFOL  = "#CEC8C1",
  MFOL = "#D8A19A",
  MOL  = "#C17970"
)
COL_STAGE_DOMINANT <- c(
  OPC_dominant  = unname(COL_STAGE["OPC"]),
  COP_dominant  = unname(COL_STAGE["COP"]),
  NFOL_dominant  = unname(COL_STAGE["NFOL"]),
  MFOL_dominant = unname(COL_STAGE["MFOL"]),
  MOL_dominant  = unname(COL_STAGE["MOL"]),
  OL_dominant   = unname(COL_STAGE["MOL"])
)
COL_EARLY_LATE <- c(
  "Early (OPC+COP)"     = "#7D98B8",
  "Intermediate (NFO)"  = unname(COL_STAGE["NFO"]),
  "Late (MFOL+MOL)"     = COL_CONTROL
)
COL_CLUSTER_LABEL <- c(
  C3_HealthyMyelin = "#90A6BF",
  C5_RootZone      = "#8E8983",
  C7_TNMyelin      = COL_CONTROL,
  other            = COL_NEUTRAL_PALE
)
COL_CATEGORY <- c(
  OPC_lineage   = COL_DISEASE,
  NFO_premyelin = "#9AA9BA",
  Mature_myelin = COL_CONTROL,
  Neuronal      = COL_NEUTRAL_DARK
)
diverging_palette <- function(n = 100) {
  colorRampPalette(unname(COL_DIVERGING))(n)
}
seq_blue_palette <- function(n = 100) {
  colorRampPalette(c(COL_NEUTRAL_LIGHT, "#D9E2EC", "#AFC0D3", COL_DISEASE))(n)
}

seq_red_palette <- function(n = 100) {
  colorRampPalette(c(COL_NEUTRAL_LIGHT, "#E6D6D1", "#D5A6A0", COL_CONTROL))(n)
}

# cluster_TREZ ####
base <- "."
mat_dir <- file.path(base, "Matrix")
out_dir <- file.path(base, "results_paper_TREZ")
for (sub in c("A_clustering", "B_mapping", "C_TREZ_analysis")) {
  dir.create(file.path(out_dir, sub), showWarnings = FALSE, recursive = TRUE)}
set.seed(42)
options(future.globals.maxSize = 64 * 1024^3)

required <- c(
  "Seurat", 
  "MAST", 
  "UCell", 
  "clusterProfiler",
  "org.Rn.eg.db", 
  "enrichplot", 
  "EnhancedVolcano",
  "dplyr", 
  "tidyr", 
  "ggplot2", 
  "patchwork", 
  "ggpubr")
missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  cat("The following packages are missing; please install them first:\n"); print(missing)
  cat("\nInstallation command:\n")
  cat('if (!require("BiocManager")) install.packages("BiocManager")\n')
  cat('BiocManager::install(c(', paste0('"', missing, '"', collapse = ", "),'))\n', sep = "") 
  stop("Install missing dependencies first, then continue.")}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
  library(UCell)
  library(MAST)
  library(clusterProfiler)
  library(org.Rn.eg.db)
  library(enrichplot)
  library(EnhancedVolcano)})

load_sample <- function(name) {
  path <- file.path(mat_dir, name, "Gene/filtered")
  cat("load:", name, "from", path, "\n")
  m <- ReadMtx(
    mtx = file.path(path, "matrix.mtx"),
    cells = file.path(path, "barcodes.tsv"),
    features = file.path(path, "features.tsv"))
  m <- sanitize_feature_names(m)
  s <- CreateSeuratObject(m, project = name, min.cells = 3, min.features = 200)
  s$sample <- name
  s$group  <- ifelse(name == "HC", "HC", "TN")
  s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^Mt-")
  s}
samples <- lapply(c("HC", "TN1", "TN2", "TN3"), load_sample)
obj <- merge(samples[[1]], samples[-1], add.cell.ids = c("HC", "TN1", "TN2", "TN3"))
cat("Number of spots after merging:", ncol(obj), "\n")
cat("Group distribution:\n"); print(table(obj$sample, obj$group))

obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 20)
cat("Number of QC spot counts:", ncol(obj), "\n\n")

obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = TRUE)
obj <- RunPCA(obj, npcs = 30, verbose = TRUE)
obj <- RunUMAP(obj, dims = 1:30, verbose = TRUE, umap.method = "uwot", metric = "cosine")
obj <- FindNeighbors(obj, dims = 1:30, verbose = TRUE)
res_grid <- c(0.15, 0.18, 0.20, 0.22, 0.25, 0.28, 0.30)
res_log <- sapply(res_grid, function(r) {
  tmp <- FindClusters(obj, resolution = r, verbose = FALSE)
  length(unique(Idents(tmp)))})
res_df <- data.frame(resolution = res_grid, n_clusters = res_log)
print(res_df)
write.csv(res_df, file.path(out_dir, "A_clustering", "resolution_sweep.csv"), row.names = FALSE)

pca_emb <- Embeddings(obj, "pca")[, 1:30]
set.seed(42)
km <- kmeans(pca_emb, centers = 10, nstart = 25, iter.max = 500, algorithm = "MacQueen")
obj$seurat_clusters <- factor(km$cluster, levels = 1:10)
Idents(obj) <- "seurat_clusters"
cat("K-means convergence state:", km$ifault, "(0 = Normal convergence)\n")
cat("Final number of clusters:", length(unique(Idents(obj))), "\n")
cat("Number of spots in each cluster:\n"); print(table(Idents(obj)))
saveRDS(obj, file.path(out_dir, "A_clustering", "obj_10clusters.rds"))

myelin_mk <- c("Mobp", "Sema4d", "Mog", "Mlc1", "Agt", "Thy1", "Mt2A", "S100a1", "Atp1a3", "Gja1")
root_mk <- c("S100a8", "S100a9", "Alox15", "Slc4a1", "Ibsp", "RatNP-3b", "Lgals5", "Ahsp")
myelin_mk_present <- intersect(myelin_mk, rownames(obj))
root_mk_present   <- intersect(root_mk,   rownames(obj))
cat("Original myelin gene (", length(myelin_mk_present), "/", length(myelin_mk), "): ", paste(myelin_mk_present, collapse = ", "), "\n", sep = "")
cat("Original root gene (", length(root_mk_present), "/", length(root_mk), "): ", paste(root_mk_present, collapse = ", "), "\n\n", sep = "")

p_dot <- DotPlot(obj, features = c(myelin_mk_present, root_mk_present)) +
  RotatedAxis() +
  ggtitle("Paper markers — identify C3 (HC-myelin), C5 (root zone), C7 (TN-myelin)")
safe_ggsave(file.path(out_dir, "B_mapping", "dotplot_paper_markers.pdf"), p_dot, width = 10, height = 6)

comp <- as.data.frame.matrix(table(Cluster = Idents(obj), Group = obj$group))
comp$total  <- rowSums(comp)
comp$pct_HC <- round(comp$HC / comp$total * 100, 1)
comp$pct_TN <- round(comp$TN / comp$total * 100, 1)
comp$Cluster <- rownames(comp)
comp <- comp[, c("Cluster", "HC", "TN", "total", "pct_HC", "pct_TN")]
write.csv(comp, file.path(out_dir, "B_mapping", "cluster_composition.csv"), row.names = FALSE)
cat("Composition of each cluster HC vs TN:\n"); print(comp)
paper_C3 <- "3"
paper_C5 <- "9"
paper_C7 <- "8"

if (any(c(paper_C3, paper_C5, paper_C7) == c("X", "Y", "Z"))) {
  cat("Pause: Please open the following two files and fill in paper_C3 / paper_C5 / paper_C7:\n")
  cat(" 1)", file.path(out_dir, "B_mapping", "dotplot_paper_markers.pdf"), "\n")
  cat(" 2)", file.path(out_dir, "B_mapping", "cluster_composition.csv"), "\n")
  cat("Then rerun the script (or continue from Part B).\n")
  stop("You need to manually map the cluster ID before continuing.", call. = FALSE)}
stopifnot(all(c(paper_C3, paper_C5, paper_C7) %in% as.character(levels(Idents(obj)))))
map_df <- data.frame(
  paper_label = c("C3_HealthyMyelin", "C5_RootZone", "C7_TNMyelin"),
  seurat_id   = c(paper_C3, paper_C5, paper_C7))
write.csv(map_df, file.path(out_dir, "B_mapping", "cluster_label_map.csv"), row.names = FALSE)
cat("\n  Cluster mapping:\n"); print(map_df)
label_map <- setNames(map_df$paper_label, map_df$seurat_id)
obj$paper_label <- unname(label_map[as.character(Idents(obj))])
obj$paper_label[is.na(obj$paper_label)] <- "other"
print(table(obj$paper_label))

obj_trez <- subset(obj, subset = paper_label %in% 
                     c("C3_HealthyMyelin", "C5_RootZone", "C7_TNMyelin"))
cat("Total TREZ spots:", ncol(obj_trez), "\n")

counts <- obj_trez@meta.data %>%
  dplyr::group_by(sample, group, paper_label) %>%
  dplyr::summarise(n = n(), .groups = "drop")
write.csv(counts, file.path(out_dir, "C_TREZ_analysis", "TREZ_spot_counts.csv"), row.names = FALSE)
cat("Sample distribution:\n"); print(counts)
min_by_sample <- obj_trez@meta.data %>% 
  dplyr::group_by(sample) %>% summarise(n = n())
if (any(min_by_sample$n < 20)) {warning("The sample has fewer than 20 TREZ spots, so the statistics may be unstable.")}

DefaultAssay(obj_trez) <- "RNA"
obj_trez <- NormalizeData(obj_trez, verbose = FALSE)
obj_trez <- ScaleData(obj_trez, verbose = FALSE)

opc_strict <- intersect(c("Olig1", "Olig2", "Gpr17", "Neu4", "Pcdh15", "Sox6", "Nkx2-2", "Vcan", "Bcan"), rownames(obj_trez))
ol_strict <- intersect(c("Mog", "Cnp", "Enpp6", "Myrf", "Ugt8", "Fa2h", "Mag", "Ermn", "Aspa"), rownames(obj_trez))
opc_classic <- intersect(c("Pdgfra", "Cspg4"), rownames(obj_trez))
volcano_label_sets <- list(
  OPC_lineage = c(
    "Pdgfra", "Cspg4", "Olig1", "Olig2", "Gpr17", "Neu4",
    "Pcdh15", "Nkx2-2", "Sox6", "Sox10", "Vcan", "Bcan",
    "Bmp4", "Sox4", "Fyn"),
  NFO_premyelin = c(
    "Itpr2", "Enpp6", "Tcf7l2", "Tns3", "Gjb1", "Plp1",
    "Opalin", "Cnp", "Mal", "Ctps1"),
  Mature_myelin = c(
    "Mog", "Mag", "Mbp", "Mobp", "Klk6", "Ermn", "Myrf",
    "Ugt8", "Fa2h", "Aspa", "Apod", "Ptgds", "Hapln2",
    "Gjc2", "Mlc1", "Tf", "S100b", "Agpat4", "Cdkn1c",
    "Cldn11", "Pllp", "Ndrg1"))
volcano_marker_df <- data.frame(
  gene = unlist(volcano_label_sets, use.names = FALSE),
  category = rep(names(volcano_label_sets), lengths(volcano_label_sets)),
  stringsAsFactors = FALSE)
volcano_marker_df <- volcano_marker_df[!duplicated(volcano_marker_df$gene), ]
cat("OPC_strict (", length(opc_strict), "): ", paste(opc_strict, collapse = ", "), "\n", sep = "")
cat("OL (", length(ol_strict), "): ", paste(ol_strict, collapse = ", "), "\n", sep = "")
cat("OPC_classic (", length(opc_classic), "): ", paste(opc_classic, collapse = ", "), "\n", sep = "")

feat_list <- list(
  OPC_strict  = opc_strict,
  OL = ol_strict,
  OPC_classic = opc_classic)
feat_list <- feat_list[sapply(feat_list, length) > 0]
obj_trez <- AddModuleScore_UCell(
  obj_trez, features = feat_list,
  assay = "RNA", slot = "data", name = "",
  BPPARAM = ucell_bpparam())
obj_trez$OPC_OL_ratio <- obj_trez$OPC_strict / (obj_trez$OL + 0.01)

m <- obj_trez@meta.data
metrics <- c("OPC_strict", "OL", "OPC_OL_ratio")
if ("OPC_classic" %in% names(m)) metrics <- c(metrics, "OPC_classic")
stats_df <- data.frame(
  metric = metrics,
  HC_median = sapply(metrics, function(g) median(m[[g]][m$group == "HC"])),
  TN_median = sapply(metrics, function(g) median(m[[g]][m$group == "TN"])),
  p_value_spot_level = sapply(metrics, function(g) safe_wilcox_p(m[[g]], m$group)))
stats_df$log2FC <- log2((stats_df$TN_median + 0.01) / (stats_df$HC_median + 0.01))
rownames(stats_df) <- NULL
write.csv(stats_df, file.path(out_dir, "C_TREZ_analysis", "OPC_OL_UCell_stats.csv"), row.names = FALSE)
cat("\n Key statistics:\n"); print(stats_df)

cols_g <- COL_GROUP[c("HC", "TN")]
make_main_violin <- function(feature, plot_title) {
  plot_df <- data.frame(
    group = factor(m$group, levels = c("HC", "TN")),
    score = m[[feature]])
  plot_df$group_x <- as.numeric(plot_df$group)
  y_values <- plot_df$score
  y_values <- y_values[is.finite(y_values)]
  y_range <- if (length(y_values) > 0) range(y_values) else c(0, 1)
  y_span <- diff(y_range)
  if (!is.finite(y_span) || y_span == 0) y_span <- 0.1
  p_val <- stats_df$p_value_spot_level[match(feature, stats_df$metric)]
  line_y <- y_range[2] + 0.22 * y_span
  label_y <- y_range[2] + 0.32 * y_span
  panel_top <- y_range[2] + 0.48 * y_span
  panel_bottom <- y_range[1] - 0.03 * y_span
  p_df <- data.frame(
    x = 1.08,
    xend = 1.92,
    xmid = 1.5,
    line_y = line_y,
    label_y = label_y,
    p_label = format_p_label(p_val))
  ggplot(plot_df, aes(x = group_x, y = score, fill = group)) +
    geom_violin(
      aes(group = group),
      trim = TRUE, scale = "width",
      color = "grey30", linewidth = 0.35,
      width = 0.72, na.rm = TRUE) +
    geom_segment(
      data = p_df,
      aes(x = x, xend = xend, y = line_y, yend = line_y),
      inherit.aes = FALSE, linewidth = 0.45) +
    geom_text(
      data = p_df,
      aes(x = xmid, y = label_y, label = p_label),
      inherit.aes = FALSE, size = 3.4, fontface = "bold") + 
    scale_fill_manual(values = cols_g, guide = "none") +
    scale_x_continuous(breaks = c(1, 2), labels = c("HC", "TN"), expand = expansion(add = 0.35)) +
    labs(x = NULL, y = feature, title = plot_title) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank()) +
    coord_cartesian(ylim = c(panel_bottom, panel_top)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)))}
p1 <- make_main_violin("OPC_strict", "OPC score (UCell, strict)")
p2 <- make_main_violin("OL", "Mature OL score (UCell)")
p3 <- make_main_violin("OPC_OL_ratio", "OPC/OL score ratio")
safe_ggsave(file.path(out_dir, "C_TREZ_analysis", "violin_main.pdf"), p1 | p2 | p3, width = 12, height = 4.5)

p_bysample <- m %>% 
  dplyr::select(sample, group, OPC_strict, OL, OPC_OL_ratio) %>% 
  tidyr::pivot_longer(cols = c(OPC_strict, OL, OPC_OL_ratio), names_to = "metric") %>%
  ggplot(aes(x = sample, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.3) +
  scale_fill_manual(values = cols_g) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  ggtitle("TREZ per-sample breakdown (HC / TN1 / TN2 / TN3)")
safe_ggsave(file.path(out_dir, "C_TREZ_analysis", "boxplot_per_sample.pdf"), p_bysample, width = 12, height = 4)

key_gene_targets <- c("Olig1", "Olig2", "Gpr17", "Mog", "Cnp", "Enpp6", "Myrf")
key_genes <- resolve_feature_names(key_gene_targets, rownames(obj_trez[["RNA"]]))
if (length(attr(key_genes, "missing_features")) > 0) {
  cat("Key genes not found: ",paste(attr(key_genes, "missing_features"), collapse = ", "), "\n", sep = "")}
if (length(key_genes) > 0) {
  p_keys <- make_hc_tn_gene_violin(
    obj_trez,
    features = key_gene_targets,
    group_col = "group",
    assay = "RNA",
    layer = "data",
    cols = cols_g,
    ncol = 3,
    title = "Key OPC / OL genes in TREZ: HC vs TN",
    subtitle = "Short line marks the HC-TN comparison; P values are descriptive spot-level Wilcoxon tests")
  write.csv(attr(p_keys, "p_values"), file.path(out_dir, "C_TREZ_analysis", "violin_key_gene_pvalues.csv"), row.names = FALSE)
  safe_ggsave(file.path(out_dir, "C_TREZ_analysis", "violin_key_genes.pdf"), p_keys, width = 11, height = max(5.2, 2.45 * ceiling(length(key_genes) / 3)))}

all_markers <- c(opc_strict, ol_strict)
if (length(all_markers) > 0) {
  p_dot2 <- DotPlot(
    obj_trez, features = all_markers,
    group.by = "group", scale = FALSE,
    cols = c(COL_DIVERGING["low"], COL_DIVERGING["high"])) + 
    RotatedAxis() +
    ggtitle("OPC / OL markers in TREZ: HC vs TN")
  safe_ggsave(file.path(out_dir, "C_TREZ_analysis", "dotplot_OPC_OL.pdf"),p_dot2, width = 10, height = 5)}

DefaultAssay(obj_trez) <- "SCT"
obj_trez <- PrepSCTFindMarkers(obj_trez, verbose = FALSE)
Idents(obj_trez) <- "group"
deg <- suppressMessages(quiet_known_package_warnings(FindMarkers(
  obj_trez, ident.1 = "TN", ident.2 = "HC",
  test.use = "MAST", assay = "SCT",
  logfc.threshold = 0, min.pct = 0.05,
  verbose = FALSE)))
deg$gene <- rownames(deg)
write.csv(deg, file.path(out_dir, "C_TREZ_analysis", "DEG_MAST_all.csv"), row.names = FALSE)
deg_sig <- subset(deg, p_val < 0.05 & abs(avg_log2FC) >= log2(1.2))
write.csv(deg_sig, file.path(out_dir, "C_TREZ_analysis", "DEG_MAST_sig.csv"), row.names = FALSE)
cat("DEG total:", nrow(deg)," spot-level flagged genes (p_val<0.05 & |log2FC|>=log2(1.2)):", nrow(deg_sig), "\n")

deg_markers <- deg[deg$gene %in% c(opc_strict, ol_strict, opc_classic), ]
deg_markers$type <- case_when(
  deg_markers$gene %in% opc_strict ~ "OPC_strict",
  deg_markers$gene %in% ol_strict  ~ "OL",
  deg_markers$gene %in% opc_classic ~ "OPC_classic",
  TRUE ~ "other")
deg_markers <- deg_markers %>% arrange(type, desc(avg_log2FC))
write.csv(deg_markers, file.path(out_dir, "C_TREZ_analysis", "DEG_MAST_markers_only.csv"), row.names = FALSE)
print(deg_markers[, c("gene", "type", "avg_log2FC", "p_val_adj", "pct.1", "pct.2")])

deg$p_val_plot <- floor_zero_pvalues(deg$p_val)
volcano_labels <- deg_sig$gene[deg_sig$gene %in% volcano_marker_df$gene]
volcano_missing_labels <- setdiff(volcano_marker_df$gene, deg$gene)
volcano_labeled_df <- deg[match(volcano_labels, deg$gene), c("gene", "avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")]
volcano_labeled_df$category <- volcano_marker_df$category[match(volcano_labeled_df$gene, volcano_marker_df$gene)]
volcano_labeled_df$direction <- ifelse(volcano_labeled_df$avg_log2FC > 0, "TN_up", "TN_down")
volcano_labeled_df$is_differential <- TRUE
volcano_labeled_df$label_source <- "differential_pval_0.05_fc_1.2"
volcano_labeled_df <- volcano_labeled_df %>%
  dplyr::arrange(factor(category, levels = names(volcano_label_sets)), desc(avg_log2FC), gene)
write.csv(volcano_labeled_df, file.path(out_dir, "C_TREZ_analysis", "volcano_TREZ_labeled_genes.csv"), row.names = FALSE)
cat("Differential volcano labels (OPC/NFO/myelin genes, p_val<0.05 & |log2FC|>=log2(1.2)): ", length(volcano_labels), "\n", sep = "")
if (length(volcano_missing_labels) > 0) {cat("  Volcano label genes not in DEG table: ", paste(volcano_missing_labels, collapse = ", "), "\n", sep = "")}
p_vol <- quiet_known_package_warnings(
  EnhancedVolcano(
    deg, lab = deg$gene,
    selectLab = volcano_labels,
    x = "avg_log2FC", y = "p_val_plot",
    pCutoff = 0.05, FCcutoff = log2(1.2),
    title = "TREZ (C3+C5+C7): TN vs HC",
    subtitle = sprintf(
      "p < 0.05 and |log2FC| >= log2(1.2); OPC/NFO/myelin labels shown (n = %d)", length(volcano_labels)),
    drawConnectors = TRUE, max.overlaps = Inf,
    labSize = 2.6, pointSize = 1.25, colAlpha = 0.5,
    col = c(COL_NEUTRAL_PALE, COL_NEUTRAL_MID, COL_NEUTRAL_DARK, COL_CONTROL)))
safe_ggsave(file.path(out_dir, "C_TREZ_analysis", "volcano_TREZ.pdf"), p_vol, width = 11, height = 8.5)

up_genes <- deg_sig$gene[deg_sig$avg_log2FC > 0]
dn_genes <- deg_sig$gene[deg_sig$avg_log2FC < 0]
enrich_plots <- list()
enrich_results <- list()
for (direction in c("up", "down")) {
  genes <- if (direction == "up") up_genes else dn_genes
  if (length(genes) < 10) {
    cat("  ", direction, "Directional gene number <10, Skip enrichment\n")
    next}
  eg <- tryCatch(
    bitr(genes, "SYMBOL", "ENTREZID", OrgDb = org.Rn.eg.db)$ENTREZID,
    error = function(e) NULL, warning = function(w) NULL)
  if (is.null(eg) || length(eg) < 5) next
  go <- tryCatch(
    enrichGO(
      eg, OrgDb = org.Rn.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05),
    error = function(e) NULL)
  kk <- tryCatch(
    enrichKEGG(eg, organism = "rno", qvalueCutoff = 0.05),
    error = function(e) NULL)
  enrich_results[[paste0("GO_", direction)]] <- go
  enrich_results[[paste0("KEGG_", direction)]] <- kk
  if (!is.null(go) && nrow(as.data.frame(go)) > 0) {
    enrich_plots[[paste0("GO_", direction)]] <- dotplot(go, showCategory = 15) +
      ggtitle(paste("GO BP -", direction))}
  if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
    enrich_plots[[paste0("KEGG_", direction)]] <- dotplot(kk, showCategory = 15) +
      ggtitle(paste("KEGG -", direction))}}
if (length(enrich_plots) > 0) {
  safe_ggsave(
    file.path(out_dir, "C_TREZ_analysis", "GO_KEGG_TREZ.pdf"),
    wrap_plots(enrich_plots, ncol = 2),
    width = 14, height = 4 * ceiling(length(enrich_plots) / 2))
  saveRDS(enrich_results,
          file.path(out_dir, "C_TREZ_analysis", "enrich_results.rds"))} else {
            cat("The enrichment result is empty, skip GO_KEGG_TREZ.pdf\n")}
saveRDS(obj_trez, file.path(out_dir, "C_TREZ_analysis", "obj_TREZ.rds"))

safe_ggsave(
  file.path(out_dir, "A_clustering", "umap_clusters.pdf"), 
  DimPlot(obj, label = TRUE) + ggtitle("10 clusters (K-means)"),
  width = 7, height = 6)

safe_ggsave(
  file.path(out_dir, "A_clustering", "umap_split_group.pdf"),
  DimPlot(obj, split.by = "group", label = TRUE),
  width = 12, height = 6)

obj$is_trez <- ifelse(
  obj$paper_label %in% c("C3_HealthyMyelin", "C5_RootZone", "C7_TNMyelin"),
  obj$paper_label, "other")
safe_ggsave(
  file.path(out_dir, "A_clustering", "umap_TREZ_highlight.pdf"), 
  DimPlot(obj, group.by = "is_trez",
          cols = COL_CLUSTER_LABEL[c("C3_HealthyMyelin", "C5_RootZone", "C7_TNMyelin", "other")]) +
    ggtitle("TREZ region (C3 + C5 + C7)"), width = 7, height = 6)
comp_long <- as.data.frame(table(Cluster = Idents(obj), Group = obj$group))
safe_ggsave(
  file.path(out_dir, "A_clustering", "composition_stack.pdf"),
  ggplot(
    comp_long, aes(Cluster, Freq, fill = Group)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = cols_g) +
    theme_bw() + ylab("Proportion") +
    ggtitle("Per-cluster HC vs TN composition"),
  width = 8, height = 5)

base <- "."
in_obj  <- file.path(base, "results_paper_TREZ/C_TREZ_analysis/obj_TREZ.rds")
out_dir <- file.path(base, "results_OPC_lineage_gated")

for (sub in c("01_gating", "02_scores", "03_DEG", "04_main_figure")) {
  dir.create(file.path(out_dir, sub), showWarnings = FALSE, recursive = TRUE)}
set.seed(42)
options(future.globals.maxSize = 64 * 1024^3)

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(MAST)
  library(EnhancedVolcano)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
  library(ggrepel)
  library(svglite)})

cols_g <- COL_GROUP[c("HC", "TN")]
stopifnot(file.exists(in_obj))
obj <- readRDS(in_obj)
signature_preflight <- function(obj, sig_list) {
  genes_obj <- rownames(obj)
  bind_rows(lapply(names(sig_list), function(nm) {
    sig <- unique(sig_list[[nm]])
    found <- intersect(sig, genes_obj)
    missing <- setdiff(sig, genes_obj)
    data.frame(
      signature = nm,
      n_total = length(sig),
      n_found = length(found),
      missing_rate = round(length(missing) / length(sig), 4),
      found_genes = paste(found, collapse = ", "),
      missing_genes = paste(missing, collapse = ", "),
      stringsAsFactors = FALSE)}))}

DefaultAssay(obj) <- "RNA"
lineage_gate_genes <- c("Olig1", "Olig2")
stopifnot(all(lineage_gate_genes %in% rownames(obj)))

gate_df <- fetch_seurat_data(obj, vars = lineage_gate_genes, layer = "counts")
olig1_pos <- gate_df[, "Olig1"] >= 1
olig2_pos <- gate_df[, "Olig2"] >= 1
names(olig1_pos) <- rownames(gate_df)
names(olig2_pos) <- rownames(gate_df)
lineage_pos <- olig1_pos | olig2_pos
obj$Olig1_pos <- olig1_pos
obj$Olig2_pos <- olig2_pos
obj$lineage_pos <- lineage_pos
gate_stats <- obj@meta.data %>%
  dplyr::group_by(sample, group) %>%
  dplyr::summarise(
    total_spots = n(),
    Olig1_pos = sum(Olig1_pos),
    Olig2_pos = sum(Olig2_pos),
    lineage_pos = sum(lineage_pos),
    lineage_pct = round(100 * sum(lineage_pos) / n(), 2),
    .groups = "drop")
write.csv(gate_stats, file.path(out_dir, "01_gating", "gating_summary_per_sample.csv"), row.names = FALSE)

gate_by_grp <- obj@meta.data %>%
  group_by(group) %>%
  summarise(total = n(), lineage = sum(lineage_pos),
            lineage_pct = 100 * lineage / total, .groups = "drop")

p_gate_bar <- ggplot(gate_stats, aes(sample, lineage_pct, fill = group)) +
  geom_col() +
  scale_fill_manual(values = cols_g) +
  geom_text(aes(label = sprintf("%d/%d", lineage_pos, total_spots)),
            vjust = -0.3, size = 3) +
  theme_bw() + ylab("Lineage-positive spots (%)") +
  ggtitle("Lineage-positive spots (Olig1+ OR Olig2+)"); p_gate_bar
safe_ggsave(file.path(out_dir, "01_gating", "gating_barplot.pdf"), p_gate_bar, width = 7, height = 5)
obj_gated <- subset(obj, subset = lineage_pos == TRUE)
check_groups_present(obj_gated@meta.data)
cat(sprintf(
  "\nCompleted gating: kept %d spots (HC %d + TN %d)\n",
  ncol(obj_gated),
  sum(obj_gated$group == "HC"),
  sum(obj_gated$group == "TN")))

DefaultAssay(obj_gated) <- "RNA"
obj_gated <- NormalizeData(obj_gated, verbose = FALSE)

opc_pure <- c("Neu4", "Pcdh15", "Nkx2-2", "Gpr17", "Pdgfra", "Cspg4")
ol_pure  <- c("Mog", "Mag", "Mbp", "Mobp", "Plp1", "Cnp")
opc_pure_p <- intersect(opc_pure, rownames(obj_gated))
ol_pure_p  <- intersect(ol_pure,  rownames(obj_gated))
sig_qc <- signature_preflight(obj_gated, list(
  gate = lineage_gate_genes,
  OPC_pure = opc_pure,
  OL_pure = ol_pure))
write.csv(sig_qc, file.path(out_dir, "02_scores", "signature_qc_gated.csv"), row.names = FALSE)
cat("signature pre-flight:\n"); print(sig_qc[, c("signature", "n_total", "n_found", "missing_rate")])
cat("OPC_pure (", length(opc_pure_p), "): ", paste(opc_pure_p, collapse = ", "), "\n", sep = "")
cat("OL_pure  (", length(ol_pure_p), "): ", paste(ol_pure_p,  collapse = ", "), "\n", sep = "")
if (length(opc_pure_p) < 2 || length(ol_pure_p) < 2) {
  stop("Too few detected genes in OPC_pure or OL_pure signature after gating.")}

obj_gated <- AddModuleScore_UCell(
  obj_gated,
  features = list(
    OPC_pure = opc_pure_p,
    OL_pure  = ol_pure_p),
  assay = "RNA",
  slot = "counts",
  maxRank = 1500,
  missing_genes = "skip",
  name = "",
  BPPARAM = ucell_bpparam())
obj_gated$OPC_OL_ratio_pure <- obj_gated$OPC_pure / (obj_gated$OL_pure + 0.01)

m <- obj_gated@meta.data
metrics <- c("OPC_pure", "OL_pure", "OPC_OL_ratio_pure")
stats_df <- data.frame(
  metric = metrics,
  HC_median = sapply(metrics, function(g) median(m[[g]][m$group == "HC"])),
  TN_median = sapply(metrics, function(g) median(m[[g]][m$group == "TN"])),
  HC_mean = sapply(metrics, function(g) mean(m[[g]][m$group == "HC"])),
  TN_mean = sapply(metrics, function(g) mean(m[[g]][m$group == "TN"])),
  p_wilcox_spot_level = sapply(metrics, function(g) safe_wilcox_p(m[[g]], m$group)))
stats_df$log2FC_median <- log2((stats_df$TN_median + 0.01) / (stats_df$HC_median + 0.01))
stats_df$log2FC_mean <- log2((stats_df$TN_mean + 0.01) / (stats_df$HC_mean + 0.01))
write.csv(stats_df, file.path(out_dir, "02_scores", "UCell_stats_gated.csv"), row.names = FALSE)
print(stats_df)

get_metric_p_label <- function(metric) {
  p <- stats_df$p_wilcox_spot_level[stats_df$metric == metric]
  format_p_label(p[1])}
make_gated_score_violin <- function(
    metric, title, y_label,
    pseudo_log = FALSE, show_nonzero_points = FALSE) {
  plot_df <- data.frame(
    group = factor(as.character(m$group), levels = c("HC", "TN")),
    score = as.numeric(m[[metric]]),
    stringsAsFactors = FALSE)
  plot_df <- plot_df[!is.na(plot_df$group) & is.finite(plot_df$score), , drop = FALSE]
  point_df <- if (show_nonzero_points) {
    plot_df[plot_df$score > 0, , drop = FALSE]
  } else {
    plot_df[FALSE, , drop = FALSE]}
  y_range <- range(plot_df$score, na.rm = TRUE)
  y_span <- diff(y_range)
  if (!is.finite(y_span) || y_span <= 0) y_span <- max(abs(y_range), 0.1)
  annot_df <- data.frame(
    x = 1.10,
    xend = 1.90,
    xmid = 1.50,
    line_y = y_range[2] + 0.35 * y_span,
    label_y = y_range[2] + 0.47 * y_span,
    panel_top = y_range[2] + 0.68 * y_span,
    p_label = get_metric_p_label(metric),
    stringsAsFactors = FALSE)
  p <- ggplot(plot_df, aes(x = group, y = score, fill = group)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.78,
                color = "grey25", linewidth = 0.35, width = 0.72,
                na.rm = TRUE) +
    scale_fill_manual(values = cols_g, drop = FALSE, name = NULL) +
    labs(x = "Identity", y = y_label, title = title) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5))
  p <- p +
    geom_blank(data = annot_df, aes(x = xmid, y = panel_top),
               inherit.aes = FALSE) +
    geom_segment(data = annot_df,
                 aes(x = x, xend = xend, y = line_y, yend = line_y),
                 inherit.aes = FALSE, linewidth = 0.45) +
    geom_text(data = annot_df, aes(x = xmid, y = label_y, label = p_label),
              inherit.aes = FALSE, size = 3.1, fontface = "bold")
  if (pseudo_log) {
    p <- p + scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 0.001),
      expand = expansion(mult = c(0.02, 0.04)))} else {
        p <- p + scale_y_continuous(expand = expansion(mult = c(0.02, 0.04)))}
  p}
p1 <- make_gated_score_violin(
  "OPC_pure", "OPC_pure (lineage-gated)",
  "UCell score (pseudo-log scale)",
  pseudo_log = TRUE, show_nonzero_points = TRUE)
p2 <- make_gated_score_violin(
  "OL_pure", "OL_pure (lineage-gated)",
  "UCell score",
  pseudo_log = FALSE, show_nonzero_points = TRUE)
p3 <- make_gated_score_violin(
  "OPC_OL_ratio_pure", "OPC/OL score ratio",
  "OPC/OL ratio (pseudo-log scale)",
  pseudo_log = TRUE, show_nonzero_points = TRUE)
p_violin <- p1 | p2 | p3; p_violin
safe_ggsave(file.path(out_dir, "02_scores", "violin_gated_three.pdf"), p_violin, width = 12, height = 4.5)
ggsave(
  filename = file.path(out_dir, "02_scores", "violin_gated_three.svg"),
  plot = p_violin,
  device = svglite,  
  width = 12, 
  height = 4.5,
  dpi = 300,
  system_fonts = list(
    sans = "Arial", 
    serif = "Times New Roman"),)

sample_levels <- m %>%
  distinct(sample, group) %>%
  arrange(factor(group, levels = c("HC", "TN")), sample) %>%
  pull(sample); sample_levels

make_per_sample_gated_boxplot <- function(
    metric, title, y_label,       
    pseudo_log = FALSE,
    highlight_nonzero = FALSE) {
  plot_df <- data.frame(
    sample = factor(as.character(m$sample), levels = sample_levels),
    group = factor(as.character(m$group), levels = c("HC", "TN")),
    score = as.numeric(m[[metric]]),
    stringsAsFactors = FALSE)
  plot_df <- plot_df[!is.na(plot_df$sample) & !is.na(plot_df$group) &
                       is.finite(plot_df$score), , drop = FALSE]
  plot_df$sample_x <- as.numeric(plot_df$sample)
  point_df <- if (highlight_nonzero) {
    plot_df[plot_df$score > 0, , drop = FALSE]} else {plot_df}
  y_range <- range(plot_df$score, na.rm = TRUE)
  y_span <- diff(y_range)
  if (!is.finite(y_span) || y_span <= 0) {
    y_span <- max(abs(y_range), 0.1, na.rm = TRUE)}
  annot_df <- data.frame(
    x = 1.15,
    xend = length(sample_levels) - 0.15,
    xmid = (1.15 + length(sample_levels) - 0.15) / 2,
    line_y = y_range[2] + 0.24 * y_span,
    label_y = y_range[2] + 0.43 * y_span,
    panel_top = y_range[2] + 0.78 * y_span,
    p_label = get_metric_p_label(metric),
    stringsAsFactors = FALSE)
  p <- ggplot(plot_df, aes(x = sample_x, y = score, fill = group)) +
    geom_boxplot(aes(group = sample), width = 0.55, outlier.shape = NA,
                 alpha = 0.58, color = "grey25", linewidth = 0.32,
                 na.rm = TRUE) +
    geom_point(data = point_df, aes(color = group),
               position = position_jitter(width = 0.10, height = 0, seed = 42),
               size = if (highlight_nonzero) 0.50 else 0.28,
               alpha = if (highlight_nonzero) 0.55 else 0.20,
               show.legend = FALSE, na.rm = TRUE) +
    geom_blank(data = annot_df, aes(x = xmid, y = panel_top),
               inherit.aes = FALSE) +
    geom_segment(data = annot_df,
                 aes(x = x, xend = xend, y = line_y, yend = line_y),
                 inherit.aes = FALSE, linewidth = 0.45) +
    geom_text(data = annot_df,
              aes(x = xmid, y = label_y, label = p_label),
              inherit.aes = FALSE, size = 3.0, fontface = "bold",
              vjust = 0) +
    scale_x_continuous(breaks = seq_along(sample_levels),
                       labels = sample_levels,
                       limits = c(0.55, length(sample_levels) + 0.45)) +
    scale_fill_manual(values = cols_g, drop = FALSE, name = NULL) +
    scale_color_manual(values = cols_g, drop = FALSE, name = NULL) +
    labs(x = NULL, y = y_label, title = title) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5))
  if (pseudo_log) {
    p <- p + scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 0.001),
      expand = expansion(mult = c(0.03, 0.04)))} else {
        p <- p + scale_y_continuous(expand = expansion(mult = c(0.03, 0.04)))}
  p}

p_bysample <- make_per_sample_gated_boxplot(
  "OPC_pure", "OPC_pure", "UCell score (pseudo-log scale)",
  pseudo_log = TRUE, highlight_nonzero = TRUE) | make_per_sample_gated_boxplot(
    "OL_pure", "OL_pure", "UCell score",
    pseudo_log = FALSE, highlight_nonzero = FALSE) | make_per_sample_gated_boxplot(
      "OPC_OL_ratio_pure", "OPC/OL score ratio",
      "OPC/OL ratio (pseudo-log scale)", 
      pseudo_log = TRUE, highlight_nonzero = TRUE)
p_bysample <- p_bysample +
  plot_annotation(title = "Per-sample (lineage-gated): HC / TN1 / TN2 / TN3") &
  theme(plot.title = element_text(face = "bold", hjust = 0.5)); p_bysample
safe_ggsave(file.path(out_dir, "02_scores", "boxplot_per_sample_gated.pdf"), p_bysample, width = 12, height = 4)

sample_scores <- m %>%
  group_by(sample, group) %>%
  summarise(
    n_spots = n(),
    OPC_pure = mean(OPC_pure),
    OL_pure = mean(OL_pure),
    OPC_OL_ratio_pure = mean(OPC_OL_ratio_pure),
    .groups = "drop"); sample_scores
write.csv(sample_scores, file.path(out_dir, "02_scores", "sample_level_scores_gated.csv"), row.names = FALSE)

p_samp <- sample_scores %>%
  dplyr::select(sample, group, OPC_pure, OL_pure, OPC_OL_ratio_pure) %>%
  tidyr::pivot_longer(-c(sample, group), names_to = "metric", values_to = "score") %>%
  dplyr::mutate(metric = factor(metric, levels = c("OPC_pure","OL_pure","OPC_OL_ratio_pure"))) %>%
  ggplot(aes(group, score, fill = group)) +
  geom_bar(stat = "summary", fun = "mean", alpha = 0.5, width = 0.6) +
  geom_jitter(width = 0.1, size = 3, aes(color = sample)) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = cols_g) +
  theme_bw() +
  ggtitle("Sample-level (n=1 HC, n=3 TN) — lineage-gated"); p_samp
safe_ggsave(file.path(out_dir, "02_scores", "sample_level_barplot_gated.pdf"), p_samp, width = 12, height = 4)

DefaultAssay(obj_gated) <- "SCT"
obj_gated <- PrepSCTFindMarkers(obj_gated, verbose = FALSE)
Idents(obj_gated) <- "group"
deg <- suppressMessages(quiet_known_package_warnings(FindMarkers(
  obj_gated, ident.1 = "TN", ident.2 = "HC",
  test.use = "MAST", assay = "SCT",
  logfc.threshold = 0, min.pct = 0.05,
  verbose = FALSE)))
deg$gene <- rownames(deg)
write.csv(deg, file.path(out_dir, "03_DEG", "DEG_MAST_gated_all.csv"), row.names = FALSE)
deg_sig <- subset(deg, p_val_adj < 0.001 & abs(avg_log2FC) >= log2(1.5))
write.csv(deg_sig, file.path(out_dir, "03_DEG", "DEG_MAST_gated_sig.csv"), row.names = FALSE)
deg_focus <- deg[deg$gene %in% c(opc_pure_p, ol_pure_p), ]
deg_focus$set <- ifelse(deg_focus$gene %in% opc_pure_p, "OPC_pure", "OL_pure")
deg_focus <- deg_focus[order(deg_focus$set, deg_focus$p_val_adj), ]; deg_focus
write.csv(deg_focus, file.path(out_dir, "03_DEG", "DEG_focus_OPC_OL.csv"), row.names = FALSE)

deg$p_val_adj_plot <- floor_zero_pvalues(deg$p_val_adj)
volcano_labels <- select_volcano_labels(
  deg,
  p_col = "p_val_adj_plot",
  fc_col = "avg_log2FC",
  candidate_genes = c(opc_pure_p, ol_pure_p),
  max_labels = 4,
  max_per_direction = 0)
p_vol <- quiet_known_package_warnings(EnhancedVolcano(
  deg, lab = deg$gene,
  selectLab = volcano_labels,
  x = "avg_log2FC", y = "p_val_adj_plot",
  pCutoff  = 0.001, FCcutoff = log2(1.5),
  title = "Lineage-gated TREZ: TN vs HC",
  subtitle = sprintf("n = %d spots; compact marker labels shown", ncol(obj_gated)),
  drawConnectors = TRUE, max.overlaps = 15,
  labSize = 3.1, pointSize = 1.25, colAlpha = 0.5,
  col = c(COL_NEUTRAL_PALE, COL_NEUTRAL_MID, COL_NEUTRAL_DARK, COL_CONTROL))); p_vol
safe_ggsave(file.path(out_dir, "03_DEG", "volcano_gated.pdf"), p_vol, width = 9, height = 7.2)

DefaultAssay(obj_gated) <- "RNA"
p_dot <- DotPlot(
  obj_gated, features = c(opc_pure_p, ol_pure_p),
  group.by = "group", scale = FALSE,
  cols = c(COL_DIVERGING["low"], COL_DIVERGING["high"])) +
  RotatedAxis() +
  ggtitle("OPC/OL markers (lineage-gated) - HC vs TN"); p_dot
safe_ggsave(file.path(out_dir, "03_DEG", "dotplot_gated.pdf"), p_dot, width = 9, height = 4)
key_gene_targets <- unique(c(lineage_gate_genes, opc_pure, ol_pure))
key_genes <- resolve_feature_names(key_gene_targets, rownames(obj_gated[["RNA"]]))
if (length(attr(key_genes, "missing_features")) > 0) {
  cat("Key genes not found: ", paste(attr(key_genes, "missing_features"), collapse = ", "), "\n", sep = "")}

if (length(key_genes) > 0) {
  p_keys <- make_hc_tn_gene_violin(
    obj_gated,
    features = key_gene_targets,
    group_col = "group",
    assay = "RNA",
    layer = "data",
    cols = cols_g,
    ncol = 3,
    title = "Key OPC / OL genes after lineage gating: HC vs TN",
    subtitle = "Olig1/Olig2 are shown as gate-marker QC; P values are descriptive spot-level Wilcoxon tests")
  write.csv(attr(p_keys, "p_values"), file.path(out_dir, "03_DEG", "violin_key_gene_pvalues_gated.csv"), row.names = FALSE)
  safe_ggsave(file.path(out_dir, "03_DEG", "violin_key_genes_gated.pdf"), p_keys, width = 12, height = max(5.2, 2.45 * ceiling(length(key_genes) / 3)))}

p_main <- (p_gate_bar | p_violin) / (p_samp | p_vol) +
  plot_annotation(
    title = "Lineage-gated OPC/OL score summary (Olig1+ OR Olig2+)",
    subtitle = sprintf("Gated: %d / %d spots (%.1f%%) | spot-level displays are descriptive",
                       ncol(obj_gated), ncol(obj),
                       100 * ncol(obj_gated) / ncol(obj)),
    tag_levels = "A"); p_main
safe_ggsave(file.path(out_dir, "04_main_figure", "main_summary_gated.pdf"), p_main, width = 18, height = 12)
saveRDS(obj_gated, file.path(out_dir, "obj_gated.rds"))

base <- "."
in_obj <- file.path(base, "results_OPC_lineage_gated/obj_gated.rds")
out_dir <- file.path(base, "results_OPC_volcano_enhance")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(42)
options(future.globals.maxSize = 64 * 1024^3)

suppressPackageStartupMessages({
  library(Seurat)
  library(EnhancedVolcano)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
  library(ggrepel)})
stopifnot(file.exists(in_obj))
obj_gated <- readRDS(in_obj)
check_groups_present(obj_gated@meta.data)

opc_pure <- c("Neu4", "Pcdh15", "Nkx2-2", "Gpr17", "Pdgfra", "Cspg4")
ol_pure  <- c("Mog", "Mag", "Mbp", "Mobp", "Plp1", "Cnp")
opc_pure_p <- intersect(opc_pure, rownames(obj_gated))
ol_pure_p  <- intersect(ol_pure,  rownames(obj_gated))
cat("OPC_pure (", length(opc_pure_p), "): ",
    paste(opc_pure_p, collapse = ", "), "\n", sep = "")
cat("OL_pure  (", length(ol_pure_p), "): ",
    paste(ol_pure_p,  collapse = ", "), "\n\n", sep = "")
if (length(opc_pure_p) < 2 || length(ol_pure_p) < 2) {
  stop("Too few detected genes in OPC_pure or OL_pure signature for focused plots.")}

DefaultAssay(obj_gated) <- "RNA"
Idents(obj_gated) <- "group"

deg_wx <- suppressMessages(FindMarkers(
  obj_gated, ident.1 = "TN", ident.2 = "HC",
  test.use = "wilcox",
  logfc.threshold = 0, min.pct = 0.05,
  verbose = FALSE))
deg_wx$gene <- rownames(deg_wx)
write.csv(deg_wx, file.path(out_dir, "DEG_wilcox_all.csv"), row.names = FALSE)

deg_focus_wx <- deg_wx[deg_wx$gene %in% c(opc_pure_p, ol_pure_p), ]
deg_focus_wx$set <- ifelse(deg_focus_wx$gene %in% opc_pure_p, "OPC_pure", "OL_pure")
deg_focus_wx <- deg_focus_wx[order(deg_focus_wx$set, deg_focus_wx$p_val_adj), ]
write.csv(deg_focus_wx, file.path(out_dir, "DEG_wilcox_OPC_OL_focus.csv"), row.names = FALSE)
print(deg_focus_wx)

mast_path <- file.path(base, "results_OPC_lineage_gated/03_DEG/DEG_focus_OPC_OL.csv")
if (file.exists(mast_path)) {
  deg_focus_mast <- read.csv(mast_path)
  cmp <- data.frame(
    gene = deg_focus_wx$gene,
    set = deg_focus_wx$set,
    log2FC_wilcox = deg_focus_wx$avg_log2FC,
    padj_wilcox = deg_focus_wx$p_val_adj) %>%
    dplyr::left_join(data.frame(
      gene = deg_focus_mast$gene,
      log2FC_mast = deg_focus_mast$avg_log2FC,
      padj_mast = deg_focus_mast$p_val_adj), by = "gene") %>%
    dplyr::arrange(set, padj_wilcox)
  cmp$sig_wilcox <- ifelse(
    cmp$padj_wilcox < 0.001, "***",
    ifelse(cmp$padj_wilcox < 0.01,  "**",
           ifelse(cmp$padj_wilcox < 0.05,  "*",  "ns")))
  cmp$sig_mast <- ifelse(
    cmp$padj_mast < 0.001, "***", 
    ifelse(cmp$padj_mast < 0.01,  "**",
           ifelse(cmp$padj_mast < 0.05,  "*",   "ns")))
  write.csv(cmp, file.path(out_dir, "MAST_vs_Wilcoxon_comparison.csv"), row.names = FALSE)
  cat("\nMAST vs Wilcoxon:\n"); print(cmp)}

deg_wx$p_val_adj_plot <- floor_zero_pvalues(deg_wx$p_val_adj)
volcano_labels_wx <- select_volcano_labels(
  deg_wx,
  p_col = "p_val_adj_plot",
  fc_col = "avg_log2FC",
  candidate_genes = c(opc_pure_p, ol_pure_p),
  max_labels = 4,
  max_per_direction = 0)
p_vol_wx <- quiet_known_package_warnings(EnhancedVolcano(
  deg_wx, lab = deg_wx$gene,
  selectLab = volcano_labels_wx,
  x = "avg_log2FC", y = "p_val_adj_plot",
  pCutoff = 0.05, FCcutoff = log2(1.2),
  title = "Lineage-gated TREZ: TN vs HC (Wilcoxon)",
  subtitle = sprintf("n = %d spots; compact marker labels shown", ncol(obj_gated)),
  drawConnectors = TRUE, max.overlaps = 15,
  labSize = 3.1, pointSize = 1.25, colAlpha = 0.5,
  col = c(COL_NEUTRAL_PALE, COL_NEUTRAL_MID, COL_NEUTRAL_DARK, COL_CONTROL))); p_vol_wx
safe_ggsave(file.path(out_dir, "volcano_wilcox.pdf"), p_vol_wx, width = 9, height = 7.2)

focus_genes <- c(opc_pure_p, ol_pure_p)
expr_df <- fetch_seurat_data(obj_gated, vars = c(focus_genes, "group"), layer = "data")
boot_log2fc <- function(gene, df, n_boot = 1000) {
  e_hc <- df[[gene]][df$group == "HC"]
  e_tn <- df[[gene]][df$group == "TN"]
  est  <- log2((mean(e_tn) + 0.01) / (mean(e_hc) + 0.01))
  boots <- replicate(n_boot, {
    h <- mean(sample(e_hc, replace = TRUE))
    t <- mean(sample(e_tn, replace = TRUE))
    log2((t + 0.01) / (h + 0.01))})
  data.frame(
    gene = gene,
    est = est,
    lower = as.numeric(quantile(boots, 0.025)),
    upper = as.numeric(quantile(boots, 0.975)))}
forest_df <- do.call(rbind, lapply(focus_genes, boot_log2fc, df = expr_df))
forest_df$set <- ifelse(forest_df$gene %in% opc_pure_p, "OPC_pure", "OL_pure")
forest_df$set <- factor(forest_df$set, levels = c("OL_pure", "OPC_pure"))
forest_df$direction <- ifelse(
  forest_df$lower > 0, "TN up; interval excludes 0", 
  ifelse(forest_df$upper < 0, "TN down; interval excludes 0", "interval includes 0"))
forest_df <- forest_df %>% arrange(set, est)
forest_df$gene <- factor(forest_df$gene, levels = forest_df$gene); forest_df
write.csv(forest_df, file.path(out_dir, "forest_data.csv"), row.names = FALSE)

p_forest <- ggplot(forest_df, aes(x = est, y = gene, color = direction)) +
  geom_vline(xintercept = 0, lty = 2, color = "grey40") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.25, orientation = "y") +
  geom_point(size = 3) +
  facet_grid(set ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(
    values = c(
      "TN up; interval excludes 0" = COL_CONTROL,
      "TN down; interval excludes 0" = COL_DISEASE,
      "interval includes 0" = COL_NEUTRAL_DARK)) +
  labs(
    x = "log2 Fold Change (TN / HC, 95% bootstrap CI)", y = NULL,
    title = "OPC / OL gene-level forest plot (lineage-gated)",
    subtitle = "Point = log2FC; bars = 95% spot-level bootstrap CI (1000 resamples)") +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold"),
    legend.position = "top"); p_forest
safe_ggsave(file.path(out_dir, "forest_OPC_OL.pdf"), p_forest, width = 8, height = 6)

bar_df <- deg_focus_wx %>%
  dplyr::mutate(
    set = factor(set, levels = c("OL_pure", "OPC_pure")),
    sig = case_when(
      p_val_adj < 0.001 ~ "***",
      p_val_adj < 0.01  ~ "**",
      p_val_adj < 0.05  ~ "*",
      TRUE ~ "ns"),
    label = paste0(gene, " ", sig)) %>%
  dplyr::arrange(set, avg_log2FC); bar_df
bar_df$gene <- factor(bar_df$gene, levels = bar_df$gene)

p_bar <- ggplot(bar_df, aes(x = avg_log2FC, y = gene, fill = avg_log2FC > 0)) +
  geom_col(alpha = 0.85) +
  geom_vline(xintercept = 0, lty = 1, color = "black") +
  facet_grid(set ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(
    values = COL_SIGN,
    labels = c("TRUE" = "TN up", "FALSE" = "TN down"),
    name   = NULL) +
  labs(
    x = "Wilcoxon avg_log2FC (TN vs HC)", y = NULL,
    title = "OPC / OL direction (Wilcoxon, lineage-gated)",
    subtitle = "Direction and effect size shown; inferential scope remains spot-level descriptive") +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold"),
    legend.position = "top"); p_bar
safe_ggsave(file.path(out_dir, "barplot_log2FC.pdf"), p_bar, width = 8, height = 6)

p_main <- p_vol_wx / (p_forest | p_bar) +
  plot_annotation(
    title = "OPC up / OL down in TREZ - enhanced visualization",
    subtitle = "Wilcoxon DEG on lineage-gated spots; spot-level output is descriptive",
    tag_levels = "A"); p_main
safe_ggsave(file.path(out_dir, "main_enhanced_volcano.pdf"), p_main, width = 16, height = 14)

# stage_ucell ####
base <- "."
in_obj <- file.path(base, "results_OPC_lineage_gated/obj_gated.rds")
out_dir <- file.path(base, "results_OPC_stage")
for (sub in c("01_stage_scores", "02_stats", "03_plots")) {
  dir.create(file.path(out_dir, sub), showWarnings = FALSE, recursive = TRUE)}
set.seed(42)
options(future.globals.maxSize = 64 * 1024^3)

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
  library(ggrepel)})

stopifnot(file.exists(in_obj))
obj <- readRDS(in_obj)
check_groups_present(obj@meta.data)
signature_preflight <- function(obj, sig_list) {
  genes_obj <- rownames(obj)
  dplyr::bind_rows(lapply(names(sig_list), function(nm) {
    sig <- unique(sig_list[[nm]])
    found <- intersect(sig, genes_obj)
    missing <- setdiff(sig, genes_obj)
    data.frame(
      signature = nm,
      n_total = length(sig),
      n_found = length(found),
      missing_rate = round(length(missing) / length(sig), 4),
      found_genes = paste(found, collapse = ", "),
      missing_genes = paste(missing, collapse = ", "),
      stringsAsFactors = FALSE)}))}

stage_sets <- list(
  OPC = c("Pdgfra", "Cspg4", "Olig1", "Olig2"),
  COP = c("Bmp4", "Neu4", "Sox6", "Gpr17", "Fyn"),
  NFOL = c("Itpr2", "Enpp6", "Tcf7l2", "Bcas1", "Casr"),
  MFOL = c("Cnp", "Mal", "Plp1", "Opalin"),
  MOL = c("Mog", "Mag", "Mbp", "Mobp", "Klk6", "Ermn"))

stage_qc <- signature_preflight(obj, stage_sets)
write.csv(stage_qc, file.path(out_dir, "01_stage_scores", "stage_signature_qc.csv"), row.names = FALSE)
stage_sets_p <- lapply(stage_sets, function(g) intersect(g, rownames(obj)))
for (s in names(stage_sets_p)) {
  cat(sprintf("  %-5s (%d): %s\n", s, length(stage_sets_p[[s]]), paste(stage_sets_p[[s]], collapse = ", ")))}

min_gene_check <- sapply(stage_sets_p, length)
if (any(min_gene_check < 2)) {
  warning(
    "Stage-specific gene sets < 2, UCell might be unstable: ", 
    paste(names(min_gene_check)[min_gene_check < 2], collapse = ", "))}
if (any(stage_qc$missing_rate > 0.2)) {
  warning(
    "If the gene set dropout rate exceeds 20% at certain stages, please clearly disclose this in the Methods or Supplementary Information.")}

DefaultAssay(obj) <- "RNA"
ucell_max_rank <- 1500L
obj <- AddModuleScore_UCell(
  obj,
  features = stage_sets_p,
  maxRank = ucell_max_rank,
  slot = "counts",
  missing_genes = "skip",
  assay = "RNA",
  name = "", 
  BPPARAM = ucell_bpparam())
stage_order <- c("OPC", "COP", "NFOL", "MFOL", "MOL")
eps_main <- 1e-2
eps_grid <- c(1e-3, 1e-2, 1e-1)
obj$MaturationIndex <- log2(
  (obj$MFOL + obj$MOL + eps_main) /
    (obj$OPC  + obj$COP + eps_main))
obj$EarlyLateRatio <- log2(
  (obj$OPC + obj$COP + eps_main) /
    (obj$NFOL + obj$MFOL + obj$MOL + eps_main))
mi_sensitivity <- bind_rows(lapply(eps_grid, function(eps) {
  data.frame(
    sample = obj$sample,
    group = obj$group,
    eps = eps,
    MaturationIndex = log2(
      (obj$MFOL + obj$MOL + eps) /
        (obj$OPC  + obj$COP + eps)))}))
write.csv(
  mi_sensitivity,
  file.path(out_dir, "02_stats", "maturation_index_eps_sensitivity_per_spot.csv"),
  row.names = FALSE)

stage_margin_cutoff <- 0.05
stage_assignment <- assign_dominant_stage(
  obj@meta.data, stage_order,
  margin_cutoff = stage_margin_cutoff,
  use_z = TRUE)
obj$dominant_stage <- stage_assignment$dominant_stage
obj$stage_margin <- stage_assignment$stage_margin
obj$stage_ambiguous <- stage_assignment$stage_ambiguous
write.csv(stage_assignment$score_matrix, file.path(out_dir, "01_stage_scores", "stage_score_zmatrix_for_classification.csv"))
score_cols <- c("OPC", "COP", "NFOL", "MFOL", "MOL", "MaturationIndex", "EarlyLateRatio")
qc_cols <- c("stage_margin")
score_df <- obj@meta.data %>%
  dplyr::select(sample, group, dominant_stage, stage_ambiguous, all_of(score_cols), all_of(qc_cols))
write.csv(score_df, file.path(out_dir, "01_stage_scores", "per_spot_stage_scores.csv"))

m <- obj@meta.data
stats_df <- data.frame(
  stage = score_cols,
  HC_median = sapply(score_cols, function(g) median(m[[g]][m$group == "HC"])),
  TN_median = sapply(score_cols, function(g) median(m[[g]][m$group == "TN"])),
  HC_mean = sapply(score_cols, function(g) mean(m[[g]][m$group == "HC"])),
  TN_mean = sapply(score_cols, function(g) mean(m[[g]][m$group == "TN"])),
  p_wilcox_spot_level = sapply(score_cols, function(g) safe_wilcox_p(m[[g]], m$group)))
stats_df$delta_mean <- stats_df$TN_mean - stats_df$HC_mean
stats_df$direction <- ifelse(stats_df$delta_mean > 0, "TN up",
                             ifelse(stats_df$delta_mean < 0, "TN down", "-")); stats_df
write.csv(stats_df, file.path(out_dir, "02_stats", "stage_stats_HC_vs_TN.csv"), row.names = FALSE)
sample_scores <- m %>%
  group_by(sample, group) %>%
  summarise(
    across(all_of(score_cols), mean), 
    stage_ambiguous_pct = 100 * mean(stage_ambiguous, na.rm = TRUE),
    n_spots = n(), .groups = "drop"); sample_scores
write.csv(sample_scores, file.path(out_dir, "02_stats", "stage_sample_level.csv"), row.names = FALSE)

cols_g <- COL_GROUP[c("HC", "TN")]
m_long <- m %>%
  dplyr::select(sample, group, all_of(stage_order)) %>%
  tidyr::pivot_longer(cols = all_of(stage_order), names_to = "stage", values_to = "score") %>%
  dplyr::mutate(stage = factor(stage, levels = stage_order), stage_x = as.numeric(stage)); m_long
p_label_one <- function(p) {
  if (!is.finite(p)) return("p=NA")
  if (p < 1e-4) return("p<0.0001")
  paste0("p=", sprintf("%.4f", p))}
score_range <- range(m_long$score, na.rm = TRUE)
score_span <- diff(score_range)
if (!is.finite(score_span) || score_span == 0) score_span <- 0.1
annot_line_y <- score_range[2] + 0.24 * score_span
annot_label_y <- score_range[2] + 0.32 * score_span
p_annot <- m_long %>%
  dplyr::group_by(stage) %>%
  dplyr::summarise(
    stage_x = dplyr::first(stage_x),
    line_y = annot_line_y,
    label_y = annot_label_y,
    .groups = "drop") %>%
  dplyr::left_join(
    stats_df %>%
      dplyr::filter(stage %in% stage_order) %>%
      transmute(
        stage = factor(stage, levels = stage_order),
        p_label = vapply(p_wilcox_spot_level, p_label_one, character(1))),
    by = "stage"); p_annot

head(m_long, n = 10)
p_violin <- ggplot2::ggplot(
  data = m_long, 
  mapping = aes(x = group, y = score, fill = group)) + 
  ggplot2::geom_violin(trim = FALSE, alpha = 0.6, width = 1, linewidth = 1) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.1, linewidth = 1) +
  ggplot2::scale_fill_manual(values = c("#5F7FA8", "#C97A74"), name = NULL) +
  ggpubr::stat_compare_means(
    comparisons = list(c("HC", "TN")), 
    method = "wilcox.test", 
    tip.length = 0.03, 
    bracket.size = 0.3,
    label.y = 1.1, 
    label = "p.format", # p.format
    size = 12) +
  ggplot2::facet_wrap(~ stage , nrow = 1) + 
  ggplot2::labs(
    x = "Differentiation stage (early -> late)",
    y = "UCell score",
    title = "OL Lineage 5 Stages: HC vs TN",
    subtitle = "Merged stage-level violin plot;\nP values are descriptive spot-level Wilcoxon tests",
    fill = NULL) + 
  ggplot2::coord_cartesian(ylim = c(-0.05, 1.2)) +
  ggplot2::scale_y_continuous(labels = scales::comma, breaks = seq(0, 1.2, 0.2)) + 
  ggplot2::theme_classic(base_size = 10) + 
  ggplot2::theme(
    panel.grid = element_blank(), 
    strip.background = element_rect(fill = "#FFFFFFFF", color = "#FFFFFFFF"), 
    strip.text = element_text(size = 35, face = "bold"),
    axis.title.x = element_text(size = 30, face = "bold"), 
    axis.title.y = element_text(size = 30, face = "bold"), 
    axis.text.x = element_text(size = 25, face = "bold"), 
    axis.text.y = element_text(size = 25, face = "bold"), 
    legend.text = element_text(size = 30, face = "bold"),
    legend.title = element_text(size = 35, face = "bold"), 
    panel.spacing.x = unit(0, "lines"), 
    axis.line.x = element_line(linewidth = 1, colour = "black"), 
    legend.position = "right",
    legend.key.height = unit(5, "lines"),
    legend.key.width  = unit(5, "lines"),
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 35, face = "bold", hjust = 0.5)); p_violin
ggsave(
  plot = p_violin, 
  filename = "results_OPC_stage/03_plots/p_violin.pdf", 
  width = 35, height = 14, dpi = 300, units = "in", limitsize = FALSE)

dom_comp <- m %>%
  dplyr::group_by(sample, group, dominant_stage) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::group_by(sample) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup(); dom_comp
write.csv(dom_comp, file.path(out_dir, "02_stats", "dominant_stage_composition.csv"), row.names = FALSE)

dom_group <- dom_comp %>%
  dplyr::group_by(group, dominant_stage) %>%
  dplyr::summarise(n = sum(n), .groups = "drop") %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    group = factor(group, levels = c("HC", "TN")),
    dominant_stage = factor(dominant_stage, levels = stage_order)); dom_group
write.csv(dom_group, file.path(out_dir, "02_stats", "dominant_stage_composition_group_pooled.csv"), row.names = FALSE)
head(dom_group)
p_traj <- ggplot2::ggplot(
  data = dom_group, 
  mapping = aes(
    x = group, 
    y = pct, 
    fill = dominant_stage, 
    stratum = dominant_stage, 
    alluvium = dominant_stage)) + 
  ggplot2::scale_fill_manual(values = c("#6B89AEFF", "#93A7BEFF", "#CEC8C1FF", "#D8A19AFF", "#C17970FF" )) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    y = "Spot-level % (dominant stage per spot)",
    x = "Group",
    title = "Spot-level stage classification",
    subtitle = "Pooled by group; spot-level descriptive") + 
  ggplot2::theme(
    panel.grid = element_blank(), 
    axis.text.y = element_blank(), 
    axis.title.x = element_text(size = 14, face = "bold"), 
    axis.title.y = element_text(size = 14, face = "bold"), 
    axis.text.x = element_text(size = 12, face = "bold"), 
    legend.text = element_text(size = 12, face = "bold"), 
    legend.title = element_text(size = 14, face = "bold"), 
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12, face = "bold"),
    legend.key.height = unit(1.5, "lines"),
    legend.key.width  = unit(1.5, "lines")) + 
  ggplot2::geom_col(width = 0.5, color = NA, linewidth = 1) + 
  ggalluvial::geom_flow(width = 0.5, alpha = 0.25, knot.pos = 0.3, color = "#FFFFFFFF", linewidth = 0.5) + 
  ggalluvial::geom_alluvium(width = 0.5, alpha = 1, knot.pos = 0.3, fill = NA, color = "#FFFFFFFF", linewidth = 0.5); p_traj
ggsave(
  plot = p_traj, 
  filename = "results_OPC_stage/03_plots/p_traj.pdf", 
  width = 8, height = 8, dpi = 300, units = "in", limitsize = FALSE)
saveRDS(obj, file.path(out_dir, "obj_staged.rds"))
