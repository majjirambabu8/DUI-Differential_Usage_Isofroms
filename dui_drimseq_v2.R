#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(argparse); library(DRIMSeq); library(data.table); library(GenomicFeatures)
  library(tximport); library(stageR); library(dplyr); library(reshape2); library(ggplot2)
  library(ggbeeswarm); library(rtracklayer); library(RColorBrewer); library(tidyr)
  library(BiocParallel); library(matrixStats)
})

# Utility functions
standardize_names <- function(names) make.names(names)
setup_dirs <- function(output_dir, contrast_name) { 
  dirs <- list(
    main = file.path(output_dir, contrast_name),
    plots = file.path(output_dir, contrast_name, "plots"),
    rds = file.path(output_dir, contrast_name, "rds"),
    tables = file.path(output_dir, contrast_name, "tables")
  )
  sapply(dirs, function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE))
  return(dirs)
}

setup_parallel <- function(threads = 1, verbose = TRUE) {
  if (threads <= 1) {
    if (verbose) cat("Using serial processing\n")
    return(BiocParallel::SerialParam())
  }
  
  max_cores <- parallel::detectCores()
  if (!is.na(max_cores) && threads > max_cores) threads <- max_cores
  if (verbose) cat("Setting up", threads, "threads\n")
  
  tryCatch({
    bp_param <- if (.Platform$OS.type == "unix") {
      BiocParallel::MulticoreParam(workers = threads, progressbar = verbose)
    } else {
      BiocParallel::SnowParam(workers = threads, progressbar = verbose, type = "SOCK")
    }
    BiocParallel::bplapply(1:2, function(x) x^2, BPPARAM = bp_param)
    BiocParallel::register(bp_param)
    return(bp_param)
  }, error = function(e) {
    warning("Parallel setup failed, using serial: ", e$message)
    return(BiocParallel::SerialParam())
  })
}

save_plot <- function(plot, path, width = 8, height = 6, dpi = 300) {
  ggsave(filename = path, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
}

# Build a robust transcript -> gene mapping from a GTF.
# Returns a data.frame with three columns: transcript_id, gene_id (Ensembl), gene_name (symbol).
# Strips Ensembl version suffixes so it works whether quant.sf uses ENSMUST00000193812
# or ENSMUST00000193812.2, since tximport is called with ignoreTxVersion=TRUE.
build_tx2gene <- function(gtf_path, verbose = TRUE) {
  if (verbose) cat("Building tx2gene from GTF:", gtf_path, "\n")
  gtf <- as.data.frame(rtracklayer::import(gtf_path, format = "gtf"))
  
  if (!"type" %in% colnames(gtf)) stop("GTF missing 'type' column from rtracklayer import")
  tx_rows <- gtf[gtf$type == "transcript", , drop = FALSE]
  if (nrow(tx_rows) == 0) {
    # Fall back to deriving from exons if no transcript-level features exist
    tx_rows <- gtf[gtf$type == "exon", , drop = FALSE]
  }
  
  if (!all(c("transcript_id", "gene_id") %in% colnames(tx_rows))) {
    stop("GTF must contain transcript_id and gene_id attributes")
  }
  
  tx2gene <- unique(data.frame(
    transcript_id = as.character(tx_rows$transcript_id),
    gene_id       = as.character(tx_rows$gene_id),
    gene_name     = if ("gene_name" %in% colnames(tx_rows)) as.character(tx_rows$gene_name)
                    else as.character(tx_rows$gene_id),
    stringsAsFactors = FALSE
  ))
  
  # Drop rows with missing transcript_id (shouldn't happen for valid GTF)
  tx2gene <- tx2gene[!is.na(tx2gene$transcript_id) & nzchar(tx2gene$transcript_id), , drop = FALSE]
  
  # Fill missing gene_name with gene_id to avoid NA groups downstream
  na_name <- is.na(tx2gene$gene_name) | !nzchar(tx2gene$gene_name)
  tx2gene$gene_name[na_name] <- tx2gene$gene_id[na_name]
  
  if (verbose) {
    cat("tx2gene rows:", nrow(tx2gene),
        "| unique transcripts:", length(unique(tx2gene$transcript_id)),
        "| unique genes:", length(unique(tx2gene$gene_id)), "\n")
  }
  return(tx2gene)
}

# Original stacked barplot function (unchanged)
create_stacked_barplot <- function(prop_matrix, gene_name, sample_metadata, output_path = NULL, 
                                 plot_width = 16, plot_height = 8, plot_dpi = 300) {
  gene_data <- filter(prop_matrix, gene_id == gene_name)
  if (nrow(gene_data) < 2) {
    message("Gene ", gene_name, " has fewer than 2 transcripts, skipping...")
    return(FALSE)
  }
  
  groups <- unique(sample_metadata$group)
  if (length(groups) != 2) stop("Need exactly 2 groups for comparison")
  
  # Calculate group averages
  group1_samples <- standardize_names(sample_metadata$sample_id[sample_metadata$group == groups[1]])
  group2_samples <- standardize_names(sample_metadata$sample_id[sample_metadata$group == groups[2]])
  
  existing_samples1 <- intersect(group1_samples, colnames(gene_data))
  existing_samples2 <- intersect(group2_samples, colnames(gene_data))
  
  if (length(existing_samples1) == 0 || length(existing_samples2) == 0) {
    warning("No samples found for one or both groups for gene ", gene_name)
    return(FALSE)
  }
  
  base_data <- gene_data %>%
    mutate(
      Control_strength = rowMeans(select(., all_of(existing_samples1)), na.rm = TRUE) * 100,
      Treatment_strength = rowMeans(select(., all_of(existing_samples2)), na.rm = TRUE) * 100
    ) %>%
    select(isoformID = feature_id, Control_strength, Treatment_strength) %>%
    arrange(desc(Control_strength))
  
  plot_data <- gather(base_data, key = "Condition", value = "Percentage", -isoformID)
  plot_data$isoformID_short <- substr(plot_data$isoformID, 1, 15)
  
  legend_data <- plot_data %>%
    spread(key = Condition, value = Percentage, fill = 0) %>%
    mutate(
      legend_label = paste0(isoformID_short, " (", groups[1], ": ", 
                           round(Control_strength, 1), "%, ", groups[2], ": ", 
                           round(Treatment_strength, 1), "%)")
    )
  
  legend_labels <- setNames(legend_data$legend_label, legend_data$isoformID)
  plot_data$isoformID <- factor(plot_data$isoformID, levels = base_data$isoformID)
  plot_data$Condition <- factor(plot_data$Condition, levels = c("Control_strength", "Treatment_strength"))
  
  # Colors
  n_isoforms <- length(unique(plot_data$isoformID))
  if (n_isoforms <= 8) {
    colors <- brewer.pal(max(n_isoforms, 3), "Dark2")
  } else if (n_isoforms <= 12) {
    colors <- brewer.pal(max(n_isoforms, 3), "Set1")
  } else {
    colors <- rainbow(n_isoforms, s = 0.8, v = 0.7)
  }
  names(colors) <- levels(plot_data$isoformID)
  
  p <- ggplot(plot_data, aes(x = Condition, y = Percentage, fill = isoformID)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6, 
             color = "white", size = 1.0) +
    labs(
      title = paste("Isoform Usage:", gene_name),
      subtitle = paste(groups[1], "vs", groups[2], "- DRIMSeq Analysis"),
      y = "Percentage (%)",
      x = "Condition"
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100), breaks = seq(0, 100, 25)) +
    scale_x_discrete(breaks = c("Control_strength", "Treatment_strength"), 
                    labels = groups) +
    scale_fill_manual(values = colors, name = "Isoform (Percentages)", 
                     labels = legend_labels) +
    coord_cartesian(ylim = c(0, 100)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "white", color = "black", size = 1),
      plot.background = element_rect(fill = "white"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", size = 0.5),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
      axis.title = element_text(size = 12, face = "bold", color = "black"),
      axis.text = element_text(size = 11, color = "black"),
      legend.position = "right",
      legend.text = element_text(size = 9, color = "black"),
      legend.title = element_text(size = 11, face = "bold", color = "black"),
      legend.background = element_rect(fill = "white", color = "grey80"),
      plot.margin = unit(c(10, 10, 10, 10), "points")
    ) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE))
  
  if (!is.null(output_path)) save_plot(p, output_path, plot_width, plot_height, plot_dpi)
  return(p)
}

# Load Salmon data using GTF-derived tx2gene for grouping.
# tx2gene must contain columns: transcript_id, gene_id (Ensembl), gene_name (symbol).
# Output gene_id column = Ensembl gene_id (used for DRIMSeq grouping).
# Output gene_name column = symbol (kept as annotation; written alongside).
load_salmon_data <- function(salmon_dir, metadata, sample_id_col, tx2gene) {
  cat("Loading Salmon data...\n")
  
  # Debug: show what we're working with
  cat("Sample name mapping (first 5):\n")
  for(i in 1:min(5, nrow(metadata))) {
    cat(sprintf("  Original: %s -> Standardized: %s\n", 
                metadata$sample_id_original[i], 
                metadata$sample_id[i]))
  }
  
  # Use ORIGINAL names for file paths (with hyphens)
  files <- file.path(salmon_dir, metadata$sample_id_original, "quant.sf")
  # Use STANDARDIZED names for the file names attribute (with dots for R)
  names(files) <- metadata$sample_id
  
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    cat("\n=== DEBUGGING INFO ===\n")
    cat("Missing files:\n")
    cat(paste(missing, collapse = "\n"), "\n")
    
    cat("\nExpected directory structure:\n")
    for(i in 1:min(5, nrow(metadata))) {
      expected_dir <- file.path(salmon_dir, metadata$sample_id_original[i])
      exists_status <- if(dir.exists(expected_dir)) "✓" else "✗"
      cat(sprintf("  %s %s\n", exists_status, expected_dir))
    }
    
    stop("Missing Salmon quantification files - check paths above")
  }
  
  # tximport at transcript level. ignoreTxVersion=TRUE strips ".N" suffixes from
  # both the quant.sf Name column and the tx2gene transcript_id column for matching,
  # so the index can have been built with or without versioned IDs.
  txi <- tximport(files, type = "salmon", txOut = TRUE,
                 tx2gene = tx2gene[, c("transcript_id", "gene_id", "gene_name")],
                 ignoreTxVersion = TRUE, countsFromAbundance = "no")
  
  counts <- round(txi$counts)
  colnames(counts) <- standardize_names(colnames(counts))
  
  # Build per-transcript annotation by joining quant rownames against tx2gene.
  # Match is on un-versioned transcript_id (same normalization tximport used).
  strip_version <- function(x) sub("\\..*$", "", x)
  quant_tx_raw      <- rownames(counts)
  quant_tx_unversioned <- strip_version(quant_tx_raw)
  tx2gene_unversioned  <- strip_version(tx2gene$transcript_id)
  
  match_idx <- match(quant_tx_unversioned, tx2gene_unversioned)
  tx_data <- data.frame(
    feature_id = quant_tx_raw,
    gene_id    = tx2gene$gene_id[match_idx],
    gene_name  = tx2gene$gene_name[match_idx],
    row.names  = quant_tx_raw,
    stringsAsFactors = FALSE
  )
  
  n_unmatched <- sum(is.na(tx_data$gene_id))
  if (n_unmatched > 0) {
    cat("Note:", n_unmatched, "of", nrow(tx_data),
        "quantified transcripts had no match in the GTF and will be dropped.\n")
  }
  
  # Drop transcripts that didn't map to any gene
  valid <- !is.na(tx_data$gene_id)
  counts <- counts[valid, , drop = FALSE]
  tx_data <- tx_data[valid, , drop = FALSE]
  
  if (nrow(counts) == 0) {
    stop("No transcripts could be mapped to genes via GTF tx2gene. ",
         "Check that the GTF matches the Salmon index transcriptome.")
  }
  
  # Keep only genes with >= 2 transcripts (DRIMSeq requirement for DTU)
  gene_counts <- table(tx_data$gene_id)
  multi_tx <- tx_data$gene_id %in% names(gene_counts)[gene_counts >= 2]
  counts <- counts[multi_tx, , drop = FALSE]
  tx_data <- tx_data[multi_tx, , drop = FALSE]
  
  cat("Loaded", nrow(counts), "transcripts from",
      length(unique(tx_data$gene_id)), "genes (multi-isoform only)\n")
  return(list(counts = counts, tx_data = tx_data))
}

# Proportion SD filter
smallProportionSD <- function(d, filter = 0.1) {
  cts <- as.matrix(subset(counts(d), select = -c(gene_id, feature_id)))
  gene.cts <- rowsum(cts, counts(d)$gene_id)
  total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
  props <- cts/total.cts
  sqrt(matrixStats::rowVars(props)) < filter
}

# Main DRIMSeq analysis
drimseq_analysis <- function(counts, design, contrast, output_dirs, contrast_name, args, samps, bp_param) {
  set.seed(args$seed)
  cat("Running DRIMSeq for:", contrast, "\n")
  
  samps$sample_id <- standardize_names(samps$sample_id)
  
  # Create and filter DRIMSeq object
  dds <- dmDSdata(counts = counts, samples = samps)
  n.small <- min(table(samps$group))
  
  d <- dmFilter(dds,
    min_samps_feature_expr = n.small, min_feature_expr = args$min_feature_expr,
    min_samps_feature_prop = n.small, min_feature_prop = args$min_feature_prop,
    min_samps_gene_expr = nrow(samps), min_gene_expr = args$min_gene_expr,
    run_gene_twice = FALSE
  )
  
  cat("Genes with multiple isoforms after filtering:\n")
  print(table(table(counts(d)$gene_id)))
  
  # DRIMSeq workflow
  cat("Estimating precision parameters...\n")
  d <- dmPrecision(d, design = design, mean_expression = TRUE, prec_adjust = TRUE,
                  prec_subset = args$prec_subset, prec_tol = args$prec_tol,
                  prec_init = args$prec_init, prec_span = args$prec_span,
                  one_way = TRUE, verbose = args$verbose, BPPARAM = bp_param)
  
  cat("Fitting regression models...\n")
  d <- dmFit(d, design = design, one_way = TRUE, bb_model = TRUE,
            prop_tol = args$prop_tol, coef_tol = args$coef_tol,
            verbose = args$verbose, add_uniform = FALSE, BPPARAM = bp_param)
  
  cat("Performing hypothesis testing...\n")
  ref_contrast <- paste0('group', contrast)
  d <- dmTest(d, coef = ref_contrast, BPPARAM = bp_param)
  
  proportion.matrix <- proportions(d)
  res.g <- DRIMSeq::results(d, level = "gene")
  res.t <- DRIMSeq::results(d, level = "feature")
  
  # Save core results
  cat("Saving results to:", output_dirs$tables, "\n")
  saveRDS(d, file.path(output_dirs$rds, paste0(contrast_name, ".rds")))
  write.table(counts, file.path(output_dirs$tables, paste0(contrast_name, "_counts.txt")), 
              sep = "\t", quote = FALSE, row.names = TRUE)
  write.table(proportion.matrix, file.path(output_dirs$tables, paste0(contrast_name, "_proportion_matrix.txt")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(res.g, file.path(output_dirs$tables, paste0(contrast_name, "_gene_results.txt")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(res.t, file.path(output_dirs$tables, paste0(contrast_name, "_transcript_results.txt")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  return(d)
}

# Merge results for final output. Adds gene_name (symbol) annotation alongside gene_id.
merge_results <- function(df_count, df_proportion, df_res, tx_annotation = NULL) {
  common_ids <- intersect(intersect(row.names(df_count), df_proportion$feature_id), df_res$feature_id)
  cat("Common feature_ids found:", length(common_ids), "\n")
  if (length(common_ids) == 0) stop("No common feature_ids found")
  
  df_count <- df_count[common_ids, ]
  df_proportion <- df_proportion[df_proportion$feature_id %in% common_ids, ]
  df_res <- df_res[df_res$feature_id %in% common_ids, ]
  
  samples <- setdiff(names(df_proportion), c("feature_id", "gene_id"))
  
  # Build a gene_id -> gene_name lookup from tx_annotation (if provided)
  gene_name_lookup <- NULL
  if (!is.null(tx_annotation) && all(c("gene_id", "gene_name") %in% colnames(tx_annotation))) {
    gn <- unique(tx_annotation[, c("gene_id", "gene_name")])
    # If a gene_id maps to multiple gene_names (rare), keep the first
    gn <- gn[!duplicated(gn$gene_id), , drop = FALSE]
    gene_name_lookup <- setNames(gn$gene_name, gn$gene_id)
  }
  
  header <- c('gene_id', 'gene_name', 'n_transcripts', 'transcript_ids', 'p_values', 'adj_pvalues',
              paste0(samples, '_count'), paste0(samples, '_proportion'))
  
  results <- lapply(unique(df_res$gene_id), function(g) {
    df_res_g <- df_res[df_res$gene_id == g, ][order(df_res[df_res$gene_id == g, ]$feature_id), ]
    df_count_g <- df_count[df_res_g$feature_id, ]
    df_prop_g <- df_proportion[df_proportion$gene_id == g, ]
    gname <- if (!is.null(gene_name_lookup) && g %in% names(gene_name_lookup)) gene_name_lookup[[g]] else NA
    
    c(g, gname, nrow(df_count_g), paste(rownames(df_count_g), collapse = ','),
      paste(df_res_g$pvalue, collapse = ','), paste(df_res_g$adj_pvalue, collapse = ','),
      sapply(samples, function(s) paste(df_count_g[[s]], collapse = ',')),
      sapply(samples, function(s) paste(df_prop_g[[s]], collapse = ',')))
  })
  
  final_df <- data.frame(do.call(rbind, results), stringsAsFactors = FALSE)
  names(final_df) <- header
  return(final_df)
}

# Generate individual gene plots
generate_gene_plots <- function(gene_id, d, df_proportion, plots_dir, contrast_name, args) {
  tryCatch({
    prop_plot_dir <- file.path(plots_dir, "proportion_plots")
    if (!dir.exists(prop_plot_dir)) dir.create(prop_plot_dir, recursive = TRUE)
    
    # Boxplot (original)
    plotProportions(d, gene_id = gene_id, group_variable = "group")
    save_plot(last_plot(), file.path(prop_plot_dir, paste0(gene_id, "_proportions.png")),
              args$plot_width, args$plot_height, args$plot_dpi)
    
    # Lineplot (original)
    plotProportions(d, gene_id = gene_id, group_variable = "group", plot_type = "lineplot")
    save_plot(last_plot(), file.path(prop_plot_dir, paste0(gene_id, "_proportions_lineplot.png")),
              args$plot_width, args$plot_height, args$plot_dpi)
    
    # Ribbonplot (original)
    plotProportions(d, gene_id = gene_id, group_variable = "group", plot_type = "ribbonplot")
    save_plot(last_plot(), file.path(prop_plot_dir, paste0(gene_id, "_proportions_ribbonplot.png")),
              args$plot_width, args$plot_height, args$plot_dpi)
    
    # Stacked barplot (original)
    create_stacked_barplot(
      prop_matrix = df_proportion,
      gene_name = gene_id,
      sample_metadata = d@samples,
      output_path = file.path(prop_plot_dir, paste0(gene_id, "_stacked_barplot.png")),
      plot_width = args$plot_width * 2,
      plot_height = args$plot_height,
      plot_dpi = args$plot_dpi
    )
    
    return(paste("Completed plots for", gene_id))
  }, error = function(e) {
    return(paste("Error plotting", gene_id, ":", e$message))
  })
}

# Generate all plots and results
generate_plots <- function(d, output_dirs, contrast_name, args, tx_annotation = NULL) {
  cat("Generating diagnostic plots...\n")
  
  # Generate original diagnostic plots with new naming
  plotData(d)
  save_plot(last_plot(), file.path(output_dirs$plots, paste0(contrast_name, ".png")), 
           args$plot_width, args$plot_height, args$plot_dpi)
  
  plotPrecision(d)
  save_plot(last_plot(), file.path(output_dirs$plots, paste0(contrast_name, "_precision.png")), 
           args$plot_width, args$plot_height, args$plot_dpi)
  
  plotPValues(d)
  save_plot(last_plot(), file.path(output_dirs$plots, paste0(contrast_name, "_pvalues.png")), 
           args$plot_width, args$plot_height, args$plot_dpi)
  
  res.g <- DRIMSeq::results(d, level = "gene")
  res.t <- DRIMSeq::results(d, level = "feature")
  
  cat("Transcripts with adjusted p-value < 0.05 before filtering:\n")
  print(table(res.t$adj_pvalue < 0.05))
  
  # Apply SD filter if specified
  if (args$prop_sd_filter > 0) {
    cat("Applying proportion SD filter with threshold", args$prop_sd_filter, "...\n")
    filt <- smallProportionSD(d, filter = args$prop_sd_filter)
    res.t$pvalue[filt] <- 1
    res.t$adj_pvalue[filt] <- 1
    cat("Transcripts with adjusted p-value < 0.05 after filtering:\n")
    print(table(res.t$adj_pvalue < 0.05))
  }
  
  res.t <- res.t[order(res.t$pvalue), ]
  
  # Save results
  write.table(res.t, file.path(output_dirs$tables, paste0(contrast_name, "_filtered.txt")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  df_final <- merge_results(as.data.frame(d@counts@unlistData), proportions(d), res.t,
                            tx_annotation = tx_annotation)
  write.table(df_final, file.path(output_dirs$tables, paste0(contrast_name, "_all.txt")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Generate proportion plots for top_k genes
  cat("Generating proportion plots for top", args$top_k, "genes...\n")
  genes_to_plot <- unique(res.t$gene_id[1:min(args$top_k, nrow(res.t))])
  
  plot_results <- lapply(genes_to_plot, 
                        FUN = generate_gene_plots,
                        d = d,
                        df_proportion = proportions(d),
                        plots_dir = output_dirs$plots,
                        contrast_name = contrast_name,
                        args = args)
  
  for (result in plot_results) cat(result, "\n")
  
  # Summary file
  summary_text <- c(
    paste("=== DRIMSeq Analysis Summary ==="),
    paste("Contrast:", contrast_name), paste("Analysis Date:", Sys.time()), "",
    paste("=== Data Overview ==="),
    paste("Total samples:", nrow(d@samples)),
    paste("Total genes analyzed:", length(unique(counts(d)$gene_id))),
    paste("Total transcripts analyzed:", nrow(counts(d))), "",
    paste("=== Filtering Parameters ==="),
    paste("Min feature expression:", args$min_feature_expr),
    paste("Min feature proportion:", args$min_feature_prop),
    paste("Min gene expression:", args$min_gene_expr),
    paste("Proportion SD filter:", args$prop_sd_filter), "",
    paste("=== Results Summary ==="),
    paste("Genes with p-value < 0.05:", sum(res.g$pvalue < 0.05, na.rm = TRUE)),
    paste("Genes with adjusted p-value < 0.05:", sum(res.g$adj_pvalue < 0.05, na.rm = TRUE)),
    paste("Transcripts with p-value < 0.05:", sum(res.t$pvalue < 0.05, na.rm = TRUE)),
    paste("Transcripts with adjusted p-value < 0.05:", sum(res.t$adj_pvalue < 0.05, na.rm = TRUE))
  )
  writeLines(summary_text, file.path(output_dirs$main, paste0(contrast_name, "_summary.txt")))
  
  return(list(res.g = res.g, res.t = res.t))
}

# Argument parser
create_parser <- function() {
  parser <- ArgumentParser(description = "DRIMSeq Differential Isoform Usage Analysis")
  
  # Core arguments
  parser$add_argument("-g", "--gtf", required = TRUE, help = "GTF annotation file")
  parser$add_argument("-s", "--salmon-dir", required = TRUE, help = "Salmon quantification directory")
  parser$add_argument("-m", "--metadata", required = TRUE, help = "Metadata file (tab-separated)")
  parser$add_argument("-o", "--output-dir", required = TRUE, help = "Output directory")
  
  # Sample configuration
  parser$add_argument("-i", "--sample-id-col", default = "sample", help = "Sample ID column [sample]")
  parser$add_argument("-c", "--contrast-col", default = "genotype", help = "Group column [genotype]")
  parser$add_argument("-r", "--control-contrast", default = "CT", help = "Control group [CT]")
  
  # Parameters
  parser$add_argument("--min-feature-expr", type = "integer", default = 1, help = "Min feature expression [1]")
  parser$add_argument("--min-feature-prop", type = "double", default = 0.1, help = "Min feature proportion [0.1]")
  parser$add_argument("--min-gene-expr", type = "integer", default = 1, help = "Min gene expression [1]")
  parser$add_argument("--prop-sd-filter", type = "double", default = 0.1, help = "Proportion SD filter [0.1]")
  parser$add_argument("--prec-subset", type = "double", default = 0.1, help = "Precision subset [0.1]")
  parser$add_argument("--prec-tol", type = "integer", default = 10, help = "Precision tolerance [10]")
  parser$add_argument("--prec-init", type = "integer", default = 100, help = "Initial precision [100]")
  parser$add_argument("--prec-span", type = "double", default = 0.1, help = "Precision span [0.1]")
  parser$add_argument("--prop-tol", type = "double", default = 1e-12, help = "Proportion tolerance [1e-12]")
  parser$add_argument("--coef-tol", type = "double", default = 1e-12, help = "Coefficient tolerance [1e-12]")
  
  # Visualization
  parser$add_argument("-k", "--top-k", type = "integer", default = 10, help = "Top genes to plot [10]")
  parser$add_argument("--plot-width", type = "double", default = 8, help = "Plot width [8]")
  parser$add_argument("--plot-height", type = "double", default = 6, help = "Plot height [6]")
  parser$add_argument("--plot-dpi", type = "integer", default = 300, help = "Plot DPI [300]")
  
  # Processing
  parser$add_argument("-t", "--threads", type = "integer", default = 1, help = "Threads [1]")
  parser$add_argument("-p", "--parallel-contrasts", action = "store_true", help = "Parallel contrasts")
  parser$add_argument("-e", "--seed", type = "integer", default = 12345, help = "Random seed [12345]")
  parser$add_argument("-v", "--verbose", type = "integer", default = 1, help = "Verbosity [1]")
  
  return(parser)
}

# Main function - FIXED VERSION (Option B: group by Ensembl gene_id, annotate with gene_name)
main <- function() {
  args <- create_parser()$parse_args()
  
  if (args$verbose > 0) {
    cat("=== DRIMSeq Analysis ===\n")
    for (name in names(args)) {
      cat(sprintf("%-20s: %s\n", name, args[[name]]))
    }
    cat("===================================\n\n")
  }
  
  # Validate inputs
  for (path in c(args$gtf, args$metadata, args$salmon_dir)) {
    if (!file.exists(path) && !dir.exists(path)) stop("Path not found: ", path)
  }
  if (!dir.exists(args$output_dir)) dir.create(args$output_dir, recursive = TRUE)
  
  # Load data - preserve original sample names IMMEDIATELY
  metadata_raw <- read.table(args$metadata, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Extract unique samples with their true original names (before any R processing)
  sample_col_index <- which(colnames(metadata_raw) == args$sample_id_col)
  group_col_index <- which(colnames(metadata_raw) == args$contrast_col)
  
  if (length(sample_col_index) == 0) stop("Sample column '", args$sample_id_col, "' not found")
  if (length(group_col_index) == 0) stop("Group column '", args$contrast_col, "' not found")
  
  # Get unique combinations of sample and group
  unique_samples <- unique(metadata_raw[, c(sample_col_index, group_col_index)])
  colnames(unique_samples) <- c("sample_id_original", "group")
  
  # Create the processed metadata
  metadata <- data.frame(
    sample_id_original = unique_samples$sample_id_original,  # TRUE originals (with hyphens)
    sample_id = standardize_names(unique_samples$sample_id_original),  # R-compatible (with dots)
    group = unique_samples$group,
    stringsAsFactors = FALSE
  )
  
  if (args$verbose > 0) {
    cat("Loaded metadata for", nrow(metadata), "samples\n")
    cat("Groups found:", paste(unique(metadata$group), collapse = ", "), "\n")
  }
  
  # Build tx2gene from GTF. Group by Ensembl gene_id; keep gene_name as annotation.
  tx2gene <- build_tx2gene(args$gtf, verbose = args$verbose > 0)
  
  bp_param <- setup_parallel(args$threads, args$verbose > 0)
  count_data <- load_salmon_data(args$salmon_dir, metadata, "sample_id", tx2gene)
  
  all_counts <- cbind(count_data$tx_data[, c("gene_id", "feature_id")], count_data$counts)
  # gene_id column = Ensembl gene_id (used by DRIMSeq for grouping)
  # gene_name annotation is preserved in count_data$tx_data and passed into merge_results below
  colnames(all_counts)[3:ncol(all_counts)] <- standardize_names(colnames(all_counts)[3:ncol(all_counts)])
  
  # Run analysis for each contrast
  contrasts <- setdiff(unique(metadata$group), args$control_contrast)
  if (length(contrasts) == 0) stop("No contrasts found")
  
  analyze_contrast <- function(contrast) {
    contrast_name <- paste0(args$control_contrast, "_vs_", contrast)
    if (args$verbose > 0) cat("\n--- Analyzing contrast:", contrast_name, "---\n")
    
    output_dirs <- setup_dirs(args$output_dir, contrast_name)
    
    samps <- metadata[metadata$group %in% c(args$control_contrast, contrast), ]
    samps$group <- relevel(factor(samps$group), ref = args$control_contrast)
    design <- model.matrix(~group, data = samps)
    counts <- all_counts[, c('gene_id', 'feature_id', samps$sample_id)]
    
    nested_bp <- if (args$parallel_contrasts && !inherits(bp_param, "SerialParam")) SerialParam() else bp_param
    d <- drimseq_analysis(counts, design, contrast, output_dirs, contrast_name, args, samps, nested_bp)
    generate_plots(d, output_dirs, contrast_name, args, tx_annotation = count_data$tx_data)
    
    if (args$verbose > 0) cat("Completed:", contrast_name, "\n")
    return(list(contrast = contrast, success = TRUE))
  }
  
  # Run contrasts
  if (args$parallel_contrasts && length(contrasts) > 1 && args$threads > 1) {
    bp_contrasts <- if (.Platform$OS.type == "unix") {
      MulticoreParam(workers = min(args$threads, length(contrasts)))
    } else {
      SnowParam(workers = min(args$threads, length(contrasts)))
    }
    results <- bplapply(contrasts, analyze_contrast, BPPARAM = bp_contrasts)
    if (inherits(bp_contrasts, "SnowParam")) try(bpstop(bp_contrasts), silent = TRUE)
  } else {
    results <- lapply(contrasts, analyze_contrast)
  }
  
  # Cleanup
  if (!inherits(bp_param, "SerialParam") && inherits(bp_param, "SnowParam")) {
    try(if (bpisup(bp_param)) bpstop(bp_param), silent = TRUE)
  }
  
  cat("\n=== Analysis Complete ===\n")
  cat("All results saved to:", args$output_dir, "\n")
  cat("Total contrasts analyzed:", length(contrasts), "\n")
  
  # Force cleanup but DON'T call quit here
  gc()  # Garbage collection
  
  # Return success status
  return(invisible(0))
}

# FIXED execution wrapper
if (!interactive()) {
  # Flush output before starting
  flush.console()
  
  # Capture any warnings to prevent broken pipe
  options(warn = 1)  # Print warnings as they occur
  
  # Run main and capture result
  result <- tryCatch(
    {
      suppressWarnings({
        main()
      })
      0  # Success
    },
    error = function(e) {
      message("Error: ", e$message)
      1  # Failure
    }
  )
  
  # Ensure all output is flushed
  flush.console()
  if (exists("sink.number") && sink.number() > 0) {
    sink()  # Close any open sinks
  }
  
  # Clean exit with proper status
  quit(save = "no", status = result, runLast = FALSE)
}