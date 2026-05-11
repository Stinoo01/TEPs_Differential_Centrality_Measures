# ============================================================
# functions.R
# ============================================================

# ----------------------------
# 1) Caricamento librerie
# ----------------------------
load_libraries <- function(pkgs) {
  suppressPackageStartupMessages({
    for (p in pkgs) {
      if (!requireNamespace(p, quietly = TRUE)) {
        stop(sprintf("Pacchetto mancante: %s. Installalo prima.", p))
      }
      library(p, character.only = TRUE)
    }
  })
  invisible(TRUE)
}

# Lista pacchetti (puoi aggiungere/togliere)
get_default_packages <- function() {
  c(
    "clusterProfiler", "org.Hs.eg.db", "stringr", "dplyr", "pathview",
    "igraph", "survival", "survminer", "TCGAbiolinks", "network",
    "DESeq2", "DT", "ggplot2", "NetworkToolbox", "psych", "GGally",
    "sna", "BiocGenerics", "SummarizedExperiment", "reshape2",
    "vegan", "edgeR", "ggrepel","ReactomePA"
  )
}

# ----------------------------
# 2) Pre-processing (NO I/O)
#    - sistemi nomi colonne
#    - selezioni patologia vs HC
#    - costruisci DESeq2 + DEGs + annot
# ----------------------------
preprocess_data <- function(counts, colData, genes_info,
                            condition = "GBM",
                            min_row_sum = 100,
                            padj_cutoff = 0.01,
                            lfc_cutoff = 0.5) {
  
  # sample_name pulito
  colData$sample_name <- sapply(colData$source_name_ch1, function(x) gsub("-", ".", x))
  
  # rimuovi "X" iniziale nei nomi colonne conte (tipico di read.table)
  colnames(counts) <- sapply(colnames(counts), function(x) sub("^X", "", x))
  
  # subset sample malati e sani
  condition.samples <- colData[colData$cancer.type.ch1 == condition, ]$sample_name
  healthy.samples   <- colData[colData$cancer.type.ch1 == "HC", ]$sample_name
  
  all.samples <- c(condition.samples, healthy.samples)
  
  # filtra counts
  counts <- counts[, intersect(all.samples, colnames(counts)), drop = FALSE]
  
  # colData allineato a counts
  colData_sub <- colData[colData$sample_name %in% colnames(counts), c("sample_name", "cancer.type.ch1")]
  colData_sub <- colData_sub[match(colnames(counts), colData_sub$sample_name), ]
  stopifnot(all(colData_sub$sample_name == colnames(counts)))
  
  group <- ifelse(colData_sub$cancer.type.ch1 == condition, "Cancer", "Healthy")
  group <- factor(group)
  group <- relevel(group, ref = "Healthy")
  
  colData_DESeq <- data.frame(row.names = colnames(counts), group = group)
  
  # DESeq2
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData_DESeq, design = ~ group)
  dds <- dds[rowSums(counts(dds)) > min_row_sum, ]
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("group", "Cancer", "Healthy"))
  res <- as.data.frame(res)
  res <- res[complete.cases(res), ]
  
  res_universe <- res
  res_sig <- res[res$padj < padj_cutoff, , drop = FALSE]
  
  res_up   <- res_sig[res_sig$log2FoldChange >  lfc_cutoff, , drop = FALSE]
  res_down <- res_sig[res_sig$log2FoldChange < -lfc_cutoff, , drop = FALSE]
  
  # annot (uniformo nomi)
  annot <- genes_info[, c("EnsemblGeneID", "Symbol", "GeneType")]
  colnames(annot) <- c("ENSEMBL", "SYMBOL", "TYPE")
  
  # normalizzazione
  vsd <- vst(dds, blind = TRUE)
  norm_mat <- assay(vsd)
  
  list(
    counts = counts,
    colData_sub = colData_sub,
    group = group,
    dds = dds,
    res = res,
    res_universe = res_universe,
    res_sig = res_sig,
    res_up = res_up,
    res_down = res_down,
    annot = annot,
    norm_mat = norm_mat
  )
}



# ----------------------------
# 3) Grafi + metriche + top genes
# ----------------------------
build_coexp_graph <- function(expr_mat, thr = 0.7, signed = FALSE) {
  gene_sd <- apply(expr_mat, 1, sd, na.rm = TRUE)
  expr_mat <- expr_mat[gene_sd > 0, , drop = FALSE]
  
  cor_mat <- cor(t(expr_mat), method = "spearman", use = "pairwise.complete.obs")
  diag(cor_mat) <- 0
  
  if (signed) {
    cor_mat[cor_mat < thr] <- 0
  } else {
    cor_mat[abs(cor_mat) < thr] <- 0
    cor_mat <- abs(cor_mat)  # <--- rende i pesi positivi
  }
  
  igraph::graph_from_adjacency_matrix(
    cor_mat, mode = "undirected", weighted = TRUE, diag = FALSE
  )
}


node_stats <- function(g) {
  data.frame(
    ENSEMBL = igraph::V(g)$name,
    degree = igraph::degree(g),
    
    # IGNORA I PESI (fondamentale quando ci sono pesi negativi)
    betweenness = igraph::betweenness(g, normalized = TRUE, weights = NA),
    
    # Anche questa può usare i pesi se presenti: meglio ignorarli
    eigen_centrality = igraph::eigen_centrality(g, weights = NA)$vector,
    
    stringsAsFactors = FALSE
  )
}


get_top_ensembl <- function(df, metric = c("degree", "betweenness"), top = 50) {
  metric <- match.arg(metric)
  delta_col <- paste0("delta_", metric)
  df2 <- df[!is.na(df[[delta_col]]), , drop = FALSE]
  df2 <- df2[order(abs(df2[[delta_col]]), decreasing = TRUE), , drop = FALSE]
  head(df2$ENSEMBL, top)
}





run_network_analysis <- function(norm_mat, group, degs_ensembl, annot,
                                 thr = 0.7, signed = TRUE, top_n = 50,
                                 select = c("intersection", "union", "degree", "betweenness")) {
  select <- match.arg(select)
  
  # --- solo DEGs presenti nella matrice ---
  degs_ensembl <- intersect(degs_ensembl, rownames(norm_mat))
  norm_counts_degs <- norm_mat[degs_ensembl, , drop = FALSE]
  
  # --- split per gruppo ---
  expr_cancer  <- norm_counts_degs[, group == "Cancer",  drop = FALSE]
  expr_healthy <- norm_counts_degs[, group == "Healthy", drop = FALSE]
  
  # --- grafi co-espressione ---
  g_cancer  <- build_coexp_graph(expr_cancer,  thr = thr, signed = signed)
  g_healthy <- build_coexp_graph(expr_healthy, thr = thr, signed = signed)
  
  # --- statistiche nodi ---
  stats_cancer  <- node_stats(g_cancer)
  stats_healthy <- node_stats(g_healthy)
  
  stats_compare <- merge(
    stats_cancer, stats_healthy,
    by = "ENSEMBL",
    suffixes = c("_cancer", "_healthy"),
    all = FALSE
  )
  
  stats_compare$delta_degree      <- stats_compare$degree_cancer - stats_compare$degree_healthy
  stats_compare$delta_betweenness <- stats_compare$betweenness_cancer - stats_compare$betweenness_healthy
  
  # --- aggiungi annotazione ---
  stats_annot <- merge(stats_compare, annot, by = "ENSEMBL", all.x = TRUE)
  stats_annot$RNA_class <- ifelse(stats_annot$TYPE == "protein-coding", "mRNA", "ncRNA")
  
  # --- top genes per metrica ---
  top_deg <- get_top_ensembl(stats_annot, metric = "degree", top = top_n)
  top_bet <- get_top_ensembl(stats_annot, metric = "betweenness", top = top_n)
  
  # --- selezione finale ---
  top_selected <- switch(
    select,
    intersection = intersect(top_deg, top_bet),
    union        = union(top_deg, top_bet),
    degree       = top_deg,
    betweenness  = top_bet
  )
  
  # --- tabella top (con SYMBOL) ---
  top_df <- data.frame(ENSEMBL = unique(top_selected), stringsAsFactors = FALSE)
  top_df <- merge(top_df, annot[, c("ENSEMBL", "SYMBOL")], by = "ENSEMBL", all.x = TRUE)
  
  list(
    g_cancer = g_cancer,
    g_healthy = g_healthy,
    stats_annot = stats_annot,
    top_ensembl = unique(top_selected),
    top_df = top_df,
    top_degree = top_deg,          # utile per debug / report
    top_betweenness = top_bet,     # utile per debug / report
    selection = select
  )
}




# ----------------------------
# 4) Plotting (curva rank)
#    ritorna lista di ggplot
# ----------------------------
make_delta_plots <- function(stats_annot,
                             n_labels = 25,
                             highlight_ensembl = NULL,
                             title_prefix = NULL) {
  
  # helper: cap simmetrico per asse y
  sym_cap <- function(x) {
    m <- max(abs(x), na.rm = TRUE)
    c(-m, m)
  }
  
  df_plot <- stats_annot %>%
    dplyr::mutate(
      ENSEMBL_clean = sub("\\..*$", "", ENSEMBL),
      gene_label = dplyr::if_else(!is.na(SYMBOL) & SYMBOL != "", SYMBOL, ENSEMBL_clean),
      abs_delta_degree = abs(delta_degree),
      abs_delta_betweenness = abs(delta_betweenness),
      RNA_class = dplyr::if_else(is.na(RNA_class), "unknown", RNA_class),
      sign_degree = dplyr::if_else(delta_degree >= 0, "↑ Cancer", "↑ Healthy"),
      sign_bet    = dplyr::if_else(delta_betweenness >= 0, "↑ Cancer", "↑ Healthy"),
      is_highlight = if (!is.null(highlight_ensembl)) ENSEMBL %in% highlight_ensembl else FALSE
    ) %>%
    dplyr::filter(!is.na(delta_degree), !is.na(delta_betweenness)) %>%
    dplyr::distinct(ENSEMBL_clean, .keep_all = TRUE)
  
  # ranking
  df_deg <- df_plot %>% dplyr::arrange(dplyr::desc(abs_delta_degree)) %>% dplyr::mutate(rank = dplyr::row_number())
  df_bet <- df_plot %>% dplyr::arrange(dplyr::desc(abs_delta_betweenness)) %>% dplyr::mutate(rank = dplyr::row_number())
  
  # common theme
  base_theme <- ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 11),
      legend.position = "right"
    )
  
  # label data (top N by abs metric; include highlights if provided)
  lab_deg <- df_deg %>% dplyr::slice_head(n = n_labels)
  lab_bet <- df_bet %>% dplyr::slice_head(n = n_labels)
  
  if (!is.null(highlight_ensembl)) {
    lab_deg <- dplyr::bind_rows(lab_deg, df_deg %>% dplyr::filter(is_highlight)) %>% dplyr::distinct(ENSEMBL, .keep_all = TRUE)
    lab_bet <- dplyr::bind_rows(lab_bet, df_bet %>% dplyr::filter(is_highlight)) %>% dplyr::distinct(ENSEMBL, .keep_all = TRUE)
  }
  
  # ----------------------------
  # Δdegree (signed)
  # ----------------------------
  ttl1 <- paste(c(title_prefix, "Δ degree (Cancer − Healthy)"), collapse = " — ")
  
  p_deg_signed <- ggplot2::ggplot(df_deg, ggplot2::aes(rank, delta_degree)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed", alpha = 0.7) +
    ggplot2::geom_line(linewidth = 0.5, alpha = 0.8) +
    ggplot2::geom_point(
      ggplot2::aes(color = sign_degree, shape = RNA_class, alpha = is_highlight),
      size = 2
    ) +
    ggrepel::geom_label_repel(
      data = lab_deg,
      ggplot2::aes(label = gene_label),
      size = 3,
      label.size = 0.2,
      label.r = grid::unit(0.12, "lines"),
      min.segment.length = 0,
      max.overlaps = Inf
    ) +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.55, `TRUE` = 1), guide = "none") +
    ggplot2::coord_cartesian(ylim = sym_cap(df_deg$delta_degree)) +
    ggplot2::labs(
      title = ttl1,
      subtitle = sprintf("Ranked by |Δdegree| (top labels = %d)", n_labels),
      x = "Rank (descending |Δdegree|)",
      y = "Δ degree",
      color = "Direction"
    ) +
    base_theme
  
  # ----------------------------
  # |Δdegree| (absolute)
  # ----------------------------
  ttl2 <- paste(c(title_prefix, "|Δ degree|"), collapse = " — ")
  
  p_deg_abs <- ggplot2::ggplot(df_deg, ggplot2::aes(rank, abs_delta_degree)) +
    ggplot2::geom_line(linewidth = 0.6, alpha = 0.85) +
    ggplot2::geom_point(
      ggplot2::aes(shape = RNA_class, alpha = is_highlight),
      size = 2
    ) +
    ggrepel::geom_label_repel(
      data = lab_deg,
      ggplot2::aes(label = gene_label),
      size = 3,
      label.size = 0.2,
      label.r = grid::unit(0.12, "lines"),
      min.segment.length = 0,
      max.overlaps = Inf
    ) +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.55, `TRUE` = 1), guide = "none") +
    ggplot2::labs(
      title = ttl2,
      subtitle = sprintf("Ranked by |Δdegree| (top labels = %d)", n_labels),
      x = "Rank (descending |Δdegree|)",
      y = "|Δ degree|"
    ) +
    base_theme
  
  # ----------------------------
  # Δbetweenness (signed)
  # ----------------------------
  ttl3 <- paste(c(title_prefix, "Δ betweenness (Cancer − Healthy)"), collapse = " — ")
  
  p_bet_signed <- ggplot2::ggplot(df_bet, ggplot2::aes(rank, delta_betweenness)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed", alpha = 0.7) +
    ggplot2::geom_line(linewidth = 0.5, alpha = 0.8) +
    ggplot2::geom_point(
      ggplot2::aes(color = sign_bet, shape = RNA_class, alpha = is_highlight),
      size = 2
    ) +
    ggrepel::geom_label_repel(
      data = lab_bet,
      ggplot2::aes(label = gene_label),
      size = 3,
      label.size = 0.2,
      label.r = grid::unit(0.12, "lines"),
      min.segment.length = 0,
      max.overlaps = Inf
    ) +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.55, `TRUE` = 1), guide = "none") +
    ggplot2::coord_cartesian(ylim = sym_cap(df_bet$delta_betweenness)) +
    ggplot2::labs(
      title = ttl3,
      subtitle = sprintf("Ranked by |Δbetweenness| (top labels = %d)", n_labels),
      x = "Rank (descending |Δbetweenness|)",
      y = "Δ betweenness",
      color = "Direction"
    ) +
    base_theme
  
  # ----------------------------
  # |Δbetweenness| (absolute)
  # ----------------------------
  ttl4 <- paste(c(title_prefix, "|Δ betweenness|"), collapse = " — ")
  
  p_bet_abs <- ggplot2::ggplot(df_bet, ggplot2::aes(rank, abs_delta_betweenness)) +
    ggplot2::geom_line(linewidth = 0.6, alpha = 0.85) +
    ggplot2::geom_point(
      ggplot2::aes(shape = RNA_class, alpha = is_highlight),
      size = 2
    ) +
    ggrepel::geom_label_repel(
      data = lab_bet,
      ggplot2::aes(label = gene_label),
      size = 3,
      label.size = 0.2,
      label.r = grid::unit(0.12, "lines"),
      min.segment.length = 0,
      max.overlaps = Inf
    ) +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.55, `TRUE` = 1), guide = "none") +
    ggplot2::labs(
      title = ttl4,
      subtitle = sprintf("Ranked by |Δbetweenness| (top labels = %d)", n_labels),
      x = "Rank (descending |Δbetweenness|)",
      y = "|Δ betweenness|"
    ) +
    base_theme
  
  list(
    p_deg_signed = p_deg_signed,
    p_deg_abs = p_deg_abs,
    p_bet_signed = p_bet_signed,
    p_bet_abs = p_bet_abs
  )
}





# ----------------------------
# 5) Enrichment (GO BP/MF/CC)
# ----------------------------
run_enrichment_all <- function(gene_ens, universe_ens,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.1,
                               pAdjustMethod = "BH",
                               OrgDb = org.Hs.eg.db,
                               readable = TRUE) {
  
  # ENSEMBL -> ENTREZ
  gene_entrez <- bitr(gene_ens, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
  universe_entrez <- bitr(universe_ens, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID
  
  gene_entrez <- unique(na.omit(gene_entrez))
  universe_entrez <- unique(na.omit(universe_entrez))
  
  # ---- GO ----
  go_fun <- function(ont) {
    obj <- enrichGO(
      gene = gene_ens,
      universe = universe_ens,
      OrgDb = OrgDb,
      keyType = "ENSEMBL",
      ont = ont,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      readable = readable
    )
    list(
      object = obj,
      table  = as.data.frame(obj)
    )
  }
  
  # ---- KEGG ----
  kegg_obj <- enrichKEGG(
    gene = gene_entrez,
    universe = universe_entrez,
    organism = "hsa",
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod
  )
  
  # ---- Reactome ----
  react_obj <- enrichPathway(
    gene = gene_entrez,
    universe = universe_entrez,
    organism = "human",
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    readable = TRUE
  )
  
  list(
    GO_BP = go_fun("BP"),
    GO_MF = go_fun("MF"),
    GO_CC = go_fun("CC"),
    KEGG = list(object = kegg_obj, table = as.data.frame(kegg_obj)),
    Reactome = list(object = react_obj, table = as.data.frame(react_obj))
  )
}











get_top_genes <- function(stats_df,
                          metric = c("degree", "betweenness"),
                          n = 100,
                          mode = c("abs", "up", "down")) {
  
  metric <- match.arg(metric)
  mode   <- match.arg(mode)
  
  delta_col <- paste0("delta_", metric)
  
  df <- stats_df %>%
    dplyr::filter(!is.na(.data[[delta_col]]))
  
  df <- switch(
    mode,
    abs  = df %>% dplyr::arrange(dplyr::desc(abs(.data[[delta_col]]))),
    up   = df %>% dplyr::arrange(dplyr::desc(.data[[delta_col]])),
    down = df %>% dplyr::arrange(.data[[delta_col]])
  )
  
  df %>%
    dplyr::slice_head(n = n) %>%
    dplyr::pull(ENSEMBL) %>%
    unique()
}


run_intersection_enrichment <- function(stats_df,
                                        n = 100,
                                        mode = c("abs", "up", "down"),
                                        universe_ens,
                                        OrgDb = org.Hs.eg.db) {
  
  mode <- match.arg(mode)
  
  top_deg <- get_top_genes(stats_df, metric = "degree", n = n, mode = mode)
  top_bet <- get_top_genes(stats_df, metric = "betweenness", n = n, mode = mode)
  
  genes_int <- intersect(top_deg, top_bet)
  
  enrich <- run_enrichment_all(
    gene_ens = genes_int,
    universe_ens = universe_ens,
    OrgDb = OrgDb
  )
  
  list(
    mode = mode,
    n = n,
    top_degree = top_deg,
    top_betweenness = top_bet,
    intersection = genes_int,
    intersection_size = length(genes_int),
    
    GO_BP = enrich$GO_BP,
    GO_MF = enrich$GO_MF,
    GO_CC = enrich$GO_CC,
    KEGG = enrich$KEGG,
    Reactome = enrich$Reactome
  )
}




run_permanova_permdisp_sets <- function(
    dds = NULL,
    vsd = NULL,
    group,
    degs_use,
    genes_network,
    vst_blind = TRUE,
    dist_method = "euclidean",
    n_iter = 100,
    n_perm_disp = 999,
    seed = 123
) {
  stopifnot(!is.null(vsd) || !is.null(dds))
  stopifnot(length(group) >= 2)
  
  if (!is.null(dds) && is.null(vsd)) {
    vsd <- DESeq2::vst(dds, blind = vst_blind)
  }
  
  mat_all <- SummarizedExperiment::assay(vsd)
  
  # assicuro che i set siano presenti nella matrice
  degs_use      <- intersect(degs_use, rownames(mat_all))
  genes_network <- intersect(genes_network, rownames(mat_all))
  
  prep_mat <- function(gene_vec, mat_all) {
    gene_vec <- intersect(gene_vec, rownames(mat_all))
    stopifnot(length(gene_vec) >= 2)
    
    m <- mat_all[gene_vec, , drop = FALSE]
    
    keep <- apply(m, 1, sd, na.rm = TRUE) > 0
    m <- m[keep, , drop = FALSE]
    stopifnot(nrow(m) >= 2)
    
    m <- t(scale(t(m)))
    
    m <- as.matrix(m)
    storage.mode(m) <- "double"
    m <- m[complete.cases(m), , drop = FALSE]
    m <- m[is.finite(rowSums(m)), , drop = FALSE]
    stopifnot(nrow(m) >= 2)
    
    m
  }
  
  get_permanova_stats <- function(mat, group, gene_set_name) {
    d <- dist(t(mat), method = dist_method)
    ad <- vegan::adonis2(d ~ group)
    
    data.frame(
      Gene_set = gene_set_name,
      N_genes  = nrow(mat),
      R2       = ad$R2[1],
      F_stat   = ad$F[1],
      p_value  = ad$`Pr(>F)`[1]
    )
  }
  
  get_permdisp_stats <- function(mat, group, gene_set_name) {
    d <- dist(t(mat), method = dist_method)
    bd <- vegan::betadisper(d, group)
    pt <- vegan::permutest(bd, permutations = n_perm_disp)
    
    data.frame(
      Gene_set = gene_set_name,
      N_genes  = nrow(mat),
      F_stat   = unname(pt$tab[1, "F"]),
      p_value  = unname(pt$tab[1, "Pr(>F)"])
    )
  }
  
  permanova_random <- function(n_genes, mat_all, group, n_iter, exclude_genes = NULL) {
    set.seed(seed)
    
    pool <- rownames(mat_all)
    if (!is.null(exclude_genes)) {
      exclude_genes <- intersect(exclude_genes, pool)
      pool <- setdiff(pool, exclude_genes)
    }
    stopifnot(length(pool) >= n_genes)
    
    res <- replicate(n_iter, {
      genes_rand <- sample(pool, n_genes)
      mat_rand <- prep_mat(genes_rand, mat_all)
      
      ad <- vegan::adonis2(dist(t(mat_rand), method = dist_method) ~ group)
      c(R2 = ad$R2[1], F = ad$F[1], p = ad$`Pr(>F)`[1])
    })
    t(res)
  }
  
  permdisp_random <- function(n_genes, mat_all, group, n_iter, exclude_genes = NULL) {
    set.seed(seed)
    
    pool <- rownames(mat_all)
    if (!is.null(exclude_genes)) {
      exclude_genes <- intersect(exclude_genes, pool)
      pool <- setdiff(pool, exclude_genes)
    }
    stopifnot(length(pool) >= n_genes)
    
    res <- replicate(n_iter, {
      genes_rand <- sample(pool, n_genes)
      mat_rand <- prep_mat(genes_rand, mat_all)
      
      d <- dist(t(mat_rand), method = dist_method)
      bd <- vegan::betadisper(d, group)
      pt <- vegan::permutest(bd, permutations = n_perm_disp)
      
      c(F = unname(pt$tab[1, "F"]), p = unname(pt$tab[1, "Pr(>F)"]))
    })
    t(res)
  }
  
  # ---- prepara matrici ----
  mat_degs <- prep_mat(degs_use, mat_all)
  mat_net  <- prep_mat(genes_network, mat_all)
  
  # ---- PERMANOVA ----
  perm_degs <- get_permanova_stats(mat_degs, group, "DEGs")
  perm_net  <- get_permanova_stats(mat_net, group, sprintf("Network genes (%d)", nrow(mat_net)))
  
  # ESCLUDO i genes_network dai random
  rand_perm <- permanova_random(
    n_genes = nrow(mat_net),
    mat_all = mat_all,
    group = group,
    n_iter = n_iter,
    exclude_genes = genes_network
  )
  
  perm_rand <- data.frame(
    Gene_set = sprintf("Random genes (%d, mean ± sd)", nrow(mat_net)),
    N_genes  = nrow(mat_net),
    R2       = paste0(round(mean(rand_perm[, "R2"]), 3), " ± ", round(sd(rand_perm[, "R2"]), 3)),
    F_stat   = paste0(round(mean(rand_perm[, "F"]),  2), " ± ", round(sd(rand_perm[, "F"]),  2)),
    p_value  = format.pval(median(rand_perm[, "p"]), eps = 1e-3)
  )
  
  perm_table <- rbind(perm_degs, perm_net, perm_rand)
  
  # ---- PERMDISP ----
  disp_degs <- get_permdisp_stats(mat_degs, group, "DEGs")
  disp_net  <- get_permdisp_stats(mat_net, group, sprintf("Network genes (%d)", nrow(mat_net)))
  
  # ESCLUDO i genes_network dai random
  rand_disp <- permdisp_random(
    n_genes = nrow(mat_net),
    mat_all = mat_all,
    group = group,
    n_iter = n_iter,
    exclude_genes = genes_network
  )
  
  disp_rand <- data.frame(
    Gene_set = sprintf("Random genes (%d, mean ± sd)", nrow(mat_net)),
    N_genes  = nrow(mat_net),
    F_stat   = paste0(round(mean(rand_disp[, "F"]), 2), " ± ", round(sd(rand_disp[, "F"]), 2)),
    p_value  = format.pval(median(rand_disp[, "p"]), eps = 1e-3)
  )
  
  disp_table <- rbind(disp_degs, disp_net, disp_rand)
  
  # ---- tabella finale unita (solo DEGs + Net; random resta separata perché ha stringhe) ----
  final_table <- merge(
    perm_table[perm_table$Gene_set != perm_rand$Gene_set, ],
    disp_table[disp_table$Gene_set != disp_rand$Gene_set, ],
    by = c("Gene_set", "N_genes"),
    suffixes = c("_PERMANOVA", "_PERMDISP")
  )
  
  list(
    vsd = vsd,
    mat_degs = mat_degs,
    mat_net  = mat_net,
    permanova_table = perm_table,
    permdisp_table  = disp_table,
    final_table     = final_table,
    random_permanova = rand_perm,
    random_permdisp  = rand_disp
  )
}





make_geneset_heatmap <- function(vsd = NULL,
                                 expr_mat = NULL,
                                 top_df,          # top_all_df (ENSEMBL + SYMBOL)
                                 group,
                                 main = NULL,
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,                 # <-- pazienti con dendrogramma
                                 cluster_cols_within_group = FALSE,   # <-- opzionale
                                 scale_rows = TRUE,
                                 show_colnames = FALSE,
                                 fontsize_row = 7,
                                 min_genes = 5) {
  
  stopifnot(!is.null(expr_mat) || !is.null(vsd))
  stopifnot(!is.null(top_df), !is.null(group))
  
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Pacchetto mancante: pheatmap")
  }
  
  if (is.null(expr_mat)) {
    expr_mat <- SummarizedExperiment::assay(vsd)
  }
  
  # --- geni e label dal top_df ---
  genes_ens <- unique(top_df$ENSEMBL)
  label_map <- top_df$SYMBOL
  label_map[label_map == "" | is.na(label_map)] <- top_df$ENSEMBL[label_map == "" | is.na(label_map)]
  names(label_map) <- top_df$ENSEMBL
  
  genes_present <- intersect(genes_ens, rownames(expr_mat))
  if (length(genes_present) < min_genes) {
    warning("Troppi pochi geni per la heatmap.")
    return(invisible(NULL))
  }
  
  expr_sub <- expr_mat[genes_present, , drop = FALSE]
  labels <- label_map[rownames(expr_sub)]
  
  # --- collassa duplicati di label (media) ---
  expr_sum <- rowsum(expr_sub, group = labels, reorder = FALSE)
  label_counts <- as.numeric(table(labels)[rownames(expr_sum)])
  expr_lab <- expr_sum / label_counts
  
  # --- annotazione colonne ---
  ann <- data.frame(group = group)
  rownames(ann) <- colnames(expr_lab)
  
  # --- ordine colonne (se vuoi raggruppare Cancer/Healthy) ---
  gaps_col <- NULL
  
  if (cluster_cols_within_group) {
    # ordina per gruppo, ma poi clusterizza dentro gruppo passando cluster_cols=TRUE
    ord <- order(group)
    expr_lab <- expr_lab[, ord, drop = FALSE]
    ann <- ann[ord, , drop = FALSE]
    
    # separatore visivo tra gruppi
    rle_g <- rle(as.character(group[ord]))
    gaps_col <- cumsum(rle_g$lengths)[-length(rle_g$lengths)]
    
    # cluster pazienti entro gruppo:
    # usiamo un hclust globale ma con ordine iniziale per gruppo (pheatmap lo ricalcola comunque),
    # quindi l'effetto più robusto è: fare due hclust separati e incollarli:
    make_block_hclust <- function(mat, grp) {
      ord2 <- integer(0)
      for (lvl in unique(grp)) {
        idx <- which(grp == lvl)
        if (length(idx) >= 2) {
          hc <- hclust(dist(t(mat[, idx, drop = FALSE])))
          ord2 <- c(ord2, idx[hc$order])
        } else {
          ord2 <- c(ord2, idx)
        }
      }
      ord2
    }
    
    ord2 <- make_block_hclust(expr_lab, ann$group)
    expr_lab <- expr_lab[, ord2, drop = FALSE]
    ann <- ann[ord2, , drop = FALSE]
    
    # aggiorna gaps dopo riordino
    rle_g2 <- rle(as.character(ann$group))
    gaps_col <- cumsum(rle_g2$lengths)[-length(rle_g2$lengths)]
    
    cluster_cols <- FALSE  # già ordinati con clustering "within-group"
  } else {
    # se NON fai within-group:
    # se cluster_cols=TRUE lasciamo clustering libero,
    # se cluster_cols=FALSE mettiamo ordine per gruppo e gap
    if (!cluster_cols) {
      ord <- order(group)
      expr_lab <- expr_lab[, ord, drop = FALSE]
      ann <- ann[ord, , drop = FALSE]
      
      rle_g <- rle(as.character(group[ord]))
      gaps_col <- cumsum(rle_g$lengths)[-length(rle_g$lengths)]
    }
  }
  
  # --- scaling per riga ---
  mat <- expr_lab
  if (scale_rows) {
    mat <- t(scale(t(mat)))
    mat <- mat[complete.cases(mat), , drop = FALSE]
  }
  
  if (is.null(main)) {
    main <- sprintf("Expression heatmap — network genes (n=%d)", nrow(mat))
  }
  
  pheatmap::pheatmap(
    mat,
    annotation_col = ann,
    show_colnames = show_colnames,
    fontsize_row = fontsize_row,
    main = main,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,   # <-- dendrogramma pazienti se TRUE
    gaps_col = gaps_col,
    legend = TRUE
  )
}



# ============================================================
# SALVATAGGIO NETWORK IN PDF (6 PDF: cancer/healthy x 3 viste)
# - tuoi geni selezionati: ROSSO (più grandi)
# - altri geni: GRIGIO (più piccoli)
# - label: SYMBOL se disponibile, altrimenti ENSEMBL (solo selected)
# - layout "fit to page" per ridurre spazio bianco
# ============================================================
save_network_views_pdf <- function(
    g,
    selected_genes,
    prefix,
    outdir = ".",
    top_all_df,
    seed = 123,
    pdf_width = 11,
    pdf_height = 8.5,
    v_size_selected = 7,
    v_size_other = 3.5,
    edge_width = 0.35,
    label_cex = 0.6,
    label_only_selected = TRUE,
    layout_each_view = TRUE,
    fr_niter = 2000
) {
  stopifnot(igraph::is_igraph(g))
  if (is.null(igraph::V(g)$name)) stop("Il grafo deve avere V(g)$name (ENSEMBL).")
  
  if (!all(c("ENSEMBL", "SYMBOL") %in% colnames(top_all_df))) {
    stop("top_all_df deve contenere le colonne: ENSEMBL, SYMBOL")
  }
  
  # ---- mapping ENSEMBL -> SYMBOL ----
  df_map <- top_all_df[, c("ENSEMBL", "SYMBOL")]
  df_map <- df_map[!is.na(df_map$ENSEMBL) & df_map$ENSEMBL != "", , drop = FALSE]
  df_map <- df_map[order(df_map$ENSEMBL, is.na(df_map$SYMBOL)), , drop = FALSE]
  df_map <- df_map[!duplicated(df_map$ENSEMBL), , drop = FALSE]
  symbol_map <- setNames(df_map$SYMBOL, df_map$ENSEMBL)
  
  # ---- selected presenti nel grafo ----
  selected_genes <- unique(intersect(selected_genes, igraph::V(g)$name))
  if (length(selected_genes) == 0) stop("Nessun gene selezionato presente nel grafo.")
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # ---- layout helper ----
  fit_layout <- function(g_for_layout) {
    set.seed(seed)
    lay <- igraph::layout_with_fr(g_for_layout, niter = fr_niter)
    igraph::norm_coords(lay, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
  }
  
  lay_full_ref <- fit_layout(g)
  
  plot_pdf <- function(g_plot, layout_mat, file, selected_in_plot, title_txt) {
    all_names <- igraph::V(g_plot)$name
    is_sel <- all_names %in% selected_in_plot
    
    v_col <- ifelse(is_sel, "red", "grey70")
    v_size <- ifelse(is_sel, v_size_selected, v_size_other)
    e_col <- grDevices::adjustcolor("grey80", alpha.f = 0.35)
    
    labels_mapped <- unname(symbol_map[all_names])
    labels_final <- ifelse(is.na(labels_mapped) | labels_mapped == "",
                           all_names,
                           labels_mapped)
    v_lab <- if (label_only_selected) ifelse(is_sel, labels_final, NA) else labels_final
    
    grDevices::pdf(file = file, width = pdf_width, height = pdf_height, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)
    
    par(mar = c(0, 0, 2, 0))
    igraph::plot.igraph(
      g_plot,
      layout = layout_mat,
      rescale = FALSE,
      vertex.color = v_col,
      vertex.frame.color = NA,
      vertex.size = v_size,
      vertex.label = v_lab,
      vertex.label.cex = label_cex,
      vertex.label.color = "black",
      edge.color = e_col,
      edge.width = edge_width,
      main = title_txt
    )
    invisible(TRUE)
  }
  
  # -------- 1) FULL NETWORK --------
  f1 <- file.path(outdir, sprintf("%s_full_network.pdf", prefix))
  lay1 <- if (layout_each_view) fit_layout(g) else lay_full_ref
  plot_pdf(g, lay1, f1, selected_genes, paste0(prefix, " - full network"))
  
  # -------- 2) SELECTED ONLY --------
  g_sel <- igraph::induced_subgraph(g, vids = selected_genes)
  f2 <- file.path(outdir, sprintf("%s_selected_only.pdf", prefix))
  
  lay2 <- if (layout_each_view) {
    fit_layout(g_sel)
  } else {
    igraph::norm_coords(
      lay_full_ref[match(igraph::V(g_sel)$name, igraph::V(g)$name), , drop = FALSE],
      xmin = -1, xmax = 1, ymin = -1, ymax = 1
    )
  }
  
  plot_pdf(g_sel, lay2, f2, selected_genes, paste0(prefix, " - selected only"))
  
  # -------- 3) SELECTED + 1st NEIGHBORS --------
  neigh_list <- igraph::neighborhood(g, order = 1, nodes = selected_genes)
  neigh_nodes <- unique(unlist(neigh_list))
  
  n_selected <- length(selected_genes)
  n_total_1st <- length(neigh_nodes)
  n_neighbors_only <- n_total_1st - n_selected
  
  message(
    sprintf(
      "[%s] Selected + 1st neighbors -> selected: %d | 1st neighbors: %d | total nodes: %d",
      prefix,
      n_selected,
      n_neighbors_only,
      n_total_1st
    )
  )
  
  g_1st <- igraph::induced_subgraph(g, vids = neigh_nodes)
  f3 <- file.path(outdir, sprintf("%s_selected_plus_1st_neighbors.pdf", prefix))
  
  lay3 <- if (layout_each_view) {
    fit_layout(g_1st)
  } else {
    igraph::norm_coords(
      lay_full_ref[match(igraph::V(g_1st)$name, igraph::V(g)$name), , drop = FALSE],
      xmin = -1, xmax = 1, ymin = -1, ymax = 1
    )
  }
  
  plot_pdf(g_1st, lay3, f3, selected_genes,
           paste0(prefix, " - selected + 1st neighbors"))
  
  invisible(list(full = f1, selected_only = f2, selected_plus_1st = f3))
}


plot_gene_boxplot <- function(gene,
                              norm_mat,
                              colData_sub,
                              annot,
                              group_col = "cancer.type.ch1",
                              condition_label = "Cancer") {
  
  # --- ENSEMBL -> SYMBOL se serve ---
  if (gene %in% annot$SYMBOL) {
    gene_ens <- annot$ENSEMBL[annot$SYMBOL == gene][1]
  } else {
    gene_ens <- gene
  }
  
  stopifnot(gene_ens %in% rownames(norm_mat))
  
  # --- long format ---
  df <- data.frame(
    sample = colnames(norm_mat),
    expression = as.numeric(norm_mat[gene_ens, ]),
    condition = colData_sub[[group_col]]
  )
  
  df$condition <- ifelse(df$condition == condition_label, "Cancer", "Healthy")
  df$condition <- factor(df$condition, levels = c("Healthy", "Cancer"))
  
  ggplot(df, aes(condition, expression, fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
    labs(
      title = paste0("Expression of ", gene),
      y = "VST expression",
      x = ""
    ) +
    theme_minimal(base_size = 13) +
    scale_fill_manual(values = c("Healthy" = "#4DBBD5", "Cancer" = "#E64B35")) +
    theme(legend.position = "none")
}





