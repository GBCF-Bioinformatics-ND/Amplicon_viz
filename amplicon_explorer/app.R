# Created by: Bharat Mishra
# Edited by: Elizabeth Brooks
# Modified: 21 Sept 2026

required_pkgs <- c("shiny", "phyloseq", "ggplot2", "vegan")
optional_pkgs <- c("biomformat")

missing_required <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0) {
  stop(
    paste0(
      "Missing required packages: ",
      paste(missing_required, collapse = ", "),
      "\nInstall CRAN packages with install.packages(...).",
      "\nInstall Bioconductor packages with BiocManager::install(...)."
    ),
    call. = FALSE
  )
}

missing_optional <- optional_pkgs[!vapply(optional_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_optional) > 0) {
  message(
    "Optional packages not installed (some features may be unavailable): ",
    paste(missing_optional, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(shiny)
  library(phyloseq)
  library(ggplot2)
  library(vegan)
})

options(shiny.maxRequestSize = 500 * 1024^2)

rank_names_default <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

default_paths <- list(
  otu = "data/otu_table.csv",
  tax = "data/taxonomy_table.csv",
  meta = "data/sample_metadata.csv",
  biom = "data/centrifuge_reports.biom"
)


read_table_with_rownames <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(x) < 2) {
    stop("Input table has fewer than 2 columns.")
  }

  first_name <- colnames(x)[1]
  if (is.na(first_name) || first_name %in% c("", "X", "...1")) {
    rownames(x) <- as.character(x[[1]])
    x <- x[, -1, drop = FALSE]
  }
  x
}

normalize_taxonomy <- function(tax_df) {
  tax <- as.data.frame(tax_df, stringsAsFactors = FALSE)
  tax[] <- lapply(tax, function(col) {
    value <- as.character(col)
    value <- trimws(value)
    value[value %in% c("", "NA", "N/A")] <- NA_character_
    value
  })

  if (all(grepl("^Rank", colnames(tax), ignore.case = TRUE)) && ncol(tax) <= length(rank_names_default)) {
    colnames(tax) <- rank_names_default[seq_len(ncol(tax))]
  }

  tax
}

normalize_tax_label <- function(x) {
  y <- trimws(as.character(x))
  y <- gsub("^[A-Za-z]__", "", y)
  tolower(y)
}

sanitize_count_table <- function(df, table_label = "Count table") {
  raw <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  raw_mat <- as.matrix(raw)

  raw_chr <- matrix(as.character(raw_mat), nrow = nrow(raw_mat), ncol = ncol(raw_mat))
  raw_chr <- trimws(raw_chr)

  missing_like <- is.na(raw_chr) | raw_chr %in% c("", "NA", "N/A", "na", "n/a", "NULL", "null")

  num_mat <- suppressWarnings(matrix(as.numeric(raw_chr), nrow = nrow(raw_chr), ncol = ncol(raw_chr)))
  invalid <- is.na(num_mat) & !missing_like
  if (any(invalid)) {
    invalid_idx <- which(invalid, arr.ind = TRUE)
    max_show <- min(5, nrow(invalid_idx))
    preview <- apply(invalid_idx[seq_len(max_show), , drop = FALSE], 1, function(i) {
      sprintf("row '%s', column '%s'", rownames(raw_mat)[i[1]], colnames(raw_mat)[i[2]])
    })
    stop(
      sprintf(
        "%s contains non-numeric values outside recognized missing markers at: %s",
        table_label,
        paste(preview, collapse = "; ")
      )
    )
  }

  n_missing <- sum(missing_like, na.rm = TRUE)
  if (n_missing > 0) {
    warning(sprintf("%s: replacing %d missing values with 0.", table_label, n_missing))
  }

  num_mat[is.na(num_mat)] <- 0
  colnames(num_mat) <- colnames(raw_mat)
  rownames(num_mat) <- rownames(raw_mat)

  as.data.frame(num_mat, check.names = FALSE)
}

ensure_taxonomy_columns <- function(ps) {
  tt <- tax_table(ps, errorIfNULL = FALSE)
  if (is.null(tt)) {
    fallback_tax <- matrix(
      taxa_names(ps),
      ncol = 1,
      dimnames = list(taxa_names(ps), "FeatureID")
    )
    tax_table(ps) <- tax_table(fallback_tax)
    return(ps)
  }

  if (is.null(colnames(tt)) || any(colnames(tt) == "")) {
    colnames(tt) <- paste0("Rank", seq_len(ncol(tt)))
    tax_table(ps) <- tt
  }

  ps
}

attach_or_build_metadata <- function(ps, meta_path = NULL) {
  if (!is.null(meta_path) && nzchar(meta_path)) {
    meta <- read.csv(meta_path, check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(meta) < 2) {
      stop("Metadata must include a sample ID column plus at least one metadata column.")
    }

    normalize_sample_id <- function(x) {
      y <- trimws(as.character(x))
      y <- gsub("_kreport$", "", y, ignore.case = TRUE)
      y <- gsub("\\.fastq(\\.gz)?$", "", y, ignore.case = TRUE)
      y <- gsub("\\.fq(\\.gz)?$", "", y, ignore.case = TRUE)
      y <- gsub("[[:space:]]+", "", y)
      y
    }

    meta_ids_raw <- as.character(meta[[1]])
    rownames(meta) <- meta_ids_raw
    meta <- meta[, -1, drop = FALSE]

    biom_ids_raw <- sample_names(ps)
    common_raw <- intersect(biom_ids_raw, rownames(meta))

    if (length(common_raw) < 2) {
      biom_norm <- normalize_sample_id(biom_ids_raw)
      meta_norm <- normalize_sample_id(rownames(meta))

      biom_lookup <- setNames(biom_ids_raw, biom_norm)
      meta_lookup <- setNames(rownames(meta), meta_norm)
      common_norm <- intersect(names(biom_lookup), names(meta_lookup))

      if (length(common_norm) >= 2) {
        common_biom <- unname(biom_lookup[common_norm])
        common_meta <- unname(meta_lookup[common_norm])

        ps <- prune_samples(common_biom, ps)
        meta_sub <- meta[common_meta, , drop = FALSE]
        rownames(meta_sub) <- sample_names(ps)
        sample_data(ps) <- sample_data(meta_sub)
        return(ps)
      }

      warning(
        "Could not match BIOM and metadata sample IDs (even after normalization). ",
        "Falling back to BIOM-only sample metadata."
      )
      fallback_meta <- data.frame(
        SampleID = sample_names(ps),
        row.names = sample_names(ps),
        stringsAsFactors = FALSE
      )
      sample_data(ps) <- sample_data(fallback_meta)
      return(ps)
    }

    ps <- prune_samples(common_raw, ps)
    sample_data(ps) <- sample_data(meta[common_raw, , drop = FALSE])
  } else {
    fallback_meta <- data.frame(
      SampleID = sample_names(ps),
      row.names = sample_names(ps),
      stringsAsFactors = FALSE
    )
    sample_data(ps) <- sample_data(fallback_meta)
  }

  ps
}

postprocess_ps <- function(ps, assay_type, keep_kingdom) {
  ps <- ensure_taxonomy_columns(ps)

  if (!identical(keep_kingdom, "All")) {
    tt <- tax_table(ps, errorIfNULL = FALSE)
    if (!is.null(tt)) {
      tt_mat <- as(tt, "matrix")
      rank_candidates <- c("Kingdom", "kingdom", "Superkingdom", "superkingdom", "Domain", "domain")
      rank_col <- rank_candidates[rank_candidates %in% colnames(tt_mat)]

      if (length(rank_col) > 0) {
        rank_col <- rank_col[1]
        kingdom_vals <- normalize_tax_label(tt_mat[, rank_col])
        target <- normalize_tax_label(keep_kingdom)

        keep_idx <- !is.na(kingdom_vals) & kingdom_vals != "" &
          (kingdom_vals == target | grepl(target, kingdom_vals, fixed = TRUE))

        if (any(keep_idx)) {
          ps <- prune_taxa(taxa_names(ps)[keep_idx], ps)
        } else {
          warning(
            sprintf(
              "Kingdom filter '%s' found no matching taxa in column '%s'. Keeping all taxa.",
              keep_kingdom,
              rank_col
            )
          )
        }
      } else {
        warning("No Kingdom/Domain taxonomy column found. Kingdom filter was skipped.")
      }
    }
  }

  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  sample_data(ps)$AssayType <- assay_type
  ps
}

build_phyloseq_from_csv <- function(otu_path, tax_path, meta_path, assay_type, keep_kingdom) {
  otu <- read_table_with_rownames(otu_path)
  tax <- normalize_taxonomy(read_table_with_rownames(tax_path))
  meta <- read.csv(meta_path, check.names = FALSE, stringsAsFactors = FALSE)

  if (ncol(meta) < 2) {
    stop("Metadata must include a sample ID column plus at least one metadata column.")
  }

  rownames(meta) <- as.character(meta[[1]])
  meta <- meta[, -1, drop = FALSE]

  otu <- sanitize_count_table(otu, table_label = "OTU table")

  common <- Reduce(intersect, list(colnames(otu), rownames(meta)))
  if (length(common) < 2) {
    stop("Fewer than 2 overlapping samples between OTU table and metadata.")
  }

  otu <- otu[, common, drop = FALSE]
  meta <- meta[common, , drop = FALSE]

  common_taxa <- intersect(rownames(otu), rownames(tax))
  if (length(common_taxa) < 2) {
    stop("Fewer than 2 overlapping taxa between OTU and taxonomy tables.")
  }

  otu <- otu[common_taxa, , drop = FALSE]
  tax <- tax[common_taxa, , drop = FALSE]

  ps <- phyloseq(
    otu_table(as.matrix(otu), taxa_are_rows = TRUE),
    tax_table(as.matrix(tax)),
    sample_data(meta)
  )

  postprocess_ps(ps, assay_type, keep_kingdom)
}

build_phyloseq_from_biom <- function(biom_path, meta_path, assay_type, keep_kingdom) {
  if (!requireNamespace("biomformat", quietly = TRUE)) {
    stop("Package 'biomformat' is required for BIOM imports. Please install it first.")
  }

  ps <- import_biom(biom_path)
  ps <- attach_or_build_metadata(ps, meta_path)
  postprocess_ps(ps, assay_type, keep_kingdom)
}

parse_qiime2_taxonomy <- function(taxon_values) {
  ranks <- rank_names_default
  out <- matrix(NA_character_, nrow = length(taxon_values), ncol = length(ranks))
  colnames(out) <- ranks

  for (i in seq_along(taxon_values)) {
    val <- trimws(as.character(taxon_values[i]))
    if (is.na(val) || !nzchar(val)) {
      next
    }

    parts <- trimws(unlist(strsplit(val, ";", fixed = TRUE)))
    parts <- parts[nzchar(parts)]
    if (length(parts) == 0) {
      next
    }

    parts <- gsub("^[A-Za-z]__", "", parts)
    parts[parts == ""] <- NA_character_
    n_use <- min(length(parts), length(ranks))
    out[i, seq_len(n_use)] <- parts[seq_len(n_use)]
  }

  as.data.frame(out, stringsAsFactors = FALSE)
}

read_qiime2_feature_table <- function(feature_path) {
  lines <- readLines(feature_path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  lines <- lines[!grepl("^#\\s*Constructed", lines)]
  lines <- lines[!grepl("^#\\s*OTU table", lines)]

  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(lines, tmp)

  ft <- read.delim(tmp, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
  if (ncol(ft) < 2) {
    stop("QIIME2 feature table appears malformed.")
  }

  # QIIME2 exports often include a row declaring per-column data types.
  if (nrow(ft) > 0 && identical(as.character(ft[[1]][1]), "#q2:types")) {
    ft <- ft[-1, , drop = FALSE]
  }

  rownames(ft) <- as.character(ft[[1]])
  ft <- ft[, -1, drop = FALSE]

  tax_col <- NULL
  if ("Taxon" %in% colnames(ft)) {
    tax_col <- as.character(ft$Taxon)
  } else if ("taxonomy" %in% colnames(ft)) {
    tax_col <- as.character(ft$taxonomy)
  } else if ("Consensus.Lineage" %in% colnames(ft)) {
    tax_col <- as.character(ft$Consensus.Lineage)
  }

  confidence_col <- NULL
  if ("Confidence" %in% colnames(ft)) {
    confidence_col <- as.character(ft$Confidence)
  } else if ("confidence" %in% colnames(ft)) {
    confidence_col <- as.character(ft$confidence)
  }

  tax_df <- NULL
  if (!is.null(tax_col)) {
    tax_df <- parse_qiime2_taxonomy(tax_col)
    tax_df <- cbind(
      Taxon = tax_col,
      taxonomy = tax_col,
      tax_df,
      stringsAsFactors = FALSE
    )
    if (!is.null(confidence_col)) {
      tax_df$Confidence <- confidence_col
    }
    rownames(tax_df) <- rownames(ft)
  }

  # Drop common non-count annotation columns.
  annotation_cols <- c("Taxon", "taxonomy", "Consensus.Lineage", "Confidence", "confidence")
  ft <- ft[, setdiff(colnames(ft), annotation_cols), drop = FALSE]

  # Keep only columns that are fully numeric-like (sample count columns).
  is_numeric_like_col <- vapply(ft, function(col) {
    vals <- trimws(as.character(col))
    vals <- vals[!is.na(vals)]
    vals <- vals[!vals %in% c("", "NA", "N/A", "na", "n/a", "NULL", "null")]
    if (length(vals) == 0) {
      return(FALSE)
    }
    all(grepl("^-?[0-9]+(\\.[0-9]+)?$", vals))
  }, logical(1))

  ft <- ft[, is_numeric_like_col, drop = FALSE]
  if (ncol(ft) < 2) {
    stop("QIIME2 feature table has fewer than 2 numeric sample columns after filtering.")
  }

  ft <- sanitize_count_table(ft, table_label = "QIIME2 feature table")

  list(otu = ft, tax = tax_df)
}

read_metadata_flexible <- function(meta_path) {
  ext <- tolower(tools::file_ext(meta_path))
  if (ext == "csv") {
    md <- read.csv(meta_path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    md <- read.delim(meta_path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
  }

  if (ncol(md) < 2) {
    stop("Metadata must include a sample ID column plus at least one metadata column.")
  }

  first_col <- as.character(md[[1]])
  if (length(first_col) > 0 && identical(first_col[1], "#q2:types")) {
    md <- md[-1, , drop = FALSE]
  }

  rownames(md) <- as.character(md[[1]])
  md <- md[, -1, drop = FALSE]
  md
}

build_phyloseq_from_qiime2 <- function(feature_path, meta_path, assay_type, keep_kingdom) {
  qiime2_data <- read_qiime2_feature_table(feature_path)
  otu <- qiime2_data$otu
  tax <- qiime2_data$tax
  meta <- read_metadata_flexible(meta_path)

  common <- intersect(colnames(otu), rownames(meta))
  if (length(common) < 2) {
    stop("Fewer than 2 overlapping samples between QIIME2 feature table and metadata.")
  }

  otu <- otu[, common, drop = FALSE]
  meta <- meta[common, , drop = FALSE]

  if (is.null(tax)) {
    tax <- data.frame(
      FeatureID = rownames(otu),
      row.names = rownames(otu),
      stringsAsFactors = FALSE
    )
  } else {
    tax <- tax[rownames(otu), , drop = FALSE]
  }

  tax <- normalize_taxonomy(tax)

  ps <- phyloseq(
    otu_table(as.matrix(otu), taxa_are_rows = TRUE),
    tax_table(as.matrix(tax)),
    sample_data(meta)
  )

  postprocess_ps(ps, assay_type, keep_kingdom)
}

analysis_help_block <- function(title, bullets) {
  wellPanel(
    tags$h4(title),
    tags$ul(lapply(bullets, tags$li))
  )
}

app_ui <- fluidPage(
  titlePanel("Amplicon Explorer: 16S / ITS"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      radioButtons(
        "assay_type",
        "Assay type",
        choices = c("16S", "ITS"),
        selected = "16S",
        inline = TRUE
      ),
      selectInput(
        "keep_kingdom",
        "Kingdom filter",
        choices = c("All", "Bacteria", "Fungi"),
        selected = "All"
      ),
      tags$hr(),
      radioButtons(
        "data_mode",
        "Input mode",
        choices = c("CSV tables" = "csv", "BIOM" = "biom", "QIIME2 feature table" = "qiime2"),
        selected = "csv"
      ),
      checkboxInput("use_defaults", "Use workspace default files", value = FALSE),
      conditionalPanel(
        condition = "input.data_mode == 'csv'",
        fileInput("otu_file", "OTU/ASV count table (CSV)", accept = ".csv"),
        fileInput("tax_file", "Taxonomy table (CSV)", accept = ".csv"),
        fileInput("meta_file", "Sample metadata (CSV)", accept = ".csv")
      ),
      conditionalPanel(
        condition = "input.data_mode == 'biom'",
        fileInput("biom_file", "Feature table (BIOM)", accept = c(".biom", ".json")),
        fileInput("biom_meta_file", "Sample metadata (CSV, optional)", accept = ".csv")
      ),
      conditionalPanel(
        condition = "input.data_mode == 'qiime2'",
        fileInput("qiime2_feature_file", "QIIME2 feature table (txt/tsv)", accept = c(".txt", ".tsv")),
        fileInput("qiime2_meta_file", "Sample metadata (csv/txt/tsv)", accept = c(".csv", ".txt", ".tsv"))
      ),
      actionButton("load_data", "Load / Reload data", class = "btn-primary"),
      tags$hr(),
      uiOutput("group_var_ui"),
      uiOutput("facet_var_ui"),
      uiOutput("tax_rank_ui"),
      sliderInput("top_n", "Top taxa to show", min = 5, max = 30, value = 10)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel(
          "Guide",
          tags$br(),
          tags$h3("How To Use This App"),
          tags$p("This app explores amplicon count tables, taxonomy, and sample metadata to describe within-sample diversity, between-sample community structure, taxonomic composition, and prevalent taxa."),
          analysis_help_block(
            "Recommended workflow",
            list(
              "Load an OTU/ASV table, taxonomy table, and metadata file, or load a BIOM file with optional metadata.",
              "Check the Data summary tab for sample counts, feature counts, sequencing depth, and taxonomy availability.",
              "Use Alpha diversity to compare richness and evenness within individual samples.",
              "Use Beta ordination, clustering, and PERMANOVA to examine whole-community differences between samples or groups.",
              "Use Taxa composition to inspect relative-abundance patterns at a selected taxonomic rank.",
              "Use Core microbiome to find taxa that meet a chosen prevalence and abundance threshold.",
              "Treat plots as exploratory evidence and interpret them with replication, effect size, and statistical results."
            )
          ),
          analysis_help_block(
            "Important concepts",
            list(
              "A taxon is a biological group such as a phylum, genus, or species.",
              "Relative abundance means counts are converted to proportions within each sample, so samples can be compared even if sequencing depth differs.",
              "Alpha diversity describes diversity within one sample.",
              "Beta diversity describes differences between samples.",
              "An ordination plot places similar samples close together and dissimilar samples farther apart.",
              "Relative abundance is the proportion of reads assigned to each taxon within a sample; it is not the same as absolute cell abundance.",
              "A p-value measures evidence against a null hypothesis, not the size or biological importance of an effect.",
              "PERMANOVA can reflect differences in group centroids or within-group dispersion, so inspect ordination and replication together."
            )
          ),
          analysis_help_block(
            "Good practice for interpretation",
            list(
              "Look for consistent patterns across multiple plots instead of relying on one result.",
              "Interpret statistics together with effect size, group separation, dispersion, and sample size.",
              "Be careful when groups have very different numbers of samples.",
              "A statistically significant result does not always mean a biologically large effect.",
              "Taxonomic labels can be incomplete; an unassigned feature is still included in abundance calculations."
            )
          )
        ),
        tabPanel(
          "Data summary",
          tags$br(),
          analysis_help_block(
            "What this tab shows",
            list(
              "The number of samples and features loaded after matching and filtering.",
              "The total number of reads and sequencing depth represented in the current object.",
              "The taxonomy ranks available for labeling composition plots.",
              "A preview of the sample metadata used to color, group, and facet plots.",
              "A taxonomy preview showing the labels associated with feature IDs."
            )
          ),
          verbatimTextOutput("summary_text"),
          tags$h4("Sample metadata preview"),
          tableOutput("sample_table"),
          downloadButton("download_sample_table", "Download sample table"),
          tags$br(),
          tags$br(),
          tags$h4("Taxonomy preview"),
          tableOutput("taxonomy_table"),
          downloadButton("download_taxonomy_table", "Download taxonomy table")
        ),
        tabPanel(
          "Alpha diversity",
          tags$br(),
          analysis_help_block(
            "How to read alpha diversity",
            list(
              "Observed richness counts how many taxa were detected in each sample.",
              "Shannon diversity increases when a sample has both many taxa and a more even distribution among them.",
              "Simpson diversity gives more weight to dominant taxa.",
              "Each point is one sample. The boxplot summarizes the group.",
              "Differences between groups suggest different within-sample diversity; use the points and group sizes to assess consistency."
            )
          ),
          selectInput("alpha_measure", "Alpha metric", choices = c("Observed", "Shannon", "Simpson")),
          plotOutput("alpha_plot", height = 420),
          downloadButton("download_alpha_plot", "Download alpha plot"),
          downloadButton("download_alpha_table", "Download alpha table")
        ),
        tabPanel(
          "Beta ordination",
          tags$br(),
          analysis_help_block(
            "How to read ordination",
            list(
              "Each point is one sample; the axes summarize distances among samples rather than individual taxa.",
              "Samples that are close together have more similar community composition under the selected distance.",
              "Samples that are far apart have more different community composition, but axis values are not direct abundance measurements.",
              "Bray-Curtis uses abundance, Jaccard uses presence/absence, and Euclidean uses the numeric transformed table.",
              "Group separation is exploratory; assess it alongside PERMANOVA and within-group spread."
            )
          ),
          fluidRow(
            column(
              4,
              selectInput("distance_method", "Distance", choices = c("bray", "jaccard", "euclidean"), selected = "bray")
            ),
            column(
              4,
              selectInput("ordination_method", "Ordination", choices = c("NMDS", "PCoA"), selected = "NMDS")
            )
          ),
          plotOutput("ordination_plot", height = 460),
          downloadButton("download_ordination_plot", "Download ordination plot"),
          downloadButton("download_ordination_table", "Download ordination table")
        ),
        tabPanel(
          "Taxa composition",
          tags$br(),
          analysis_help_block(
            "How to read taxa composition",
            list(
              "Each bar is one sample.",
              "Colors represent taxa at the selected rank, such as phylum or genus.",
              "Bar height segments show relative abundance, not absolute count.",
              "This plot is useful for identifying dominant taxa and broad relative-abundance shifts.",
              "The 'Other' category groups lower-abundance taxa so the plot stays readable; it does not mean those taxa are absent."
            )
          ),
          plotOutput("taxa_plot", height = 520),
          downloadButton("download_taxa_plot", "Download taxa plot"),
          downloadButton("download_taxa_table", "Download taxa table")
        ),
        tabPanel(
          "Clustering",
          tags$br(),
          analysis_help_block(
            "How to read clustering",
            list(
              "Samples connected by shorter branches are more similar under the selected distance.",
              "Samples that join only near the top of the tree are less similar to the rest of the cluster.",
              "Clustering can reveal sample groupings or potential outliers, but the tree is dependent on the distance and linkage method.",
              "Terminal labels are colored by the selected metadata grouping variable; the legend identifies each group.",
              "This is exploratory and should be interpreted together with ordination, metadata, and replication."
            )
          ),
          plotOutput("cluster_plot", height = 520),
          downloadButton("download_cluster_plot", "Download clustering plot"),
          downloadButton("download_cluster_table", "Download clustering table")
        ),
        tabPanel(
          "PERMANOVA",
          tags$br(),
          analysis_help_block(
            "How to read PERMANOVA",
            list(
              "PERMANOVA tests whether overall community composition differs among groups.",
              "The R2 value estimates how much of the variation is explained by the grouping variable.",
              "A small p-value suggests group centroids differ more than expected under the permutation null model.",
              "R2 is an effect-size-like measure of explained variation, not the percent of samples classified correctly.",
              "PERMANOVA can also be sensitive to unequal within-group dispersion, so inspect ordination spread and group sizes."
            )
          ),
          verbatimTextOutput("permanova_text"),
          downloadButton("download_permanova_table", "Download PERMANOVA table")
        ),
        tabPanel(
          "Core microbiome",
          tags$br(),
          analysis_help_block(
            "How to read the core microbiome",
            list(
              "The core microbiome is the set of taxa found across many samples at or above a chosen abundance threshold.",
              "Prevalence threshold asks: in what fraction of samples must a taxon appear?",
              "Minimum abundance asks: how abundant must a taxon be before we count it as present?",
              "The heatmap shows abundance patterns for taxa that pass both thresholds; it does not prove that they are biologically essential.",
              "Changing the taxonomic rank changes how features are aggregated and labeled in the core summary."
            )
          ),
          uiOutput("core_rank_ui"),
          sliderInput("core_prev", "Prevalence threshold", min = 0.5, max = 1, value = 0.8, step = 0.05),
          sliderInput("core_abund", "Minimum abundance (%)", min = 0.01, max = 1, value = 0.1, step = 0.01),
          plotOutput("core_heatmap", height = 520),
          downloadButton("download_core_plot", "Download core heatmap"),
          downloadButton("download_core_table", "Download core table")
        ),
      )
    )
  )
)

app_server <- function(input, output, session) {
  require_tab <- function(tab_name) {
    req(input$tabs == tab_name)
  }

  get_inputs <- reactive({
    if (isTRUE(input$use_defaults)) {
      if (identical(input$data_mode, "csv")) {
        required_defaults <- unlist(default_paths[c("otu", "tax", "meta")])
        missing_defaults <- names(required_defaults)[!file.exists(required_defaults)]
        if (length(missing_defaults) > 0) {
          stop(sprintf("Default files missing for CSV mode: %s", paste(missing_defaults, collapse = ", ")))
        }
        return(list(mode = "csv", otu = default_paths$otu, tax = default_paths$tax, meta = default_paths$meta))
      }

      if (!file.exists(default_paths$biom)) {
        stop("Default BIOM file missing: centrifuge_reports.biom")
      }

      biom_meta_path <- if (file.exists(default_paths$meta)) default_paths$meta else ""
      return(list(mode = "biom", biom = default_paths$biom, biom_meta = biom_meta_path))
    }

    if (identical(input$data_mode, "csv")) {
      req(input$otu_file, input$tax_file, input$meta_file)
      return(list(mode = "csv", otu = input$otu_file$datapath, tax = input$tax_file$datapath, meta = input$meta_file$datapath))
    }

    if (identical(input$data_mode, "qiime2")) {
      req(input$qiime2_feature_file, input$qiime2_meta_file)
      return(list(mode = "qiime2", feature = input$qiime2_feature_file$datapath, meta = input$qiime2_meta_file$datapath))
    }

    req(input$biom_file)
    list(
      mode = "biom",
      biom = input$biom_file$datapath,
      biom_meta = if (!is.null(input$biom_meta_file)) input$biom_meta_file$datapath else ""
    )
  })

  ps_obj <- eventReactive(input$load_data, {
    data_inputs <- get_inputs()
    if (identical(data_inputs$mode, "csv")) {
      return(
        build_phyloseq_from_csv(
          otu_path = data_inputs$otu,
          tax_path = data_inputs$tax,
          meta_path = data_inputs$meta,
          assay_type = input$assay_type,
          keep_kingdom = input$keep_kingdom
        )
      )
    }

    if (identical(data_inputs$mode, "qiime2")) {
      return(
        build_phyloseq_from_qiime2(
          feature_path = data_inputs$feature,
          meta_path = data_inputs$meta,
          assay_type = input$assay_type,
          keep_kingdom = input$keep_kingdom
        )
      )
    }

    build_phyloseq_from_biom(
      biom_path = data_inputs$biom,
      meta_path = data_inputs$biom_meta,
      assay_type = input$assay_type,
      keep_kingdom = input$keep_kingdom
    )
  }, ignoreInit = FALSE)

  ps_rel <- reactive({
    transform_sample_counts(ps_obj(), function(x) x / sum(x))
  })

  dist_obj <- reactive({
    phyloseq::distance(ps_rel(), method = input$distance_method)
  })

  otu_sample_matrix <- reactive({
    ps <- ps_obj()
    mat <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) {
      mat <- t(mat)
    }
    mat
  })

  sample_table_df <- reactive({
    as(sample_data(ps_obj()), "data.frame")
  })

  taxonomy_table_df <- reactive({
    ps <- ps_obj()
    tt <- tax_table(ps, errorIfNULL = FALSE)
    validate(need(!is.null(tt), "No taxonomy table available in this dataset."))

    tax_df <- as.data.frame(tt, stringsAsFactors = FALSE)
    tax_df$FeatureID <- rownames(tax_df)

    preferred_cols <- c("FeatureID", "Taxon", "taxonomy", "Confidence", "confidence")
    present_preferred <- intersect(preferred_cols, colnames(tax_df))

    if (length(present_preferred) > 1) {
      out <- tax_df[, unique(present_preferred), drop = FALSE]
    } else {
      other_cols <- setdiff(colnames(tax_df), "FeatureID")
      keep_other <- head(other_cols, 6)
      out <- tax_df[, c("FeatureID", keep_other), drop = FALSE]
    }

    out
  })

  alpha_df <- reactive({
    req(input$group_var)
    ps <- ps_obj()
    df <- estimate_richness(ps, measures = input$alpha_measure)
    df$Sample <- rownames(df)

    meta <- as(sample_data(ps), "data.frame")
    meta$Sample <- rownames(meta)
    out <- merge(df, meta, by = "Sample")
    keep <- !is.na(out[[input$group_var]]) & trimws(as.character(out[[input$group_var]])) != ""
    out[keep, , drop = FALSE]
  })

  alpha_plot_obj <- reactive({
    df <- alpha_df()
    ggplot(df, aes_string(x = input$group_var, y = input$alpha_measure, color = input$group_var)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.25) +
      geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
      theme_bw(base_size = 13) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
      labs(x = input$group_var, y = input$alpha_measure)
  })

  ordination_result <- reactive({
    req(input$group_var)
    ps_rel_local <- ps_rel()
    meta_df <- as(sample_data(ps_rel_local), "data.frame")
    keep <- !is.na(meta_df[[input$group_var]]) & trimws(as.character(meta_df[[input$group_var]])) != ""
    ps_rel_local <- prune_samples(rownames(meta_df)[keep], ps_rel_local)
    validate(need(nsamples(ps_rel_local) >= 2, "Need at least 2 samples with non-missing group values for ordination."))

    if (identical(input$ordination_method, "NMDS")) {
      ord <- ordinate(ps_rel_local, method = "NMDS", distance = input$distance_method)
    } else {
      ord <- ordinate(ps_rel_local, method = "PCoA", distance = input$distance_method)
    }
    list(phyloseq = ps_rel_local, ordination = ord)
  })

  ordination_df <- reactive({
    result <- ordination_result()
    plot_ordination(
      result$phyloseq,
      result$ordination,
      color = input$group_var,
      justDF = TRUE
    )
  })

  ordination_plot_obj <- reactive({
    result <- ordination_result()
    plot_ordination(result$phyloseq, result$ordination, color = input$group_var) +
      geom_point(size = 3, alpha = 0.9) +
      theme_bw(base_size = 13)
  })

  taxa_df <- reactive({
    req(input$tax_rank)
    ps_rank <- tax_glom(ps_rel(), taxrank = input$tax_rank)
    df <- psmelt(ps_rank)
    df$Abundance <- df$Abundance * 100
    df[[input$tax_rank]] <- as.character(df[[input$tax_rank]])
    df[[input$tax_rank]][is.na(df[[input$tax_rank]]) | df[[input$tax_rank]] == ""] <- "Unassigned"
    top_taxa <- names(sort(tapply(df$Abundance, df[[input$tax_rank]], sum), decreasing = TRUE))[1:input$top_n]
    df[[input$tax_rank]][!df[[input$tax_rank]] %in% top_taxa] <- "Other"
    df
  })

  taxa_plot_obj <- reactive({
    df <- taxa_df()
    p <- ggplot(df, aes_string(x = "Sample", y = "Abundance", fill = input$tax_rank)) +
      geom_bar(stat = "identity") +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
      labs(y = "Relative abundance (%)", fill = input$tax_rank)

    if (!is.null(input$facet_var) && input$facet_var != "None") {
      p <- p + facet_wrap(stats::as.formula(paste("~", input$facet_var)), scales = "free_x")
    }
    p
  })

  cluster_obj <- reactive({
    hclust(dist_obj(), method = "average")
  })

  draw_cluster_plot <- function() {
    req(input$group_var)
    hc <- cluster_obj()
    meta_df <- as(sample_data(ps_rel()), "data.frame")
    req(input$group_var %in% colnames(meta_df))

    group_values <- trimws(as.character(meta_df[hc$labels, input$group_var]))
    group_values[is.na(group_values) | group_values == ""] <- "Missing"
    group_levels <- sort(unique(group_values))
    group_colors <- setNames(grDevices::hcl.colors(length(group_levels), "Dark 3"), group_levels)
    label_colors <- unname(group_colors[group_values[hc$order]])
    ordered_labels <- hc$labels[hc$order]

    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mar = c(10, 4, 4, 2) + 0.1)
    plot(
      hc,
      labels = FALSE,
      main = paste("UPGMA clustering (", input$distance_method, ")", sep = ""),
      xlab = "",
      sub = "",
      axes = FALSE
    )
    axis(2)
    box()
    text(
      x = seq_along(ordered_labels),
      y = par("usr")[3],
      labels = ordered_labels,
      col = label_colors,
      srt = 90,
      adj = 1,
      xpd = NA,
      cex = 0.8
    )
    legend(
      "topright",
      legend = names(group_colors),
      col = unname(group_colors),
      pch = 15,
      title = input$group_var,
      bty = "n",
      cex = 0.8
    )
  }

  cluster_table <- reactive({
    hc <- cluster_obj()
    data.frame(
      Sample = hc$labels[hc$order],
      ClusterOrder = seq_along(hc$order),
      stringsAsFactors = FALSE
    )
  })

  permanova_result <- reactive({
    req(input$group_var)
    ps_rel_local <- ps_rel()
    meta_df <- as(sample_data(ps_rel_local), "data.frame")
    req(input$group_var %in% colnames(meta_df))

    keep <- !is.na(meta_df[[input$group_var]]) & trimws(as.character(meta_df[[input$group_var]])) != ""
    validate(need(sum(keep) >= 2, "Need at least 2 samples with non-missing group values for PERMANOVA."))

    ps_sub <- prune_samples(rownames(meta_df)[keep], ps_rel_local)
    meta_sub <- as(sample_data(ps_sub), "data.frame")
    meta_sub$GroupVar <- as.factor(as.character(meta_sub[[input$group_var]]))
    validate(need(length(unique(meta_sub$GroupVar)) >= 2, "PERMANOVA needs at least 2 groups after removing missing values."))

    dist_sub <- phyloseq::distance(ps_sub, method = input$distance_method)
    vegan::adonis2(dist_sub ~ GroupVar, data = meta_sub)
  })

  permanova_table <- reactive({
    ad <- permanova_result()
    out <- as.data.frame(ad)
    out$Term <- rownames(out)
    rownames(out) <- NULL
    out
  })

  observe({
    ps <- ps_obj()
    meta_df <- as(sample_data(ps), "data.frame")
    meta_cols <- colnames(meta_df)
    rank_cols <- rank_names(ps)

    # Prefer a grouping variable with repeated levels for downstream comparisons.
    group_candidates <- meta_cols[vapply(meta_cols, function(col) {
      vals <- trimws(as.character(meta_df[[col]]))
      vals <- vals[!is.na(vals) & vals != ""]
      if (length(vals) < 2) {
        return(FALSE)
      }
      lev_counts <- table(vals)
      length(lev_counts) >= 2 && any(lev_counts >= 2)
    }, logical(1))]

    preferred_group <- if (length(group_candidates) > 0) group_candidates[1] else meta_cols[1]

    updateSelectInput(session, "alpha_measure", selected = "Shannon")

    output$group_var_ui <- renderUI({
      selectInput("group_var", "Group/color variable", choices = meta_cols, selected = preferred_group)
    })

    output$facet_var_ui <- renderUI({
      selectInput("facet_var", "Facet variable (taxa plot)", choices = c("None", meta_cols), selected = "None")
    })

    output$tax_rank_ui <- renderUI({
      selectInput("tax_rank", "Taxonomic rank", choices = rank_cols, selected = if ("Genus" %in% rank_cols) "Genus" else rank_cols[1])
    })

    output$core_rank_ui <- renderUI({
      selectInput("core_tax_rank", "Core microbiome label rank", choices = rank_cols, selected = if ("Genus" %in% rank_cols) "Genus" else rank_cols[1])
    })

    output$da_var_ui <- renderUI({
      selectInput("da_var", "Model variable", choices = meta_cols, selected = meta_cols[1])
    })

  })

  observe({
    req(input$da_var)
    meta_df <- as(sample_data(ps_obj()), "data.frame")
    req(input$da_var %in% colnames(meta_df))

    levs <- unique(as.character(meta_df[[input$da_var]]))
    levs <- levs[!is.na(levs) & levs != ""]
    if (length(levs) < 2) {
      output$da_ref_level_ui <- renderUI({
        helpText("Selected variable has fewer than 2 non-empty levels.")
      })
      output$da_comp_level_ui <- renderUI({
        helpText("")
      })
      return()
    }

    output$da_ref_level_ui <- renderUI({
      selectInput("da_ref_level", "Reference level", choices = levs, selected = levs[1])
    })
    output$da_comp_level_ui <- renderUI({
      selectInput("da_comp_level", "Comparison level", choices = levs, selected = levs[2])
    })
  })


  output$summary_text <- renderPrint({
    ps <- ps_obj()
    cat("Amplicon dataset loaded successfully\n")
    cat(sprintf("Samples: %d\n", nsamples(ps)))
    cat(sprintf("Features: %d\n", ntaxa(ps)))
    cat(sprintf("Assay type: %s\n", unique(sample_data(ps)$AssayType)))
    cat(sprintf("Total reads: %s\n", format(sum(sample_sums(ps)), big.mark = ",")))
    cat("\nAvailable taxonomy ranks:\n")
    print(rank_names(ps))
  })

  output$sample_table <- renderTable({
    sample_table_df()
  }, rownames = TRUE)

  output$taxonomy_table <- renderTable({
    head(taxonomy_table_df(), 25)
  }, rownames = FALSE)

  output$alpha_plot <- renderPlot({
    require_tab("Alpha diversity")
    alpha_plot_obj()
  })

  output$ordination_plot <- renderPlot({
    require_tab("Beta ordination")
    ordination_plot_obj()
  })

  output$taxa_plot <- renderPlot({
    require_tab("Taxa composition")
    taxa_plot_obj()
  })

  output$cluster_plot <- renderPlot({
    require_tab("Clustering")
    draw_cluster_plot()
  })

  output$permanova_text <- renderPrint({
    require_tab("PERMANOVA")
    print(permanova_result())
  })

  output$core_heatmap <- renderPlot({
    require_tab("Core microbiome")
    req(input$core_tax_rank)
    ps_rel_local <- ps_rel()
    otu <- as(otu_table(ps_rel_local), "matrix")
    if (!taxa_are_rows(ps_rel_local)) {
      otu <- t(otu)
    }

    prev <- rowMeans(otu > (input$core_abund / 100))
    core_taxa <- names(prev[prev >= input$core_prev])

    validate(need(length(core_taxa) > 0, "No core taxa at current thresholds."))

    core_ps <- prune_taxa(core_taxa, ps_rel_local)
    core_df <- psmelt(core_ps)
    core_df$AbundancePct <- core_df$Abundance * 100

    if (!(input$core_tax_rank %in% colnames(core_df))) {
      core_df$TaxaLabel <- as.character(core_df$OTU)
    } else {
      core_df$TaxaLabel <- as.character(core_df[[input$core_tax_rank]])
      core_df$TaxaLabel[is.na(core_df$TaxaLabel) | core_df$TaxaLabel == ""] <- "Unassigned"
    }

    core_df <- aggregate(AbundancePct ~ Sample + TaxaLabel, data = core_df, FUN = sum)
    core_df <- core_df[order(core_df$TaxaLabel, core_df$Sample), , drop = FALSE]

    ggplot(core_df, aes(x = Sample, y = TaxaLabel, fill = AbundancePct)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "steelblue") +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
      labs(
        title = paste("Core microbiome heatmap:", input$core_tax_rank),
        x = "Sample",
        y = input$core_tax_rank,
        fill = "Abundance (%)"
      )
  })

  core_table <- reactive({
    req(input$core_tax_rank)
    ps_rel_local <- ps_rel()
    otu <- as(otu_table(ps_rel_local), "matrix")
    if (!taxa_are_rows(ps_rel_local)) {
      otu <- t(otu)
    }

    prev <- rowMeans(otu > (input$core_abund / 100))
    core_taxa <- names(prev[prev >= input$core_prev])
    validate(need(length(core_taxa) > 0, "No core taxa at current thresholds."))

    core_ps <- prune_taxa(core_taxa, ps_rel_local)
    core_df <- psmelt(core_ps)
    core_df$AbundancePct <- core_df$Abundance * 100

    if (!(input$core_tax_rank %in% colnames(core_df))) {
      core_df$TaxaLabel <- as.character(core_df$OTU)
    } else {
      core_df$TaxaLabel <- as.character(core_df[[input$core_tax_rank]])
      core_df$TaxaLabel[is.na(core_df$TaxaLabel) | core_df$TaxaLabel == ""] <- "Unassigned"
    }

    out <- aggregate(AbundancePct ~ Sample + TaxaLabel, data = core_df, FUN = sum)
    out[order(out$TaxaLabel, out$Sample), , drop = FALSE]
  })

  output$download_sample_table <- downloadHandler(
    filename = function() paste0("sample_metadata_", Sys.Date(), ".csv"),
    content = function(file) write.csv(sample_table_df(), file, row.names = TRUE)
  )

  output$download_taxonomy_table <- downloadHandler(
    filename = function() paste0("taxonomy_table_", Sys.Date(), ".csv"),
    content = function(file) write.csv(taxonomy_table_df(), file, row.names = FALSE)
  )

  output$download_alpha_plot <- downloadHandler(
    filename = function() paste0("alpha_diversity_", Sys.Date(), ".png"),
    content = function(file) ggsave(file, plot = alpha_plot_obj(), width = 9, height = 5, dpi = 300)
  )

  output$download_alpha_table <- downloadHandler(
    filename = function() paste0("alpha_diversity_table_", Sys.Date(), ".csv"),
    content = function(file) write.csv(alpha_df(), file, row.names = FALSE)
  )

  output$download_ordination_plot <- downloadHandler(
    filename = function() paste0("ordination_", Sys.Date(), ".png"),
    content = function(file) ggsave(file, plot = ordination_plot_obj(), width = 9, height = 5, dpi = 300)
  )

  output$download_ordination_table <- downloadHandler(
    filename = function() paste0("ordination_table_", Sys.Date(), ".csv"),
    content = function(file) write.csv(ordination_df(), file, row.names = FALSE)
  )

  output$download_taxa_plot <- downloadHandler(
    filename = function() paste0("taxa_composition_", Sys.Date(), ".png"),
    content = function(file) ggsave(file, plot = taxa_plot_obj(), width = 11, height = 6, dpi = 300)
  )

  output$download_taxa_table <- downloadHandler(
    filename = function() paste0("taxa_composition_table_", Sys.Date(), ".csv"),
    content = function(file) write.csv(taxa_df(), file, row.names = FALSE)
  )

  output$download_cluster_plot <- downloadHandler(
    filename = function() paste0("clustering_", Sys.Date(), ".png"),
    content = function(file) {
      png(file, width = 1200, height = 700, res = 120)
      draw_cluster_plot()
      dev.off()
    }
  )

  output$download_cluster_table <- downloadHandler(
    filename = function() paste0("clustering_table_", Sys.Date(), ".csv"),
    content = function(file) write.csv(cluster_table(), file, row.names = FALSE)
  )

  output$download_permanova_table <- downloadHandler(
    filename = function() paste0("permanova_", Sys.Date(), ".csv"),
    content = function(file) write.csv(permanova_table(), file, row.names = FALSE)
  )

  output$download_core_plot <- downloadHandler(
    filename = function() paste0("core_microbiome_", Sys.Date(), ".png"),
    content = function(file) {
      g <- ggplot(core_table(), aes(x = Sample, y = TaxaLabel, fill = AbundancePct)) +
        geom_tile() +
        scale_fill_gradient(low = "white", high = "steelblue") +
        theme_bw(base_size = 11) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        labs(
          title = paste("Core microbiome heatmap:", input$core_tax_rank),
          x = "Sample",
          y = input$core_tax_rank,
          fill = "Abundance (%)"
        )
      ggsave(file, plot = g, width = 11, height = 6, dpi = 300)
    }
  )

  output$download_core_table <- downloadHandler(
    filename = function() paste0("core_microbiome_table_", Sys.Date(), ".csv"),
    content = function(file) write.csv(core_table(), file, row.names = FALSE)
  )
}

shinyApp(ui = app_ui, server = app_server)
