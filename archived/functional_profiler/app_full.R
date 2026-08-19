# Created by: Bharat Mishra
# Edited by: Elizabeth Brooks
# Modified: 15 June 2026

# install any missing packages
options(repos = c(CRAN = "https://cloud.r-project.org/"))
biocList <- c("DESeq2", "phyloseq", "ALDEx2", "edgeR", "limma")
packageList <- c("shiny", "ggplot2", "vegan", "ggpicrust2", "RcppEigen", "RcppParallel")
newBioc <- biocList[!(biocList %in% installed.packages()[,"Package"])]
newPackages <- packageList[!(packageList %in% installed.packages()[,"Package"])]
if(length(newBioc)){
  install.packages("BiocManager")
  BiocManager::install(newBioc)
}
if(length(newPackages)){
  install.packages(newPackages)
}
suppressPackageStartupMessages({
  library(shiny)
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(ggpicrust2)
})

options(shiny.maxRequestSize = 500 * 1024^2)

rank_names_default <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

default_paths <- list(
  otu = "data/otu_table.csv",
  tax = "data/taxonomy_table.csv",
  meta = "data/sample_metadata.csv",
  biom = "data/centrifuge_reports.biom"
)

ggpicrust_cache_dir <- file.path(getwd(), "ggpicrust2_cache")
ko_reference_rds_path <- file.path(ggpicrust_cache_dir, "ko_reference.rds")

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
  titlePanel("Functional Profiler: DESeq2 & PICRUSt2"),
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
      checkboxInput("use_defaults", "Use workspace default files", value = TRUE),
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
      uiOutput("group_var_ui")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel(
          "Guide",
          tags$br(),
          tags$h3("How To Use This App"),
          tags$p("This app is designed for exploratory microbiome analysis of 16S or ITS data. Each tab answers a different biological question. Start with the data summary, then move from simple description to more advanced comparison tests."),
          analysis_help_block(
            "Recommended workflow",
            list(
              "Load your count table, taxonomy table, and sample metadata, or use a BIOM file.",
              "Check the Data summary tab first to make sure the number of samples and taxa look reasonable.",
              "Use Alpha diversity to ask whether some groups have more within-sample diversity than others.",
              "Use Beta ordination and Clustering to ask whether whole communities differ among groups.",
              "Use Taxa composition to see which taxa dominate samples or groups.",
              "Use PERMANOVA and Differential abundance only after checking the exploratory plots, because statistical output is easier to interpret when you already understand the data pattern.",
              "Use Core microbiome to identify taxa that are common across many samples."
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
              "A p-value is evidence against a null hypothesis, but biological interpretation should not rely on p-values alone.",
              "Adjusted p-values control for multiple testing and are more appropriate when many taxa or pathways are tested at once."
            )
          ),
          analysis_help_block(
            "Good practice for interpretation",
            list(
              "Look for consistent patterns across multiple plots instead of relying on one result.",
              "Always interpret statistics together with effect size, group separation, and sample size.",
              "Be careful when groups have very different numbers of samples.",
              "A statistically significant result does not always mean a biologically large effect.",
              "Taxonomic labels can be incomplete; an 'Unassigned' label does not necessarily mean the feature is unimportant."
            )
          )
        ),
        tabPanel(
          "Data summary",
          tags$br(),
          analysis_help_block(
            "What this tab shows",
            list(
              "The total number of samples and taxa in the current dataset.",
              "The total number of reads after filtering.",
              "The taxonomic ranks available in your taxonomy table.",
              "A preview of your sample metadata, which is used to color or group many of the plots.",
              "A taxonomy preview table with Taxon/taxonomy and Confidence columns when available (for QIIME2 and compatible inputs)."
            )
          ),
          verbatimTextOutput("summary_text"),
          tags$h4("Sample metadata preview"),
          tableOutput("sample_table"),
          downloadButton("download_sample_table", "Download sample table"),
          tags$br(),
          tags$br(),
          tags$h4("Taxonomy preview (includes Taxon and Confidence when available)"),
          tableOutput("taxonomy_table"),
          downloadButton("download_taxonomy_table", "Download taxonomy table")
        ),
        tabPanel(
          "Differential abundance",
          tags$br(),
          analysis_help_block(
            "How to read differential abundance",
            list(
              "This analysis tests whether taxa differ in abundance between two groups.",
              "Log2 fold change shows the direction and size of the difference between the comparison group and the reference group.",
              "Positive values mean higher abundance in the comparison group; negative values mean higher abundance in the reference group.",
              "The volcano plot combines effect size and statistical evidence.",
              "Focus on taxa with both meaningful fold change and low adjusted p-value."
            )
          ),
          uiOutput("da_var_ui"),
          fluidRow(
            column(6, uiOutput("da_ref_level_ui")),
            column(6, uiOutput("da_comp_level_ui"))
          ),
          actionButton("run_deseq", "Run DESeq2"),
          tags$br(),
          tags$br(),
          verbatimTextOutput("deseq_status"),
          tableOutput("deseq_table"),
          plotOutput("volcano_plot", height = 440),
          downloadButton("download_deseq_table", "Download DESeq2 table"),
          downloadButton("download_volcano_plot", "Download volcano plot")
        ),
        tabPanel(
          "ggpicrust2",
          tags$br(),
          analysis_help_block(
            "How to use ggpicrust2 in this app",
            list(
              "Upload a PICRUSt2 functional abundance table (for example, pred_metagenome_unstrat_descrip.tsv).",
              "Use the same sample IDs as your metadata so the app can match columns to samples.",
              "Choose a grouping variable, then run differential functional analysis.",
              "The table shows annotated differential results, and plots show effect patterns and clustering."
            )
          ),
          fileInput("picrust_file", "PICRUSt2 abundance table (txt/tsv)", accept = c(".txt", ".tsv")),
          uiOutput("picrust_group_var_ui"),
          uiOutput("picrust_reference_ui"),
          uiOutput("picrust_contrast_ui"),
          fluidRow(
            column(
              4,
              selectInput(
                "picrust_daa_method",
                "DAA method",
                choices = c("LinDA", "ALDEx2", "DESeq2", "edgeR"),
                selected = "LinDA"
              )
            ),
            column(
              4,
              selectInput(
                "picrust_pathway",
                "Pathway type",
                choices = c("KO", "MetaCyc", "EC"),
                selected = "KO"
              )
            ),
            column(
              4,
              checkboxInput("picrust_ko_to_kegg", "Convert KO to KEGG pathways", value = TRUE)
            )
          ),
          actionButton("run_picrust", "Run ggpicrust2", class = "btn-primary"),
          tags$br(),
          tags$br(),
          verbatimTextOutput("picrust_status"),
          tableOutput("picrust_table"),
          downloadButton("download_picrust_table", "Download ggpicrust2 table"),
          tags$br(),
          tags$br(),
          plotOutput("picrust_errorbar_plot", height = 500),
          downloadButton("download_picrust_errorbar_plot", "Download errorbar plot"),
          tags$br(),
          tags$br(),
          plotOutput("picrust_pca_plot", height = 500),
          downloadButton("download_picrust_pca_plot", "Download PCA plot"),
          tags$br(),
          tags$br(),
          plotOutput("picrust_heatmap_plot", height = 520),
          downloadButton("download_picrust_heatmap_plot", "Download heatmap plot")
        )
      )
    )
  )
)

app_server <- function(input, output, session) {
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

  get_picrust_file <- reactive({
    if (!is.null(input$picrust_file) && nzchar(input$picrust_file$datapath)) {
      return(input$picrust_file$datapath)
    }

    default_picrust <- "qiimeandpicrust_oyester/pred_metagenome_unstrat_descrip.tsv"
    if (isTRUE(input$use_defaults) && file.exists(default_picrust)) {
      return(default_picrust)
    }

    stop("Upload a PICRUSt2 abundance table in the ggpicrust2 tab.")
  })

  read_picrust_abundance <- function(path) {
    ab <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(ab) < 3) {
      stop("PICRUSt2 table appears malformed: expected feature column plus at least 2 sample columns.")
    }

    feature_candidates <- c("function", "#NAME", "pathway", "KO", "feature")
    feature_col <- feature_candidates[feature_candidates %in% colnames(ab)]
    feature_col <- if (length(feature_col) > 0) feature_col[1] else colnames(ab)[1]

    ab[[feature_col]] <- trimws(as.character(ab[[feature_col]]))
    ab <- ab[!is.na(ab[[feature_col]]) & ab[[feature_col]] != "", , drop = FALSE]
    rownames(ab) <- make.unique(ab[[feature_col]])
    ab <- ab[, setdiff(colnames(ab), feature_col), drop = FALSE]

    is_numeric_like_col <- vapply(ab, function(col) {
      vals <- trimws(as.character(col))
      vals <- vals[!is.na(vals)]
      vals <- vals[!vals %in% c("", "NA", "N/A", "na", "n/a", "NULL", "null")]
      if (length(vals) == 0) {
        return(FALSE)
      }
      all(grepl("^-?[0-9]+(\\.[0-9]+)?$", vals))
    }, logical(1))

    ab <- ab[, is_numeric_like_col, drop = FALSE]
    if (ncol(ab) < 2) {
      stop("PICRUSt2 table has fewer than 2 numeric sample columns after filtering.")
    }

    sanitize_count_table(ab, table_label = "PICRUSt2 table")
  }

  has_ko_to_kegg_reference <- function() {
    out <- tryCatch(
      {
        data("ko_to_kegg_reference", package = "ggpicrust2", envir = environment())
        exists("ko_to_kegg_reference", inherits = TRUE)
      },
      error = function(e) FALSE
    )
    isTRUE(out)
  }

  ensure_ko_reference_rds <- function(cache_path = ko_reference_rds_path) {
    if (file.exists(cache_path)) {
      return(cache_path)
    }

    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)

    ns <- asNamespace("ggpicrust2")
    ko_ref <- NULL

    if (exists("ko_reference", envir = ns, inherits = FALSE)) {
      ko_ref <- get("ko_reference", envir = ns, inherits = FALSE)
    } else if (exists("KO_reference", envir = ns, inherits = FALSE)) {
      ko_ref <- get("KO_reference", envir = ns, inherits = FALSE)
    } else {
      ref_env <- new.env(parent = emptyenv())
      loaded <- tryCatch(
        {
          utils::data("ko_reference", package = "ggpicrust2", envir = ref_env)
          exists("ko_reference", envir = ref_env, inherits = FALSE)
        },
        error = function(e) FALSE
      )

      if (!loaded) {
        loaded <- tryCatch(
          {
            utils::data("KO_reference", package = "ggpicrust2", envir = ref_env)
            exists("KO_reference", envir = ref_env, inherits = FALSE)
          },
          error = function(e) FALSE
        )
      }

      if (!loaded) {
        stop(
          "Could not load KO reference from ggpicrust2 (neither 'ko_reference' nor 'KO_reference'). ",
          "Please reinstall/update ggpicrust2."
        )
      }

      ko_ref_name <- if (exists("ko_reference", envir = ref_env, inherits = FALSE)) "ko_reference" else "KO_reference"
      ko_ref <- get(ko_ref_name, envir = ref_env, inherits = FALSE)
    }

    saveRDS(ko_ref, cache_path)
    cache_path
  }

  ensure_ko_reference_loaded <- function(cache_path = ko_reference_rds_path) {
    cache_file <- ensure_ko_reference_rds(cache_path)
    ko_ref <- readRDS(cache_file)

    if (!is.data.frame(ko_ref) || nrow(ko_ref) == 0) {
      stop("Cached ko_reference.rds is invalid or empty.")
    }

    ns <- asNamespace("ggpicrust2")
    has_binding <- exists("ko_reference", envir = ns, inherits = FALSE)
    if (!has_binding) {
      assign("ko_reference", ko_ref, envir = ns)
    }

    has_caps_binding <- exists("KO_reference", envir = ns, inherits = FALSE)
    if (!has_caps_binding) {
      assign("KO_reference", ko_ref, envir = ns)
    }

    invisible(cache_file)
  }

  make_notice_plot <- function(label) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = label, size = 5) +
      xlim(0, 1) +
      ylim(0, 1) +
      theme_void()
  }



  ordination_obj <- reactive({
    ps_rel_local <- ps_rel()
    if (identical(input$ordination_method, "NMDS")) {
      ordinate(ps_rel_local, method = "NMDS", distance = input$distance_method)
    } else {
      ordinate(ps_rel_local, method = "PCoA", distance = input$distance_method)
    }
  })

  ordination_df <- reactive({
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
    plot_ordination(ps_rel_local, ord, color = input$group_var, justDF = TRUE)
  })

  ordination_plot_obj <- reactive({
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
    plot_ordination(ps_rel_local, ord, color = input$group_var) +
      geom_point(size = 3, alpha = 0.9) +
      theme_bw(base_size = 13)
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

    output$picrust_group_var_ui <- renderUI({
      selectInput("picrust_group_var", "ggpicrust2 group variable", choices = meta_cols, selected = preferred_group)
    })

    output$picrust_reference_ui <- renderUI({
      req(input$picrust_group_var)
      levs <- unique(trimws(as.character(meta_df[[input$picrust_group_var]])))
      levs <- levs[!is.na(levs) & levs != ""]
      choices <- c("Auto (first level)" = "", levs)
      selectInput("picrust_reference", "ggpicrust2 reference level (optional)", choices = choices, selected = "")
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
    cat("Amplicon object loaded successfully\n")
    cat(sprintf("Samples: %d\n", nsamples(ps)))
    cat(sprintf("Taxa: %d\n", ntaxa(ps)))
    cat(sprintf("Assay type: %s\n", unique(sample_data(ps)$AssayType)))
    cat(sprintf("Total reads: %s\n", format(sum(sample_sums(ps)), big.mark = ",")))
    cat("\nTaxonomic ranks:\n")
    print(rank_names(ps))
  })

  output$sample_table <- renderTable({
    sample_table_df()
  }, rownames = TRUE)

  output$taxonomy_table <- renderTable({
    head(taxonomy_table_df(), 25)
  }, rownames = FALSE)






  deseq_results <- eventReactive(input$run_deseq, {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      stop("Package DESeq2 is not installed.")
    }

    req(input$da_var, input$da_ref_level, input$da_comp_level)
    if (identical(input$da_ref_level, input$da_comp_level)) {
      stop("Reference and comparison levels must be different.")
    }

    ps <- ps_obj()
    meta_df <- as(sample_data(ps), "data.frame")
    req(input$da_var %in% colnames(meta_df))

    keep_samples <- rownames(meta_df)[as.character(meta_df[[input$da_var]]) %in% c(input$da_ref_level, input$da_comp_level)]
    ps_sub <- prune_samples(keep_samples, ps)
    meta_sub <- as(sample_data(ps_sub), "data.frame")
    meta_sub[[input$da_var]] <- factor(as.character(meta_sub[[input$da_var]]), levels = c(input$da_ref_level, input$da_comp_level))
    sample_data(ps_sub) <- sample_data(meta_sub)

    dds <- phyloseq_to_deseq2(ps_sub, stats::as.formula(paste("~", input$da_var)))
    dds <- DESeq2::DESeq(dds, fitType = "parametric", quiet = TRUE)
    res <- DESeq2::results(dds, contrast = c(input$da_var, input$da_comp_level, input$da_ref_level))

    res_df <- as.data.frame(res)
    res_df$Taxon <- rownames(res_df)

    if (!is.null(tax_table(ps_sub, errorIfNULL = FALSE))) {
      tax_df <- as.data.frame(tax_table(ps_sub), stringsAsFactors = FALSE)
      tax_df$Taxon <- rownames(tax_df)
      res_df <- merge(res_df, tax_df, by = "Taxon", all.x = TRUE, sort = FALSE)
    }

    res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
    res_df
  })

  picrust_base_results <- eventReactive(input$run_picrust, {
    if (!requireNamespace("ggpicrust2", quietly = TRUE)) {
      stop("Package ggpicrust2 is not installed. Install with install.packages('ggpicrust2').")
    }

    ko_ref_cache_used <- tryCatch(
      {
        ensure_ko_reference_loaded()
        TRUE
      },
      error = function(e) {
        FALSE
      }
    )

    req(input$picrust_group_var)
    picrust_path <- get_picrust_file()
    abundance <- read_picrust_abundance(picrust_path)

    meta <- as(sample_data(ps_obj()), "data.frame")
    meta$sample_name <- rownames(meta)
    req(input$picrust_group_var %in% colnames(meta))

    colnames(abundance) <- trimws(colnames(abundance))
    common_samples <- intersect(colnames(abundance), meta$sample_name)
    if (length(common_samples) < 2) {
      stop("Fewer than 2 overlapping sample IDs between PICRUSt2 table and loaded metadata.")
    }

    abundance <- abundance[, common_samples, drop = FALSE]
    meta <- meta[match(common_samples, meta$sample_name), , drop = FALSE]

    group_vals <- trimws(as.character(meta[[input$picrust_group_var]]))
    keep <- !is.na(group_vals) & group_vals != ""
    abundance <- abundance[, keep, drop = FALSE]
    meta <- meta[keep, , drop = FALSE]

    if (ncol(abundance) < 2) {
      stop("Need at least 2 samples with non-missing group values for ggpicrust2.")
    }

    meta_gg <- data.frame(
      sample_name = meta$sample_name,
      GroupVar = as.character(meta[[input$picrust_group_var]]),
      stringsAsFactors = FALSE
    )

    requested_reference <- NULL
    if (!is.null(input$picrust_reference) && nzchar(input$picrust_reference)) {
      requested_reference <- input$picrust_reference
    }

    if (!is.null(requested_reference) && !(requested_reference %in% unique(meta_gg$GroupVar))) {
      stop("Selected ggpicrust2 reference level is not present after sample filtering.")
    }

    if (length(unique(meta_gg$GroupVar)) < 2) {
      stop("Need at least 2 groups after removing missing values for ggpicrust2.")
    }

    grp_counts <- table(meta_gg$GroupVar)
    if (!any(grp_counts >= 2)) {
      stop(
        paste0(
          "Selected grouping variable has no replicated groups for ggpicrust2. ",
          "Choose a variable where at least one group has >= 2 samples."
        )
      )
    }

    analysis_abundance <- abundance
    notes <- character(0)
    effective_ko_to_kegg <- identical(input$picrust_pathway, "KO") && isTRUE(input$picrust_ko_to_kegg)
    if (effective_ko_to_kegg && !has_ko_to_kegg_reference()) {
      notes <- c(
        notes,
        "Reference 'ko_to_kegg_reference' not found; running without KO-to-KEGG conversion."
      )
      effective_ko_to_kegg <- FALSE
    }

    if (isTRUE(ko_ref_cache_used)) {
      notes <- c(notes, paste0("Loaded local KO reference cache: ", ko_reference_rds_path))
    } else {
      notes <- c(
        notes,
        "Could not load local KO reference cache; ggpicrust2 may require package reinstallation if KO annotation fails."
      )
    }

    if (effective_ko_to_kegg) {
      ko_ids <- trimws(rownames(abundance))
      ko_like <- grepl("^(ko:)?K[0-9]{5}$", ko_ids)
      if (!all(ko_like)) {
        stop("KO to KEGG conversion requires KO IDs as row names (K##### or ko:K#####).")
      }

      ko_input <- data.frame(
        ko_ids,
        abundance,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      colnames(ko_input)[1] <- "function"

      analysis_abundance <- ggpicrust2::ko2kegg_abundance(data = ko_input)
    }

    daa_results <- ggpicrust2::pathway_daa(
      abundance = analysis_abundance,
      metadata = meta_gg,
      group = "GroupVar",
      daa_method = input$picrust_daa_method,
      p_adjust_method = "BH",
      reference = requested_reference
    )

    if (is.null(daa_results) || nrow(daa_results) == 0) {
      stop(
        paste0(
          "ggpicrust2 returned no differential results (empty table). ",
          "Try a grouping variable with stronger replication or switch DAA method."
        )
      )
    }

    annotated_results <- ggpicrust2::pathway_annotation(
      pathway = input$picrust_pathway,
      daa_results_df = daa_results,
      ko_to_kegg = effective_ko_to_kegg
    )

    if (is.null(annotated_results) || nrow(annotated_results) == 0) {
      stop(
        paste0(
          "ggpicrust2 annotation produced an empty table. ",
          "Check grouping variable and pathway settings."
        )
      )
    }

    available_contrasts <- character(0)
    best_contrast <- NULL
    if (all(c("group1", "group2") %in% colnames(annotated_results))) {
      contrast_key <- paste(annotated_results$group1, annotated_results$group2, sep = "__vs__")
      split_idx <- split(seq_len(nrow(annotated_results)), contrast_key)
      contrast_scores <- vapply(split_idx, function(idx) {
        pvals <- annotated_results$p_adjust[idx]
        pvals <- pvals[!is.na(pvals)]
        if (length(pvals) == 0) {
          return(Inf)
        }
        min(pvals)
      }, numeric(1))

      available_contrasts <- names(split_idx)
      best_key <- names(which.min(contrast_scores))[1]
      if (!is.na(best_key) && nzchar(best_key)) {
        best_contrast <- best_key
      }
    }

    pca_plot <- tryCatch(
      ggpicrust2::pathway_pca(
        abundance = analysis_abundance,
        metadata = meta_gg,
        group = "GroupVar"
      ),
      error = function(e) {
        notes <<- c(notes, paste0("PCA plot unavailable: ", conditionMessage(e)))
        make_notice_plot("PCA plot unavailable for current settings")
      }
    )

    list(
      abundance = analysis_abundance,
      metadata = meta_gg,
      daa = daa_results,
      annotated = annotated_results,
      pca = pca_plot,
      n_samples = ncol(analysis_abundance),
      n_features = nrow(analysis_abundance),
      notes = notes,
      reference = requested_reference,
      effective_ko_to_kegg = effective_ko_to_kegg,
      available_contrasts = available_contrasts,
      default_contrast = best_contrast
    )
  })

  picrust_results <- reactive({
    req(input$run_picrust > 0)
    base <- picrust_base_results()

    analysis_abundance <- base$abundance
    meta_gg <- base$metadata
    annotated_results <- base$annotated
    notes <- base$notes

    plot_results <- annotated_results
    contrast_label <- NULL
    if (length(base$available_contrasts) > 0 && all(c("group1", "group2") %in% colnames(plot_results))) {
      selected_contrast <- input$picrust_contrast
      if (is.null(selected_contrast) || !nzchar(selected_contrast) || !(selected_contrast %in% base$available_contrasts)) {
        selected_contrast <- base$default_contrast
      }

      if (!is.null(selected_contrast) && nzchar(selected_contrast)) {
        contrast_key <- paste(plot_results$group1, plot_results$group2, sep = "__vs__")
        plot_results <- plot_results[contrast_key == selected_contrast, , drop = FALSE]
        contrast_label <- selected_contrast
      }

      if (length(base$available_contrasts) > 1 && !is.null(contrast_label)) {
        if (!is.null(input$picrust_contrast) && nzchar(input$picrust_contrast)) {
          notes <- c(notes, paste0("Manual contrast selected: ", gsub("__vs__", " vs ", contrast_label), "."))
        } else {
          notes <- c(notes, paste0("Multiple group contrasts detected; plotting contrast: ", gsub("__vs__", " vs ", contrast_label), "."))
        }
      }
    }

    x_lab_value <- if (base$effective_ko_to_kegg) "pathway_name" else "description"
    order_value <- if (base$effective_ko_to_kegg) "pathway_class" else "group"

    errorbar_plot <- tryCatch(
      {
        ggpicrust2::pathway_errorbar(
          abundance = analysis_abundance,
          daa_results_df = plot_results,
          Group = meta_gg$GroupVar,
          ko_to_kegg = base$effective_ko_to_kegg,
          p_values_threshold = 0.05,
          order = order_value,
          select = NULL,
          p_value_bar = TRUE,
          colors = NULL,
          x_lab = x_lab_value
        )
      },
      error = function(e) {
        notes <<- c(notes, paste0("Errorbar plot at FDR<0.05 unavailable: ", conditionMessage(e)))
        tryCatch(
          ggpicrust2::pathway_errorbar(
            abundance = analysis_abundance,
            daa_results_df = plot_results,
            Group = meta_gg$GroupVar,
            ko_to_kegg = base$effective_ko_to_kegg,
            p_values_threshold = 1,
            order = order_value,
            select = NULL,
            p_value_bar = TRUE,
            colors = NULL,
            x_lab = x_lab_value
          ),
          error = function(e2) {
            notes <<- c(notes, paste0("Fallback errorbar plot unavailable: ", conditionMessage(e2)))
            make_notice_plot("No errorbar plot available for current settings")
          }
        )
      }
    )

    sig_features <- plot_results$feature[
      !is.na(plot_results$p_adjust) & plot_results$p_adjust < 0.05
    ]
    sig_features <- intersect(sig_features, rownames(analysis_abundance))

    heatmap_plot <- NULL
    if (length(sig_features) > 0) {
      heatmap_plot <- tryCatch(
        ggpicrust2::pathway_heatmap(
          abundance = analysis_abundance[sig_features, , drop = FALSE],
          metadata = meta_gg,
          group = "GroupVar"
        ),
        error = function(e) {
          notes <<- c(notes, paste0("Heatmap unavailable: ", conditionMessage(e)))
          NULL
        }
      )
    } else {
      notes <- c(notes, "No significant pathways at FDR < 0.05; heatmap not generated.")
    }

    list(
      daa = base$daa,
      annotated = base$annotated,
      errorbar = errorbar_plot,
      pca = base$pca,
      heatmap = heatmap_plot,
      sig_features = sig_features,
      n_samples = base$n_samples,
      n_features = base$n_features,
      notes = notes,
      contrast = contrast_label,
      reference = base$reference,
      available_contrasts = base$available_contrasts,
      default_contrast = base$default_contrast
    )
  })

  output$picrust_contrast_ui <- renderUI({
    if (is.null(input$run_picrust) || input$run_picrust < 1) {
      return(NULL)
    }

    base <- picrust_base_results()
    if (length(base$available_contrasts) <= 1) {
      return(NULL)
    }

    contrast_labels <- gsub("__vs__", " vs ", base$available_contrasts)
    choices <- c("Auto (most significant contrast)" = "", stats::setNames(base$available_contrasts, contrast_labels))

    selected_val <- ""
    if (!is.null(input$picrust_contrast) && input$picrust_contrast %in% base$available_contrasts) {
      selected_val <- input$picrust_contrast
    }

    selectInput("picrust_contrast", "Contrast to plot", choices = choices, selected = selected_val)
  })

  output$deseq_status <- renderPrint({
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      cat("DESeq2 package is not installed. Install with BiocManager::install('DESeq2').\n")
      return(invisible(NULL))
    }
    cat("Ready. Choose variable and levels, then click Run DESeq2.\n")
  })

  output$picrust_status <- renderPrint({
    if (!requireNamespace("ggpicrust2", quietly = TRUE)) {
      cat("ggpicrust2 package is not installed. Install with install.packages('ggpicrust2').\n")
      return(invisible(NULL))
    }

    if (is.null(input$run_picrust) || input$run_picrust < 1) {
      cat("Ready. Upload a PICRUSt2 table in the ggpicrust2 tab and click Run ggpicrust2.\n")
      return(invisible(NULL))
    }

    res <- picrust_results()
    cat("ggpicrust2 run completed.\n")
    cat(sprintf("Samples used: %d\n", res$n_samples))
    cat(sprintf("Features used: %d\n", res$n_features))
    cat(sprintf("DAA rows: %d\n", nrow(res$daa)))
    cat(sprintf("Annotated rows: %d\n", nrow(res$annotated)))
    cat(sprintf("Significant pathways (FDR < 0.05): %d\n", length(res$sig_features)))
    if (!is.null(res$reference)) {
      cat(sprintf("Reference level: %s\n", res$reference))
    }
    if (!is.null(res$contrast)) {
      cat(sprintf("Plotted contrast: %s\n", gsub("__vs__", " vs ", res$contrast)))
    }

    if (length(res$notes) > 0) {
      cat("\nNotes:\n")
      cat(paste0("- ", res$notes, collapse = "\n"), "\n")
    }
  })

  output$deseq_table <- renderTable({
    req(input$run_deseq > 0)
    head(deseq_results(), 25)
  }, rownames = FALSE)

  output$volcano_plot <- renderPlot({
    req(input$run_deseq > 0)
    res_df <- deseq_results()
    req(nrow(res_df) > 0)

    res_df$significant <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05, "FDR < 0.05", "NS")

    ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c("FDR < 0.05" = "red", "NS" = "grey50")) +
      theme_bw(base_size = 12) +
      labs(
        title = "Differential abundance volcano plot",
        x = paste0("log2 fold change (", input$da_comp_level, " vs ", input$da_ref_level, ")"),
        y = "-log10 adjusted p-value",
        color = "Significance"
      )
  })

  output$picrust_table <- renderTable({
    req(input$run_picrust > 0)
    head(picrust_results()$annotated, 25)
  }, rownames = FALSE)

  output$picrust_errorbar_plot <- renderPlot({
    req(input$run_picrust > 0)
    print(picrust_results()$errorbar)
  })

  output$picrust_pca_plot <- renderPlot({
    req(input$run_picrust > 0)
    print(picrust_results()$pca)
  })

  output$picrust_heatmap_plot <- renderPlot({
    req(input$run_picrust > 0)
    validate(need(!is.null(picrust_results()$heatmap), "No significant pathways at FDR < 0.05 to plot heatmap."))
    print(picrust_results()$heatmap)
  })



  volcano_plot_obj <- reactive({
    req(input$run_deseq > 0)
    res_df <- deseq_results()
    req(nrow(res_df) > 0)
    res_df$significant <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05, "FDR < 0.05", "NS")

    ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c("FDR < 0.05" = "red", "NS" = "grey50")) +
      theme_bw(base_size = 12) +
      labs(
        title = "Differential abundance volcano plot",
        x = paste0("log2 fold change (", input$da_comp_level, " vs ", input$da_ref_level, ")"),
        y = "-log10 adjusted p-value",
        color = "Significance"
      )
  })

  output$volcano_plot <- renderPlot({
    volcano_plot_obj()
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

  output$download_rarefaction_plot <- downloadHandler(
    filename = function() paste0("rarefaction_", Sys.Date(), ".png"),
    content = function(file) {
      png(file, width = 1200, height = 800, res = 120)
      draw_rarefaction_plot()
      dev.off()
    }
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
      plot(cluster_obj(), main = "UPGMA clustering (selected distance)", xlab = "", sub = "")
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

  output$download_deseq_table <- downloadHandler(
    filename = function() paste0("deseq2_results_", Sys.Date(), ".csv"),
    content = function(file) write.csv(deseq_results(), file, row.names = FALSE)
  )

  output$download_volcano_plot <- downloadHandler(
    filename = function() paste0("volcano_", Sys.Date(), ".png"),
    content = function(file) ggsave(file, plot = volcano_plot_obj(), width = 8, height = 5, dpi = 300)
  )

  output$download_picrust_table <- downloadHandler(
    filename = function() paste0("ggpicrust2_results_", Sys.Date(), ".csv"),
    content = function(file) write.csv(picrust_results()$annotated, file, row.names = FALSE)
  )

  output$download_picrust_errorbar_plot <- downloadHandler(
    filename = function() paste0("ggpicrust2_errorbar_", Sys.Date(), ".png"),
    content = function(file) {
      png(file, width = 1400, height = 900, res = 140)
      print(picrust_results()$errorbar)
      dev.off()
    }
  )

  output$download_picrust_pca_plot <- downloadHandler(
    filename = function() paste0("ggpicrust2_pca_", Sys.Date(), ".png"),
    content = function(file) {
      png(file, width = 1200, height = 900, res = 140)
      print(picrust_results()$pca)
      dev.off()
    }
  )

  output$download_picrust_heatmap_plot <- downloadHandler(
    filename = function() paste0("ggpicrust2_heatmap_", Sys.Date(), ".png"),
    content = function(file) {
      validate(need(!is.null(picrust_results()$heatmap), "No significant pathways at FDR < 0.05 to save heatmap."))
      png(file, width = 1400, height = 1000, res = 140)
      print(picrust_results()$heatmap)
      dev.off()
    }
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
