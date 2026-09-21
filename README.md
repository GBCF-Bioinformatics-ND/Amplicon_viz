# Amplicon Analysis & Functional Profiling Suite

This project contains two specialized Shiny applications for microbiome analysis.

## Directory Structure

```
Amplicon_viz/
├── amplicon_explorer/          # Taxonomic diversity and community structure
├── functional_profiler/        # PICRUSt2 functional pathway analysis
├── data/                       # Shared input datasets
├── ggpicrust2_cache/           # Cached KO reference data
├── Amplicon_viz.Rproj          # RStudio project file
└── README.md
```

## Quick Start

From the project root:

```r
shiny::runApp("amplicon_explorer")
shiny::runApp("functional_profiler")
```

## Amplicon Explorer

The explorer focuses on taxonomic diversity and community structure.

Analyses include:

- Data summary, sample metadata, and taxonomy previews
- Alpha diversity: Observed, Shannon, and Simpson
- Beta ordination: NMDS or PCoA with Bray-Curtis, Jaccard, or Euclidean distance
- Taxa composition with top-taxa grouping
- Hierarchical clustering
- PERMANOVA
- Core microbiome prevalence and abundance filtering

It accepts CSV tables, BIOM files, and QIIME2 feature tables. Rarefaction and
taxonomic differential-abundance analysis are not included.

## Functional Profiler

The functional profiler analyzes PICRUSt2 predicted functions with ggpicrust2.

Supported pathway types:

- KO (KEGG Ortholog)
- MetaCyc
- EC (Enzyme Commission)

DAA methods include LinDA, ALDEx2, and edgeR. `edgeR` is selected by default.
KO-to-KEGG conversion is optional and disabled by default.

The default maximum analysis size is 2,000 features. Errorbar and heatmap plots
use the top 10 DAA-ranked features by default; both settings can be adjusted in
the app.

## Shared Data Requirements

- OTU/ASV table: feature IDs in the first column and sample IDs in the remaining columns
- Taxonomy table: feature IDs matching the OTU/ASV table
- Sample metadata: sample IDs matching the count-table columns plus grouping variables
- PICRUSt2 abundance table: functional IDs and numeric sample-abundance columns

Sample IDs must match exactly after normalization where applicable. Grouping
variables should contain at least two groups, with replication preferred for DAA.

## Shared Resources

- `data/otu_table.csv`
- `data/taxonomy_table.csv`
- `data/sample_metadata.csv`
- `data/centrifuge_reports.biom`
- `ggpicrust2_cache/ko_reference.rds`

## Troubleshooting

**Fewer than two overlapping samples**: Check sample IDs across the count table,
metadata, and PICRUSt2 table.

**No results from ggpicrust2**: Choose a grouping variable with at least two
groups and replicated samples, then try a different DAA method or feature filter.

**Missing package**: Install CRAN packages with `install.packages()` and
Bioconductor packages with `BiocManager::install()`.

**Created**: April 2026  
**Framework**: R Shiny  
