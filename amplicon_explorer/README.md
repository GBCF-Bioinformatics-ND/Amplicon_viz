# Amplicon Explorer: 16S / ITS Analysis App

**Exploratory microbiome analysis focused on taxonomic diversity and community structure**

## Launch

```bash
shiny::runApp()
# or from terminal:
# Rscript -e "shiny::runApp()"
```

## Analyses Available

### Guide Tab
- How to use this app
- Recommended analysis workflow
- Important concepts explained
- Good practices for interpretation

### Data Summary
- Total samples and taxa in dataset
- Total sequencing reads
- Available taxonomic ranks
- Sample metadata preview
- Taxonomy preview with confidence scores

### Alpha Diversity
- **Observed**: Count of taxa detected
- **Shannon**: Accounts for both richness and evenness
- **Simpson**: Emphasizes dominant taxa
- Box plots + individual points per sample
- Tests whether within-sample diversity differs between groups

### Beta Ordination
- **NMDS** (default) or **PCoA**
- **Bray-Curtis**: Uses abundance
- **Jaccard**: Presence/absence only
- **Euclidean**: Standard distance
- Visualizes community similarity between samples

### Taxa Composition
- Stacked bar charts
- Colors = taxa at chosen rank (Phylum, Genus, etc.)
- Top N taxa shown; others grouped as "Other"
- Faceting by second grouping variable available
- Great for spotting dominant taxa and group patterns

### Clustering
- UPGMA dendrograms
- Shows hierarchical similarity between samples
- Samples joining low = more similar
- Helps identify natural groupings or outliers

### PERMANOVA
- Tests overall community composition differences
- R² = fraction of variance explained by grouping
- p-value = evidence against null hypothesis
- Complement with ordination plots for interpretation

### Core Microbiome
- Taxa found across many samples (prevalence filter)
- Heatmap shows which core taxa appear where
- Adjustable abundance threshold
- Useful for identifying "keystone" taxa

## Sidebar Controls

- **Assay type**: 16S or ITS data
- **Kingdom filter**: All, Bacteria, or Fungi
- **Input mode**: 
  - CSV tables (OTU, taxonomy, metadata)
  - BIOM file
  - QIIME2 feature table
- **Use defaults**: Optional; loads pre-configured files from workspace when selected
- **Group/color variable**: Which metadata column to use for coloring plots
- **Facet variable**: Optional second grouping for taxa plot
- **Taxonomic rank**: Which level to show (Phylum, Genus, etc.)
- **Top taxa**: How many taxa to label (rest become "Other")

## Data Format Requirements

### OTU Table (CSV)
```
,Sample1,Sample2,Sample3
OTU1,100,50,75
OTU2,200,150,100
```
- First column = OTU/feature IDs
- Other columns = samples
- Values = read counts (integers)

### Taxonomy Table (CSV)
```
,Kingdom,Phylum,Class,Order,Family,Genus,Species
OTU1,Bacteria,Firmicutes,Clostridia,Clostridiales,Lachnospiraceae,Roseburia,
OTU2,Bacteria,Bacteroidetes,Bacteroidia,Bacteroidales,Bacteroidaceae,Bacteroides,
```
- First column = OTU/feature IDs (must match OTU table)
- Columns = taxonomic ranks
- Values = taxon names or empty

### Sample Metadata (CSV)
```
SampleID,Location,Treatment,DaysSick
Sample1,Gut,Placebo,3
Sample2,Gut,Antibiotic,1
Sample3,Skin,Placebo,5
```
- First column = sample IDs (must match OTU table)
- Other columns = any metadata variables
- Values = groups, continuous, or categorical

### BIOM File
- Standard Biological Observation Matrix format
- Can include taxonomy and sample metadata
- Optionally paired with metadata CSV

### QIIME2 Feature Table
- TSV format exported from QIIME2
- Includes optional Taxon column (semicolon-delimited)
- Paired with metadata file

## Interpretation Tips

### Alpha Diversity
- High Shannon = many taxa, evenly distributed
- Low Shannon = few taxa or one dominates
- Compare boxplot shapes across groups

### Beta Ordination
- Tight clusters = homogeneous communities
- Separate clouds = distinct community types
- Overlap = communities share composition

### Taxa Composition
- Look for consistent patterns across replicates
- Sudden shifts suggest important drivers
- Rare taxa in "Other" category = not worth highlighting

### PERMANOVA + Ordination
- **PERMANOVA p < 0.05 + clear separation** = strong evidence
- **PERMANOVA p < 0.05 + high overlap** = significant but small effect
- **PERMANOVA p > 0.05 + clear separation** = high variance within groups

### Core Microbiome
- Taxa appearing across most samples = "true core"
- Taxa in few samples = environment-specific
- Combine with abundance threshold to find keystone taxa

## Key Concepts

**Taxon**: A biological group (Phylum, Genus, Species, etc.)

**Relative Abundance**: Proportions within each sample; allows comparison despite different sequencing depths

**Alpha Diversity**: Diversity within ONE sample

**Beta Diversity**: Diversity BETWEEN samples (community differences)

**Ordination**: 2D plot arranging similar samples close together

**OTU/ASV**: Operational Taxonomic Unit / Amplicon Sequence Variant (the unit being counted)

**Taxonomy**: Classification of each OTU/ASV

**PERMANOVA**: Permutational MANOVA; tests if groups differ in overall composition

## Downloads

All plots and tables can be downloaded as:
- **Plots**: PNG files (300 dpi, publication-ready)
- **Tables**: CSV files (can open in Excel)

## Troubleshooting

**Error: "Fewer than 2 overlapping samples"**
- Check that column names in OTU table match row names in metadata
- Sample IDs are case-sensitive!

**Error: "Kingdom filter found no matching taxa"**
- Check Kingdom column name and values
- May need "Superkingdom" or "Domain" instead

**Plots not appearing**
- Load data first (Data Summary tab should show counts)
- Select a grouping variable with ≥2 non-empty groups
- Ensure at least 2 samples remain after filtering

**Very slow rendering**
- Large datasets (>10,000 taxa) may be slow
- Consider filtering to top abundant taxa first
- Contact instructor if persistent

## Citations

This app uses:
- **phyloseq**: McMurdie & Holmes (2013) PLoS Comp Biol
- **vegan**: Oksanen et al., R package
- **ggplot2**: Wickham (2009) Springer

---

**Last Updated**: April 2026
