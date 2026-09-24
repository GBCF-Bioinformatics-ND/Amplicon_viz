# Amplicon Explorer: 16S / ITS Analysis App

**Exploratory amplicon analysis focused on diversity, community structure, and taxonomy**

## Launch

```bash
shiny::runApp()
# or from terminal:
# Rscript -e "shiny::runApp()"
```

## Analyses Available

### Guide Tab
- Current workflow for count tables, taxonomy, and metadata
- Explanation of each analysis and figure
- Guidance for interpreting effect size, replication, and p-values

### Data Summary
- Total samples and features in dataset
- Total sequencing reads
- Available taxonomic ranks
- Sample metadata preview
- Taxonomy preview with confidence scores

### Alpha Diversity
- **Observed**: Count of taxa detected
- **Shannon**: Accounts for both richness and evenness
- **Simpson**: Emphasizes dominant taxa
- Box plots + individual points per sample
- Helps compare within-sample diversity between groups; it does not identify which taxa cause a difference

### Beta Ordination
- **NMDS** (default) or **PCoA**
- **Bray-Curtis**: Uses abundance
- **Jaccard**: Presence/absence only
- **Euclidean**: Standard distance
- Visualizes community similarity between samples; axis values are not individual taxa or direct abundance values

### Taxa Composition
- Stacked bar charts
- Colors = taxa at chosen rank (Phylum, Genus, etc.)
- Top N taxa shown; others grouped as "Other"
- Faceting by second grouping variable available
- Useful for spotting dominant taxa and broad relative-abundance patterns, not absolute abundance

### Clustering
- UPGMA dendrograms
- Shows hierarchical similarity between samples
- Samples joining low = more similar
- Helps identify sample groupings or outliers; results depend on distance and linkage choices

### PERMANOVA
- Tests overall community composition differences
- R² = fraction of variance explained by grouping
- p-value = evidence against the permutation null model
- Interpret R2 as explained variation, not classification accuracy
- Complement with ordination plots and checks for unequal within-group dispersion

### Core Microbiome
- Taxa found across many samples (prevalence filter)
- Heatmap shows which core taxa appear where
- Adjustable abundance threshold
- Does not prove that a taxon is biologically essential or a "keystone" taxon

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
- Compare group medians, spread, individual points, and sample sizes

### Beta Ordination
- Tight clusters = homogeneous communities
- Separate clouds = distinct community types
- Overlap = communities share composition or within-group variation is high
- Use PERMANOVA to test group differences; do not infer significance from visual separation alone

### Taxa Composition
- Look for consistent patterns across replicates
- Shifts suggest differences worth investigating, but do not identify a causal driver
- Rare taxa grouped as "Other" remain part of the total relative abundance

### PERMANOVA + Ordination
- **Small p-value + meaningful R2 + consistent separation** = stronger evidence of a group-associated community difference
- **Small p-value + low R2** = statistically detectable but potentially modest effect
- **Visual separation without a small p-value** = exploratory pattern requiring more replication or a closer look at within-group variation
- **Different within-group spread** = interpret PERMANOVA cautiously because dispersion can influence the result

### Core Microbiome
- Taxa appearing across most samples meet the selected prevalence definition of core
- Taxa in fewer samples may be condition-specific or low-prevalence
- Combine prevalence with the abundance threshold and inspect the underlying table

## Key Concepts

**Feature**: An OTU or ASV counted in the input table

**Taxon**: A biological group assigned to one or more features (Phylum, Genus, Species, etc.)

**Relative Abundance**: Proportions within each sample; allows comparison despite different sequencing depths

**Alpha Diversity**: Diversity within ONE sample

**Beta Diversity**: Diversity BETWEEN samples (community differences)

**Ordination**: A 2D representation of sample-to-sample distances; proximity is meaningful, while axis values are usually not directly interpretable

**OTU/ASV**: Operational Taxonomic Unit / Amplicon Sequence Variant, the feature being counted

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

**Last Updated**: September 2026
