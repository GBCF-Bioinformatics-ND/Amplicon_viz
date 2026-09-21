# Functional Profiler: PICRUSt2

**Functional differential analysis and pathway visualization**

## Launch

```bash
shiny::runApp()
# or from terminal:
# Rscript -e "shiny::runApp()"
```

## Analyses Available

### Data Summary
- Total samples and taxa in phyloseq object
- Total sequencing reads
- Available taxonomic ranks
- Sample metadata preview
- Taxonomy information

### Functional Profiler (ggpicrust2)
- Analyzes metabolic/functional pathways from PICRUSt2 predictions
- **Pathway Types**:
  - **KO** (KEGG Ortholog): Individual genes/proteins
  - **MetaCyc**: Metabolic pathways
  - **EC** (Enzyme Commission): Enzyme functions
- **DAA Methods** (pathway-level differential abundance analysis):
  - LinDA (default, conservative)
  - ALDEx2 (compositional data specialist)
  - edgeR (default; RNA-seq inspired)
- **Features**:
  - PCA to visualize functional profiles
  - Errorbar plots and heatmaps of top DAA-ranked pathways
  - KO-to-KEGG pathway annotation
- **Output**: Annotated results table + publication-ready plots

## Sidebar Controls

- **Assay type**: 16S or ITS data
- **Kingdom filter**: All, Bacteria, or Fungi  
- **Input mode**:
  - CSV tables
  - BIOM file
  - QIIME2 feature table
- **Use defaults**: Optional; loads pre-configured files when selected
- **Maximum features for analysis**: Defaults to 2,000 to limit memory use
- **Top features for errorbar and heatmap**: Defaults to 10 DAA-ranked features

## Step-by-Step Workflow

### 1. Load Count Data
1. Choose input mode (CSV, BIOM, or QIIME2)
2. Upload OTU/taxonomy/metadata files
3. Click "Load / Reload data"
4. Check Data Summary tab for row/column counts

### 2. ggpicrust2 Analysis
1. Go to "Functional Profiler" tab
2. Upload PICRUSt2 abundance table (pred_metagenome_unstrat_descrip.tsv or similar)
3. Select grouping variable (same metadata column)
4. (Optional) Select reference level
5. Choose pathway type: KO (default), MetaCyc, or EC
6. Choose DAA method (default: LinDA)
7. Choose a DAA method; `edgeR` is selected by default
8. (Optional) Enable KO-to-KEGG conversion for KEGG pathway names
9. Set the maximum analysis features and the number of plotted top features
10. Click "Run ggpicrust2"
9. Review results:
   - Table shows annotated differential results
   - Errorbar plot shows significant pathways with effect size
   - PCA shows functional profile separation
   - Heatmap shows significant pathway abundance patterns
10. Download desired tables/plots

## Data Format Requirements

### OTU Table (CSV)
```
,Sample1,Sample2,Sample3
OTU1,100,50,75
OTU2,200,150,100
```
- First column = feature IDs
- Other columns = samples
- Values = integer counts

### Taxonomy Table (CSV)
- First column = feature IDs (matching OTU table)
- Other columns = Kingdom, Phylum, Class, Order, Family, Genus, Species
- Can include Taxon and Confidence columns for QIIME2 data

### Sample Metadata (CSV)
```
SampleID,Location,Treatment,Timepoint
Sample1,Gut,Antibiotic,Day1
Sample2,Gut,Placebo,Day1
Sample3,Gut,Antibiotic,Day7
```
- First column = sample IDs (matching OTU table)
- Other columns = grouping variables, continuous variables, etc.

### PICRUSt2 Abundance Table (TSV/TXT)
**Example structure**:
```
function	Sample1	Sample2	Sample3	description
K00001	150	200	180	pyrophosphate phospho-hydrolase
K00002	100	80	120	phosphoglycerate kinase
```
- First column: Feature IDs (KO IDs, MetaCyc pathways, or EC numbers)
- Middle columns: Sample abundances (numeric, decimal OK)
- Last column (optional): Annotations/descriptions
- Row names should be KO/MetaCyc/EC IDs

**Obtaining PICRUSt2 output**:
```bash
# QIIME2 plugin command
qiime picrust2 full-pipeline ...
# Outputs pred_metagenome_unstrat_descrip.tsv, etc.
```

## Interpretation Guide

### ggpicrust2 Results

#### Errorbar Plot
- X-axis: Pathway name or ID
- Y-axis: Group
- Error bars: 95% confidence interval
- Color: Effect direction and magnitude
- **Look for**: Non-overlapping error bars between groups

#### PCA Plot
- Each point = one sample
- Closeness = functional profile similarity
- Separation = distinct functional profiles between groups

#### Heatmap
- Rows: Selected top DAA-ranked pathways
- Columns: Samples
- Color intensity: Abundance
- Hierarchical clustering groups similar samples/pathways

#### Statistical Table
- Feature/Pathway ID
- Group1/Group2: Abundance in each group
- log2FoldChange: Effect size
- p_value, p_adjust: Statistical significance
- Description: Pathway annotation (if available)

## Key Concepts

**Differential Abundance Analysis (DAA)**: Testing whether predicted pathways differ in abundance between groups

**DAA (Differential Abundance Analysis)**: Statistical method for this test

**ggpicrust2**: Applies DAA methods to functional pathways

**KO (KEGG Ortholog)**: Functional gene/protein identifier

**PICRUSt2**: "Prediction of Metabolic Intermediate Genes by Functional Tools"
- Predicts metabolic genes/pathways from 16S/ITS data

**Fold-Change**: Ratio of abundances between groups

## Troubleshooting

**"fewer than 2 overlapping sample IDs"**
- Check that PICRUSt2 table column names exactly match metadata sample IDs
- Case-sensitive!
- Remove whitespace if present

**"ggpicrust2 returned no differential results"**
- May need stronger biological effect or different DAA method
- Try more samples or higher replication
- Consider filtering to more abundant pathways

**"Selected reference level is not present"**
- The reference level you chose doesn't exist in filtered data
- Check for empty values in grouping variable
- Use "Auto" to let the app choose

**"KO to KEGG conversion requires KO IDs"**
- Your table doesn't have proper KO ID format (K##### or ko:K#####)
- Check PICRUSt2 output file; may need different prediction level
- Disable KO-to-KEGG conversion to proceed without KEGG annotations

**Plots not rendering**
- Ensure data is loaded (check Data Summary tab)
- PICRUSt2 file must have numeric abundance columns
- Sample IDs must match between OTU and PICRUSt2 tables

## Advanced Features

### Multiple Contrasts
- If ≥3 groups, app auto-detects contrasts
- "Contrast to plot" selector appears
- Shows most significant contrast by default
- Manual selection available

### Statistical Methods
- **LinDA**: Conservative; recommended for small sample sizes
- **ALDEx2**: Compositional specialist; good for 16S data
- **edgeR**: Default; fast and suitable for larger feature tables

### KO-to-KEGG Conversion
- Maps individual KO IDs → functional modules/pathways
- Provides more interpretable pathway names
- Only works with KO input
- Requires ko_to_kegg_reference data (auto-loaded)

## Output Quality

All plots saved as:
- **Format**: PNG, 300 dpi
- **Size**: Publication-ready
- **Resolution**: High quality for papers/presentations

All tables saved as:
- **Format**: CSV
- **Compatibility**: Excel, R, Python, etc.

## Citations

This app uses:
- **ggpicrust2**: Liang et al. (2023)
- **PICRUSt2**: Douglas et al. (2020) mBio

For DAA methods:
- **LinDA**: Zeller et al. (2021) bioRxiv
- **ALDEx2**: Fernandes et al. (2014) PLoS ONE
- **edgeR**: Robinson et al. (2010) Bioinformatics

---

**Last Updated**: April 2026
