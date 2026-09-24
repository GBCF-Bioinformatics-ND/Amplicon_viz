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
- Total samples and PICRUSt2 features
- Sample metadata preview
- PICRUSt2 feature ID preview

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

## Input Controls

- **PICRUSt2 abundance table**: `pred_metagenome_unstrat_descrip.tsv` or similar
- **Sample metadata**: CSV with a sample ID column and grouping variables
- **Use defaults**: Optional; loads the workspace default PICRUSt2 table and metadata
- **Maximum features for analysis**: Defaults to 2,000 to limit memory use
- **Top features for errorbar and heatmap**: Defaults to 10 DAA-ranked features

## Step-by-Step Workflow

### 1. Load PICRUSt2 Data
1. On the Guide tab, upload the PICRUSt2 abundance table and sample metadata
2. Click "Load / Reload data"
3. Check the Data summary tab for sample and feature counts

### 2. ggpicrust2 Analysis
1. Go to the "ggpicrust2" tab
2. Confirm the PICRUSt2 abundance table and metadata were loaded on the Guide tab
3. Select grouping variable (same metadata column)
4. (Optional) Select reference level
5. Choose pathway type: KO (default), MetaCyc, or EC
6. Choose DAA method (default: LinDA)
7. Choose a DAA method; `edgeR` is selected by default
8. (Optional) Enable KO-to-KEGG conversion for KEGG pathway names
9. Set the maximum analysis features and the number of plotted top features
10. Click "Run ggpicrust2"
11. Review results:
   - Table shows annotated differential results
  - Errorbar plot compares selected functions across groups
  - PCA shows similarity and separation among sample functional profiles
  - Heatmap shows abundance patterns for top DAA-ranked functions
12. Download desired tables/plots

## Data Format Requirements

### Sample Metadata (CSV)
```
SampleID,Location,Treatment,Timepoint
Sample1,Gut,Antibiotic,Day1
Sample2,Gut,Placebo,Day1
Sample3,Gut,Antibiotic,Day7
```
- First column = sample IDs matching the PICRUSt2 sample columns
- Other columns = grouping variables and other sample annotations

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
- Compares the selected top functions across metadata groups.
- The point and interval summarize the estimated group-level effect and uncertainty.
- Larger separation and less overlap suggest stronger differences, but exact significance comes from the results table.

#### PCA Plot
- Each point represents one sample based on its complete predicted functional profile.
- Nearby points have similar profiles; separated points have more dissimilar profiles.
- Group separation is exploratory and should be checked against the DAA statistics.

#### Heatmap
- Rows are the top DAA-ranked functions and columns are samples.
- Color intensity represents abundance after the plotting transformation.
- Clustering groups samples or functions with similar abundance patterns.
- A visible pattern is not, by itself, evidence of statistical significance.

#### Statistical Table
- Feature or pathway ID and description
- Group comparisons and estimated effect size
- Raw p-value and BH-adjusted p-value (FDR)
- Interpret FDR together with effect size, replication, and the plots

## Key Concepts

**Differential Abundance Analysis (DAA)**: Testing whether predicted functions differ in abundance between groups.

**Adjusted p-value (FDR)**: A multiple-testing-adjusted measure of evidence. It helps control false discoveries when many functions are tested.

**ggpicrust2**: Applies DAA methods to functional pathways

**KO (KEGG Ortholog)**: Functional gene/protein identifier

**PICRUSt2**: "Prediction of Metabolic Intermediate Genes by Functional Tools"
- Predicts metabolic genes/pathways from 16S/ITS data

**Effect size**: The estimated magnitude and direction of the difference between groups. Statistical significance and effect size answer different questions.

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
