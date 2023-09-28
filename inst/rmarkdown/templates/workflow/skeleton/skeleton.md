---
title: "Analysis Workflow"
author: "Alexandra Cassano"
date: "2023-09-26"
output:
  html_document: 
    df_print: tibble
    highlight: tango
    theme: united
    number_sections: yes
    keep_md: yes
    toc: yes
    toc_float: 
        collapsed: true
  pdf_document:
    number_sections: yes
    keep_tex: yes
---



# Introduction

This code should be usable to analyze RNA-seq data.
Please read through the document before making any changes.
This is a template and will indicate where you need to make changes for your specific experiment.

This R markdown file is broken down into "code chunks" that can be run individually.
Downstream chunks rely on upstream chunks so running code out of order is not recommended.
In R Studio, you can run a chunk by hitting the "play button (▶️)" on the top right of the chunk.
You can run all upstream chunks by clicking the button that has a down arrow pointing to a green bar.

*Built with R version 4.3.1. Most recent template changes made on September 26, 2023*

# Prepare work-space

First you have to set up the "work-space" so that your computer knows what code it will need to import and where all of the raw data is.

## Install Packages

This segment of code should only be run if the packages are not already installed on whichever computer this is being run on.
In order to have this chunk run, change `eval=FALSE` to `eval=TRUE` in the chunk header.
If you're running chunk by chunk be sure to run this chunk the first time but do not run this code subsequent times.

[**Do not edit this chunk (just hit the run chunk button if needed)**]{style="color: red"}


```r
install.packages(c("BiocManager",
                   "tidyverse",
                   "gplots",
                   "ggplot2",
                   "ggrepel",
                   "stringi",
                   "pheatmap",
                   "stats",
                   "RColorBrewer"))

BiocManager::install(c("DESeq2",
                       "GO.db",
                       "GOstats",
                       "pathview",
                       "gage",
                       "uwot",
                       "gageData",
                       "GenomicRanges",
                       "Repitools",
                       "clusterProfiler",
                       "DOSE",
                       "pathview",
                       "treemap",
                       "biomaRt",
                       "ashr"))
BiocManager::install("enrichplot")
```

## Load packages

The first step is to load the R "packages" needed for all of the code below.
These will be used for data manipulation, analysis, and visualization.

[**Do not edit this chunk**]{style="color: red"}


```r
library(devtools)
library(roxygen2)
library(tidyverse)
library(utils) 
library(stringi)
library(stringr)
library(knitr)
library(rmarkdown)
library(biomaRt) 
library(DESeq2) 
library(ashr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(uwot)
library(ggVennDiagram)
library(pheatmap)
library(tinytex)
library(stats)
library(enrichplot)
library(clusterProfiler)
library(DOSE)

if (!require(devtools)) { install.packages('devtools') }
devtools::install_github('JamesJoly/DGSEA')
devtools::install_github('amcassano/RNAseqAnalysisPackage') 
library(DGSEA)
library(analyzeRNA) #arguably the most important line of code, this imports all the code I wrote
library(org.Mm.eg.db)
library(ggridges)
library(ggpubr)
library(treemap)
library(pathview)
library(SPIA)
```

## Set working directory

This chunk will set the 'current working directory'.
You must change this to reflect what folder you'll be working in.
The working directory is all where all outputs will be saved.
All files being imported should also be in this folder.

**Your working directory should contain:**

1.  **CSV file containing all of the raw read counts obtained from the fastq files**
2.  **CSV file containing the meta-data for your experiment**
3.  **CSV files with any gene sets of interest**

[**Change the working directory variable (change what is in the quotes) You will need to use the complete path. The path should end with the name of the folder, you should omit the trailing '/'**]{style="color: red"}


```r
# set working directory
cwd <- "parentfolder/childfolder/yourfolder"
setwd(cwd)
```

## Data import & preparation

### Import CSV files

The `raw_counts` variable is set as the csv file containing the output of feature counts after using STAR to align reads to the reference genome.
The csv file must be "cleaned up" by removing unneeded columns such as *Chr*, *start*, *end* prior to loading into R.

The `meta` variable is set as the csv file containing metadata about the current experiment.
At a minimum, it must contain sample names and conditions/groups.
**The condition column must be named "Group" & the sample ID column must be named "SampleID".** You can also include other information such as sac day, cell type, or anything else you might want to keep track of, let me know if you do though so I can help with some minor code adjustments.

Samples not of interest (in a combined experiment, irrelevant groups) should be removed from both raw counts and metadata.

The samples **must** be in the same order in the `meta` file as they are in `raw_counts`.

Pre-filtering is applied to the raw data to remove any genes that have a total count number across samples less than the user defined minimum.
Additionally, genes that have a total count that is less than a user defined value, when the single max count is subtracted are removed.
This is done to filter out genes that have a high read count for just one sample, but are not expressed in other samples.
These filtering steps minimize noise and remove outliers.

*Avoid using spaces or other white space in the names of your files, use '\_' or CamelCase to separate words instead. Your file names must exactly match what you type below. Be sure to include the `.csv` file extension.*

[**Change`countsCSV` and `metadataCSV` file names**]{style="color: red"}


```r
# import files
countsCSV <- "yourcountsfile.csv"
countsCSV <- paste(cwd, countsCSV, sep = "/")
raw_counts <- read.csv(countsCSV, header = TRUE, sep = ",")

metadataCSV <- "yourmetadatafile.csv"
metadataCSV <- paste(cwd, metadataCSV, sep = "/")
meta <- read.csv(metadataCSV, header = TRUE, sep = ",")
```

### CSV Cleanup

This code will accomplish a few things:

1.  Set the row names as the GeneID in the counts file
2.  Set groups as the data type `factor` : Replace whats in quotes with the names of your groups (add or remove lines as necessary). The contents of the "" must exactly match the group as it is written in your metadata csv file. The order they appear is the order your groups will appear in plots and legends so fill in accordingly
3.  Apply pre-filtering to your data to retain only rows that have more than minimum counts across all samples and have a total count that is less than a user defined value, when the single max count is subtracted are removed. The more samples you have, the higher the number should be, the fewer samples you have the lower the number should be. It makes sense to set the number 300-1000. This is done to filter out genes that have a high read count for just one sample, but are not expressed in other samples. This minimizes noise and removes outliers.

**Things to change in this chunk:** Group names, count minimum, number of samples


```r
# renames each row as the GeneID and removes the Geneid column
rownames(raw_counts) <- raw_counts$Geneid
raw_counts <-  dplyr::select(raw_counts, -Geneid)

# set groups as factors
meta$Group <- factor(meta$Group,
                         levels = c("group1",
                                    "group2",
                                    "group3",
                                    "group4",
                                    "group5",
                                    "group6"))

# pre-filtering
count_minimum <- 1000 
numberofsamples <- 20
raw_counts$row_max <- apply(raw_counts[,1:numberofsamples], 1, max)
raw_counts$row_sum <- apply(raw_counts[,1:numberofsamples], 1, sum)
raw_counts <- raw_counts %>%
  dplyr::filter((raw_counts$row_sum - raw_counts$row_max) > count_minimum) %>%
  dplyr::select(-c(row_max, row_sum))

# set the sample ID as the row name for the metadata 
rownames(meta) <- meta$SampleID
meta <- meta %>% dplyr::select(-c(SampleID))
```

### Rename groups (optional)

If you want to rename any groups from how they appear in your csv files (i.e. 'rej' -\> 'Acute Rejection') use the code below.
If this is needed change `eval=FALSE` to `eval=TRUE` .
Delete or add lines as needed.
Replace `group1`, `group2` etc with the name of the group you want to change.
Replace the text in quotes (`"new name 1"`) with what you want the new group name to be.
Fill in the factoring list so that the new group names replace any old ones, include unchanged group names as well


```r
# rename groups 
meta$Group <- recode_factor(meta$Group, group1 = "new name 1")
meta$Group <- recode_factor(meta$Group, group2 = "New name 2")
meta$Group <- recode_factor(meta$Group, group3 = "new name 3")

# refactor groups
meta$Group <- factor(meta$Group,
                         levels = c("new name 1",
                                    "New name 2",
                                    "new name 3",
                                    "group4",
                                    "group5",
                                    "group6"))
```

# Set colors, shapes, labels for plots

This code will set global variables for plotting such as the color palette, labels, and shapes.
The aesthetics dataframe is created for use in creating consistent plots throughout.
One of the package functions has the different options for the colors and shapes in it.

For colors you can choose from black, blue, orange, green, purple, red, yellow, brown, white, or grey.
Fill the `""` with your choice, add or remove rows as necessary but make sure each row ends in a `,` except the last row.
If you want a solid symbol make the fill and outline the same color.
If you would like an open symbol set the fill as white.

For shapes you can choose from 'circle', 'crossed circle', 'diamond', 'crossed diamond', 'square', 'crossed square', 'triangle up', 'triangle down', 'x', 'plus', 'asterik'.
Fill in the `""` (same concept as outline color)

For each group/condition you have you will specify the shape of symbol as well as the outline and fill colors desired.
Each thing should be in the same order of the groups, the code below will also print out your options and the order of your groups.


```r
# get the names of the groups to assign to the colors/shapes
sample_labels <- (meta$Group)
group_labels <- unique(sample_labels) %>% sort() %>% print()

analyzeRNA::getPalettes()
```

**Things to change in this chunk:** the colors and shapes.
Make sure that the number of arguments provided match the number of conditions you have.
Add or remove lines from the `anno_colors` in the below format for the number of different conditions you have.


```r
plot_aes <- analyzeRNA::set_aes(group_labels, 
                    c("shape1", "shape2", "shape3"),
                    c("outline1", "outline2", "outlin3"),
                    c("fill1", "fill2", "fill3"))

anno_colors <- list(
  Group = c(
    as.character(plot_aes$labels[1]) = as.character(plot_aes$outlineID[1]),
    as.character(plot_aes$labels[2]) = as.character(plot_aes$outlineID[2]),
    as.character(plot_aes$labels[3]) = as.character(plot_aes$outlineID[3]),
    as.character(plot_aes$labels[4]) = as.character(plot_aes$outlineID[4]),
    as.character(plot_aes$labels[5]) = as.character(plot_aes$outlineID[5]),
    as.character(plot_aes$labels[6]) = as.character(plot_aes$outlineID[6])
  )
)             
```

# DESeq2 package

Now that you've got everything set up, its time to start analysis!
The DESeq2 package has functions to find differentially expressed genes, statistics associated with this, and normalized counts.

First a DESeq data object is created using `raw_counts` and `meta`, specifying that `group` is the variable of interest in the metadata table.
The `DESeq` function is run and from this output, normalized results (rlog transformed, and normalized counts).
Differential gene expression requires pairwise comparisons, however analysis using the normalized counts (UMAP, PCA, Sample Distances) use all samples at once.

## Create DESeq object and run DESeq

This code will construct the `DESeqDataSet` which takes the raw counts and metadata as inputs.
`design = ~Group` means that Group is your experimental variable, this variable expresses how the counts for each gene depend on the variables in metadata.
If the experimental design is more complex than just group (i.e. group and cell type) then you can combine variables (i.e. `~ Group + CellType`) (but lmk if this is the case as we might need to make some slight edits to the code).

DESeq uses a generalized linear model and outputs a DEseq2 object, we will use this to get normalized counts, p-values, and log-2 fold changes.

**Only change this chunk if you need to change the experiment design**


```r
dseq <- analyzeRNA::getDESeq(raw_counts, meta, ~Group)
```

## Regularized Log Transformation

This code will apply a regularized log transformation to the data.
This transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size.
- - The `rlogtransformation` produces a similar variance stabilizing effect as `varianceStabilizingTransformation`, though rlog is more robust in the case when the size factors vary widely.
The transformation is useful when checking for outliers or as input for machine learning techniques such as clustering or linear discriminant analysis.
rlog takes as input a `DESeqDataSet` and returns a `RangedSummarizedExperiment` object.
The rlog transformation takes more time to calculate than the VST function but it is more robust.

`blind = FALSE` must be included in the code.
This determines whether to blind the transformation to the experimental design.
`blind=TRUE` should be used for comparing samples in an manner unbiased by prior information on samples (for QC).
`blind=FALSE` should be used for transforming data for downstream analysis, where the full use of the design information should be made.
If many of genes have large differences in counts due to the experimental design, it is important to set `blind=FALSE` for downstream analysis.

**Nothing needs to be changed in this chunk**


```r
rlog_norm <- DESeq2::rlog(dseq, blind = FALSE)
rlog_df <- as.data.frame(SummarizedExperiment::assay(rlog_norm))
```

## Pairwise comparisons

Multiple pairwise comparisons can be made to determine the differentially expressed genes for any 2 given groups.
Pairwise comparisons between 2 specified groups are saved as objects for use further on.

You can choose to use a Log Fold Change shrinking algorithm.
This mitigates the effects of any single sample on the log fold change calculations, minimizing noise.
Add `log2FCshrink = TRUE` to the function call if you want to use this.
The calculation used is the "ashr" algorithm [@stephens2016].

**To change:** Add the desired comparisons, change the names of the numerator and denominator and the names of the objects you're storing the results in.


```r
resultsNames(dseq)

one_vs_two <- analyzeRNA::pairwiseDEGresults(numerator = "comparisongroup",
                                denominator = "baselinegroup",
                                deseq_obj = dseq)
                                
two_vs_three <- analyzeRNA::pairwiseDEGresults(numerator = "othergroup",
                                denominator = "anotherone",
                                deseq_obj = dseq)
```

# Annotate Results

The initial `raw_counts` table, and therefore all DESeq results derived from it, use the *Ensembl ID* to identify each gene.
However, this is not a useful identifier for most humans.
The BiomaRt package can be used to match the *Ensembl ID* with it's corresponding *MGI Symbol*, *MGI Description*, (and optionally things like the *Gene Biotype*, and *Entrez (NCBI) ID*).
The genemap can be altered to omit or include different attributes as needed.
The genemap is then used along with a `deseq results object` in the `annotate_biomart()` function to add these descriptors to the results.

Any data frame that has *Ensembl ID*s in a column called `GeneID` can be used in this annotation function.

## Make genemap

**Do not change** (unless you want to add more things like gene biotype)


```r
allgenes <- tibble::rownames_to_column(rlog_df, var = "GeneID")
allgenes <- select(allgenes, GeneID)
gmap <- analyzeRNA::createGenemap(allgenes,  GMattributes = c("ensembl_gene_id", "mgi_symbol", "mgi_description"))
```

## Annotate pairwise comparisons (DEGs)

**To change:** the names of your pairwise comparisons, add so that all of your comparisons get annotated


```r
one_vs_two <- analyzeRNA::annotate_biomart(one_vs_two, gmap) 
two_vs_three <- analyzeRNA::annotate_biomart(two_vs_three, gmap) 

rlog_df <- tibble::rownames_to_column(rlog_df, var = "GeneID")
rlog_df <- analyzeRNA::annotate_biomart(rlog_df, gmap)
```

# Filter & sort significant DEGs

The list of genes returned by the `results` function of the DESeq2 package contains all the genes above a 0.1 P-value.
This is a rather generous filtering, especially with so many genes and so many pairwise comparisons being done.
The results also include some genes that have uncalculated P-values or Log2 Fold Changes due to the presence of a sample that has been detected as "an extreme outlier".
There are also results with uncalculated adjusted P-values if there is a very low mean normalized count.

The `signif_deg` function will filter out these genes and only keep DEGs that meet a user specified adjusted P-value threshold and log 2 fold change threshold.
In addition, it will add a `Change Direction` column to indicate if a gene is up-regulated or down-regulated in the experimental group compared to the baseline group (as specified when retrieving results above).
`Change Metric` is also added to each gene.
This will calculate $Change\space Metric = \log_2(Fold \space Change) *-\log(adj.\space P \space value)$.

**To change:** threshold values (if desired), and adjust the results data frames being generated


```r
# set thresholds for adjusted P values and log2 fold changes in determining which genes to keep as significant
padj_threshold <- 0.01
l2fc_threshold <- 1

# create ranked and filtered list of dif expr genes
ranked_OneTwo <- analyzeRNA::signif_deg(one_vs_two, padj_threshold, l2fc_threshold)
ranked_TwoThree <- analyzeRNA::signif_deg(two_vs_three, padj_threshold, l2fc_threshold)
```

## Export lists of significantly DEGs

The `export_genelist` function will save the output of the `signif_deg` function as a csv, exporting it to the current working directory.
The option to export just those genes that are up-regulated, just those that are down-regulated, or all genes is included

**To Change:** the exports and the file names, feel free to export just up/down/both and not all three


```r
analyzeRNA::export_genelist(ranked_OneTwo, filename = "1vs2", direction = "up")
analyzeRNA::export_genelist(ranked_OneTwo, filename = "1vs2", direction = "down")
analyzeRNA::export_genelist(ranked_OneTwo, filename = "1vs2")

analyzeRNA::export_genelist(ranked_TwoThree, filename = "2vs3", direction = "up")
analyzeRNA::export_genelist(ranked_TwoThree, filename = "2vs3", direction = "down")
analyzeRNA::export_genelist(ranked_TwoThree, filename = "2vs3")
```

# Data visualizations

For the sheer amount of data in an RNAseq experiment, one of the most useful tools is data visualizations.
This section includes dimensionality reduction visualization techniques as well visualizations of pairwise comparisons and visualizations of expression of specific genes of interest.
Visualizations that use all of the samples across conditions will use the rlog normalized data.
Visualizations showing DEGs between two conditions will use the annotated results of pairwise comparisons.

## Sample to Sample Distances

Plotting sample distances will illustrate the overall similarity between samples.
This can be compared to what is expected by the experiment design.
Rlog transformed data will be used for this to ensure that sample distance isn't overly due to certain highly expressed genes.

**nothing needs to be changed here**


```r
analyzeRNA::sample_distances(rlog_norm, anno_colors, meta)
```

## Principal Component Analysis

Principal component analysis will illustrate the differences between samples as a whole, reducing the number of dimensions.
PCA is most useful for a smaller number of samples (i.e. what we're working with in bulk RNA seq experiments, but not in single cell sequencing experiments).

PC1 will be the x axis and differences in PC1 are responsible for the most variance in the data set.
The percent of variation accounted for in PC1 will be specified in the axis label.
PC2 is the y axis and PC2 is responsible for the second most variance.
The percent of variation here will also be included in the axis label.

The major limitation of PC analysis is that you cannot say for certain what factor or combinations of factors are contributing to the variance.
However, you do get very clear insight into overall similarity and differences between samples.
Furthermore, inclusion of a group of samples that are vastly different than the rest will limit the information you can glean from a PCA plot.
If the difference between naive samples and any treatment condition is more significant than any other differences, then you will not be able to see other differences present.
To further interrogate your samples I recommend removing samples like naive from your data set and running the PCA analysis on that.
This will let you see what groupings emerge without something like naive driving \~80% of the variance.

For a PCA plot to really give good information, you would hope to see the principal components account for \~50-60% of the total variance.
Ideally, PC1 would account for \~30+% of the variance.

**Change the title and if you want to label the dots change it to TRUE**


```r
# run the function
analyzeRNA::pca_analysis(dseq_transform = rlog_norm,
             plot_aes = plot_aes,
             plot_title = "PCA",
             label_samples = FALSE)
```

## UMAP

UMAP is another dimensionality reduction technique.
UMAP is better suited for a much higher number of samples (single cell vs bulk).
Odds are you won't really use UMAP for anything here but I already wrote the code so I'm not getting rid of it.
UMAP, unlike PCA does not provide you with the numerical amount of variance accounted for.

*Warning:* Every time the UMAP function is (re)run, the resulting plot will be different.
This is due to "randomness" inherent in the UMAP algorithm.
Either take steps to force reproducible results, or be cautious when re-running UMAP function.
The result of the UMAP function is saved as an object so the plot can be displayed without re-calculating UMAP.
It does not incorporate "true randomness" since computers are not capable of this.
Functions achieve "pseudorandomness" by setting a number called a seed which changes every time the function is run.
This number is generated using the system clock of the computer.
The results can be forced to be reproducible by adding `batch = TRUE` to the umap() call, and by adding `set.seed(x)` where `x` is an integer you choose.

**you can change the seed around, play with the neighbors and m Dist, change the title and toggle labels**


```r
set.seed(5)
rlogUMAP <- umap_analysis(dseq_transform = rlog_norm,
                          neighbors = 3,
                          m_dist = 0.5,
                          cond_list = sample_labels,
                          plot_aes = plot_aes,
                          plot_title = "rlog Normalized UMAP",
                          batch = TRUE,
                          label_samples = FALSE)

rlogUMAP
```

## Plot counts

-   The counts for specific genes can be plotted by group using the normalized counts.
-   This can be done for specific genes of interest.
-   Using the `plot_genecounts` function, the data can be saved to an object for plotting with ggplot
-   For this function, you will provide the MGI symbol for the gene of interest.
    -   The function also takes as input the annotated data frame of rlog annotated counts, the metadata data frame, and the plot aesthetics.
    -   Optional arguments include
        -   a list of comparisons to get p values for. The format for this is \`comparisons = list(c("condition1", "condition2"), c("condition2", "condition3")) add as many comparisons as desired. The order of comparisons in the list will correspond to the order in which the pvalues are plotted on the graph.
        -   changing the y axis label using `yaxis = "my custom axis label"`. it defaults to "rLog Normalized Reads"
        -   labeling the samples using `label_samples = TRUE`, it defaults to false
        -   you can also set it so the Y axis is forced to go to 0. By default the Y axis will zoom in around the counts. You can change this by including `yaxistozero = TRUE`
        -   you can change the location of the kruskal wallis statistic that will be included using `kw_hjust = yournumberhere, kw_vjust = yournumberhere`
        -   by default the x axis will split your samples according to the Group variable, if there is a different name for it (i.e. condition) you can include `grouping_var = "Condition"`
-   **change the mgi symbol names & add as you wish. edit the optional variables as you wish and add comparisons as desired, the below are just examples**


```r
statcomps <- list(c("Naive", "Tolerant"), c("Naive", "Rejection"), c("Tolerant", "Rejection"))

plot_genecounts("Pdcd1", rlog_df, meta, plot_aes)
plot_genecounts("Tox", rlog_df, meta, plot_aes, comparisons = statcomps, yaxistozero = TRUE)
```

## Volcano Plot

-   A volcano plot will show all of the genes in the pairwise comparison results

    -   the X axis on the volcano plot is the Log2 Fold Change for each gene
    -   the Y axis is the -Log10 Adjusted P value for each gene
    -   a horizontal and vertical line will be added to the plot to indicate the threshold for significance for fold change and p value respectively

-   The adjusted P value threshold defaults to 0.01.
    to change that, add `pval_cutoff = 0.005` to the function call (use any number you desire, thats just an example)

-   the log2 fold change threshold defaults to 2.
    to change that add `l2fc_cutoff = 1` to the function call (again use any number you want)

-   the plot will label significant genes with their MGI symbol (not all genes may be labeled if there are too many to plot without crowding the plot)

-   the function takes as input:

    -   results data frame for pairwise comparison, this should be annotated already, use the full data frame, not the one that has been filtered for significant genes already
    -   plot title string
    -   2 strings indicating the names of the conditions in the pairwise comparison. they should be in order so that the numerator for the DEG results is listed before the denominator
    -   *optional* change the threshold values as specified above
    -   *optional*

    

## Gene Set Enrichment Analysis

The dotplot shows the number of genes associated with the first 50 terms (size) and the p-adjusted values for these terms (color).
This plot displays the top 50 genes by gene ratio (\# genes related to GO term / total number of sig genes), not p-adjusted value.
-



## Heatmaps

This code chunk will create a heatmap of the most significantly changed genes across all samples and their expression in each sample.
This visualization can easily show the effects of treatments on different genes.
It also will cluster genes.
Gene sets gotten from GSEA



### Shared Genes



### VennDiagram


