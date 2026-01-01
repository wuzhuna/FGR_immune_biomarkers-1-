# FGR_immune_biomarkers

##Study Overview
This code is used to reproduce the analytical workflow described in the paper Screening of Key Immune Biomarkers Reveals GAL and F2R as Promising Candidates for Diagnosing Fetal Growth Restriction.

##Runtime Environment
R Version: 4.1.3
Required R Packages: limma, ggplot2, clusterProfiler, glmnet, pROC, e1071, sva, rms, etc.
Detailed package versions are provided in sessionInfo.txt.

##File Descriptions
DEG_analysis.R: Screening and visualization of differentially expressed genes
venn_analysis.R: intersecting IR-DEGs between IRGs and DEGs
GO_analysis.R: GO analysis
KEGG_analysis.R: KEGG analysis
PPI_analysis: TRING database(https://string-db.org/)
hubGene_analysis:Cytoscape 3.10.0
heatmap_analysis.R:expression patterns of the top 10 hub IR-DEGs in the FGR
vol_analysis.R:the distribution of the top 10 hub IR-DEGs
corrr_analysis.R:correlations between 10 hub genes.
scatter_analysis.R:some highly correlated IR-DEGs is provided
lasso_analysis.R:LASSO regression
SVM_analysis.R:SVM-RFE algorithms
venn_analysis.R:The intersection of LASSO and SVM-RFE
circlize_analysis.R:Chromosomal positions of the key IR-DEGs are presented
PCA_analysis.R:illustrates the distribution of samples based on the expression profiles of the 3 key IR-DEGs
boxplot_analysis.R:Comparative expression levels of three crucial IR-DEGs
ROC_analysis.R:validate the efficacy of three crucial IR-DEGs in predicting FGR
Nomo_analysis.R:predicting the risk of FGR is constructed based on three IR-DEGs
CIBERSORT_analysis.R:the relative abundances of 22 distinct immune cell types
barplot_analysis.R:the differentiation in ratios of 22 immune cell types between FGR and AGA
cor_analysis.R:illustrate the correlative relationships between two IR-DEGs and infiltrating immune cells

##Running Procedures
Install the required R packages.
Run each script in sequence.
Output results include figures and tables.

##Citation
Please cite this paper if you use this code.

##Contact Information
For any questions, please contact: [wuzhuna@aliyun.com]
