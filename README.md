# CellFindPy

## **Description:**
CellFindPy performs recursive community detection with automatic gamma/resolution parameter optimization based on biologically meaningful parameters such as gene log2 fold expression difference and false discovery rate. Initial community partitioning produces parent communities which are recursively partitioned. A child community is partitioned from the parent community if it passes a defined genetic test. However, the child community is only retained if it is robust to bootstrapping (resampling with replacement). If a child community is retained, it is recursively partitioned and so on. After all communities have been partitioned, a final genetic test ensures that all communities are unique. With synthetic sc-RNA-seq data composed of known ground truth communities, CellFindPy demonstrates superior performance at detecting cell populations than Leiden or Louvain algorithms (mean adjusted Rand Index improvement = ~15%). With a real sc-RNA-seq dataset of over 1.2 million cells across 23 replicates, CellFindPy partitions up to 7 times as many small populations (population size <1% of total sample) while maintaining excellent replicability that is comparable to leiden and louvain replicability. Additionally, CellFindPy enables biologists to partition cell populations using biologically meaningful parameters in place of abstract parameters.

## **How to use:**
1) Download the package into a directory of your choice.
2) From the command line, cd into the directory containing your Anndata .h5ad file.
3) Type `python /path/to/CellFindPy/CellFindPy.py -i file.h5ad -o output_folder`

### **Optional Parameters:**
Optional parameters can be used by advanced users or if the default parameters are not optimal for your dataset.
\
For more information, type `python /path/to/CellFindPy/CellFindPy.py -h`
\
-i -------- Provide Anndata file name.\
-o -------- Provide output folder.\
-s -------- True or False: Do you want to produce gene expression excel sheets?\
-c -------- leiden or louvain community detection?\
-b -------- True or False: Do you want scanpy.rank_genes_groups to adjust for batch?\
-bk ------- If you want scanpy.rank_genes_groups to adjust for batch, what is the batch key in Anndata.obs?\
-mcnts ---- Genes with less counts than this threshold will be filtered out.\
-mcells --- Genes found in less than this number of cells will be filtered out.\
-nboot ---- For community stability (resampling with replacmenet), how many bootstraps?\
-nbf ------ For community stability (resampling with replacmenet), Minimum fraction of cells consistently clustering together?\
-lg2FC ---- log2(value) threshold for average expression difference of all significant genes.\
-p -------- False Discovery Rate.\
-msg ------ Minimum number of significant genes for each community.\
-ncom ----- Initial number of communities to opitmize resolution to?\
-cfrac ---- Maximum size of largest community during resolution optimization?\
-ng ------- Number of genes from each community to make UMAP projections for?

## **Output:**
CellfindPy aims to identify all possible communities within a dataset based on user defined biological parameters. For larger datasets, this can result in a very large number of communities. To facilitate ease of interpretation and since the CellFindPy communities are hierarchically related, CellFindPy community labels are produced at 3 levels (branch points) of increasing granularity/depth in addition to a results file containing all communities. For interpretation purposes, a nested folder structure is constructed containing excel files containing useful genetic information (e.g. top expressed genes and top differentially expressed genes) and comparisons for each community against other communities at the same level. Additionally, to facilitate ease of comparison with other covariates within the Anndata structure, UMAP projections of all Anndata.obs covariates including community assignments are generated and stored in `/figures`. Within `/figures`, UMAP projections of the top n differentially expressed genes are generated for community (default 20 genes). 

The genetic excel sheet contains the following:
1) Top 100 expressed genes for each community.
2) Mean expression of each gene within each community (non-significant values are zeroed).
3) Top 100 differentially expressed genes for each community (some communities may have less than 100 genes).
4) Differential expression of each gene within each community versus all others at the same level (non-significant values are zeroed).
5) Differential expression of each gene within each community versus a specific community at the same level (non-significant values are zeroed).

## **Requirements:**
anndata==0.6.22.post1\
leidenalg==0.7.0\
numpy==1.19.1\
pandas==1.1.2\
scanpy==1.5.1
