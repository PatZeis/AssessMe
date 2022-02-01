# AssessME
A cluster assessment tool for preprocessing and clustering optimisation

# Description
A tool for assessment and comparison of cluster partitions based on different:filtering, feature selection, normalization, batch correction, imputation, clustering algorithms

# Vignette

## Install required packages
``` r
install.packages("devtools")
library(devtools)
install_github("PatZeis/AssessMe")
```

``` r
library(AssessME)
```

#### run Seurat with increasing resolution on small intestinal epithelial cells from Haber et al.
``` r
x <- read.delim("~/Downloads/GSE92332_atlas_UMIcounts.txt")
entero <- CreateSeuratObject(counts = x, project = "10x", min.cells = 3, min.features = 200)
entero <- NormalizeData(entero, normalization.method = "RC", scale.factor = 10000)
entero <- FindVariableFeatures(entero, selection.method = "vst", nfeatures = 3000)
features <- Seurat::VariableFeatures(entero)
entero <- ScaleData(entero, features = features)
entero <- RunPCA(entero, features = features, npcs = 100)
entero <- FindNeighbors(entero, dims = 1:100)
resolution <- c(1:10)
for (i in resolution)  {
  entero  <- FindClusters(entero , resolution = i)
}
```



#### run AssessME for the different resolution based cluster partitions
``` r
res <- colnames(entero@meta.data)[c(4,6:length(colnames(entero@meta.data)))]

for(i in 1:length(res)) {
  if (i == 1) {
    assess_seurat_entero <- cluster_assessment( seuratobject=entero,givepart = res[i], give2ndfiff=F, Entro_med=F, diptest=F, run_enriched=T, bwidth=F, critmass=F, gooutlier=T)
  }
  else {
    assess_seurat_entero <- cluster_assessment(assessment_list = assess_seurat_entero, seuratobject=entero,givepart = res[i], give2ndfiff=F, Entro_med=F, diptest=F, run_enriched=T, bwidth=F, critmass=F, gooutlier=T)

  }
}
```

#### plot results of AssessME and identify optimal resolution for cell type identification
``` r
opti_resolution_plot(assess_seurat_entero, f1_thr = 0.5, max_leng = 3, lcol = "red", resolu = T)
```

[Identifying resolution parameter sweet spot. -click to see image-](images/opti_resolution.png)
``` r
goblet <- c("Tff3", "Manf", "Ccl9")
middle <- c("Sox4", "Dll1")
A <- c("Serpina1c", "Ghrl")
SILA <- c("Cck", "Gal")
SILP <- c("Pyy", "Gcg")
SAKD <- c("Sst", "Iapp")
ECReg4 <- c("Reg4", "Afp")
marker_entero <- c(goblet, middle, A, SILA, SILP  , SAKD, ECReg4)
binaryclass_scatterplot(assess_seurat_entero$Sres.2, assess_seurat_entero$Sres.6, out_doub_dif = marker_entero, maplot = T, logmean = "log2")
```

[Identifying unresolved cell types by increasing resolution. -click to see image-](images/entero_scatter_F1.png)

#### preprocessing of Tusi et al., count data of hemaptopoietic progenitors
``` r
filt_norm_counts <- read.csv("GSM2388072_basal_bone_marrow.filtered_normalized_counts.csv", row.names = 1)
filt_norm_counts_matrix <- filt_norm_counts[,-(1:4)]
normalized_data_tusi <- filt_norm_counts_matrix[filt_norm_counts$seq_run_id == "seq_run4",]
raw_counts <- read.csv("GSM2388072_basal_bone_marrow.raw_umifm_counts.csv", row.names = 1)
raw_counts_matrix <- raw_counts[,-(1:4)]
rawdata <- t(raw_counts[rownames(normalized_data_tusi), colnames(normalized_data_tusi)])
```


#### run Seurat with library size normalisations and increasing resolution
``` r
require(Seurat)
tusi_sc <- CreateSeuratObject(counts = rawdata, project = "celseq", min.cells = 1, min.features = 100)
tusi_sc <- NormalizeData(tusi_sc, normalization.method = "RC", scale.factor = 10000)
tusi_sc <- FindVariableFeatures(tusi_sc, selection.method = "vst", nfeatures = 3000)
features <- Seurat::VariableFeatures(tusi_sc)
tusi_sc <- ScaleData(tusi_sc, features = features)
tusi_sc <- RunPCA(tusi_sc, features = features, npcs = 100)
tusi_sc <- FindNeighbors(tusi_sc, dims = 1:100)
resolution <- c(1:10)
for (i in resolution)  {
  tusi_sc <- FindClusters(tusi_sc, resolution = i)
}
```

#### run AssessME for the different cluster partitions of the seurat object
``` r
res <- colnames(tusi_sc@meta.data)[c(4,6:length(colnames(tusi_sc@meta.data)))]
for(i in 1:length(res)) {
  if (i == 1) {
    assess_seuratRC <- cluster_assessment( seuratobject=tusi_sc,givepart = res[i], give2ndfiff=F, Entro_med=F, diptest=F, run_enriched=T, bwidth=F, critmass=F, gooutlier=T)
  }
  else {
    assess_seuratRC <- cluster_assessment(assessment_list = assess_seuratRC,  seuratobject=tusi_sc,givepart = res[i], give2ndfiff=F, Entro_med=F, diptest=F, run_enriched=T, bwidth=F, critmass=F, gooutlier=T)

  }
}
```

#### run Seurat with log - library size normalisations and increasing resolution
``` r
tusi_sc_log <- CreateSeuratObject(counts = rawdata, project = "celseq", min.cells = 1, min.features = 100)
tusi_sc_log <- NormalizeData(tusi_sc_log , normalization.method = "LogNormalize", scale.factor = 10000)
tusi_sc_log <- ScaleData(tusi_sc_log , features = features)
tusi_sc_log <- RunPCA(tusi_sc_log , features = features, npcs = 100)
tusi_sc_log <- FindNeighbors(tusi_sc_log, dims = 1:100)
for (i in resolution)  {
  tusi_sc_log  <- FindClusters(tusi_sc_log, resolution = i)
}
tusi_sc_log[["RNA"]]@var.features <- features
```

#### run AssessME for the different cluster partions using log-library size normalisation
``` r
res <- colnames(tusi_sc@meta.data)[c(4,6:length(colnames(tusi_sc@meta.data)))]
for(i in 1:length(res)) {
  if (i == 1) {
    assess_seuratLog <- cluster_assessment( seuratobject=tusi_sc_log,givepart = res[i], give2ndfiff=F, Entro_med=F, diptest=F, run_enriched=T, bwidth=F, critmass=F, gooutlier=T)
  }
  else {
    assess_seuratLog <- cluster_assessment(assessment_list = assess_seuratLog, seuratobject=tusi_sc_log,givepart = res[i], give2ndfiff=F, Entro_med=F, diptest=F, run_enriched=T, bwidth=F, critmass=F, gooutlier=T)
  }
}
```

#### run Seurat with SCTransform and increasing resolution
``` r
tusi_sc_SCT <- CreateSeuratObject(counts = rawdata, project = "celseq", min.cells = 1, min.features = 100)
tusi_sc_SCT <- SCTransform(tusi_sc_SCT, verbose = FALSE)
tusi_sc_SCT <- RunPCA(tusi_sc_SCT, features = features, npcs = 100)
tusi_sc_SCT <- FindNeighbors(tusi_sc_SCT, dims = 1:100)
resolution <- c(1:10)
  for (i in resolution)  {
    tusi_sc_SCT <- FindClusters(tusi_sc_SCT, resolution = i)
  }
tusi_sc_SCT[["RNA"]]@var.features <- features
```

#### run AssessME for the different cluster partions using SCTransform normalisation
``` r
res <- colnames(tusi_sc@meta.data)[c(6,8:length(colnames(tusi_sc@meta.data)))]
  for(i in 1:length(res)) {
    if (i == 1) {
      assess_seuratSCT <- cluster_assessment( seuratobject=tusi_sc,givepart = res[i], give2ndfiff=F, Entro_med=F, diptest=F, run_enriched=T, bwidth=F, critmass=F, gooutlier=T) }
    else {
      assess_seuratSCT <- cluster_assessment(assessment_list = assess_seuratSCT, seuratobject=tusi_sc,givepart = res[i], give2ndfiff=F, Entro_med=F, diptest=F, run_enriched=T, bwidth=F, critmass=F, gooutlier=T)
    }
}
```

#### combine list of assessments and Plot results
``` r
assessment_list_sc <- list(seuratLog=assess_seuratLog, seuratRC=assess_seuratRC, seuratSCT=assess_seuratSCT)
preprocess <- opti_preprocess_plot(assessment_list_sc, lcol = c("red", "blue", "green"))
```

[Assessment of different normalisation methods. -click to see image-](images/opti_preprocess.png)


#### in-depth analysis RC and log-RC normalisation with the same resolution/No. clusters partitions
``` r
relevant hematopoietic cell type marker
marker_erythroid <- c("Hbb.bt","Hba.a2","Hba.a1","Alas2","Bpgm")
marker_basophil_mastcell <- c("Lmo4","Ifitm1","Ly6e","Srgn")
marker_megakaryo <- c("Pf4","Itga2b","Vwf","Pbx1","Mef2c")
marker_mpp <- c("Hlf","Gcnt2")
marker_lymphocytes <- c("Cd79a","Igll1","Vpreb3","Vpreb1","Lef1")
marker_dendritic <- c("H2.Aa","Cd74","H2.Eb1","H2.Ab1","Cst3")
marker_monocytes <- c("Csf1r","Ly6c2","Ccr2")
marker_granu_neutro <- c("Lcn2","S100a8","Ltf","Lyz2","S100a9")
marker <- c(marker_monocytes, marker_granu_neutro, marker_dendritic, marker_lymphocytes, marker_mpp, marker_megakaryo, marker_basophil_mastcell, marker_erythroid)
```

``` r
binaryclass_scatterplot(assessment_list_sc$seuratRC$Sres.2, assessment_list_sc$seuratLog$Sres.2, out_doub_dif = marker, maplot = T, logmean = "log2")
```
[Library size normalisation leads to better resolution of hematopoietic progenitors. -click to see image-](images/tusi_scatter_F1.png)

``` r
binaryclass_scatterplot(assessment_list_sc$seuratRC$Sres.2, assessment_list_sc$seuratLog$Sres.2, out_doub_dif = marker, maplot = T, logmean = "log2", toplot = "Entropy")
```

[Library size normalisation leads to better resolution of hematopoietic progenitors. -click to see image-](images/tusi_scatter_entropy.png)



#### Test robustness of clusters
``` r
accuracy_list <- accuracy(giveassessment = assess_seurat_entero$Sres.1, data = entero@assays$RNA@counts, ntree = 100, crossvali = 50)
accuracy_list <- accuracy(accuracy_list = accuracy_list,giveassessment = assess_seurat_entero$Sres.6, data = entero@assays$RNA@counts, ntree = 100, crossvali = 50)
```

#### Plot and compare robustness of different cluster partitions
``` r
accuracy_plot(accuracy_list)
```
[Comparison of accuracy of re-classification of cells of the different cluster partitions tested. -click to see image-](images/compare_accuracy.png)

### Calculate robustness of different cluster partitions and plot on a HPC system
#### "assess_seurat_entero.Rda" = saved list of assessments, "entero.RDa" saved Seurat object
``` bash
zeis@example:~$ mkdir res1
zeis@example:~$ cd res1
zeis@example:~/res1$ for i in {1..50}; do sbatch ~/AssessME_accuracy_hpc_bash_script.sh -a ~/assess_seurat_entero.Rda -j ~/entero.RDa -t $i -i 1; done
zeis@example:~/res1$ cd ..
zeis@example:~$ cd res6
zeis@example:~/res6$ for i in {1..50}; do sbatch ~/AssessME_accuracy_hpc_bash_script.sh -a ~/assess_seurat_entero.Rda -j ~/entero.RDa -t $i -i 6; done
zeis@example:~/res6$ cd ..
zeis@example:~$ Rscript --vanilla ~/AssessME_accuracy_plot_hpc_script.R ~/res1/ ~/res6
```
