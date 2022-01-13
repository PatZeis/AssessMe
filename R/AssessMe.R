#' @useDynLib AssessME, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RaceID
#' @import diptest
#' @import multimode
#' @import Rcpp
#' @import Matrix
#' @import RcppProgress
#' @import Matrix.utils
#' @import future.apply
#' @import future
#' @import ggpubr
#' @import Seurat


#' @title AssessME - a cluster assessment tool for preprocessing and clustering optimisation
#' @description tool for assessment and comparison of cluster partitions based on different:filtering, feature selection, normalization, batch correction, imputation, clustering algorithms
#' @param assessment_list list, with named slots for different assessments, to which new assessment is added. Default is \code{NULL}.
#' @param seuratobject Seurat object as input for assessment:  derives UMI count object, normalized count object, cluster partition and variable features from Seurat Object. Default = \code{NULL}.
#' @param seurat_assay if \code{seuratobject}, name of Seurat assay to retrieve slots. Default =”RNA”
#' @param seurat_lib_size logical. If \code{FALSE} performs library size normalization of UMI counts slot of \code{seuratobject} and overwrites normalized data slot within assessment object. Default = \code{FALSE}.
#' @param do.features logical. If \code{TRUE} performs feature selection and derives \code{var_feat_len} number of top variable genes. Default = \code{TRUE}.
#' @param var_feat_len number of top variable genes used for cluster assessment, if \code{var_feat_len} not equivalent of the length of "var.features" slot of \code{seuratobject}, derive top \code{var_feat_len} number of feature genes using Seurat’s variance stabilization method, requires \code{seuratobject} and \code{do.features} needs to be set \code{TRUE}. Default = \code{NULL}.
#' @param RaceIDobject RaceID object as input for assessment: derives UMI count data of cells passing filtering criteria, normalized data, cluster partition, feature genes, background noise model describing the expression variance of genes as a function of their mean and RaceID filtering criteria. Default = \code{NULL}.
#' @param RaceID_cl_table metadata data frame for a RaceID object in similar form as meta.data slot of a Seurat object with rows as cells and columns as e.g. different cluster partitions. Default = \code{NULL}.
#' @param ScanpyobjectFullpath full path to scanpy object in h5ad format, which is converted to Seurat object from which UMI counts, cluster partition and feature genes are derived. Using UMI count data and scale factor, library size normalization is performed and scaled using the scale factor.
#' @param scanpy_clust either “leiden” or “louvain”, derives cluster partition of either Leiden or Louvain clustering. Default=”leiden”.
#' @param scanpyscalefactor integer number with which relative cell counts are scaled to equal transcript counts. Default = 10,000.
#' @param rawdata UMI count expression data with genes as rows and cells as columns. Default = \code{NULL}.
#' @param ndata normalized expression data with genes as rows and cells as columns. Default = \code{NULL}.
#' @param norm performs library size normalization on provided rawdata argument. Default = \code{TRUE}.
#' @param givepart clustering partition. Either a vector of integer cluster number for each cell in the same order as UMI count table or normalized count table for RaceIDobject; or a character string representing a column name of Seurat metadata data frame of a Seurat object or similar metadata frame, \code{RaceID_cl_table},for a RaceID object. Default = \code{NULL}.
#' @param givefeatures gene vector to perform assessment. Default = \code{NULL}.
#' @param minexpr minimum required transcript count of a gene across evaluated cells. Genes not passing criteria are filtered out. Default 5. If \code{RaceIDobject}, \code{minexpr} derived from  \code{RaceIDobject}. Relevant for deriving feature genes if \code{gene.domain} and calculating fit of dependency of mean on variance.
#' @param CGenes gene vector for genes to exclude from feature selection. Only relevant if \code{seuratobject} \code{&} \code{RaceIDobject} \code{&} \code{ScanpyobjectFullpath} = \code{NULL} and \code{rawdata} is given. Default = \code{NULL}.
#' @param ccor integer value of correlation coefficient used as threshold for determining genes correlated to genes in \code{CGenes}. Only genes correlating less than \code{ccor} to all genes in \code{CGenes} are retained for analysis. Default = 0.65.
#' @param fselectRace logical. If \code{True}, performs RaceID feature selection, only  if \code{seuratobject} \code{&} \code{RaceIDobject} \code{&} \code{ScanpyobjectFullpath} \code{&} \code{givefeatures} = \code{NULL}. Default = \code{False}.
#' @param fselectSeurat logical. If \code{True},performs Seurat variance stabilization feature selection and derives \code{var_feat_len} number of top variable genes, only if \code{seuratobject} \code{&} \code{RaceIDobject} \code{&} \code{ScanpyobjectFullpath} \code{&} \code{givefeatures} = \code{NULL}. Default = \code{False}.
#' @param givebatch vector indicating batch information for cells; must have the same length and order as cluster partition. Default = \code{NULL}.
#' @param individualbatch individual batch name, element of \code{givebatch}, to perform assessment on. Default = \code{NULL}.
#' @param gene.domain logical. If \code{TRUE}, assess all genes with at least \code{minexpr} in one cell.
#' @param PCA_QA logical. If \code{TRUE}, derives first two principal components and the top \code{PCAnum} number of genes with highest or lowest loadings. Default = \code{False}.
#' @param PCAnum integer value, number of genes to be derived with top highest and top lowest loadings for the first two principal components. Default = 10.
#' @param run_cutoff logical. If \code{TRUE} calculate per gene cutoff, representing true label utilized for F1 score, entropy and enrichment of gene per cluster calculation. Default = \code{T}.
#' @param logical. If \code{TRUE} than cutoff for true label is x>0. Default = \code{False}.
#' @param cutoff either “mean” or “median”, utilizes either per gene average expression within clusters or per gene median expression within clusters to calculate the true label cutoff. The Cutoff is calculated per gene by selecting the cluster with highest average or median expression and averaging this mean, with the mean or median of the remaining clusters.
#' @param cutoffmax logical. If \code{TRUE}, then per gene cutoff is the average expression of the cluster with highest average expression across clusters. Default = \code{False}.
#' @param clustsize integer value, threshold of minimum number of cells a cluster should have to be included in the assessment.
#' @param binaclassi either “F1Score”, “Cohenkappa”, “MCC” or NULL. Statistical analysis for binary classification. F1Score, Cohenkappa or Matthews correlation coefficient (MCC). If \code{NULL} then computation is skipped. Default = “F1Score”.
#' @param Entro_tresh logical. If \code{TRUE}, calculate per gene entropy, utilizing the derived per gene cutoff as true-label, to assess label distribution across clusters. Default = \code{TRUE}, requires \code{run_cuoff}.
#' @param Entro_med logical. If \code{TRUE}, calculate per gene median expression per cluster and fraction of individual median of summed medians across clusters, which is used to calculate per gene entropy. Default = \code{F}, requires \code{run_cuoff}.
#' @param run_enriched logical. If \code{TRUE}, run enrichment analysis using fisher.test. Using cutoff, expression per gene is binarized across cells. Cells have either 1 or 0 expression. Expression is summed within clusters and enrichment per cluster is calculated for each gene using fisher.test. If cluster has enrichment for a gene( p-value < 0.05), the value per gene of a cluster is set to 1. In order to speed up computation, for each gene, fraction of positive cells within a cluster is ordered in decreasing order and enrichment is tested iterativelly along that order. If enrichment p-value of 3 clusters (flag count) is not significant, the remaining clusters are expected to be not enriched. Cluster with less cells than the number of average cells per cluster do not increase the flag count.
#' @param give2ndfiff logical. If \code{TRUE}, run differential expression analysis between every cluster and its closest cluster(s) based on highest number of co-enriched genes, for genes which are shared enriched in these clusters. If more than one cluster share the same number of co-enriched genes, differential expression of co-enriched genes is performed for all co-enriched clusters. Default = \code{T}. Co-enriched clusters can represent cell states of the same cell types.
#' @param diffexp either “nbino” or “wilcox”. Performs differential expression analysis between cells of clusters with highest number of co-enriched genes for these co-enriched genes based on Wilcoxon test or negative binomial distribution test utilizing global gene mean-variance dependence. Default = “nbino”.
#' @param vfit function of the background noise model describing the expression variance as a function of the mean expression. Input can be utilized for differential expression analysis between co-enriched genes and identification of outlier gene-expression within cluster in outlier analysis. Default = \code{NULL}.
#' @param gooutlier  logical. If \code{TRUE}, performs outlier identification based on cluster partition and identifies outlier gene expression within clusters.
#' @param individualfit logical. If \code{TRUE}, background noise model, required to infer outlier expression, is fitted for each cluster separately, default = \code{F}.
#' @param outminc integer value, minimal transcript count of a gene to be included in the background fit.
#' @param probthr integer value, outlier probability threshold for genes to exhibit outlier expression within a cluster. Probability is computed from a negative binomial background model of expression in a cluster.
#' @param diptest logical. If \code{T}, performs dip.test function from the diptest package to test for unimodality of gene expression (enriched genes) within clusters by computing Hartigans’ dip statistics per gene. Calculation is performed only on expression values with at least \code{minexpr}. As calculating is performed on library size normalized and rescaled data, \code{minexpr} is rescaled basd on scalefactor divided by \code{mintotal}. Expression is only tested, if a given fraction of a cluster, \code{unifrac}, exhibits minimal expression of rescaled \code{minexpr} or the sample size equals at least \code{clustsize}.Default = \code{T}.
#' @param bwidth  logical. If \code{T}, performs Silverman’s critical bandwidth method to test for unimodality of gene expression (enriched genes) within clusters. Calculation is performed only on expression values with at least \code{minexpr}. As calculating is performed on library size normalized and rescaled data, \code{minexpr} is rescaled based on scalefactor divided by \code{mintotal}. Expression is only tested, if a given fraction of a cluster, \code{unifrac}, exhibits minimal expression of rescaled \code{minexpr} or the sample size equals at least \code{clustsize}.Default = \code{T}.
#' @param critmass logical. If \code{T}, performs Ameijeiras-Alonsos’s method to test for unimodality of gene expression (enriched genes) within clusters. Calculation is performed only on expression values with at least \code{minexpr}. As calculating is performed on library size normalized and rescaled data, \code{minexpr} is rescaled based on scalefactor divided by \code{mintotal}. Expression is only tested, if a given fraction of a cluster, \code{unifrac}, exhibits minimal expression of rescaled \code{minexpr} or the sample size equals at least \code{clustsize}.Default = \code{T}.
#' @param mintotal minimal number of transcripts cells are expected to have, to calculate expression cutoff. Default = 3000
#' @param unifrac fraction of cluster required to exhibit at least scaled \code{minexpr} that gene is tested for unimodality. Default = 0.1.
#' @param logmodetest  logical. If \code{T}, performs log transformation before testing unimodality of gene expression. Default = \code{F}.
#' @param b_bw number of replicates used for Silverman’s critical bandwith test, default = 25.
#' @param n_bw number of equally spaced points at which density is estimated, for Silverman’s critical bandwith test, default = 128.
#' @param b_ACR number of replicates used for Ameijeiras-Alonsos’s unimodality test, default = 100.
#' @param n_ACR number of equally spaced points at which density is estimated, for Ameijeiras-Alonsos’s unimodality test, default = 1024.
#' @param batch_entrop logical. If \code{T}, calculate the entropy of batches across cluster. Default = \code{F}.
#' @param set.name set name for individual assessment within output of list of assessments. Default = \code{NULL} and name is given in the following way: if \code{seuratobject}, name is selected from metadata columns equal to Idents(), or character string given as input for givepart or character string of object name of numeric cluster partition. If \code{RaceIDobject}, name is given by character string given as input for givepart, character string of the object name of the number cluster partition or “Vdefault”.
#' @param rawdata_null logical. If \code{TRUE}, do not store UMI count table in output of assessment, default = T
#' @return List of assessments, with a named slot per assessment. Individual assessments represent a list with the following slots:
#'   \item{rawdata}{Raw expression data matrix/UMI count matrix derived from input objects, with cells as columns and genes as rows in sparse matrix format.}
#'   \item{rowmean}{mean expression of assessed features.}
#'   \item{part}{vector containing cluster partition derived from input objects.}
#'   \item{clustsize}{threshold of minimum number of cells in a cluster used for assessment.}
#'   \item{features}{vector of feature genes derived from object, used to compute its cluster partition.}
#'   \item{assessed_features}{vector of features assessed through assess me function, can differ from \code{features} when \code{var_feat_len} argument differs from length of object derived features or different set of genes given as argument with \code{givefeatures}}
#'   \item{PCA}{data.frame with 4 columns, indicating top PCAnum genes with: highest loadings for PC1, lowest loadings for PC1, highest loadings for PC2 and lowest loadings for PC2.}
#'   \item{max_cl}{vector indicating for assessed features which cluster exhibits highest mean expression.}
#'   \item{cutoff}{vector indicating calculated numeric cutoff for assessed features.}
#'   \item{f1_score}{vector indicating f1_score or alternative statistical analysis for binary classification, for the assessed features.}
#'   \item{Entropy_tresh}{vector indicating Entropy per assessed feature, calculated based on the per gene cutoff.}
#'   \item{Entropy_median}{ector indicating Entropy per assessed feature, calculated based on per gene median expression per cluster and fraction of individual medians of summed median across clusters.}
#'   \item{cluster}{vector indicating assessed clusters.}
#'   \item{enriched_features}{number of enriched features per cluster.}
#'   \item{enriched_feature_list}{list with a vector per cluster of enriched features.}
#'   \item{unique_features}{number of uniquely enriched features per cluster.}
#'   \item{unique_feature_list}{list with a vector per cluster of uniquely enriched features.}
#'   \item{second_cluster}{data.frame with rows representing a cluster and its closest clusters based on co-enriched genes and columns representing: "frac_shared_to_clos_cluster” = number of co-enriched genes,“rel_frac_shared_to_clos”: fraction of co-enriched genes of enriched genes,“frac_diff_of_shared_features “: number of differential genes of co-enriched genes,“rel_frac_diff_of_shared_to_clos”: fraction of differential genes of co-enriched genes}
#'   \item{list_2ndShared}{list with data.frame for every cluster with rows as enriched genes of a cluster and columns representing binary classification for enrichment (1= enriched, 0 = not enriched) of  a cluster and its  most similar clusters based on co-enriched genes.}
#'   \item{shared2ndgenes}{list with vector for every cluster of enriched genes with co-enrichment in closest clusters.}
#'   \item{list_2nd_diff}{list with vector for every cluster of co-enriched genes with differential expression to co-enriched clusters.}
#'   \item{outliertab}{data.frame indicating number of outlier cells per cluster with 1, 2 or 3 outlier genes. Rows representing cluster and columns representing number of cells with 1, 2 or 3 outlier genes.}
#'   \item{outlier_genes}{list with vector for every clusters indicating outlier genes.}
#'   \item{nonunimodal_list}{list with data.frame per cluster with rows representing enriched gene per cluster and columns p.value of dip.test and p.value after multiple testing correction with Bonferroni and BH method.}
#'   \item{nonunimodaltab}{data.frame indicating number of genes per cluster with non-unimodal expression before and after multiple-testing correction.}
#'   \item{bandwidth_list}{list with vector for every cluster indicating genes with non-unimodal expression derived from Silverman’s critical bandwith test.}
#'   \item{masstest_list}{list with vectors for every cluster indicating gene with non-unimodal expression based on Ameijeiras-Alonsos’s method to test for unimodality.}
#'   \item{batch_entropy}{entropy of batches across clusters}
#' @examples
#' entero <- CreateSeuratObject(counts = x, project = "10x", min.cells = 3, min.features = 200)
#' entero <- NormalizeData(entero, normalization.method = "RC", scale.factor = 10000)
#' entero <- FindVariableFeatures(entero, selection.method = "vst", nfeatures = 3000)
#' features <- Seurat::VariableFeatures(entero)
#' entero <- ScaleData(entero, features = features)
#' entero <- RunPCA(entero, features = features, npcs = 100)
#' entero <- FindNeighbors(entero, dims = 1:100)
#' resolution <- c(1:10)
#' for (i in resolution)  { entero  <- FindClusters(entero , resolution = i) }
#' res <- colnames(entero[[]])[c(4,6:length(colnames(entero[[]])))]
#' for (i in 1:length(res)) {if (i == 1) { assess_seuratRC <- cluster_assessment( seuratobject=entero,givepart = res[i], give2ndfiff=F, Entro_med=F, diptest=F, run_enriched=T, bwidth=F, critmass=F, gooutlier=T) } else { assess_seuratRC <- cluster_assessment(assessment_list = assess_seuratRC, seuratobject=entero,givepart = res[i], give2ndfiff=F, Entro_med=F, diptest=F, run_enriched=T, bwidth=F, critmass=F, gooutlier=T) }}
#' @export
cluster_assessment <- function(assessment_list=NULL, seuratobject =NULL, seurat_assay="RNA", seurat_lib_size=F,do.features=T,var_feat_len=NULL, RaceIDobject=NULL, RaceID_cl_table=NULL,ScanpyobjectFullpath=NULL, scanpy_clust = "leiden", scanpyscalefactor=10000,rawdata=NULL,ndata=NULL, norm=T, givepart=NULL,givefeatures=NULL,minexpr=5,CGenes=NULL,ccor=0.65,fselectRace=F,fselectSeurat=F,givebatch=NULL,individualbatch=NULL,gene.domain=F,PCA_QA=F,PCAnum=10,run_cutoff=T ,f1Z = F,cutoff="mean",cutoffmax=F,clustsize=10, binaclassi="F1Score", Entro_tresh=T, Entro_med=T,run_enriched=T,give2ndfiff=T,diffexp="nbino",vfit=NULL,gooutlier=T,individualfit=F, outminc=5, probthr=0.01,diptest=T, bwidth=T, critmass=T,mintotal=3000,unifrac=0.1,  logmodetest=F,b_bw=25,n_bw=128,b_ACR=100,n_ACR=1024,  batch_entropy=F,   set.name=NULL,rawdata_null=T){ # mintotal=3000, mingene=200) not used as cluster partition should be there {
  start_time_assessment_complete <- Sys.time()
  if (!cutoff %in% c("mean", "median")) { stop("cutoff must be either mean or median")}
  if (sum(c(!is.null(seuratobject), !is.null(RaceIDobject), !is.null(ScanpyobjectFullpath))) > 1) { stop("please set only one object: either Seurat or RaceID or Scanpy")}
  if (class(givepart) =="character") {
    if (is.null(seuratobject) && is.null(RaceID_cl_table)) { stop("either seurat object or RaceID metadata dataframe/cluster dataframe(RaceID_cl_table can't be NULL) has to be provided") }
  }
  cat(paste("run_extraction", "\n"))
  start_time <- Sys.time()
  features <- NULL
  g <- NULL
  if ( !is.null(seuratobject)) {
    if (! seurat_assay %in% names(seuratobject)) { stop("pleaes set seurat assay to e.g. 'RNA'")}
    else {
      rawdata <- seuratobject[[seurat_assay]]@counts  ###  e.g. seuratobject@assays$RNA
      if (seurat_lib_size==T) {
        ndata <- t(t(rawdata)/colSums(rawdata))
      }
      else {
        ndata <- seuratobject[[seurat_assay]]@data}
      if (is.null(givepart)) {
        part <- as.numeric(seuratobject$seurat_clusters)} ### needs to be corrected later for 0-notation used by seurat
      else {
        if (class(givepart) =="character") {part <- as.numeric(seuratobject@meta.data[,givepart]) }
        else {part <- as.numeric(givepart)}

      }
      names(part) <- colnames(seuratobject)
      features <- seuratobject[[seurat_assay]]@var.features
      if (do.features==T) {
        if (eval(parse(text = paste(c("length(", paste0(as.character(substitute(seuratobject)), "@assays$", seurat_assay, "@var.features"), ")"), collapse = ""))) == 0 || !is.null(var_feat_len) && length(features) != var_feat_len) {

          if (!inherits(x = seuratobject, "dgMatrix")) {
            seuratobject2 <- seuratobject[[seurat_assay]]@counts
          }
          if (!inherits(x = seuratobject2, what = "dgCMatrix")) {
            seuratobject2 <- Matrix(as.matrix(seuratobject2), sparse = T)
          }
          clip.max <- sqrt(x = ncol(x = seuratobject2))
          hvf.info <- data.frame(mean = rowMeans(x = seuratobject2))
          hvf.info$variance <- SparseRowVar2(mat = seuratobject2, mu = hvf.info$mean,
                                             display_progress = T)
          hvf.info$variance.expected <- 0
          hvf.info$variance.standardized <- 0
          not.const <- hvf.info$variance > 0
          fit <- loess(formula = log10(x = variance) ~ log10(x = mean),
                       data = hvf.info[not.const, ], span = 0.3)
          hvf.info$variance.expected[not.const] <- 10^fit$fitted
          hvf.info$variance.standardized <- SparseRowVarStd(mat = seuratobject2,
                                                            mu = hvf.info$mean, sd = sqrt(hvf.info$variance.expected),
                                                            vmax = clip.max, display_progress = T)
          colnames(x = hvf.info) <- paste0("vst.", colnames(x = hvf.info))
          hvf.info <- hvf.info[which(x = hvf.info[, 1, drop = TRUE] !=
                                       0), ]  ### mean expression >0
          hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized,
                                     decreasing = TRUE), , drop = FALSE]
          feature_genes <- head(rownames(hvf.info), var_feat_len)
        }
        else {feature_genes <- features }
      }
      else {
        feature_genes <- features
      }
    }
  }
  else if ( !is.null(RaceIDobject)) {
    rawdata <- RaceIDobject@expdata
    ndata <- RaceIDobject@ndata
    rawdata <- rawdata[,colnames(ndata)]
    if (is.null(givepart)) {
      part <- RaceIDobject@cpart
      if ( length(part) == 0) { part <- RaceIDobject@cluster$kpart }
    } ### needs to be corrected later for 0-notation used by seurat
    else {
      if (class(givepart) =="character") {
        if(is.null(RaceID_cl_table)) { stop("please include table with different resolutions")}
        else {
          part <- as.numeric(RaceID_cl_table[,givepart])
          names(part) <- rownames(RaceID_cl_table)
        }
      }
      else{
        part <- givepart }
    }
    feature_genes <- RaceIDobject@cluster$features
    features <- RaceIDobject@cluster$features
    vfit <- RaceIDobject@background$vfit ###continue integration
    minexpr <- RaceIDobject@filterpar$minexpr
    mintotal <- RaceIDobject@filterpar$mintotal
  }
  else if (!is.null(ScanpyobjectFullpath)) {
    Convert(ScanpyobjectFullpath, dest = "h5seurat", overwrite = TRUE)
    scanpytoseurat <- sub("h5ad", "h5seurat", ScanpyobjectFullpath)
    scanpy_object <- LoadH5Seurat(scanpytoseurat)
    rawdata <- scanpy_object@assays$RNA@counts
    ndata <- t(t(rawdata)/colSums(rawdata) * scanpyscalefactor)

    if ( ! scanpy_clust %in% c("leiden", "louvain") ) { stop("set either leiden or louvain ")}
    else { if  (scanpy_clust == "leiden") { part <- as.numeric(scanpy_object$leiden)}
      if (scanpy_clust == "louvain") { part <- as.numeric(scanpy_object$louvain)}
    }
    names(part) <- colnames(rawdata)
    feature_genes <- rownames(scanpy_object@assays$RNA@scale.data)
  }
  else {
    if (is.null(rawdata)) { stop("please provide rawcounts as rawdata")}
    if (is.null(givepart)) {stop("please provide cluster partition")}
    if ( norm==T) {
      cs <- apply(rawdata, 2, sum)
      ndata <- t(t(as.matrix(rawdata))/cs)
    }
    else {
      if ( is.null(ndata)) stop("provide normalised data as ndata")
    }
    rawdata <- Matrix(as.matrix(rawdata), sparse = T)
    ndata <- Matrix(as.matrix(ndata), sparse = T)
    if (is.null(givefeatures)){  ### feature selection only RaceID
      if (sum( as.numeric(fselectRace) + as.numeric(fselectSeurat)) == 2) { stop("either feature selection RaceID or Seurat")}
      if (sum( as.numeric(fselectRace) + as.numeric(fselectSeurat)) == 1) {
        fg <- rownames(rawdata)[apply(rawdata, 1, max , na.rm=T) >= minexpr]
        if (!is.null(CGenes)) {
          CGenes <- CGenes[CGenes %in% fg]
          h <- NULL
          if (length(CGenes) == 1) {
            k <- cor(as.matrix(t(ndata[fg, ])),
                     as.matrix(ndata[CGenes, ]))
          }
          else {
            k <- cor(as.matrix(t(ndata[fg, ])),
                     as.matrix(t(ndata[CGenes, ])))

          }
          h <- apply(abs(k), 1, max, na.rm = TRUE) < ccor
          h[is.na(h)] <- TRUE
          if (!is.null(h))
            genes <- fg[h]
        }
        else {
          genes <- fg
        }
        if (fselectRace) {

          bg <- fitbackground(as.matrix(ndata * min(apply(rawdata, 2, sum)))[genes, ] + 0.1)
          feature_genes <- bg$n
        }
        else if (fselectSeurat){

            if (is.null(var_feat_len)) stop("var_feat_len has to be set")
            if (!inherits(x = rawdata, "dgMatrix")) {
              seuratobject2 <- Matrix(as.matrix(rawdata), sparse = T)
            }
            if (!inherits(x = seuratobject2, what = "dgCMatrix")) {
              seuratobject2 <- rawdata
            }
            clip.max <- sqrt(x = ncol(x = seuratobject2))
            hvf.info <- data.frame(mean = rowMeans(x = seuratobject2))
            hvf.info$variance <- SparseRowVar2(mat = seuratobject2, mu = hvf.info$mean,
                                               display_progress = T)
            hvf.info$variance.expected <- 0
            hvf.info$variance.standardized <- 0
            not.const <- hvf.info$variance > 0
            fit <- loess(formula = log10(x = variance) ~ log10(x = mean),
                         data = hvf.info[not.const, ], span = 0.3)
            hvf.info$variance.expected[not.const] <- 10^fit$fitted
            hvf.info$variance.standardized <- SparseRowVarStd(mat = seuratobject2,
                                                              mu = hvf.info$mean, sd = sqrt(hvf.info$variance.expected),
                                                              vmax = clip.max, display_progress = T)
            colnames(x = hvf.info) <- paste0("vst.", colnames(x = hvf.info))
            hvf.info <- hvf.info[which(x = hvf.info[, 1, drop = TRUE] !=
                                         0), ]  ### mean expression >0
            hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized,
                                       decreasing = TRUE), , drop = FALSE]
            feature_genes <- head(rownames(hvf.info), var_feat_len)


        }
        else{
          stop("either fselectRace or fselectSeurat has to be set TRUE")
        }
      }
      else { stop("either fselect with RaceID or Seurat")}
    }
    else { feature_genes <- givefeatures}
  }
  end_time <- Sys.time()
  cat(paste("done_extraction", "time", difftime(end_time, start_time, units="secs"), "s", "\n", sep = " "))
  ### individual batch
  if (!is.null(givebatch) & !is.null(individualbatch) ) {
    index_batch <- which(givebatch == individualbatch)
    if ( length(index_batch) == 0) { stop("individual batch has to be element of givebatch ")}
    else {
      rawdata <- rawdata[,index_batch]
      ndata <- ndata[,index_batch]
      part <- part[index_batch]


    }
  }

  if (gene.domain==T) {
    if (!is.null(seuratobject) | !is.null(RaceIDobject) | !is.null(ScanpyobjectFullpath)) {
      g <- apply(rawdata, 1, max , na.rm=T) >= minexpr
      feature_genes <- rownames(rawdata)[g]
    }
    if (f1Z == F) { ### if cutoff is NOT expression = x>0
      feature_data <- ndata[feature_genes,]
    }
    else {
      feature_data <- rawdata[feature_genes,]
    }
  }
  else {
    if (is.null(givefeatures) & length(feature_genes) == 0 ) {
      stop("features need to provided or calculated")}
    else{
      if (!is.null(givefeatures)) {
        feature_genes <- rownames(rawdata)[rownames(rawdata) %in% givefeatures]
        givefeatures <- feature_genes
      }
      else{
        feature_genes <- feature_genes
      }
    }
    if (f1Z == F) {
      feature_data <- ndata[rownames(ndata) %in% feature_genes,]
    }
    else {
      feature_data <- rawdata[rownames(rawdata) %in% feature_genes,]
    }
  }

  if (PCA_QA){
    cat(paste("run_PCA", "\n"))
    pca_an <- prcomp(t(as.matrix(rawdata[feature_genes,])), center = T, scale. = T, rank. = 2)
    features_pc1 <- feature_genes[order(pca_an$rotation[,1], decreasing = T)]
    features_pc2 <- feature_genes[order(pca_an$rotation[,2], decreasing = T)]
    pctab <- cbind(pc1_head=head(features_pc1,PCAnum), pc1_tail=tail(features_pc1,PCAnum)[PCAnum:1], pc2_head=head(features_pc2, PCAnum), pc2_tail=tail(features_pc2,PCAnum)[PCAnum:1])
    cat(paste("done_PCA", "\n"))
  }
  else {pctab <- NULL}

  if (run_cutoff) {
    cat(paste("run_cutoff_table", "\n"))
    start_time <- Sys.time()

    cat(paste("run_aggregate_mean/median_cluster", "\n"))
    start_time2 <- Sys.time()

    if (cutoff == "mean") {
      agg <- aggregate.Matrix(t(feature_data), list(part), fun = "sum")[as.numeric(which(table(part) >= clustsize)),]
      agg <- agg/as.numeric(table(part)[as.numeric(which(table(part) >= clustsize))])
    }

    if ( cutoff == "median") {
      agg <- aggregate(t(as.matrix(feature_data)), list(part), median)[as.numeric(names(table(part)[table(part) >= clustsize])),]
      agg <- agg[,-1]
    }
    end_time <- Sys.time()
    cat(paste("done_aggregate","time",difftime(end_time, start_time2, units="secs"),"s", "\n", sep = " "))

    max_cl <- apply(agg, 2, function(x) { which(x == max(x))})

    ### if multiple cluster share the same maximal mean expression for a gene, kick that gene out of the feature list
    if (is.list(max_cl)) {
      kick_out <- lapply(max_cl, length)
      kick_out2 <- Reduce(append, kick_out)
      kick_int <- which(! kick_out2==1)
      agg <- agg[,-kick_int]
      feature_data <- feature_data[-kick_int,]
      feature_genes <- feature_genes[-kick_int]
      max_cl <- rownames(agg)[apply(agg, 2, function(x) { which(x == max(x))})]
    }
    else {
      max_cl <- rownames(agg)[max_cl]
    }
    names(max_cl) <- colnames(agg)

    if ( !cutoff %in% c("mean", "median")) stop("cutoff either mean or median")

    cat(paste("calculate_cutoff_per_gene", "\n"))
    start_time2 <- Sys.time()
    if ( cutoff == "mean") {
      cutoff <- apply(agg, 2, function(x) {
        maxi <- max(x)
        maxi_index <- as.numeric(which(x == max(x)))
        minor <- mean(x[-maxi_index])
        cutoff <- mean(c(maxi, minor))
        return(cutoff)
      })}
    else if (cutoff == "median") {
      cutoff <- apply(agg, 2, function(x) {
        maxi <- max(x)
        maxi_index <- as.numeric(which(x == max(x)))
        minor <- median(x[-maxi_index])
        cutoff <- mean(c(maxi, minor))
        return(cutoff)
      })}
    else if (cutoffmax) {
      cutoff <- apply(agg, 2, function(x) {
        maxi <- max(x)
        return(maxi)
      })}

    end_time <- Sys.time()
    cat(paste("done_calculte_cutoff","time", difftime(end_time, start_time2, units="secs"),"s","\n"))

    cat(paste("build_cutoff_data", "\n"))
    start_time2 <- Sys.time()

    position <- 1:nrow(feature_data)
    if (f1Z==T){
      cutoff_data_table <- feature_data > 0
      cutoff_data_table <- apply(cutoff_data_table, 2, as.numeric)
    }
    else{
      cutoff_data_table <- feature_data > cutoff
      cutoff_data_table <- apply(cutoff_data_table, 2, as.numeric)
    }

    rownames(cutoff_data_table) <- rownames(feature_data)
    ### if gene counts in all cells are above 0
    over0 <- as.numeric(which(rowSums(cutoff_data_table) == ncol(cutoff_data_table)))
    if (length(over0) > 0 ) {
      cutoff_data_table <- cutoff_data_table[-over0,]
      feature_genes <- feature_genes[-over0]
      max_cl <- max_cl[-over0]
      cutoff <- cutoff[-over0]
    }

    end_time <- Sys.time()
    cat(paste("done_build_cutoff_data","time",difftime(end_time, start_time2, units="secs"),"s", "\n"))
    end_time <- Sys.time()
    cat(paste("done_cutoff_table", "time",difftime(end_time, start_time, units="secs"),"s","\n"))
  }
  else{
    max_cl <- NULL
  }

  if (!is.null(binaclassi)) {
    if( run_cutoff ==F) { stop("F1Score/Cohenkappa/MCC calculations requires prior cutoff calculations")}
    if (!binaclassi %in% c("F1Score", "Cohenkappa", "MCC" )) { stop("statistical analysis for binary classification must be either F1Score or Cohenkappa or MCC")}
    cat(paste(binaclassi, "\n"))
    start_time <- Sys.time()

    getAUC <- function(position, data, part, max_part, F1=T, Kap=F, MCC=F) {
      if (sum(F1, Kap, MCC) != 1) { stop("set one, and only one of F1, Kap or MCC, True")}
      predictions <- as.numeric(part == max_part[position])
      labels <- as.numeric(data[position,])
      df <- table(predictions, labels)
      TP <- df[2,2]
      TN <- df[1,1]
      FP <- df[2,1]
      FN <- df[1,2]
      if (F1) {
        P <- TP/(TP+FP)
        R <- TP/(TP+FN)
        F1 <- 2*P*R/(P+R)

      }
      else if (Kap) {
        kappa <- 2*(TP*TN-FN*FP)/((TP+FP)*(FP+TN)+(TP+FN)*(FN+TN))
      }
      else if (MCC) {
        MCC <- mcc(predictions, labels)
      }
    }
    #plan(multisession)
    position2 <- 1:length(max_cl)
    if (binaclassi == "F1Score") {f1_score <- sapply(position2, getAUC, cutoff_data_table, part, as.numeric(max_cl))}
    if (binaclassi == "Cohenkappa") {f1_score <- sapply(position2, getAUC, cutoff_data_table, part, as.numeric(max_cl), F1=F, Kap=T, MCC=F)}
    if (binaclassi == "MCC") {f1_score <- sapply(position2, getAUC, cutoff_data_table, part, as.numeric(max_cl), F1=F,Kap=F, MCC=T)}

    names(f1_score) <- rownames(cutoff_data_table)
    end_time <- Sys.time()
    cat(paste("done_f1","time",difftime(end_time, start_time, units="secs"),"s","\n"))}
  else {f1_score <- NULL}

  ### entropy cutoff
  if (Entro_tresh) {
    if( run_cutoff ==F) { stop("Entropy calculations requires prior cutoff calculations")}
    cat(paste("run_entropy", "\n"))
    start_time <- Sys.time()
    getEnropy <- function(position, data, part) {
      uni_part <- sort(unique(part))
      Entropy <- c()
      labels <- data[position,]
      for ( i in 1:length(uni_part)) {
        pos_clus <- sum(data[position, part == uni_part[i]])
        p <- pos_clus/sum(labels) ### how many of the positive are in the cluster what about how many are positive in the cluster or
        Entropy <- c(Entropy, p*log(p, base = length(uni_part)))
      }

      return(sum(Entropy[!is.na(Entropy)])*(-1))
    }
    position2 <- 1:length(max_cl)

    Entropy_thresh <- sapply(position2, getEnropy, cutoff_data_table, part )
    names(Entropy_thresh) <- rownames(cutoff_data_table)
    end_time <- Sys.time()
    cat(paste("done_entropy_tresh","time", difftime(end_time, start_time, units="secs"),"s","\n"))
  }
  else {Entropy_thresh <- NULL }
  #
  if (Entro_med) {
    if( run_cutoff ==F) { stop("Entropy calculations requires prior cutoff calculations")}
    start_time <- Sys.time()
    getEntropymed <- function(position, data, part) {
      uni_part <- sort(unique(part))
      median_clus <- c()
      Entropy <- c()
      expression <- as.numeric(data[position,])
      for ( i in 1:length(uni_part)) {
        median_clus <- c(median_clus, median(expression[part == uni_part[[i]]]))
      }
      p <- median_clus/sum(median_clus)
      for ( i in 1:length(p)) {
        Entropy <- c(Entropy, p[i]*log(p[i], base = length(uni_part)))
      }
      return(sum(Entropy[!is.na(Entropy)]) *(-1))
    }
    Entropy_median <- sapply(position, getEntropymed, feature_data, part )
    names(Entropy_median) <- rownames(cutoff_data_table)
    end_time <- Sys.time()
    cat(paste("done_entropy_med", "time",difftime(end_time, start_time, units="secs"),"s","\n"))}
  else {
    Entropy_median <- NULL
  }
  ### sum feature genes per cluster # this sounds as we can work with it
  if (run_enriched) {
    if( run_cutoff ==F) { stop("Enrichment calculations requires prior cutoff calculations")}
    cat(paste("run_aggregate_cutoff", "\n"))

    start_time <- Sys.time()

    cutoff_data_table2 <- Matrix(cutoff_data_table, sparse = T)
    agg2 <- aggregate.Matrix(t(cutoff_data_table2), list(part), fun = "sum")[as.numeric(which(table(part) >= clustsize)),]
    end_time <- Sys.time()

    end_time <- Sys.time()
    cat(paste("done_aggregate_cutoff", "time",difftime(end_time, start_time, units="secs"),"s","\n"))
    cat(paste("run_enriched_features", "\n"))
    start_time <- Sys.time()

    #### calculated for every gene the enrichment p-value per cluster based on number of cells in a cluster above the cutoff versus all positive cells

    calculate_enrichment <- function(gene,pval, part, clustsize) {
      table_part <- table(part)[table(part) >= clustsize]
      pop <- sum(table_part)
      b <- sum(gene)
      fraction <- gene/as.numeric(table_part)
      ordering <- order(fraction, decreasing = T)
      pvalues <- rep(0, length(gene))
      flag <- 0
      for (i in 1:length(pvalues)){
        if (flag < 3) { #### length of the flag can be optional depending on the differentce in cluster sizes
          pvalue_row <- fisher.test(matrix(c(pop-as.numeric(table_part)[ ordering[i]], as.numeric(table_part)[ordering[i]], b - gene[ordering[i]], gene[ordering[i]]), ncol = 2), alternative = "greater")
          pv <- pvalue_row$p.value
          if (pv < pval) {
            pvalues[ordering[i]] <- 1}
          else if ( as.numeric(table_part[ordering[i]]) > mean(as.numeric(table_part))){
            flag <- flag+1}
        }
        else { break;}

      }
      return(pvalues)
    }
    plan(multisession)
    pval_data_table <-  future_apply(agg2, 2, calculate_enrichment, 0.05, part,  clustsize)
    rownames(pval_data_table) <- rownames(agg2)
    pval_data_table <- t(pval_data_table)

    end_time <- Sys.time()
    cat(paste("done_enriched_features", "time",difftime(end_time, start_time, units="secs"),"s","\n"))
    enriched_features <- colSums(pval_data_table) #### features genes enriched in cluster
    enriched_feature_list <- apply(pval_data_table, 2, function(x) {
      rownames(pval_data_table)[x == 1]
    })
    names(enriched_feature_list) <- paste("cl", rownames(agg2))
    ### unique feature genes
    unique_features <- colSums(pval_data_table[apply(pval_data_table, 1, sum) == 1,]) ### feature genes enriched in only 1 cl

    ### unique features
    unique_pval_table <- pval_data_table[apply(pval_data_table, 1, sum) == 1,]

    unique_feature_list <- apply(unique_pval_table,2, function(x) {
      rownames(unique_pval_table)[x == 1]
    })
    names(unique_feature_list) <- paste("cl", rownames(agg2))

  }
  else {
    unique_feature_list <- NULL
    unique_features <- NULL
    enriched_feature_list <- NULL
    enriched_features <- NULL
  }

  if(give2ndfiff) {
    if( run_cutoff ==F) { stop("2nd differential calculations requires prior cutoff calculations")}
    if( run_enriched ==F) { stop("2nd differential calculations requires prior enrichment calculations")}
    if (!diffexp %in% c("nbino", "wilcox")) { stop("diffexp must be either set to wilcox or nbino")}
    start_time <- Sys.time()
    cat(paste("run_2nd_diff", "\n"))
    give_2nd_fraction <- function(dataX, diffexp, vfit, gene.domain, part, rawdata, minexpr, ndata, g) {
      pos <- apply(dataX, 2, function(x){
        indices <- which(x==1)   ## indices for enriched feature genes for a cluster
      })
      if ( !is.null(seuratobject) && seurat_lib_size == F){
        ndata <-  t(t(rawdata)/colSums(rawdata))
        #fdata <- ndata * median(colSums(rawdata)) ## redundant
      }
      if (diffexp=="nbino") {
        if ( !is.null(RaceIDobject) | sum(colSums(as.matrix(ndata))) == ncol(ndata) ) {
          fdata <- ndata * median(colSums(rawdata))
        }
        else {
          fdata <- ndata
        }

        vf <- function(x, fit) 2^(coef(fit)[1] + log2(x) *
                                    coef(fit)[2] + coef(fit)[3] * log2(x)^2)
        sf <- function(x, vf, fit) x^2/(max(x + 1e-06, vf(x, fit)) -
                                          x)
        if ( is.null(vfit )) {
          if (gene.domain==F & is.null(g)) {
            g <- apply(rawdata, 1, max , na.rm=T) >= minexpr
          }
          bg <- RaceID:::fitbackground(as.matrix(rawdata[g,names(part)]))
          fit <- bg$fit ### for background model
          vfit <- fit

        }
      }
      fraction_2nd <- c()
      fraction_2nd_rel <- c()
      diffexp_2nd <- c()
      diffexp_2nd_rel <- c()
      cluster_2nd <- c()
      list_2nd <- list()
      list_2nd_diff <- list()
      for ( i in 1:ncol(dataX)) {
        data_pos <- dataX[as.numeric(pos[[i]]),]   ### subset for positive genes for the p-value cluster table
        data_pos_cs <- apply(data_pos[,-i], 2, sum) ### sum of positive genes for all other clusters to find the 2nd highest overlap
        max_overlap <- max(data_pos_cs) ### max overlap
        data_pos_2nd <- names(data_pos_cs[which(data_pos_cs == max(data_pos_cs))]) ### most similar cluster based on enriched genes
        cat(paste(data_pos_2nd, "\n", sep = ""))
        if ( length(data_pos_2nd) > 1) { cluster_2nd <- c(cluster_2nd, paste("cl", paste(data_pos_2nd, collapse = "&"), sep=""))}
        else { cluster_2nd <- c(cluster_2nd, paste("cl", data_pos_2nd, sep = ""))}


        common_genes <- rownames(data_pos)[apply(data_pos[,data_pos_2nd, drop=F], 1, max) > 0] ## also considers if more than two cluster share differential expressed genes
        fraction <- length(common_genes)

        fraction_2nd_rel <- c(fraction_2nd_rel, fraction/nrow(data_pos)) ### relative fraction of shared genes between 2nd closes
        fraction_2nd <- c(fraction_2nd, fraction) ### number of shared genes
        data_2nd <- data_pos[,c(colnames(data_pos)[i], data_pos_2nd)] ## subsetting for iterating cluster and most similar
        list_2nd[[i]] <- data_2nd
        names(list_2nd)[i] <- paste("cl", colnames(dataX)[i])
        if ( length(data_pos_2nd) > 1) {
          if ( diffexp == "wilcox") {
            pval_list <- list()
            cat(paste("shared_enriched_genes", "\n", sep=""))
            for ( k in 1:length(data_pos_2nd)) {
              wilcox_pval <- c()
              cat(paste(data_pos_2nd[k], "\n", sep=""))
              for ( n in 1:length(common_genes)) { ### to compare whether there is a signficant difference between the counts of the two, although enriched
                #cat(paste0(common_genes[n], "\n"))
                neg_nums <- as.numeric(ndata[common_genes[n],names(part)[part == as.numeric(colnames(data_2nd)[1])]])
                pos_nums <- as.numeric(ndata[common_genes[n],names(part)[part == as.numeric(colnames(data_2nd)[k+1])]])
                comp <- wilcox.test(neg_nums, pos_nums, alternative = "two.sided")
                wilcox_pval <- c(wilcox_pval, comp$p.value)
              }
              pval_list[[k]] <- wilcox_pval
            }
            wilcox_pval <- do.call(rbind, pval_list)
            diffexp_2nd_rel <- c(diffexp_2nd_rel, sum(apply(wilcox_pval, 2, min) < 0.05)/fraction_2nd[i])   ### relative fraction of shared genes differential expressed
            diffexp_2nd <- c(diffexp_2nd, sum(apply(wilcox_pval, 2, min) < 0.05))
            list_2nd_diff[[i]] <- common_genes[apply(wilcox_pval, 2, min) < 0.05]
            names(list_2nd_diff)[i] <- paste("cl", colnames(dataX)[i])
          }
          if (diffexp == "nbino") {
            pval_list <- list()
            m <- list()
            m[[1]] <- rowMeans(fdata[common_genes,names(part)[part == as.numeric(colnames(data_2nd)[1])]])
            #v[[1]] <- matrixStats::rowVars(ndata[common_genes,names(part)[part == as.numeric(colnames(data_2nd)[1])]])
            for ( k in 1:length(data_pos_2nd)) {
              cat(paste(data_pos_2nd[k], "\n", sep=""))
              m[[2]] <- rowMeans(fdata[common_genes,names(part)[part == as.numeric(colnames(data_2nd)[k+1])]])
              #v[[2]] <- matrixStats::rowVars(ndata[common_genes,names(part)[part == as.numeric(colnames(data_2nd)[k+1])]])
              pv <- apply(data.frame(m[[1]], m[[2]]), 1, function(x) {
                p12 <- (dnbinom(0:round(x[1] * sum(part == as.numeric(colnames(data_2nd)[1])) + x[2] *
                                          sum(part == as.numeric(colnames(data_2nd)[k+1])), 0), mu = mean(x) * sum(part == as.numeric(colnames(data_2nd)[1])), size = sum(part == as.numeric(colnames(data_2nd)[1])) *
                                  sf(mean(x), vf, vfit))) * (dnbinom(round(x[1] * sum(part == as.numeric(colnames(data_2nd)[1])) +
                                                                             x[2] * sum(part == as.numeric(colnames(data_2nd)[k+1])), 0):0, mu = mean(x) * sum(part == as.numeric(colnames(data_2nd)[k+1])),
                                                                     size = sum(part == as.numeric(colnames(data_2nd)[k+1])) * sf(mean(x), vf, vfit)))
                if (sum(p12) == 0)
                  0
                else sum(p12[p12 <= p12[round(x[1] * sum(part == as.numeric(colnames(data_2nd)[1])), 0) +
                                          1]])/(sum(p12))
              })
              pval_list[[k]] <- pv

            }
            wilcox_pval <- do.call(rbind, pval_list)
            diffexp_2nd_rel <- c(diffexp_2nd_rel, sum(apply(wilcox_pval, 2, min) < 0.05)/fraction_2nd[i])   ### relative fraction of shared genes differential expressed
            diffexp_2nd <- c(diffexp_2nd, sum(apply(wilcox_pval, 2, min) < 0.05))
            list_2nd_diff[[i]] <- common_genes[apply(wilcox_pval, 2, min) < 0.05]
            names(list_2nd_diff)[i] <- paste("cl", colnames(dataX)[i])
          }

        }


        else {
          if ( diffexp == "wilcox") {
            wilcox_pval <- c()
            for ( n in 1:length(common_genes)) { ### to compare whether there is a signficant difference between the counts of the two, although enriched
              #cat(paste0(common_genes[n], "\n"))
              neg_nums <- as.numeric(ndata[common_genes[n],names(part)[part == as.numeric(colnames(data_2nd)[1])]])# + rnorm(length(neg_nums))/1e+9
              pos_nums <- as.numeric(ndata[common_genes[n],names(part)[part == as.numeric(colnames(data_2nd)[2])]])# + rnorm(length(pos_nums))/1e+9
              comp <- wilcox.test(neg_nums, pos_nums, alternative = "two.sided")
              wilcox_pval <- c(wilcox_pval, comp$p.value)}
            diffexp_2nd_rel <- c(diffexp_2nd_rel, sum(wilcox_pval < 0.05)/fraction_2nd[i])   ### relative fraction of shared genes differential expressed
            diffexp_2nd <- c(diffexp_2nd, sum(wilcox_pval < 0.05)) ## number of shared genes differential
            list_2nd_diff[[i]] <- common_genes[wilcox_pval < 0.05]
            names(list_2nd_diff)[i] <- paste("cl", colnames(dataX)[i])
          }
          if (diffexp == "nbino") {
            cat(paste(data_pos_2nd, "\n", sep=""))
            m <- list()
            m[[1]] <- rowMeans(fdata[common_genes,names(part)[part == as.numeric(colnames(data_2nd)[1])]])
            m[[2]] <- rowMeans(fdata[common_genes,names(part)[part == as.numeric(colnames(data_2nd)[2])]])
            pv <- apply(data.frame(m[[1]], m[[2]]), 1, function(x) {
              p12 <- (dnbinom(0:round(x[1] * sum(part == as.numeric(colnames(data_2nd)[1])) + x[2] *
                                        sum(part == as.numeric(colnames(data_2nd)[2])), 0), mu = mean(x) * sum(part == as.numeric(colnames(data_2nd)[1])), size = sum(part == as.numeric(colnames(data_2nd)[1])) *
                                sf(mean(x), vf, vfit))) * (dnbinom(round(x[1] * sum(part == as.numeric(colnames(data_2nd)[1])) +
                                                                           x[2] * sum(part == as.numeric(colnames(data_2nd)[2])), 0):0, mu = mean(x) * sum(part == as.numeric(colnames(data_2nd)[2])),
                                                                   size = sum(part == as.numeric(colnames(data_2nd)[2])) * sf(mean(x), vf, vfit)))
              if (sum(p12) == 0)
                0
              else sum(p12[p12 <= p12[round(x[1] * sum(part == as.numeric(colnames(data_2nd)[1])), 0) +
                                        1]])/(sum(p12))
            })
            wilcox_pval <- pv
            diffexp_2nd_rel <- c(diffexp_2nd_rel, sum(wilcox_pval < 0.05)/fraction_2nd[i])   ### relative fraction of shared genes differential expressed
            diffexp_2nd <- c(diffexp_2nd, sum(wilcox_pval < 0.05)) ## number of shared genes differential
            list_2nd_diff[[i]] <- common_genes[wilcox_pval < 0.05]
            names(list_2nd_diff)[i] <- paste("cl", colnames(dataX)[i])

          }

        }
      }
      names(fraction_2nd) <- colnames(pval_data_table) ### how many genes are shared between the 2nd closest
      names(fraction_2nd_rel) <- colnames(pval_data_table)
      names(diffexp_2nd) <- colnames(pval_data_table) ### how manny of those shared are differential expressed
      names(diffexp_2nd_rel) <- colnames(pval_data_table)
      dataY <- cbind(frac_shared_to_clos_cluster=fraction_2nd, rel_frac_shared_to_clos=fraction_2nd_rel, frac_diff_of_shared_features=diffexp_2nd, rel_frac_diff_of_shared_to_clos=diffexp_2nd_rel)
      rownames(dataY) <- paste(colnames(dataX), " w ", cluster_2nd, sep = "")
      list2nd <- list(dataY=data.frame(dataY), list_2nd=list_2nd, list_2nd_diff=list_2nd_diff)
      return(list2nd)
    }

    second_cluster <- give_2nd_fraction(pval_data_table, diffexp = diffexp, part=part, vfit=vfit, gene.domain=gene.domain, rawdata=rawdata, minexpr=minexpr, ndata=ndata, g=g )

    shared2ndgenes <- lapply(second_cluster$list_2nd, function(x) {
      rownames(x)[rowSums(x) >= 2]
    })
    end_time <- Sys.time()
    cat(paste("done_2nd_cluster", "time",difftime(end_time, start_time, units="secs"),"s","\n"))
  }
  else {
    shared2ndgenes <- NULL ### for every cluster, list of enriched genes which are shared between most similar(1 or more) cluster
    second_cluster <- list()
    second_cluster[["list_2nd"]] <- NULL ### for every cluster table of enriched genes together with cluster which share most of these enriched genes
    second_cluster[["list_2nd_diff"]] <- NULL ### differential expressed genes of shared enriched genes between closest clusters
    second_cluster[["dataY"]] <- NULL   ###  total number of enriched genes shared with closest clusters,and fraction of shared enriched relative to number of enriched genes
  }                                     #### total number of differential expressed of shared,

  ### outlier cells which probthr to chose? Default is 10^-3 for RaceID, nearest neighbors with link probability < 10^-2 are discarded
  if ( gooutlier) {
    start_time <- Sys.time()
    cat(paste("run_outlier", "\n"))
    lvar  <- function(x,fit) 2**(coef(fit)[1] + log2(x)*coef(fit)[2] + coef(fit)[3] * log2(x)**2)
    lsize <- function(x,lvar,fit) x**2/(max(x + 1e-6,lvar(x,fit)) - x)
    cluster <- as.numeric(names(table(part)[table(part) >= clustsize]))
    outliergenes <- function(rawdata,data, part, clustsize, feature_genes2, vfit, outminc,minexpr,pval_outlg, individualfit, g) {
      outlier_list <- list()
      cluster <- as.numeric(names(table(part)[table(part) >= clustsize]))
      if (!is.null(g) && outminc == minexpr) {
        figenes <- g
      }
      else {
        figenes <- apply(rawdata, 1, max) >= outminc
      }
      if (individualfit==F) {

        if ( is.null(vfit )) {
          ### use RaceID function to provide the fit
          bg <- RaceID:::fitbackground(as.matrix(rawdata[figenes,names(part)]))
          fit <- bg$fit ### for background model

        }
        else {
          fit <- vfit
        }

      }
      else {
        plan(multisession)
        fit <- future_lapply(seq_along(cluster), function(x) {
          #cat(paste("fit_cluster_", cluster[x], "\n",sep = "" ))
          bg <- RaceID:::fitbackground(as.matrix(rawdata[figenes,names(part)[part == cluster[x]]]))
          vfit <- bg$fit
          vfit
        })


      }
      for (n in 1:length(cluster)) {
        cat(paste("cluster",n, "\n"))
        if (sum(part == cluster[n]) == 1) {
          cprobs <- append(cprobs, 0.5)
          names(cprobs)[length(cprobs)] <- names(part)[part ==  cluster[n]]
          next
        }
        set <- part == cluster[n]
        fdata <- rawdata[feature_genes2, as.logical(set)]
        cl <- paste("cl", n, sep = ".")
        if (individualfit==T) {
          z <- apply(fdata, 1, function(x) {
            pr <- pnbinom(round(x, 0), mu = mean(x), size = lsize(mean(x),
                                                                  lvar, fit[[n]])) ### p-value for expression of each gene
            apply(cbind(pr, 1 - pr), 1, min)
          })
          z <- t(z)
        }
        else{
          z <- apply(fdata, 1, function(x) {
            pr <- pnbinom(round(x, 0), mu = mean(x), size = lsize(mean(x),
                                                                  lvar, fit)) ### p-value for expression of each gene
            apply(cbind(pr, 1 - pr), 1, min)
          })
          z <- t(z)
        }
        z[is.na(z)] <- 1
        cp <- apply(z, 2, function(x) {
          y <- p.adjust(x, method = "BH")
          y1 <- as.numeric(y < pval_outlg)
          return(y1)
        })
        rownames(cp) <- feature_genes2
        outlier_list[[n]] <- cp
        names(outlier_list)[n] <- paste0("cl", cluster[n])
        ### cells of cluster with number of outlier genes
      }
      return(outlier_list)
    }
    outlier <- outliergenes(rawdata = rawdata,  part = part, clustsize = clustsize, outminc = outminc, minexpr = minexpr,feature_genes2 = feature_genes, vfit = vfit, pval_outlg = probthr, individualfit=individualfit,g=g)

    outliernumb <- lapply(outlier, function(x) {
      x <- x[rownames(x) %in% feature_genes,]
      y <- colSums(x)
      z <- c(outlg1=sum(y >= 1), outlg2=sum(y >= 2), outlg3=sum(y>=3))
      return(z)
    })
    outliertab <- do.call(rbind, outliernumb)
    outlier_genes <- lapply(outlier, function(x) {
      x <- x[rownames(x) %in% feature_genes,]
      rownames(x)[rowSums(x) > 0]
    })
    names(outlier_genes) <- rownames(outliertab)

    end_time <- Sys.time()
    cat(paste("done_outlier", "time",difftime(end_time, start_time, units="secs"),"s","\n"))
  }


  else {
    outliertab <- NULL
    outlier_genes <- NULL
    cluster <- NULL
  }

  if ( diptest || bwidth || critmass ) {

    if (!is.null(seuratobject)) {
      scale_factor <- seuratobject@commands$NormalizeData.RNA$scale.factor
      if (!is.null(scale_factor)) {
        unimod_tab <- t(t(rawdata) / colSums(rawdata) * seuratobject@commands$NormalizeData.RNA$scale.factor)}
      else {
        scale_factor <- median(colSums(rawdata))
        unimod_tab <- t(t(rawdata) / colSums(rawdata) * scale_factor)

      }

    }
    else if ( !is.null(RaceIDobject)) {
      unimod_tab <- ndata * mintotal; scale_factor <- mintotal
    }
    else if ( !is.null(ScanpyobjectFullpath)) {
      unimod_tab <- ndata
      scale_factor <- scanpyscalefactor
    }
    else { scale_factor <-  median(apply(rawdata, 2, sum));
    unimod_tab <- ndata *  scale_factor }
    uni_part <- as.numeric(names(table(part)))[table(part) >= clustsize]
    cat(paste("done_unimod_tab", "\n"))
  }


  if ( diptest) {
    cat(paste("run_dip", "\n"))
    pval_list <- list()


    for ( i in 1:length(uni_part)) {
      cat(paste("pval_cal",uni_part[i],"_",names(enriched_feature_list)[i], "\n"))

      unimod_tab2 <- unimod_tab[enriched_feature_list[[i]],part == uni_part[i]]

      pval <- apply(unimod_tab2, 1, function(x) {
        to_test <- x[x >= minexpr*(scale_factor/mintotal)]
        if(logmodetest) { to_test <- log(to_test)}
        if ( length(to_test) >= unifrac*sum(part == uni_part[i]) && length(to_test) >= clustsize) {
          #if ( length(to_test) >= clustsize ) {
          diptest <- dip.test(to_test)
          return(diptest$p.value)}
        else { diptest <- 1; return(diptest) }
      })
      pval_list[[paste("cl", uni_part[i], sep = "")]] <- cbind(p.val=pval, BH.pval=p.adjust(pval, method = "BH"), BF.pval=p.adjust(pval, method = "bonferroni"))
    }
    nonunimodalnumb <- lapply(pval_list, function(x) {
      return(apply(x < 0.05, 2, sum))
    })
    nonunimodaltab <- data.frame(do.call(rbind, nonunimodalnumb))
    cat(paste("done_dip", "\n"))
  }
  else{
    pval_list <- NULL
    nonunimodaltab <- NULL
  }

  ##bw
  if (bwidth) {
    cat(paste("run_bandwith_test", "\n"))
    pval_list2 <- list()
    for ( i in 1:length(uni_part)) {
      cat(paste("pval_cal_band","cl",uni_part[i], "\n"))

      unimod_tab2 <- unimod_tab[feature_genes,part == uni_part[i]]
      pval2 <- apply(unimod_tab2, 1, function(x) {
        to_test <- x[x >= minexpr*(scale_factor/mintotal)]
        if(logmodetest) { to_test <- log(to_test)}
        if ( length(to_test) >= unifrac*sum(part == uni_part[i]) && length(to_test) >= clustsize) {
          #if ( length(to_test) >= clustsize) {
          diptest <- modetest(to_test,mod0=1,method="SI",B=b_bw,submethod=1,n=n_bw,tol=1) ### retest effect of tolerance on runtime
          return(diptest$p.value)}
        else { diptest <- 1; return(diptest) }
      })
      pval_list2[[paste("cl", uni_part[i], sep = "")]] <- pval2
    }
    pval_list2 <- lapply(pval_list2, function(x) {
      x[x < 0.05]
    })
    cat(paste("done_bandwith_test", "\n"))
  }
  else { pval_list2 <- NULL}

  #### critmass
  if ( critmass) {
    cat(paste("run_critmass_test", "\n"))
    pval_list3 <- list()
    for ( i in 1:length(uni_part)) {
      cat(paste("critical mass","cl",uni_part[i], "\n"))
      pval_vec <- c()
      x = 0
      repeat {
        #unimod_tab2 <- unimod_tab[enriched_feature_list[[i]],part == uni_part[i]]
        unimod_tab2 <- unimod_tab[feature_genes,part == uni_part[i]]
        pval <- apply(unimod_tab2, 1, function(x) {
          to_test <- x[x >= minexpr*(scale_factor/mintotal)]
          if(logmodetest) { to_test <- log(to_test)}
          if ( length(to_test) >= unifrac*sum(part == uni_part[i]) && length(to_test) >= clustsize) {
            #diptest <- modetest(to_test,mod0=1,method="SI",B=25,submethod=1,n=128,tol=1)
            diptest <- modetest(to_test,mod0=1,method="ACR",B=b_ACR, n = n_ACR)
            return(diptest$p.value)}
          else { diptest <- 1; return(diptest) }
        })
        pval <- pval[pval < 0.05]
        pval_vec <- c(pval_vec, pval)
        #end_time <- Sys.time()
        x = x+1
        #end_time - start_time
        if (x == 5) {break}
      }
      pval_list3[[paste("cl", uni_part[i], sep = "")]] <- names(table(names(pval_vec))[table(names(pval_vec)) ==5])
    }
    cat(paste("done_critmass_test", "\n"))
  }
  else { pval_list3 <- NULL}
  if (batch_entropy ==T) {
    cat(paste0("run_batch_entropy", "\n"))
    if(is.null(givebatch)) {stop("please provide batch as vector")}
    batch_entropy2 <- function(part, batch, clustsize=clustsize) {
      #if ( contri_entropy==F & batch_entropy==F || contri_entropy==T & batch_entropy==T ) { stop("either batch entropy need to be set")}
      if (!identical(sub("\\_.+", "", names(part)), batch)) {
        names(part) <- paste(batch, names(part), sep = "_")
      }
      cpart <- part
      if ( is.factor(part)) { cpart <- as.numeric(as.character(cpart)); names(cpart) <- names(part)}
      cluster <- as.numeric(names(table(cpart)[table(cpart) >= clustsize]))
      cpart <- cpart[cpart %in% cluster]
      batch <- sub("\\_.+", "", names(cpart))
      batch_len <- unique(batch)
      iels <- cbind(cluster=cpart, batch=batch)
      rownames(iels) <- names(cpart)
      iels <- data.frame(iels)

      counts <- as.matrix(table(iels$batch ,iels$cluster))
      ### normalize batch per cluster count, by cluster size
      clust_norm_counts <- t(t(counts)/apply(counts,2,sum))
      ### fraction of normalized counts per cluster
      rel_over_clus_counts <- clust_norm_counts/rowSums(clust_norm_counts)

      ### batch entropy across cluster
      entropy2 <- apply(rel_over_clus_counts, 1, function(x) {
        ent <- c()
        for ( i in x) {
          ent <- c(ent, i*log(i, base = length(x)))
        }
        sum(ent[!is.na(ent)])
      })
      entropy2 <- entropy2 *(-1)
      return(entropy2)

    }
    batch_entr <- batch_entropy2(part, batch = givebatch, clustsize = clustsize )
    if (!is.null(individualbatch)){
      batch_entr <- batch_entr[individualbatch]}
    cat(paste0("done_run_batch_entropy", "\n"))
  }
  else {
    batch_entr <- NULL
  }
  end_time_assessment_complete <- Sys.time()
  cat(paste("done assessment", "time", difftime(end_time_assessment_complete, start_time_assessment_complete, units="secs"), "s", "\n", sep = " "))
  mu <- apply(rawdata[feature_genes,], 1, mean)
  if (rawdata_null==T) {
    rawdata <- NULL
  }
  assessment <- list(rawdata=rawdata,rowmean=mu,part=part, clustsize=clustsize, features=features, assessed_features=feature_genes,PCA=pctab,max_cl=max_cl,cutoff=cutoff,f1_score=f1_score, Entropy_thresh=Entropy_thresh, Entropy_median=Entropy_median,cluster=cluster,enriched_features=enriched_features, enriched_feature_list=enriched_feature_list, unique_features=unique_features, unique_feature_list=unique_feature_list,second_cluster=second_cluster$dataY,list2ndShared=second_cluster$list_2nd, shared2ndgenes=shared2ndgenes,list_2nd_diff=second_cluster$list_2nd_diff,outliertab=outliertab, outlier_genes=outlier_genes,  diptest=pval_list, nonunimodal=nonunimodaltab, nonunimodal_list=pval_list,bandwidth_list=pval_list2, masstest_list=pval_list3,  batch_entropy=batch_entr )
  if (is.null(assessment_list)) {
    assessment_list <- list()
  }
  if(!is.null(set.name)) { assessment_list[[set.name]] <- assessment}
  else {
    if (!is.null(seuratobject)) {
      if(is.null(givepart)) {
        ident_seu <- names(which(apply(seuratobject@meta.data, 2, function(x){
          identical(as.numeric(table(x)), as.numeric(table(Seurat::Idents(seuratobject))))
        }))[1])
        list_name <- paste0("S", sub(".+\\_", "", ident_seu ))
      }
      else{
        if (class(givepart) =="character") {
          list_name <- paste0("S", sub('\\"+',"", sub(".+\\_","",givepart)))
        }
        else{
          list_name <- sub(".+\\$", "", deparse(substitute(givepart)))
          list_name <- paste0("S", sub(".+\\_", "", list_name ))}
      }
      assessment_list[[list_name]] <- assessment
    }
    if (!is.null(RaceIDobject)) {
      if(is.null(givepart)) {
        assessment_list[["Vdefault"]] <- assessment
      }
      else {
        if (class(givepart) =="character") {
          list_name <- paste0("V", sub('cl.',"", sub("res.","res",givepart)))
        }
        else {
          #list_name <- sub(".+\\$", "", deparse(substitute(givepart)))
          list_name <- sub("\\$.+","", deparse(substitute(givepart)))
          list_name <- sub("cl_", "", list_name)
          list_name <- paste0("V", sub("\\_.+", "", list_name ))
          #list_name <- paste0("V", sub(".+\\_", "", list_name ))
        }
        assessment_list[[list_name]] <- assessment
      }
    }
  }
  return(assessment_list)


}

#' @import ggpubr
#' @import dplyr
#' @import ggplot2
#' @title Plot differences in f1-score or entropy of individual genes
#' @description  This function serves to explore differences in f1-score or entropy of individual genes between different assessed cluster partitions. Genes with large differences between the two assessments are highlighted.
#' @param assessment1 assessment of a cluster partition for which either f1_score or Entropy computation has been performed based on a derived per gene cutoff.
#' @param assessment2 assessment of a second partition derived for the same data.
#' @param toplot plot either “F1” or “Entropy”. Differences of F1 or Entropy calculations of genes between the assessed cluster partitions. Default = "F1".
#' @param signifi logical. if \code{T}, calculate standard deviation of absolute differences between F1/Entropy gene-pairs of two assessments and multiplies this standard deviation by 1.64485361, representing the z-score with significance level based on the normal distribution function. The expected difference of a pair equals 0. Genes with 1.64485361 standard deviations(sds) difference between the two assessments are considered significant and are highlighted. Default = \code{F} and all genes with differences of one standard deviation are highlighted.
#' @param out_doub_dif vector of gene names of user specified gene sets to highlight. Default = \code{NULL} and genes with differences larger than 1 sd or 1.64485361 sds are highlighted.
#' @param xlabs character string for x-axis label indicating assessment 1. Default = \code{NULL} and label is derived from object name of assessment 1.
#' @param ylabs character string for y-axis label indicating assessment 2. Default = \code{NULL} and label is derived from object name of assessment 2.
#' @param maplot logical. If \code{T}, plot the F1 or Entropy per gene difference of assessment 1 and assessment 2 against the mean expression of the gene.
#' @param logmean either “log2”, “log10” or NULL. If “log2” or “log10” performs log transformation of the mean expression per gene. Default = \code{NULL}.
#' @examples
#' binaryclass_scatterplot(assessmentRC$Sres.1, assessmentRC$Sres.5, maplot = T, logmean = "log2", signifi = T)
#' @export
binaryclass_scatterplot <- function(assessment1, assessment2, maplot=F,toplot= "F1",signifi=F,out_doub_dif=NULL,  xlabs=NULL, ylabs=NULL, logmean=NULL){
  if ( ! toplot %in% c("F1", "Entropy")) { stop("set to plot either to F1 or Entropy = T") }
  else if (toplot == "F1") {
    f1_1 <- assessment1$f1_score
    f1_2 <- assessment2$f1_score}
  else {
    f1_1 <- assessment1$Entropy_thresh
    f1_2 <- assessment2$Entropy_thresh
  }
  if (is.null(f1_1) || is.null(f1_2)) { stop(paste("at least one of the assessments'", toplot, "slots are NULL and thus not calculated"))}

  int_F1 <- intersect(names(f1_1), names(f1_2))
  f1_1 <- f1_1[int_F1]
  f1_2 <- f1_2[int_F1]

  require(RColorBrewer)
  set1 <- brewer.pal(4, "Set1")
  sdfco <- sd(abs(f1_2-f1_1))
  if (signifi==T) {
    sdfco = sdfco*1.64485361
  }
  data <- data.frame(f1_1=f1_1, f1_2=f1_2, diff=f1_1-f1_2)
  sig <- rep(3, nrow(data))
  sig[which( f1_1 - (f1_2 + sdfco) > 0) ] <- 1
  sig[which(f1_2 - (f1_1 + sdfco) > 0)] <- 2
  if ( maplot) {
    if (!is.null(logmean)) {
      if( ! logmean %in% c("log2", "log10")) { stop("perform either log2 or log10 transformation of the mean")}
      else if ( logmean == "log10") {
        mu <- log10(assessment1$rowmean[int_F1]) }
      else{
        mu <- log2(assessment1$rowmean[int_F1])
      }
    }
    else{
      mu <- assessment1$rowmean[int_F1]
    }
    data <- data.frame(name=names(f1_1), diff = f1_1-f1_2, mean=mu,  sig=sig)
  }
  else{
    data <- data.frame(name=names(f1_1), f1_1=f1_1, f1_2=f1_2, sig=sig)}
  data$sig <- as.factor(data$sig)
  .lev <- ggpubr:::.levels(data$sig) %>% as.numeric()
  palette = c("#B31B21", "#1465AC", "darkgray")
  palette <- palette[.lev]
  if (!is.null(out_doub_dif)) {
    labs_data <- data[data$sig == 1 | data$sig == 2 ,]
    labs_data <- labs_data[labs_data$name %in% out_doub_dif,]
  }
  else{
    labs_data <- data[data$sig == 1 | data$sig == 2,]
  }
  if( toplot == "F1") {
    new.levels <- c(paste0("1stUp: ", "F1"), paste0("2ndUp: ",
                                                    "F1"), "NS") %>% .[.lev]}
  else{
    new.levels <- c(paste0("1stUp: ", "entropy"), paste0("2ndUp: ",
                                                         "entropy"), "NS") %>% .[.lev]
  }
  data$sig <- factor(data$sig, labels = new.levels)

  if(maplot) {
    p <- ggplot(data, aes(x = mean, y = diff)) + geom_point(aes(color = sig))
  }
  else {
    p <- ggplot(data, aes(x = f1_1, y = f1_2)) + geom_point(aes(color = sig)) }
  p <- p + ggrepel::geom_text_repel(data = labs_data, mapping = aes(label = name),
                                    box.padding = unit(0.35, "lines"), point.padding = unit(0.3,
                                                                                            "lines"), force = 1, fontface = "plain",
                                    size = 2, color = "black", max.overlaps = 200)
  if (toplot == "F1") { main <- "F1_score"}
  else {main <- "Entropy"}
  if (maplot) {
    if (!is.null(logmean)) {
      p <- p + labs(x = paste(logmean, "Mean expression"), y = paste("/\\",toplot, "of" , paste(strsplit(deparse(substitute(assessment1)), "\\$" )[[1]][-1], collapse = "_"), "-", paste(strsplit(deparse(substitute(assessment2)), "\\$" )[[1]][-1], collapse = "_")), title = main, color = "") +
        geom_hline(yintercept = 0, color = "black") + geom_hline(yintercept = 0+sdfco, color = "red", linetype="dashed") +
        geom_hline(yintercept = 0-sdfco, color = "red", linetype="dashed")
    }
    else{
      p <- p + labs(x = "Mean expression", y = paste("/\\",toplot, "of" , paste(strsplit(deparse(substitute(assessment1)), "\\$" )[[1]][-1], collapse = "_"), "-", paste(strsplit(deparse(substitute(assessment1)), "\\$" )[[1]][-1], collapse = "_")), title = main, color = "") +
        geom_hline(yintercept = 0, color = "black") + geom_hline(yintercept = 0+sdfco, color = "red", linetype="dashed") +
        geom_hline(yintercept = 0-sdfco, color = "red", linetype="dashed") }
  }
  else {
    if (is.null(xlabs) && is.null(ylabs)) {
      p <- p + labs(x = paste(strsplit(deparse(substitute(assessment1)), "\\$" )[[1]][-1], collapse = "_"), y = paste(strsplit(deparse(substitute(assessment2)), "\\$" )[[1]][-1], collapse = "_"), title = main, color = "") +
        geom_abline(intercept = 0, slope = 1, color = "black") + geom_abline(intercept = 0+sdfco, slope = 1, color = "red", linetype="dashed") +
        geom_abline(intercept = 0-sdfco, slope = 1, color = "red", linetype="dashed") }
    else{
      p <- p + labs(x = xlabs, y = ylabs, title = main, color = "") +
        geom_abline(intercept = 0, slope = 1, color = "black") + geom_abline(intercept = 0+sdfco, slope = 1, color = "red", linetype="dashed") +
        geom_abline(intercept = 0-sdfco, slope = 1, color = "red", linetype="dashed")

    }
  }
  p <- ggpar(p, palette = palette, ggtheme = theme_classic())
  p
}

#' @title Plot function to determine optimal resolution parameter of community detection algorithms after cluster assessment.
#' @description  After cluster assessment, this function serves to  optimize community detection clustering resolution parameter by highlighting resolution with maximal number of clusters with a user defined threshold and saturation point at which number increase of resolution only linearly decreases average number of detected outlier cells across clusters.
#' @param assessment_list list of assessments for partitions of e.g. increasing resolution parameters of community detection methods.
#' @param cex  numeric, graphical parameter indicating the amount by which the line connecting data points should be scaled. Default = 1.
#' @param f1_thr numeric, threshold used to calculate how many clusters have at least 1 gene with F1-score above this threshold for different cluster partitions assessed. Default = 0.5.
#' @param max_leng numeric, calculation of number of clusters with at least \code{max_leng} genes with minimal F1-score of \code{f1_thr} for the different cluster partitions assessed. Default = 3.
#' @param lcol color used for highlighting the line connecting the data points. Default = “red”.
#' @param resolu logical. If \code{T}. Calculates and highlights saturation point when number of average outlier cells decreases linear with increasing cluster number using two approaches:1) linear model of dependency of average number of outlier cells on the number of different clusters of the different assessed resolutions/cluster partitions. Resolution with the largest negative distance to the fit is highlighted as saturation point. 2) calculates saturation point using an elbow criterion. In addition, highlights the resolution with the maximal number of clusters with 1 gene and \code{max_leng} genes with F1-scores >= \code{f1_thr}. Default = \code{T}.
#' @param sat2 logical. If \code{T}, resolution is fulfilling saturation criterion, if one of the next 3 resolutions also fulfills the saturation criterion.
#' @return plot with 6 graphs, plotting information about cluster partition against number of clusters and number of clusters against average F1-score, average Entropy, average No. of outlier genes across clusters, number of clusters with F1-scores >= \code{f1_thr} and number of clusters with \code{max_leng} genes with F1-score >= \code{f1_thr}.
#'   \item{output_tab}{data.frame with with different resolutions/cluster partitions as rows and the following columns: “No.cluster” = number of assessed clusters, “mean_F1” = mean F1_Score across genes, “mean_Entropy” = mean Entropy across genes, “mean_No.cell_outlg1” = mean number of cells with 1 outlier gene expression across clusters, “No._cluster_F1_>=_f1thr”= number of clusters with at least 1 genes with F1 >= f1thr, “No._cluster_max_leng _genes_w._F1_>=_ f1thr “= number of clusters with at least max_leng genes with F1 >= f1thr, “F1_max_genes” = for every resolution, highest rank genes based on F1-score for clusters if F1>= f1thr.}
#' @examples
#' opti_resolution_plot(assess_seuratRC, f1_thr = 0.5, max_leng = 3, lcol = "red", resolu = T)
#' @export
opti_resolution_plot <- function(assesment_list, cex=1, f1_thr=0.5, max_leng=3, lcol="red", resolu=T, sat2=F) {

  one_putput <- lapply(assesment_list, function(x) {
    output <- c(No.cluster=length(x$cluster),mean_F1=mean(x$f1_score), mean_Entropy=mean(x$Entropy_thresh),mean_No.cell_outlg1=mean(x$outliertab[,1]))
  })
  output_tab <- data.frame(do.call(rbind, one_putput))
  F1_output <- lapply(assesment_list, function(x){
    F1_max_clus <- list()
    max_cl <- x$max_cl
    f1_score <- x$f1_score
    cluster <- x$cluster
    for ( i in 1:length(cluster) ){
      F1_max_clus[[i]] <- f1_score[which(as.numeric(max_cl) == cluster[i])]
    }
    names(F1_max_clus) <- paste0("cl",cluster)
    F1_max_clus <- lapply(F1_max_clus, function(x) {
      x[order(x, decreasing = T)]
    })
    return(F1_max_clus)
  })
  ### for every clustering the max F1 score per cluster
  F1_max <- lapply(F1_output, function(x) {
    y <- lapply(x, function(z){
      z[1]
    })
    y <- unlist(y)
    y[is.na(y)] <- 0
    y
  })

  F1_max_tab <- lapply(F1_max, function(x, f1_thr){
    x[x >= f1_thr]
  },f1_thr=f1_thr)
  F1_max_tab <- unlist(lapply(F1_max_tab, function(x){
    y <- names(x)
    y <- paste(y, collapse = ",")
  }))


  F1_max_genes <- strsplit(F1_max_tab, split = ",")
  F1_max_genes <- unlist(lapply(F1_max_genes, function(x) {
    x <- gsub("cl\\d+\\.","",x)
    x <- paste(x, collapse = ",")
  }))

  F1_max_thr_abs <- unlist(lapply(F1_max, function(x, f1_thr){
    sum(x >= f1_thr)
  },f1_thr=f1_thr))

  output_tab <- cbind(output_tab, F1_max_thr_abs)
  ### mean of top length()
  F1_max_cl <- lapply(F1_output, function(x, lengi){
    y <- lapply(x, function(z,lengi){
      z[lengi]
    }, lengi=lengi)
    y <- unlist(y)
    y[is.na(y)] <- 0
    y
  }, lengi=max_leng)

  F1_max_cl_thr_abs <- unlist(lapply(F1_max_cl, function(x, f1_thr){
    sum(x >= f1_thr)
  },f1_thr=f1_thr))

  output_tab <- cbind(output_tab, F1max_g_absthr=F1_max_cl_thr_abs)
  colnames(output_tab)[5] <- paste("No._cluster_F1_>=", f1_thr, sep = "_")
  colnames(output_tab)[6] <- paste("No._cluster", max_leng,"genes_w._F1_>=", f1_thr, sep = "_")

  maxi <- apply(output_tab, 2, max)
  mini <- apply(output_tab, 2, min)
  mini <- sapply(mini, function(x) {x-0.1*x})
  maxi <- sapply(maxi, function(x) {x+0.1*x})
  clus_range <- maxi[1]
  par(mfrow=c(2,3))
  for ( n in 1:ncol(output_tab)) {
    if (n == 1) {
      plot(output_tab[,n], type = "l", ylim = c(0,clus_range),ylab = colnames(output_tab)[n], xaxt="n", xlab="", col=lcol, lwd=2)
      axis(1, at=1:length(assesment_list), labels = rownames(output_tab), las=2, cex.axis=cex )
    }
    else{
      maxis <- maxi[n]
      minis <- mini[n]
      plot(x=output_tab[,1], type = "l", y=output_tab[,n],xlim = c(0,clus_range), ylim = c(minis, maxis),ylab = colnames(output_tab)[n], xlab="No. cluster", col=lcol, lwd=2)
      if(resolu==T){
        if ( n == 4) {
          outlg <- as.numeric(output_tab[,n])
          fit <- lm(outlg~as.numeric(output_tab[,1]))
          resi <- which(residuals(fit) == min(residuals(fit)))
          resi <- as.numeric(resi)
          points(output_tab[,1][resi], output_tab[,n][resi],cex = 1.5,col = "green", pch=20)
          text(output_tab[,1][resi], output_tab[,n][resi], labels = sub("S", "", rownames(output_tab)[resi]), cex = 1)
          ### commented out for simplicity
          outlier <- outlg/sum(outlg)
          y <- outlier[-length(outlier)] - outlier[-1]
          mm <- numeric(length(y))
          nn <- numeric(length(y))
          for (k in 1:length(y)) {
            mm[k] <- mean(y[k:length(y)])
            nn[k] <- sqrt(var(y[k:length(y)]))
          }
          ind <- which(y - (mm + nn) < 0)
          if (sat2==T) {
            for (p in ind) {
              if (sum(p:(p + 3) %in% ind) >= 2) {
                ind <- p
                break
              }
            }
            points(output_tab[,1][ind], output_tab[,n][ind],cex = 1.5,col = "blue", pch=20)
            text(output_tab[,1][ind], output_tab[,n][ind], labels = sub("S", "", rownames(output_tab)[ind]), cex = 1)
          }
          else{
            ind <- as.numeric(ind[1])
            points(output_tab[,1][ind], output_tab[,n][ind],cex = 1.5,col = "blue", pch=20)
            text(output_tab[,1][ind], output_tab[,n][ind], labels = sub("S", "", rownames(output_tab)[ind]), cex = 1)
          }
        }
        if ( n %in% c(5,6)) {
          f1f1max <- as.numeric(which(output_tab[,n] == max(output_tab[,n])))
          points(output_tab[,1][f1f1max], output_tab[,n][f1f1max],cex = 1.5,col = "green", pch=20)
          text(output_tab[,1][f1f1max], output_tab[,n][f1f1max], labels = sub("S", "", rownames(output_tab)[f1f1max]), cex = 1)
        }

      }

    }
  }
  return(cbind(output_tab,F1_max_tab=F1_max_tab, F1_max_genes=F1_max_genes ))
}

#' @title Plot function to determine optimal pre-processing method.
#' @description  After cluster assessment, this function serves to identify optimal pre-processing method, independent of the number of clusters, plotting the following criteria: average F1-Score, average Entropy, average number of outlier genes across cluster and average number of enriched features across clusters.
#' @param assessment_list2 list of assessment lists exhibiting different assessments, for example: different normalization methods and for each increasing resolution parameters.
#' @param cex = numeric, graphical parameter indicating the amount by which the line connecting the data points should be scaled. Default = 1
#' @param f1_thr numeric, threshold used to calculate how many clusters have at least 1 gene with F1-score above this threshold for different cluster partitions assessed. Default = 0.5.
#' @param max_leng numeric, calculation of number of clusters with at least \code{max_leng} genes with minimal F1-score of \code{f1_thr} for the different cluster partitions assessed. Default = 3.
#' @param lcol = vector of colors used for highlighting slots of list of assessments, for each list of lists of assessment one color: e.g. list of resolution optimizations for e.g. normalization A, and another color for list of resolution optimization for e.g. normalization B.
#' @param map = logical. If \code{T}, then legend is shown. Default = \code{T}.
#' @param leg = logical. If \code{T}, then the legend is shown. Default = \code{T}.
#' @return plot with 6 graphs, plotting information about cluster partition against number of clusters and number of assessed genes, as well as plotting number of clusters against average F1-score, average Entropy, average number of enriched features assessed and average No. of outlier genes across clusters.
#'   \item{output_list}{with with data.frame for every list of assessments, with with different resolutions/cluster partitions as rows and the following columns: “No.cluster” = number of assessed clusters, “mean_F1” = mean F1_Score across genes, “mean_Entropy” = mean Entropy across genes, “mean_No.cell_outlg1” = mean number of cells with 1 outlier gene expression across clusters, "mean_enriched_features" = mean numer of enriched feautres across clusters of assessed features and "assessed_features" = Number of assessed features.}
#' @examples
#' opti_preprocess_plot(assessment_list, lcol = c("red", "blue", "green"))
#' @export
opti_preprocess_plot <- function(assesment_list2, cex=1, lcol=c("red"), map=T, leg=T) {
  if (class(assesment_list2[[1]][[1]]) != "list") {
    naming <- deparse(substitute(assesment_list2))
    assesment_list2 <- list(assesment_list2)
    names(assesment_list2) <- naming
  }
  output_list <- list()
  for (i in 1:length(assesment_list2)) {
    assesment_list <- assesment_list2[[i]]

    one_putput <- lapply(assesment_list, function(x) {
      output <- c(No.cluster=length(x$cluster),mean_F1=mean(x$f1_score), mean_Entropy=mean(x$Entropy_thresh),assessed_features=length(x$assessed_features),mean_enriched_features=mean(x$enriched_features),mean_No.cell_outlg1=mean(x$outliertab[,1]))
    })
    output_tab <- data.frame(do.call(rbind, one_putput))
    output_list[[i]] <- output_tab
  }
  names(output_list) <- names(assesment_list2)
  maxi <- lapply(output_list, function(x) {
    apply(x, 2, max)
  })
  maxi <- do.call(rbind, maxi)
  maxi <- apply(maxi, 2, max)
  mini <- lapply(output_list, function(x) {
    apply(x, 2, min)
  })
  mini <- do.call(rbind, mini)
  mini <- apply(mini, 2, min)
  mini <- sapply(mini, function(x) {x-0.1*x})
  maxi <- sapply(maxi, function(x) {x+0.1*x})
  clus_range <- maxi[1]
  par(mfrow=c(2,3))
  for ( n in 1:ncol(output_list[[1]])) {
    if (n == 1 ) {
      for ( i in 1:length(output_list)) {
        if ( i == 1) {
          if ( map == T) {
            plot(x = 1:nrow(output_list[[1]]),y=output_list[[i]][,n], type = "l", ylim = c(0,clus_range),ylab = colnames(output_list[[i]])[n], xaxt="n", xlab="", col=lcol[i], lwd=2)
            axis(1, at=1:length(assesment_list2[[i]]), labels = rownames(output_list[[i]]), las=2, cex.axis=cex ) }
          else {
            plot(x = 1:nrow(output_list[[1]]),y=output_list[[i]][,n], xlab=NA, ylab=NA, axes = FALSE, cex=0, pch=20)
          }
          if (length(output_list) > 1) {
            if ( leg == T) {
              legend("topright", names(assesment_list2), fill = lcol, cex=1, bty="n")}
          }
        }
        else{
          if ( map == T) {
            lines(x = 1:nrow(output_list[[1]]), y=output_list[[i]][,n],  col=lcol[i], lwd=2)
          }
        }
      }

    }
    else if ( n == 4) {
      for ( i in 1:length(output_list)) {
        if ( i == 1) {
          maxis <- maxi[n]
          minis <- mini[n]
          if ( map == T) {
            plot(x=1:nrow(output_list[[1]]),y=output_list[[i]][,n], type = "l",ylab = colnames(output_list[[i]])[n], xaxt="n", xlab="", col=lcol[i], lwd=2, ylim = c(minis, maxis))
            axis(1, at=1:length(assesment_list2[[i]]), labels = rownames(output_list[[i]]), las=2, cex.axis=cex )
          }
          else {
            plot(x=1:nrow(output_list[[1]]),y=output_list[[i]][,n], xlab=NA, ylab=NA, axes = FALSE, cex=0, pch=20)
          }
          if (length(output_list) > 1) {
            if ( leg == T) {
              legend("topright", names(assesment_list2), fill = lcol, cex=1, bty="n")}
          }
        }
        else{
          if ( map == T) {
            lines(x = 1:nrow(output_list[[1]]), y=output_list[[i]][,n],  col=lcol[i], lwd=2)}
        }
      }

    }
    else{
      for ( i in 1:length(output_list)) {
        if ( i == 1) {
          maxis <- maxi[n]
          minis <- mini[n]
          if ( map == T) {
            plot(x=output_list[[i]][,1], type = "l", y=output_list[[i]][,n],xlim = c(0,clus_range), ylim = c(minis, maxis),ylab = colnames(output_list[[i]])[n], xlab="No. cluster", col=lcol[i], lwd=2)
          }
          else {
            plot(x=output_list[[i]][,1], y=output_list[[i]][,n],xlab=NA, ylab=NA, axes = FALSE, cex=0, pch=20)
          }
          if (length(output_list) > 1) {
            if ( leg == T) {
              legend("topright", names(assesment_list2), fill = lcol, cex=1, bty="n")}
          }

        }
        else{
          if ( map == T) {
            lines(x=output_list[[i]][,1], y=output_list[[i]][,n], col=lcol[i], lwd=2)   }
        }
      }
    }
  }
  return(output_list)
}


