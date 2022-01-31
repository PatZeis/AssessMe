#!/package/R-4.0.3/bin/Rscript --vanilla
library(AssessME)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2 ) {
    stop("please provide at least 2 paths to nfold cross validations of assessments ")
} 

leng_args <- length(args)
combined_list <- list()
for (i in 1:leng_args) {
  files <- list.files(path = as.character(args[i]), pattern = ".Rda", full.names=T)
  for ( n in 1:length(files)) {
    if ( n == 1) {
      accuracy <- unlist(readRDS(files[n]))
    }
    else{
      accuracy2 <- unlist(readRDS(files[n]))
      accuracy <- rbind(accuracy, accuracy2)
    }
  }
  combined_list[[i]] <- accuracy
  file_nam <- strsplit(as.character(args[i]), split = "/")
  file_nam <- file_nam[[1]][length(file_nam[[1]])]
  names(combined_list)[i] <- file_nam
}

pdf(file.path(".", "accuracy.pdf"))
accuracy_plot_hpc(combined_list)
dev.off()
