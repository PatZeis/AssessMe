#!/package/R-4.0.3/bin/Rscript --vanilla
library(AssessME)
library(optparse)

option_list <- list(
  make_option(c("-a", "--assessment"), type = "character", default = NULL, help = "name of assessment", metavar = "character"),
  make_option(c("-j", "--obj"), type= "character", default = NULL, help = "name of object", metavar = "character"),
  make_option(c("-c", "--crossvali"), type = "integer", default = 50, help = "how many crossvalidation to perform", metavar = "integer"),
  make_option(c("-i", "--assessindex"), type = "integer", default = 1, help = "index of assessment if assessment is list of assessmens", metavar = "integer"),
  make_option(c("-r", "--rawdata"), type = "character", default=NULL, help = "set countdata", metavar = "character"),
  make_option(c("-p", "--sampling"), type = "integer", default = 70, help = "set fraction you want to sample of each cluster for cross validation", metavar = "integer"),
  make_option(c("-t", "--iterator"), type = "integer", default = 1, help = "loop iterator", metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!is.null(opt$assessment)) {
  assessment_input = readRDS(opt$assessment) }
if (!is.null(opt$obj)) {
  giveobject = readRDS(opt$obj)}
csvali <- opt$crossvali
assess_index = opt$assessindex
if (!is.null(opt$rawdata)) {
  raw <- readRDS(opt$rawdata)}
sample2 <- opt$sampling
iterat <- opt$iterator

tes_generate <- generate_input_accuracy_hpc(assessment=assessment_input, object=giveobject, crossvali=csvali, assessment_no = assess_index, rawdata = raw, tosample=sample2 )
give_accuracy <- accuracy_hpc(generate = tes_generate, crossvali = iterat, ntree = 100)
saveRDS(give_accuracy, paste("accuracy_", iterat, "_assessment", ".Rda", sep = ""))
