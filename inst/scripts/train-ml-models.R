# Train the five models shipped in inst/extdata of huGSVcalleR
# This script is a record of that training. It is not called by the package
# It needs SNPSimulatoR (https://github.com/TalvinStrannik/SNPSimulatoR)
# Use that package to simulate another table and
# train a replacement model
#
# GenerateSNPSeqGamma draws read counts from gamma distributions
# The distribution parameters were fitted to NA12878 pile-ups and to
# HG001 variant calls on human chromosome 22 (GRCh37)
# The function writes <seqname>.csv and SNPs_<seqname>.csv

library(SNPSimulatoR)

set.seed(144)
seqname <- "Test"
total_sites <- 20000
modelled <- GenerateSNPSeqGamma(
    seqname = seqname,
    nucleotide_quantity = total_sites,
    snp_quantity = total_sites / 2
)
target_file <- paste0("SNPs_", seqname, ".csv")
train_n <- floor(0.8 * nrow(modelled))
train_rows <- sample(seq_len(nrow(modelled)), size = train_n)
train_data <- modelled[train_rows, ]
methods <- c("rpart", "SVM", "log_reg", "rand_forest", "xgboost")
MakeMLModels(
    data = train_data,
    TargetData = target_file,
    models = methods,
    features = c("x1", "x2", "x3", "x4")
)

# MakeMLModels writes xgboost.rds
# huGSVcalleR loads inst/extdata/xgboost.json
xgb_model <- readRDS("xgboost.rds")
xgboost::xgb.save(xgb_model, "xgboost.json")
