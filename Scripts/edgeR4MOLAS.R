#!/usr/bin/env Rscript
library(optparse, quietly = TRUE)
library(edgeR, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(ggrepel, quietly = TRUE)

# Setting arguments ----
option_list = list(
  make_option(c("-a", "--analysis_type"), type="character", default = NULL,
              help="select 'normalize' or 'statistical_test'", metavar = "character"),
  make_option(c("-r", "--feature_count_table"), type="character", default="InputData/FeatureCount.csv", 
              help="dataset file.csv path [default= %default]", metavar="character"),
  make_option(c("-m", "--metadata"), type = "character", default="InputData/metadata.csv",
              help="metadata file.csv path [default= %default]", metavar = "character"),
  make_option(c("-o", "--output_folder"), type="character", default="OutputData", 
              help="output folder name [default= %default]", metavar="character"),
  make_option(c("-c", "--comparison"), type = "character", default = NULL,
              help="the name of compared groups, which should be 'GroupA,GroupB'", metavar = "character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser, args = c("-a", "statistical_test", "-c", "DMSO,Fostamatinib-treated"))

# Check arguments
if(is.null(opt$analysis_type)){
  stop("Choose analysis type: 'normalize' or 'statistical_test'", call.=FALSE)
} else if (is.null(opt$output_folder)) {
  stop("Set output folder path", call.=FALSE)
} 
if(opt$analysis_type == "normalize"){
  if(is.null(opt$feature_count_table) | is.null(opt$metadata)){
    stop("'feature_count_table' and 'metadata' is required", call.=FALSE)
  }
}
if (opt$analysis_type == "statistical_test"){
  if(is.null(opt$feature_count_table ) | is.null(opt$metadata)){
    stop("'feature_count_table' and 'metadata' is required", call.=FALSE)
  } else if (is.null(opt$comparison)){
    print_help(opt_parser)
    stop("At least one comparison should be supplied", call.=FALSE)
  } else {
    opt$comparison_processed <- unlist(strsplit(opt$comparison, ","))
    opt$statistical_method <- "GLM_LRT" 
  }
}

# Normalization ----
# Load data
FeatureCount <- read.csv(opt$feature_count_table, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
metadata <- read.csv(opt$metadata, stringsAsFactors = FALSE, check.names = FALSE)

# Check data format and sample completeness
CheckSampleCompleteness <- identical(colnames(FeatureCount), metadata$SampleID)
if(!CheckSampleCompleteness){stop("Samples in Feature count table are not identical to samples in metadata table")}
CheckFeatureCountTableFormat <- any(c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique") %in% rownames(FeatureCount))
if(CheckFeatureCountTableFormat){stop("Please remove special counters from htseq-count output")}

# Function
edgeR.Filter.Normalize <- function(Express, Group, NormMethod){
  # Argument:
  #   Express: edgeR.input(count table, not fpkm or other transformed number), data.frame
  #   Group: sample classification, factor
  #   NormMethod: normalizing method, character
  
  y <- DGEList(counts = Express, group = Group)
  
  # Filtering
  y_keep <- filterByExpr(y)
  print(paste0("True expressing gene number: ", sum(y_keep == TRUE)))
  y <- y[y_keep, , keep.lib.size=FALSE]
  
  Group_name <- paste(levels(Group), table(Group), sep = "*", collapse = ",")
  
  # Normalization
  y <- calcNormFactors(y, method = NormMethod)
  write.csv(cpm(y), paste0("OutputData/FeatureCount_filtered_normalized(", NormMethod, ")(", Group_name, ").csv"), quote = FALSE)
  
  return(y)
}
glmLRTestPipe <- function(Target, Coefficient=NULL, Contrast=NULL, P.value, Log.fold.change=NULL){
  # Argument:
  #   Target: DGEList object, already calculated the dispersion
  #   Coefficient/Contrast: experiment design and chosed comparison/like c(-1,1,0) = B compare to A, choose one of them.
  #   P.value: p-value threshold, numeric
  #   Log.fold.change: the numeric scaler of the log2-fold change threshold above which differential expression is to be considered.
  
  # Transform args.Contrast format
  Contrast_null <- rep(0, length(levels(Target$samples$group)))
  Contrast_null[levels(Target$samples$group) == Contrast[1]] <- 1
  Contrast_null[levels(Target$samples$group) == Contrast[2]] <- -1
  Contrast <- Contrast_null
  
  # Model Test
  fit <- glmFit(Target)
  if (is.null(Log.fold.change)) {
    TestResult <- glmLRT(fit, coef = Coefficient, contrast = Contrast)
  } else if (is.numeric(Log.fold.change)) {
    TestResult <- glmTreat(fit, coef = Coefficient, contrast = Contrast, lfc = Log.fold.change)
  }
  
  print(topTags(TestResult))
  TestResult.is.de <- decideTests(TestResult, p.value = P.value)
  print(summary(TestResult.is.de))
  FDR <- topTags(TestResult, n = nrow(TestResult), sort.by = "none")$table$FDR
  Output_file <- cbind(TestResult$table, FDR, TestResult.is.de)
  
  # File name
  SelectedContrast <- if(!is.null(Coefficient)){TestResult$comparison} else{paste0(Contrast, "*", colnames(design), collapse = "_")}
  file_name <- if(is.null(Log.fold.change)) {
    paste(opt$output_folder, "/glmLRTest(",SelectedContrast,")(FDR=",P.value,").csv", sep = "")
  } else if (is.numeric(Log.fold.change)) {
    paste(opt$output_folder, "/glmLRTest(",SelectedContrast,")(LogFoldChange_threshold=",Log.fold.change,")(FDR=",P.value,").csv", sep = "")
  }
  
  write.csv(Output_file, file = file_name)
  
  return(TestResult)
}
exactTestPipe <- function(Target, Pair, P.value){
  # Argument:
  #   Target: DGEList object, already calculated the dispersion
  #   Pair: Specified comparison, e.g. c("A","B") is compare A and B
  #   P.value: p-value threshold, numeric
  
  # Transform args.Pair format
  Pair <- c(Pair[1], Pair[2])
  
  et <- exactTest(object = Target, pair = Pair)
  print(topTags(et))
  et.is.de <- decideTests(et, p.value = P.value)
  print(summary(et.is.de))
  Output_file <- cbind(et$table, et.is.de)
  
  write.csv(Output_file, paste0(opt$output_folder, "/exactTest(",et$comparison[2],"v.s.",et$comparison[1],")(FDR=",P.value,").csv"))
  
  return(et)  
}
Export_MDS_plot <- function(metadata, filename, width, height){
  
  require(ggplot2, quietly = TRUE)
  require(ggrepel, quietly = TRUE)
  
  MDS_data <- plotMDS(edgeR_preprocessed, top = 500, plot = FALSE)
  MDS_plot_data <- data.frame(x = MDS_data[['x']], y = MDS_data[['y']], label = metadata[,1], group = metadata[,2])
  
  p <- ggplot(MDS_plot_data, aes(x = x, y = y, color = group)) + 
    geom_point() + 
    geom_text_repel(aes(label = label), show.legend = FALSE) + 
    theme_bw() + 
    xlab("Leading logFC dim 1") + ylab("Leading logFC dim 2")
  
  ggsave(filename = filename, plot = p, device = "png", width = width, height = height, dpi = 300, units = "cm")
  
}
Export_BCV_plot <- function(filename, width, height){
  png(filename = filename, width = width, height = height, units = "cm", res = 300)
  plotBCV(edgeR_preprocessed_disp)
  dev.off()
}

# Gene filtration and Sample normalization
edgeR_preprocessed <- edgeR.Filter.Normalize(Express = FeatureCount,
                                             Group = factor(metadata$Group, levels = unique(metadata$Group)),
                                             NormMethod = "TMM")

# Whether continue to statistical test
if (opt$analysis_type == "normalize"){quit(save = "no")}

# Statistical test ----
# Construct MDS plot
Export_MDS_plot(metadata = metadata, filename = paste0(opt$output_folder, "/MDS_plot.png"), width = 15, height = 12)

# Set experiment design matrix
design <- model.matrix(~0+group, data = edgeR_preprocessed$samples)

# Estimate dispersion
edgeR_preprocessed_disp <- estimateDisp(edgeR_preprocessed, design, robust = TRUE)
invisible(capture.output(Export_BCV_plot(filename = paste0(opt$output_folder, "/biological_coefficient_of_variation.png"), width = 20, height = 15)))

# statistical calculation
edgeR_glmLRTest <- glmLRTestPipe(Target = edgeR_preprocessed_disp, Coefficient = NULL, Contrast = opt$comparison_processed, P.value = 0.05, Log.fold.change = NULL)

