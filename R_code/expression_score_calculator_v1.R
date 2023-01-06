library(readr)
library(preprocessCore)

# follow methods from Ayers, JCI 2017;Damotte, JTM 2019
## quantile normalize gene expression, then average the normalized expression for genes in the signature

###---data import---------------------------------------------------------------
# parameters file
parameters <- data.frame(read_csv("~/Documents/BFX_proj/expression_score_calculator/_input/Chow_PNAS_2020.csv"))
rownames(parameters) <- parameters[, "description"]

# checkpoints
checkpoint_file <- parameters["checkpoint_file", "value"]
if("checkpoint.csv" %in% list.files(parameters["output_folder", "value"])){
  checkpoint <- data.frame(read_csv(checkpoint_file))
} else {
  checkpoint <- data.frame(matrix(ncol = 2, nrow = 0))
  names(checkpoint) <- c("file", "date_created")
  checkpoint[1, ] <- c("checkpoint.csv", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  write.csv(checkpoint, checkpoint_file, row.names = F)
}

###---quantile normalization----------------------------------------------------
input_counts <- data.frame(read.csv(parameters["counts_file", "value"])) # load counts table
counts <- as.matrix(input_counts[, 2:ncol(input_counts)]) # reformat to matrix
quant_counts <- log2(normalize.quantiles(counts) + 1) # normalize
row.names(quant_counts) <- input_counts$gene
colnames(quant_counts) <- colnames(counts)

metadata <- data.frame(read.csv(parameters["metadata_file", "value"], row.names = 1))
rownames(metadata) <- str_replace_all(rownames(metadata), "-", ".")
metadata <- metadata[colnames(quant_counts), ]

###---evaluate scores-----------------------------------------------------------
score_ <- c("IFNG", "PRF1", "GZMB")
score_ <- score_[score_ %in% rownames(quant_counts)]
sample_score_ <- colMeans(quant_counts[score_, ])

score_plot <- metadata
score_plot[, "score"] <- sample_score_

ggplot(score_plot, aes(treatment, score)) +
  geom_boxplot()

