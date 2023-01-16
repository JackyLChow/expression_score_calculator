library(readr)
library(stringr)
library(preprocessCore)

# follow methods described in Ayers, JCI 2017;Damotte, JTM 2019
## quantile normalize gene expression, then average the normalized expression for genes in the signature

###---data import---------------------------------------------------------------
# parameters file
## complete counts normalization and batch corrections before running
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
input_counts <- readRDS(parameters["counts_file", "value"]) # load normalized and transformed counts
quant_counts <- log2(normalize.quantiles(counts) + 1) # normalize; quantile normalization should be robust to counts type
row.names(quant_counts) <- rownames(input_counts)
colnames(quant_counts) <- colnames(counts)

metadata <- data.frame(read.csv(parameters["metadata_file", "value"], row.names = 1))
rownames(metadata) <- str_replace_all(rownames(metadata), "-", ".")
metadata <- metadata[colnames(quant_counts), ]

###---calculate scores----------------------------------------------------------
scores <- data.frame(read.csv(parameters["signatures_file", "value"])) # load scores file
row.names(scores) <- scores[, parameters["signature_name_column", "value"]]

missing_score_genes <- c() # initialize vector for signature genes not in expression data
sample_scores <- data.frame() # score calculator outputs

for(i in rownames(scores)){
  #i <- rownames(scores)[1]
  # extract and curate gene signatures
  score_ <- scores[i, parameters["signature_genes_column", "value"]] # extract genes in score
  score_ <- unlist(str_split(score_, "\\, |\\,|; |\\;")) # split into vector
  missing_score_genes <- c(missing_score_genes, score_[!score_ %in% rownames(quant_counts)]) # append missing genes
  score_ <- score_[score_ %in% rownames(quant_counts)] # subset to genes present in expression data
  # calculate score
  sample_score_ <- data.frame(score_name = i,
                              sample = colnames(quant_counts),
                              score = colMeans(quant_counts[score_, ]),
                              metadata,
                              row.names = NULL)
  
  sample_scores <- rbind(sample_scores, sample_score_)
}

saveRDS(sample_scores, paste0(parameters["output_folder", "value"], "signature_scores.rds"))


