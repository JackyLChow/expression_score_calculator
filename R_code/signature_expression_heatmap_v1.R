library(readr)
library(stringr)
library(spatstat)
library(ComplexHeatmap)
library(circlize)
library(pals)

###---data import---------------------------------------------------------------
# parameters file
## complete counts normalization and batch corrections before running
parameters <- data.frame(read_csv("~/Documents/BFX_proj/expression_score_calculator/_input/TCGA_PanCancer.csv"))
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

###---load data-----------------------------------------------------------------
input_counts <- readRDS(parameters["counts_file", "value"]) # load normalized/transformed counts
metadata <- readRDS(parameters["metadata_file", "value"])

#foo <- sample(1:10065, 200)
#input_counts <- input_counts[, foo]
#metadata <- metadata[foo, ]

scores <- data.frame(read.csv(parameters["signatures_file", "value"])) # load scores file
row.names(scores) <- scores[, parameters["signature_name_column", "value"]]

###---scores heatmap------------------------------------------------------------
score_annotation <- data.frame(row.names = unique(unlist(str_split(scores$Genes, "\\, |\\,|; |\\;")))) # extract genes from all scores
score_annotation <- score_annotation[rownames(score_annotation) %in% rownames(input_counts), ] # filter out genes not in counts
for(i in unique(scores[, parameters["signature_name_column", "value"]])){
  score_annotation[, i] <- ifelse(rownames(score_annotation) %in%
                                             unlist(str_split(scores[i, parameters["signature_genes_column", "value"]], "\\, |\\,|; |\\;")),
                                  1, 0)
}

heatmap_counts <- input_counts[rownames(score_annotation), ]
heatmap_counts <- t(scale(t(heatmap_counts)))

top_annotation <- HeatmapAnnotation(comparison = metadata[, parameters["metadata_comparison_column", "value"]],
                                    col = list(comparison = setNames(
                                      c(alphabet(26), alphabet2(26))[1:length(unique(metadata[, parameters["metadata_comparison_column", "value"]]))],
                                      c(unique(metadata[, parameters["metadata_comparison_column", "value"]])))))

set.seed(415); left_annotation <- Heatmap(as.matrix(score_annotation),
                                          col = colorRamp2(c(0, 3), c("white", "black")),
                                          border = "grey20",
                                          cluster_rows = F,
                                          cluster_columns = T,
                                          show_row_names = F,
                                          show_row_dend = F,
                                          show_heatmap_legend = F,
                                          width = unit(0.4 * ncol(score_annotation), "cm"),
                                          height = unit(20, "cm"))

set.seed(415); signature_heatmap <- Heatmap(heatmap_counts,
                                            col = colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow")),
                                            top_annotation = top_annotation,
                                            show_row_names = F,
                                            show_column_names = F,
                                            use_raster = F,
                                            name = "z-score",
                                            height = unit(20, "cm"))

jpeg(paste0(parameters["output_folder", "value"], "signatures_heatmap.jpg"), width = 1250, height = 1000)
draw(left_annotation + signature_heatmap, main_heatmap = "z-score")
dev.off()
