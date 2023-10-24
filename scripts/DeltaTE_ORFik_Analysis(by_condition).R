# Argument parsing:
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 4) {
  stop("Four arguments are required: 
  path to rna counts file, 
  path to ribo counts file,
  path to sample info file, 
  and if you want to do a batch correction [T or F]")
}


rna_counts_file <- args[1]
ribo_counts_file <- args[2]
sample_info_file <- args[3]
batch_b <- args[4]

#rna_counts_file <- "/home/gionmattia/Desktop/CANADATA_NEXTFLOW/merged_input.tsv"
#ribo_counts_file <- "/home/gionmattia/Desktop/CANADATA_NEXTFLOW/merged_poly.tsv"
#sample_info_file <- "/home/gionmattia/Desktop/CANADATA_NEXTFLOW/sample_info.txt"
#batch_b <- F


# Load the libraries

library(ORFik)
library(data.table)
library(DESeq2)

# Load the functions missing from ORFik 1.18

DTEG_model_results <- function(ddsMat_rna, ddsMat_ribo, ddsMat_te,
                               target.contrast, pairs, p.value = 0.05,
                               complex.categories) {
  # Do result analysis: per contrast selected
  dt.between <- data.table()
  for(i in pairs) {
    name <- paste("Comparison:", i[1], "vs", i[2])
    message(name)
    # Results
    current.contrast <- c(target.contrast, i[1], i[2])
    res_te <- results(ddsMat_te, contrast = current.contrast)
    
    res_ribo <- results(ddsMat_ribo, contrast = current.contrast)
    suppressMessages(res_ribo <- lfcShrink(ddsMat_ribo, contrast=current.contrast,
                                           res=res_ribo, type = "normal"))
    
    res_rna <- results(ddsMat_rna, contrast = current.contrast)
    suppressMessages(res_rna <- lfcShrink(ddsMat_rna, contrast = current.contrast,
                                          res = res_rna, type = "normal"))
    
    ## The differential regulation groupings (padj is padjusted)
    both <- which(res_te$padj < p.value & res_ribo$padj < p.value & res_rna$padj < p.value)
    
    ## The 4 classes of genes
    # Forwarded are non significant in TE, diagonal line
    forwarded <- rownames(res_te)[which(res_te$padj > p.value & res_ribo$padj < p.value & res_rna$padj < p.value)]
    # These two are the X and Y axis
    exclusive.translation <- rownames(res_te)[which(res_te$padj < p.value & res_ribo$padj < p.value & res_rna$padj > p.value)]
    exclusive.expression <- rownames(res_te)[which(res_te$padj > p.value & res_ribo$padj > p.value & res_rna$padj < p.value)]
    ## These are the remaining groups
    # Also called mRNA abundance
    intensified <- rownames(res_te)[both[which(res_te[both, 2]*res_rna[both, 2] > 0)]]
    # Also called inverse mRNA abundance
    inverse <- rownames(res_te)[both[which(res_te[both, 2]*res_rna[both, 2] < 0)]]
    # Stable protein output
    buffered <- c(rownames(res_te)[which(res_te$padj < p.value & res_ribo$padj > p.value & res_rna$padj < p.value)])
    
    n <- rownames(res_te)
    Regulation <- rep("No change", nrow(res_te))
    Regulation[n %in% forwarded] <- "Forwarded" # Old Buffering
    Regulation[n %in% buffered] <- "Buffering"
    Regulation[n %in% inverse] <- "Inverse"
    Regulation[n %in% exclusive.translation] <- "Translation"
    Regulation[n %in% exclusive.expression] <- "Expression"
    Regulation[n %in% intensified] <- "mRNA abundance"
    
    if (!complex.categories) {
      Regulation[n %in% exclusive.expression] <- "Buffering"
      Regulation[n %in% inverse] <- "Buffering"
      Regulation[n %in% forwarded] <- "Buffering"
    }
    print(table(Regulation))
    
    
    dt.between <-
      rbindlist(list(dt.between,
                     data.table(contrast = name,
                                Regulation = Regulation,
                                id = rownames(ddsMat_rna),
                                rna = res_rna$log2FoldChange,
                                rfp = res_ribo$log2FoldChange,
                                te = res_te$log2FoldChange,
                                rna.padj = res_rna$padj,
                                rfp.padj = res_ribo$padj,
                                te.padj = res_te$padj,
                                rna.meanCounts = res_rna$baseMean,
                                rfp.meanCounts = res_ribo$baseMean
                     )))
  }
  regulation_levels <-  c("No change", "Translation", "Buffering",
                          "mRNA abundance", "Expression", "Forwarded",
                          "Inverse")
  dt.between[, Regulation :=
               factor(Regulation, levels = regulation_levels, ordered = TRUE)]
  return(dt.between)
}


DEG_DESeq <- function(counts, main.design, message = "Creating DESeq model:") {
  message(message)
  message("----------------------")
  ddsMat <- DESeqDataSet(se = counts, design = main.design)
  ddsMat <- DESeq(ddsMat)
  message("----------------------")
  return(ddsMat)
}

# Load the files

ribo_counts = setDT(read.delim(ribo_counts_file))
colnames(ribo_counts) <- gsub("\\.", "-", colnames(ribo_counts))
rna_counts = setDT(read.delim(rna_counts_file))
colnames(rna_counts) <- gsub("\\.", "-", colnames(rna_counts))
sample_info = setDT(read.delim(sample_info_file))

# Good check for your data (use it)
stopifnot(all(ribo_counts$Name==rna_counts$Name))

Gene_names = ribo_counts$Name

ribo_counts$Name <- NULL
rna_counts$Name <- NULL

# This code creates a subdirectory to store the output in the same working directory.
output_dir <- "./ORFik_output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# This prepares the ribo_counts
anyNA(sample_info[SeqType == "RIBO",]$SampleID)
anyNA(colnames(ribo_counts))

ribo_counts <- ribo_counts[,data.table::chmatch(sample_info[SeqType == "RIBO",]$SampleID, colnames(ribo_counts)), with = F]
ribo_counts <- as.matrix(ribo_counts)
#ribo_counts <- round(ribo_counts)
mode(ribo_counts) <- "integer"
RPF_counts <- SummarizedExperiment(assays = ribo_counts, colData = sample_info[sample_info$SeqType=="RIBO",])
rownames(RPF_counts) <- Gene_names


# This prepares the rna_counts

rna_counts <- rna_counts[,data.table::chmatch(sample_info[SeqType == "RNA",]$SampleID, colnames(rna_counts)), with = F]
rna_counts <- as.matrix(rna_counts)
#rna_counts <- round(rna_counts)
mode(rna_counts) <- "integer"
RNA_counts <- SummarizedExperiment(assays = rna_counts, colData = sample_info[sample_info$SeqType=="RNA",])
rownames(RNA_counts) <- Gene_names


# Correct for batch effect (optional)

batch.effect <- batch_b
be <- ifelse(batch.effect, "replicate + ", "")




sample_info_per_libtype <- copy(sample_info)
sample_info_per_libtype$SeqType = NULL

design <- colnames(colData(RPF_counts))
design <- design[-c(1,3,4)]
design

##### START THE ANALYSIS USING CONTRAST: CELL ####
target.contrast <- "Cell"

te.design <- as.formula(paste0("~ SeqType + ", be,
                               paste(design, collapse = " + "),
                               "+ SeqType:", target.contrast))

te.design

main.design <- as.formula(paste0("~ ", be, paste(design, collapse = " + ")))
message("----------------------")
message("Full exper. design: ", main.design)
message("Interaction design: ", te.design)
message("Target -- contrast: ", target.contrast)
message("----------------------")

# TE count table (cbind of both)

se <- cbind(assay(RPF_counts), assay(RNA_counts))

colData <- rbind(colData(RPF_counts), colData(RNA_counts))
combined_se <- SummarizedExperiment(list(counts = se),colData = colData)


## DESeq models (total: 3)


# TE
message("----------------------")
ddsMat_te <- DEG_DESeq(combined_se, te.design, "Model 1/3: TE")

# Ribo
ddsMat_ribo <- DEG_DESeq(RPF_counts, main.design, "Model 2/3: Ribo-seq")

# RNA
ddsMat_rna <- DEG_DESeq(RNA_counts, main.design, "Model 3/3: RNA-seq")

# Do result analysis: per contrast selected

cell <- sample_info$Cell
pairs <- ORFik:::combn.pairs(unique(cell))

dt <- DTEG_model_results(ddsMat_rna, ddsMat_ribo, ddsMat_te,
                         target.contrast, pairs, p.value = 0.05,
                         complex.categories = F)


###### ANALYSIS ONLY FOR CONTROL ######

colData(combined_se)
combined_se[, colData(combined_se)$Condition=="ctrl"]
assay(combined_se[, colData(combined_se)$Condition=="ctrl"])
combined_se_ctrl <- combined_se[, colData(combined_se)$Condition=="ctrl"]
RPF_counts_ctrl <- RPF_counts[, colData(RPF_counts)$Condition=="ctrl"]
RNA_counts_ctrl <-RNA_counts[, colData(RNA_counts)$Condition=="ctrl"]


design_ctrl <- design[2]
te.design.ctrl <- as.formula(paste0("~ SeqType + ", be,
                                    paste(design_ctrl, collapse = " + "),
                                    "+ SeqType:", target.contrast))

te.design.ctrl

main.design.ctrl <- as.formula(paste0("~ ", be, paste(design_ctrl, collapse = " + ")))
message("----------------------")
message("Full exper. design: ", main.design.ctrl)
message("Interaction design: ", te.design.ctrl)
message("Target -- contrast: ", target.contrast)
message("----------------------")


# TE
message("----------------------")
ddsMat_te <- DEG_DESeq(combined_se_ctrl, te.design.ctrl, "Model 1/3: TE")

# Ribo
ddsMat_ribo <- DEG_DESeq(RPF_counts_ctrl, main.design.ctrl, "Model 2/3: Ribo-seq")

# RNA
ddsMat_rna <- DEG_DESeq(RNA_counts_ctrl, main.design.ctrl, "Model 3/3: RNA-seq")

# Do result analysis: per contrast selected
cell <- sample_info$Cell
pairs <- ORFik:::combn.pairs(unique(cell))

dt <- DTEG_model_results(ddsMat_rna, ddsMat_ribo, ddsMat_te,
                         target.contrast, pairs, p.value = 0.05,
                         complex.categories = F)

## Plotting ##

output.dir <- "./ORFik_output"
#plot.title <- "Condition: Adherence - VC vs BP"
baseline <- pairs[[1]][2]
case <- pairs[[1]][1]
plot.title <- paste("Condition: ctrl |", case, "(target) vs", baseline, "(baseline)")
plot <- DTEG.plot(dt, output.dir, p.value = 0.05, plot.title, width = 6, height = 6,
                  relative.name = "ctrl_subset_contrast_cell.pdf")

file_name <- paste("./ORFik_output/ctrl_subset", case, "vs", baseline)
file_name <- gsub(" ", "_", file_name)
fwrite(dt, file = file_name)


###### ANALYSIS ONLY FOR ANOI ######

colData(combined_se)
combined_se[, colData(combined_se)$Condition=="anoi"]
assay(combined_se[, colData(combined_se)$Condition=="anoi"])
combined_se_anoi <- combined_se[, colData(combined_se)$Condition=="anoi"]
RPF_counts_anoi <- RPF_counts[, colData(RPF_counts)$Condition=="anoi"]
RNA_counts_anoi <-RNA_counts[, colData(RNA_counts)$Condition=="anoi"]


design_anoi <- design[2]
te.design.anoi <- as.formula(paste0("~ SeqType + ", be,
                                    paste(design_anoi, collapse = " + "),
                                    "+ SeqType:", target.contrast))

te.design.anoi

main.design.anoi <- as.formula(paste0("~ ", be, paste(design_anoi, collapse = " + ")))
message("----------------------")
message("Full exper. design: ", main.design.anoi)
message("Interaction design: ", te.design.anoi)
message("Target -- contrast: ", target.contrast)
message("----------------------")


# TE
message("----------------------")
ddsMat_te <- DEG_DESeq(combined_se_anoi, te.design.anoi, "Model 1/3: TE")

# Ribo
ddsMat_ribo <- DEG_DESeq(RPF_counts_anoi, main.design.anoi, "Model 2/3: Ribo-seq")

# RNA
ddsMat_rna <- DEG_DESeq(RNA_counts_anoi, main.design.anoi, "Model 3/3: RNA-seq")

# Do result analysis: per contrast selected
cell <- sample_info$Cell
pairs <- ORFik:::combn.pairs(unique(cell))

dt <- DTEG_model_results(ddsMat_rna, ddsMat_ribo, ddsMat_te,
                         target.contrast, pairs, p.value = 0.05,
                         complex.categories = F)

## Plotting ##
output.dir <- "./ORFik_output"
baseline <- pairs[[1]][2]
case <- pairs[[1]][1]
plot.title <- paste("Condition: anoi |", case, "(target) vs", baseline, "(baseline)")
plot <- DTEG.plot(dt, output.dir, p.value = 0.05, plot.title, width = 6, height = 6,
                  relative.name = "anoi_subset_contrast_cell.pdf")

file_name <- paste("./ORFik_output/anoi_subset", case, "vs", baseline)
file_name <- gsub(" ", "_", file_name)
fwrite(dt, file = file_name)


