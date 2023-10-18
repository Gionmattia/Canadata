# Arg parser. Needs: path to files, annotation file, output path

args <- commandArgs(trailingOnly = TRUE)
path_to_file <- args[1]
path_to_annotation_file <- args[2]

# tximport is used to convert from Transcripts to Genes,
# GenomicFeatures creates the conversion table from transcriptID to geneID

library(tximport)
library(GenomicFeatures)

txdb <- makeTxDbFromGFF(path_to_annotation_file)

# Get the vector with all transcripts
transcripts <- keys(txdb, keytype = "TXNAME")

# Creates the tx2gene dataframe, by selecting the info from the txdb
tx2gene <- select(txdb, keys = transcripts, keytype = "TXNAME", columns = "GENEID") # nolint
tx2gene <- subset(tx2gene, select = c("TXNAME", "GENEID"))

# Define a function to remove the transcript versions
remove_string <- function(x) {
  gsub("\\..*", "", x)
}

# Remove the transcript version from tx2gene
tx2gene[, 1] <- apply(tx2gene[, 1, drop = FALSE], 1, remove_string)

# call tximport on the single file to get the output (it's a complex list)
output <- tximport(path_to_file, type = "salmon",
                    tx2gene = tx2gene, ignoreTxVersion = T)

# Obtain a df out of output
output_df <- data.frame(output$counts)

# Rename the column accordingly to sample name  # Still keeps those weird nts...
sample_name <- sub("_quant\\.sf$", "", file)
colnames(output_df)[1] <- sample_name

# Creates the new name and saves everything as a tsv in the working directory.
output_file_name <- paste0(path_to_output, sample_name, "_gene_counts.tsv")
write.table(output_df, output_file_name, sep = "\t", row.names = TRUE)

# OPERATIONS TO RETRIEVE THE UNCONVERTED IDs
# Read the original file
data <- read.csv(path_to_file, header = TRUE, sep = "\t")
data_mod <- data
# Removes the transcript ID version
data_mod[,1] <- apply(data[, 1, drop = FALSE], 1, remove_string)
# Retrieves missing transcripts
missing_transcripts <- data_mod$Name[!(data_mod$Name %in% tx2gene$TXNAME)]
# Creates a dataframe with the info of the missed transcripts for the file
missing_data <- data.frame()
for (root in missing_transcripts) {
  matching_rows <- grepl(root, data$Name)
  missing_data <- rbind(missing_data, data[matching_rows, ])

  # Assembles the file, once given the right name and saves it
  colnames(missing_data)[1] <- sample_name
  missing_data_name <- paste0(sample_name, "_missing_data.tsv")
  write.table(missing_data, missing_data_name, sep = "\t", row.names = TRUE)
  }


  
# I have decided to use an "agnostic" tx2gene table, in which I manually removed 
# the transcript version from the new annotation I used.
# This because the "ignoreTxVersion"of tximport wasn't working.

# I used the GENCODE release VM22 (released 06.2019 and based on the genome assembly GRCm38.p6) since
# It's the one that contains most of the transcript IDs identified by the salmon index adopted
# during the counting step

# The option to ignore the version is still not available. Nice. 