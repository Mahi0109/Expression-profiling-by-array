#WRITE FIST LIBRARY SELECTION AND RUN ALL IMPORTANT LIBRARIES
library(org.Hs.eg.db)
library(AnnotationDbi)
library(biomaRt)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
BiocManager::install("biomaRt")
BiocManager::install("EnsDb.Hsapiens.v86")

library(GEOquery)
#first upload a file and setup or untar 
getGEOSuppFiles("GSE325008")
file.info("C:/Users/Mahima/Downloads/GSE325008_RAW.tar")
utils::untar(
  "C:/Users/Mahima/Downloads/GSE325008_RAW.tar",
  list = TRUE,
  tar = "internal"
)
tools::file_ext("C:/Users/Mahima/Downloads/GSE325008_RAW.tar")

#after extract all file 
dir.create("GSE325008_RAW", showWarnings = FALSE)

utils::untar(
  "C:/Users/Mahima/Downloads/GSE325008_RAW.tar",
  exdir = "C:/Users/Mahima/Downloads/GSE325008_RAW",
  tar = "internal"
)
# method 1 biomaRr read one untar with library to see all raw file and extact

library(data.table)

file <- "C:/Users/Mahima/Downloads/GSE325008_RAW/GSM9590859_HPNE.txt.gz"

raw <- fread(file)
colnames(raw)
con <- gzfile("C:/Users/Mahima/Downloads/GSE325008_RAW/GSM9590859_HPNE.txt.gz", "rt")


raw <- fread(
  "C:/Users/Mahima/Downloads/GSE325008_RAW/GSM9590859_HPNE.txt.gz",
  sep = "\t",
  header = FALSE,
  fill = TRUE
)
gene_table <- raw[, .(V13, V14, V15, V21)]
gene_table <- gene_table

genes <- unique(gene_table$GeneSymbol)
head(genes)
length(genes)

gene_table <- raw[, c("ProbeName","GeneName","SystematicName","gProcessedSignal")]
gene_table <- raw[, c("V13","V14","V15","V21")]
# 
colnames(gene_table) <- c(
  "ProbeName",
  "GeneName",
  "SystematicName",
  "gProcessedSignal"
)
gene_table <- gene_table[
  !(GeneName %in% c(
    "GE_BrightCorner",
    "DarkCorner"
  )),
]
gene_table <- gene_table[-c(1:10), ]

genes <- unique(gene_table$GeneName)

length(genes)
head(genes, 20)
# methde 2 and 3 combine 

library(org.Hs.eg.db)
library(AnnotationDbi)

ensembl_ids <- mapIds(
  org.Hs.eg.db,
  keys = genes,
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"
)
head(ensembl_ids)

annotation <- data.frame(
  GeneName = names(ensembl_ids),
  EnsemblID = as.character(ensembl_ids)
)

sum(!is.na(annotation$EnsemblID))   # mapped
sum(is.na(annotation$EnsemblID))    # unmapped

library(biomaRt)

mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

na_genes <- annotation$GeneName[is.na(annotation$EnsemblID)]

extra <- getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id"
  ),
  filters = "hgnc_symbol",
  values = na_genes,
  mart = mart
)

head(extra)

enst_ids <- grep("^ENST", genes, value = TRUE)

enst_map <- getBM(
  attributes = c(
    "ensembl_transcript_id",
    "ensembl_gene_id"
  ),
  filters = "ensembl_transcript_id",
  values = enst_ids,
  mart = mart
)

head(enst_map)

sum(!is.na(annotation$EnsemblID))
sum(is.na(annotation$EnsemblID))

# download maine file to create GO and KEGG enrichmentpart
write.csv(
  extra,
  "GSE325008_expression_with_ENSG.csv",
  row.names = FALSE
)
