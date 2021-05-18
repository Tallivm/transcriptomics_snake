library(DESeq2)
library(EnhancedVolcano)

#file_name <- '../fcounts/fcount_result.txt'
#out_dir <- 'DE_results'
#prep_names <- c('Collibri', 'KAPA')

file_name <- snakemake@input[[1]]
out_dir <- snakemake@params[[1]]
prep_names <- snakemake@params[[2]]

dir.create(out_dir)

featurecounts <- read.delim(file_name, header=TRUE, skip=1)
column.from <- which(colnames(featurecounts)=="Length") + 1
geneid <- featurecounts[, 1]
featurecounts <- featurecounts[, c(column.from:ncol(featurecounts))]

prep_methods <- c(sub("(samtools_sorted.)([^_]*)_(.*)", "\\2", colnames(featurecounts)))

for (prep_name in prep_names) {
  sel <- featurecounts[, prep_methods==prep_name]
  types <- c(sub("(.*)(HBR|UHRR)(.*)", "\\2", colnames(sel)))
  types <- data.frame(condition=types, stringsAsFactors=TRUE)
  dds <- DESeqDataSetFromMatrix(countData=sel,
                                colData=types,
                                design=~condition)
  rownames(dds) <- geneid

  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]

  dds <- DESeq(dds)
  res <- results(dds, contrast=c("condition","UHRR","HBR"))
  DE <- res[res$padj < 0.1, ]

  write.csv(as.data.frame(DE),
            file=paste(out_dir, "/",prep_name, "_DE_results.csv", sep=''))

  pdf(paste(out_dir, "/", prep_name, "_DE_volcano.pdf", sep=''))
  plot(EnhancedVolcano(res,
                       title=paste("UHRR vs HBR", "-", prep_name),
                       lab=rownames(res),
                       x='log2FoldChange',
                       y='pvalue'))
  dev.off()
}
