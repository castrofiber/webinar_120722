library(tidyverse)
library(tximport)
library(DESeq2)
library(clusterProfiler)


setwd("~/Sandbox/Blastim_webinar/webinar_snake/results/Kallisto/") # change accordingly

sample_names <- c(paste0("WT", 1:3), paste0("GM", 1:3))

run_names <- c("SRR16481365",
               "SRR16481366",
               "SRR16481367",
               "SRR16481371",
               "SRR16481372",
               "SRR16481373")

kallisto_dirs <- paste0("./", run_names)

samples <- data.frame(sample = sample_names,
                      condition = c(rep("WT", 3), rep("GM", 3)),
                      path = kallisto_dirs)

files <- file.path(samples$path, "abundance.h5")
txi <- tximport(files, type = 'kallisto', txOut = T)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

dds <- ddsTxi[rowSums(counts(ddsTxi)) >= 10,]

dds$condition <- relevel(dds$condition, ref = "WT")

# Differential expression

dds <- DESeq(dds)

res <- results(dds, name = "condition_GM_vs_WT")
res[order(res$pvalue),]
sum(res$padj < 0.05, na.rm=TRUE)

library("pheatmap")
ntd <- normTransform(dds)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,"condition"])
pheatmap(assay(ntd)[select,], cluster_rows=T, show_rownames=T,
         cluster_cols=T, annotation_col=df)

plotPCA(ntd)

sign_results <- res %>%
  as.data.frame %>%
  rownames_to_column("gene_name") %>%
  mutate(ens_id = str_replace(gene_name, "_mRNA", "")) %>%
  filter(padj < .05)


sign_up <- sign_results %>% filter(log2FoldChange > 0)
sign_dw <- sign_results %>% filter(log2FoldChange < 0)

# Enrichment

KEGG_enrich <- enrichKEGG(sign_dw$ens_id, organism = "sce")

barplot(KEGG_enrich)
cnetplot(KEGG_enrich)
