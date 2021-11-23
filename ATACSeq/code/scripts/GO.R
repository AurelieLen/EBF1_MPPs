library("clusterProfiler")
library("org.Mm.eg.db")
library("msigdbr")

# MPP3 WT

WT_MPP3 <- read.delim('../../data/ATAC/output/MPP3_KOvsWT_edgeR.tsv')
WT_MPP3 <- WT_MPP3[WT_MPP3$gene_id != "",]
WT_MPP3 <- WT_MPP3[WT_MPP3$FDR < 0.1,]
WT_MPP3 <- WT_MPP3[WT_MPP3$logFC < 0,]
WT_MPP3_genes <- WT_MPP3$logFC
names(WT_MPP3_genes) <- WT_MPP3$gene_id
WT_MPP3_GO <- enrichGO(gene  = names(WT_MPP3_genes),
                           OrgDb         = org.Mm.eg.db,
                           keyType       = 'ENSEMBL',
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.5,
                           qvalueCutoff  = 0.5)

write.table(WT_MPP3_GO , file = "../../data/ATAC/GO/WT_MPP3_GO.tsv", row.names=FALSE, quote=FALSE, sep='\t')

# MPP3 KO

KO_MPP3 <- read.delim('../../data/ATAC/output/MPP3_KOvsWT_edgeR.tsv')
KO_MPP3 <- KO_MPP3[KO_MPP3$gene_id != "",]
KO_MPP3 <- KO_MPP3[KO_MPP3$FDR < 0.1,]
KO_MPP3 <- KO_MPP3[KO_MPP3$logFC > 0,]
KO_MPP3_genes <- KO_MPP3$logFC
names(KO_MPP3_genes) <- KO_MPP3$gene_id
KO_MPP3_GO <- enrichGO(gene  = names(KO_MPP3_genes),
                           OrgDb         = org.Mm.eg.db,
                           keyType       = 'ENSEMBL',
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.5,
                           qvalueCutoff  = 0.5)
write.table(KO_MPP3_GO , file = "../../data/ATAC/GO/KO_MPP3_GO.tsv", row.names=FALSE, quote=FALSE, sep='\t')

########

# MPP4 WT

WT_MPP4 <- read.delim('../../data/ATAC/output/MPP4_KOvsWT_edgeR.tsv')
WT_MPP4 <- WT_MPP4[WT_MPP4$gene_id != "",]
WT_MPP4 <- WT_MPP4[WT_MPP4$FDR < 0.1,]
WT_MPP4 <- WT_MPP4[WT_MPP4$logFC < 0,]
WT_MPP4_genes <- WT_MPP4$logFC
names(WT_MPP4_genes) <- WT_MPP4$gene_id
WT_MPP4_GO <- enrichGO(gene  = names(WT_MPP4_genes),
                           OrgDb         = org.Mm.eg.db,
                           keyType       = 'ENSEMBL',
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.5,
                           qvalueCutoff  = 0.5)
write.table(WT_MPP4_GO , file = "../../data/ATAC/GO/WT_MPP4_GO.tsv", row.names=FALSE, quote=FALSE, sep='\t')

# MPP4 KO

KO_MPP4 <- read.delim('../../data/ATAC/output/MPP4_KOvsWT_edgeR.tsv')
KO_MPP4 <- KO_MPP4[KO_MPP4$gene_id != "",]
KO_MPP4 <- KO_MPP4[KO_MPP4$FDR < 0.1,]
KO_MPP4 <- KO_MPP4[KO_MPP4$logFC > 0,]
KO_MPP4_genes <- KO_MPP4$logFC
names(KO_MPP4_genes) <- KO_MPP4$gene_id
KO_MPP4_GO <- enrichGO(gene  = names(KO_MPP4_genes),
                           OrgDb         = org.Mm.eg.db,
                           keyType       = 'ENSEMBL',
                           ont           = "all",
                           pAdjustMethod = "BH",,
                           pvalueCutoff  = 0.5,
                           qvalueCutoff  = 0.5)
write.table(KO_MPP4_GO , file = "../../data/ATAC/GO/KO_MPP4_GO.tsv", row.names=FALSE, quote=FALSE, sep='\t')
