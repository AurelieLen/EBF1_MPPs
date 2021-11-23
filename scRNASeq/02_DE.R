# ---
# jupyter:
#   jupytext:
#     formats: ipynb,R:percent
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# %%
library(data.table)
library(ggplot2)
library(dplyr)
library(SingleCellExperiment)
library(anndata)
library(scater)
library(edgeR)
source('utils/funs.R')

#ggplot theme
theme_paper = theme_bw(base_size = 18, base_family = "",
                         base_line_size = 0.5, base_rect_size = 0.5)
theme_set(theme_paper)

#output directories
procdata_dir = './procdata/02script/'
dir.create(procdata_dir, showWarnings = FALSE)

figures_dir = './figures/02script/'
dir.create(figures_dir, showWarnings = FALSE)

depbdir_b1b2 = paste0(procdata_dir, 'clusterDEpb_b1b2/')
dir.create(depbdir_b1b2, showWarnings = FALSE)

# %% [markdown] tags=[]
# ## Loading and preparing data

# %%
data = read_h5ad('./procdata/01script/aur4_comb_proc.h5ad')


sce <- SingleCellExperiment(assays = List(counts = t(data$layers[['counts']])),
                            colData=data$obs,
                            rowData=data$var,
                            reducedDims = data$obsm)

colData(sce)$condition = factor(colData(sce)$condition, levels = c('WT', 'KO'))
colData(sce)$scaled_n_genes = scale(colData(sce)$n_genes)
sce = logNormCounts(sce)
g = plotReducedDim(sce, dimred = "X_umap_2d", colour_by = "leiden") + scale_color_discrete()
g

# %%
#Removing all clusters which have < 20 cells in combination: leiden-condition-sex
z = table(colData(sce)$leiden, colData(sce)$condition, colData(sce)$experiment, colData(sce)$sex)
z = data.table(z)
colnames(z) = c('leiden', 'condition', 'sex', 'experiment', 'count')
toexclude = unique(z[count< 20, leiden])
# toinclude = unique(z$leiden)[!unique(z$leiden) %in% toexclude]
sce = sce[,!as.character(colData(sce)$leiden) %in% toexclude]

colData(sce)$leiden = droplevels(colData(sce)$leiden)

# %% [markdown] tags=[]
# ## DE by cluster (pseudobulks)

# %%
colData(sce)$mouse = paste(colData(sce)$condition, colData(sce)$sex, colData(sce)$experiment, sep = '_')
pb = make_pseudobulks(sce, columns = c('leiden', "mouse"))

#Adding back metadata which is necessary for DE testing
colData(pb)$leiden = factor(gsub('([0-9]*)_(.*)', '\\1', colnames(pb)))
colData(pb)$condition = factor(gsub('([0-9]*)_(KO|WT)_(.*)', '\\2', colnames(pb)), levels = c('WT', 'KO'))
colData(pb)$experiment = gsub('(.*)(b[0-9]$)', '\\2', colnames(pb))
colData(pb)$sex = gsub('([0-9]*)_(KO|WT)_(female|male)_(.*)', '\\3', colnames(pb))


# %%
deb1b2_pb = edgeR_percluster(pb,
                      csvsave = paste0(depbdir_b1b2, 'DEb1b2'),
                      pdfsave = paste0(depbdir_b1b2, 'DEb1b2'),
                      model_formula = '~condition + sex + experiment',
                      expr_tr = 0,
                      cell_fraction = 0,
                      meanexpr_tr = 5,
                            FDR_tr = 0.1,
                            logFC_tr = 0.2)

# %%
pb = logNormCounts(pb)

# %% [markdown]
# Making a DE summary

# %%
for (i in names(deb1b2_pb)){
    deb1b2_pb[[i]]$cluster = i
}
deb1b2_pb_all = do.call(rbind, deb1b2_pb)

deb1b2_pb_allsig = sig(deb1b2_pb_all, FDR_tr = 0.1, logFC_tr = 0.2)

deb12_sum = deb1b2_pb_allsig[, .(experiment = 'b1b2',
                   UP = sum(logFC > 0),
                   DOWN = sum(logFC < 0),
                   DEno = .N),
                 by = 'cluster']
fwrite(deb12_sum, paste0(depbdir_b1b2, 'DEsum.csv'))
deb12_sum

# %%
fwrite(deb1b2_pb_all, paste0(depbdir_b1b2, 'DEallclusters.csv'))
fwrite(deb1b2_pb_allsig, paste0(depbdir_b1b2, 'DEallclusters_sig.csv'))

# %% [markdown]
# For each gene the number of clusters in which it is differentially expressed

# %%
x = deb1b2_pb_allsig[,.N, by = genesymbol][order(N, decreasing = TRUE),]
fwrite(x, paste0(depbdir_b1b2, 'DEsumpergene.csv'))
x

# %%
varslot = data$var
varslot$index = row.names(varslot)
fwrite(varslot, paste0(procdata_dir, '/varslot.csv'))

# %% [markdown]
# ## Plotting Lymphoid gene expression

# %%
pb = logNormCounts(pb)

# %%
lygenes = c('Igkc', 'Iglc1', 'Jchain', 'Ighm', 'Iglc1', 'Iglc2', 'Iglc3')
dt = as.data.table(get_flatgenes(pb, genes= lygenes, cellfeatures=c("condition", "leiden")))
dtmean = dt[,.(mean_logn = mean(logn)), by=.(genesymbol,condition)]

setkey(dt, leiden)
dt = dt[leiden %in% c('0', '1', '2', '4', '7', '14'),]

# %%
options(repr.plot.width=18, repr.plot.height=5)
g = ggplot(dt, aes(x = leiden, y =logn, color = condition, group=condition)) + 
    geom_point(position = position_dodge(width = 0.7), alpha = 0.7) + 
    stat_summary(fun="mean", geom="point", size=10,
             shape=95,
             position = position_dodge(width = 0.7)) + 
    facet_wrap(~genesymbol)
ggsave(paste0(depbdir_b1b2, 'Lygenes_expr.pdf'), g, width=15, height=7)
g
options(repr.plot.width=4, repr.plot.height=4)

# %% [markdown]
# ## DE on LSK middle clusters - 2, 1, 0, 9

# %% [markdown]
# Comparison to capture most of intermediate progenitors (overlapping with the MPP3 population) - KO vs WT.

# %%
sce_mid = sce[,colData(sce)$leiden %in% c('2', '1', '0', '9')]
g = plotReducedDim(sce_mid, dimred = "X_umap_2d", colour_by = "leiden") + scale_color_discrete()
g

# %%
colData(sce_mid)$mouse = paste(colData(sce_mid)$condition, colData(sce_mid)$sex, colData(sce_mid)$experiment, sep = '_')
pb_mid = make_pseudobulks(sce_mid, columns = c("mouse"))

#Adding back metadata which is necessary for DE testing
colData(pb_mid)$condition = factor(gsub('(KO|WT)_(.*)', '\\1', colnames(pb_mid)), levels = c('WT', 'KO'))
colData(pb_mid)$experiment = gsub('(.*)(b[0-9]$)', '\\2', colnames(pb_mid))
colData(pb_mid)$sex = gsub('(KO|WT)_(female|male)_(.*)', '\\2', colnames(pb_mid))

colData(pb_mid)$population = as.factor("cl0129")

# %%
depbdir_b1b2_cl2109 = paste0(procdata_dir, 'clusterDEpb_b1b2_cl2109/')
dir.create(depbdir_b1b2_cl2109, showWarnings = FALSE)

deb1b2_pb_mid = edgeR_percluster(pb_mid,
                                 cluster_column = 'population',
                                csvsave = paste0(depbdir_b1b2_cl2109, 'DEb1b2'),
                              pdfsave = paste0(depbdir_b1b2_cl2109, 'DEb1b2'),
                              model_formula = '~condition + sex + experiment',
                              expr_tr = 0,
                              cell_fraction = 0,
                              meanexpr_tr = 5,
                                FDR_tr = 0.1,
                                logFC_tr = 0.2)

# %%
pb_mid = logNormCounts(pb_mid)

# %% [markdown]
# Plotting all the lymphoid genes found DE

# %%
dt = as.data.table(get_flatgenes(pb_mid, genes= c('Igkc', 'Iglc1', 'Iglc2', 'Iglc3', 'Igha', 'Ighg1', 'Ighg2b', 'Igkv14-126', 'Ly6d', 'Jchain',
                               'Ighg2c'), cellfeatures=c("condition")))
dtmean = dt[,.(mean_logn = mean(logn)), by=.(genesymbol,condition)]

# %%
options(repr.plot.width=15, repr.plot.height=5)
g = ggplot(dt, aes(x = genesymbol, y =logn, color = condition, group=condition)) + 
geom_point(position = position_dodge(width = 0.7), alpha = 0.7) + 
stat_summary(fun="mean", geom="point", size=12,
             shape=95,
             position = position_dodge(width = 0.7))
ggsave(paste0(depbdir_b1b2_cl2109, 'Lygenes_expr.pdf'), g, width=15, height=5)
g
options(repr.plot.width=4, repr.plot.height=4)

# %% [markdown]
# ## Cell cycle score analysis

# %%
options(repr.plot.width=20, repr.plot.height=6)
x = as.data.table(colData(sce))

xsum = x[, .(meanS = mean(S_score), meanG2M = mean(G2M_score)), by = .(leiden, mouse, condition, experiment, sex)]
clusters = c('0', '1', '2', '4', '7', '14')
g1 = ggplot(xsum[leiden %in% clusters,], aes(x = leiden, y = meanS, color = condition)) + 
       geom_point(position = position_dodge(width = 0.7), alpha = 0.7, size = 4) + 
        stat_summary(fun="mean", geom="point", size=18,
             shape=95,
             position = position_dodge(width = 0.7))
print(g1)
ggsave(paste0(figures_dir, 'S_scores_percluster.pdf'), g1, width = 14)

g2 = ggplot(xsum[leiden %in% clusters,], aes(x = leiden, y = meanG2M, color = condition)) + 
       geom_point(position = position_dodge(width = 0.7), alpha = 0.7, size = 4) + 
        stat_summary(fun="mean", geom="point", size=18,
             shape=95,
             position = position_dodge(width = 0.7))
print(g2)
ggsave(paste0(figures_dir, 'G2M_scores_percluster.pdf'), g2, width = 14)
options(repr.plot.width=4, repr.plot.height=4)

# %%
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1, "group"]
  newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
      1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

# %%
options(repr.plot.width=20, repr.plot.height=6)

g1 = ggplot(x[leiden %in% clusters,], aes(x = leiden, y = S_score, fill = condition)) + 
       geom_split_violin(draw_quantiles = 0.5)

g2 = ggplot(x[leiden %in% clusters,], aes(x = leiden, y = G2M_score, fill = condition)) + 
       geom_split_violin(draw_quantiles = c(0.5))

ggsave(paste0(figures_dir, 'S_scores_violin.pdf'), g1, width = 14)
ggsave(paste0(figures_dir, 'G2M_scores_violin.pdf'), g2, width = 14)

g1
g2
options(repr.plot.width=4, repr.plot.height=4)
