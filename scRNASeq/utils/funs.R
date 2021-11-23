#General functions, should be useful for any analysis with minimal modifications

#' Calling significant genes on output from edgeR_DE function
sig = function(x, FDR_tr = 0.05, logFC_tr = 0.5){

    x = x[x$FDR < (FDR_tr) & abs(x$logFC) > (logFC_tr),]
  return(x)
}

#' Filtering based on expresion and/or fraction of cells/sampls expressing the gene
gene_filter = function(x, expr_tr = 0, cell_fraction = 0, meanexpr_tr = 0){
  n = ncol(x)
  above_tr = apply(x, 2, function(x) x >= expr_tr)
  expr_no = rowSums(above_tr)

  keep1 = expr_no/n >= cell_fraction
  print(paste0("Genes expressed (expr_tr >=", expr_tr, ") in ", cell_fraction, " of cells: ", sum(keep1)))
# Filtering by mean expr
  keep2 = rowMeans(x) >= meanexpr_tr
  print(paste0("Genes with mean expression >= ", meanexpr_tr, ": ", sum(keep2)))
  
  return(keep1 & keep2)
}

#' Running edgeR DE between conditions for each cluster
edgeR_percluster = function(sce,
                            cluster_column = 'leiden',
                            csvsave = 'DE',
                            pdfsave = NULL,
                            model_formula = '~condition + sex + scaled_n_genes',
                            coef = 'conditionKO',
                            expr_tr = 1,
                            cell_fraction = 0.075,
                            meanexpr_tr = 0,
                           FDR_tr = 0.05,
                           logFC_tr = 0.5){

  #Specifying clusters to compare (taking levels from the factor (converting to factor first if not a factor))
  clusters = colData(sce)[,cluster_column]
  if (is.factor(clusters)) {
    clusters = levels(colData(sce)[,cluster_column])
  }
  else {
    clusters = factor(colData(sce)[,cluster_column])
    clusters = levels(clusters)
  }

  delist = list() #output DE list
  if (!is.null(save)) pdf(paste0(pdfsave, '_DEfigures.pdf')) #Where to save plots
  for (i in clusters){
    print(paste0('DE analysis of cluster: ', i))
    sceI = sce[, as.character(colData(sce)[,cluster_column]) == i] #Subsetting for cells in relevant clusters
    sceI = sceI[rowSums(assays(sceI)$counts) > 0,] #Removing all 0 rows
    #Printing cell numbers to check
    print(table(colData(sceI)$sex, colData(sceI)$condition, droplevels(colData(sceI)[,cluster_column])))
    
    de = edgeR_DE(sceI,
                  plot_title = paste0('cluster', i),
                  model_formula = model_formula,
                  coef = coef,
                  expr_tr = expr_tr,
                  cell_fraction = cell_fraction,
                  meanexpr_tr = meanexpr_tr)

    write.csv(de, paste0(csvsave, '_DE_cluster', i, '.csv'))
    desig = de[de$FDR < FDR_tr & abs(de$logFC) > logFC_tr,]
    write.csv(desig, paste0(csvsave, '_DEsig_cluster', i, '.csv'))

    delist[[paste0('cluster', i)]]  = de
  }
  if (!is.null(save)) dev.off()
  return(delist)
}

#' EdgeR function for DE testing, general purpose, allows passing model formula and choice of the coefficient to test
#' Not implemented the contrasts
edgeR_DE = function(sce,  plot_title ='',
                    model_formula = '~condition + sex + scaled_n_genes',
                    coef = 'conditionKO',
                    expr_tr = 1,
                    cell_fraction = 0.075,
                    meanexpr_tr = 0){
  require("edgeR", quietly = TRUE)
  #' in edger expression can be obtained with; y$counts, cpm(y)

  y <- DGEList(counts=assays(sce)$counts, group = colData(sce)$condition)
  y <- calcNormFactors(y)

  #'Filtering low expressed genes
  #'Soneson Robinson paper 2018 shows that filtering as follows:
  #' After filtering, retaining only genes with an estimated expression above 1 TPM in more than 25% of the cells.
  #' reduces the FPRs in edgeR, MAST does not need this filtering

  #'Based on Dal Molin, Baruzzo, and Di Camillo 2017
  #'edgeR has higher recall but lower precision than MAST
  keep = gene_filter(y$counts, expr_tr = expr_tr, cell_fraction = cell_fraction, meanexpr_tr = meanexpr_tr)
  y <- y[keep, , keep.lib.sizes=TRUE]
  print(dim(y$counts))

  design = model.matrix(as.formula(model_formula), colData(sce))

  y <- estimateDisp(y, design, robust = TRUE)
  plotBCV(y, main = plot_title)
  fit <- glmQLFit(y, design, robust=TRUE)
  plotQLDisp(fit, main = plot_title) #Plotting the Q dispersions by mean counts

  lrt <- glmQLFTest(fit, coef = coef)

  destat = as.vector(decideTestsDGE(lrt, p.value = 0.05, lfc = 0.5))
  print(table(destat))
  print(plotMD(lrt, status = destat, main = paste0(plot_title, ' ', coef)))

  res = topTags(lrt, n = Inf, p.value = 1)$table
  ## ressig = res[res$FDR < 0.01 & abs(res$logFC) > 0.5, ]
  ## cpmsort = edgeR::cpm(y)[row.names(res),]

  colnames(res) = c("logFC", "logCPM", "F", "pvalue", "FDR")
  res$genesymbol = row.names(res)
  return(as.data.table(res))
}

#' Creating pseudobulk RNA-Seq profiles by summing up the counts
make_pseudobulks = function(sce, columns = c('leiden', 'mouse')){
  #' aggregate by cluster-sample
  require(Matrix.utils)
  require(SingleCellExperiment)
  z = t(assays(sce)$counts)
  groups <- colData(sce)[, columns]
  #column matrices
  pb <- aggregate.Matrix(as(z, 'CsparseMatrix'), groupings = groups, fun = "sum") 
  #'Reverting back to row-matrices just in case it causes any problems with SingleCellExperiment
  pb = as(t(pb), 'RsparseMatrix')
  pb <- SingleCellExperiment(assays = pb)
  names(assays(pb)) = 'counts'
  return(pb)
}

#' Plot multiple ggplots function
#'
#' Function which takes multiple ggpltos and plots them in a single window. Ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects).
#' #' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#' Source: Cookbook R website: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#' @param cols:   Number of columns in layout
#' @param layout: A matrix specifying the layout. If present, 'cols' is ignored.
#' @export
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' Convenience functions for extracting gene expression 
#' with annotation in a long format
get_flatgenes = function(sce, genes, cellfeatures = c("dpt_pseudotime", "condition")){

  require(data.table)

  df = assays(sce)$logcounts[genes, , drop = FALSE]
  df = t(df)
  df = as.matrix(df)
  dfm = reshape2::melt(df)
  colnames(dfm) = c("cellid", "genesymbol", "logn")
  dfm$cellid = as.character(dfm$cellid)
  dfm$genesymbol = as.character(dfm$genesymbol)

  toadd = colData(sce)[match(dfm$cellid, row.names(colData(sce))), cellfeatures, drop = FALSE]
  dfm = cbind(dfm, toadd)
  return(dfm)

}