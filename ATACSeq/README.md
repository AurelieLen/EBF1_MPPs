## EBF1 in HSPC/MPPs

This folder contains the code, scripts and environments used for the ATACSeq data.
Single scripts reside in code/scripts, while analysis notebooks (jupyter) sit under code/books.
Note that the notebooks refer to the single scripts where appropriate (relative path), and also refer to an LFS folder (Lenaerts_et_al_ATAC_LFS.tar) (available on zenodo: 10.5281/zenodo.5720286).

A number of conda environments are used:
  - ATACofthesnake.yml: main environment for differential ATAC analysis (Note, pip install from the repository. https://github.com/WardDeb/ATACofthesnake)
  - condaEnv.yml : main environment for ATAC-seq analysis
  - GO.yml : environment for GO analysis (ATAC-seq).
  - snakePipes.yml: environment for snakePipes (ATAC-seq & RNA-seq).
  - TOBIASenv.yml : environment for footprinting analysis (ATAC-seq).

Notebooks:
  - ATAC.ipynb:  
    - preprocessing ( trimming, alignment )
    - differential accessibility
    - basic plots (correlations, genomeTracks, heatMaps, GSEA/GO, RNA correlation)
  - ATAC_motifs.ipynb:  
    - motif enrichments
    - motif clustering
    - co-occurence analysis
    - aggregation plots (footprinting)
  - Public_ChIP.ipynb:  
    - preprocessing (trimming, alignment) of public ChIP data (GSE107242 , GSE146128 , GSE59636)
    - enhancer prediction
    - ChIP profiles in ATAC peaks
    - genomeTrack
