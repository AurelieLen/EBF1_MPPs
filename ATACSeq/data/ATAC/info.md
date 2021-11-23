## Preprocessing

Samples:

  - MPP3_KO1  
  - MPP3_KO2  
  - MPP3_WT1  
  - MPP3_WT2  
  - MPP4_KO1  
  - MPP4_KO2  
  - MPP4_WT1  
  - MPP4_WT2  

 1.  Reads are trimmed for nextera transposase adapter before snakePipes:

```
cutadapt -j 20 --minimum-length 25 -A CTGTCTCTTATA -a CTGTCTCTTATA -G GATGTGTATAAGAG -g CGATGTGTATAAGAG -o fqTrim/${sample}_R1.fastq.gz -p fqTrim/${sample}_R2.fastq.gz ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz
```

 2.  Run snakePipes for DNA-mapping (version 2.5.0) on fqTrim folder

```
DNA-mapping --trim --fastqc --mapq 3 --insertSizeMax 2000 -j 10 --dedup -i ./fqTrim/ -o DNA_mapq3_I2000_dedup mm10
```
Move filtered bam files (deduplicated) into 'bams' folder

 3.  run ATACofthesnake

```
ATAC --bamDir bams/ --outDir AOS_NFR --blackList blacklist.bed --genes genes.gtf --genomeFasta genome.fa --sampleSheet ss.tsv --motifs motifs.meme --downStream --fragSizeMax 150 --genomeSize 1.87e9 --mergeBam --condaPrefix /localenv/deboutte/anaconda/miniconda3
```

 - blacklist: akundaje blacklist (http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz)
 - genes.gtf: GRCm38 ensemble release 91
 - genome.fa: GRCm38
 - motifs: HOCOMOCOv11 CORE mouse motifs

 4.  output files (under Output folder)

 - peaks:  

   - MPP3_KO.bed: open sites in MPP3 KO
   - MPP3_WT.bed: open sites in MPP3 WT
   - MPP3_peaks.bed: union of all peaks in MPP3 samples
   - MPP4_KO.bed: open sites in MPP4 KO
   - MPP4_WT.bed: open sites in MPP4 WT
   - MPP4_peaks.bed: union of all peaks in MPP4 samples

 - motif enrichment results:

tsv files of motif enrichment results (part of ATACofthesnake).

   - MPP3_WT_ame.tsv 
   - MPP3_KO_ame.tsv 
   - MPP4_WT_ame.tsv 
   - MPP4_KO_ame.tsv 




MPP3_KO.bed
MPP3_peaks.bed
MPP3_WT.bed
MPP4_KO.bed
MPP4_peaks.bed
MPP4_WT.bed
