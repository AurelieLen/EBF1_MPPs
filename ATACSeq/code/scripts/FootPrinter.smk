import os
import itertools


# -c "SlurmEasy -l clusterlogs -t {threads} -n {rule}" 
# vars:
# set these specifically for what run you need.
bamDir = 'input/bams'
peaks = 'input/peaks/peaks.bed'
blackList = 'input/blacklist.bed'
genome = 'input/genome.fa'
motifs = 'input/motifs.meme'


inBams = ['bam1_rep1.bam','bam1_rep2.bam','bam2_rep1.bam','bam2_rep2.bam']
mergedOut = {
    'bam1': ['bam1_rep1.bam', 'bam1_rep2.bam'],
    'bam2': ['bam2_rep1.bam', 'bam2_rep2.bam']
}

rule all:
    input:
        ## merge bam files, run ATACorrect and Footprintscores.
        # expand('mergedBam/{merge}.bam', merge = list(mergedOut.keys()) ),
        # expand('ATACorrect/{merge}_corrected.bw', merge = list(mergedOut.keys()) ),
        expand('ATACorrect/{merge}_footprints.bw', merge = list(mergedOut.keys()) ),
        ## cluster Motifs.
        'motifClus/motif_consensus_motifs.meme',
        'bindetect_output/bindetect_results.xlsx',
        expand('subsample/{merge}.bam', merge = list(mergedOut.keys()) ),
        expand('pyDNAse_{merge}/{merge}.WellingtonFootprints.wig', merge = list(mergedOut.keys()) )
    
rule mergeBam:
    input:
        lambda wildcards: expand( os.path.join(bamDir, '{bamFile}'), bamFile = mergedOut[wildcards.merge] )
    output:
        'mergedBam/{merge}.bam'
    log:
        out = 'logs/sambamba_{merge}.out',
        err = 'logs/sambamba_{merge}.err'
    threads: 15
    shell:'''
    sambamba merge -t 15 {output} {input} > {log.out} 2> {log.err}
    '''

rule ATACorrect:
    input:
        'mergedBam/{merge}.bam'
    output:
        'ATACorrect/{merge}_corrected.bw'
    log:
        out = 'logs/ATACorrect_{merge}.out',
        err = 'logs/ATACorrect_{merge}.err'
    params:
        genome = genome,
        peaks = peaks,
        bl = blackList
    threads: 15
    shell:'''
    TOBIAS ATACorrect --bam {input} --genome {params.genome} \\
    --peaks {params.peaks} --blacklist {params.bl} --outdir ATACorrect \\
    --cores {threads} > {log.out} 2> {log.err}
    '''

rule FPScore:
    input:
        'ATACorrect/{merge}_corrected.bw'
    output:
        'ATACorrect/{merge}_footprints.bw'
    threads: 15
    log:
        out = 'logs/FootprintScores_{merge}.out',
        err = 'logs/FootprintScores_{merge}.err'
    params:
        peaks = peaks
    shell:'''
    TOBIAS FootprintScores --signal {input} --regions {params.peaks} \\
    --output {output} --cores {threads} > {log.out} 2> {log.err}
    '''

rule clusterMots:
    input: motifs
    output:
        'motifClus/motif_consensus_motifs.meme'
    log:
        out = 'logs/clusterMots.out',
        err = 'logs/clusterMots.err'
    shell:'''
    TOBIAS ClusterMotifs -m {input} -t 0.4 -a meme \\
    --dist_method 'seqcor' --clust_method 'complete' \\
    -p motif -o motifClus > {log.out} 2> {log.err}
    '''

rule BINDetect:
    input:
        mots =  'motifClus/motif_consensus_motifs.meme',
        signal1 = 'ATACorrect/' + list(mergedOut.keys())[0] + '_footprints.bw',
        signal2 = 'ATACorrect/' + list(mergedOut.keys())[1] + '_footprints.bw'
    output:
        'bindetect_output/bindetect_results.xlsx'
    params:
        peaks = peaks,
        genome = genome,
        condNames = list(mergedOut.keys())[0] + ' ' + list(mergedOut.keys())[1]
    log:
        out = 'logs/BINDetect.out',
        err = 'logs/BINDetect.err'
    threads: 15
    shell:'''
    TOBIAS BINDetect --signals {input.signal1} {input.signal2} \\
    --peaks {params.peaks} --genome {params.genome} --motifs {input.mots} \\
    --cond-names {params.condNames} --cores {threads} > {log.out} 2> {log.err}
    '''

rule subSample:
    input:
        'mergedBam/{merge}.bam'
    output:
        'subsample/{merge}.bam'
    log:
        out = 'logs/subSample_{merge}.out',
        err = 'logs/subSample_{merge}.out'
    threads: 15
    shell:'''
    count=$(samtools idxstats {input} | cut -f3 | awk 'BEGIN {{total=0}} {{total += $1}} END {{print total}}')
    FRAC=$(bc -l <<< 20000000/$count)
    echo $FRAC > {log.out}
    sambamba view -f bam -t {threads} -s $FRAC {input} -o {output} >> {log.out} 2> {log.err}
    '''

rule wellington:
    input:
        'subsample/{merge}.bam'
    output:
        'pyDNAse_{merge}/{merge}.WellingtonFootprints.wig'
    log:
        out = 'logs/pyDNAse_{merge}.out',
        err = 'logs/pyDNAse_{merge}.err'
    params:
        peaks = peaks,
        outDir = lambda wildcards: 'pyDNAse_' + wildcards.merge,
        prefix = lambda wildcards: wildcards.merge
    threads: 15
    shell:'''
    wellington_footprints.py -b -A -p 15 -o {params.prefix} {params.peaks} {input} {params.outDir} > {log.out} 2> {log.err}
    '''
