import pybedtools
import pandas as pd
import argparse
from rich import print
import sys
import os
import shutil
from multiprocessing import Pool
import subprocess
from itertools import combinations
import glob
import numpy as np
import seaborn as sns

# code based on Rscript from Pierre Cauchy.


def parseArguments():
    parser = argparse.ArgumentParser(description='coMotif bootstrapping.')
    parser.add_argument('--iBed', type=str, required=True, help='BED file containing the sample sites.')
    parser.add_argument('--bgBed', type=str, required=True, help='BED file containing the background sites.')
    parser.add_argument('--motifs', type=str, required=True, help='motifs in meme format. (at least 2).')
    parser.add_argument('--fasta', type=str, required=True, help='Genome fasta file.')
    parser.add_argument('--reps', type=int, default=1000, required=False, help='Number of bootstraps.')
    parser.add_argument('--diskHammer', action='store_true', default=False, required=False, help='Force continue if reps > 5000')
    parser.add_argument('--threads', type=int, default=4, required=False, help='Force continue if reps > 5000')
    parser.add_argument('--tmpdir', default='/tmp', required=False, type=str, help='Set a tmp dir.')
    parser.add_argument('--prefix', default='', required=False, type=str, help='prefix for output file names.')
    args = parser.parse_args()
    return args

# Checks:
# bgBed should be >>> iBed. (100:1 ?).
# motifs > 1.
# pybedtools puts the files on disk, so for large number of reps this could get nasty.

def bedNums(iBed, bgBed):
    if (len(bgBed)/len(iBed)) < 5:
        return False
    else:
        return True

def motLen(motFile):
    with open(motFile) as f:
        motCount = 0
        for line in f:
            if 'letter-probability' in line:
                motCount += 1
    if motCount > 1:
        return True
    else:
        return False

def diskHammer(reps, diskHammer):
    if reps > 5000 and diskHammer:
        return True, "[red]Reps more then 5000 but diskhammer invoked. Moving on.[/red]"
    elif reps > 5000 and not diskHammer:
        return False, "[red]Reps more then 5000 but diskhammer NOT invoked. Stopping.[/red]"
    else:
        return True, "[green]The number of reps seems feasible. Moving on.[/green]"

def bedSamples(bed, n, samples):
    bedLis = []
    for i in range(1,samples+1):
        bedLis.append(bed.random_subset(n=n))
    return bedLis

def runFimo(argLis):
    motifFile = argLis[0]
    fnaFile = argLis[1]
    outDir = fnaFile + '_FIMO'
    logFile = 'logs/' + fnaFile.split('/')[-1] + '_' + motifFile.split('/')[-1] + '.out'
    errFile = 'logs/' + fnaFile.split('/')[-1] + '_' + motifFile.split('/')[-1] + '.err'
    fimoCMD = ['fimo', '--o',outDir, '--max-stored-scores', "10000000", motifFile, fnaFile]
    with open(logFile, 'w') as f, open(errFile, 'w') as e:
        subprocess.run(fimoCMD, stdout=f, stderr=e)

def getMots(memeFile):
    motLis = []
    with open(memeFile) as f:
        for line in f:
            if line.startswith('MOTIF'):
                motLis.append(line.strip().replace("MOTIF ",""))
    return motLis

def pairwiseIntersect(argLis):
    tsvFile = argLis[0]
    fullMotifs = argLis[1]
    motsList = [i.split(' ')[0] for i in fullMotifs]
    motifDic = {}
    for smallMot, largeMot in zip(motsList, fullMotifs):
        motifDic[smallMot] = largeMot
    tmpDF = pd.read_csv( tsvFile, sep='\t', comment='#', header=0, skip_blank_lines=True )
    # Filter by pvalue.
    tmpDF = tmpDF[tmpDF['p-value'] < 1e-3]
    motSitesDic = {}
    for hit in motsList:
        if hit in list(tmpDF['motif_id']):
            hitSites = []
            for index,row in tmpDF[tmpDF['motif_id'] == hit][['sequence_name','start','stop']].iterrows():
                peakLis = row['sequence_name'].replace(':','_').replace('-','_').split('_')
                start = row['start']
                stop = row['start']
                hitSites.append( (peakLis[0], int(peakLis[1]) + start -25, int(peakLis[1]) + stop + 25) )
            motSitesDic[motifDic[hit]] = pybedtools.BedTool(hitSites)
        else:
            # if no motif hit was found, just create an empty bed.
            motSitesDic[motifDic[hit]] = pybedtools.BedTool("", from_string=True)
    lenLis = []
    for comb in combinations(fullMotifs, 2):
        bed1 = motSitesDic[comb[0]]
        bed2 = motSitesDic[comb[1]]
        intLen = len(bed1.intersect(bed2))
        lenLis.append([comb[0], comb[1], intLen ])
        lenLis.append([comb[1], comb[0], intLen ])
    mat = pd.DataFrame(0, index=fullMotifs, columns=fullMotifs)
    for comb in lenLis:
        mat.loc[comb[0],comb[1]] = comb[2]
    pybedtools.helpers.cleanup()
    return mat

def main():
    args = parseArguments()
    pybedtools.helpers.set_tempdir(args.tmpdir)
    # Populate vars.
    iBed = pybedtools.BedTool(args.iBed)
    bgBed = pybedtools.BedTool(args.bgBed)
    
    # Checks.
    if bedNums(iBed, bgBed):
        print("[green]Number of sites well balanced. Going on.[/green]")
    else:
        print("[red]Number of sites in background is too close to number of sites in iBed. Stopping. [/red]")
        sys.exit()
    if motLen(args.motifs):
        print("[green]More then 1 motif found. Going on.[/green]")
    else:
        print("[red] only 1 or 0 motifs found. You sure the file is correctly formatted ? [meme format]. Stopping. [/red]")
        sys.exit()
    diskStat, retStr = diskHammer(args.reps, args.diskHammer)
    if diskStat:
        print(retStr)
    else:
        print(retStr)
        sys.exit()
    
    # set random subsamples.

    print("[green]Building bed files. This takes a while.[/green]")
    bgSamples = bedSamples(bgBed, len(iBed), args.reps)

    # Write out fastas.
    if os.path.exists('tmpFasta'):
        print("[red] tmpFasta folder found. Removing.[/red]")
        shutil.rmtree('tmpFasta')
    os.mkdir('tmpFasta')
    counter = 0
    iBed.sequence(fi=args.fasta, fo='tmpFasta/samples.fasta')
    fnas = ['tmpFasta/samples.fasta']
    for sample in bgSamples:
        outStr = 'tmpFasta/tmp_' + str(counter) + '.fna'
        fnas.append(outStr)
        sample.sequence(fi=args.fasta, fo=outStr)
        counter += 1
    
    # run fimos.
    if not os.path.exists('logs'):
        os.mkdir('logs')
    argLists = []
    for fnaFile in fnas:
        argLists.append([ args.motifs, fnaFile ])
    p = Pool(args.threads)
    p.map(runFimo, argLists)

    # Combinations.
    # parseFimo
    argLis = ["tmpFasta/samples.fasta_FIMO/fimo.tsv", getMots(args.motifs)]
    samples_IntersectMat = pairwiseIntersect(argLis)
    print(samples_IntersectMat)

    # backGrounds.
    bgLists = []
    argLis = []
    for tsv in glob.glob('tmpFasta/tmp*FIMO/fimo.tsv'):
        argLis.append( [tsv, getMots(args.motifs) ] )
    p = Pool(args.threads)
    bgList = p.map(pairwiseIntersect, argLis)
    arrayList = []
    for df in bgList:
        arrayList.append(df.to_numpy())
    meanArray = np.mean( np.array(arrayList), axis=0 )
    sdArray = np.std( np.array(arrayList), axis=0 )
    with open("bgArr.txt",'w') as f:
        for i in np.array(arrayList):
            np.savetxt(f, i, header=','.join(list(samples_IntersectMat.columns)) )
    with open("countArr.txt",'w') as f:
        for i in samples_IntersectMat.to_numpy():
            np.savetxt(f, i, header=','.join(list(samples_IntersectMat.columns)) )

    #Zscores
    z_scores = ( (samples_IntersectMat.to_numpy() + 1 ) - (meanArray+ 1) ) / (sdArray + 1)
    z_scoreDF = pd.DataFrame(z_scores)
    z_scoreDF.columns = samples_IntersectMat.columns
    z_scoreDF.index = samples_IntersectMat.index

    z_scoreDF.to_csv(args.prefix + 'Zscores.tsv', sep='\t')

    g = sns.clustermap(z_scoreDF)
    g.savefig(args.prefix + 'Zscores.png')

    # populate argLis for the backgrounds.

    # fimoFiles = glob.glob('tmpFasta/tmp*FIMO/fimo.tsv')

    # bgMats = []
    # for tsv in fimoFiles:
    #     bed = tsv.split('/')[1].replace('FIMO','')
    #     print(bed)
    #     motSitesDic = parseFimoTSV(tsv)
    #     bgMats = pairwiseIntersect(motSitesDic)
    # print(bgMats)



if __name__ == "__main__":
    main()

# Get sampled bedFiles from bgBed.
