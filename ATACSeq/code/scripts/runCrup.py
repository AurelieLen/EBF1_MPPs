import os
import subprocess
import pandas as pd

# all the bamFiles that need to be included in the analysis need to sit under bamLinks. 
# Results get written by default to 'mergeBAM'

# Steps: 
# 1. merge bam files
# 2. CRUP normalize
# 3. CRUP predict

bamDir = 'bamLinks'
outDir = 'mergeBAM'
if not os.path.exists(outDir): 
    os.mkdir(outDir) 
baseNames = [
    'H3K27ac_CLP',
    'H3K27ac_GMP',
    'H3K27ac_MPP3',
    'H3K27ac_MPP4',
    'H3K4me1_CLP',
    'H3K4me1_GMP',
    'H3K4me1_MPP3',
    'H3K4me1_MPP4',
    'H3K4me3_CLP',
    'H3K4me3_GMP',
    'H3K4me3_MPP3',
    'H3K4me3_MPP4']

baseRep = {}
for base in baseNames:
    tempLis = []
    for bamFile in os.listdir(bamDir):
        if base in bamFile:
            tempLis.append(os.path.join(bamDir, bamFile))
    baseRep[base] = tempLis

for base in baseRep:
    bamOut = os.path.join(outDir, base + '.bam')
    if not os.path.exists(bamOut):
        mergeCMD = ['samtools','merge','-@','20',bamOut] + baseRep[base]
        print(' '.join(mergeCMD))
        subprocess.run(mergeCMD)

# Define celltypes. 
cellTypes = [] 
for i in baseRep:
    if i.split('_')[1] not in cellTypes:
        cellTypes.append(i.split('_')[1])

print("Working with cellTypes {}".format(','.join(cellTypes)))

for cell in cellTypes:
    marks = ['H3K27ac','H3K4me1','H3K4me3']
    with open(cell, 'w') as f:
        f.write('feature\tbam_file\n')
        for mark in marks:
            for bamFile in os.listdir(outDir):
                if 'bai' not in bamFile and 'bam' in bamFile and mark in bamFile and cell in bamFile:
                    f.write(mark + '\t' + os.path.join(outDir, bamFile) + '\n')
    # run CRUP normalization
    crupN = ['Rscript','CRUP.R','-N','-I','-x','10','-i',cell, '-g', 'mm10','-s','single']
    print("Running crup Norm on {}".format(cell))
    subprocess.run(crupN)
    # Ship matrix.
    os.mkdir(cell + '_CRUP')
    os.rename(cell + '.data_matrix.rds',cell + '_CRUP/' + cell + '.data_matrix.rds')
    os.chdir(cell + '_CRUP') 
    print('Running crup pred on {}'.format(cell))
    crupP = ['Rscript','CRUP.R','-P','-x','10', '-m',cell+'.data_matrix.rds','-c','/data/manke/group/deboutte/repo/CRUP/DATA/CLASSIFIER/']
    subprocess.run(crupP)
    os.chdir('../')
