import argparse
import os

parser = argparse.ArgumentParser(description='read fimo tsv file and create BED file centered around regions.')

parser.add_argument('-i', type=str, required=True, help='input fimo TSV.')
parser.add_argument('-o', type=str, required=False, default='./' ,help='output directory')
parser.add_argument('--cond', type=str, required=False, default='fimoparse', help='condition of fimo file [fimoparse].')
parser.add_argument('--pval', type=float, required=False, default=0.05, help='p-value cutoff for inclusion. default=0.05.')
parser.add_argument('--qval', type=float, required=False, default=1, help='q-value cutoff for inclusion. default=1 (no qvalue filtering).')
parser.add_argument('--USfasta', type=bool, required=False, default=True, help='expects fasta headers to be underscore delimited [e.g. chr1_200_300], if set to False, expects chr1:200-300')
args = parser.parse_args()

def retBedfromHead(col, delims):
    if not delims:
        chrom = col.split(':')[0]
        start = col.split(':')[1].split('-')[0]
        stop = col.split('-')[1]
    else:
        chrom = col.split('_')[0]
        start = col.split('_')[1]
        stop = col.split('_')[2]
    return chrom, int(start), int(stop)

with open(args.i) as fimoOut:
    bedFormat = {}
    trackDic = {}
    for line in fimoOut:
        if not line.startswith('motif') and not line.startswith('#') and len(line.strip().split()) > 0:
            pval = float(line.strip().split()[7])
            qval = float(line.strip().split()[8])
            if pval <= args.pval and qval <= args.qval:
                motID= line.strip().split()[0]
                if motID not in trackDic:
                    trackDic[motID] = []
                seqName = line.strip().split()[2]
                chrStr, seqStart, seqStop = retBedfromHead(seqName, args.USfasta)
                motStart = int(line.strip().split()[3])
                motStop = int(line.strip().split()[4])
                trackStr = seqName + str(motStart) + str(motStop)
                if trackStr not in trackDic[motID]:
                    if motID not in bedFormat:
                        bedFormat[motID] = []
                    bedFormat[motID].append([str(chrStr), str((seqStart+motStart)-0 ),str((seqStart+motStop) + 1)])
                    trackDic[motID].append(trackStr)

for i in bedFormat:
    outFile = os.path.join(args.o, i + "." + args.cond + ".bed")
    with open(outFile, 'w') as f:
        f.writelines('\t'.join(entry) + '\n' for entry in bedFormat[i])
