import argparse

parser = argparse.ArgumentParser(
        description='Fetch motifs based on ame tsv results',
        prog='parseAmeresults.',
        )

parser.add_argument(
        '-i',
        type=str,
        required=True,
        help='input file (ame.tsv)'
        )

parser.add_argument(
        '--pval',
        type=float,
        required=False,
        default=None,
        help='p-value to filter on.'
        )

parser.add_argument(
        '--motifs',
        type=str,
        required=True,
        help='Path to motif file.'
        )

parser.add_argument(
        '-o',
        type=str,
        required=True,
        help='output file (without extension).'
        )

args = parser.parse_args()

keepers = []
with open(args.i) as f:
    for line in f:
        if not line.startswith('rank') and not line.startswith('#'):
            if len(line.strip().split()) > 1:
                ID = line.strip().split()[2]
                pval = float(line.strip().split()[5])
                Eval = float(line.strip().split()[7])
                if args.pval:
                    if pval < args.pval:
                        keepers.append(ID)
                else:
                    keepers.append(ID)

# parse JasparDB meme.
motNes = []

with open(args.motifs) as f:
    headerStatus = True
    motKeep = False
    for line in f:
        if line.startswith('MOTIF'):
            headerStatus = False
            if line.strip().split(' ')[1] in keepers:
                motNes.append(line.strip())
                motKeep = True
            else:
                motKeep = False
        elif headerStatus == True:
            motNes.append(line.strip())
        elif motKeep == True:
            motNes.append(line.strip())

with open(args.o, 'w') as f:
    for item in motNes:
        f.write("%s\n" % item)

