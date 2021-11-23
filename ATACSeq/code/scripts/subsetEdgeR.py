import argparse

parser = argparse.ArgumentParser(description='subset edgeR based on bed.')

parser.add_argument('--iBed', type=str, help='bed', required=True)
parser.add_argument('--table', type=str, help='edgeR tsv', required=True)

args = parser.parse_args()

bedL = []
with open(args.iBed) as f:
    for line in f:
        bedL.append('_'.join(line.strip().split()))
bedSet = set(bedL)

with open(args.table) as f:
    for line in f:
        if line.strip().split()[0] in bedSet:
            print(line.strip())
