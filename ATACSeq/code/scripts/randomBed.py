import argparse
import random
import os

parser = argparse.ArgumentParser(description='read BED file and create random subsets.')

parser.add_argument('-i', type=str, required=True, help='input bed file.')
parser.add_argument('-o', type=str, required=False, default='./' ,help='output directory')
parser.add_argument('-l', type=int, required=True, help='Number of entries within 1 subsampled bed file.')
parser.add_argument('-n', type=int, required=True, help='How many subsamplings should be done.')

args = parser.parse_args()

os.mkdir(args.o)

#Read bed
bedLis = []
with open(args.i) as f:
    for line in f:
        bedLis.append(line.strip())

# subsample
for i in range(1,args.n+1):
    sampleStr = "subsam_" + str(i) + '.bed'
    subsam = random.sample(bedLis, args.l)
    with open(args.o + '/' + sampleStr, 'w') as f:
        for entr in subsam:
            f.write(entr + '\n')
    if i%10 == 0:
        print("{} written.".format(i))
