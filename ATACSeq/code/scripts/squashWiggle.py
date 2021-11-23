import argparse
import sys

parser = argparse.ArgumentParser(description='go over wiggle file and replace negative with 0.')

parser.add_argument('-i', type=str, required=True, help='wiggle')
args = parser.parse_args()

with open(args.i) as f:
    for line in f:
        if not line.startswith('fixedStep'):
            stripLine = line.strip().split()
            try:
                if float(stripLine[3]) < 0:
                    print("{}\t{}\t{}\t{}".format(stripLine[0], stripLine[1],stripLine[2], 0.000000))
                else:
                    print("{}\t{}\t{}\t{}".format(stripLine[0], stripLine[1], stripLine[2], stripLine[3]))
            except:
                print("error with {}".format(line.strip()), file=sys.stderr)
