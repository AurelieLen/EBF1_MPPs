import argparse
import yaml

parser = argparse.ArgumentParser(
        description='rename motif file (meme) based on a cluster.yaml file (TOBIAS ClusterMotifs output.)',
        prog='parseAmeresults.',
        )

parser.add_argument(
        '-i',
        type=str,
        required=True,
        help='input file.meme'
        )

parser.add_argument(
        '--yaml',
        type=str,
        required=True,
        help='yaml file.'
        )

args = parser.parse_args()

with open(args.yaml) as f:
    clusDic = yaml.safe_load(f)

with open(args.i) as f:
    for line in f:
        if line.startswith('MOTIF'):
            clusName = line.strip().split(' ')[1]
            if clusName in clusDic:
                clusDesc = ",".join( [i.replace(' ','-') for i in clusDic[clusName]] )
                print("MOTIF {} {}".format(clusName, clusDesc.replace('(...)','')))
            else:
                fixStat = False
                for clus in clusDic:
                    for entry in clusDic[clus]:
                        if clusName in entry:
                            clusDesc =  ",".join( [i.replace(' ','-') for i in clusDic[clus]] )
                            print("MOTIF {} {}".format(clus, clusDesc.replace('(...)', '')))
                            fixStat = True
                if fixStat == False:
                    print("ERROR WITH: {}".format(clusName))
        else:
            print( line.strip() )
