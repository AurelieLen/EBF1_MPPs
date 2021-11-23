import yaml
import os
import pandas as pd
import numpy as np
import pdfkit

# Since we clustered motifs on 2 levels (1 condition specific, 1 over all conditions), things start to get messy.
# In this script we parse over all the clustering steps, retaining information on the 'base level' e.g. motifs as they are in jaspar.
# It creates a rather excessive file, where every cluster has a consensus motif, the individual motifs that are clustered and expression values.
# The idea is to make it more clear on how the manual name for a cluster was chosen.

yamlFile = "../../data/ATAC/coMotifEnrichment/enrichedMots_cluster/enrichedMots_clusters.yml"
exprTSV = '../../data/RNA/MPP3_KOvsWT.tsv'

with open(yamlFile) as f:
    motDic = yaml.safe_load(f)

# Single yamls.
singleYaml = [
    '../../data/ATAC/motifCluster/MPP3_WT_clus/MPP3_WT_clus_clusters.yml',
    '../../data/ATAC/motifCluster/MPP4_WT_clus/MPP4_WT_clus_clusters.yml',
    '../../data/ATAC/motifCluster/MPP3_KO_clus/MPP3_KO_clus_clusters.yml',
    '../../data/ATAC/motifCluster/MPP4_KO_clus/MPP4_KO_clus_clusters.yml'
]

singleMotDics = {}
for i in singleYaml:
    name = i.split('/')[-1].replace("_clus_clusters.yml", "").replace("_","")
    with open(i) as f:
        singleMotDics[name] = yaml.safe_load(f)

expr = pd.read_csv(exprTSV, index_col=0, sep='\t')
expr.columns = ['baseMean','log2FoldChange','lfcSE','pvalue','padj']
expr['symbol'] = expr.index


# Populate the final dic with all information.
finDic = {}
sampleDic = {'MPP3WT':'MPP3_WT','MPP3KO':'MPP3_KO','MPP4WT':'MPP4_WT','MPP4KO':'MPP4_KO'}
for cluster in motDic:
    finDic[cluster] = []
    for motif in motDic[cluster]:
        # These motifs are actually clusters (or singletons), from our first round of clustering.
        # Where do they come from ?
        sample = motif.split('_')[0]
        if 'Cluster' in motif.split('_')[1]:
            sampleMotif = motif.split(' ')[0].replace(sample + "_", '')
            for singlemot in singleMotDics[sample][sampleMotif]:
                singlemotPNG = '../../LFS/genomeData/subsetMots/' + singlemot.split(' ')[0] + '.png'
                finDic[cluster].append( [sample, singlemot, singlemotPNG] )
        else:
            sampleMotif = motif.split('_')[1].split(' ')[0]
            for singleclus in singleMotDics[sample]:
                for singlemot in singleMotDics[sample][singleclus]:
                    if singlemot.split(' ')[0] == sampleMotif:
                        singlemotPNG = '../../LFS/genomeData/subsetMots/' + singlemot.split(' ')[0] + '.png'
                        finDic[cluster].append( [sample, singlemot, singlemotPNG] )

finExprDic = {}
for cluster in finDic:
    finExprDic[cluster] = []
    for lis in finDic[cluster]:
        sample = lis[0]
        mot = lis[1]
        gene = mot.split(' ')[1]
        if '::' in gene:
            exprs = []
            for factor in gene.split('::'):
                if factor.lower().capitalize() not in list(expr['symbol']):
                    #print("{} not found, appending 0.".format(factor.lower().capitalize()))
                    exprs.append(float(0))
                else:
                    exprs.append(np.mean(float(expr[expr['symbol'] == factor.lower().capitalize()]['baseMean']) ))
        elif '(' in gene:
            exprs = []
            factor = gene.split('(')[0].lower().capitalize()
            if factor not in list(expr['symbol']):
                #print("{} not found, appending 0.".format(factor))
                exprs.append(float(0))
            else:
                exprs.append(float(np.mean(expr[expr['symbol'] == factor]['baseMean'])))
        else:
            exprs = []
            factor = gene.lower().capitalize()
            if factor not in list(expr['symbol']):
                #print("{} not found, appending 0.".format(factor))
                exprs.append(float(0))
            else:
                exprs.append(float(np.mean(expr[expr['symbol'] == factor]['baseMean'])))
        finExprDic[cluster].append( [sample, mot,lis[2], np.round( np.mean(exprs), 2) ] )


finExprDicRep = {}
for cluster in finExprDic:
    max = 0
    for lis in finExprDic[cluster]:
        if lis[3] > max:
            max = lis[3]
    finExprDicRep[cluster] = []
    for lis in finExprDic[cluster]:
        if lis[3] == max:
            finExprDicRep[cluster].append( [lis[0], lis[1], lis[2],lis[3],'rep'] )
        else:
            finExprDicRep[cluster].append( [lis[0], lis[1], lis[2],lis[3], 'norep'] )


# Fetch the actual names the consensus motifs had in meme file:
memeNames = []
with open('../../data/ATAC/coMotifEnrichment/enrichedMots_cluster/enrichedMots_consensus_motifs.meme') as f:
    for line in f:
        if line.strip().startswith('MOTIF'):
            memeNames.append(line.strip().replace("MOTIF ", ""))
cluster2memeDic = {}
for i in memeNames:
    if i.startswith('Cluster'):
        if i.split(' ')[0] in finExprDicRep:
            cluster2memeDic[i.split(' ')[0]] = i
        else:
            print("This can't be right. Double check {}".format(i))
    else:
        for key in motDic:
            for hit in motDic[key]:
                if hit.split(' ')[0] == i.split(' ')[0]:
                    cluster2memeDic[key] = i
finExprDicRepMeme = {}
for cluster in finExprDicRep:
    finExprDicRepMeme[cluster] = []
    for lis in finExprDicRep[cluster]:
        finExprDicRepMeme[cluster].append( [lis[0], lis[1], lis[2], lis[3], lis[4], cluster2memeDic[cluster] ] )

manualNameDic = {   
    "MPP3WT_MA0739.1 MPP3WT_Cl7": "Hic1",
    "Cluster_7 MPP3WT_Cl4,MPP4KO_Cl5,MPP4KO_Cl7(...)": "Kruppel-like",
    "MPP3WT_Cluster_5 MPP3WT_Cl5": "PKNOX1",
    "MPP3WT_MA0513.1 MPP3WT_Cl6": "SMAD",
    "MPP3WT_MA1581.1 MPP3WT_Cl10": "ZBTB6",
    "MPP3WT_Cluster_11 MPP3WT_Cl11": "Zincfinger",
    "Cluster_1 MPP4KO_Cl6,MPP4KO_Cl3": "EGR1",
    "MPP4WT_Cluster_5 MPP4WT_Cl5": "Homeobox2",
    "MPP4WT_MA0041.1 MPP4WT_Cl14": "Foxd3",
    "MPP4WT_MA0850.1 MPP4WT_Cl15": "FOXP3",
    "MPP4WT_MA0109.1 MPP4WT_Cl10": "HLTF",
    "Cluster_18 MPP4WT_Cl19,MPP4WT_Cl18": "POU2F1",
    "MPP4WT_Cluster_11 MPP4WT_Cl11": "MEF2C",
    "MPP4WT_MA0676.1 MPP4WT_Cl8": "NR2E1",
    "MPP4WT_MA0687.1 MPP4WT_Cl1": "SPIC",
    "MPP4WT_Cluster_16 MPP4WT_Cl16": "TCF7L2",
    "MPP4WT_MA1585.1 MPP4WT_Cl2" : "ZKSCAN1",
    "MPP3KO_MA1536.1 MPP3KO_Cl17":"NR2C2",
    "MPP3KO_MA0506.1 MPP3KO_Cl14": "NRF1",
    "Cluster_16 MPP3KO_Cl12,MPP4WT_Cl17": "SOX15" ,
    "MPP3KO_MA0067.1 MPP3KO_Cl23": "Pax2",
    "Cluster_12 MPP3KO_Cl7,MPP4WT_Cl21,MPP4WT_Cl20": "TBP",
    "MPP3KO_MA1616.1 MPP3KO_Cl3": "Prdm15",
    "Cluster_29 MPP3KO_Cl19,MPP4WT_Cl3": "RFX1",
    "MPP3KO_Cluster_2 MPP3KO_Cl2": "SIX",
    "Cluster_4 MPP3WT_Cl14,MPP3WT_Cl3,MPP4KO_Cl1": "TCF3",
    "Cluster_5 MPP3WT_Cl8,MPP4KO_Cl4": "EBF1",
    "Cluster_22 MPP3WT_Cl1,MPP3WT_Cl2": "ETS-related",
    "Cluster_23 MPP3KO_Cl1,MPP4WT_Cl12" : "STAT/IRF",
    "Cluster_33 MPP3KO_Cl13,MPP3WT_Cl13,MPP3WT_Cl12" : "TAl1::TCF3",
    "Cluster_26 MPP3KO_Cl5,MPP3KO_Cl4,MPP4WT_Cl4" : "NFAT",
    "Cluster_39 MPP3KO_Cl20,MPP3WT_Cl9" : "NKX",
    "Cluster_14 MPP3KO_Cl6,MPP4WT_Cl22": "Arid5a",
    "Cluster_3 MPP3KO_Cl21,MPP4KO_Cl2": "Arnt/Hif",
    "Cluster_30 MPP3KO_Cl22,MPP3KO_Cl16,MPP3KO_Cl15": "CEBP/FOS/JUN",
    "Cluster_13 MPP3KO_Cl9,MPP4KO_Cl7,MPP4WT_Cl6": "KLF",
    "MPP3KO_Cluster_18 MPP3KO_Cl18": "E2F1",
    "Cluster_9 MPP3KO_Cl11,MPP4WT_Cl13": "FOXP1",
    "MPP3KO_Cluster_10 MPP3KO_Cl10": "Gfi",
    "Cluster_25 MPP3KO_Cl8,MPP4WT_Cl9": "Homeobox1"
    }

finExprDicRepMemeName = {}
for cluster in finExprDicRepMeme:
    finExprDicRepMemeName[cluster] = []
    for lis in finExprDicRepMeme[cluster]:
        finExprDicRepMemeName[cluster].append( [lis[0], lis[1], lis[2], lis[3], lis[4], lis[5], manualNameDic[ lis[5] ] ] )


# Find location of consensus images.
clusterConsDic = {}
for cluster in finExprDicRepMemeName:
    if len(motDic[cluster]) > 1:
        checkPath = os.path.join("../../data/ATAC/coMotifEnrichment/enrichedMots_cluster/consensus_motifs_img/" + cluster + "_consensus.png")
        if os.path.exists(checkPath):
            clusterConsDic[cluster] = checkPath
        else:
            print("Big problem with {}".format(checkPath) )
    else:
        checkPath = os.path.join("../../data/ATAC/coMotifEnrichment/enrichedMots_cluster/consensus_motifs_img/" + motDic[cluster][0].split(' ')[0] + "_consensus.png")
        if os.path.exists(checkPath):
            clusterConsDic[cluster] = checkPath
        else:
            print("Big problem with {}".format(checkPath))

finExprDicRepMemeNameCons = {}
for cluster in finExprDicRepMemeName:
    finExprDicRepMemeNameCons[ (cluster, clusterConsDic[cluster]) ] = []
    for lis in finExprDicRepMemeName[cluster]:
        finExprDicRepMemeNameCons[ (cluster, clusterConsDic[cluster]) ].append( [lis[0], lis[1], lis[2], lis[3], lis[4], lis[5], lis[6]] )

print(finExprDicRepMemeNameCons)
# Our final finExprDicRepMemeNameCons dictionary contains:
# keys: tuple(cluster, consensus img). The cluster is as defined by the highest level of clustering e.g.: 
# cluster level 1 = motifs enriched in MPP3WT, MPP3KO, MPP4WT, MPP4KO
# cluster level 2 = clustered motifs from level 1.
# the nested list associated to the dict contain (in order):
# sample (enriched motifs comes from e.g. MPP4KO)
# motif name (accession as it is in JASPAR database)
# single consensus path.
# expression
# rep/norep (string saying if it contains the highest counts of all members)
# memeName (string as it is named in the meme file.)
# manualName (how we named the cluster in the end)



message = """<html>
<head></head>
<body><p>MPP4WT!</p></body>
<table>
<tr>
<th>Cluster</th>
<th>Consensus</th>
<th>clusterName</th>
<th>baseMean</th>
<th>manualName</th>
<th>PWM</th>
<th>memeName</th>
</tr>
"""

# Iterate over the clusters.
for cluster in finExprDicRepMemeNameCons:
        for lis in finExprDicRepMemeNameCons[cluster]:
                if 'MPP4WT' in lis[5]:
                        message += "<tr>"
                        message += "<td>{}</td>".format( cluster[0] )
                        message += "<td><img src={}></td>".format( cluster[1] )
                        message += "<td>{}</td>".format( lis[5] )
                        message += "<td>{}</td>".format( lis[3] )
                        message += "<td>{}</td>".format( lis[6] )
                        message += "<td><img src={}></td>".format( lis[2] )
                        message += "<td>{}</td>".format( lis[1] )
                        message += "</tr>"
message += """
</table>
</html>"""

with open('MPP4WT.html', 'w') as f:
        f.write(message)
