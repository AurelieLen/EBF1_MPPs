bamFiles = []
with open('names') as f:
	for line in f:
		bamFiles.append(line.strip())

bamDir = '../../../../snakePipes_hoxb8HeptadChIP/DNA_trim_dedup_map3/filtered_bam/'
ctrl = '../../../../snakePipes_hoxb8HeptadChIP/DNA_trim_dedup_map3/filtered_bam/IgG_ctrl.filtered.bam'
blackList = '../../blacklist.bed'

rule all:
	input:
		expand('{sample}.log2.rpkm.bw', sample=bamFiles)

rule bamCompare:
	input:
		b1 = bamDir + '{sample}.filtered.bam',
		b2 = ctrl
	output:
		'{sample}.log2.rpkm.bw'
	params:
		bl = blackList
	threads: 10
	shell:'''
	bamCompare -b1 {input.b1} -b2 {input.b2} -o {output} --operation log2 -bs 10 --blackListFileName {params.bl} --normalizeUsing RPKM --scaleFactorsMethod None
	'''
	
