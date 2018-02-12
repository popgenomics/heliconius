#!/home/roux/bin/pypy
# to run on the 'nereis' server: /gepv/home2/croux/heliconius/data

import sys

contig = sys.argv[1]

def parseDNA(x):
	res = {}
	names = []
	seq = []
	infile = open(x, "r")
	for i in infile:
		if i[0] == ">":
			names.append(i.strip()[1:])
		else:
			seq.append(i.strip())
	infile.close()
	res['id'] = names
	res['seq'] = seq
	return(res)

def findSNP(ama, chi, flo, mal, ros, txn ,vul, zel, pos):
	# returns a string with: refA / altA / {nRef_allele / nAlt_allele / nN} x {ama, chi, flo, mal, ros, txn, vul, zel}
	amapos = []
	chipos = []
	flopos = []
	malpos = []
	rospos = []
	txnpos = []
	vulpos = []
	zelpos = []
	tot = []
	for ind in ama:
		tmp = ind[pos]
		amapos.append(tmp)
		tot.append(tmp)
	for ind in chi:
		tmp = ind[pos]
		chipos.append(tmp)
		tot.append(tmp)
	for ind in flo:
		tmp = ind[pos]
		flopos.append(tmp)
		tot.append(tmp)
	for ind in mal:
		tmp = ind[pos]
		malpos.append(tmp)
		tot.append(tmp)
	for ind in ros:
		tmp = ind[pos]
		rospos.append(tmp)
		tot.append(tmp)
	for ind in txn:
		tmp = ind[pos]
		txnpos.append(tmp)
		tot.append(tmp)
	for ind in vul:
		tmp = ind[pos]
		vulpos.append(tmp)
		tot.append(tmp)
	for ind in zel:
		tmp = ind[pos]
		zelpos.append(tmp)
		tot.append(tmp)
	alleles = list(set(tot))
	if 'N' in alleles:
		alleles.pop(alleles.index('N'))
	nAlleles = len(alleles)
	
	if nAlleles == 2:
		refA = alleles[0]
		altA = alleles[1]
		res = "{0}\t{1}\t".format(refA, altA)
		res += "{0}\t{1}\t{2}\t".format(amapos.count(refA), amapos.count(altA), amapos.count('N'))
		res += "{0}\t{1}\t{2}\t".format(chipos.count(refA), chipos.count(altA), chipos.count('N'))
		res += "{0}\t{1}\t{2}\t".format(flopos.count(refA), flopos.count(altA), flopos.count('N'))
		res += "{0}\t{1}\t{2}\t".format(malpos.count(refA), malpos.count(altA), malpos.count('N'))
		res += "{0}\t{1}\t{2}\t".format(rospos.count(refA), rospos.count(altA), rospos.count('N'))
		res += "{0}\t{1}\t{2}\t".format(txnpos.count(refA), txnpos.count(altA), txnpos.count('N'))
		res += "{0}\t{1}\t{2}\t".format(vulpos.count(refA), vulpos.count(altA), vulpos.count('N'))
		res += "{0}\t{1}\t{2}\n".format(zelpos.count(refA), zelpos.count(altA), zelpos.count('N'))
		return(res)
	else:
		if nAlleles == 0:
			res = "N\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n"
			return(res)
                if nAlleles == 1:
                        refA = alleles[0]
                        res = "{0}\t{1}\t".format(alleles[0], 'N')
                        res += "{0}\t{1}\t{2}\t".format(amapos.count(refA), 0, amapos.count('N'))
                        res += "{0}\t{1}\t{2}\t".format(chipos.count(refA), 0, chipos.count('N'))
                        res += "{0}\t{1}\t{2}\t".format(flopos.count(refA), 0, flopos.count('N'))
                        res += "{0}\t{1}\t{2}\t".format(malpos.count(refA), 0, malpos.count('N'))
                        res += "{0}\t{1}\t{2}\t".format(rospos.count(refA), 0, rospos.count('N'))
                        res += "{0}\t{1}\t{2}\t".format(txnpos.count(refA), 0, txnpos.count('N'))
                        res += "{0}\t{1}\t{2}\t".format(vulpos.count(refA), 0, vulpos.count('N'))
                        res += "{0}\t{1}\t{2}\n".format(zelpos.count(refA), 0, zelpos.count('N'))
                        return(res)
		if nAlleles > 2:
			res = "multiallelic\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n"
			return(res)


ama = parseDNA("/home/roux/scratch/heliconius/data/ama/1_contigs/contig_{0}_ama.fasta".format(contig))['seq']
chi = parseDNA("/home/roux/scratch/heliconius/data/chi/1_contigs/contig_{0}_chi.fasta".format(contig))['seq']
flo = parseDNA("/home/roux/scratch/heliconius/data/flo/1_contigs/contig_{0}_flo.fasta".format(contig))['seq']
mal = parseDNA("/home/roux/scratch/heliconius/data/mal/1_contigs/contig_{0}_mal.fasta".format(contig))['seq']
ros = parseDNA("/home/roux/scratch/heliconius/data/ros/1_contigs/contig_{0}_ros.fasta".format(contig))['seq']
txn = parseDNA("/home/roux/scratch/heliconius/data/txn/1_contigs/contig_{0}_txn.fasta".format(contig))['seq']
vul = parseDNA("/home/roux/scratch/heliconius/data/vul/1_contigs/contig_{0}_vul.fasta".format(contig))['seq']
zel = parseDNA("/home/roux/scratch/heliconius/data/zel/1_contigs/contig_{0}_zel.fasta".format(contig))['seq']


# read the annotation file producted from the gff
annoFile = open('/home/roux/scratch/heliconius/data/annotations/annotations_{0}.txt'.format(contig), 'r')

res = annoFile.readline().strip() + "\trefA\taltA"
for i in ['ama', 'chi', 'flo', 'mal', 'ros', 'txn', 'vul', 'zel']:
	res += "\tnRefA_{0}\tnAltA_{0}\tnN_{0}".format(i)
res += "\n"

# write the output
outfile = open("{0}_freq.txt".format(contig), "w")
outfile.write(res)
outfile.close()

outfile = open("{0}_freq.txt".format(contig), "a")
cnt = 0
for i in annoFile:
	res = i.strip() + "\t" + findSNP(ama, chi, flo, mal, ros, txn ,vul, zel, cnt)
	outfile.write(res)
	cnt += 1

annoFile.close()
outfile.close()

