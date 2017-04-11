#!/usr/bin/pypy
#from Bio.Seq import Seq
from time import gmtime, strftime
import sys
import os

if len(sys.argv) != 7:
	print("\nvcf_filter.py extract CDS (i.e: concatenated exons) from a genome, using informations")
	print("get from VCF (about coverage) and gff-file (positions of exons).")
	print("\nIt outputs one Fasta alignment per CDS.")
	print("\n\033[1;33m Example: ./vcf_filter.py Hmel2.fa Hmel2.gff subVCF_ama_0.vcf 3 corr_table_ama.txt ama\033[0m\n")
	print("\t\targ1 = name of the reference genome in Fasta (used to produce VCF + GFF files)")
	print("\t\targ2 = name of the GFF file, corresponding to the reference genome in Fasta")
	print("\t\targ3 = name of the \033[1;33msub-VCF file\033[0m produced after splitting the VCF using https://github.com/popgenomics/heliconius/blob/master/scripts/split_vcf.cpp")
	print("\t\targ4 = integer, corresponding to the minimum number of reads for calling ONE ALLELE at ONE GENE COPY")
	print("\t\targ5 = name of the correspondance table produced by https://github.com/popgenomics/heliconius/blob/master/scripts/split_vcf.cpp")
	print("\t\targ6 = single word describing the studied species")
	if(len(sys.argv)<7):
		sys.exit("\n\033[1;31m ERROR: 6 arguments are required: {0} missing\033[0m\n".format(7-len(sys.argv)))
	if(len(sys.argv)>7):
		sys.exit("\n\033[1;31m ERROR: 6 arguments are required: {0} too much\033[0m\n".format(len(sys.argv)-7))


genomeFileName = sys.argv[1]
gffFileName = sys.argv[2]
vcfFileName = sys.argv[3]
minCov = int(sys.argv[4])
corrTableFileName = sys.argv[5]
species = sys.argv[6]


testFile = os.path.isfile(genomeFileName)
if testFile == False:
	sys.exit("\n\033[1;31m ERROR: Reference genome file '{0}' is not found\033[0m\n".format(genomeFileName))


testFile = os.path.isfile(gffFileName)
if testFile == False:
	sys.exit("\n\033[1;31m ERROR: GFF file '{0}' is not found\033[0m\n".format(gffFileName))


testFile = os.path.isfile(vcfFileName)
if testFile == False:
	sys.exit("\n\033[1;31m ERROR: VCF file '{0}' is not found\033[0m\n".format(vcfFileName))


testFile = os.path.isfile(corrTableFileName)
if testFile == False:
	sys.exit("\n\033[1;31m ERROR: Correspondance table '{0}' is not found\033[0m\n".format(corrTableFileName))


#genomeFileName = "/home/roux/Documents/ABCheliconius/Hmel2.fa"
#gffFileName = "/home/roux/Documents/ABCheliconius/Hmel2.gff"
#vcfFileName = "/home/roux/Documents/ABCheliconius/test_ama10.Hmel2.bwa.default.HC.vcf"
#minCov = 3
#corrTableFileName = ""
#species = "mel"


# parse fasta
def fasta2dic(f):
	fasta = open(f).readlines()
	seqName = [x.split(" ")[0].rstrip().replace('>','') for x in fasta if x[0] == '>']
	seq = ''.join([x.rstrip() if x[0]!='>' else '@' for x in fasta])[1:].split('@')
	res = {}
	for i in range(len(seq)):
		res[seqName[i]] = seq[i]
	return (res)


def writeContig(contigName, alignment):
	outfile = open("{0}.fas".format(contigName), "w")
	for i in alignment:
		outfile.write(">{0}_allele1\n{1}\n".format(i, alignment[i]['allele1']))
		outfile.write(">{0}_allele2\n{1}\n".format(i, alignment[i]['allele2']))
	outfile.close()


def reverseComp(seq):
	seq2 = ""
	for i in seq[::-1]:
		if i in ["A", "T", "G", "C"]:
			if i == "A":
				seq2 += "T"
			if i == "T":
				seq2 += "A"
			if i == "G":
				seq2 += "C"
			if i == "C":
				seq2 += "G"
		else:
			seq2 += "N"
	return(seq2)


def extractGene(contigName, nGenes, alignment, gff, genome, species):
	for i in range(1, nGenes+1, 1):
		gene = gff[contigName + "_" + str(i)]
		if gene['nExons'] > 0:
			gene['min'].sort()
			gene['max'].sort()
			alignment2 = {}
			for j in alignment:
				if j not in alignment2:
					alignment2[j] = {}
				for allele in ["allele1", "allele2"]:
					cds = ""
					for k in range(gene['nExons']):
						cds = cds + alignment[j][allele][gene['min'][k]:(gene['max'][k]+1)]
					if gene['strand'] == '-':
						cds = reverseComp(cds)
					alignment2[j][allele] = cds
			res = ""
			for j in alignment2:
				nameA = contigName + "_" + str(i) + "|" + species + "|" + j + "|allele1" 
				nameB = contigName + "_" + str(i) + "|" + species + "|" + j + "|allele2" 
				res += ">{0}\n{1}\n>{2}\n{3}\n".format(nameA, alignment2[j]['allele1'], nameB, alignment2[j]['allele2'])
			outfile = open(species + "_" + contigName + "_" + str(i) + ".fas", "w")
			outfile.write(res)
			outfile.close()


def test_gene(x, genome):
	# x = entry of gff dictionnary
	# test the existance of codon stops
	x['min'].sort()
	x['max'].sort()
	contig = x['contig']
	gene = ""
	for i in range(x['nExons']):
		gene += genome[contig][x['min'][i]:(x['max'][i]+1)]
	if x['strand'] == '-':
		gene = reverseComp(gene)
	nStop = Seq(gene[:-3]).translate().count("*")
	print("{0} stop codons over {1} nucleotides".format(nStop, len(gene)))
	return(nStop)


#genomeFileName = sys.argv[1]
#gffFileName = sys.argv[2]
#vcfFileName = sys.argv[3]
#minCov = int(sys.argv[4])
#corrTableFileName = sys.argv[5]

surveyedContig = [] # list of surveyed contigs corresponding to the one present in the subVCF-file.
infile = open(corrTableFileName, "r")
for i in infile:
	i = i.strip().split("\t")
	if i[1] == vcfFileName:
		surveyedContig.append(i[0])
infile.close()


if len(surveyedContig) == 0:
	sys.exit("\n\033[1;31m ERROR: No contigs are present in {0}\033[0m\n".format(vcfFileName))


##############
bases = ['A', 'T', 'G', 'C']

### Genome ###
genome = fasta2dic(genomeFileName)


### GFF ###

gff = {} # contains all genes and informations about them
list_of_contigs = [] # contains the list of contigs 
nGenes_per_contigs = [] # contains the number of genes per contigs

geneID = 0
input = open(gffFileName, "r")
for i in input:
	if i[0:2] != "##":
		i = i.strip().split("\t")
		contig = i[0]
		if contig in surveyedContig:	
			if contig not in list_of_contigs:
				list_of_contigs.append(contig)
				nGenes_per_contigs.append(0)
				geneID = 0
			if i[2] == "gene":
				geneID += 1
				geneName = "{0}_{1}".format(contig, geneID)
			if geneName not in gff:
				nGenes_per_contigs[len(nGenes_per_contigs)-1] += 1
				gff[geneName] = {}
				gff[geneName]['min'] = []
				gff[geneName]['max'] = []
				gff[geneName]['strand'] = ''
				gff[geneName]['nExons'] = 0
			#if i[2] == "exon":
			if i[2] == "CDS":
				gff[geneName]['min'].append(int(i[3]) - 1) # convert gff unit (1-based) into python unit (0-based)
				gff[geneName]['max'].append(int(i[4]) - 1)
				gff[geneName]['strand'] = i[6]
				gff[geneName]['nExons'] += 1
				gff[geneName]['contig'] = contig
input.close()


# simply get the header and positions of columns
input = open(vcfFileName, "r")
for i in input:
	i = i.strip()
	if i[0:2] != "##":
		if i[0] == "#":
			i = i.split("\t")
			header = i
			pos_col = [ j for j in range(len(header)) if header[j] == 'POS' ][0]
			reference_col = [ j for j in range(len(header)) if header[j] == 'REF' ][0]
			alternative_col = [ j for j in range(len(header)) if header[j] == 'ALT' ][0]
			format_col = [ j for j in range(len(header)) if header[j] == 'FORMAT' ][0]
			break
input.close()


# get individuals present in the header
individuals = [] # list of headers
indiv_col = []
for i in range(len(header)):
	if i>8: # TO CHANGE IF REQUIRED: indicates columns beyond/after the FORMAT column in VCF file
		individuals.append(header[i])
		indiv_col.append(i)


# contigs found in vcf:
contigs_in_vcf = []
input = open(vcfFileName, "r")
for i in input:
	if i[0] != "#":
		i = i.strip().split("\t")
		if i[0] not in contigs_in_vcf:
			contigs_in_vcf.append(i[0])
input.close()
 

### LOOP OVER CONTIGS:
	### lOOP OVER CONTIG_GENES
#list_of_contigs = [] # contains the list of contigs 
#nGenes_per_contigs = [] # contains the number of genes per contigs

for i in range(len(list_of_contigs)):
	print("start of treatment of contig {0} over {1}: {2} {3})".format(i, len(list_of_contigs), list_of_contigs[i], strftime("%Y-%m-%d %H:%M:%S", gmtime())))
	contigName = list_of_contigs[i]
	if contigName in contigs_in_vcf:
		contig = genome[contigName]
		# prepare a matrix of sequences for individuals
		alignment = {}
		for j in individuals:
			alignment[j] = {}
			alignment[j]["allele1"] = contig 
			alignment[j]["allele2"] = contig 
		# get informations from vcf
		test = 0
		input = open(vcfFileName, "r")
		for j in input:
			if j[0]!='#':
				DP_presence = 0
				AD_presence = 0
				j = j.strip().split("\t")
				if j[0] != contigName and test==1: # stop reading vcf if all lines related to 'contigName' had been readen
					break
				if j[0] == contigName:
					test = 1
					pos = int(j[pos_col]) - 1 # converts the 'VCF' format (1-based) into 'python' format (0-based)
					ref = j[reference_col]
					alt = j[alternative_col]
					code = j[format_col].split(":") # informations contained per individual in columns
					if "DP" in code:
						DP = code.index("DP")
						DP_presence = 1
					if "AD" in code:
						AD = code.index("AD") 
						AD_presence = 1
					for k in range(len(indiv_col)): # loop over individuals
						ind_col_k = indiv_col[k] 
						ind_name_k = individuals[k]
						tmp = j[ind_col_k].split(":")
						if len(tmp) != len(code):
							alignment[individuals[k]]['allele1'] =  alignment[individuals[k]]['allele1'][:pos] + "N" + alignment[individuals[k]]['allele1'][(pos+1):]
							alignment[individuals[k]]['allele2'] =  alignment[individuals[k]]['allele2'][:pos] + "N" + alignment[individuals[k]]['allele2'][(pos+1):]
						if len(tmp) == len(code):
							################################
							# if ref has only one base, and if no alternative allele was found
							if len(ref)==1 and alt==".":
								if(tmp[DP] == "."):
									DP_k = 0
								else:
									DP_k = int(tmp[DP])
								if DP_k < (2 * minCov): # if not enough reads{both alleles are 'N'}else{keep the reference's allele}
									alignment[individuals[k]]['allele1'] =  alignment[individuals[k]]['allele1'][:pos] + "N" + alignment[individuals[k]]['allele1'][(pos+1):]
									alignment[individuals[k]]['allele2'] =  alignment[individuals[k]]['allele2'][:pos] + "N" + alignment[individuals[k]]['allele2'][(pos+1):]
							################################
							# if one ref has one base only, and if a A, T, G or C was found as an alternative allele
							if len(ref)==1 and ref in bases and alt in bases:
								AD_k = [ int(l) for l in tmp[AD].split(',') ]
								if tmp[DP] != ".":
									DP_k = int(tmp[DP])
								else:
									DP_k = sum(AD_k)
								if AD_k[0] < minCov and AD_k[1] < minCov:
									alignment[individuals[k]]['allele1'] =  alignment[individuals[k]]['allele1'][:pos] + "N" + alignment[individuals[k]]['allele1'][(pos+1):]
									alignment[individuals[k]]['allele2'] =  alignment[individuals[k]]['allele2'][:pos] + "N" + alignment[individuals[k]]['allele2'][(pos+1):]
#								if AD_k[0] >= minCov and AD_k[0] < (2 * minCov):
								if AD_k[0] >= minCov:
									alignment[individuals[k]]['allele1'] =  alignment[individuals[k]]['allele1'][:pos] + ref + alignment[individuals[k]]['allele1'][(pos+1):]
									if AD_k[1] < minCov:
										if AD_k[0] < (2 * minCov):
											alignment[individuals[k]]['allele2'] =  alignment[individuals[k]]['allele2'][:pos] + "N" + alignment[individuals[k]]['allele2'][(pos+1):]
										else:
											alignment[individuals[k]]['allele2'] =  alignment[individuals[k]]['allele2'][:pos] + ref + alignment[individuals[k]]['allele2'][(pos+1):]
									if AD_k[1] >= minCov:
										alignment[individuals[k]]['allele2'] =  alignment[individuals[k]]['allele2'][:pos] + alt + alignment[individuals[k]]['allele2'][(pos+1):]
								if AD_k[0] < minCov:
									if AD_k[1] >= minCov and AD_k[1] < (2 * minCov):
										alignment[individuals[k]]['allele1'] =  alignment[individuals[k]]['allele1'][:pos] + "N" + alignment[individuals[k]]['allele1'][(pos+1):]
										alignment[individuals[k]]['allele2'] =  alignment[individuals[k]]['allele2'][:pos] + alt + alignment[individuals[k]]['allele2'][(pos+1):]
									if AD_k[1] >= (2 * minCov):
										alignment[individuals[k]]['allele1'] =  alignment[individuals[k]]['allele1'][:pos] + alt + alignment[individuals[k]]['allele1'][(pos+1):]
										alignment[individuals[k]]['allele2'] =  alignment[individuals[k]]['allele2'][:pos] + alt + alignment[individuals[k]]['allele2'][(pos+1):]
		input.close()
		#writeContig(contigName, alignment # write the alignment of individuals for all contigs
		# deal with GFF file:
		nGenes = nGenes_per_contigs[i]
		extractGene(contigName, nGenes, alignment, gff, genome, species)
#		print("end of treatment of contig {0} over {1}: {2} {3})".format(i, len(list_of_contigs), list_of_contigs[i], strftime("%Y-%m-%d %H:%M:%S", gmtime())))





