#!/usr/bin/pypy
import sys
import os

# check the arguments
if len(sys.argv) != 5:
	print("\n\tfasta2ABC.py produces a msLike output file from two alignments (one per species)")
	print("\n\033[1;33m\tExample: ./fasta2ABC.py ama_Hmel200001_1.fas chi_Hmel200001_1.fas 0.75 10\033[0m\n")
	print("\t\targ1 =\tname of the fasta file containing the alignment for species A")
	print("\t\targ2 =\tname of the fasta file containing the alignment for species B")
	print("\t\targ3 =\tvalue in [0-1]. Corresponds to a threshold of %N above which a sequence is rejected")
	print("\t\targ4 =\tinteger, corresponding to the minimum number of retained sequences (after rejection).\n\t\t\tif not enough sequences are retained, the loci is excluded from the analysis")
	if(len(sys.argv)<5):
		sys.exit("\n\033[1;31m ERROR: 4 arguments are required: {0} missing\033[0m\n".format(5-len(sys.argv)))
	if(len(sys.argv)>5):
		sys.exit("\n\033[1;31m ERROR: 4 arguments are required: {0} too much\033[0m\n".format(len(sys.argv)-5))

seqA = sys.argv[1] # alignment of sequences from species A
seqB = sys.argv[2] # alignment of sequences from species B
threshold_N = float(sys.argv[3]) # if an allele has %N > threshold_N --> sequence is rejected
nMin = int(sys.argv[4]) # minimum number of individuals within a species


test = os.path.isfile(seqA)
if test == False:
	sys.exit("\n\033[1;31m ERROR: alignement '{0}' of species A is not found\033[0m\n".format(seqA))

test = os.path.isfile(seqB)
if test == False:
	sys.exit("\n\033[1;31m ERROR: alignement '{0}' of species B is not found\033[0m\n".format(seqB))


def coloredSeq(seq):
	# print sequences with the standard color code
	seq = seq.replace("A", '\x1b[5;30;41m' + 'A' + '\x1b[0m')
	seq = seq.replace("T", '\x1b[5;30;44m' + 'T' + '\x1b[0m')
	seq = seq.replace("G", '\x1b[6;30;43m' + 'G' + '\x1b[0m')
	seq = seq.replace("C", '\x1b[5;30;42m' + 'C' + '\x1b[0m')
	return(seq)


def trunc2triplets(size):
	# trunc a value to its closest and smaller multiple of 3
	size = int(size)
	for i in range(2):
		if size%3 != 0:
			size -= 1
	return(size)


# nN = number of non-synonymous sites in the codon i: example for CGG -> nN = 2/3 + 3/3 + 0/3
# nS = number of synonymous sites in the codon i: example for CGG -> n> = 1/3 + 0/3 + 3/3
codonTable = {'AAA': {'aa': 'K', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAC': {'aa': 'N', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAG': {'aa': 'K', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAT': {'aa': 'N', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'ACA': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACC': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACG': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACT': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'AGA': {'aa': 'R', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'AGC': {'aa': 'S', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AGG': {'aa': 'R', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'AGT': {'aa': 'S', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'ATA': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'ATC': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'ATG': {'aa': 'M', 'nN': 3.0, 'nS': 0.0}, 'ATT': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'CAA': {'aa': 'Q', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAC': {'aa': 'H', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAG': {'aa': 'Q', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAT': {'aa': 'H', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CCA': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCC': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCG': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCT': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CGA': {'aa': 'R', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CGC': {'aa': 'R', 'nN': 2.0, 'nS': 1.0}, 'CGG': {'aa': 'R', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CGT': {'aa': 'R', 'nN': 2.0, 'nS': 1.0}, 'CTA': {'aa': 'L', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CTC': {'aa': 'L', 'nN': 2.0, 'nS': 1.0}, 'CTG': {'aa': 'L', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CTT': {'aa': 'L', 'nN': 2.0, 'nS': 1.0}, 'GAA': {'aa': 'E', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAC': {'aa': 'D', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAG': {'aa': 'E', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAT': {'aa': 'D', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GCA': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCC': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCG': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCT': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GGA': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGC': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGG': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGT': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GTA': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTC': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTG': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTT': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'TAC': {'aa': 'Y', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TAT': {'aa': 'Y', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TCA': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCC': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCG': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCT': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TGC': {'aa': 'C', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TGG': {'aa': 'W', 'nN': 3.0, 'nS': 0.0}, 'TGT': {'aa': 'C', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TTA': {'aa': 'L', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'TTC': {'aa': 'F', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TTG': {'aa': 'L', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'TTT': {'aa': 'F', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}}


def fasta2dic(fastaFile):
	fasta = open(fastaFile).readlines()
	seqName = [x.split(" ")[0].rstrip().replace('>','') for x in fasta if x[0] == '>']
	seq = ''.join([x.rstrip() if x[0]!='>' else '@' for x in fasta])[1:].split('@')
	res = {}
	for i in range(len(seq)):
		res[seqName[i]] = seq[i]
	return (res)



alignA = fasta2dic(seqA)
alignB = fasta2dic(seqB)


L = len(alignA[alignA.keys()[0]]) - 3 # remove the last 3 bases to excluse final stop codon
L = trunc2triplets(L) # convert the remaining length into a multiple of 3


interspe = []
nA = 0
nB = 0
for i in alignA:
	seq = alignA[i][0:L]
	propN = seq.count("N")/(1.0 * L)
	if propN < threshold_N:
		nA += 1
		interspe.append(seq)

for i in alignB:
	seq = alignB[i][0:L]
	propN = seq.count("N")/(1.0 * L)
	if propN < threshold_N:
		nB += 1
		interspe.append(seq)


if nA < nMin:
	if nB >= nMin:
		sys.exit("Species A has less than {0} sequences. Species B is ok".format(nMin))
	if nB < nMin:
		sys.exit("Species A and B have less than {0} sequences".format(nMin))
if nB < nMin:
	sys.exit("Species B have less than {0} sequences. Species A is ok". format(nMin))


#for i in interspe:
#	print(">{0}\n{1}".format(i, coloredSeq(interspe[i][0:120])))



nSites = 0 # total number of synonymous sites within the sequence, computed using codonTable
nSynSegSite = 0 # number of synonymous segregating sites among the nSites
positions = [] # list of synonymous polymorphic positions: doesn't correspond to the SNP position, but to the first codon position
msStyle = [] # contains the msStyle format
for ind in range(nA):
	msStyle.append([])
for ind in range(nB):
	msStyle.append([])

# loop over codons:

for pos in range(L)[::3]:
	alignmentOfCodons = [] # set of codons in the alignment, starting at the position 'pos1'
	# loop over individuals:
	# get all codons in the alignment
	for ind in range(nA + nB):
		pos1 = interspe[ind][pos]
		pos2 = interspe[ind][pos + 1]
		pos3 = interspe[ind][pos + 2]
		base = pos1 + pos2 + pos3 
		alignmentOfCodons.append(base)
	
	polyMcodons = list(set(alignmentOfCodons)) # list of codons found in the alignment
	nCodons = 0
	nCodons = len(polyMcodons)
	testN = False # False if no codon with 'N'; True if a 'N' is found in at least one codon for one individual
	testStopCodon = False # False if no stop codon was found; True if a stop codon was found
	for i in polyMcodons: # loop to test for some 'N'
		if 'N' in i:
			testN = True
		if i not in codonTable:
			testStopCodon = True
	
	# if: 1) a maximum of 2 polymorphic codons, and, 2) no codon with 'N', and, 3) all codons effectively code for an amino acid
	if nCodons <= 2 and testN==False and testStopCodon==False: 
		nSites_pos = 0.0
		for i in alignmentOfCodons:
			nSites_pos += codonTable[i]['nS']
		nSites += nSites_pos/len(alignmentOfCodons)
		
		# if two codons --> there is a polymorphism
		if nCodons == 2:
			alignmentOfAminoAcids = []
			for i in alignmentOfCodons:
				alignmentOfAminoAcids.append(codonTable[i]['aa'])
			setOfAminoAcids = list(set(alignmentOfAminoAcids))
			
			# if two codons but one amino acids --> synonymous polymorphism
			if len(setOfAminoAcids) == 1:
				nSynSegSite += 1
				positions.append(pos) # positions: list of first codon position of polymorphic synonymous codons
				ancestralAllele = polyMcodons[0] # in absence of outgroup --> the ancestral allele is the first in the alignement
				derivedAllele = polyMcodons[1] # without outgroup --> the derived allele is the one who is not the first...
				for i in range(nA + nB):
					if alignmentOfCodons[i] == ancestralAllele:
						msStyle[i].append('0')
					if alignmentOfCodons[i] == derivedAllele:
						msStyle[i].append('1')



# print some global informations
#print("# Length = {0}".format(L))
#print("# of synonymous positions = {0}".format(nSites))
#print("\t".join( [ str(i) for i in positions ] ))

# core of the output files names
geneName = seqA
if "/" in geneName:
	geneName.split("/")[-1:]
geneName = geneName.split(".")[0][4:]


geneName = seqA.split(".")[0][4:]


# ms_like output files
locus_ms = "./msnsam tbs 20 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ej tbs 1 2 -eN tbs tbs\n3579 27011 59243\n\n"
locus_ms = locus_ms + "//" + "\n"
locus_ms = locus_ms + "segsites: {0}\n".format(int(nSynSegSite))
if nSynSegSite != 0:
	locus_ms += "positions: {0}\n".format( " ".join([ str(round((1.0*i)/L, 4)) for i in positions ]))
	for i in msStyle:
		locus_ms = locus_ms + "".join( [ str(j) for j in i ] ) + "\n"

outfile = open(geneName + ".ms", "w")
outfile.write(locus_ms)
outfile.close()

# informations about locus
res = "locusName\tL_including_N\tLsyno\tnSynSegSite\tnsamA\tnsamB\n"
res += "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(geneName, L, nSites, nSynSegSite, nA, nB)

outfile = open(geneName + "_info.txt", "w")
outfile.write(res)
outfile.close()
 
