#!/usr/bin/python

from Bio.SeqIO import parse
import sys
import os


if len(sys.argv) != 2:
	print("\n\tPhases a fasta alignment into a new fasta file")
	print("\n\tCalls fastPHASE from Paul Scheet and Matthew Stephens (http://scheet.org/software.html)\n")
	print("\033[1;33m\tExample: ./fasta2phase.py ama_Hmel200102_1.fas \033[0m")
	sys.exit("\n\033[1;31m\tNeeds one argument: input fasta file\033[0m\n")


fastaFileName = sys.argv[1]


testFile = os.path.isfile(fastaFileName)
if testFile == False:
	sys.exit("\n\033[1;31m ERROR: The fasta alignment '{0}' is not found\033[0m\n".format(fastaFileName))


# get the gene's name
geneName = fastaFileName.split("/")[-1:]
geneName = geneName[0].split(".")[0]


# seq = dictionary; seq['individual']['allele1', 'allele2'] = sequence
seq = {}
individuals = []
infile = parse(fastaFileName, "fasta")


fullNames = []
for i in infile:
	fullNames.append(i.id)
	ind = i.id.split("|")[2]
	if ind not in seq:
		individuals.append(ind)
		seq[ind] = {}
	if "allele1" in i.id:
		seq[ind]['allele1'] = str(i.seq)
		seq[ind]['name_allele1'] = i.id
	if "allele2" in i.id:
		seq[ind]['allele2'] = str(i.seq)
		seq[ind]['name_allele2'] = i.id


L = len(i.seq)


# seq2 = fictionary; seq2['individual']['allele1', 'allele2'] = "polymorphic alleles"
seq2 = {}
for i in individuals:
	if i not in seq2:
		seq2[i] = {}
	seq2[i]['allele1'] = ""
	seq2[i]['allele2'] = ""


nBiallSite = 0
positions = []
for pos in range(L):
	alignment = []
	for ind in individuals:
		alignment.append(seq[ind]['allele1'][pos])
		alignment.append(seq[ind]['allele2'][pos])
	alignment_noN = [ i for i in alignment if i!='N' ]
	setOfBases = list(set(alignment_noN))
	if len(setOfBases) == 2:
		nBiallSite += 1
		positions.append(pos)
		for ind in individuals:
			allele1 = seq[ind]['allele1'][pos]
			if allele1 == 'N':
				allele1='?'
			allele2 = seq[ind]['allele2'][pos]
			if allele2 == 'N':
				allele2='?'
			seq2[ind]['allele1'] += allele1
			seq2[ind]['allele2'] += allele2



# If there is at least 2 polymorphic positions
if len(positions) >= 2:
	res = "{0}\n{1}\n".format(len(individuals), nBiallSite)
	res += "P {0}\n".format(" ".join([ str(i) for i in positions] ))
	for ind in individuals:
		res += "# {0}\n".format(ind)
		res += "{0}\n{1}\n".format(seq2[ind]['allele1'], seq2[ind]['allele2'])


	fastPhase_input = geneName + "_fastPhase.inp"
	outfile = open(fastPhase_input, "w")
	outfile.write(res)
	outfile.close()


	# run fastPHASE
	commandLine = "fastPHASE -o{0} -T10 -KL6 -KU12 -Ki2 -Ks50 -Km1000 -Kp.05 {1} >/dev/null 2>&1".format("outputPhase_" + geneName, fastPhase_input)
	testCommand = os.system(commandLine)
	if testCommand != 0:
		sys.exit("\n\033[1;31m\tERROR: fastPHASE failed of being executed for locus {0}\033[0m\n".format(geneName))


	# rm outputPhase_ama_Hmel200102_1_origchars
	testCommand = os.system("rm {0}".format(fastPhase_input))
	if testCommand != 0:
		sys.exit("\n\033[1;31m\tERROR: cannot delete the file {0}\033[0m\n".format(fastPhase_input))


	# rm outputPhase_ama_Hmel200102_1_origchars
	targetFile = "outputPhase_" + geneName + "_origchars"
	testCommand = os.system("rm {0}".format(targetFile))
	if testCommand != 0:
		sys.exit("\n\033[1;31m\tERROR: cannot delete the file {0}\033[0m\n".format(targetFile))


	# rm outputPhase_ama_Hmel200102_1_finallikelihoods
	targetFile = "outputPhase_" + geneName + "_finallikelihoods"
	testCommand = os.system("rm {0}".format(targetFile))
	if testCommand != 0:
		sys.exit("\n\033[1;31m\tERROR: cannot delete the file {0}\033[0m\n".format(targetFile))


	# get fastPHASE's output file
	outputFileName = "outputPhase_" + geneName + "_hapguess_switch.out"

	output = {}
	test = 0
	infile = open(outputFileName, "r")
	for i in infile:
		i = i.strip()
		if i == "BEGIN GENOTYPES":
			test = 1
			continue
		if test == 1 and i != "END GENOTYPES":
			if i[0] == "#":
				allele = 0
				gene = i.split(" ")[1]
				output[gene] = {}
				continue
			if i[0] != "#":
				allele += 1
				output[gene]['allele' + str(allele) ] = i.replace(' ', '')
	infile.close()

	testCommand = os.system("rm {0}".format(outputFileName))
	if testCommand != 0:
		sys.exit("\n\033[1;31m\tERROR: cannot delete the file {0}\033[0m\n".format(outputFileName))


	# put 'N' in data
	for ind in seq2:
		for allele in ['allele1', 'allele2']:
			for pos in range(len(seq2[ind][allele])):
				base = seq2[ind][allele][pos]
				if base == '?':
					output[ind]['allele1'] = output[ind]['allele1'][:pos] + 'N' + output[ind]['allele1'][(pos+1):]
					output[ind]['allele2'] = output[ind]['allele2'][:pos] + 'N' + output[ind]['allele2'][(pos+1):]
		for i in positions:
			seq[ind]['allele1'] = seq[ind]['allele1'][:pos] + output[ind]['allele1'][pos]  + seq[ind]['allele1'][(pos+1):]
			seq[ind]['allele2'] = seq[ind]['allele2'][:pos] + output[ind]['allele2'][pos]  + seq[ind]['allele2'][(pos+1):]
# END OF BLOCK: "If there is at least 2 polymorphic positions 

res = ""
for i in seq:
	res += ">{0}\n{1}\n".format(seq[i]['name_allele1'], seq[i]['allele1'])
	res += ">{0}\n{1}\n".format(seq[i]['name_allele2'], seq[i]['allele2'])


outfile = open("PHASED_{0}.fas".format(geneName), "w")
outfile.write(res)
outfile.close()

