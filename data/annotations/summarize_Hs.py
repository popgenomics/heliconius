#!/gepv/home2/croux/bin/pypy

from os import listdir

# populations
populations = ['ama', 'chi', 'flo', 'mal', 'ros', 'txn', 'vul', 'zel']
pos = ['ig', 'intron', 'miRNA', 'pos_1', 'pos_2', 'pos_3', 'pseudogenic_exon', 'rRNA', 'tRNA', 'gene']
rejected_contigs = ["Hmel219003", "Hmel201002", "Hmel200084", "Hmel206006"]

# compute the mean of a list 'x':
def mean(x):
	sum_tot = 0.0
	cnt = 0
	for i in x:
		cnt += 1
		sum_tot += i
	if cnt > 0:
		return(sum_tot/cnt)
	else:
		return(-9)

# treatFile = function getting informations in files
def treatFile(fileName, HS, populations):
	infile = open(fileName, "r")
	header = infile.readline().strip().split("\t")
	index_pos_FOR = header.index("pos_FOR")
	index_pos_REV = header.index("pos_REV")
	index_pop = {}
	# get column index for HS
	for i in populations:
		index_pop[i] = header.index("HS_{0}".format(i))

	for i in infile:
		i = i.strip().split("\t")
		pos_FOR = i[index_pos_FOR]	
		pos_REV = i[index_pos_REV]
		if pos_FOR == 'ig':
			for pop in populations:
				HS[pop][pos_REV].append(float(i[index_pop[pop]]))
		else:
			if pos_REV == 'ig':
				for pop in populations:
					HS[pop][pos_FOR].append(float(i[index_pop[pop]]))
	infile.close()

# test whether the contig is into a rejected list
def test_contig(FstFile, list_rejected):
        test = 0
        for i in list_rejected:
                if i in FstFile:
                        test += 1
        if test == 0:
                return(1)
        else:
                return(0)


def write_output(pos, populations, HS, chromoZ):
	# print the output
	output = "population"
	for i in pos:
		output += "\t{0}".format(i)
	output += "\n"
	
	for pop in populations:
		output += "{0}".format(pop)
		for pos_tmp in pos:
			output += "\t{0:.5f}".format(mean(HS[pop][pos_tmp]))
		output += "\n"

	outfile = open("table_polymorphism_{0}.txt".format(chromoZ), "w")
	outfile.write(output)
	outfile.close()


def initiate_HS(populations):
	HS = {}
	for i in populations:
		if i not in HS:
			HS[i] = {}
			# ig intron miRNA pos_1 pos_2 pos_3 pseudogenic_exon rRNA tRNA
			HS[i]['ig'] = []
			HS[i]['intron'] = []
			HS[i]['miRNA'] = []
			HS[i]['pos_1'] = []
			HS[i]['pos_2'] = []
			HS[i]['pos_3'] = []
			HS[i]['pseudogenic_exon'] = []
			HS[i]['rRNA'] = []
			HS[i]['tRNA'] = []
			HS[i]['gene'] = []
	return(HS)

# get the table corresponding to auto vs gonosomes
autosomes = [ i for i in listdir('./') if "Hmel" in i and "Hmel221" not in i ]
sexChro = [ i for i in listdir('./') if "Hmel" in i and "Hmel221" in i ]

# analyse the autosomes
HS = initiate_HS(populations)
for i in autosomes:
	test = test_contig(i, rejected_contigs)
	if test == 1:
                print(i)
		treatFile(i, HS, populations)
        else:
                print("reject: {0}".format(i))
write_output(pos, populations, HS, "autosome")


# analyse the sexChro 
HS = initiate_HS(populations)
for i in sexChro:
	test = test_contig(i, rejected_contigs)
	if test == 1:
                print(i)
		treatFile(i, HS, populations)
        else:
                print("reject: {0}".format(i))
write_output(pos, populations, HS, "sexChro")



