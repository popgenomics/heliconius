#!/gepv/home2/croux/bin/pypy

from os import listdir

# populations
populations = ['ama', 'chi', 'flo', 'mal', 'ros', 'txn', 'vul', 'zel'] # list of populations
pos = ['ig', 'intron', 'miRNA', 'pos_1', 'pos_2', 'pos_3', 'pseudogenic_exon', 'rRNA', 'tRNA', 'gene'] # list of positions
stats = ['same', 'Sf', 'Ss', 'SxA', 'SxB', 'FST'] # list of statistics

# compute the mean of a list 'x':
def mean(x):
	sum_tot = 0.0
	cnt = 0
	for i in x:
		if i != -9:
			cnt += 1
			sum_tot += i
	if cnt > 0:
		return(sum_tot/cnt)
	else:
		return(-9)

# treatFile = function getting informations in files
def treatFile(fileName, div, populations, stats, pos):
	# div div[pair][pos: ig, intron, ... ][stat: FST, sxA, sxB, same, ... ]
	infile = open(fileName, "r")
	header = infile.readline().strip().split("\t")
	index_pos_FOR = header.index("pos_FOR")
	index_pos_REV = header.index("pos_REV")
	index_pair = {}
	pairs = []
	# get column index for div
	for pop_A in range(len(populations)-1):
		for pop_B in range(pop_A+1, len(populations), 1):
			pair = "{0}_{1}".format(populations[pop_A], populations[pop_B])
			pairs.append(pair)
			index_pair[pair] = {} # index_pair[pair][position]['FST' or 'site']
			index_pair[pair] = {}
			index_pair[pair]['FST'] = header.index("FST_{0}".format(pair))
			index_pair[pair]['site'] = header.index("site_{0}".format(pair))

	for i in infile:
		i = i.strip().split("\t")
		pos_FOR = i[index_pos_FOR]	
		pos_REV = i[index_pos_REV]
		if pos_FOR == 'ig':
			for pair in pairs:
				fst = float(i[index_pair[pair]['FST']])
				if fst != -9:
					div[pair][pos_REV]['FST'].append(fst) # append the FST
				div[pair][pos_REV][i[index_pair[pair]['site']]] += 1 # add +1 to the site
		else:
			if pos_REV == 'ig':
				for pair in pairs:
					fst = float(i[index_pair[pair]['FST']])
					if fst != -9:
						div[pair][pos_FOR]['FST'].append(fst) # append the FST
					div[pair][pos_FOR][i[index_pair[pair]['site']]] += 1 # add +1 to the site
	infile.close()


# get the table corresponding to auto vs gonosomes
autosomes = [ i for i in listdir('./') if "Hmel" in i and "Hmel221" not in i ]
sexChro = [ i for i in listdir('./') if "Hmel" in i and "Hmel221" in i ]

# store the results
def initiate_div(populations, pos, stats):
	div = {} # div[pair][pos: ig, intron, ... ][stat: FST, sxA, sxB, same, ... ]
	for pop_A in range(len(populations)-1):
		for pop_B in range(pop_A+1, len(populations), 1):
			pair = "{0}_{1}".format(populations[pop_A], populations[pop_B])
			if pair not in div:
				div[pair] = {}
				# ig intron miRNA pos_1 pos_2 pos_3 pseudogenic_exon rRNA tRNA
				for pos_i in pos:
					div[pair][pos_i] = {}
					for stat_i in stats:
						if stat_i == 'FST':
							div[pair][pos_i][stat_i] = [] # to compute the mean
						else:
							div[pair][pos_i][stat_i] = 0 # to compute the total number
	return(div)


def write_output(stats, pos, populations, div, chromoZ, contig):
	for stat_i in stats:
		output = "population"
		for pos_i in pos:
			output += "\t{0}".format(pos_i)
		output += "\n"
		
		for pop_A in range(len(populations)-1):
			for pop_B in range(pop_A+1, len(populations), 1):
				pair = "{0}_{1}".format(populations[pop_A], populations[pop_B])
				output += pair
				for pos_i in pos:
					if stat_i == "FST":
						output += "\t{0:.5f}".format(mean(div[pair][pos_i][stat_i]))
					else:
						output += "\t{0}".format(div[pair][pos_i][stat_i])
				output += "\n"
		
		# write into a file
#		outfile = open("table_divergence_autosomes_{0}_{1}.txt".format(stat_i, contig), "w")
		outfile = open("table_divergence_{0}_{1}_{2}.txt".format(chromoZ, stat_i, contig), "w")
		outfile.write(output)
		outfile.close()

# read the files:
#for i in autosomes:
#for i in sexChro:
#i="annot_Hmel207002_Fst.txt"
#i="annot_Hmel218018_Fst.txt"
#autosomes = ["annot_Hmel217001_Fst.txt", "annot_Hmel219014_Fst.txt", "annot_Hmel220003_Fst.txt", "annot_Hmel201011_Fst.txt", "annot_Hmel212013_Fst.txt", "annot_Hmel209007_Fst.txt", "annot_Hmel206009_Fst.txt", "annot_Hmel211004_Fst.txt", "annot_Hmel217004_Fst.txt", "annot_Hmel202006_Fst.txt", "annot_Hmel210011_Fst.txt", "annot_Hmel220005_Fst.txt", "annot_Hmel212001_Fst.txt", "annot_Hmel207002_Fst.txt", "annot_Hmel219003_Fst.txt", "annot_Hmel210004_Fst.txt", "annot_Hmel201009_Fst.txt", "annot_Hmel218018_Fst.txt"]

print("\nSexChro:")
for i in sexChro:
	div = initiate_div(populations, pos, stats)
	treatFile(i, div, populations, stats, pos)
	write_output(stats, pos, populations, div, "sexChro", i.split("_")[1])
	print(i)


print("\n\nAutosomes:")
for i in autosomes:
	div = initiate_div(populations, pos, stats)
	treatFile(i, div, populations, stats, pos)
	write_output(stats, pos, populations, div, "autosome", i.split("_")[1])
	print(i)

