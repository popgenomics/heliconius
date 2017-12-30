#!/usr/bin/python
# -*-coding:Latin-1 -*

from numpy.random import uniform
from numpy.random import beta
from numpy.random import choice 
# uniform(min, max, nSample)
# beta(a, b, nSample)

# msms line:
# msms nsam nIter*nLoci -s tbs -r tbs tbs
# -n 1 10 -n 2 tbs -n 3 tbs -n 4 tbs -n 5 tbs -n 6 tbs -n 7 tbs -n 8 tbs
# -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 1 4 tbs -m 4 1 tbs -m 1 5 tbs -m 5 1 tbs -m 1 6 tbs -m 6 1 tbs -m 1 7 tbs -m 7 1 tbs -m 1 8 tbs -m 8 1 tbs
# -m 2 3 tbs -m 3 2 tbs -m 2 4 tbs -m 4 2 tbs -m 2 5 tbs -m 5 2 tbs -m 2 6 tbs -m 6 2 tbs -m 2 7 tbs -m 7 2 tbs -m 2 8 tbs -m 8 2 tbs
# -m 3 4 tbs -m 4 3 tbs -m 3 5 tbs -m 5 3 tbs -m 3 6 tbs -m 6 3 tbs -m 3 7 tbs -m 7 3 tbs -m 3 8 tbs -m 8 3 tbs
# -m 4 5 tbs -m 5 4 tbs -m 4 6 tbs -m 6 4 tbs -m 4 7 tbs -m 7 4 tbs -m 4 8 tbs -m 8 4 tbs
# -m 5 6 tbs -m 6 5 tbs -m 5 7 tbs -m 7 5 tbs -m 5 8 tbs -m 8 5 tbs
# -m 6 7 tbs -m 7 6 tbs -m 6 8 tbs -m 8 6 tbs
# -m 7 8 tbs -m 8 7 tbs


# prior distribution
Tau_prior = [0, 2.0]
tau_scale = [1, 10]

N_scale = [0, 2]
ref_N = 10
shape_beta = [0, 10]

M_prior = [0, 0.1]
list_M = [0, 0.1, 10]

# n multilocus iterations
nIter = 10

# read informations about loci from "bpfile"
infile = open("bpfile", "r") # nSNPs L n1 n2 n3 n4 n5 n6 n7 n8 rho
tmp = infile.readline() # header
tmp = infile.readline().strip().split("\t")
nSNPs = [ int(i) for i in tmp ]

tmp = infile.readline().strip().split("\t")
L = [ int(i) for i in tmp ]

tmp = infile.readline().strip().split("\t")
n1sam = [ int(i) for i in tmp ]

tmp = infile.readline().strip().split("\t")
n2sam = [ int(i) for i in tmp ]

tmp = infile.readline().strip().split("\t")
n3sam = [ int(i) for i in tmp ]

tmp = infile.readline().strip().split("\t")
n4sam = [ int(i) for i in tmp ]

tmp = infile.readline().strip().split("\t")
n5sam = [ int(i) for i in tmp ]

tmp = infile.readline().strip().split("\t")
n6sam = [ int(i) for i in tmp ]

tmp = infile.readline().strip().split("\t")
n7sam = [ int(i) for i in tmp ]

tmp = infile.readline().strip().split("\t")
n8sam = [ int(i) for i in tmp ]

tmp = infile.readline().strip().split("\t")
rho = [ float(i) for i in tmp ]

infile.close()

nLoci = len(nSNPs)

nSamTot_multilocus = [ n1sam[i]+n2sam[i]+n3sam[i]+n4sam[i]+n5sam[i]+n6sam[i]+n7sam[i]+n8sam[i] for i in range(nLoci) ] # total number of samples for each of the nLoci loci

# labels migration matrix
migration_matrix_labels = "" 
for i in range(1,8,1):
	for j in range(i+1,9,1):
		migration_matrix_labels += "m_{0}_{1}\tm_{1}_{0}\t".format(i,j)

# priorfile 
priorfile = "N1\tN2\tN3\tN4\tN5\tN6\tN7\tN8\t"
priorfile += "T_1_2\tT_3_4\tT_5_6\tT_7_8\t"
priorfile += "T_12_34\tT_56_78\t"
priorfile += "T_1234_5678\t"
priorfile += "{0}\n".format(migration_matrix_labels.strip())


# loop over multilocus iterations
for iter_i in range(nIter):
	# TIMES
	tau_1 = uniform(Tau_prior[0], Tau_prior[1]) # between 1 and 2
	tau_2 = uniform(Tau_prior[0], Tau_prior[1]) # between 3 and 4
	tau_3 = uniform(Tau_prior[0], Tau_prior[1]) # between 5 and 6
	tau_4 = uniform(Tau_prior[0], Tau_prior[1]) # between 7 and 8
	tau_5 = uniform(tau_scale[0], tau_scale[1]) * max([tau_1, tau_2]) # between 1_2 and 3_4
	tau_6 = uniform(tau_scale[0], tau_scale[1]) * max([tau_3, tau_4]) # between 5_6 and 7_8
	tau_7 = uniform(tau_scale[0], tau_scale[1]) * max([tau_5, tau_6]) # between 1234 and 5678


	# EFFECTIVE SIZES	
	shape1_N = uniform(shape_beta[0], shape_beta[1])
	shape2_N = uniform(shape_beta[0], shape_beta[1])
	
	bf = beta(shape1_N, shape2_N, nLoci)
	
	# current populations
	N1 = ref_N # pop 1
	N2 = ref_N * uniform(N_scale[0], N_scale[1]) # pop 2
	N3 = ref_N * uniform(N_scale[0], N_scale[1]) # pop 3
	N4 = ref_N * uniform(N_scale[0], N_scale[1]) # pop 4
	N5 = ref_N * uniform(N_scale[0], N_scale[1]) # pop 5
	N6 = ref_N * uniform(N_scale[0], N_scale[1]) # pop 6
	N7 = ref_N * uniform(N_scale[0], N_scale[1]) # pop 7
	N8 = ref_N * uniform(N_scale[0], N_scale[1]) # pop 8
	
	N1_multilocus = [ ref_N * i for i in bf ] # pop 1
	N2_multilocus = [ N2 * i for i in bf ] # pop 2	
	N3_multilocus = [ N3 * i for i in bf ] # pop 3	
	N4_multilocus = [ N4 * i for i in bf ] # pop 4	
	N5_multilocus = [ N5 * i for i in bf ] # pop 5	
	N6_multilocus = [ N6 * i for i in bf ] # pop 6	
	N7_multilocus = [ N7 * i for i in bf ] # pop 7	
	N8_multilocus = [ N8 * i for i in bf ] # pop 8	
	
	# ancestral populations
	N1a = (N1+N2)/2.0 # ancestor between 1 and 2
	N3a = (N3+N4)/2.0 # ancestor between 3 and 4
	N1b = (N1a+N3a)/2.0 # ancestor between 12 and 34
	
	N5a = (N5+N6)/2.0 # ancestor between 1 and 2
	N7a = (N7+N8)/2.0 # ancestor between 3 and 4
	N5b = (N5a+N7a)/2.0 # ancestor between 12 and 34
	
	N1c = (N1b+N5b)/2.0 # final ancestor

	N1a_multilocus = [ N1a * i for i in bf ] # ancestor between 1 and 2
	N3a_multilocus = [ N3a * i for i in bf ] # ancestor between 3 and 4
	N1b_multilocus = [ N1b * i for i in bf ] # ancestor between 12 and 34
	N5a_multilocus = [ N5a * i for i in bf ] # ancestor between 5 and 6
	N7a_multilocus = [ N7a * i for i in bf ] # ancestor between 7 and 8
	N5b_multilocus = [ N5b * i for i in bf ] # ancestor between 56 and 78
	N1c_multilocus = [ N1c * i for i in bf ] # ancestor between 1234 and 5678
	
	# MIGRATION RATES
	migration_matrix = []
	for i in range(1,8,1):
		for j in range(i+1,9,1):
			migration_matrix.append(choice(list_M))
			migration_matrix.append(choice(list_M))
			#migration_matrix.append(uniform(M_prior[0], M_prior[1]))
			#migration_matrix.append(uniform(M_prior[0], M_prior[1]))
			#res+="-m {0} {1} tbs -m {1} {0} tbs ".format(i,j)

	# priorfile
	priorfile_tmp = "\t".join([ str(i) for i in [N1, N2, N3, N4, N5, N6, N7, N8, tau_1, tau_2, tau_3, tau_4, tau_5, tau_6, tau_7] ])
	priorfile_tmp += "\t"
	priorfile_tmp += "\t".join([ str(i) for i in migration_matrix ])
	priorfile += "{0}\n".format(priorfile_tmp)

	# print values for simulator	
	for locus_i in range(nLoci):
		# msms tbs nIter.nLoci -S tbs -r tbs tbs -I 8 tbs tbs tbs tbs tbs tbs tbs tbs 0
		tmp = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t".format(nSamTot_multilocus[locus_i], nSNPs[locus_i], rho[locus_i], L[locus_i], n1sam[locus_i], n2sam[locus_i], n3sam[locus_i], n4sam[locus_i], n5sam[locus_i], n6sam[locus_i], n7sam[locus_i], n8sam[locus_i])
		
		# -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 1 4 tbs -m 4 1 tbs -m 1 5 tbs -m 5 1 tbs -m 1 6 tbs -m 6 1 tbs -m 1 7 tbs -m 7 1 tbs -m 1 8 tbs -m 8 1 tbs -m 2 3 tbs -m 3 2 tbs -m 2 4 tbs -m 4 2 tbs -m 2 5 tbs -m 5 2 tbs -m 2 6 tbs -m 6 2 tbs -m 2 7 tbs -m 7 2 tbs -m 2 8 tbs -m 8 2 tbs -m 3 4 tbs -m 4 3 tbs -m 3 5 tbs -m 5 3 tbs -m 3 6 tbs -m 6 3 tbs -m 3 7 tbs -m 7 3 tbs -m 3 8 tbs -m 8 3 tbs -m 4 5 tbs -m 5 4 tbs -m 4 6 tbs -m 6 4 tbs -m 4 7 tbs -m 7 4 tbs -m 4 8 tbs -m 8 4 tbs -m 5 6 tbs -m 6 5 tbs -m 5 7 tbs -m 7 5 tbs -m 5 8 tbs -m 8 5 tbs -m 6 7 tbs -m 7 6 tbs -m 6 8 tbs -m 8 6 tbs -m 7 8 tbs -m 8 7 tbs
		tmp += "{0}\t".format("\t".join([ str(i) for i in migration_matrix ]))
		
		# -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -n 5 tbs -n 6 tbs -n 7 tbs -n 8 tbs
		tmp += "{0}\t".format("\t".join([ str(i) for i in [N1_multilocus[locus_i], N2_multilocus[locus_i], N3_multilocus[locus_i], N4_multilocus[locus_i], N5_multilocus[locus_i], N6_multilocus[locus_i], N7_multilocus[locus_i], N8_multilocus[locus_i] ] ]))
		
		# -ej tbs 2 1 -en tbs 1 tbs
		tmp += "{0}\t{0}\t{1}\t".format(tau_1, N1a_multilocus[locus_i])
		
		# -ej tbs 4 3 -en tbs 3 tbs
		tmp += "{0}\t{0}\t{1}\t".format(tau_2, N3a_multilocus[locus_i])
		
		# -ej tbs 6 5 -en tbs 5 tbs
		tmp += "{0}\t{0}\t{1}\t".format(tau_3, N5a_multilocus[locus_i])
		
		# -ej tbs 8 7 -en tbs 7 tbs
		tmp += "{0}\t{0}\t{1}\t".format(tau_4, N7a_multilocus[locus_i])
		
		# -ej tbs 3 1 -en tbs 1 tbs 
		tmp += "{0}\t{0}\t{1}\t".format(tau_5, N1b_multilocus[locus_i])
		
		# -ej tbs 7 5 -en tbs 5 tbs
		tmp += "{0}\t{0}\t{1}\t".format(tau_6, N5b_multilocus[locus_i])

		# -ej tbs 5 1 -en tbs 1 tbs
		tmp += "{0}\t{0}\t{1}".format(tau_7, N1c_multilocus[locus_i])
		
		print(tmp)

outfile = open("priorfile.txt", "w")
outfile.write(priorfile)
outfile.close()

# msms tbs 40 -s tbs -r tbs tbs -I 8 tbs tbs tbs tbs tbs tbs tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 1 4 tbs -m 4 1 tbs -m 1 5 tbs -m 5 1 tbs -m 1 6 tbs -m 6 1 tbs -m 1 7 tbs -m 7 1 tbs -m 1 8 tbs -m 8 1 tbs -m 2 3 tbs -m 3 2 tbs -m 2 4 tbs -m 4 2 tbs -m 2 5 tbs -m 5 2 tbs -m 2 6 tbs -m 6 2 tbs -m 2 7 tbs -m 7 2 tbs -m 2 8 tbs -m 8 2 tbs -m 3 4 tbs -m 4 3 tbs -m 3 5 tbs -m 5 3 tbs -m 3 6 tbs -m 6 3 tbs -m 3 7 tbs -m 7 3 tbs -m 3 8 tbs -m 8 3 tbs -m 4 5 tbs -m 5 4 tbs -m 4 6 tbs -m 6 4 tbs -m 4 7 tbs -m 7 4 tbs -m 4 8 tbs -m 8 4 tbs -m 5 6 tbs -m 6 5 tbs -m 5 7 tbs -m 7 5 tbs -m 5 8 tbs -m 8 5 tbs -m 6 7 tbs -m 7 6 tbs -m 6 8 tbs -m 8 6 tbs -m 7 8 tbs -m 8 7 tbs -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -n 5 tbs -n 6 tbs -n 7 tbs -n 8 tbs -ej tbs 2 1 -en tbs 1 tbs -ej tbs 4 3 -en tbs 3 tbs -ej tbs 6 5 -en tbs 5 tbs -ej tbs 8 7 -en tbs 7 tbs -ej tbs 3 1 -en tbs 1 tbs -ej tbs 7 5 -en tbs 5 tbs -ej tbs 5 1 -en tbs 1 tbs

# msms tbs nIter.nLoci -S tbs -r tbs tbs -I 8 tbs tbs tbs tbs tbs tbs tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 1 4 tbs -m 4 1 tbs -m 1 5 tbs -m 5 1 tbs -m 1 6 tbs -m 6 1 tbs -m 1 7 tbs -m 7 1 tbs -m 1 8 tbs -m 8 1 tbs -m 2 3 tbs -m 3 2 tbs -m 2 4 tbs -m 4 2 tbs -m 2 5 tbs -m 5 2 tbs -m 2 6 tbs -m 6 2 tbs -m 2 7 tbs -m 7 2 tbs -m 2 8 tbs -m 8 2 tbs -m 3 4 tbs -m 4 3 tbs -m 3 5 tbs -m 5 3 tbs -m 3 6 tbs -m 6 3 tbs -m 3 7 tbs -m 7 3 tbs -m 3 8 tbs -m 8 3 tbs -m 4 5 tbs -m 5 4 tbs -m 4 6 tbs -m 6 4 tbs -m 4 7 tbs -m 7 4 tbs -m 4 8 tbs -m 8 4 tbs -m 5 6 tbs -m 6 5 tbs -m 5 7 tbs -m 7 5 tbs -m 5 8 tbs -m 8 5 tbs -m 6 7 tbs -m 7 6 tbs -m 6 8 tbs -m 8 6 tbs -m 7 8 tbs -m 8 7 tbs -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -n 5 tbs -n 6 tbs -n 7 tbs -n 8 tbs -ej tbs 2 1 -en tbs 1 tbs -ej tbs 4 3 -en tbs 3 tbs -ej tbs 6 5 -en tbs 5 tbs -ej tbs 8 7 -en tbs 7 tbs -ej tbs 3 1 -en tbs 1 tbs -ej tbs 7 5 -en tbs 5 tbs -ej tbs 5 1 -en tbs 1 tbs

