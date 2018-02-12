#!/usr/bin/Rscript
# script that takes the Hmel*_Fst.txt files as inputs and produces 2 columns for FOR and REV strands with positions:
# intron / ig / pseudo / pos1 / pos2 / pos3


input = commandArgs()[6]
x = read.table(input, h=T)

pos_FOR = as.character(x$feature_strandFOR)
pos_REV = as.character(x$feature_strandREV)

genes_FOR = table(x$geneName_strandFOR)
genes_FOR = names(genes_FOR[grep("gene",names(genes_FOR))]) 

genes_REV = table(x$geneName_strandREV)
genes_REV = names(genes_REV[grep("gene",names(genes_REV))]) 


for(i in genes_FOR){
	y = grep(i, x$geneName_strandFOR)
	cnt = -1
	for(j in 1:length(y)){
		pos = y[j]
		if( is.na(pos_FOR[pos]) == FALSE ){
			if( pos_FOR[pos] == "CDS" ){
				cnt = cnt + 1
				pos_FOR[pos] = paste("pos_", cnt%%3+1, sep="")
			}
			if( pos_FOR[pos] == "mRNA" ){
				pos_FOR[pos] = "intron"
			}
		}
	}
}


for(i in genes_REV){
	y = grep(i, x$geneName_strandREV)
	cnt = -1
	for(j in length(y):1){
		pos = y[j]
		if( is.na(pos_REV[pos]) == FALSE ){
			if( pos_REV[pos] == "CDS" ){
				cnt = cnt + 1
				pos_REV[pos] = paste("pos_", cnt%%3+1, sep="")
			}
			if( pos_REV[pos] == "mRNA" ){
				pos_REV[pos] = "intron"
			}
		}
	}
}

x = cbind(x[,1:4], pos_FOR, x[,5:6], pos_REV, x[,7:96])

x = na.omit(x)

output = paste("annot", input, sep="_")
write.table(x, output, col.names=T, row.names=F, sep="\t", quote=F)

