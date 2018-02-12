# reads the gff file and gets the annotation
for i in $(cat Hmel2.gff | grep "##" | grep Hmel | awk '{print $2}'); do
	sbatch --nodes=1 --ntasks-per-node=1 --time=2:00:00 -J ${i} --wrap="./extract.py  Hmel2.gff ${i}"
	done

# counts the number of alleles in different populations
for i in $(cat Hmel2.gff | grep "##" | grep Hmel | awk '{print $2}'); do
        sbatch --nodes=1 --ntasks-per-node=1 --time=2:00:00 -J ${i} --wrap="./count_alleles.py ${i}"
done

# compute some statistics at each position: Hs, Fst, {SxA, SxB, Ss, Sf}
for i in $(cat Hmel2.gff | grep "##" | grep Hmel | awk '{print $2}'); do
        sbatch --nodes=1 --ntasks-per-node=1 --time=2:00:00 -J ${i} --wrap="./Fst.py ${i}"
done

#3 compute Fst and annotate all individual positions:
for i in $(cat Hmel2.gff | grep "##" | grep Hmel | awk '{print $2}'); do
        echo "./Fst.py ${i}; ./positions.R ${i}_Fst.txt" >>bsub.txt; done
done

#4 summarize diversity at pos1, pos2, pos3, Ig, pseudo, tRNA, miRNA, intron for all of the 8 populations and auto versus sexChro
./summarize_Div.py # global divergence
./summarize_Div_tmp.py # divergence per contig
./summarize_Hs.py # global diversity



