# splits large VCF files in subVCF:
bsub -q dee-hugemem <<<"/scratch/ul/monthly/croux/ABCh/heliconius/scripts/split_vcf /scratch/ul/monthly/croux/ABCh/ama10.Hmel2.bwa.default.HC.vcf ama 20"
bsub -q dee-hugemem <<<"/scratch/ul/monthly/croux/ABCh/heliconius/scripts/split_vcf /scratch/ul/monthly/croux/ABCh/chi10.Hmel2.bwa.default.HC.vcf chi 20"
bsub -q dee-hugemem <<<"/scratch/ul/monthly/croux/ABCh/heliconius/scripts/split_vcf /scratch/ul/monthly/croux/ABCh/flo10.Hmel2.bwa.default.HC.vcf flo 20"
bsub -q dee-hugemem <<<"/scratch/ul/monthly/croux/ABCh/heliconius/scripts/split_vcf /scratch/ul/monthly/croux/ABCh/mal10.Hmel2.bwa.default.HC.vcf mal 20"
bsub -q dee-hugemem <<<"/scratch/ul/monthly/croux/ABCh/heliconius/scripts/split_vcf /scratch/ul/monthly/croux/ABCh/ros10.Hmel2.bwa.default.HC.vcf ros 20"
bsub -q dee-hugemem <<<"/scratch/ul/monthly/croux/ABCh/heliconius/scripts/split_vcf /scratch/ul/monthly/croux/ABCh/txn10.Hmel2.bwa.default.HC.vcf txn 20"
bsub -q dee-hugemem <<<"/scratch/ul/monthly/croux/ABCh/heliconius/scripts/split_vcf /scratch/ul/monthly/croux/ABCh/vul10.Hmel2.bwa.default.HC.vcf vul 20"
bsub -q dee-hugemem <<<"/scratch/ul/monthly/croux/ABCh/heliconius/scripts/split_vcf /scratch/ul/monthly/croux/ABCh/zel10.Hmel2.bwa.default.HC.vcf zel 20"


# produces alignments of CDS in fasta:
cd /scratch/ul/monthly/croux/ABCh/ama/fasta
for species in ama chi flo mal ros txn vul zel; do
        cd /scratch/ul/monthly/croux/ABCh/${species}/fasta
        for i in $(ls ../subVCF_${species}_*); do bsub -q dee-hugemem -R "rusage[mem=8000]" -M 8000000 <<<"/home/croux/work/heliconius/scripts/vcf_filter.py ../../Hmel2.fa ../../Hmel2.gff ${i} 3 ../corr_table_${species}.txt ${species}"; done
done


# converts fasta in inputs for ABC
# phasing
for species in ama chi flo mal ros txn vul zel; do
	mkdir /scratch/ul/monthly/croux/ABCh/${species}/vcf
	mv /scratch/ul/monthly/croux/ABCh/${species}/*.vcf /scratch/ul/monthly/croux/ABCh/${species}/vcf/
done

for species in ama chi flo mal ros txn vul zel; do
	mkdir /scratch/ul/monthly/croux/ABCh/${species}/fasta/unphased
	mkdir /scratch/ul/monthly/croux/ABCh/${species}/fasta/phased
	mv /scratch/ul/monthly/croux/ABCh/${species}/fasta/*.fas /scratch/ul/monthly/croux/ABCh/${species}/fasta/unphased/
done

cnt=0
for species in ama chi flo mal ros txn vul zel; do
	echo $species
	cd /scratch/ul/monthly/croux/ABCh/${species}/fasta/phased
	for locus in $(ls /scratch/ul/monthly/croux/ABCh/${species}/fasta/unphased/*.fas); do
		cnt=$((cnt+1))
		bsub -q dee-hugemem <<< "fasta2phase.py ${locus}";
		if (( cnt % 100 == 0 )); then sleep 1; fi
	done
	sleep 5
done
 

