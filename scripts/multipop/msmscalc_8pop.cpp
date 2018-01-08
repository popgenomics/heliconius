#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>

void get_sample_sizes(std::vector <unsigned int> & n1, std::vector <unsigned int> & n2, std::vector <unsigned int> & n3, std::vector <unsigned int> & n4, std::vector <unsigned int> & n5, std::vector <unsigned int> & n6, std::vector <unsigned int> & n7, std::vector <unsigned int> & n8, std::vector <unsigned int> & nTot, std::vector <size_t> & nSNPs, unsigned int & nLoci);

void compute_stats(const size_t nSNP, const unsigned int nsam1, const unsigned int nsam2, const unsigned int nsam3, const unsigned int nsam4, const unsigned int nsam5, const unsigned int nsam6, const unsigned int nsam7, const unsigned int nsam8, const std::vector <std::string> & block_haplotypes, std::vector <unsigned int> & sx_1, std::vector <unsigned int> & sx_2, std::vector <unsigned int> & sx_3, std::vector <unsigned int> & sx_4, std::vector <unsigned int> & sx_5, std::vector <unsigned int> & sx_6, std::vector <unsigned int> & sx_7, std::vector <unsigned int> & sx_8, std::vector <unsigned int> & sx_12, std::vector <unsigned int> & sx_34, std::vector <unsigned int> & sx_56, std::vector <unsigned int> & sx_78, std::vector <unsigned int> & sx_1234, std::vector <unsigned int> & sx_5678, std::vector <unsigned int> & sf_1_2, std::vector <unsigned int> & sf_3_4, std::vector <unsigned int> & sf_5_6, std::vector <unsigned int> & sf_7_8, std::vector <unsigned int> & sf_12_34, std::vector <unsigned int> & sf_56_78, std::vector <unsigned int> & sf_1234_5678, std::vector <unsigned int> & ss_1_2, std::vector <unsigned int> & ss_3_4, std::vector <unsigned int> & ss_5_6, std::vector <unsigned int> & ss_7_8, std::vector <unsigned int> & ss_12_34, std::vector <unsigned int> & ss_56_78, std::vector <unsigned int> & ss_1234_5678);

void get_catFreq(const float & freq, const size_t & nsam, float & pop, unsigned int & cat);

void print_stats(const unsigned int & nLoci, const size_t & nDataset, std::vector <unsigned int> & sx_1, std::vector <unsigned int> & sx_2, std::vector <unsigned int> & sx_3, std::vector <unsigned int> & sx_4, std::vector <unsigned int> & sx_5, std::vector <unsigned int> & sx_6, std::vector <unsigned int> & sx_7, std::vector <unsigned int> & sx_8, std::vector <unsigned int> & sx_12, std::vector <unsigned int> & sx_34, std::vector <unsigned int> & sx_56, std::vector <unsigned int> & sx_78, std::vector <unsigned int> & sx_1234, std::vector <unsigned int> & sx_5678, std::vector <unsigned int> & sf_1_2, std::vector <unsigned int> & sf_3_4, std::vector <unsigned int> & sf_5_6, std::vector <unsigned int> & sf_7_8, std::vector <unsigned int> & sf_12_34, std::vector <unsigned int> & sf_56_78, std::vector <unsigned int> & sf_1234_5678, std::vector <unsigned int> & ss_1_2, std::vector <unsigned int> & ss_3_4, std::vector <unsigned int> & ss_5_6, std::vector <unsigned int> & ss_7_8, std::vector <unsigned int> & ss_12_34, std::vector <unsigned int> & ss_56_78, std::vector <unsigned int> & ss_1234_5678);

void mean_deviation(const std::vector <unsigned int> & tableau, const unsigned int & nLoci, float & mean, float & sd);

// g++ msmscalc_8pop.cpp -std=c++17 -O3 -o msmscalc_multipop
int main(int argc, char* argv[]){
	size_t i(0);
	unsigned int locus = 0;
	unsigned int nLoci = 0;
	
	std::vector <size_t> nSNPs;
	std::vector <unsigned int> n1;
	std::vector <unsigned int> n2;
	std::vector <unsigned int> n3;
	std::vector <unsigned int> n4;
	std::vector <unsigned int> n5;
	std::vector <unsigned int> n6;
	std::vector <unsigned int> n7;
	std::vector <unsigned int> n8;
	std::vector <unsigned int> nTot;

	get_sample_sizes(n1, n2, n3, n4, n5, n6, n7, n8, nTot, nSNPs, nLoci);

	for( i=0; i<nLoci; ++i ){
		nTot.push_back(n1[i]+n2[i]+n3[i]+n4[i]+n5[i]+n6[i]+n7[i]+n8[i]);
	}
	
	// summary statistics
	std::vector <unsigned int> sx_1;
	std::vector <unsigned int> sx_2;
	std::vector <unsigned int> sx_3;
	std::vector <unsigned int> sx_4;
	std::vector <unsigned int> sx_5;
	std::vector <unsigned int> sx_6;
	std::vector <unsigned int> sx_7;
	std::vector <unsigned int> sx_8;
	std::vector <unsigned int> sx_12;
	std::vector <unsigned int> sx_34;
	std::vector <unsigned int> sx_56;
	std::vector <unsigned int> sx_78;
	std::vector <unsigned int> sx_1234;
	std::vector <unsigned int> sx_5678;
	
	std::vector <unsigned int> sf_1_2;
	std::vector <unsigned int> sf_3_4;
	std::vector <unsigned int> sf_5_6;
	std::vector <unsigned int> sf_7_8;
	std::vector <unsigned int> sf_12_34;
	std::vector <unsigned int> sf_56_78;
	std::vector <unsigned int> sf_1234_5678;

	std::vector <unsigned int> ss_1_2;
	std::vector <unsigned int> ss_3_4;
	std::vector <unsigned int> ss_5_6;
	std::vector <unsigned int> ss_7_8;
	std::vector <unsigned int> ss_12_34;
	std::vector <unsigned int> ss_56_78;
	std::vector <unsigned int> ss_1234_5678;
	
	
	size_t nDataset(0); // count the number of simulated datasets over the whole msmsFile
	unsigned int test(0);
	unsigned int nHaplotype(0);

	std::vector <std::string> block_haplotypes;
	// read the msms output from the stdin
	std::string line; // contains the msms line
	for(line; std::getline(std::cin, line);){
		
		// new loci
		// 'segsites' line
		if(line[0]=='s'){ // if the line starts by a 's', then we expect: 'segsites: int'
			test=1;
			continue; // stop reading the 'segites' line
		} // end of: "if the line contains the string: 'segsites'


		/* If the line is not empty, and doesn't start by '/', nor 'segsites', nor 'positions'*/
		if(line!="" && line[0]!='/' && line[0]!='s' && line[0]!='p' && test==1){
			++nHaplotype;
			block_haplotypes.push_back( line );
//			std::cout << "dataset : " << nDataset << "\tlocus : " << locus << "\tline : " << nHaplotype << ", " << line << std::endl;

			if( nHaplotype == nTot[locus] ){ // if the haplotype is the last of a monolocus block
				nHaplotype = 0;
				
				// treat the alignement of the locus i
				//compute_stats();
				compute_stats(nSNPs[locus], n1[locus], n2[locus], n3[locus], n4[locus], n5[locus], n6[locus], n7[locus], n8[locus], block_haplotypes, sx_1, sx_2, sx_3, sx_4, sx_5, sx_6, sx_7, sx_8, sx_12, sx_34, sx_56, sx_78, sx_1234, sx_5678, sf_1_2, sf_3_4, sf_5_6, sf_7_8, sf_12_34, sf_56_78, sf_1234_5678, ss_1_2, ss_3_4, ss_5_6, ss_7_8, ss_12_34, ss_56_78, ss_1234_5678 );
//				std::cout << "POUET: " << block_haplotypes.size() << std::endl;
				
				if( locus == nLoci-1 ){
					++nDataset;
					print_stats(nLoci, nDataset, sx_1, sx_2, sx_3, sx_4, sx_5, sx_6, sx_7, sx_8, sx_12, sx_34, sx_56, sx_78, sx_1234, sx_5678, sf_1_2, sf_3_4, sf_5_6, sf_7_8, sf_12_34, sf_56_78, sf_1234_5678, ss_1_2, ss_3_4, ss_5_6, ss_7_8, ss_12_34, ss_56_78, ss_1234_5678);
					locus = 0;
				}else{
					++locus;
				}
				
				block_haplotypes.clear();
			}

		}
		
	}
	
/*	for(i=0; i<sx_1.size(); ++i){
		std::cout << sx_1[i] << " ";
	}
	std::cout << std::endl;*/
	return 0;
}


void get_sample_sizes(std::vector <unsigned int> & n1, std::vector <unsigned int> & n2, std::vector <unsigned int> & n3, std::vector <unsigned int> & n4, std::vector <unsigned int> & n5, std::vector <unsigned int> & n6, std::vector <unsigned int> & n7, std::vector <unsigned int> & n8, std::vector <unsigned int> & nTot, std::vector <size_t> & nSNPs, unsigned int & nLoci){
	// read the bpfile to get some informations about sample sizes
	//ex :
	/*
	# nSNPs L n1 n2 n3 n4 n5 n6 n7 n8 rho
	10      20      30      40
	100     200     300     400
	4       4       4       4
	4       4       4       4
	4       4       4       4
	4       4       4       4
	4       4       4       4
	4       4       4       4
	4       4       4       4
	4       4       4       4
	5       5       5       5
	*/

	size_t i(0);	
	size_t count_line(0);
	std::string line;
	std::string word;
	std::ifstream infile("bpfile");
	if( infile ){
		while (std::getline(infile, line)){
			++count_line;
			std::istringstream iss(line);
			
			while( std::getline( iss, word, '\t' ) ){
				if( count_line == 2 ){ nSNPs.push_back( std::stoi(word) ); }
				if( count_line == 4 ){ n1.push_back( std::stoi(word) ); }
				if( count_line == 5 ){ n2.push_back( std::stoi(word) ); }
				if( count_line == 6 ){ n3.push_back( std::stoi(word) ); }
				if( count_line == 7 ){ n4.push_back( std::stoi(word) ); }
				if( count_line == 8 ){ n5.push_back( std::stoi(word) ); }
				if( count_line == 9 ){ n6.push_back( std::stoi(word) ); }
				if( count_line == 10 ){ n7.push_back( std::stoi(word) ); }
				if( count_line == 11 ){ n8.push_back( std::stoi(word) ); }
			}
		}
	
		for( i=0; i<n1.size(); ++i ){
			nTot.push_back(n1[i]+n2[i]+n3[i]+n4[i]+n5[i]+n6[i]+n7[i]+n8[i]);
//			std::cout << "total number individuals for loci #" << i << " : " << nTot[i] << " (" << nSNPs[i] << " SNPs)" << std::endl;
		}
	
	nLoci = nSNPs.size();

	}else{
		std::cerr << "Couldn't open bpfile for reading\n";
	}
}

void compute_stats(const size_t nSNP, const unsigned int nsam1, const unsigned int nsam2, const unsigned int nsam3, const unsigned int nsam4, const unsigned int nsam5, const unsigned int nsam6, const unsigned int nsam7, const unsigned int nsam8, const std::vector <std::string> & block_haplotypes, std::vector <unsigned int> & sx_1, std::vector <unsigned int> & sx_2, std::vector <unsigned int> & sx_3, std::vector <unsigned int> & sx_4, std::vector <unsigned int> & sx_5, std::vector <unsigned int> & sx_6, std::vector <unsigned int> & sx_7, std::vector <unsigned int> & sx_8, std::vector <unsigned int> & sx_12, std::vector <unsigned int> & sx_34, std::vector <unsigned int> & sx_56, std::vector <unsigned int> & sx_78, std::vector <unsigned int> & sx_1234, std::vector <unsigned int> & sx_5678, std::vector <unsigned int> & sf_1_2, std::vector <unsigned int> & sf_3_4, std::vector <unsigned int> & sf_5_6, std::vector <unsigned int> & sf_7_8, std::vector <unsigned int> & sf_12_34, std::vector <unsigned int> & sf_56_78, std::vector <unsigned int> & sf_1234_5678, std::vector <unsigned int> & ss_1_2, std::vector <unsigned int> & ss_3_4, std::vector <unsigned int> & ss_5_6, std::vector <unsigned int> & ss_7_8, std::vector <unsigned int> & ss_12_34, std::vector <unsigned int> & ss_56_78, std::vector <unsigned int> & ss_1234_5678){
	size_t pos(0);
	size_t ind(0);

	// 0, 0.5, 1 for monomorphic 0, polymorphism, monomorphic 1	
	float pop_1(0.0);
	float pop_2(0.0);
	float pop_3(0.0);
	float pop_4(0.0);
	float pop_5(0.0);
	float pop_6(0.0);
	float pop_7(0.0);
	float pop_8(0.0);
	
	// 0 for monomorphic, 1 for polymorphic
	unsigned int cat_1(0.0);
	unsigned int cat_2(0.0);
	unsigned int cat_3(0.0);
	unsigned int cat_4(0.0);
	unsigned int cat_5(0.0);
	unsigned int cat_6(0.0);
	unsigned int cat_7(0.0);
	unsigned int cat_8(0.0);

	std::vector <float> freq_1;
	std::vector <float> freq_2;
	std::vector <float> freq_3;
	std::vector <float> freq_4;
	std::vector <float> freq_5;
	std::vector <float> freq_6;
	std::vector <float> freq_7;
	std::vector <float> freq_8;

	// results
	unsigned int sx_1_tmp(0);
	unsigned int sx_2_tmp(0);
	unsigned int sx_3_tmp(0);
	unsigned int sx_4_tmp(0);
	unsigned int sx_5_tmp(0);
	unsigned int sx_6_tmp(0);
	unsigned int sx_7_tmp(0);
	unsigned int sx_8_tmp(0);
	unsigned int sx_12_tmp(0);
	unsigned int sx_34_tmp(0);
	unsigned int sx_56_tmp(0);
	unsigned int sx_78_tmp(0);
	unsigned int sx_1234_tmp(0);
	unsigned int sx_5678_tmp(0);

	unsigned int sf_1_2_tmp(0);
	unsigned int sf_3_4_tmp(0);
	unsigned int sf_5_6_tmp(0);
	unsigned int sf_7_8_tmp(0);
	unsigned int sf_12_34_tmp(0);
	unsigned int sf_56_78_tmp(0);
	unsigned int sf_1234_5678_tmp(0);

	unsigned int ss_1_2_tmp(0);
	unsigned int ss_3_4_tmp(0);
	unsigned int ss_5_6_tmp(0);
	unsigned int ss_7_8_tmp(0);
	unsigned int ss_12_34_tmp(0);
	unsigned int ss_56_78_tmp(0);
	unsigned int ss_1234_5678_tmp(0);

	for( pos=0; pos<nSNP; ++pos){ // loop over SNPs
		freq_1.push_back(0.0);
		freq_2.push_back(0.0);
		freq_3.push_back(0.0);
		freq_4.push_back(0.0);
		freq_5.push_back(0.0);
		freq_6.push_back(0.0);
		freq_7.push_back(0.0);
		freq_8.push_back(0.0);
		
		for( ind=0; ind<nsam1; ++ind ){ // loop over individuals
			if( block_haplotypes[ind][pos]=='1'){ ++freq_1[pos];}
		}
		
		for( ind=nsam1; ind<(nsam1+nsam2); ++ind ){
			if( block_haplotypes[ind][pos]=='1'){ ++freq_2[pos];}
		}
		
		for( ind=nsam1+nsam2; ind<(nsam1+nsam2+nsam3); ++ind ){
			if( block_haplotypes[ind][pos]=='1'){ ++freq_3[pos];}
		}
		
		for( ind=nsam1+nsam2+nsam3; ind<(nsam1+nsam2+nsam3+nsam4); ++ind ){
			if( block_haplotypes[ind][pos]=='1'){ ++freq_4[pos];}
		}
		
		for( ind=nsam1+nsam2+nsam3+nsam4; ind<(nsam1+nsam2+nsam3+nsam4+nsam5); ++ind ){
			if( block_haplotypes[ind][pos]=='1'){ ++freq_5[pos];}
		}
		
		for( ind=nsam1+nsam2+nsam3+nsam4+nsam5; ind<(nsam1+nsam2+nsam3+nsam4+nsam5+nsam6); ++ind ){
			if( block_haplotypes[ind][pos]=='1'){ ++freq_6[pos];}
		}
		
		for( ind=nsam1+nsam2+nsam3+nsam4+nsam5+nsam6; ind<(nsam1+nsam2+nsam3+nsam4+nsam5+nsam6+nsam7); ++ind ){
			if( block_haplotypes[ind][pos]=='1'){ ++freq_7[pos];}
		}
		
		for( ind=nsam1+nsam2+nsam3+nsam4+nsam5+nsam6+nsam7; ind<(nsam1+nsam2+nsam3+nsam4+nsam5+nsam6+nsam7+nsam8); ++ind ){
			if( block_haplotypes[ind][pos]=='1'){ ++freq_8[pos];}
		} // end of loop over individuals
	} // end of loop over SNPs

	for( pos=0; pos<nSNP; ++pos ){
		get_catFreq(freq_1[pos], nsam1, pop_1, cat_1);
		get_catFreq(freq_2[pos], nsam2, pop_2, cat_2);
		get_catFreq(freq_3[pos], nsam3, pop_3, cat_3);
		get_catFreq(freq_4[pos], nsam4, pop_4, cat_4);
		get_catFreq(freq_5[pos], nsam5, pop_5, cat_5);
		get_catFreq(freq_6[pos], nsam6, pop_6, cat_6);
		get_catFreq(freq_7[pos], nsam7, pop_7, cat_7);
		get_catFreq(freq_8[pos], nsam8, pop_8, cat_8);
		
		// sx
		// sx_1, sx_2, sx_3, sx_4, sx_5, sx_6, sx_7, sx_8, sx_12, sx_34, sx_56, sx_78, sx_1234, sx_5678
		if( cat_1 == 1 && cat_2+cat_3+cat_4+cat_5+cat_6+cat_7+cat_8 == 0){
			++sx_1_tmp; 
		}
		if( cat_2 == 1 && cat_1+cat_3+cat_4+cat_5+cat_6+cat_7+cat_8 == 0){
			++sx_2_tmp; 
		}
		if( cat_3 == 1 && cat_1+cat_2+cat_4+cat_5+cat_6+cat_7+cat_8 == 0){
			++sx_3_tmp; 
		}
		if( cat_4 == 1 && cat_1+cat_2+cat_3+cat_5+cat_6+cat_7+cat_8 == 0){
			++sx_4_tmp; 
		}
		if( cat_5 == 1 && cat_1+cat_2+cat_3+cat_4+cat_6+cat_7+cat_8 == 0){
			++sx_5_tmp; 
		}
		if( cat_6 == 1 && cat_1+cat_2+cat_3+cat_4+cat_5+cat_7+cat_8 == 0){
			++sx_6_tmp; 
		}
		if( cat_7 == 1 && cat_1+cat_2+cat_3+cat_4+cat_5+cat_6+cat_8 == 0){
			++sx_7_tmp; 
		}
		if( cat_8 == 1 && cat_1+cat_2+cat_3+cat_4+cat_5+cat_6+cat_7 == 0){
			++sx_8_tmp; 
		}

		if( cat_1 == 1 && cat_2 == 1 && cat_3+cat_4+cat_5+cat_6+cat_7+cat_8 == 0){
			++sx_12_tmp; 
		}

		if( cat_3 == 1 && cat_4 == 1 && cat_1+cat_2+cat_5+cat_6+cat_7+cat_8 == 0){
			++sx_34_tmp; 
		}

		if( cat_5 == 1 && cat_6 == 1 && cat_1+cat_2+cat_3+cat_4+cat_7+cat_8 == 0){
			++sx_56_tmp; 
		}

		if( cat_7 == 1 && cat_8 == 1 && cat_1+cat_2+cat_3+cat_4+cat_5+cat_6 == 0){
			++sx_78_tmp; 
		}

		if( cat_1 == 1 && cat_2 == 1 && cat_3 == 1 && cat_4 == 1 && cat_5+cat_6+cat_7+cat_8 == 0){
			++sx_1234_tmp; 
		}

		if( cat_5 == 1 && cat_6 == 1 && cat_7 == 1 && cat_8 == 1 && cat_1+cat_2+cat_3+cat_4 == 0){
			++sx_5678_tmp; 
		}

		// sf
		// sf_1_2, sf_3_4, sf_5_6, sf_7_8, sf_12_34, sf_56_78, sf_1234_5678
		if( cat_1 == 0 && cat_2 == 0 ){
			if( pop_1 != pop_2 ){
				++sf_1_2_tmp;
			}
		}

		if( cat_3 == 0 && cat_4 == 0 ){
			if( pop_3 != pop_4 ){
				++sf_3_4_tmp;
			}
		}
		
		if( cat_5 == 0 && cat_6 == 0 ){
			if( pop_5 != pop_6 ){
				++sf_5_6_tmp;
			}
		}

		if( cat_7 == 0 && cat_8 == 0 ){
			if( pop_7 != pop_8 ){
				++sf_7_8_tmp;
			}
		}

		if( cat_1 == 0 && cat_2 == 0 && cat_3 == 0 && cat_4 == 0 ){
			if( pop_1 == pop_2 ){
				if( pop_3 == pop_4 ){
					if( pop_1 != pop_3){
						++sf_12_34_tmp;
					}
				}
			}
		}

		if( cat_5 == 0 && cat_6 == 0 && cat_7 == 0 && cat_8 == 0 ){
			if( pop_5 == pop_6 ){
				if( pop_7 == pop_8 ){
					if( pop_5 != pop_7){
						++sf_56_78_tmp;
					}
				}
			}
		}

		if( cat_1 + cat_2 + cat_3 + cat_4 + cat_5 + cat_6 + cat_7 + cat_8 == 0){
			if( pop_1 + pop_2 + pop_3 + pop_4 == 0 ){ // if populations 1:4 fixed the same allele '0'
				if( pop_5 + pop_6 + pop_7 + pop_8 == 4){ // if populations 5:8 fixed the same allele '1'
					++sf_1234_5678_tmp;
				}
			}else if( pop_1 + pop_2 + pop_3 + pop_4 == 4){
				if( pop_5 + pop_6 + pop_7 + pop_8 == 0){ // if populations 5:8 fixed the same allele '9'
					++sf_1234_5678_tmp;
				}
			}
		}
		
		// ++ss
		// ++ss_1_2, ++ss_3_4, ++ss_5_6, ++ss_7_8, ++ss_12_34, ++ss_56_78, ++ss_1234_5678
		if( cat_1 + cat_2 == 2){
			++ss_1_2_tmp;
		}
		
		if( cat_3 + cat_4 == 2){
			++ss_3_4_tmp;
		}
		
		if( cat_5 + cat_6 == 2){
			++ss_5_6_tmp;
		}
	
		if( cat_7 + cat_8 == 2){	
			++ss_7_8_tmp;
		}

		if( cat_1 + cat_2 + cat_3 + cat_4 == 4){			
			++ss_12_34_tmp;
		}
		
		if( cat_5 + cat_6 + cat_7 + cat_8 == 4){
			++ss_56_78_tmp;
		}

		if( cat_1 + cat_2 + cat_3 + cat_4 + cat_5 + cat_6 + cat_7 + cat_8 == 8){
			++ss_1234_5678_tmp;
		}
	}
	sx_1.push_back(sx_1_tmp);
	sx_2.push_back(sx_2_tmp);
	sx_3.push_back(sx_3_tmp);
	sx_4.push_back(sx_4_tmp);
	sx_5.push_back(sx_5_tmp);
	sx_6.push_back(sx_6_tmp);
	sx_7.push_back(sx_7_tmp);
	sx_8.push_back(sx_8_tmp);
	sx_12.push_back(sx_12_tmp);
	sx_34.push_back(sx_34_tmp);
	sx_56.push_back(sx_56_tmp);
	sx_78.push_back(sx_78_tmp);
	sx_1234.push_back(sx_1234_tmp);
	sx_5678.push_back(sx_5678_tmp);
	
	sf_1_2.push_back(sf_1_2_tmp);
	sf_3_4.push_back(sf_3_4_tmp);
	sf_5_6.push_back(sf_5_6_tmp);
	sf_7_8.push_back(sf_7_8_tmp);
	sf_12_34.push_back(sf_12_34_tmp);
	sf_56_78.push_back(sf_56_78_tmp);
	sf_1234_5678.push_back(sf_1234_5678_tmp);
	
	ss_1_2.push_back(ss_1_2_tmp);
	ss_3_4.push_back(ss_3_4_tmp);
	ss_5_6.push_back(ss_5_6_tmp);
	ss_7_8.push_back(ss_7_8_tmp);
	ss_12_34.push_back(ss_12_34_tmp);
	ss_56_78.push_back(ss_56_78_tmp);
	ss_1234_5678.push_back(ss_1234_5678_tmp);
}


void get_catFreq(const float & freq, const size_t & nsam, float & pop, unsigned int & cat){
	if( freq == 0.0 ){
		cat = 0;
		pop = 0.0;
	}else if( freq == nsam ){
		cat = 0;
		pop = 1.0;
	}else{
		cat = 1;
		pop = 0.5;
	}
}

void print_stats(const unsigned int & nLoci, const size_t & nDataset, std::vector <unsigned int> & sx_1, std::vector <unsigned int> & sx_2, std::vector <unsigned int> & sx_3, std::vector <unsigned int> & sx_4, std::vector <unsigned int> & sx_5, std::vector <unsigned int> & sx_6, std::vector <unsigned int> & sx_7, std::vector <unsigned int> & sx_8, std::vector <unsigned int> & sx_12, std::vector <unsigned int> & sx_34, std::vector <unsigned int> & sx_56, std::vector <unsigned int> & sx_78, std::vector <unsigned int> & sx_1234, std::vector <unsigned int> & sx_5678, std::vector <unsigned int> & sf_1_2, std::vector <unsigned int> & sf_3_4, std::vector <unsigned int> & sf_5_6, std::vector <unsigned int> & sf_7_8, std::vector <unsigned int> & sf_12_34, std::vector <unsigned int> & sf_56_78, std::vector <unsigned int> & sf_1234_5678, std::vector <unsigned int> & ss_1_2, std::vector <unsigned int> & ss_3_4, std::vector <unsigned int> & ss_5_6, std::vector <unsigned int> & ss_7_8, std::vector <unsigned int> & ss_12_34, std::vector <unsigned int> & ss_56_78, std::vector <unsigned int> & ss_1234_5678){
	float mean_sx_1(0.0); float sd_sx_1(0.0);
	float mean_sx_2(0.0); float sd_sx_2(0.0);
	float mean_sx_3(0.0); float sd_sx_3(0.0);
	float mean_sx_4(0.0); float sd_sx_4(0.0);
	float mean_sx_5(0.0); float sd_sx_5(0.0);
	float mean_sx_6(0.0); float sd_sx_6(0.0);
	float mean_sx_7(0.0); float sd_sx_7(0.0);
	float mean_sx_8(0.0); float sd_sx_8(0.0);
	float mean_sx_12(0.0); float sd_sx_12(0.0);
	float mean_sx_34(0.0); float sd_sx_34(0.0);
	float mean_sx_56(0.0); float sd_sx_56(0.0);
	float mean_sx_78(0.0); float sd_sx_78(0.0);
	float mean_sx_1234(0.0); float sd_sx_1234(0.0);
	float mean_sx_5678(0.0); float sd_sx_5678(0.0);

	float mean_sf_1_2(0.0); float sd_sf_1_2(0.0);
	float mean_sf_3_4(0.0); float sd_sf_3_4(0.0);
	float mean_sf_5_6(0.0); float sd_sf_5_6(0.0);
	float mean_sf_7_8(0.0); float sd_sf_7_8(0.0);
	float mean_sf_12_34(0.0); float sd_sf_12_34(0.0);
	float mean_sf_56_78(0.0); float sd_sf_56_78(0.0);
	float mean_sf_1234_5678(0.0); float sd_sf_1234_5678(0.0);

	float mean_ss_1_2(0.0); float sd_ss_1_2(0.0);
	float mean_ss_3_4(0.0); float sd_ss_3_4(0.0);
	float mean_ss_5_6(0.0); float sd_ss_5_6(0.0);
	float mean_ss_7_8(0.0); float sd_ss_7_8(0.0);
	float mean_ss_12_34(0.0); float sd_ss_12_34(0.0);
	float mean_ss_56_78(0.0); float sd_ss_56_78(0.0);
	float mean_ss_1234_5678(0.0); float sd_ss_1234_5678(0.0);


	mean_deviation(sx_1, nLoci, mean_sx_1, sd_sx_1);
	mean_deviation(sx_2, nLoci, mean_sx_2, sd_sx_2);
	mean_deviation(sx_3, nLoci, mean_sx_3, sd_sx_3);
	mean_deviation(sx_4, nLoci, mean_sx_4, sd_sx_4);
	mean_deviation(sx_5, nLoci, mean_sx_5, sd_sx_5);
	mean_deviation(sx_6, nLoci, mean_sx_6, sd_sx_6);
	mean_deviation(sx_7, nLoci, mean_sx_7, sd_sx_7);
	mean_deviation(sx_8, nLoci, mean_sx_8, sd_sx_8);
	mean_deviation(sx_12, nLoci, mean_sx_12, sd_sx_12);
	mean_deviation(sx_34, nLoci, mean_sx_34, sd_sx_34);
	mean_deviation(sx_56, nLoci, mean_sx_56, sd_sx_56);
	mean_deviation(sx_78, nLoci, mean_sx_78, sd_sx_78);
	mean_deviation(sx_1234, nLoci, mean_sx_1234, sd_sx_1234);
	mean_deviation(sx_5678, nLoci, mean_sx_5678, sd_sx_5678);
	
	mean_deviation(sf_1_2, nLoci, mean_sf_1_2, sd_sf_1_2);
	mean_deviation(sf_3_4, nLoci, mean_sf_3_4, sd_sf_3_4);
	mean_deviation(sf_5_6, nLoci, mean_sf_5_6, sd_sf_5_6);
	mean_deviation(sf_7_8, nLoci, mean_sf_7_8, sd_sf_7_8);
	mean_deviation(sf_12_34, nLoci, mean_sf_12_34, sd_sf_12_34);
	mean_deviation(sf_56_78, nLoci, mean_sf_56_78, sd_sf_56_78);
	mean_deviation(sf_1234_5678, nLoci, mean_sf_1234_5678, sd_sf_1234_5678);
	
	mean_deviation(ss_1_2, nLoci, mean_ss_1_2, sd_ss_1_2);
	mean_deviation(ss_3_4, nLoci, mean_ss_3_4, sd_ss_3_4);
	mean_deviation(ss_5_6, nLoci, mean_ss_5_6, sd_ss_5_6);
	mean_deviation(ss_7_8, nLoci, mean_ss_7_8, sd_ss_7_8);
	mean_deviation(ss_12_34, nLoci, mean_ss_12_34, sd_ss_12_34);
	mean_deviation(ss_56_78, nLoci, mean_ss_56_78, sd_ss_56_78);
	mean_deviation(ss_1234_5678, nLoci, mean_ss_1234_5678, sd_ss_1234_5678);
	
	sx_1.clear();
	sx_2.clear();
	sx_3.clear();
	sx_4.clear();
	sx_5.clear();
	sx_6.clear();
	sx_7.clear();
	sx_8.clear();
	sx_12.clear();
	sx_34.clear();
	sx_56.clear();
	sx_78.clear();
	sx_1234.clear();
	sx_5678.clear();

	sf_1_2.clear();
	sf_3_4.clear();
	sf_5_6.clear();
	sf_7_8.clear();
	sf_12_34.clear();
	sf_56_78.clear();
	sf_1234_5678.clear();

	ss_1_2.clear();
	ss_3_4.clear();
	ss_5_6.clear();
	ss_7_8.clear();
	ss_12_34.clear();
	ss_56_78.clear();
	ss_1234_5678.clear();

	std::ofstream ofs;
	if( nDataset == 1 ){
		ofs.open ("ABCstat_8pop.txt", std::ofstream::out);

		ofs << "sx_1_avg\tsx_1_std\t";
		ofs << "sx_2_avg\tsx_2_std\t";
		ofs << "sx_3_avg\tsx_3_std\t";
		ofs << "sx_4_avg\tsx_4_std\t";
		ofs << "sx_5_avg\tsx_5_std\t";
		ofs << "sx_6_avg\tsx_6_std\t";
		ofs << "sx_7_avg\tsx_7_std\t";
		ofs << "sx_8_avg\tsx_8_std\t";
		ofs << "sx_12_avg\tsx_12_std\t";
		ofs << "sx_34_avg\tsx_34_std\t";
		ofs << "sx_56_avg\tsx_56_std\t";
		ofs << "sx_78_avg\tsx_78_std\t";
		ofs << "sx_1234_avg\tx_1234_std\t";
		ofs << "sx_5678_avg\tx_5678_std\t";
		
		ofs << "sf_1_2_avg\tsf_1_2_std\t";
		ofs << "sf_3_4_avg\tsf_3_4_std\t";
		ofs << "sf_5_6_avg\tsf_5_6_std\t";
		ofs << "sf_7_8_avg\tsf_7_8_std\t";
		ofs << "sf_12_34_avg\tsf_12_34_std\t";
		ofs << "sf_56_78_avg\tsf_56_78_std\t";
		ofs << "sf_1234_5678_avg\tsf_1234_5678_std\t";
		
		ofs << "ss_1_2_avg\tss_1_2_std\t";
		ofs << "ss_3_4_avg\tss_3_4_std\t";
		ofs << "ss_5_6_avg\tss_5_6_std\t";
		ofs << "ss_7_8_avg\tss_7_8_std\t";
		ofs << "ss_12_34_avg\tss_12_34_std\t";
		ofs << "ss_56_78_avg\tss_56_78_std\t";
		ofs << "ss_1234_5678_avg\tss_1234_5678_std" << std::endl;
		
		ofs.close();
	}
	
	ofs.open ("ABCstat_8pop.txt", std::ofstream::out | std::ofstream::app);
	ofs << mean_sx_1 << "\t" << sd_sx_1 << "\t";
	ofs << mean_sx_2 << "\t" << sd_sx_2 << "\t";
	ofs << mean_sx_3 << "\t" << sd_sx_3 << "\t";
	ofs << mean_sx_4 << "\t" << sd_sx_4 << "\t";
	ofs << mean_sx_5 << "\t" << sd_sx_5 << "\t";
	ofs << mean_sx_6 << "\t" << sd_sx_6 << "\t";
	ofs << mean_sx_7 << "\t" << sd_sx_7 << "\t";
	ofs << mean_sx_8 << "\t" << sd_sx_8 << "\t";
	ofs << mean_sx_12 << "\t" << sd_sx_12 << "\t";
	ofs << mean_sx_34 << "\t" << sd_sx_34 << "\t";
	ofs << mean_sx_56 << "\t" << sd_sx_56 << "\t";
	ofs << mean_sx_78 << "\t" << sd_sx_78 << "\t";
	ofs << mean_sx_1234 << "\t" << sd_sx_1234 << "\t";
	ofs << mean_sx_5678 << "\t" << sd_sx_5678 << "\t";
	
	ofs << mean_sf_1_2 << "\t" << sd_sf_1_2 << "\t";
	ofs << mean_sf_3_4 << "\t" << sd_sf_3_4 << "\t";
	ofs << mean_sf_5_6 << "\t" << sd_sf_5_6 << "\t";
	ofs << mean_sf_7_8 << "\t" << sd_sf_7_8 << "\t";
	ofs << mean_sf_12_34 << "\t" << sd_sf_12_34 << "\t";
	ofs << mean_sf_56_78 << "\t" << sd_sf_56_78 << "\t";
	ofs << mean_sf_1234_5678 << "\t" << sd_sf_1234_5678 << "\t";
	
	ofs << mean_ss_1_2 << "\t" << sd_ss_1_2 << "\t";
	ofs << mean_ss_3_4 << "\t" << sd_ss_3_4 << "\t";
	ofs << mean_ss_5_6 << "\t" << sd_ss_5_6 << "\t";
	ofs << mean_ss_7_8 << "\t" << sd_ss_7_8 << "\t";
	ofs << mean_ss_12_34 << "\t" << sd_ss_12_34 << "\t";
	ofs << mean_ss_56_78 << "\t" << sd_ss_56_78 << "\t";
	ofs << mean_ss_1234_5678 << "\t" << sd_ss_1234_5678 << std::endl;
	ofs.close();
}

void mean_deviation(const std::vector <unsigned int> & tableau, const unsigned int & nLoci, float & mean, float & sd){
	size_t i(0);
	
	// compute the mean
	float mean_tmp(0.0);
	for(i=0; i<nLoci; ++i){
		mean_tmp += tableau[i];
	}
	mean_tmp /= (1.0 * nLoci);

	// compute the sd
	float sd_tmp(0.0);
	for(i=0; i<nLoci; ++i){
		sd_tmp += pow(tableau[i] - mean_tmp, 2);
	}
	mean = mean_tmp;
	sd = sqrt(sd_tmp/(1.0 * nLoci));
}

