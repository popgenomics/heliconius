#include <sstream>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

int main(int argc, char* argv[]){
	const std::string vcfFileSPA(argv[1]); // vcf's name for species A

	const unsigned int minCov(10);

	std::string contigTMP;
	std::string alleleRefTMP;
	std::string alleleAltTMP;
	std::string formatTMP; // contains format of individual informations in vcf
	std::string positionTMP;
	
	std::vector<std::string> contig;
	std::vector<std::string> position;
	std::vector<std::string> alleleRef;
	std::vector<std::string> format;

	// species A
	std::vector<std::string> alleleAlt_spA;
	std::vector<unsigned int> nAlleleRef_spA;
	std::vector<unsigned int> nAlleleAlt_spA;
	
	// species B 
	std::vector<std::string> alleleAlt_spB;
	std::vector<unsigned int> nAlleleRef_spB;
	std::vector<unsigned int> nAlleleAlt_spB;

	std::string line;
	std::string word; // used to parse a line
	std::string word2; // used to parse a complex word 
	std::string word3; // used to parse AD (allele depth)

	int i(0); // used in loops
	int j(0); // used in sub-loops
	int k(0); // used in sub-sub-loops
	int test_AD(0); // =0: no alternative allele; =1: alternative allele
	int nReads(0);
	int nReadsRef(0);
	int nReadsAlt(0);


	//int DP(0); // indicates the position of the read depth
	int AD(0); // indicates the position of the allelic depths (ref,alt)
	int nRefAll_spA_TMP(0); // number of individuals with the reference allele in species A for a given position
	int nRefAll_spB_TMP(0); // number of individuals with the reference allele in species B for a given position
	int nAltAll_spA_TMP(0); // number of individuals with the alternative allele in species A for a given position
	int nAltAll_spB_TMP(0); // number of individuals with the alternative allele in species B for a given position

	size_t found(0); // used to test the presence of some expression in a string using std::string::find

	// read vcf file for species A
	std::ifstream vcfA(vcfFileSPA.c_str(), std::ios::in);
	if(vcfA){
		while(getline(vcfA, line)){
			if(line[0]!='#' & line[1]!='#'){ // start of loop over positions
				i = -1;
				std::istringstream iss(line);
				while(std::getline(iss, word, '\t') ){ // start of loop along a vcf line
					++i;
					if( i==0 ){ // column = contig
						contigTMP = word;
						continue;
					}
					if( i==1 ){ // column = position
						positionTMP = word;
						continue;
					}
					if( i==3 ){ // column = reference allele
						alleleRefTMP = word;
						continue;
					}
					if( i==4 ){ // column = alternative allele 
						alleleAltTMP = word;
						continue;
					}
					if( i==8 ){ // column = alternative allele 
						formatTMP = word;

						nRefAll_spA_TMP = 0;
						nAltAll_spA_TMP = 0;
						
						test_AD = 0;
						found = 0;
						found = formatTMP.find("AD");
						if(found != std::string::npos){
							test_AD = 1; // there is an alternative allele
//							std::cout << "Pos: " << positionTMP << "\t2 alternative alleles: " << alleleRefTMP << " and " << alleleAltTMP << std::endl;
						}else{
							test_AD = 0; // it's only a matter of reference allele
						}
						continue;
					}
					
					if( alleleRefTMP.size()==1 && alleleAltTMP.size()==1 ){ // if ref and alt alleles are not a nucleotide -> alleleTMP='.' and nAll=0 
						if( i>8 ){ // start of loop over individuals
							nReads = 0;
							j=-1;
							std::istringstream iss2(word);
							while(std::getline(iss2, word2, ':') ){ // start of loop within individuals
								++j;

								// if no alternative allele:
								if( test_AD==0 ){
									if( j==1 ){ // gets the information about coverage
										nReads = atoi(word2.c_str());
										if( nReads<minCov ){
											nRefAll_spA_TMP+=0;
										}

										if( nReads>=minCov && nReads<(2*minCov) ){ // only call one allele if: minCov <= nReads < 2minCov
											nRefAll_spA_TMP+=1;
										}

										if( nReads>=(2*minCov) ){ // call two alleles if: nReads >= 2*minCov
											nRefAll_spA_TMP+=2;
										}
									}
								}

								// if there is an alternative allele:
								if( test_AD==1 ){
									if( j==1 ){ // gets the information about coverage
										nReadsRef = 0;
										nReadsAlt = 0;
										k=-1;
										std::istringstream iss3(word2);
										while(std::getline(iss3, word3, ',') ){ // loop over alleles depth
											++k;
											if( k==0 ){ // ref allele: in case of alleleRefTMP.size()==1 && alleleAltTMP.size()==1
												nReadsRef = atoi(word3.c_str());
											}
											if( k==1 ){ // alt allele: in case of alleleRefTMP.size()==1 && alleleAltTMP.size()==1
												nReadsAlt = atoi(word3.c_str());
											}
										}
									// start of treatment of nReads
									if( nReadsRef<minCov && nReadsAlt<minCov ){
											nRefAll_spA_TMP+=0;
											nAltAll_spA_TMP+=0;
										}
									if( nReadsRef>=minCov && nReadsAlt>=minCov ){
											nRefAll_spA_TMP+=1;
											nAltAll_spA_TMP+=1;
										}
									if( nReadsRef<minCov && nReadsAlt>=minCov ){
											nRefAll_spA_TMP+=0;
											if( nReadsAlt<(2*minCov) ){
												nAltAll_spA_TMP+=1;
											}else{
												nAltAll_spA_TMP+=2;
											}
										}
									if( nReadsRef>=minCov && nReadsAlt<minCov ){
											nAltAll_spA_TMP+=0;
											if( nReadsRef<(2*minCov) ){
												nRefAll_spA_TMP+=1;
											}else{
												nRefAll_spA_TMP+=2;
											}
										}
									// end of treatment of nReads
									}
								}
								
							} // end of loop within individuals
						} // end of loop over individuals
					} // end of loop of the condition: if a alleleRefTMP.size()==1 and allelleAltTMP.size()==1 (only consider biallelic variation, no indel)
				if( alleleRefTMP.size()!=1 & alleleAltTMP.size()!=1 ){ // if ref or alt allele are: neither biallelic, or contain indel 
					alleleRefTMP="N";
					alleleAltTMP=".";
					nRefAll_spA_TMP=0;
					nAltAll_spA_TMP=0;
				}

				} // end of loop along a vcf line
//			std::cout << positionTMP << std::endl;
			contig.push_back(contigTMP);
			position.push_back(positionTMP);
			alleleRef.push_back(alleleRefTMP);
			alleleAlt_spA.push_back(alleleAltTMP);
			format.push_back(formatTMP);
			nAlleleRef_spA.push_back(nRefAll_spA_TMP);
			nAlleleAlt_spA.push_back(nAltAll_spA_TMP);
			} // end of loop over positions
		} // end of VCF file

		vcfA.close();
/*		std::cout << contig.size() << " contigs" << std::endl;
		std::cout << position.size() << " positions" << std::endl;
		std::cout << alleleRef.size() << " ref alleles" << std::endl;
		std::cout << alleleAlt_spA.size() << " alt alleles" << std::endl;
		std::cout << format.size() << " formats" << std::endl;
		std::cout << nAlleleRef_spA.size() << " counts of ref alleles" << std::endl;
		std::cout << nAlleleAlt_spA.size() << " counts of alt alleles" << std::endl;*/

		std::ofstream outputFlux("output.txt", std::ios::out);
		if(outputFlux){
			outputFlux << "contig\tposition\tformat\tallele_Ref\tallele_Alt_spA\tnAllRef_spA\tnAllAlt_spA\n";
			for(i=0; i<contig.size(); ++i){
				outputFlux << contig[i] << "\t" << position[i] << "\t" << format[i] << "\t" << alleleRef[i] << "\t" << alleleAlt_spA[i] << "\t" << nAlleleRef_spA[i] << "\t" << nAlleleAlt_spA[i] << std::endl;
			}
		}else{
			std::cerr <<  "ERROR: cannot oppen output.xt" << std::endl;
			exit(0);
		}

	}else{
		std::cerr << "ERROR: file " << vcfFileSPA << " is not found" << std::endl;
	}
	return(0);
}

