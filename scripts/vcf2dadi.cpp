#include <sstream>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>


// g++ -std=c++17 vcf2dadi.cpp -O3 -o vcf2dadi


void treatVCFfile(std::string const vcfFileSPA, const unsigned int minCov, std::vector<std::string> & alleleAlt_spA, std::vector<unsigned int> & nAlleleRef_spA, std::vector<unsigned int> & nAlleleAlt_spA, std::vector<std::string> & contig, std::vector<std::string> & position, std::vector<std::string> & alleleRef, std::vector<std::string> & format);

int main(int argc, char* argv[]){
	if(argc != 8){
		std::cerr << std::endl;
		std::cerr << "\033[1;31m ERROR: 7 arguments are needed\033[0m\n";
		std::cerr << "\t arg1 = name of the VCF file for species A" << std::endl;
		std::cerr << "\t arg2 = name of the VCF file for species B" << std::endl;
		std::cerr << "\t arg3 = minimum number of reads (>= nReads) to call an allelic state" << std::endl;
		std::cerr << "\t arg4 = minimum number of called allelic state (>= nRef+nAlt) in a species to keep the SNP" << std::endl;
		std::cerr << "\t arg5 = name of species A (ex: ama)" << std::endl;
		std::cerr << "\t arg6 = name of species B (ex: chi)" << std::endl;
		std::cerr << "\t arg7 = label to specify the output files" << std::endl << std::endl;
		std::cerr << "\t camille.roux.1983@gmail.com; 2017" << std::endl << std::endl;
		std::cerr << "\033[1;33m example: ./vcf2dadi subVCF_ama_105.vcf subVCF_chi_105.vcf 10 6 ama chi 105\033[0m\n" << std::endl << std::endl;
		exit(0);
	}

	const std::string vcfFileSPA(argv[1]); // vcf's name for species A
	const std::string vcfFileSPB(argv[2]); // vcf's name for species B
	const unsigned int minCov(atoi(argv[3])); // minimum number of reads to call an allelic state
	const unsigned int minAllState(atoi(argv[4])); // minimum number of allelic states to keep a SNP
	const std::string speciesA(argv[5]); // name of species A
	const std::string speciesB(argv[6]); // name of species B
	const std::string labelOutput(argv[7]); // label put in output file name for being recognized

	size_t i(0);

	// species A
	std::vector<std::string> contig_spA; // vector of contig names
	std::vector<std::string> position_spA; // vector of positions in a contig
	std::vector<std::string> alleleRef_spA; // vector of alleles of the reference genome. INDELS + MULTIALLELISM ARE CURRENTLY MASKED
	std::vector<std::string> alleleAlt_spA; // vector of alternative alleles. INDELS + MULTIALLELISM ARE CURRENTLY MASKED
	std::vector<std::string> format_spA; // vector of VCF code describing the format of read depths
	std::vector<unsigned int> nAlleleRef_spA; // vector of number of copies of the reference allele
	std::vector<unsigned int> nAlleleAlt_spA; // vector of number of copies of the alternative allele
	
	// read vcf file for species A
	treatVCFfile(vcfFileSPA, minCov, alleleAlt_spA, nAlleleRef_spA, nAlleleAlt_spA, contig_spA, position_spA, alleleRef_spA, format_spA);

	// species B 
	std::vector<std::string> contig_spB;
	std::vector<std::string> position_spB;
	std::vector<std::string> alleleRef_spB;
	std::vector<std::string> format_spB;
	std::vector<std::string> alleleAlt_spB;
	std::vector<unsigned int> nAlleleRef_spB;
	std::vector<unsigned int> nAlleleAlt_spB;

	// read vcf file for species B
	treatVCFfile(vcfFileSPB, minCov, alleleAlt_spB, nAlleleRef_spB, nAlleleAlt_spB, contig_spB, position_spB, alleleRef_spB, format_spB);


	// write the results
		std::ofstream outputFlux("output_" + speciesA + "_" + speciesB + "_" + labelOutput + ".txt", std::ios::out);
		if(outputFlux){
			outputFlux << "Ing\tOut\tAllele1\t" << speciesA << "\t" << speciesB << "\tAllele2\t" << speciesA << "\t" << speciesB << "\tGene\tPosition\n";
//			outputFlux << "contig\tposition\tallele_Ref\tallele_Alt_spA\tnAllRef_spA\tnAllAlt_spA\tallele_Alt_spB\tnAllRef_spB\tnAllAlt_spB\n";
			for(i=0; i<contig_spA.size(); ++i){
				if(alleleRef_spA[i]!="N" && alleleRef_spA[i]!="." && alleleAlt_spA[i]!="." && alleleRef_spB[i]!="N" && alleleRef_spB[i]!="." && alleleAlt_spB[i]!="."){
					if( nAlleleAlt_spA[i]!=0 || nAlleleAlt_spB[i]!=0 ){
						if(  nAlleleRef_spA[i]+nAlleleAlt_spA[i]>=minAllState && nAlleleRef_spB[i]+nAlleleAlt_spB[i]>=minAllState ){
							if( alleleAlt_spA[i]==alleleAlt_spB[i] ){
								outputFlux << "-" << alleleRef_spA[i] << "-\t---\t" << alleleRef_spA[i] << "\t";
								outputFlux << nAlleleRef_spA[i] << "\t" << nAlleleRef_spB[i] << "\t";
								outputFlux << alleleAlt_spA[i] << "\t";
								outputFlux << nAlleleAlt_spA[i] << "\t" << nAlleleAlt_spB[i] << "\t";
								outputFlux << contig_spA[i] << "\t" << position_spA[i] << "\n";
							}
						}
					}
				}
			}
		}else{
			std::cerr <<  "ERROR: cannot oppen output.xt" << std::endl;
			exit(0);
		}

	return(0);
}

void treatVCFfile(std::string const vcfFileSPA, const unsigned int minCov, std::vector<std::string> & alleleAlt_spA, std::vector<unsigned int> & nAlleleRef_spA, std::vector<unsigned int> & nAlleleAlt_spA, std::vector<std::string> & contig, std::vector<std::string> & position, std::vector<std::string> & alleleRef, std::vector<std::string> & format){
	std::string contigTMP;
	std::string alleleRefTMP;
	std::string alleleAltTMP;
	std::string formatTMP; // contains format of individual informations in vcf
	std::string positionTMP;
	
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

	std::ifstream vcfA(vcfFileSPA.c_str(), std::ios::in);
	if(vcfA){
		while(getline(vcfA, line)){
			if(line[0]!='#'){ // start of loop over positions
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

/*		std::ofstream outputFlux("output.txt", std::ios::out);
		if(outputFlux){
			outputFlux << "contig\tposition\tformat\tallele_Ref\tallele_Alt_spA\tnAllRef_spA\tnAllAlt_spA\n";
			for(i=0; i<contig.size(); ++i){
				outputFlux << contig[i] << "\t" << position[i] << "\t" << format[i] << "\t" << alleleRef[i] << "\t" << alleleAlt_spA[i] << "\t" << nAlleleRef_spA[i] << "\t" << nAlleleAlt_spA[i] << std::endl;
			}
		}else{
			std::cerr <<  "ERROR: cannot oppen output.xt" << std::endl;
			exit(0);
		}
*/
	}else{
		std::cerr << "ERROR: file " << vcfFileSPA << " is not found" << std::endl;
	}
}

