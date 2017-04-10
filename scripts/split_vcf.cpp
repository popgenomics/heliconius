#include <sstream>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

std::string getContigName(const std::string & line);

int main(int argc, char* argv[]){
	if(argc != 4){
		std::cerr << std::endl;
		std::cerr << "\033[1;31m ERROR: 3 arguments are needed\033[0m\n";
		std::cerr << "\t arg1 = name of the VCF file to split" << std::endl;
		std::cerr << "\t arg2 = species name to print in output files (single word: short but recognizable)" << std::endl;
		std::cerr << "\t arg3 = number of contigs to put in each sub-VCF files" << std::endl << std::endl;
		std::cerr << "\033[1;33m example: ./split_vcf ama10.Hmel2.bwa.default.HC.vcf ama 10\033[0m\n" << std::endl << std::endl;
		exit(0);
	}

	const std::string inputFile(argv[1]); // file to read
	const std::string speciesName(argv[2]); // variable used for the name of the outputFile
	const unsigned int groupSize = atoi(argv[3]); // variable used for the name of the outputFile

	std::ifstream fifo(inputFile.c_str());
	if(fifo){
		std::string line; // contains the readen line

		std::string contigName; // contains the first word of a row
		std::string previousContigName("");
		unsigned int test(0);
		unsigned int nContig(0);
		unsigned int nGroups = 0; 
	
		std::string header = "";

		// declaration of the first outputFile: subVCF_{speciesName}_0.vcf
		std::ostringstream oss;
		oss << "subVCF_" << speciesName << "_" << nGroups << ".vcf"; 
		std::string outfileName = oss.str();
		std::ofstream file;
		file.open(outfileName, std::ios::out);
	
		// declaration of the correspondance_table file: corr_table_{speciesName}.txt
		std::ostringstream oss2;
		oss2 << "corr_table_" << speciesName << ".txt"; 
		std::string outfileName2 = oss2.str();
		std::ofstream file2;
		file2.open(outfileName2, std::ios::out);
				
		size_t cnt = 0;

		while(std::getline(fifo, line)){ // read the input file

			if(line[0]=='#'){
				if(line[1]=='C'){
					header = line; // stores the header of table, which will be common for all outputFiles

					file << header << std::endl; // writes the header in the first outputFile
				}
			}else{ // if the line doesn't start by a '#'
				contigName = getContigName(line); // gets the first word of a row of the table 
				
				if(test==0){ // if it's the first row of the table:
					test=1;
					file2 << contigName << "\t" << outfileName << std::endl;
					previousContigName = contigName; 
				}else{ // it the first row of the table had been already treated
					if(contigName != previousContigName){ // if the current row has a different contigName than the previous
						++nContig; // one contig because it's a new one
						previousContigName = contigName;

						if(nContig%groupSize==0){ // every groupSize contigs:
							++nGroups;
//							std::cout << "close " << outfileName << std::endl;
							file.close(); // close the former outputFile
							
							outfileName = "";

							oss.str(std::string());
							oss << "subVCF_" << speciesName << "_" << nGroups << ".vcf";

							outfileName = oss.str();
//							std::ofstream file;
							file.open(outfileName, std::ios::out);
							file << header << std::endl; // writes the header in the new outputFile

						}
						file2 << contigName << "\t" << outfileName << std::endl;
					}
				}

				if(file){
					file << line << std::endl; // writes the current line in the outputFile
				}else{	
					std::cerr << std::endl << "ERROR: the file " << outfileName << " was not found" << std::endl;
					exit(0);
				}
			}
			
		} // end of VCF file
	}else{
		exit(0);
	}
	return(0);
}
std::string getContigName(const std::string & line){
	std::string word;
	std::istringstream iss(line);
	size_t i(0);
	while( std::getline(iss, word, '\t') ){
		if( i==0 ){ // if word is the first item of the line:
			return(word);
		}
		++i;
	}
}

