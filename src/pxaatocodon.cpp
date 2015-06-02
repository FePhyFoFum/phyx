/*
 * pxaatocodon.cpp
 *
 *  Created on: Jun 2, 2015
 *      Author: joe
 */


//Basic Program to change an amino acid alignment
//to codons using the unaligned nucleotide fasta

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <map>
using namespace std;

map <string, string> FastaToOneLine (string fasta){

        map <string, string> sequences;
        string line, dna, name_hold;
        ifstream readline;
        float count = 0;
        readline.open(fasta.c_str());
        if (readline.is_open()){
                while (getline (readline, line)){
                        if (line[0] == '>'){
                                if (count != 0){
                                    line.erase (std::remove(line.begin(), line.end(), '>'), line.end());
                                    sequences[name_hold] = dna;
                                    dna = "";

                                }else{
                                    line.erase (std::remove(line.begin(), line.end(), '>'), line.end());
                                }
                                name_hold = line;
                        }else{
                                count = 1;
                                dna += line;

                        }
                }
        }
        sequences[name_hold] = dna;

        return sequences;
}
map <string, string> ChangeToCodon(map<string, string> AminoAcid, map<string,string> Nucleotide){

	map<string, string> CodonSequences;
	map<string, string>::iterator iter;
	string AminoAcidSequence;
	string NucleotideSequence;
	string temp;
	string::iterator it;
	int j;
    for (iter = AminoAcid.begin(); iter != AminoAcid.end(); iter++){

    	if (Nucleotide.find(iter -> first) == Nucleotide.end()){

    		cout << "Only in the AA File: " << iter -> first << endl;
    	}else{

    		AminoAcidSequence = iter -> second;
    		NucleotideSequence = Nucleotide[iter -> first];
    		for(int i = 0; i < AminoAcidSequence.size(); i++){

    			j = 0;
    			if (AminoAcidSequence[i] == '-'){

    				temp += "---";
    			}else{

    				temp += NucleotideSequence[0];
    				temp += NucleotideSequence[1];
    				temp += NucleotideSequence[2];
    				NucleotideSequence.erase(0,2);
    			}
    		}
    		CodonSequences[iter -> first] = temp;
    		temp = "";
    	}
    }


	return CodonSequences;
}

int main(){

	map<string, string> AASequences;
	map<string, string> NucSequences;
	map<string, string> CodonSequences;
	map<string, string>::iterator iter;
	string AAFasta;
	string NucFasta;
	//cout << "Amino Acid File" << endl;
	//cin >> AAFasta;
	// TestFiles/AA.fa
	AAFasta = ("TestFiles/AA.fa");
	//cout << "Unaligned Nucleotide File" << endl;
	//cin >> NucFasta;
	// TestFiles/Codon.fa
	NucFasta = ("TestFiles/Codon.fa");
	AASequences = FastaToOneLine(AAFasta);
	NucSequences = FastaToOneLine(NucFasta);
	CodonSequences = ChangeToCodon(AASequences, NucSequences);
    for (iter = CodonSequences.begin(); iter != CodonSequences.end(); iter++){
    	cout << ">" << iter -> first << "\n" << iter -> second << endl;
    }
}

