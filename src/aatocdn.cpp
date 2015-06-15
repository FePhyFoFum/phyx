/*
 * aatocdn.cpp
 *
 *  Created on: Jun 15, 2015
 *      Author: joe
 */

#include "aatocdn.h"

map<string, string> aa_to_cdn::ChangeToCodon(map<string, string>& AminoAcid, map<string,string>& Nucleotide){

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


	return CodonAln;
}

//Used Twice so everything is declared in the function
map <string, string> aa_to_cdn::FastaToOneLine (string& fasta){

    map <string, string> sequences;
    string line, dna, name_hold;
    ifstream readline;
	bool round_one = true;
	readline.open(fasta.c_str());
	if (readline.is_open()){
        while (getline (readline, line)){
            if (line[0] == '>'){

                if (round_one == false){

                	line.erase(std::remove(line.begin(), line.end(), '>'), line.end());
                    sequences[name_hold] = dna;
                    dna = "";

                }else{
                    line.erase(std::remove(line.begin(), line.end(), '>'), line.end());
                }
                name_hold = line;
            }else{
                round_one = false;
                dna += line;
            }
        }
    }
    sequences[name_hold] = dna;
    return sequences;

}


aa_to_cdn::aa_to_cdn() {
	// TODO Auto-generated constructor stub

}

aa_to_cdn::~aa_to_cdn() {
	// TODO Auto-generated destructor stub
}

