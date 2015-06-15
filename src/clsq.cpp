/*
 * clsq.cpp
 *
 *  Created on: Jun 15, 2015
 *      Author: joe
 */


#include "clsq.h"

void CheckMissing(double MissingData [], string& dna){

	for (int i = 0; i < dna.size(); i++){

		if (dna[i] == 'N' || dna[i] == '-' ||  dna[i] == 'n' ||  dna[i] == 'X' ||  dna[i] == 'x'){

			MissingData[i]++;
		}
	}
}

map<string, string> clsq::FastaToOneLine (string& fasta, double& MissingAllowed){

    ifstream readline;
    double NumbOfSequences;
	string new_dna;
    bool round_one = true;
    int StillMissing = 0;
    readline.open(fasta.c_str());
    if (readline.is_open()){
        while (getline (readline, line)){
            if (line[0] == '>'){

                if (round_one == false){
                    line.erase (std::remove(line.begin(), line.end(), '>'), line.end());
                    sequences[name_hold] = dna;
                    dna = "";

                }else{
                    line.erase (std::remove(line.begin(), line.end(), '>'), line.end());
                }
                name_hold = line;
            }else{
                round_one = false;
                dna += line;
            }
        }
    }
    sequences[name_hold] = dna;
    int size = dna.size();
    double MissingData[size];
    double PercentMissingData[size];
    for(iter = sequences.begin(); iter != sequences.end(); iter++){

    	new_dna = iter -> second;
    	CheckMissing(MissingData, new_dna);
    	NumbOfSequences++;
    }
    for(iter = sequences.begin(); iter != sequences.end(); iter++){

    	string to_stay = "";
    	new_dna = iter -> second;
    	StillMissing = 0;
    	for (int i = 0; i < new_dna.size(); i++){

    		PercentMissingData[i] =  MissingData[i] / NumbOfSequences;
    		if(PercentMissingData[i] >= MissingAllowed){

    		}else{

    			to_stay += new_dna[i];
    		}

    	}
        StillMissing = 0;
		for (int j = 0; j < to_stay.size(); j++){

			if (to_stay[j] == 'N' ||to_stay[j] == '-' ||  to_stay[j] == 'n' ||  to_stay[j] == 'X' || to_stay[j] == 'x'){
				StillMissing += 1;
			}
		}
    	if (StillMissing == to_stay.size()){

    		cout << "Removed: " << iter -> first << endl;

    	}else{
    	    Trimmed[iter -> first] = to_stay;

    	}

    }
    return Trimmed;
}



clsq::clsq() {
	// TODO Auto-generated constructor stub

}

clsq::~clsq() {
	// TODO Auto-generated destructor stub
}

