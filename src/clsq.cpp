/*
 * clsq.cpp
 *
 *  Created on: Jun 15, 2015
 *      Author: joe
 */


#include "clsq.h"
#include "sequence.h"
#include "seq_reader.h"


SequenceCleaner::SequenceCleaner(istream* pios, double& missing):numTaxa (0), 
        numChar (0), missingAllowed(missing) {
    read_sequences (pios); // read in sequences on initialization
    clean_sequences ();
}

void SequenceCleaner::read_sequences (istream* pios) {
    
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*pios, retstring);
    // should put error-checking here e.g. check that sequences are all the same length
    // not doing that here; just assuming things are cool
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        sequences[seq.get_id()] = seq.get_sequence();
        if (numChar == 0) { // just getting this from an arbitrary sequence for now
            numChar = (int)seq.get_sequence().size();
        }
    }
    if (ft == 2) {
		sequences[seq.get_id()] = seq.get_sequence();
    }
    numTaxa = sequences.size();
}
 // not used
int SequenceCleaner::get_num_taxa () {
    return numTaxa;
}

// not used
map<string, string> SequenceCleaner::get_trimmed_seqs () {
    return trimmedSeqs;
}

void SequenceCleaner::write_seqs (ostream* poos) {
    //cout << "About to write sequences for " << trimmedSeqs.size() << " taxa!" << endl;
    // hmm. doesn't work...
    //for (auto& kv : trimmedSeqs) {
    //    (*poos) << ">" << kv.first << endl;
    //    (*poos) << kv.second << endl;
    //}
    for (iter = trimmedSeqs.begin(); iter != trimmedSeqs.end(); iter++) {
        (*poos) << ">" << iter->first << endl;
        (*poos) << iter->second << endl;
    }
}

void SequenceCleaner::clean_sequences () {
    
    /*
    ifstream readline;
    double NumbOfSequences = 0.0;
    string new_dna;
    bool round_one = true;
    unsigned int StillMissing = 0;
    readline.open(fasta.c_str());
    if (readline.is_open()){
        while (getline (readline, line)) {
            if (line[0] == '>') {
                if (round_one == false) {
                    line.erase (std::remove(line.begin(), line.end(), '>'), line.end());
                    sequences[name_hold] = dna;
                    dna = "";
                } else {
                    line.erase (std::remove(line.begin(), line.end(), '>'), line.end());
                }
                name_hold = line;
            } else {
                round_one = false;
                dna += line;
            }
        }
    }
    sequences[name_hold] = dna;
    int size = dna.size();
    */
    
    double MissingData[numChar];
    double PercentMissingData[numChar];
    
    string new_dna;
    unsigned int stillMissing = 0;
    
    for (iter = sequences.begin(); iter != sequences.end(); iter++) {

        new_dna = iter -> second;
        CheckMissing(MissingData, new_dna);
        //NumbOfSequences++;
    }
    for (iter = sequences.begin(); iter != sequences.end(); iter++) {

        string to_stay = "";
        new_dna = iter -> second;
        stillMissing = 0;
        for (unsigned int i = 0; i < new_dna.size(); i++) {
            PercentMissingData[i] =  MissingData[i] / (double)numTaxa;
            if (PercentMissingData[i] >= missingAllowed) {
                
                // *** something missing here? ***
                
            } else {
                to_stay += new_dna[i];
            }
        }
        stillMissing = 0;
        for (unsigned int j = 0; j < to_stay.size(); j++) {
            if (to_stay[j] == 'N' || to_stay[j] == '-' ||  to_stay[j] == 'n'
                ||  to_stay[j] == 'X' || to_stay[j] == 'x') {
                stillMissing += 1;
            }
        }
        if (stillMissing == to_stay.size()) {
            cout << "Removed: " << iter -> first << endl;
        } else {
            trimmedSeqs[iter -> first] = to_stay;
        }
    }
    //return trimmedSeqs;
}
void SequenceCleaner::CheckMissing(double MissingData [], string& dna){

    for (int i = 0; i < numChar; i++) {
        // use tolower
        if (dna[i] == 'N' || dna[i] == '-' || dna[i] == 'n' || 
            dna[i] == 'X' ||  dna[i] == 'x'){
            MissingData[i]++;
        }
    }
}

SequenceCleaner::~SequenceCleaner() {
    // TODO Auto-generated destructor stub
}

