/*
 * clsq.cpp
 *
 *  Created on: Jun 15, 2015
 *      Author: joe
 */


#include "clsq.h"
#include "sequence.h"
#include "seq_reader.h"


SequenceCleaner::SequenceCleaner(istream* pios, double& missing):num_taxa_(0), 
        num_char_(0.0), missing_allowed_(missing) {
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
        sequences_[seq.get_id()] = seq.get_sequence();
        if (num_char_ == 0.0) { // just getting this from an arbitrary (first) sequence for now
            num_char_ = (double)seq.get_sequence().size();
        }
    }
    if (ft == 2) {
        sequences_[seq.get_id()] = seq.get_sequence();
    }
    num_taxa_ = sequences_.size();
}
 // not used
int SequenceCleaner::get_num_taxa () {
    return num_taxa_;
}

// not used
map<string, string> SequenceCleaner::get_trimmed_seqs () {
    return trimmed_seqs_;
}

void SequenceCleaner::write_seqs (ostream* poos) {
    if (trimmed_seqs_.size() == 0) {
        for (iter_ = sequences_.begin(); iter_ != sequences_.end(); iter_++) {
            (*poos) << ">" << iter_->first << endl;
            (*poos) << "-" << endl;
        }
    }
    for (iter_ = trimmed_seqs_.begin(); iter_ != trimmed_seqs_.end(); iter_++) {
        (*poos) << ">" << iter_->first << endl;
        (*poos) << iter_->second << endl;
    }
}

void SequenceCleaner::clean_sequences () {
    
    double MissingData[num_char_];
    double PercentMissingData[num_char_];
    for (int i = 0; i < num_char_; i++) {
        MissingData[i] = 0.0;
    }
    string new_dna;
    unsigned int stillMissing = 0;
    
    for (iter_ = sequences_.begin(); iter_ != sequences_.end(); iter_++) {
        new_dna = iter_ -> second;
        CheckMissing(MissingData, new_dna);
        //NumbOfSequences++;
    }
    for (iter_ = sequences_.begin(); iter_ != sequences_.end(); iter_++) {

        string to_stay = "";
        new_dna = iter_ -> second;
        stillMissing = 0;
        for (unsigned int i = 0; i < new_dna.size(); i++) {
            PercentMissingData[i] =  MissingData[i] / (double)num_taxa_;
            //cout << "Position: " << i << "Amount Missing: " << MissingData[i] << " Percent Missing: " << PercentMissingData[i] << "Number of Taxa: " <<  (double)numTaxa  << " Allowed Missing: " << missingAllowed << endl;
            if (PercentMissingData[i] > missing_allowed_) {
                
                
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
            cout << "Removed: " << iter_ -> first << endl;
        } else {
            trimmed_seqs_[iter_ -> first] = to_stay;
        }
    }
}

void SequenceCleaner::CheckMissing(double MissingData [], string& dna) {

    for (int i = 0; i < num_char_; i++) {
        // use tolower
        if (dna[i] == 'N' || dna[i] == '-' || dna[i] == 'n' || 
            dna[i] == 'X' ||  dna[i] == 'x') {
            MissingData[i]++;
            //cout << "Position: " << i << " DNA: " << dna[i] <<  " Missing: " << MissingData[i] << endl;
        }
    }
}

SequenceCleaner::~SequenceCleaner() {
    // TODO Auto-generated destructor stub
}

