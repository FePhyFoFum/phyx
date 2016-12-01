/*
 * clsq.cpp
 *
 *  Created on: Jun 15, 2015
 *      Author: joe
 */

#include <iostream>

using namespace std;

#include "clsq.h"
#include "sequence.h"
#include "seq_reader.h"

SequenceCleaner::SequenceCleaner(istream* pios, double& proportion, bool& force_protein,
        bool const& verbose):num_taxa_(0), num_char_(0), required_present_(proportion) {
    //cout << MolDna << endl;
    missing_allowed_ = 1.0 - required_present_;
    is_dna_ = !force_protein;
    verbose_ = verbose;
    read_sequences (pios); // read in sequences on initialization
    clean_sequences ();
}

void SequenceCleaner::read_sequences (istream* pios) {
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*pios, retstring);
    int num_current_char = 0;
    bool first = true;
    
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        sequences_[seq.get_id()] = seq.get_sequence();
        num_current_char = seq.get_sequence().size();
        if (first) {
            num_char_ = num_current_char; // just getting this from an arbitrary (first) sequence for now
            if (is_dna_) {
                string alpha_name = seq.get_alpha_name();
                if (alpha_name == "AA") {
                    is_dna_ = false;
                    //cout << "I believe this is a protein!" << endl;
                }
            }
            first = false;
            continue;
        } else {
            if (num_current_char != num_char_) {
                cout << "Error: sequences are not all of the same length. Exiting."
                    << endl;
                exit(0);
            }
        }
    }
    if (ft == 2) {
        sequences_[seq.get_id()] = seq.get_sequence();
        num_current_char = seq.get_sequence().size();
        if (num_current_char != num_char_) {
            cout << "Error: sequences are not all of the same length. Exiting."
                << endl;
            exit(0);
        }
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
        CheckMissing(MissingData, new_dna, is_dna_);
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
        if (is_dna_) {
            for (unsigned int j = 0; j < to_stay.size(); j++) {
                if (to_stay[j] == 'N' || to_stay[j] == '-' ||  to_stay[j] == 'n'
                    ||  to_stay[j] == 'X' || to_stay[j] == 'x') {
                    stillMissing += 1;
                }
            }
        }else{
            for (unsigned int j = 0; j < to_stay.size(); j++) {
                if (to_stay[j] == '-' ||  to_stay[j] == 'X' || to_stay[j] == 'x') {
                    stillMissing += 1;
                }
            }
            
        }
        if (stillMissing == to_stay.size()) {
            if (verbose_) {
                cout << "Removed: " << iter_ -> first << endl;
            }
        } else {
            trimmed_seqs_[iter_ -> first] = to_stay;
        }
    }
}

void SequenceCleaner::CheckMissing(double MissingData [], string& dna, bool& type) {

    if (type == true){
        for (int i = 0; i < num_char_; i++) {
            if (tolower(dna[i]) == 'n' || dna[i] == '-' || tolower(dna[i]) == 'x') {
                MissingData[i]++;
                //cout << "Position: " << i << " DNA: " << dna[i] <<  " Missing: " << MissingData[i] << endl;
            }
        }
    }else{
        for (int i = 0; i < num_char_; i++) {
            if (dna[i] == '-' || tolower(dna[i]) == 'x') {
                MissingData[i]++;
                //cout << "Position: " << i << " DNA: " << dna[i] <<  " Missing: " << MissingData[i] << endl;
            }
        }
    }
}

SequenceCleaner::~SequenceCleaner() {
    // TODO Auto-generated destructor stub
}

