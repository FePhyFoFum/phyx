#include <string>
#include <map>
#include <iomanip>
#include <iostream>

using namespace std;

#include "ls_seq.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"

// for each character in the alphabet
void SeqInfo::count_chars_indiv_seq(string& seq) {

    total_.clear();
    for (unsigned int i = 0; i < seq_chars_.length(); i++) {
        total_[seq_chars_[i]] = 0.0;
    }
    for (unsigned int i = 0; i < seq.length(); i++) {
        seq[i] = toupper(seq[i]);
        // Ensure there is no weird J or whatever characters (includes '?')
        if (total_.find(seq[i]) == total_.end()) {
            if (is_protein_) {
                total_['X']++;
            } else {
                total_['N']++;
            }
        } else {
            total_[seq[i]]++;
        }
    }
}

void SeqInfo::print_stats (ostream* poos) {
    
        const char separator = ' ';
        const int nameWidth = 10;
        double divide = 0.0;
        if (is_protein_) {
            seq_type_ = "Prot";
        } else {
            seq_type_ = "Nucl";
        }
        if (finished_ == true) {
            (*poos) << "General Stats For All Sequences" << endl;
            (*poos) << "File type: " << file_type_ << endl;
            (*poos) << "Number of sequences: " << seqcount_ << endl;
            (*poos) << "Total Length of All Combined: " << concatenated_.length() << endl;
            divide = concatenated_.length();
        } else {
            (*poos) << "General Stats For " << name_ << endl;
            (*poos) << "Total Length: " << temp_seq_.length() << endl;    
            divide = temp_seq_.length();
        }
        (*poos) << "--------" << seq_type_ << " TABLE---------" << endl;
        (*poos) << seq_type_ << "\tTotal\tPercent" << endl;
        for (unsigned int i = 0; i < seq_chars_.length(); i++) {
            (*poos) << left << setw(nameWidth) << setfill(separator) << seq_chars_[i]
                << total_[seq_chars_[i]] << "\t"
                << ((total_[seq_chars_[i]] / divide)*100.0) << endl;
        }
        if (!is_protein_) {
        (*poos) << left << setw(nameWidth) << setfill(separator) << "G+C"
            << (total_['G'] + total_['C']) << "\t"
            << (((total_['G'] + total_['C']) / divide)*100.0) << endl;
    }        
    (*poos) << "--------" << seq_type_ << " TABLE---------" << endl;
}

void SeqInfo::collect_taxon_labels () {
    
}

void SeqInfo::set_alphabet () {
    if (is_protein_) {
        seq_chars_ = "ACDEFGHIKLMNPQRSTVWXY*";
    } else {
        seq_chars_ = "ACGTN-";
    }
}

SeqInfo::SeqInfo (istream* pios, bool& all, bool const& force_protein, ostream* poos) {

    //Concatenated will be used for all stats
    finished_ = false;
    seqcount_ = 0;
    output_indiv_ = (all == true) ? true : false;
    
    bool first = true;
    
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*pios, retstring);
    
    file_type_ = get_filetype_string(ft);
    
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        if (first) {
            // infer sequence type rather than setting it
            if (force_protein) {
                is_protein_ = true;
            } else {
                string alpha_name = seq.get_alpha_name();
                if (alpha_name == "AA") {
                    //cout << "I believe this is: " << alpha_name << "!" << endl;
                    is_protein_ = true;
                }
            }
            set_alphabet ();
            
            first = false;
        }
        seqcount_++;
        concatenated_ += seq.get_sequence();
        temp_seq_ = seq.get_sequence();
        name_ = seq.get_id();
        if (output_indiv_) {
            count_chars_indiv_seq(temp_seq_);
            print_stats(poos);
        }
    }
    if (ft == 2) {
        seqcount_++;
        concatenated_ += seq.get_sequence();
        temp_seq_ = seq.get_sequence();
        name_ = seq.get_id();
        if (output_indiv_) {
            count_chars_indiv_seq(temp_seq_);
            print_stats(poos);
        }
    }
    finished_ = true;
    count_chars_indiv_seq(concatenated_);
    print_stats(poos);
}
