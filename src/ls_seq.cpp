#include <string>
#include <map>
#include <iomanip>
#include <iostream>

using namespace std;

#include "ls_seq.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"

void Stats::STAT_Getter(string& seq, bool& prot) {

    total_.clear();
    for (unsigned int i = 0; i < seq_chars_.length(); i++) {
        total_[seq_chars_[i]] = 0.0;
    }
    for (unsigned int i = 0; i < seq.length(); i++) {
        seq[i] = toupper(seq[i]);
        //Ensure there is no weird J or whatever characters
        if (total_.find(seq[i]) == total_.end()) {
            if (prot == true){
                total_['X']++;
            } else {
                total_['N']++;
            }
        } else {
            total_[seq[i]]++;
        }
    }
}

void Stats::Printer (bool& prot, ostream* poos) {
    
        const char separator = ' ';
        const int nameWidth = 10;
        double divide = 0.0;
        if (prot == true){
            seq_type_ = "Prot ";
        }else{
            seq_type_ = "Nucl ";
        }
        if (finished_ == true) {
            (*poos) << "General Stats For All Sequences" << endl;
            (*poos) << "File Type: " << type_ << endl;
            (*poos) << "Number of Sequences: " << seqcount_ << endl;
            (*poos) << "Total Length of All Combined: " << concatenated_.length() << endl;
            divide = concatenated_.length();
        } else {
            (*poos) << "General Stats For " << name_ << endl;
            (*poos) << "Total Length: " << temp_seq_.length() << endl;    
            divide = temp_seq_.length();
        }
        (*poos) << "--------" << seq_type_ << "TABLE---------" << endl;
        (*poos) << seq_type_ << "\tTotal\tPercent" << endl;
        for (unsigned int i = 0; i < seq_chars_.length(); i++) {
            (*poos) << left << setw(nameWidth) << setfill(separator) << seq_chars_[i]
                << total_[seq_chars_[i]] << "\t"
                << ((total_[seq_chars_[i]] / divide)*100.0) << endl;
        }
        if (prot == false) {
        (*poos) << left << setw(nameWidth) << setfill(separator) << "G+C"
            << (total_['G'] + total_['C']) << "\t"
            << (((total_['G'] + total_['C']) / divide)*100.0) << endl;
    }        
    (*poos) << "--------" << seq_type_ << "TABLE---------" << endl;

}

Stats::Stats (istream* pios, bool& all, bool& prot, ostream* poos) {

    //Concatenated will be used for all stats
    finished_ = false;
    seqcount_ = 0;
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*pios, retstring);
    if (prot == true) {
        seq_chars_ = "ACDEFGHIKLMNPQRSTVWXY*";
    } else {
        seq_chars_ = "ACGTN-";
    }
    while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
        
        seqcount_++;
        concatenated_ += seq.get_sequence();
        temp_seq_ = seq.get_sequence();
        name_ = seq.get_id();
        if (all == true) {
            STAT_Getter(temp_seq_, prot);
            Printer(prot, poos);
        }
        if (ft == 1) {
            type_ = "Phylip";
        }
        if (ft == 0) {
            type_ = "Nexus";
        }
    }
    if (ft == 2) {
        seqcount_++;
        concatenated_ += seq.get_sequence();
        temp_seq_ = seq.get_sequence();
        name_ = seq.get_id();
        type_ = "Fasta";
        if (all == true) {
            STAT_Getter(temp_seq_, prot);
            Printer(prot, poos);
        }
    }
    finished_ = true;
    STAT_Getter(concatenated_, prot);
    Printer(prot, poos);
}
