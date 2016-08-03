#include <string>
#include <map>
#include <iomanip>
#include <iostream>
#include <algorithm>

using namespace std;

#include "ls_seq.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"

// for each character in the alphabet 'seq_chars_'
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
    if (finished_) {
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
            << ((total_[seq_chars_[i]] / divide) * 100.0) << endl;
    }
    if (!is_protein_) {
    (*poos) << left << setw(nameWidth) << setfill(separator) << "G+C"
        << (total_['G'] + total_['C']) << "\t"
        << (((total_['G'] + total_['C']) / divide) * 100.0) << endl;
    }
    (*poos) << "--------" << seq_type_ << " TABLE---------" << endl;
}

// transpose original atbel
void SeqInfo::print_stats_alt (ostream* poos) {
    
    const char separator = ' ';
    const int nameWidth = 10;
    double divide = 0.0;
    if (is_protein_) {
        seq_type_ = "Prot";
    } else {
        seq_type_ = "Nucl";
    }
    if (finished_) {
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
            << ((total_[seq_chars_[i]] / divide) * 100.0) << endl;
    }
    if (!is_protein_) {
    (*poos) << left << setw(nameWidth) << setfill(separator) << "G+C"
        << (total_['G'] + total_['C']) << "\t"
        << (((total_['G'] + total_['C']) / divide) * 100.0) << endl;
    }
    (*poos) << "--------" << seq_type_ << " TABLE---------" << endl;
}

// just grab labels, disregard the rest
void SeqInfo::collect_taxon_labels () {
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*pios_, retstring);
    while (read_next_seq_from_stream(*pios_, ft, retstring, seq)) {
        name_ = seq.get_id();
        taxon_labels_.push_back(name_);
    }
    if (ft == 2) {
        name_ = seq.get_id();
        taxon_labels_.push_back(name_);
    }
    sort(taxon_labels_.begin(), taxon_labels_.end());
}

// assumed aligned if all seqs are the same length
void SeqInfo::check_is_aligned () {
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*pios_, retstring);
    while (read_next_seq_from_stream(*pios_, ft, retstring, seq)) {
        int terp = seq.get_sequence().size();
        seq_lengths_.push_back(terp);
    }
    if (ft == 2) {
        int terp = seq.get_sequence().size();
        seq_lengths_.push_back(terp);
    }
    // check if all seqs are the same length
    if (std::adjacent_find( seq_lengths_.begin(), seq_lengths_.end(), std::not_equal_to<int>()) == seq_lengths_.end() ) {
        is_aligned_ = true;
    } else {
        is_aligned_ = false;
    }
}

void SeqInfo::get_nseqs () {
    seqcount_ = 0;
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*pios_, retstring);
    while (read_next_seq_from_stream(*pios_, ft, retstring, seq)) {
        seqcount_++;
    }
    if (ft == 2) {
        seqcount_++;
    }
    sort(taxon_labels_.begin(), taxon_labels_.end());
}

void SeqInfo::get_nchars () {
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*pios_, retstring);
    while (read_next_seq_from_stream(*pios_, ft, retstring, seq)) {
        int terp = seq.get_sequence().size();
        name_ = seq.get_id();
        taxon_labels_.push_back(name_);
        seq_lengths_.push_back(terp);
    }
    if (ft == 2) {
        int terp = seq.get_sequence().size();
        name_ = seq.get_id();
        taxon_labels_.push_back(name_);
        seq_lengths_.push_back(terp);
    }
    // check if all seqs are the same length
    if (std::adjacent_find( seq_lengths_.begin(), seq_lengths_.end(), std::not_equal_to<int>()) == seq_lengths_.end() ) {
        is_aligned_ = true;
        seq_length_ = seq_lengths_[0];
    } else {
        is_aligned_ = false;
        seq_length_ = -1;
    }
    seqcount_ = (int)seq_lengths_.size();
}

// get the longest label. for printing purposes
void SeqInfo::get_longest_taxon_label () {
    longest_tax_label_ = 0;
    for (int i = 0; i < seqcount_; i++) {
        if (taxon_labels_[i].size() > longest_tax_label_) {
            longest_tax_label_ = taxon_labels_[i].size();
        }
    }
}

void SeqInfo::set_alphabet () {
    if (is_protein_) {
        seq_chars_ = "ACDEFGHIKLMNPQRSTVWXY*";
    } else {
        seq_chars_ = "ACGTN-";
    }
}

SeqInfo::SeqInfo (istream* pios, ostream* poos, bool& indiv, bool const& force_protein) {
    // set parameters
    output_indiv_ = (indiv == true) ? true : false;
    if (force_protein) {
        is_protein_ = true;
    }
    pios_ = pios;
    poos_ = poos;
}

// return wichever property set to true
void SeqInfo::get_property (bool const& get_labels, bool const& check_aligned,
        bool const& get_nseq, bool const& get_freqs, bool const& get_nchar) {
    
    if (get_labels) {
        collect_taxon_labels();
        for (unsigned int i = 0; i < taxon_labels_.size(); i++) {
            (*poos_) << taxon_labels_[i] << endl;
        }
    } else if (check_aligned) {
        check_is_aligned();
        (*poos_) << std::boolalpha << is_aligned_ << endl;
    } else if (get_nseq) {
        get_nseqs ();
        (*poos_) << seqcount_ << endl;
    } else if (get_freqs) {
        // use original code
        
    } else if (get_nchar) {
        get_nchars ();
        if (!output_indiv_) { // single return value
            if (seq_length_ != -1) {
                (*poos_) << seq_length_ << endl;
            } else {
                // not aligned
                (*poos_) << "sequences are not aligned" << endl;
            }
        } else { // individual lengths
            get_longest_taxon_label();
            for (int i = 0; i < seqcount_; i++) {
                int diff = longest_tax_label_ - taxon_labels_[i].size();
                (*poos_) << taxon_labels_[i];
                if (diff > 0) {
                    string pad = std::string(diff, ' ');
                    (*poos_) << pad;
                }
                (*poos_) << " " << seq_lengths_[i] << endl;
            }
        }
        
    }
}

void SeqInfo::summarize () {

    //Concatenated will be used for all stats
    finished_ = false;
    seqcount_ = 0;
    
    bool first = true;
    
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*pios_, retstring);
    
    file_type_ = get_filetype_string(ft);
    
    while (read_next_seq_from_stream(*pios_, ft, retstring, seq)) {
        if (first) {
            // infer sequence type rather than setting it
            if (!is_protein_) {
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
            print_stats(poos_);
        }
    }
    if (ft == 2) {
        seqcount_++;
        concatenated_ += seq.get_sequence();
        temp_seq_ = seq.get_sequence();
        name_ = seq.get_id();
        if (output_indiv_) {
            count_chars_indiv_seq(temp_seq_);
            print_stats(poos_);
        }
    }
    finished_ = true;
    count_chars_indiv_seq(concatenated_);
    print_stats(poos_);
}


SeqInfo::SeqInfo (istream* pios, bool& indiv, bool const& force_protein, ostream* poos) {

    //Concatenated will be used for all stats
    finished_ = false;
    seqcount_ = 0;
    output_indiv_ = (indiv == true) ? true : false;
    
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
