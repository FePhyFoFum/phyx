#include <string>
#include <map>
#include <iomanip>
#include <iostream>
#include <algorithm>

using namespace std;

#include "seq_info.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"

// for each character in the alphabet 'seq_chars_'
void SeqInfo::count_chars_indiv_seq(string& seq) {
    
    seq = string_to_upper(seq);
    
    total_.clear();
    for (unsigned int i = 0; i < seq_chars_.length(); i++) {
        total_[seq_chars_[i]] = 0.0;
    }
    for (unsigned int i = 0; i < seq.length(); i++) {
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

// alternate to above, accumulate char counts across seqs
void SeqInfo::count_chars (string& seq) {
    unsigned int sum = 0;
    
    seq = string_to_upper(seq);
    
    if (output_indiv_) {
        vector <int> icounts(seq_chars_.length(), 0);
        
        for (unsigned int i = 0; i < seq_chars_.length(); i++) {
            int num = count(seq.begin(), seq.end(), seq_chars_[i]);
            char_counts_[i] += num;
            icounts[i] += num;
            sum += num;
        }
        // add invalid char counts, add to missing char count
        if (sum < seq.length()) {
            char_counts_[char_counts_.size() - 1] += (seq.length() - sum);
            icounts[icounts.size() - 1] += (seq.length() - sum);
        }
        indiv_char_counts_.push_back(icounts);
        // this is cool, but unnecessary here
        //std::transform(char_counts_.begin(), char_counts_.end(), icounts.begin(), char_counts_.begin(), std::plus<int>());
        
    } else {
        for (unsigned int i = 0; i < seq_chars_.length(); i++) {
            int num = count(seq.begin(), seq.end(), seq_chars_[i]);
            char_counts_[i] += num;
            sum += num;
        }
        // add invalid char counts, add to missing char count
        if (sum < seq.length()) {
            char_counts_[char_counts_.size() - 1] += (seq.length() - sum);
        }
    }
}

// calculate character state frequencies
void SeqInfo::calculate_freqs () {
    bool first = true;
    Sequence seq;
    string retstring;
    seqcount_ = 0;
    int ft = test_seq_filetype_stream(*pios_, retstring);
    file_type_ = get_filetype_string(ft);
    while (read_next_seq_from_stream(*pios_, ft, retstring, seq)) {
        if (first) {
            if (!is_protein_) {
                string alpha_name = seq.get_alpha_name();
                if (alpha_name == "AA") {
                    is_protein_ = true;
                }
            }
            set_alphabet ();
            first = false;
        }
        seqcount_++;
        temp_seq_ = seq.get_sequence();
        name_ = seq.get_id();
        seq_lengths_.push_back(temp_seq_.length());
        count_chars(temp_seq_);
        taxon_labels_.push_back(name_);
    }
    if (ft == 2) {
        seqcount_++;
        temp_seq_ = seq.get_sequence();
        name_ = seq.get_id();
        seq_lengths_.push_back(temp_seq_.length());
        count_chars(temp_seq_);
        taxon_labels_.push_back(name_);
    }
}

// alt to print_stats. essential difference is transposed results
void SeqInfo::return_freq_table (ostream* poos) {
    const char separator = ' ';
    const int colWidth = 10;
    if (output_indiv_) {
        // need to take into account longest_tax_label_
        get_longest_taxon_label();
        string pad = std::string(longest_tax_label_, ' ');
        // header
        (*poos) << pad << " ";
        for (unsigned int i = 0; i < seq_chars_.length(); i++) {
            (*poos) << right << setw(colWidth) << setfill(separator)
                << seq_chars_[i] << " ";
        }
        // return nchar for individual seqs
        (*poos) << right << setw(colWidth) << setfill(separator) << "Nchar" << endl;
        for (int i = 0; i < seqcount_; i++) {
            int diff = longest_tax_label_ - taxon_labels_[i].size();
            (*poos_) << taxon_labels_[i];
            if (diff > 0) {
                pad = std::string(diff, ' ');
                (*poos_) << pad;
            }
            (*poos_) << " ";
            for (unsigned int j = 0; j < seq_chars_.length(); j++) {
                (*poos) << right << setw(colWidth) << setfill(separator)
                    << (double)indiv_char_counts_[i][j] / (double)seq_lengths_[i] << " ";
            }
            (*poos) << right << setw(colWidth) << setfill(separator) << seq_lengths_[i] << endl;
        }
    } else {
        // header
        for (unsigned int i = 0; i < seq_chars_.length(); i++) {
            (*poos) << right << setw(colWidth) << setfill(separator)
                << seq_chars_[i];
            if (i != seq_chars_.length() - 1) {
                (*poos) << " ";
            }
        }
        (*poos) << endl;
        // counts
        for (unsigned int i = 0; i < seq_chars_.length(); i++) {
            (*poos) << right << setw(colWidth) << setfill(separator)
                << char_counts_[i];
            if (i != seq_chars_.length() - 1) {
                (*poos) << " ";
            }
        }
        (*poos) << endl;
        // freqs
        int total_num_chars = sum(char_counts_);
        for (unsigned int i = 0; i < seq_chars_.length(); i++) {
            (*poos) << fixed << right << setw(colWidth) << setfill(separator)
                << (double)char_counts_[i] / (double)total_num_chars;
            if (i != seq_chars_.length() - 1) {
                (*poos) << " ";
            }
        }
        (*poos) << endl;
    }
}

void SeqInfo::print_stats (ostream* poos) {
    const char separator = ' ';
    const int colWidth = 10;
    double divide = 0.0;
    if (is_protein_) {
        seq_type_ = "Prot";
    } else {
        seq_type_ = "Nucl";
    }
    
    //(*poos) << "General Stats For All Sequences" << endl;
    (*poos) << "File type: " << file_type_ << endl;
    (*poos) << "Number of sequences: " << seqcount_ << endl;
    if (std::adjacent_find( seq_lengths_.begin(), seq_lengths_.end(), std::not_equal_to<int>()) == seq_lengths_.end() ) {
        is_aligned_ = true;
    } else {
        is_aligned_ = false;
    }
    (*poos_) << "Is aligned: " << std::boolalpha << is_aligned_ << endl;
    if (is_aligned_) {
        seq_length_ = seq_lengths_[0];
        (*poos_) << "Sequence length: " << seq_length_ << endl;
    }
    //(*poos) << "Total Length of All Combined: " << concatenated_.length() << endl; // not really useful, is it?
    divide = concatenated_.length();
    
    (*poos) << "--------" << seq_type_ << " TABLE---------" << endl;
    (*poos) << left << setw(6) << setfill(separator) << seq_type_ << " "
        << setw(colWidth) << setfill(separator) << "Total" << " "
        << setw(colWidth) << setfill(separator) << "Proportion" << endl;
    for (unsigned int i = 0; i < seq_chars_.length(); i++) {
        (*poos) << left << setw(6) << setfill(separator) << seq_chars_[i] << " "
            << setw(colWidth) << setfill(separator) << total_[seq_chars_[i]] << " "
            << ((total_[seq_chars_[i]] / divide)) << endl;
    }
    if (!is_protein_) {
        (*poos) << left << setw(6) << setfill(separator) << "G+C" << " "
            << setw(colWidth) << setfill(separator) << (total_['G'] + total_['C']) << " "
            << (((total_['G'] + total_['C']) / divide)) << endl;
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

// does not currently do per-individual...
void SeqInfo::calc_missing () {
    calculate_freqs();
    // missing data are the last two characters (-N for DNA, -X for protein)
    int total_num_chars = sum(char_counts_);
    double temp = 0.0;
    for (unsigned int i = seq_chars_.length()-2; i < seq_chars_.length(); i++) {
        temp += (double)char_counts_[i] / (double)total_num_chars;
    }
    percent_missing_ = temp;
}

// get the longest label. for printing purposes
void SeqInfo::get_longest_taxon_label () {
    longest_tax_label_ = 0;
    for (int i = 0; i < seqcount_; i++) {
        if ((int)taxon_labels_[i].size() > longest_tax_label_) {
            longest_tax_label_ = taxon_labels_[i].size();
        }
    }
}

void SeqInfo::set_alphabet () {
    if (is_protein_) {
        seq_chars_ = "ACDEFGHIKLMNPQRSTVWY-X";
    } else {
        seq_chars_ = "ACGT-N";
    }
    char_counts_.resize(seq_chars_.size(), 0);
}

SeqInfo::SeqInfo (istream* pios, ostream* poos, bool& indiv, bool const& force_protein) {
    // set parameters
    output_indiv_ = (indiv == true) ? true : false;
    is_protein_ = false;
    if (force_protein) {
        is_protein_ = true;
    }
    pios_ = pios;
    poos_ = poos;
}

// return whichever property set to true
void SeqInfo::get_property (bool const& get_labels, bool const& check_aligned,
        bool const& get_nseq, bool const& get_freqs, bool const& get_nchar,
        double const& get_missing) {
    
    if (get_labels) {
        collect_taxon_labels();
        for (unsigned int i = 0; i < taxon_labels_.size(); i++) {
            (*poos_) << taxon_labels_[i] << endl;
        }
    } else if (check_aligned) {
        check_is_aligned();
        (*poos_) << std::boolalpha << is_aligned_ << endl;
    } else if (get_nseq) {
        get_nseqs();
        (*poos_) << seqcount_ << endl;
    } else if (get_freqs) {
        calculate_freqs();
        return_freq_table(poos_);
    } else if (get_missing) {
        calc_missing();
        (*poos_) << percent_missing_ << endl;
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

// TODO: get rid of the concatenated_ business
void SeqInfo::summarize () {
    // Concatenated will be used for all stats
    seqcount_ = 0;
    
    bool first = true;
    
    Sequence seq;
    string retstring;
    int ft = test_seq_filetype_stream(*pios_, retstring);
    file_type_ = get_filetype_string(ft);
    
    while (read_next_seq_from_stream(*pios_, ft, retstring, seq)) {
        if (first) {
            if (!is_protein_) { // if not forced to be protein, test
                string alpha_name = seq.get_alpha_name();
                if (alpha_name == "AA") {
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
        seq_lengths_.push_back(temp_seq_.length());
        count_chars(temp_seq_);
        taxon_labels_.push_back(name_);
    }
    if (ft == 2) {
        seqcount_++;
        concatenated_ += seq.get_sequence();
        temp_seq_ = seq.get_sequence();
        name_ = seq.get_id();
        seq_lengths_.push_back(temp_seq_.length());
        count_chars(temp_seq_);
        taxon_labels_.push_back(name_);
    }
    
    if (output_indiv_) {
        // new one
        return_freq_table(poos_);
    } else {
        count_chars_indiv_seq(concatenated_);
        print_stats(poos_);
    }
}
