#include <string>
#include <map>
#include <iomanip>
#include <iostream>
#include <algorithm>

#include "seq_info.h"
#include "utils.h"
#include "seq_utils.h"
#include "sequence.h"
#include "seq_reader.h"


SeqInfo::SeqInfo (std::istream* pios, std::ostream* poos, bool& indiv,
        const bool& force_protein):concatenated_(""), seq_chars_(""), output_indiv_(indiv),
        datatype_set_(false), is_dna_(false), is_protein_(false), is_multi_(false),
        is_binary_(false), alpha_set_(false), alpha_name_(""), seq_type_(""), gap_('-'),
        missing_('?'), num_taxa_(0), percent_missing_(0.0), is_aligned_(false),
        seq_length_(0), longest_tax_label_(0) {
    // maybe get rid of this? how often is inference wrong?
    if (force_protein) {
        is_protein_ = true;
        datatype_set_ = true;
        set_alphabet();
    }
    pios_ = pios;
    poos_ = poos;
    read_in_alignment();
}


// for each character in the alphabet 'seq_chars_'
void SeqInfo::count_chars_indiv_seq (std::string& seq) {
    seq = string_to_upper(seq);
    total_.clear(); // probably unnecessary
    for (unsigned int i = 0; i < seq_chars_.length(); i++) {
        total_[seq_chars_[i]] = 0.0;
    }
    for (unsigned int i = 0; i < seq.length(); i++) {
        // Ensure there is no weird J or whatever characters (includes '?')
        if (total_.find(seq[i]) == total_.end()) {
            total_[missing_]++;
        } else {
            total_[seq[i]]++;
        }
    }
}


// alternate to above, accumulate char counts across seqs
void SeqInfo::count_chars (std::string& seq) {
    unsigned int sum = 0;
    seq = string_to_upper(seq);
    if (output_indiv_) {
        std::vector<int> icounts(seq_chars_.length(), 0);
        for (unsigned int i = 0; i < seq_chars_.length(); i++) {
            int num = std::count(seq.begin(), seq.end(), seq_chars_[i]);
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
            int num = std::count(seq.begin(), seq.end(), seq_chars_[i]);
            char_counts_[i] += num;
            sum += num;
        }
        // add invalid char counts, add to missing char count
        if (sum < seq.length()) {
            char_counts_[char_counts_.size() - 1] += (seq.length() - sum);
        }
    }
}


// moving this outside, since so many functions need it
// read everything in, process afterwards
// this means that all seqs reside in memory. oh well...
// seqs are stored, datatype and alphabet is set
void SeqInfo::read_in_alignment () {
    Sequence seq;
    std::string retstring;
    int ft = test_seq_filetype_stream(*pios_, retstring);
    int file_ntax = 0; // ntax declared in the file itself
    file_type_ = get_filetype_string(ft);
    
    // if nexus, grab metadata
    if (file_type_ == "nexus") {
        //std::cout << "Trying to read in a Nexus alignment..." << std::endl;
        bool is_interleaved = false;
        get_nexus_alignment_properties(*pios_, file_ntax, seq_length_,
                is_interleaved, alpha_name_, seq_chars_, gap_, missing_);
        // std::cout << "alpha_name_ = " << alpha_name_ << std::endl;
        set_datatype();
        if (is_multi_ && seq_chars_ != "") {
            seq_chars_ += gap_;
            seq_chars_ += missing_;
            char_counts_.resize(seq_chars_.size(), 0);
            alpha_set_ = true;
        }
        retstring = ""; // have to set so seq_reader knows we are mid-file
        if (!is_interleaved) {
            while (read_next_seq_from_stream(*pios_, ft, retstring, seq)) {
                seqs_.push_back(seq);
            }
        } else {
            // need to read in everything at once
            seqs_ = read_interleaved_nexus(*pios_, file_ntax, seq_length_);
        }
    } else {
        bool complicated_phylip = false;
        // check if we are dealing with a complicated phylip format
        if (file_type_ == "phylip") {
            get_phylip_dimensions(retstring, file_ntax, seq_length_);
            complicated_phylip = is_complicated_phylip(*pios_, seq_length_);
        }
        if (complicated_phylip) {
            seqs_ = read_phylip(*pios_, file_ntax, seq_length_);
            if (!datatype_set_) {
                alpha_name_ = seqs_[0].get_alpha_name();
                set_datatype();
            }
        } else {
            while (read_next_seq_from_stream(*pios_, ft, retstring, seq)) {
                seqs_.push_back(seq);
                if (!datatype_set_) {
                    alpha_name_ = seq.get_alpha_name();
                    set_datatype();
                }
            }
            if (ft == 2) { // fasta has an trailing one
                seqs_.push_back(seq);
            }
        }
    }
    
    // now do the counting
    if (is_multi_ && !alpha_set_) {
        // grab all unique characters from the input string
        // here, seqs from all individuals are concatenated, so represents all sampled characters
        make_concatenated_sequence();
        set_alphabet_from_sampled_seqs(concatenated_);
    }
    if (file_ntax != 0) {
        if (file_ntax != (int)seqs_.size()) {
            std::cerr << "Error: number of taxa declared in the file ("
                << ") does not match the number read (" << seqs_.size()
                << "). Exiting." << std::endl;
            exit(1);
        }
    }
    num_taxa_ = (int)seqs_.size();
}


// calculate character state frequencies
void SeqInfo::calculate_freqs () {
    Sequence seq;
    std::string name;
    for (unsigned int i = 0; i < seqs_.size(); i++) {
        seq = seqs_[i];
        temp_seq_ = seq.get_sequence();
        name = seq.get_id();
        seq_lengths_.push_back(temp_seq_.length());
        count_chars(temp_seq_);
        taxon_labels_.push_back(name);
    }
}


// alt to print_summary_table_whole_alignment. essential difference is transposed results
void SeqInfo::return_freq_table () {
    const int colWidth = 12;
    if (output_indiv_) {
        // need to take into account longest_tax_label_
        longest_tax_label_ = get_longest_label(taxon_labels_);
        std::string pad = std::string(longest_tax_label_, ' ');
        // header
        (*poos_) << pad << " ";
        for (unsigned int i = 0; i < seq_chars_.length(); i++) {
            (*poos_) << std::right << std::setw(colWidth) << seq_chars_[i] << " ";
        }
        // return nchar for individual seqs
        (*poos_) << std::right << std::setw(colWidth) << "Nchar" << std::endl;
        for (int i = 0; i < num_taxa_; i++) {
            int diff = longest_tax_label_ - taxon_labels_[i].size();
            (*poos_) << taxon_labels_[i];
            if (diff > 0) {
                pad = std::string(diff, ' ');
                (*poos_) << pad;
            }
            (*poos_) << " ";
            for (unsigned int j = 0; j < seq_chars_.length(); j++) {
                (*poos_) << std::right << std::setw(colWidth)
                    << (double)indiv_char_counts_[i][j] / (double)seq_lengths_[i] << " ";
            }
            (*poos_) << std::right << std::setw(colWidth) << seq_lengths_[i] << std::endl;
        }
    } else {
        // header
        for (unsigned int i = 0; i < seq_chars_.length(); i++) {
            (*poos_) << std::right << std::setw(colWidth) << seq_chars_[i];
            if (i != seq_chars_.length() - 1) {
                (*poos_) << " ";
            }
        }
        (*poos_) << std::endl;
        // counts
        for (unsigned int i = 0; i < seq_chars_.length(); i++) {
            (*poos_) << std::right << std::setw(colWidth) << char_counts_[i];
            if (i != seq_chars_.length() - 1) {
                (*poos_) << " ";
            }
        }
        (*poos_) << std::endl;
        // freqs
        int total_num_chars = sum(char_counts_);
        for (unsigned int i = 0; i < seq_chars_.length(); i++) {
            (*poos_) << std::fixed << std::right << std::setw(colWidth)
                << (double)char_counts_[i] / (double)total_num_chars;
            if (i != seq_chars_.length() - 1) {
                (*poos_) << " ";
            }
        }
        (*poos_) << std::endl;
    }
}


void SeqInfo::print_summary_table_whole_alignment () {
    const int colWidth = 12;
    double total_num_chars = 0.0;
    
    //(*poos) << "General Stats For All Sequences" << std::endl;
    (*poos_) << "File type: " << file_type_ << std::endl;
    (*poos_) << "Number of sequences: " << num_taxa_ << std::endl;
    if (std::adjacent_find( seq_lengths_.begin(), seq_lengths_.end(), std::not_equal_to<int>()) == seq_lengths_.end() ) {
        is_aligned_ = true;
    } else {
        is_aligned_ = false;
    }
    (*poos_) << "Is aligned: " << std::boolalpha << is_aligned_ << std::endl;
    if (is_aligned_) {
        seq_length_ = seq_lengths_[0];
        (*poos_) << "Sequence length: " << seq_length_ << std::endl;
        total_num_chars = (double)(seq_lengths_[0] * num_taxa_);
    } else {
        total_num_chars = (double)sum(seq_lengths_);
    }
    
    (*poos_) << "--------- " << seq_type_ << " TABLE ----------" << std::endl;
    (*poos_) << std::right << std::setw(4) << seq_type_ << " "
        << std::setw(colWidth) << "Total" << " "
        << std::setw(colWidth) << "Proportion" << std::endl;
    for (unsigned int i = 0; i < seq_chars_.length(); i++) {
        (*poos_) << std::right << std::setw(4) << seq_chars_[i] << " "
            << std::setw(colWidth) << (int)total_[seq_chars_[i]] << " "
            << std::setw(colWidth) << ((total_[seq_chars_[i]] / total_num_chars)) << std::endl;
    }
    if (is_dna_) {
        (*poos_) << std::right << std::setw(4) << "G+C" << " "
            << std::setw(colWidth) << (int)(total_['G'] + total_['C']) << " "
            << std::setw(colWidth) << (((total_['G'] + total_['C']) / total_num_chars)) << std::endl;
    }
    (*poos_) << "--------- " << seq_type_ << " TABLE ----------" << std::endl;
}


void SeqInfo::make_concatenated_sequence () {
    if (concatenated_.length() == 0) {
        for (unsigned int i = 0; i < seqs_.size(); i++) {
            concatenated_ += seqs_[i].get_sequence();
        }
    }
}


// just grab labels, disregard the rest
void SeqInfo::collect_taxon_labels () {
    taxon_labels_ = collect_names(seqs_);
}


// assumed aligned if all seqs are the same length
void SeqInfo::check_is_aligned () {
    is_aligned_ = is_aligned(seqs_);
}


void SeqInfo::get_num_chars () {
    for (unsigned int i = 0; i < seqs_.size(); i++) {
        seq_lengths_.push_back(seqs_[i].get_length());
    }
    // check if all seqs are the same length
    if (std::adjacent_find( seq_lengths_.begin(), seq_lengths_.end(), std::not_equal_to<int>()) == seq_lengths_.end() ) {
        is_aligned_ = true;
        seq_length_ = seq_lengths_[0];
    } else {
        is_aligned_ = false;
        seq_length_ = -1;
    }
}


// does not currently do per-individual...
// assumes gap=- and missing=?
void SeqInfo::calc_missing () {
    // missing data are the last two characters (-N for DNA, -X for protein)
    int miss = 0;
    double temp = 0.0;
    
    calculate_freqs();
    
    if (!output_indiv_) {
        // proportion for alignment as a whole
        int total_num_chars = sum(char_counts_);
        for (unsigned int i = seq_chars_.length()-2; i < seq_chars_.length(); i++) {
            temp += (double)char_counts_[i] / (double)total_num_chars;
            miss += char_counts_[i];
        }
        percent_missing_ = temp;
        //std::cout << "total_num_chars = " << total_num_chars << std::endl;
        //std::cout << "total missing = " << miss << std::endl;
    } else {
        // per individual sequence
        // should we require seqs be aligned? stats would make more sense
        Sequence seq;
        std::string name;
        for (unsigned int i = 0; i < seqs_.size(); i++) {
            seq = seqs_[i];
            temp_seq_ = seq.get_sequence();
            name = seq.get_id();
            taxon_labels_.push_back(name);
            seq_lengths_.push_back(temp_seq_.length());
            count_chars(temp_seq_);
            miss = 0;
            for (unsigned int j = seq_chars_.length()-2; j < seq_chars_.length(); j++) {
                miss += indiv_char_counts_[i][j];
            }
            missing_counts_.push_back(miss);
        }
    }
}


void SeqInfo::return_missing () {
    if (!output_indiv_) {
        (*poos_) << percent_missing_ << std::endl;
    } else {
        const int colWidth = 12;
        // need to take into account longest_tax_label_
        longest_tax_label_ = get_longest_label(taxon_labels_);
        std::string pad = std::string(longest_tax_label_, ' ');
        // header
        (*poos_) << pad << " ";
        (*poos_) << std::right << std::setw(colWidth) << "Nchar" << " ";
        (*poos_) << std::right << std::setw(colWidth) << "Missing" << " ";
        (*poos_) << std::right << std::setw(colWidth) << "Proportion" << std::endl;
        
        // return nchar for individual seqs
        for (int i = 0; i < num_taxa_; i++) {
            int diff = longest_tax_label_ - taxon_labels_[i].size();
            (*poos_) << taxon_labels_[i];
            if (diff > 0) {
                pad = std::string(diff, ' ');
                (*poos_) << pad;
            }
            (*poos_) << " ";
            (*poos_) << std::right << std::setw(colWidth) << seq_lengths_[i] << " ";
            (*poos_) << std::right << std::setw(colWidth) << missing_counts_[i] << " ";
            (*poos_) << std::right << std::setw(colWidth)
                    << (double)missing_counts_[i] / (double)seq_lengths_[i] << std::endl;
        }
    }
    
}


// return whichever property set to true
void SeqInfo::get_property (const bool& get_labels, const bool& check_aligned,
        const bool& get_nseq, const bool& get_freqs, const bool& get_nchar,
        const double& get_missing) {
    
    if (get_labels) {
        collect_taxon_labels();
        for (unsigned int i = 0; i < taxon_labels_.size(); i++) {
            (*poos_) << taxon_labels_[i] << std::endl;
        }
    } else if (check_aligned) {
        check_is_aligned();
        (*poos_) << std::boolalpha << is_aligned_ << std::endl;
    } else if (get_nseq) {
        (*poos_) << num_taxa_ << std::endl;
    } else if (get_freqs) {
        calculate_freqs();
        return_freq_table();
    } else if (get_missing) {
        calc_missing();
        return_missing();
    } else if (get_nchar) {
        get_num_chars ();
        if (!output_indiv_) { // single return value
            if (seq_length_ != -1) {
                (*poos_) << seq_length_ << std::endl;
            } else {
                // not aligned
                (*poos_) << "sequences are not aligned" << std::endl;
            }
        } else { // individual lengths
            collect_taxon_labels();
            longest_tax_label_ = get_longest_label(taxon_labels_);
            for (int i = 0; i < num_taxa_; i++) {
                int diff = longest_tax_label_ - taxon_labels_[i].size();
                (*poos_) << taxon_labels_[i];
                if (diff > 0) {
                    std::string pad = std::string(diff, ' ');
                    (*poos_) << pad;
                }
                (*poos_) << " " << seq_lengths_[i] << std::endl;
            }
        }
    }
}


void SeqInfo::set_datatype () {
    if (alpha_name_ == "DNA") {
        is_dna_ = true;
        seq_type_ = "Nucl";
        set_alphabet();
    } else if (alpha_name_ == "AA") {
        is_protein_ = true;
        seq_type_ = "Prot";
        set_alphabet();
    } else if (alpha_name_ == "BINARY") {
        is_binary_ = true;
        seq_type_ = "Binary";
        set_alphabet();
    } else if (alpha_name_ == "MULTI") {
        is_multi_ = true;
        seq_type_ = "Multi";
        // alphabet is either 1) supplied (nexus) or 2) comes from entire alignment
    } else {
        std::cerr << "Error: cannot determine alignment type. Exiting." << std::endl;
        exit(0);
    }
    datatype_set_ = true;
}


// TODO: need to add morphology ('MULTI') data
// - may be passed in (nexus)
// - otherwise, needs to be gleaned from entire alignment
void SeqInfo::set_alphabet () {
    if (is_protein_) {
        seq_chars_ = "ACDEFGHIKLMNPQRSTVWYX";
    } else if (is_binary_) {
        seq_chars_ = "01";
    } else {
        seq_chars_ = "ACGT"; // not using ambiguity codes here
    }
    seq_chars_ += gap_;
    seq_chars_ += missing_;
    char_counts_.resize(seq_chars_.size(), 0);
    alpha_set_ = true;
}


// TODO: get rid of the concatenated_ business
void SeqInfo::summarize () {
    // a concatenated seq (i.e., across all indiv) will be used for all stats
    calculate_freqs();
    
    if (output_indiv_) {
        return_freq_table();
    } else {
        // pass the seq concatenated across all individuals
        // probably can skip a bunch of the stuff above...
        make_concatenated_sequence();
        count_chars_indiv_seq(concatenated_);
        print_summary_table_whole_alignment();
    }
}


void SeqInfo::set_alphabet_from_sampled_seqs (const std::string& seq) {
    seq_chars_ = get_alphabet_from_sequence(seq);
    // expecting order: valid, gap, missing
    // remove gap and missing (if present)
    seq_chars_.erase(std::remove(seq_chars_.begin(), seq_chars_.end(), gap_), seq_chars_.end());
    seq_chars_.erase(std::remove(seq_chars_.begin(), seq_chars_.end(), missing_), seq_chars_.end());
    
    // and append them back on
    seq_chars_ += gap_;
    seq_chars_ += missing_;
    char_counts_.resize(seq_chars_.size(), 0);
    alpha_set_ = true;
}
