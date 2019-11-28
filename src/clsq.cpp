#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <iomanip>

#include "clsq.h"
#include "sequence.h"
#include "seq_reader.h"
#include "utils.h"


SequenceCleaner::SequenceCleaner (std::istream* pios, double& prop_required, const bool& by_taxon,
        const bool& count_only, const bool& verbose):num_tax_(0), num_char_(0), num_retained_(0),
        missing_allowed_(1.0 - prop_required), by_taxon_(by_taxon), count_only_(count_only), verbose_(verbose) {
    read_in_alignment(pios);
    count_missing();
    if (!count_only_) {
        generate_cleaned_sequences();
    } else {
        
    }
}


void SequenceCleaner::read_in_alignment (std::istream* pios) {
    Sequence seq;
    std::string retstring;
    alpha_name_ = "";
    int ft = test_seq_filetype_stream(*pios, retstring);
    int file_ntax = 0; // ntax declared in the file itself
    std::string file_type = get_filetype_string(ft);
    
    if (file_type.compare("nexus") == 0) {
        bool is_interleaved = false;
        char gap, missing; // dummy variables, not used but required
        std::string symbols = ""; // another dummy
        get_nexus_alignment_properties (*pios, file_ntax, num_char_,
                is_interleaved, alpha_name_, symbols, gap, missing);
        //std::cout << "alpha_name_ = " << alpha_name_ << std::endl;
        retstring = ""; // have to set so seq_reader knows we are mid-file
        if (!is_interleaved) {
            while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
                seqs_.push_back(seq);
            }
        } else {
            seqs_ = read_interleaved_nexus(*pios, file_ntax, num_char_);
        }
    } else {
        bool complicated_phylip = false;
        // check if we are dealing with a complicated phylip format
        if (file_type.compare("phylip") == 0) {
            get_phylip_dimensions(retstring, file_ntax, num_char_);
            complicated_phylip = is_complicated_phylip(*pios, num_char_);
        }
        if (complicated_phylip) {
            seqs_ = read_phylip(*pios, file_ntax, num_char_);
            if (alpha_name_ == "") {
                alpha_name_ = seqs_[0].get_alpha_name();
            }
        } else {
            while (read_next_seq_from_stream(*pios, ft, retstring, seq)) {
                seqs_.push_back(seq);
                if (alpha_name_ == "") {
                    alpha_name_ = seq.get_alpha_name();
                }
            }
            if (ft == 2) { // fasta has an trailing one
                seqs_.push_back(seq);
            }
        }
    }
    
    // some simple error-checking
    if (file_ntax != 0) {
        if (file_ntax != (int)seqs_.size()) {
            std::cerr << "Error: number of taxa declared in the file ("
                << ") does not match the number read (" << seqs_.size()
                << "). Exiting." << std::endl;
            exit(1);
        }
    }
    num_tax_ = (int)seqs_.size();
    
    num_char_ = (int)seqs_[0].get_sequence().size();
    bool aligned = true;
    for (int i = 1; i < num_tax_; i++) {
        if ((int)seqs_[i].get_sequence().size() != num_char_) {
            aligned = false;
        }
    }
    if (!aligned) {
        std::cerr << "Sequences are not aligned. Exiting." << std::endl;
        exit(0);
    }
    set_bad_chars(); // uses alpha name
}


void SequenceCleaner::set_bad_chars () {
    if (alpha_name_ == "DNA") {
        badChars_ = "NX-?";
    } else if (alpha_name_ == "AA") {
        badChars_ = "X-?";
    } else if (alpha_name_ == "BINARY" || alpha_name_ == "MULTI") {
        badChars_ = "-?";
    } else {
        std::cout << "Don't know what kind of alignment this is :( Exiting." << std::endl;
        exit(0);
    }
}


std::vector<Sequence> SequenceCleaner::get_cleaned_seqs () {
    return cleaned_seqs_;
}


void SequenceCleaner::write_seqs (std::ostream* poos) {
    if (by_taxon_) {
        for (int i = 0; i < num_tax_; i++) {
            if (missing_per_taxon_proportion_[i] < missing_allowed_) {
                (*poos) << ">" << seqs_[i].get_id() << std::endl;
                (*poos) << seqs_[i].get_sequence() << std::endl;
            }
        }
    } else {
        for (int i = 0; i < num_tax_; i++) {
            std::string seq = cleaned_seqs_[i].get_sequence();
            if (seq.find_first_not_of(badChars_) != std::string::npos) {
                (*poos) << ">" << cleaned_seqs_[i].get_id() << std::endl;
                (*poos) << cleaned_seqs_[i].get_sequence() << std::endl;
            } else {
                if (verbose_) {
                    std::cerr << "Taxon '" << cleaned_seqs_[i].get_id()
                        << "' consists only of missing characters. Skipping." << std::endl;
                }
            }
            
        }
    }
}


void SequenceCleaner::write_stats (std::ostream* poos) {
    const char separator = ' ';
    const int colWidth = 10;
    if (by_taxon_) {
        std::string pad = "";
        int diff = 0;
        int longest = get_longest_taxon_label();
        // header
        (*poos) << "Taxon" << std::string((longest - 5), ' ');
        (*poos) << " " << std::right << std::setw(colWidth) << std::setfill(separator)
                << "Missing" << " " << std::right << std::setw(colWidth)
                << std::setfill(separator) << "Prop." << std::endl;
        (*poos) << std::string((longest + 22), '-') << std::endl;
        for (int i = 0; i < num_tax_; i++) {
            (*poos) << seqs_[i].get_id();
            diff = longest - seqs_[i].get_id().size();
            if (diff > 0) {
                pad = std::string(diff, ' ');
                (*poos) << pad;
            }
            (*poos) << " " << std::right << std::setw(colWidth) << std::setfill(separator)
                    << missing_per_taxon_[i] << " " << std::right << std::setw(colWidth)
                    << std::setfill(separator) << missing_per_taxon_proportion_[i] << std::endl;
        }
    } else {
        std::cout << "num_char_ = " << num_char_ << std::endl;
        (*poos) << std::right << std::setw(colWidth) << std::setfill(separator) << "Character" << " "
                << std::right << std::setw(colWidth) << std::setfill(separator) << "Missing" << " "
                << std::right << std::setw(colWidth) << std::setfill(separator) << "Prop."
                << std::endl;
        (*poos) << std::string(32, '-') << std::endl;
        for (int i = 0; i < num_char_; i++) {
            (*poos) << std::right << std::setw(colWidth) << std::setfill(separator) << i
                    << std::right << std::setw(colWidth) << std::setfill(separator) << missing_per_site_counts_[i]
                    << std::right << std::setw(colWidth) << std::setfill(separator) << missing_per_site_proportion_[i]
                    << std::endl;
        }
    }
}


// get the longest label. for printing purposes
int SequenceCleaner::get_longest_taxon_label () {
    int longest = 0;
    int curLength = 0;
    for (int i = 0; i < num_tax_; i++) {
        curLength = (int)seqs_[i].get_id().size();
        if (curLength > longest) {
            longest = curLength;
        }
    }
    return longest;
}


// hrm should we create a new set of seqs, or just edit existing one?
void SequenceCleaner::generate_cleaned_sequences () {
    std::string name = "";
    std::string seq_string = "";
    Sequence orig_seq;
    for (int i = 0; i < num_tax_; i++) {
        Sequence new_seq;
        orig_seq = seqs_[i];
        name = orig_seq.get_id();
        if (num_retained_ > 0) {
            seq_string = get_cleaned_seq(orig_seq.get_sequence());
        } else {
            seq_string = "-"; // for when all sites are removed
        }
        new_seq.set_id(name);
        new_seq.set_sequence(seq_string);
        cleaned_seqs_.push_back(new_seq);
    } 
}

// counts both per-site and per-taxon missing characters
void SequenceCleaner::count_missing () {
    // initialize empty vectors
    missing_per_site_counts_ = std::vector<int>(num_char_, 0);
    missing_per_site_proportion_ = std::vector<double>(num_char_, 0.0);
    missing_per_taxon_ = std::vector<int>(num_tax_, 0);
    missing_per_taxon_proportion_ = std::vector<double>(num_tax_, 0.0);
    
    std::string seq_string = "";
    
    for (int i = 0; i < num_tax_; i++) {
        seq_string = string_to_upper(seqs_[i].get_sequence());
        for (int j = 0; j < num_char_; j++) {
            if (badChars_.find(seq_string[j]) != std::string::npos) {
                missing_per_site_counts_[j]++;
                missing_per_taxon_[i]++;
                //std::cout << "  BAD CHAR (" << seq_string[j] << ")!" << std::endl;
            }
        }
        missing_per_taxon_proportion_[i] = (double)missing_per_taxon_[i] / (double)num_char_;
    }
    // get percentages
    for (int i = 0; i < num_char_; i++) {
        missing_per_site_proportion_[i] = (double)missing_per_site_counts_[i] / (double)num_tax_;
        //std::cout << i << ". missing = " << missing_per_site_counts_[i] << "("
        //        << missing_per_site_proportion_[i] << ")" << std::endl;
        if (missing_per_site_proportion_[i] <= missing_allowed_) {
            retained_sites_.push_back(i);
        }
    }
    num_retained_ = (int)retained_sites_.size();
}


std::string SequenceCleaner::get_cleaned_seq (const std::string& origseq) {
    std::string seq;
    for (unsigned int i = 0; i < retained_sites_.size(); i++) {
        if (i == 0) {
            seq = origseq[retained_sites_[i]];
        } else {
            seq += origseq[retained_sites_[i]];
        }
    }
    return seq;
}


SequenceCleaner::~SequenceCleaner() {
    // TODO Auto-generated destructor stub
}
