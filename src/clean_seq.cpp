#include <iostream>
#include <string>
#include <vector>
#include <iomanip>

#include "clean_seq.h"
#include "sequence.h"
#include "seq_reader.h"
#include "seq_utils.h"
#include "utils.h"


SequenceCleaner::SequenceCleaner (std::istream* pios, double& prop_required,
        const bool& remove_empty, const int& min_chars, const bool& by_taxon,
        const bool& by_codon, const bool& count_only, const bool& verbose):num_taxa_(0u),
        num_char_(0u), num_retained_(0u), min_chars_per_site_(min_chars),
        missing_allowed_(1.0 - prop_required), by_taxon_(by_taxon), by_codon_(by_codon),
        count_only_(count_only), verbose_(verbose), remove_empty_(remove_empty) {
    read_in_sequences(pios);
    if (min_chars_per_site_ != 0) {
        min_chars_ = true;
        if (min_chars_per_site_ > num_taxa_) {
            std::cerr << "Error: minimum characters required (" << min_chars_per_site_
                    << ") exceeds number of taxa (" << num_taxa_ << "). Exiting."
                    << std::endl;
            exit(0);
        }
        max_missing_ = num_taxa_ - min_chars_per_site_;
    }
    count_missing();
    if (!count_only_) {
        generate_cleaned_sequences();
    }
}


void SequenceCleaner::read_in_sequences (std::istream* pios) {
    seqs_ = ingest_alignment(pios, alpha_name_);
    num_taxa_ = static_cast<unsigned int>(seqs_.size());
    
    // check aligned. if not, BAIL
    bool aligned = is_aligned(seqs_);
    if (!aligned) {
        std::cerr << "Error: sequences are not aligned. Exiting." << std::endl;
        exit(0);
    }
    
    // if codons, length must be a multiple of 3
    if (by_codon_) {
        if (!is_codon_alignment(seqs_)) {
            std::cerr << "Error: sequences do not appear to be codons (i.e., length a multiple of 3). Exiting."
                << std::endl;
            exit(1);
        }
        if (alpha_name_ != "DNA") {
            std::cerr << "Error: codon alignments require DNA data type, but '"
                << alpha_name_ << "' detected. Exiting." << std::endl;
            exit(1);
        }
    }
    
    num_char_ = seqs_[0].get_length();
    set_bad_chars(); // uses alpha name
}


void SequenceCleaner::set_bad_chars () {
    if (alpha_name_ == "DNA") {
        badChars_ = "N-?";
    } else if (alpha_name_ == "AA") {
        badChars_ = "X-?";
    } else if (alpha_name_ == "BINARY" || alpha_name_ == "MULTI") {
        badChars_ = "-?";
    } else {
        std::cerr << "Error: cannot determine alignment type. Exiting." << std::endl;
        exit(0);
    }
}


std::vector<Sequence> SequenceCleaner::get_cleaned_seqs () const {
    return cleaned_seqs_;
}


void SequenceCleaner::write_seqs (std::ostream* poos) {
    if (by_taxon_) {
        for (unsigned int i = 0; i < num_taxa_; i++) {
            if (missing_per_taxon_proportion_[i] < missing_allowed_) {
                (*poos) << ">" << seqs_[i].get_id() << std::endl;
                (*poos) << seqs_[i].get_sequence() << std::endl;
            }
        }
    } else {
        for (unsigned int i = 0; i < num_taxa_; i++) {
            std::string seq = cleaned_seqs_[i].get_sequence();
            if (seq.find_first_not_of(badChars_) != std::string::npos) {
                (*poos) << ">" << cleaned_seqs_[i].get_id() << std::endl;
                (*poos) << cleaned_seqs_[i].get_sequence() << std::endl;
            } else {
                // rare case where removal of sites leaves only missing data for a taxon
                if (verbose_) {
                    std::cerr << "Taxon '" << cleaned_seqs_[i].get_id()
                        << "' consists only of missing characters. Removing." << std::endl;
                }
            }
        }
    }
}



// TODO: divide by 3 if codons
void SequenceCleaner::write_stats (std::ostream* poos) {
    const int colWidth = 12;
    if (by_taxon_) {
        if (!by_codon_) {
            (*poos) << "Length of sequences: " << num_char_ << " characters" << std::endl;
        } else {
            (*poos) << "Length of sequences: " << num_char_/3 << " codons" << std::endl;
        }
        
        std::string pad;
        unsigned int longest = get_longest_taxon_label();
        // header
        (*poos) << "Taxon" << std::string((longest - 5), ' ');
        (*poos) << " " << std::right << std::setw(colWidth)
                << "Missing" << " " << std::right << std::setw(colWidth)
                << "Prop." << std::endl;
        (*poos) << std::string((longest + 26), '-') << std::endl;
        for (unsigned int i = 0; i < num_taxa_; i++) {
            (*poos) << seqs_[i].get_id();
            auto diff = longest - seqs_[i].get_id().size();
            if (diff > 0) {
                pad = std::string(diff, ' ');
                (*poos) << pad;
            }
            if (!by_codon_) {
                (*poos) << " " << std::right << std::setw(colWidth)
                    << missing_per_taxon_[i] << " " << std::right << std::setw(colWidth)
                    << missing_per_taxon_proportion_[i] << std::endl;
            } else {
                (*poos) << " " << std::right << std::setw(colWidth)
                    << missing_per_taxon_[i]/3 << " " << std::right << std::setw(colWidth)
                    << missing_per_taxon_proportion_[i] << std::endl;
            }
            
        }
    } else {
        (*poos) << "Number of sequences: " << num_taxa_ << std::endl;
        if (!by_codon_) {
            (*poos) << std::right << std::setw(9) << "Character" << " "
                    << std::right << std::setw(colWidth) << "Missing" << " "
                    << std::right << std::setw(colWidth) << "Prop." << std::endl;
            (*poos) << std::string(35, '-') << std::endl;
            for (unsigned int i = 0; i < num_char_; i++) {
                (*poos) << std::right << std::setw(9) << i << " "
                        << std::right << std::setw(colWidth) << missing_per_site_counts_[i] << " "
                        << std::right << std::setw(colWidth) << missing_per_site_proportion_[i] << " "
                        << std::endl;
            }
        } else {
            (*poos) << std::right << std::setw(5) << "Codon" << " "
                    << std::right << std::setw(colWidth) << "Missing" << " "
                    << std::right << std::setw(colWidth) << "Prop." << std::endl;
            (*poos) << std::string(34, '-') << std::endl;
            for (unsigned int i = 0; i < (num_char_/3); i++) {
                unsigned int pos = i * 3;
                (*poos) << std::right << std::setw(5) << i << " "
                        << std::right << std::setw(colWidth) << missing_per_site_counts_[pos] << " "
                        << std::right << std::setw(colWidth) << missing_per_site_proportion_[pos] << " "
                        << std::endl;
            }
        }
    }
}


// get the longest label. for printing purposes
unsigned int SequenceCleaner::get_longest_taxon_label () {
    unsigned int longest = 0;
    for (unsigned int i = 0; i < num_taxa_; i++) {
        unsigned int curLength = static_cast<unsigned int>(seqs_[i].get_id().size());
        if (curLength > longest) {
            longest = curLength;
        }
    }
    return longest;
}


// hrm should we create a new set of seqs, or just edit existing one?
void SequenceCleaner::generate_cleaned_sequences () {
    std::string name;
    std::string seq_string;
    Sequence orig_seq;
    for (unsigned int i = 0; i < num_taxa_; i++) {
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
// TODO: need to consider codon seqs
void SequenceCleaner::count_missing () {
    // initialize empty vectors
    missing_per_site_counts_ = std::vector<int>(num_char_, 0);
    missing_per_site_proportion_ = std::vector<double>(num_char_, 0.0);
    missing_per_taxon_ = std::vector<int>(num_taxa_, 0);
    missing_per_taxon_proportion_ = std::vector<double>(num_taxa_, 0.0);
    
    std::string seq_string;
    
    if (!by_codon_) {
        for (unsigned int i = 0; i < num_taxa_; i++) {
            seq_string = string_to_upper(seqs_[i].get_sequence());
            for (unsigned int j = 0; j < num_char_; j++) {
                if (badChars_.find(seq_string[j]) != std::string::npos) {
                    missing_per_site_counts_[j]++;
                    missing_per_taxon_[i]++;
                    //std::cout << "  BAD CHAR (" << seq_string[j] << ")!" << std::endl;
                }
            }
            missing_per_taxon_proportion_[i] = static_cast<double>(missing_per_taxon_[i])
                    / static_cast<double>(num_char_);
        }
    } else {
        std::string codon;
        for (unsigned int i = 0; i < num_taxa_; i++) {
            seq_string = string_to_upper(seqs_[i].get_sequence());
            for (unsigned int j = 0; j < num_char_; j += 3) {
                codon = seq_string.substr(j, 3);
                if (codon.find_first_of(badChars_) != std::string::npos) {
                    // if any char is bad, they all are
                    // counts are still stored per site, not per codon
                    missing_per_site_counts_[j]++;
                    missing_per_site_counts_[j+1]++;
                    missing_per_site_counts_[j+2]++;
                    missing_per_taxon_[i] += 3;
                    // std::cout << "  BAD CODON (" << codon << ")!" << std::endl;
                }
            }
            missing_per_taxon_proportion_[i] = static_cast<double>(missing_per_taxon_[i])
                    / static_cast<double>(num_char_);
        }
    }
    
    // get proportions
    for (unsigned int i = 0; i < num_char_; i++) {
        if (min_chars_) {
            if (static_cast<unsigned int>(missing_per_site_counts_[i]) <= max_missing_) {
                retained_sites_.push_back(i);
            }
        } else if (remove_empty_) { // equivalent to min_chars_per_site_ == 1
            if (missing_per_site_counts_[i] != static_cast<int>(num_taxa_)) {
                retained_sites_.push_back(i);
            }
        } else {
            missing_per_site_proportion_[i] = static_cast<double>(missing_per_site_counts_[i])
                / static_cast<double>(num_taxa_);
            //std::cout << i << ". missing = " << missing_per_site_counts_[i] << "("
            //        << missing_per_site_proportion_[i] << ")" << std::endl;
            if (missing_per_site_proportion_[i] <= missing_allowed_) {
                retained_sites_.push_back(i);
            }
        }
    }
    
    num_retained_ = static_cast<unsigned int>(retained_sites_.size());
}


std::string SequenceCleaner::get_cleaned_seq (const std::string& origseq) const {
    std::string seq;
    for (unsigned int i = 0; i < num_retained_; i++) {
        seq += origseq[retained_sites_[i]];
    }
    return seq;
}


SequenceCleaner::~SequenceCleaner() = default;
