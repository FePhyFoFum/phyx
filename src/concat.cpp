#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "seq_utils.h"
#include "concat.h"


// parent constructor. where things get concatenated
SequenceConcatenater::SequenceConcatenater (const bool& toupcase):num_partitions_(0),
        num_char_(0), num_taxa_(0), toupcase_(toupcase) {
}


// constructor for individual file
SequenceConcatenater::SequenceConcatenater (std::string& seqf, const bool& toupcase):num_partitions_(0),
        num_char_(0), num_taxa_(0), toupcase_(toupcase), filename_(seqf) {
    read_sequences();
}


void SequenceConcatenater::read_sequences () {
    std::istream * pios = new std::ifstream(filename_);
    
    std::string alphaName;
    seqs_ = ingest_alignment(pios, alphaName);
    num_taxa_ = static_cast<unsigned int>(seqs_.size());
    
    if (toupcase_) {
        for (unsigned int i = 0; i < num_taxa_; i++) {
            seqs_[i].set_sequence(string_to_upper(seqs_[i].get_sequence()));
        }
    }
    
    if (!is_aligned(seqs_)) {
        std::cerr << "Error: sequences in file '" << filename_ << "' are not aligned. Exiting."
            << std::endl;
        delete pios;
        exit(1);
    }
    num_char_ = seqs_[0].get_length();
    num_partitions_ = 1;
    partition_sizes_.push_back(num_char_);
    delete pios;
}


// where stuff actually happens
void SequenceConcatenater::concatenate(SequenceConcatenater& newSeqs) {
    std::string old_filler(num_char_, '-');
    unsigned int new_seq_len = newSeqs.get_sequence_length();
    std::string new_filler(static_cast<size_t>(new_seq_len), '-');
    num_char_ += new_seq_len;
    for (unsigned int i = 0; i != num_taxa_; i++) {
        bool match_found = false;
        if (newSeqs.num_taxa_ > 0) {
            for (unsigned int j = 0; j != newSeqs.num_taxa_; j++) {
                if (seqs_[i].get_id() == newSeqs.seqs_[j].get_id()) {
                    seqs_[i].set_sequence(seqs_[i].get_sequence() + newSeqs.seqs_[j].get_sequence());
                    match_found = true;
                    // erase matched entry so it won't have to be compared against again.
                    // erase in reverse order.
                    delete_sequence(newSeqs, j);
                    break;
                }
            }
        }
        if (!match_found) { // taxon is missing from present locus.
            seqs_[i].set_sequence(seqs_[i].get_sequence() + new_filler);
        }
    }

    // now, all that should be left are the novel sequences from the new file
    if (newSeqs.num_taxa_ > 0) {
        for (unsigned int i = 0; i != newSeqs.num_taxa_; i++) {
            newSeqs.seqs_[i].set_sequence(old_filler + newSeqs.seqs_[i].get_sequence());
            seqs_.push_back(newSeqs.seqs_[i]);
            num_taxa_++;
        }
    }
    num_partitions_++;
    partition_sizes_.push_back(newSeqs.get_sequence_length());
}


unsigned int SequenceConcatenater::get_sequence_length () const {
    return num_char_;
}


unsigned int SequenceConcatenater::get_num_taxa () const {
    return num_taxa_;
}


void SequenceConcatenater::delete_sequence (SequenceConcatenater& newSeqs,
        const unsigned int& index) {
    newSeqs.seqs_.erase(newSeqs.seqs_.begin() + index);
    newSeqs.num_taxa_--;
}


Sequence SequenceConcatenater::get_sequence (const unsigned int& index) const {
    return seqs_[index];
}


// not currently used
/*
std::vector<unsigned int> SequenceConcatenater::get_partition_sizes () const {
    return partition_sizes_;
}
*/

void SequenceConcatenater::write_partition_information (const std::vector<std::string>& inputFiles,
        std::string& partfile) {
    std::ofstream outfile(partfile.c_str());
    unsigned int charIndex = 1;
    
    // need to check seq type when writing this
    // use infer_alpha / get_alpha_name
    // but: are mixed seq types allowed? prolly...
    //     - so: need to check each one
    
    for (unsigned int i = 0; i < partition_sizes_.size(); i++) {
        unsigned int stopIndex = charIndex + partition_sizes_[i] - 1;
        bool going = true;
        std::string alpha;
        unsigned long j = 0;
        while (going) {
            Sequence terp = seqs_[j];
            std::string subseq = terp.get_sequence().substr((charIndex - 1), partition_sizes_[i]);
            // check if all are the same character (presumably all N, but useful either way)
            if (subseq.find_first_not_of(subseq.front()) != std::string::npos) {
                terp.set_sequence(subseq);
                alpha = terp.get_alpha_name();
                going = false;
            }
            j++;
        }
        outfile << alpha << ", " << inputFiles[i] << " = " << charIndex << "-" << stopIndex << std::endl;
        charIndex = stopIndex + 1;
    }
    outfile.close();
}
