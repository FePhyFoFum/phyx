#include <iostream>
#include <string>
#include <vector>

#include "aa2cdn.h"
#include "sequence.h"
#include "seq_utils.h"
#include "utils.h"


struct SequenceIDListCompare {
    bool operator()(const Sequence& lhs, const Sequence& rhs) {
      return lhs.get_id() < rhs.get_id();
    }
} SequenceIDListCompare;


AAtoCDN::AAtoCDN (std::vector<Sequence> nuc_seqs, std::vector<Sequence> aa_seqs,
        const bool& remove_last):remove_last_(remove_last), nuc_seqs_(std::move(nuc_seqs)),
        aa_seqs_(std::move(aa_seqs)) {
    // set up names
    nuc_names_ = collect_names(nuc_seqs_);
    aa_names_ = collect_names(aa_seqs_);
    
    sort(aa_seqs_.begin(), aa_seqs_.end(), SequenceIDListCompare);
    sort(nuc_seqs_.begin(), nuc_seqs_.end(), SequenceIDListCompare);
    check_names();
    generate_codon_alignment();
}


// check for taxa present in only nuc or aa alignments.
// print these to cerr, and remove from the respective alignment
void AAtoCDN::check_names () {    
    std::vector<std::string> diff;
    diff = in_first_not_second(nuc_names_, aa_names_);
    if (!diff.empty()) {
        std::cerr << "The following names are present in the nucleotide alignment, but not the protein alignment:"
                << std::endl;
        for (auto & d : diff) {
            std::cerr << d << std::endl;
            // remove from alignment
            for (unsigned int j = 0; j < nuc_seqs_.size(); j++) {
                if (nuc_seqs_[j].get_id() == d) {
                    nuc_seqs_.erase(nuc_seqs_.begin()+j);
                    break;
                }
            }
        }
    }
    if (nuc_seqs_.empty()) {
        std::cerr << "Error: no names are present in both the nucleotide and protein alignments. Exiting."
                << std::endl;
        exit(1);
    }
    
    // now the other direction
    diff = in_first_not_second(aa_names_, nuc_names_);
    if (!diff.empty()) {
        std::cerr << "The following names are present in the protein alignment, but not the nucleotide alignment:"
                << std::endl;
        for (const auto & di : diff) {
            std::cerr << di << std::endl;
            // remove from alignment
            for (unsigned int j = 0; j < aa_seqs_.size(); j++) {
                if (aa_seqs_[j].get_id() == di) {
                    aa_seqs_.erase(aa_seqs_.begin()+j);
                    break;
                }
            }
        }
    }
}


// at this point, alignments should be of the same size and order
void AAtoCDN::generate_codon_alignment () {
    std::string aaseq;
    std::string nucseq;
    std::string codonseq;
    Sequence seq;
    for (unsigned int i = 0; i < aa_seqs_.size(); i++) {
        seq.set_id(aa_seqs_[i].get_id());
        aaseq = aa_seqs_[i].get_sequence();
        nucseq = nuc_seqs_[i].get_sequence();
        codonseq = "";
        auto aalen = static_cast<unsigned int>(aaseq.size());
        
        // check that seq lengths correspond
        unsigned int ncodons = static_cast<unsigned int>(nucseq.length()) / 3u;
        unsigned int naachars = aalen - static_cast<unsigned int>(std::count(aaseq.begin(), aaseq.end(), '-'));
        if (ncodons != naachars) {
            std::cerr << "Error: for taxon '" << aa_seqs_[i].get_id()
                << "' nucleotide alignment involves " << ncodons
                << " codons, but protein alignment involves " << naachars
                << " amino acids. Skipping." << std::endl;
            continue;
        }
        
        if (remove_last_) {
            aalen--;
        }
        unsigned int nuccntr = 0;
        for (unsigned int j = 0; j < aalen; j++) {
            if (aaseq[j] == '-') {
                codonseq += "---";
            } else {
                for (int k = 0; k < 3; k++) {
                    codonseq += nucseq[nuccntr];
                    nuccntr++;
                }
            }
        }
        seq.set_sequence(codonseq);
        codon_seqs_.push_back(seq);
    }
}


void AAtoCDN::write_codon_alignment (std::ostream* poos) {
    for (auto & codon_seq : codon_seqs_) {
        (*poos) << ">" << codon_seq.get_id() << std::endl;
        (*poos) << codon_seq.get_sequence() << std::endl;
    }
}


// not currently used, but available
std::vector<Sequence> AAtoCDN::get_codon_alignment () const {
    return codon_seqs_;
}
