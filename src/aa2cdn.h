#ifndef PX__AA_TO_CDN_H
#define PX__AA_TO_CDN_H

#include <string>

#include "sequence.h"


class AAtoCDN {
private:
    bool remove_last_;
    std::vector<Sequence> nuc_seqs_;
    std::vector<Sequence> aa_seqs_;
    std::vector<Sequence> codon_seqs_;
    std::vector<std::string> nuc_names_;
    std::vector<std::string> aa_names_;
    
    void check_names ();
    void generate_codon_alignment ();

public:
    AAtoCDN ();
    AAtoCDN (std::vector<Sequence> nuc_seqs, std::vector<Sequence> aa_seqs,
        const bool& remove_last);
    void write_codon_alignment (std::ostream* poos);
    std::vector<Sequence> get_codon_alignment ();
};

#endif /* PX__AA_TO_CDN_H */
