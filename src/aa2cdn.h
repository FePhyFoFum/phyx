#ifndef _AA_TO_CDN_H_
#define _AA_TO_CDN_H_

#include <string>
#include <map>

#include "sequence.h"


class AAtoCDN {
private:
    // use vector<Sequence> instead of map
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
    AAtoCDN (const std::vector<Sequence>& nuc_seqs, const std::vector<Sequence>& aa_seqs,
        const bool& remove_last);
    void write_codon_alignment (std::ostream* poos);
    std::vector<Sequence> get_codon_alignment ();
};

#endif /* _AA_TO_CDN_H_ */
