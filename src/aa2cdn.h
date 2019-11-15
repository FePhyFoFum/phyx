#ifndef _AA_TO_CDN_H_
#define _AA_TO_CDN_H_

#include <string>
#include <map>

class AAtoCDN {
private:
    std::map <std::string, std::string> codon_sequences_;
    std::map <std::string, std::string>::iterator iter_;
    std::string amino_acid_sequence_;
    std::string nucleotide_sequence_;

public:
    AAtoCDN ();
    std::map <std::string, std::string> convert_to_codons (std::map <std::string,
        std::string>& aa_sequences, std::map <std::string,
        std::string>& nuc_sequences, bool& rm_last);
};

#endif /* _AA_TO_CDN_H_ */
