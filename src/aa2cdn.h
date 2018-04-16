/*
 * aatocdn.h
 *
 *  Created on: Jun 15, 2015
 *      Author: joe
 */

#ifndef _AA_TO_CDN_H_
#define _AA_TO_CDN_H_

#include <map>

using namespace std;

class AAtoCDN {
private:
    map <string, string> codon_sequences_;
    map <string, string>::iterator iter_;
    string amino_acid_sequence_;
    string nucleotide_sequence_;

public:
    AAtoCDN();
    map <string, string> convert_to_codons (map <string, string>& aa_sequences,
    map <string, string>& nuc_sequences, bool& rm_last);

};

#endif /* _AATOCDN_H_ */
