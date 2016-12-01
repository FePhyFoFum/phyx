/*
 * clsq.h
 *
 *  Created on: Jun 15, 2015
 *      Author: joe
 */

#ifndef _CLSQ_H_
#define _CLSQ_H_

#include <map>

using namespace std;

class SequenceCleaner {
private:
    
    int num_taxa_, num_char_;
    bool is_dna_;
    bool verbose_;
    double missing_allowed_;
    double required_present_;
    string fasta_, line_, dna_, name_hold_;
    map<string, string> sequences_;
    map<string, string>::iterator iter_;
    map<string, string> trimmed_seqs_;
    void read_sequences (istream* pios);
    void CheckMissing(double MissingData [], string& dna, bool& type);
    void clean_sequences();

public:
    SequenceCleaner(istream* pios, double& proportion, bool& force_protein,
        bool const& verbose);
    int get_num_taxa (); // not used
    map<string, string> get_trimmed_seqs (); // not used
    void write_seqs (ostream* poos);
    virtual ~SequenceCleaner();
};

#endif /* _CLSQ_H_ */
