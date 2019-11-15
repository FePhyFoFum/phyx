/*
 * clsq.h
 *
 *  Created on: Jun 15, 2015
 *      Author: joe
 */

#ifndef _CLSQ_H_
#define _CLSQ_H_

#include <string>
#include <map>
#include <iostream>

class SequenceCleaner {
private:
    int num_taxa_, num_char_;
    bool is_dna_;
    bool verbose_;
    double missing_allowed_;
    double required_present_;
    std::string fasta_, line_, dna_, name_hold_;
    std::map<std::string, std::string> sequences_;
    std::map<std::string, std::string>::iterator iter_;
    std::map<std::string, std::string> trimmed_seqs_;
    void read_sequences (std::istream* pios);
    void CheckMissing(double MissingData [], std::string& dna, bool& type);
    void clean_sequences();

public:
    SequenceCleaner(std::istream* pios, double& proportion, bool& force_protein,
        const bool& verbose);
    int get_num_taxa(); // not used
    std::map<std::string, std::string> get_trimmed_seqs(); // not used
    void write_seqs(std::ostream* poos);
    virtual ~SequenceCleaner();
};

#endif /* _CLSQ_H_ */
