#ifndef PX_CLSQ_H
#define PX_CLSQ_H

#include <string>
#include <vector>
#include <iostream>

#include "sequence.h"

class SequenceCleaner {
private:
    int num_taxa_;
    int num_char_;
    int num_retained_;
    double missing_allowed_;
    
    bool by_taxon_;
    bool by_codon_;
    bool count_only_;
    bool verbose_;
    
    std::string badChars_;
    std::string alpha_name_;
    
    // refactored version
    std::vector<Sequence> seqs_;
    std::vector<Sequence> cleaned_seqs_;
    std::vector<int> missing_per_site_counts_;
    std::vector<double> missing_per_site_proportion_;
    std::vector<int> missing_per_taxon_;
    std::vector<double> missing_per_taxon_proportion_;
    std::vector<int> retained_sites_;
    
    void count_missing ();
    void generate_cleaned_sequences ();
    std::string get_cleaned_seq (const std::string& origseq);
    int get_longest_taxon_label ();
    void read_in_sequences (std::istream* pios);
    void set_bad_chars ();

public:
    SequenceCleaner (std::istream* pios, double& prop_required,
        const bool& by_taxon, const bool& by_codon, const bool& count_only,
        const bool& verbose);
    std::vector<Sequence> get_cleaned_seqs (); // not used, but available
    void write_seqs (std::ostream* poos);
    void write_stats (std::ostream* poos);
    virtual ~SequenceCleaner ();
};

#endif /* PX_CLSQ_H */
