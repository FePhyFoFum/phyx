#ifndef _LS_SEQ_H_
#define _LS_SEQ_H_

#include <string>
#include <map>
#include <vector>
#include <iostream>

#include "sequence.h"

class SeqInfo {
private:
    std::string concatenated_;
    std::string temp_seq_;
    std::string seq_chars_; // the alphabet
    std::string file_type_; //"nexus", "phylip", "fasta", "fastq"
    bool output_indiv_; // report stats for each seq

    bool datatype_set_;
    bool is_dna_;
    bool is_protein_;
    bool is_multi_;
    bool is_binary_;
    bool alpha_set_;
    
    std::string alpha_name_; // phyx seq ids: DNA, AA, BINARY, MULTI
    std::string seq_type_; // label used for output table: Prot, Nucl, Mult, Binary
    char gap_;
    char missing_;
    std::map<char, double> total_;
    int num_taxa_;
    double percent_missing_;
    
    // new stuff
    std::vector<Sequence> seqs_;
    std::vector<std::string> taxon_labels_;
    std::vector<int> seq_lengths_;
    std::vector<int> char_counts_; // length seq_chars_ (i.e. the alphabet). accumulated across all seqs
    std::vector<int> missing_counts_; // for individual seqs
    std::vector< std::vector<int> > indiv_char_counts_;
    bool is_aligned_;
    int seq_length_;
    std::istream* pios_;
    std::ostream* poos_;
    int longest_tax_label_;
    
    void read_in_alignment ();
    void collect_taxon_labels ();
    void check_is_aligned ();
    void make_concatenated_sequence ();
    void get_num_chars ();
    void set_alphabet ();
    void count_chars_indiv_seq (std::string& seq);
    void count_chars (std::string& seq);
    void print_summary_table_whole_alignment ();
    void return_freq_table ();
    void calculate_freqs ();
    void calc_missing ();
    void set_datatype ();
    void set_alphabet_from_sampled_seqs (const std::string& seq);
    void return_missing ();
    
public:
    SeqInfo (std::istream* pios, std::ostream* poos, bool& indiv, const bool& force_protein);
    void summarize ();
    void get_property (const bool& get_labels, const bool& check_aligned,
        const bool& get_nseq, const bool& get_freqs, const bool& get_nchar,
        const double& get_missing);
};

#endif /* _LS_SEQ_H_ */
