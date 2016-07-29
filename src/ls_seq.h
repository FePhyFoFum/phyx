#ifndef _LS_SEQ_H_
#define _LS_SEQ_H_

#include <map>
#include <vector>

using namespace std;

class SeqInfo {
private:
    string concatenated_;
    string temp_seq_;
    string seq_chars_; // the alphabet (DNA or AA)
    string file_type_;
    bool finished_;
    bool is_protein_;
    string seq_type_;
    string name_;
    map <char, double> total_;
    int seqcount_;
    bool output_indiv_; // report stats for each seq
    
    // new stuff
    vector <string> taxon_labels_;
    bool aligned_;
    int seq_length_;
    istream* pios_;
    ostream* poos_;
    
    void collect_taxon_labels ();
    void check_is_aligned ();
    void get_nseqs ();
    void get_nchars ();
    void set_alphabet ();
    void count_chars_indiv_seq(string& seq);
    void print_stats(ostream* poos);

public:
    SeqInfo (istream* pios, bool& indiv, bool const& force_protein, ostream* poos);
    SeqInfo (istream* pios, ostream* poos, bool& indiv, bool const& force_protein);
    void summarize ();
    void get_property (bool const& get_labels, bool const& check_aligned,
        bool const& get_nseq, bool const& get_freqs, bool const& get_nchar);
};

#endif /* _LS_SEQ_H_ */
