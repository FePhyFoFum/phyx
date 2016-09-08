#ifndef _COMP_TEST_H_
#define _COMP_TEST_H_

#include <map>
#include <vector>

using namespace std;

class CompTest {
private:
    string seq_chars_; // the alphabet (DNA or AA only at present)
    bool is_protein_;
    int seqcount_;
    int df_;
    double test_stat_;
    int total_; // all valid chars
    
    vector <string> taxon_labels_;
    vector < vector <int> > indiv_char_counts_;
    vector <int> row_totals_;
    vector <int> col_totals_;
    istream* pios_;
    ostream* poos_;
    int longest_tax_label_;
    
    void set_alphabet ();
    void count_chars(string& seq);
    void return_freq_table (ostream* poos);
    void get_longest_taxon_label ();
    void read_seqs ();
    double get_cell_value (double const& observed, double const& expected);
    double calc_chi_square ();
    double calc_chi_square_prob ();
    
public:
    CompTest (istream* pios, ostream* poos, bool const& force_protein);

};



#endif /* _COMP_TEST_H_ */