#ifndef _SEQUENCE_H_
#define _SEQUENCE_H_

#include <vector>

using namespace std;

// not sure if we need this or not
// should maybe just store the DNA or AA?
// could even store those as separate
// JWB: yes, this is clunky
typedef enum {
    DNA = 0, AA = 1, BINARY = 2, MULTI = 3, CODON = 4, NA = 5
} seqAlpha; 

extern string dnachars;
extern string protchars;

class Sequence{
private:
    string id;
    string seq;
    unsigned int length;
    bool aligned;
    vector<double> qualarr;
    string qualstr;
    seqAlpha alphabet;
    vector<double> cont_chars;
    vector<int> multistate_chars;
  
public:
    Sequence();
    Sequence(string,string,bool);
    Sequence(string,string);
    seqAlpha get_alpha();
    string get_alpha_name();
    void set_alpha(seqAlpha);
    
    void infer_alpha ();
    bool is_dna_char (char & residue);
    bool is_prot_char(char & residue);
    
    bool is_aligned();
    string get_sequence()const;
    string get_id()const ;
    unsigned int get_length();
    void add_cont_char(double num);
    double get_cont_char(int _index);
    int get_num_cont_char();
    void clear_cont_char();
    void add_multistate_char(int num);
    int get_multistate_char(int _index);
    int get_num_multistate_char();
    void set_sequence(string seq);
    void set_id(string id);
    void set_aligned(bool al);
    void set_qualstr(string &,int);
    vector<double> get_qualarr();
    double get_qualarr_mean();
    string reverse_complement();
    void perm_reverse_complement();
    string get_fasta();
    string get_fasta(bool const& uppercase);
    string get_fastq();
    string seq_to_upper ();

};
#endif /* _SEQUENCE_H_ */
