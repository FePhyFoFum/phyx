#ifndef _SEQUENCE_H_
#define _SEQUENCE_H_

#include <vector>
#include <string>

// not sure if we need this or not
// should maybe just store the DNA or AA?
// could even store those as separate
// JWB: yes, this is clunky
typedef enum {
    DNA = 0, AA = 1, BINARY = 2, MULTI = 3, CODON = 4, NA = 5
} seqAlpha; 

extern std::string dnachars;
extern std::string protchars;

class Sequence {
private:
    std::string id;
    std::string seq;
    unsigned int length;
    bool aligned;
    std::vector<double> qualarr;
    std::string qualstr;
    seqAlpha alphabet;
    std::vector<double> cont_chars;
    std::vector<int> multistate_chars;
  
public:
    Sequence();
    Sequence(std::string, std::string, bool);
    Sequence(std::string, std::string);
    seqAlpha get_alpha();
    std::string get_alpha_name();
    void set_alpha(seqAlpha);
    
    void infer_alpha();
    bool is_dna_char(char& residue);
    bool is_prot_char(char& residue);
    
    bool is_aligned();
    std::string get_sequence()const;
    std::string get_id()const;
    unsigned int get_length();
    void add_cont_char(double num);
    double get_cont_char(int _index);
    int get_num_cont_char();
    void clear_cont_char();
    void add_multistate_char(int num);
    int get_multistate_char(int _index);
    int get_num_multistate_char();
    void set_sequence(std::string seq);
    void set_id(std::string id);
    void set_aligned(bool al);
    void set_qualstr(std::string&, int);
    std::vector<double> get_qualarr();
    double get_qualarr_mean();
    std::string reverse_complement();
    void perm_reverse_complement();
    std::string get_fasta();
    std::string get_fasta(const bool& uppercase);
    std::string get_fastq();
    std::string seq_to_upper ();

};
#endif /* _SEQUENCE_H_ */
