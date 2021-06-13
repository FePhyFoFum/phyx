#ifndef PX_SEQUENCE_H
#define PX_SEQUENCE_H

#include <vector>
#include <string>


typedef enum {
    DNA = 0, AA = 1, BINARY = 2, MULTI = 3, CODON = 4, NA = 5
} seqAlpha; 

extern std::string dnachars_with_ambiguous;
extern std::string protchars;

class Sequence {
    // these member variables should have lagging underscores
private:
    std::string id_;
    std::string seq_;
    unsigned int length_;
    bool aligned_;
    std::vector<double> qualarr_;
    std::string qualstr_;
    seqAlpha alphabet_;
    std::vector<double> cont_chars_;
    std::vector<int> multistate_chars_;
  
public:
    Sequence ();
    Sequence (std::string, std::string, bool);
    Sequence (std::string, std::string);
    seqAlpha get_alpha () const;
    std::string get_alpha_name ();
    void set_alpha (seqAlpha);
    
    void infer_alpha ();
    
    bool is_aligned () const;
    std::string get_sequence ()const;
    std::string get_id ()const;
    unsigned int get_length ();
    void add_cont_char (double _num);
    double get_cont_char (int _index);
    int get_num_cont_char ();
    void clear_cont_char ();
    void add_multistate_char (int _num);
    int get_multistate_char (int _index);
    int get_num_multistate_char ();
    void set_sequence (std::string _seq);
    void set_id (std::string _id);
    void set_aligned (bool _align);
    void set_qualstr (std::string&, int);
    std::vector<double> get_qualarr () const;
    double get_qualarr_mean ();
    std::string reverse_complement ();
    void perm_reverse_complement ();
    std::string get_fasta ();
    std::string get_fasta (const bool& uppercase);
    std::string get_fastq ();
    std::string seq_to_upper ();

};
#endif /* PX_SEQUENCE_H */
