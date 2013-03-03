#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <string>
#include <vector>
using namespace std;

enum seqAlpha { DNA = 0, AA = 1, BINARY = 2, MULTI = 3 }; 

class Sequence{
private:
    string id;
    string seq;
    bool aligned;
    vector<double> qualarr;
    string qualstr;
    seqAlpha alphabet;

public:
    Sequence();
    Sequence(string,string,bool);
    Sequence(string,string);
    seqAlpha get_alpha();
    void set_alpha(seqAlpha);
    bool is_aligned();
    string get_sequence();
    string get_id();
    void set_sequence(string seq);
    void set_id(string id);
    void set_aligned(bool al);
    void set_qualstr(string &,int);
    vector<double> get_qualarr();
    double get_qualarr_mean();
    string reverse_complement();
    void perm_reverse_complement();
    string get_fasta();
    string get_fastq();
};
#endif
