#ifndef _LS_SEQ_H_
#define _LS_SEQ_H_

using namespace std;

class Stats {
private:
    string Concatenated;
    string temp_seq;
    string Molecule;
    string type;
    bool finished;
    string Mol;
    string name;
    map <char, double> Total;
    int seqcount;


public:
    Stats (istream* pios, bool& all, bool& prot, ostream* poos);
    void STAT_Getter(string& seq, bool& prot);
    void Printer(bool& prot, ostream* poos);
};

#endif /* _LS_SEQ_H_ */
