#ifndef _LS_SEQ_H_
#define _LS_SEQ_H_

using namespace std;

class Stats {
private:
    string concatenated_;
    string temp_seq_;
    string seq_chars_;
    string type_;
    bool finished_;
    string seq_type_;
    string name_;
    map <char, double> total_;
    int seqcount_;

public:
    Stats (istream* pios, bool& all, bool& prot, ostream* poos);
    void STAT_Getter(string& seq, bool& prot);
    void Printer(bool& prot, ostream* poos);
};

#endif /* _LS_SEQ_H_ */
