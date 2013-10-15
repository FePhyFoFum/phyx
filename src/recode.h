#ifndef SEQ_RECODE_H_
#define SEQ_RECODE_H_

class SequenceRecoder {

public:
    SequenceRecoder ();
    string get_recoded_seq (string const& origseq);
    void ry_recode (string &s);
    
    //~SequenceRecoder();
};

#endif /* SEQ_RECODE_H_ */
