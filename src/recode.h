#ifndef _RECODE_H_
#define _RECODE_H_

class SequenceRecoder {

public:
    SequenceRecoder ();
    string get_recoded_seq (string const& origseq);
    void ry_recode (string &s);
    
    //~SequenceRecoder();
};

#endif /* _RECODE_H_ */
