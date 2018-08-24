#ifndef _RECODE_H_
#define _RECODE_H_

#include <set>

using namespace std;

class SequenceRecoder {

public:
    string recodescheme_;
    
    static set <string> schemes_;
    
    SequenceRecoder (string & recodescheme);
    void check_valid_scheme ();
    string get_recoded_seq (string const& origseq);
    void ry_recode (string &s);
    
    //~SequenceRecoder();
};

#endif /* _RECODE_H_ */
