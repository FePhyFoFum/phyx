#ifndef _RECODE_H_
#define _RECODE_H_

#include <string>
#include <set>
#include <regex>

class SequenceRecoder {

public:
    std::string recodescheme_;
    
    // individual coding scheme
    bool R_;
    bool Y_;
    bool S_;
    bool W_;
    bool M_;
    bool K_;
    bool B_;
    bool D_;
    bool H_;
    bool V_;
    bool A_;
    bool C_;
    bool G_;
    bool T_;
    
    static std::set<char> recognized_;
    
    static regex r_;
    static regex y_;
    static regex s_;
    static regex w_;
    static regex m_;
    static regex k_;
    static regex b_;
    static regex d_;
    static regex h_;
    static regex v_;
    
    SequenceRecoder(string& recodescheme);
    void parse_scheme();
    void check_valid_scheme();
    std::string get_recoded_seq(const std::string& origseq);
    void recode_seq(std::string& s);
    
    //~SequenceRecoder();
};

#endif /* _RECODE_H_ */
