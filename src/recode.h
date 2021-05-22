#ifndef PX__RECODE_H
#define PX__RECODE_H

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
    
    static std::regex r_;
    static std::regex y_;
    static std::regex s_;
    static std::regex w_;
    static std::regex m_;
    static std::regex k_;
    static std::regex b_;
    static std::regex d_;
    static std::regex h_;
    static std::regex v_;
    
    SequenceRecoder (std::string& recodescheme);
    void parse_scheme ();
    void check_valid_scheme ();
    std::string get_recoded_seq (const std::string& origseq);
    void recode_seq (std::string& s);
    
    //~SequenceRecoder();
};

#endif /* PX__RECODE_H */
