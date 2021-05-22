#ifndef PX__TLATE_H
#define PX__TLATE_H

#include <string>
#include <map>

class TLATE {
private:
    std::map<std::string, std::string> table_;
    
public:
    TLATE (const std::string& table);
    std::string translate (std::string& dna);
    //virtual ~tlate();
};

#endif /* PX__TLATE_H */
