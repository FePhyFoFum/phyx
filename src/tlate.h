#ifndef PX_TLATE_H
#define PX_TLATE_H

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

#endif /* PX_TLATE_H */
