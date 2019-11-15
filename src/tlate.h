#ifndef _TLATE_H_
#define _TLATE_H_

#include <string>
#include <map>

class TLATE {
private:
    std::map<std::string, std::string> table_;
    
public:
    TLATE(const std::string& table);
    std::string translate(std::string& dna);
    //virtual ~tlate();
};

#endif /* _TLATE_H_ */
