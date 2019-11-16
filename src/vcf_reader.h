#ifndef _VCF_READER_H_
#define _VCF_READER_H_

#include <vector>
#include <string>
#include <iostream>

class VcfReader {
private:
    std::vector<std::string> taxa_;
    std::vector<std::string> seqs_;
    
    std::vector<std::string> get_alts (const std::string& str);
    void read_vcf (std::istream* pios);

public:
    VcfReader (std::istream* pios);
    void write_seqs (const bool& uppercase, std::ostream* poos);
};

#endif /* _VCF_READER_H_ */
