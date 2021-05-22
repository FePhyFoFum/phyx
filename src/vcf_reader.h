#ifndef PX__VCF_READER_H
#define PX__VCF_READER_H

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

#endif /* PX__VCF_READER_H */
