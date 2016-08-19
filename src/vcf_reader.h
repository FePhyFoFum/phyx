#ifndef _VCF_READER_H_
#define _VCF_READER_H_

#include <vector>
#include <string>

using namespace std;

class VcfReader {
private:
    
    vector <string> taxa_;
    vector <string> seqs_;
    
    vector <string> get_alts (const string& str);
    void read_vcf (istream* pios);

public:
    VcfReader(istream* pios);
    void write_seqs (bool const& uppercase, ostream* poos);
};








#endif /* _VCF_READER_H_ */