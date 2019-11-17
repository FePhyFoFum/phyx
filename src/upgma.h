#ifndef _UPGMA_H_
#define _UPGMA_H_

#include <string>
#include <vector>
#include <map>
#include <iostream>

class UPGMA {
private:
    std::map<std::string, std::string> sequences;
    std::map<std::string, std::string>::iterator iter;
    std::vector<std::string> names;
    std::vector< std::vector<double> > Matrix;
    void set_name_key();
    
    // additions:
    int ntax;
    int nchar;
    std::string seqfile;
    std::map<int, std::string> NameKey;
    std::string newickstring; // temporary

public:
    UPGMA ();
    UPGMA (std::istream* pios);
    std::map<std::string, std::string> FastaToOneLine (std::ifstream&);
    std::vector< std::vector<double> > BuildMatrix (std::map<std::string,
        std::string>& sequences);
    void TREEMAKE (std::vector<std::string>&, std::map <int, std::string>&,
        std::vector< std::vector<double> >&);
    std::string get_newick ();
};

#endif /* _UPGMA_H_ */
