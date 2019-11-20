#ifndef _UPGMA_H_
#define _UPGMA_H_

#include <string>
#include <vector>
#include <map>
#include <iostream>

class UPGMA {
private:
    std::map<std::string, std::string> sequences_;
    std::map<std::string, std::string>::iterator iter_;
    std::vector<std::string> names_;
    std::vector< std::vector<double> > distmatrix_;
    void set_name_key();
    
    // additions:
    int ntax_;
    int nchar_;
    std::string seqfile_;
    std::map<int, std::string> nameKey_;
    std::string newickstring_; // temporary

public:
    UPGMA (std::istream* pios);
    std::map<std::string, std::string> FastaToOneLine (std::ifstream&);
    std::vector< std::vector<double> > BuildMatrix (std::map<std::string,
        std::string>& sequences);
    void TREEMAKE (std::vector<std::string>&, std::map <int, std::string>&,
        std::vector< std::vector<double> >&);
    std::string get_newick ();
};

#endif /* _UPGMA_H_ */
