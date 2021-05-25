#ifndef PX_UPGMA_H
#define PX_UPGMA_H

#include <string>
#include <vector>
#include <iostream>

#include "sequence.h"

class Tree; // forward declaration


class UPGMA {
private:
    int num_taxa_; 
    int num_char_;
    std::string newickstring_;
    std::vector<std::string> names_;
    std::vector<Sequence> seqs_;
    std::vector< std::vector<double> > full_distmatrix_;
    Tree* tree_;
    
    void update_tree (std::string& newname, std::vector<std::string>& names,
        std::vector< std::vector<double> >& newmatrix);
    std::vector< std::vector<double> > build_matrix ();
    double get_smallest_distance (const std::vector< std::vector<double> >& dmatrix,
        unsigned long& mini1, unsigned long& mini2);
    void construct_tree ();

public:
    UPGMA (std::istream* pios);
    std::vector< std::vector<double> > get_matrix () const;
    std::string get_newick ();
};

#endif /* PX_UPGMA_H */
