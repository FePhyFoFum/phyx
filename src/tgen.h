#ifndef PX_TGEN_H
#define PX_TGEN_H

#include <string>
#include <vector>
#include <map>
#include <iostream>


class TopologyGenerator {
private:
    unsigned int num_taxa_;
    bool rooted_;
    std::string lprefix_;
    
    unsigned long int ntopos_; // number of possible topologies for n taxa. depends on rootedness
    unsigned int nedges_; // the number of edges in the final trees. depends on rootedness
    unsigned int curtax_; // the next taxon to add
    unsigned int curnode_; // the next unsigned internal node to add
    
    std::vector< std::vector< std::vector<unsigned int> > > trees_;
    std::vector<std::string> newicks_; // sweet, sweet results
    
    void initialize();
    unsigned int get_num_edges (const unsigned int& num_taxa, const unsigned int& rooted);
    std::vector< std::vector<unsigned int> > initialize_edge_matrix_unrooted (const unsigned int& n);
    std::vector< std::vector<unsigned int> > initialize_edge_matrix_rooted (const unsigned int& n);
    std::vector< std::vector< std::vector<unsigned int> > > add_taxon_unrooted (std::vector< std::vector< std::vector<unsigned int> > > edges,
        const unsigned int& taxon, const unsigned int& new_node);
    std::vector< std::vector< std::vector<unsigned int> > > add_taxon_rooted (std::vector< std::vector< std::vector<unsigned int> > > edges,
        const unsigned int& taxon, const unsigned int& new_node);
    void newick_from_tree_map (unsigned int node, std::map<unsigned int, std::vector<unsigned int> > m, std::string& tree);
    std::string edge_matrix_to_newick (const std::vector< std::vector<unsigned int> >& edges);
    void generate_trees ();
    
public:
    TopologyGenerator (const unsigned int& num_taxa, const bool& rooted, const std::string& lprefix);
    void get_newicks (std::ostream* poos);
    //~TopologyGenerator();
};

#endif /* PX_TGEN_H */
