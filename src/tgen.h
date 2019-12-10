#ifndef _TGEN_H_
#define _TGEN_H_

#include <string>
#include <vector>
#include <map>
#include <iostream>


class TopologyGenerator {
private:
    int num_taxa_;
    bool rooted_;
    std::string lprefix_;
    
    unsigned long int ntopos_; // number of possible topologies for n taxa. depends on rootedness
    int nedges_; // the number of edges in the final trees. depends on rootedness
    int curtax_; // the next taxon to add
    int curnode_; // the next internal node to add
    
    std::vector< std::vector< std::vector<int> > > trees_;
    std::vector<std::string> newicks_; // sweet, sweet results
    
    void initialize();
    int get_num_edges (const int& num_taxa, const int& rooted);
    std::vector< std::vector<int> > initialize_edge_matrix_unrooted (const int& n);
    std::vector< std::vector<int> > initialize_edge_matrix_rooted (const int& n);
    std::vector< std::vector< std::vector<int> > > add_taxon_unrooted (std::vector< std::vector< std::vector<int> > > edges,
        const int& taxon, const int& new_node);
    std::vector< std::vector< std::vector<int> > > add_taxon_rooted (std::vector< std::vector< std::vector<int> > > edges,
        const int& taxon, const int& new_node);
    void newick_from_tree_map (int node, std::map<int, std::vector<int> > m, std::string& tree);
    std::string edge_matrix_to_newick (const std::vector< std::vector<int> >& edges);
    void generate_trees ();
    
public:
    TopologyGenerator (const int& num_taxa, const bool& rooted, const std::string& lprefix);
    void get_newicks (std::ostream* poos);
    //~TopologyGenerator();
};

#endif /* _TGEN_H_ */
