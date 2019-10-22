#ifndef _TGEN_H_
#define _TGEN_H_

#include <map>

using namespace std;

class TopologyGenerator {
private:
    int ntax_;
    bool rooted_;
    string lprefix_;
    
    int ntopos_; // number of possible topologies for n taxa. depends on rootedness
    int nedges_; // the number of edges in the final trees. depends on rootedness
    int curtax_; // the next taxon to add
    int curnode_; // the next internal node to add
    
    vector < vector < vector <int> > > trees_;
    vector <string> newicks_; // sweet, sweet results
    
    void initialize ();
    int get_num_edges (int const& ntax, int const& rooted);
    vector < vector <int> > initialize_edge_matrix_unrooted (int const& n);
    vector < vector < vector <int> > > add_taxon_unrooted (vector < vector < vector <int> > > edges,
        int const& taxon, int const& new_node);
    void newick_from_tree_map (int node, map <int, vector <int>> m, string & tree);
    string edge_matrix_to_newick (vector < vector <int> > const& edges);
    void generate_trees ();
    
public:
    TopologyGenerator(int const& ntax, bool const& rooted, string const& lprefix);
    void get_newicks(ostream* poos);
    //~TopologyGenerator();
};

#endif /* _TGEN_H_ */
