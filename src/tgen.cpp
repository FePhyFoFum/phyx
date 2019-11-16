#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
//#include <ctime>
//#include <random>
#include <map>

#include "utils.h"
#include "tgen.h"
#include "tree_utils.h"

TopologyGenerator::TopologyGenerator(const int& ntax, const bool& rooted,
        const std::string& lprefix):ntax_(ntax), rooted_(rooted), lprefix_(lprefix) {
    initialize();
}


// generate starting star tree (edge matrix), set starting taxon and node
void TopologyGenerator::initialize () {
    std::vector< std::vector<int> > init_edges_;
    ntopos_ = get_num_possible_trees(ntax_, rooted_);
    nedges_ = get_num_edges(ntax_, rooted_);
    if (!rooted_) {
        init_edges_ = initialize_edge_matrix_unrooted(ntax_);
        curtax_ = 4;
    } else {
        init_edges_ = initialize_edge_matrix_rooted (ntax_);
        curtax_ = 3;
    }
    curnode_ = init_edges_[0][0] + 1; // increment largest initialized node
    trees_.push_back(init_edges_);
}


// the number of edges in the final trees
// 2n-3 for unrooted, 2n-2 for rooted
int TopologyGenerator::get_num_edges (const int& ntax, const int& rooted) {
    int nedges = 2 * ntax - 3 + (int)rooted;
    return nedges;
}


// trees are stored as edge matrices - avoids overhead of node objects
// 2 columns: ancestor node index, descendant node index
// 2n-3 rows (edges)
// first 3 edges are initialized (i.e., star tree)
std::vector< std::vector<int> > TopologyGenerator::initialize_edge_matrix_unrooted (const int& n) {
    int num_edges = 2 * n - 3;
    std::vector< std::vector<int> > edges(num_edges, std::vector<int>(2, 0));
    edges[0][0] = edges[1][0] = edges[2][0] = n + 1;
    // initialize 3 taxon star tree
    for (int i = 0; i < 3; i++) {
        edges[i][1] = i + 1;
    }
    return edges;
}


// rooted version of above
// 2 columns: ancestor node index, descendant node index
// 2n-2 rows (edges)
// first 2 edges are initialized
std::vector< std::vector<int> > TopologyGenerator::initialize_edge_matrix_rooted (const int& n) {
    int num_edges = 2 * n - 2;
    std::vector< std::vector<int> > edges(num_edges, std::vector<int>(2, 0));
    edges[0][0] = edges[1][0] = n + 1;
    // initialize 2 taxon tree
    for (int i = 0; i < 2; i++) {
        edges[i][1] = i + 1;
    }
    return edges;
}


// add next taxon to existing edge matrices
// for the nth taxon, there are 2(n-1)-3 possible attachment points
// - i.e., the edges of the previous iteration
std::vector< std::vector< std::vector<int> > > TopologyGenerator::add_taxon_unrooted (std::vector < std::vector< std::vector <int> > > edges,
        const int& taxon, const int& new_node) {
    std::vector< std::vector< std::vector<int> > > trees;
    int num_attach = 2 * (taxon - 1) - 3; // also the row where next edge will be added
    std::vector< std::vector<int> > edge;
    for (unsigned int j = 0; j < edges.size(); j++) {
        for (int i = 0; i < num_attach; i++) {
            edge = edges[j];
            int prev_node = edge[i][0];
            edge[i][0] = new_node; // dec stays the same
            // new edges
            edge[num_attach][0] = new_node;
            edge[num_attach][1] = taxon; // new tip
            
            edge[num_attach + 1][0] = prev_node;
            edge[num_attach + 1][1] = new_node;
            
            trees.push_back(edge);
        }
    }
    return trees;
}


// rooted version of above
// add next taxon to existing edge matrices
// for the nth taxon, there are 2(n-1)-2 possible attachment points
// - i.e., the edges of the previous iteration
// need to keep a copy of the first matrix for flipping in the last
std::vector< std::vector< std::vector<int> > > TopologyGenerator::add_taxon_rooted (std::vector< std::vector< std::vector<int> > > edges,
        const int& taxon, const int& new_node) {
    std::vector< std::vector< std::vector<int> > > trees;
    int num_attach = 2 * (taxon - 1) - 2; // also the row where next edge will be added
    std::vector< std::vector<int> > edge;
    
    for (unsigned int j = 0; j < edges.size(); j++) {
        for (int i = 0; i < num_attach; i++) {
            edge = edges[j];
            int prev_node = edge[i][0];
            edge[i][0] = new_node; // dec stays the same
            // new edges
            edge[num_attach][0] = new_node;
            edge[num_attach][1] = taxon; // new tip
            
            edge[num_attach + 1][0] = prev_node;
            edge[num_attach + 1][1] = new_node;
            
            trees.push_back(edge);
        }
        // now for the flipped matrix
        edge = edges[j];
        // first, swap new node in for root node (always ntax+1)
        int root = ntax_ + 1;
        for (unsigned int i = 0; i < edge.size(); i++) {
            if (edge[i][0] == root) {
                edge[i][0] = new_node;
            } else if (edge[i][0] == 0) {
                // early break: if 0, then edge has yet to be generated
                break;
            }
        }
        // second, add original root back in, attching to new taxon
        edge[num_attach][0] = root;
        edge[num_attach][1] = taxon;
        // finally, attach root above new node
        edge[num_attach + 1][0] = root;
        edge[num_attach + 1][1] = new_node;
        // and add to the tree collection
        trees.push_back(edge);
        
    }
    return trees;
}


void TopologyGenerator::generate_trees () {
    while (curtax_ <= ntax_) {
        //cout << "adding taxon: " << curtax_ << std::endl;
        if (!rooted_) {
            trees_ = add_taxon_unrooted (trees_, curtax_, curnode_);
        } else {
            trees_ = add_taxon_rooted (trees_, curtax_, curnode_);
        }
        curtax_++;
        curnode_++;
    }
    //cout << "done." << std::endl;
}


std::string TopologyGenerator::edge_matrix_to_newick (const std::vector< std::vector<int> >& edges) {
    std::string tree = "";
    //int num_edges = (int)edges.size();
    //int ntax = (int)((num_edges + 3)/2);
    int min_node = ntax_ + 1;
    int max_node = nedges_ + 1;
    
    std::map<int, std::vector<int>> m;
    
    for (int i = min_node; i <= max_node; i++) {
        std::vector<int> clade;
        for (int j = 0; j < nedges_; j++) {
            if (edges[j][0] == i) {
                clade.push_back(edges[j][1]);
            }
        }
        m[i] = clade;
    }
    newick_from_tree_map(min_node, m, tree);
    tree += ";";
    return tree;
}


// recursive, so newick string is passed by reference
void TopologyGenerator::newick_from_tree_map (int node, std::map<int, std::vector<int>> m, std::string& tree) {
    tree += "("; // every entrance to this function is a tree
    
    std::vector<int> decnodes = m[node];
    int numdec = (int)decnodes.size();
    
    for (int i = 0; i < numdec; i++) {
        int curnode = decnodes[i];
        // is curnode a clade (i.e., needs to be further processed)?
        // could also use number of tips if things are labelled reasonably
        if (m.count(curnode)) {
            newick_from_tree_map(curnode, m, tree); // recursion, baby
        } else {
            // terminal
            tree += lprefix_ + std::to_string(curnode);
        }
        if (i < (numdec - 1)) {
            tree += ",";
        }
    }
    tree += ")";
}


// this will send newick to poos
void TopologyGenerator::get_newicks (std::ostream* poos) {
    generate_trees();
    std::vector < std::vector<int> > tree;
    for (unsigned long int i = 0; i < ntopos_; i++) {
        tree = trees_[i];
        std::string newick = edge_matrix_to_newick(tree);
        (*poos) << newick << std::endl;
    }
}
