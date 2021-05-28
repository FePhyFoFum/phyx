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


TopologyGenerator::TopologyGenerator(const unsigned int& num_taxa, const bool& rooted,
        const std::string& lprefix):num_taxa_(num_taxa), rooted_(rooted), lprefix_(lprefix),
        ntopos_(0), nedges_(0), curtax_(0), curnode_(0) {
    initialize();
}


// generate starting star tree (edge matrix), set starting taxon and node
void TopologyGenerator::initialize () {
    std::vector< std::vector<unsigned int> > init_edges_;
    ntopos_ = get_num_possible_trees(num_taxa_, rooted_);
    nedges_ = get_num_edges();
    if (!rooted_) {
        init_edges_ = initialize_edge_matrix_unrooted();
        curtax_ = 4u;
    } else {
        init_edges_ = initialize_edge_matrix_rooted();
        curtax_ = 3u;
    }
    curnode_ = init_edges_[0][0] + 1; // increment largest initialized node
    trees_.push_back(init_edges_);
}


// the number of edges in the final trees
// 2n-3 for unrooted, 2n-2 for rooted
unsigned int TopologyGenerator::get_num_edges () {
    unsigned int nedges = 2 * num_taxa_ - 3 + static_cast<unsigned int>(rooted_);
    return nedges;
}


// trees are stored as edge matrices - avoids overhead of node objects
// 2 columns: ancestor node index, descendant node index
// 2n-3 rows (edges)
// first 3 edges are initialized (i.e., star tree)
std::vector< std::vector<unsigned int> > TopologyGenerator::initialize_edge_matrix_unrooted () {
    std::vector< std::vector<unsigned int> > edges(nedges_, std::vector<unsigned int>(2, 0));
    edges[0][0] = edges[1][0] = edges[2][0] = num_taxa_ + 1;
    // initialize 3 taxon star tree
    for (unsigned int i = 0; i < 3; i++) {
        edges[i][1] = i + 1;
    }
    return edges;
}


// rooted version of above
// 2 columns: ancestor node index, descendant node index
// 2n-2 rows (edges)
// first 2 edges are initialized
std::vector< std::vector<unsigned int> > TopologyGenerator::initialize_edge_matrix_rooted () {
    std::vector< std::vector<unsigned int> > edges(nedges_, std::vector<unsigned int>(2, 0));
    edges[0][0] = edges[1][0] = num_taxa_ + 1;
    // initialize 2 taxon tree
    for (unsigned int i = 0; i < 2; i++) {
        edges[i][1] = i + 1;
    }
    return edges;
}


// add next taxon to existing edge matrices
// for the nth taxon, there are 2(n-1)-3 possible attachment points
// - i.e., the edges of the previous iteration
std::vector< std::vector< std::vector<unsigned int> > > TopologyGenerator::add_taxon_unrooted (
        std::vector< std::vector< std::vector<unsigned int> > > edges,
        const unsigned int& taxon, const unsigned int& new_node) {
    std::vector< std::vector< std::vector<unsigned int> > > trees;
    unsigned int num_attach = 2 * (taxon - 1) - 3; // also the row where next edge will be added
    std::vector< std::vector<unsigned int> > edge;
    for (const auto & j : edges) {
        for (unsigned int i = 0; i < num_attach; i++) {
            edge = j;
            unsigned int prev_node = edge[i][0];
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
std::vector< std::vector< std::vector<unsigned int> > > TopologyGenerator::add_taxon_rooted (
        std::vector< std::vector< std::vector<unsigned int> > > edges,
        const unsigned int& taxon, const unsigned int& new_node) {
    std::vector< std::vector< std::vector<unsigned int> > > trees;
    unsigned int num_attach = 2 * (taxon - 1) - 2; // also the row where next edge will be added
    std::vector< std::vector<unsigned int> > edge;
    
    for (const auto & j : edges) {
        for (unsigned int i = 0; i < num_attach; i++) {
            edge = j;
            unsigned int prev_node = edge[i][0];
            edge[i][0] = new_node; // dec stays the same
            // new edges
            edge[num_attach][0] = new_node;
            edge[num_attach][1] = taxon; // new tip
            
            edge[num_attach + 1][0] = prev_node;
            edge[num_attach + 1][1] = new_node;
            
            trees.push_back(edge);
        }
        // now for the flipped matrix
        edge = j;
        // first, swap new node in for root node (always num_taxa+1)
        unsigned int root = num_taxa_ + 1;
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
    while (curtax_ <= num_taxa_) {
        //std::cout << "adding taxon: " << curtax_ << std::endl;
        if (!rooted_) {
            trees_ = add_taxon_unrooted(trees_, curtax_, curnode_);
        } else {
            trees_ = add_taxon_rooted(trees_, curtax_, curnode_);
        }
        curtax_++;
        curnode_++;
    }
    //std::cout << "done." << std::endl;
}


std::string TopologyGenerator::edge_matrix_to_newick (const std::vector< std::vector<unsigned int> >& edges) {
    std::string tree;
    //int num_edges = (int)edges.size();
    //int num_taxa = (int)((num_edges + 3)/2);
    unsigned int min_node = num_taxa_ + 1;
    unsigned int max_node = nedges_ + 1;
    
    std::map<unsigned int, std::vector<unsigned int>> m;
    
    for (unsigned int i = min_node; i <= max_node; i++) {
        std::vector<unsigned int> clade;
        for (unsigned int j = 0; j < nedges_; j++) {
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
void TopologyGenerator::newick_from_tree_map (unsigned int node,
        std::map<unsigned int, std::vector<unsigned int>> m, std::string& tree) {
    tree += "("; // every entrance to this function is a tree
    
    std::vector<unsigned int> decnodes = m[node];
    auto numdec = static_cast<unsigned int>(decnodes.size());
    
    for (unsigned int i = 0; i < numdec; i++) {
        unsigned int curnode = decnodes[i];
        // is curnode a clade (i.e., needs to be further processed)?
        // could also use number of tips if things are labelled reasonably
        if (m.count(curnode) != 0u) {
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
    std::vector< std::vector<unsigned int> > tree;
    for (unsigned long int i = 0; i < ntopos_; i++) {
        tree = trees_[i];
        std::string newick = edge_matrix_to_newick(tree);
        (*poos) << newick << std::endl;
    }
}