#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
//#include <ctime>
//#include <random>
#include <map>

using namespace std;

#include "utils.h"
#include "tgen.h"
#include "tree_utils.h"

TopologyGenerator::TopologyGenerator(int const& ntax, bool const& rooted):ntax_(ntax),
        rooted_(rooted) {
    initialize();
}

// generate starting star tree (edge matrix), set starting taxon and node
void TopologyGenerator::initialize () {
    vector < vector <int> > init_edges_;
    ntopos_ = get_num_possible_trees(ntax_, rooted_);
    nedges_ = get_num_edges(ntax_, rooted_);
    if (!rooted_) {
        init_edges_ = initialize_edge_matrix_unrooted(ntax_);
        cout << "There are " << ntopos_ << " possible unrooted topologies for "
                << ntax_ << " taxa." << endl;
        curtax_ = 4;
    } else {
        cout << "There are " << ntopos_ << " possible rooted topologies for "
                << ntax_ << " taxa." << endl;
        cout << "This will be implemented shortly." << endl;
        curtax_ = 3;
        exit(0);
    }
    curnode_ = init_edges_[0][0] + 1; // increment largest initialized node
    trees_.push_back(init_edges_);
}

// the number of edges in the final trees
// 2n-3 for unrooted, 2n-2 for rooted
int TopologyGenerator::get_num_edges (int const& ntax, int const& rooted) {
    int nedges = 2 * ntax - 3 + (int)rooted;
    return nedges;
}

// trees are stored as edge matrices - avoids overhead of node objects
// 2 columns: ancestor node index, descendant node index
// 2n-3 rows (edges)
vector < vector <int> > TopologyGenerator::initialize_edge_matrix_unrooted (int const& n) {
    int num_edges = 2 * n - 3;
    vector < vector <int> > edges(num_edges, vector<int>(2, 0));
    edges[0][0] = edges[1][0] = edges[2][0] = n + 1;
    // initialize 3 taxon star tree
    for (int i = 0; i < 3; i++) {
        edges[i][1] = i + 1;
    }
    return edges;
}

// add next taxon to existing edge matrices
// for the nth taxon, there are 2(n-1)-3 possible attachment points
// - i.e., the edges of the previous iteration
vector < vector < vector <int> > > TopologyGenerator::add_taxon_unrooted (vector < vector < vector <int> > > edges,
        int const& taxon, int const& new_node) {
    vector < vector < vector <int> > > trees;
    int num_attach = 2 * (taxon - 1) - 3; // also the row where next edge will be added
    vector < vector <int> > edge;
    for (int j = 0; j < (int)edges.size(); j++) {
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

void TopologyGenerator::generate_trees () {
    while (curtax_ <= ntax_) {
        //cout << "adding taxon: " << curtax_ << endl;
        trees_ = add_taxon_unrooted (trees_, curtax_, curnode_);
        curtax_++;
        curnode_++;
    }
    //cout << "done." << endl;
}

string TopologyGenerator::edge_matrix_to_newick (vector < vector <int> > const& edges) {
    string tree = "";
    int num_edges = (int)edges.size();
    int ntax = (int)((num_edges + 3)/2);
    int min_node = ntax + 1;
    int max_node = num_edges + 1;
    
    map <int, vector <int>> m;
    
    for (int i = min_node; i <= max_node; i++) {
        vector <int> clade;
        for (int j = 0; j < num_edges; j++) {
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
void TopologyGenerator::newick_from_tree_map (int node, map <int, vector <int>> m, string & tree) {
    tree += "("; // every entrance to this function is a tree
    
    vector <int> decnodes = m[node];
    int numdec = (int)decnodes.size();
    
    for (int i = 0; i < numdec; i++) {
        int curnode = decnodes[i];
        // is curnode a clade (i.e., needs to be further processed)?
        // could also use number of tips if things are labelled reasonably
        if (m.count(curnode)) {
            newick_from_tree_map(curnode, m, tree); // recursion, baby
        } else {
            tree += to_string(curnode);
        }
        if (i < (numdec - 1)) {
            tree += ",";
        }
    }
    tree += ")";
}

// this will send newick to poos
void TopologyGenerator::get_newicks (ostream* poos) {
    generate_trees();
    vector < vector <int> > tree;
    for (int i = 0; i < ntopos_; i++) {
        tree = trees_[i];
        string newick = edge_matrix_to_newick(tree);
        cout << newick << endl;
    }
}
