#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <getopt.h>

#include "upgma.h"
#include "utils.h"
#include "seq_utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "node.h"
#include "tree.h"
#include "tree_utils.h"


UPGMA::UPGMA (std::istream* pios):num_taxa_(0), num_char_(0), newickstring_("") {
    std::string alphaName; // not used, but required by reader
    seqs_ = ingest_alignment(pios, alphaName);
    num_taxa_ = static_cast<int>(seqs_.size());
    num_char_ = static_cast<int>(seqs_[0].get_length());
    
    // check that it is aligned (doesn't make sense otherwise)
    if (!is_aligned(seqs_)) {
        std::cerr << "Error: sequences are not aligned. Exiting." << std::endl;
        exit(0);
    }
    names_ = collect_names(seqs_);
    full_distmatrix_ = build_matrix();
}


std::vector< std::vector<double> > UPGMA::build_matrix () {
    // 1) skip self comparisons
    // 2) only calculate one half of matrix (i.e., no duplicate calcs)
    std::vector< std::vector<double> > distances(num_taxa_, std::vector<double>(num_taxa_, 0.0));
    
    double tempScore = 0.0;
    for (int i = 0; i < num_taxa_; i++) {
        std::string seq1 = seqs_[i].get_sequence();
        for (int j = (i + 1); j < num_taxa_; j++) {
            std::string seq2 = seqs_[j].get_sequence();
            // get distance
            tempScore = static_cast<double>(calc_hamming_dist(seq1, seq2));
            // put scale in terms of number of sites. original version did not do this
            tempScore /= static_cast<double>(num_char_);
            // put in both top and bottom of matrix, even though only top is used
            distances[i][j] = distances[j][i] = tempScore;
        }
    }
    
    // just for debugging
    /*
    std::cout << "\t";
    for (unsigned int i = 0; i < names_.size(); i++) {
        std::cout << names_[i] << "\t";
    }
    std::cout << std::endl;
    for (int i = 0; i < num_taxa_; i++) {
        std::cout << names_[i] << "\t";
        for (int j = 0; j < num_taxa_; j++) {
            std::cout << distances[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    */
    
    return distances;
}


// find smallest pairwise distance
// will always find this on the top half of the matrix i.e., mini1 < mini2
double UPGMA::get_smallest_distance (const std::vector< std::vector<double> >& dmatrix, int& mini1, int& mini2) {
    // super large value
    double minD = 99999999999.99;
    int numseqs = dmatrix.size();
    for (int i = 0; i < (numseqs - 1); i++) {
        int idx = std::min_element(dmatrix[i].begin() + (i + 1), dmatrix[i].end()) - dmatrix[i].begin();
        if (dmatrix[i][idx] < minD) {
            minD = dmatrix[i][idx];
            mini1 = i;
            mini2 = idx;
        }
    }
    return minD;
}


void UPGMA::construct_tree () {
    // location of minimum distance (top half)
    int ind1 = 0;
    int ind2 = 0;
    double minD = 0.0;
    
    // initialize
    std::vector< std::vector<double> > dMatrix = full_distmatrix_;
    double newHeight = 0.0;
    int numClusters = num_taxa_;
    Node * anc = nullptr; // new node, ancestor of 2 clusters
    Node * left = nullptr;
    Node * right = nullptr;
    
    // keep list of nodes left to be clustered. initially all terminal nodes
    std::vector<Node *> nodes(num_taxa_);
    for (int i = 0; i < num_taxa_; i++) {
        Node * nd = new Node();
        nd->setName(names_[i]);
        nd->setHeight(0.0);
        nodes[i] = nd;
    }
    
    while (numClusters > 1) {
        // 1. get smallest distance present in the matrix
        minD = get_smallest_distance(dMatrix, ind1, ind2);
        left = nodes[ind1];
        right = nodes[ind2];
        
        // 2. create new ancestor node
        anc = new Node();
        
        // 3. add nodes in new cluster above as children to new ancestor
        anc->addChild(*left); // addChild calls setParent
        anc->addChild(*right);
        
        // 4. compute edgelengths: half of the distance
        // edgelengths must subtract the existing height
        newHeight = 0.5 * minD;
        left->setBL(newHeight - left->getHeight());
        right->setBL(newHeight - right->getHeight());
        
        // make sure to set the height of anc for the next iteration to use
        anc->setHeight(newHeight);
        
        // 5. compute new distance matrix (1 fewer rows & columns)
        // new distances are proportional averages (size of clusters)
        // new cluster is placed first (row & column)
        std::vector<double> avdists(numClusters, 0.0);
        double Lweight = left->isExternal() ? 1.0 : static_cast<double>(left->getChildCount());
        double Rweight = right->isExternal() ? 1.0 : static_cast<double>(right->getChildCount());
        for (int i = 0; i < numClusters; i++) {
            avdists[i] = ((dMatrix[ind1][i] * Lweight) + (dMatrix[ind2][i] * Rweight)) / (Lweight + Rweight);
        }
        
        numClusters--;
        std::vector< std::vector<double> > newDistances(numClusters, std::vector<double>(numClusters, 0.0));
        
        // put in distances to new clusters first
        double tempDist = 0.0;
        int count = 0;
        for (int i = 0; i < static_cast<int>(nodes.size()); i++) {
            if (i != ind1 && i != ind2) {
                count++;
                tempDist = avdists[i];
                newDistances[0][count] = tempDist;
                newDistances[count][0] = tempDist;
            }
        }
        
        // now, fill in remaining
        int icount = 1;
        int jcount = 1;
        int ndsize = static_cast<int>(nodes.size());
        for (int i = 0; i < ndsize; i++) {
            jcount = 1;
            if (i != ind1 && i != ind2) {
                for (int j = 0; j < ndsize; j++) {
                    if (j != ind1 && j != ind2) {
                        newDistances[icount][jcount] = dMatrix[i][j];
                        newDistances[jcount][icount] = dMatrix[i][j];
                        jcount++;
                    }
                }
                icount++;
            }
        }
        
        // replace distance matrix
        dMatrix = newDistances;
        
        // 6. finally, update node vector (1 shorter). new node always goes first)
        std::vector<Node *> newNodes(numClusters);
        newNodes[0] = anc;
        int counter = 1;
        for (int i = 0; i < ndsize; i++) {
            if (i != ind1 && i != ind2) {
                newNodes[counter] = nodes[i];
                counter++;
            }
        }
        // replace node vector
        nodes = newNodes;
    }
    tree_ = new Tree(anc);
    tree_->setEdgeLengthsPresent(true); // used by newick writer
}


std::string UPGMA::get_newick () {
    if (newickstring_.empty()) {
        construct_tree();
    }
    newickstring_ = getNewickString(tree_);
    return newickstring_;
}


std::vector< std::vector<double> > UPGMA::get_matrix () const {
    return full_distmatrix_;
}
