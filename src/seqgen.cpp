/*
 * seqgen.cpp
 *
 *  Created on: Jun 23, 2015
 *      Author: joe
 */

#include <iostream>
#include <string>
#include <vector>
#include <nlopt.hpp>
#include <armadillo>
#include <random>
#include <numeric>

using namespace std;
using namespace arma;

#include "seqgen.h"
#include "utils.h"
#include "node.h"
#include "tree.h"
#include "tree_reader.h"

// TODO: do we want this order?
/*Default Rate Matrix looks like this, I don't know why but I always go A,T,C,G
 *
 *    A   T   C   G
 * A -1  .33 .33 .33
 * T .33 -1  .33 .33
 * G .33 .33  -1 .33
 * G .33 .33 .33  -1
 */


// this should enable easy changing of order, if desired

map <char, int> SequenceGenerator::nucMap = {
   {'A', 0},
   {'C', 1},
   {'G', 2},
   {'T', 3}
};
string SequenceGenerator::nucleotides = "ACGT";

/*
map <char, int> SequenceGenerator::nucMap = {
   {'A', 0},
   {'T', 1},
   {'C', 2},
   {'G', 3}
};
string SequenceGenerator::nucleotides = "ATCG";
*/
 
/* Use the P matrix probabilities and randomly draw numbers to see
 * if each individual state will undergo some type of change
 */
string SequenceGenerator::simulate_sequence (string const& anc, 
    vector < vector <double> >& QMatrix, float const& brlength) {

    std::vector<double>::iterator low;
    vector < vector <double> > PMatrix(4, vector <double>(4, 0.0));
    
    string newstring = anc; // instead of building, set size and replace
    for (int i = 0; i < seqlen; i++) {
        float RandNumb = get_uniform_random_deviate();
        int ancChar = nucMap[anc[i]];
        float brnew = brlength * site_rates[i];
        PMatrix = calculate_p_matrix(QMatrix, brnew);
        for (int i = 0; i < 4; i++) {
            // this calculates a cumulative sum
            std::partial_sum(PMatrix[i].begin(), PMatrix[i].end(), PMatrix[i].begin(), plus<double>());
        }
        low = std::lower_bound (PMatrix[ancChar].begin(), PMatrix[ancChar].end(), RandNumb);
        newstring[i] = nucleotides[low - PMatrix[ancChar].begin()];
    }
    return newstring;
}


/*
 * Calculate the Q Matrix (Substitution rate matrix)
 */
vector < vector <double> > SequenceGenerator::calculate_q_matrix () {

    vector < vector <double> > bigpi(4, vector <double>(4, 1.0));
    vector < vector <double> > t(4, vector <double>(4, 0.0));
    
    double tscale = 0.0;
    
    // doing the same looping multiple times here. simplify?
    for (unsigned int i = 0; i < rmatrix.size(); i++) {
        for (unsigned int j = 0; j < rmatrix.size(); j++) {
            if (i != j) {
                bigpi[i][j] *= basefreqs[i] * basefreqs[j] * rmatrix[i][j];
                tscale += bigpi[i][j];
            } else {
                bigpi[i][j] = 0.0;
            }
        }
    }
    for (unsigned int i = 0; i < rmatrix.size(); i++) {
        for (unsigned int j = 0; j < rmatrix.size(); j++) {
            if (i != j) {
                bigpi[i][j] /= tscale;
            } else {
                // set the diagnols to zero *** are they not set to zero above?
                bigpi[i][j] = 0.0;
            }
        }
    }
    for (unsigned int i = 0; i < rmatrix.size(); i++) {
        double diag = 0.0;
        for (unsigned int j = 0; j < rmatrix.size(); j++) {
            if (i != j) {
                diag -= bigpi[i][j];
            }
        }
        bigpi[i][i] = diag;
    }
    //Divide and Transpose
    for (unsigned int i = 0; i < rmatrix.size(); i++) {
        for (unsigned int j = 0; j < rmatrix.size(); j++) {
            bigpi[i][j] /= basefreqs[i];
        }
    }
    return bigpi;
}


/* Calculate the P Matrix (Probability Matrix)
 * Changes to armadillos format then back I don't like the way could be more
 * efficient but yeah...
 */
vector < vector <double> > SequenceGenerator::calculate_p_matrix (vector < vector <double> > QMatrix, float br) {

    vector < vector <double> > Pmatrix(4, vector <double>(4, 0.0));
    mat A = randn<mat>(4,4);
    mat B = randn<mat>(4,4);
    int count = 0;
    //Q * t moved into Matrix form for armadillo
    for (unsigned int i = 0; i < QMatrix.size(); i++) {
        for (unsigned int j = 0; j < QMatrix.size(); j++) {
            A[count] = (QMatrix[i][j] * br);
            count++;
        }
    }
   //exponentiate the matrix
   B = expmat(A);
   //cout << B << endl;
   count = 0;
   //convert the matrix back to C++ vector
   for (unsigned int i = 0; i < Pmatrix.size(); i++) {
        for (unsigned int j = 0; j < Pmatrix.size(); j++) {
            Pmatrix[i][j] = B[count];
            count++;
        }
   }
   return Pmatrix;
}


/*
 * Pre-Order traversal works
 * Calculates the JC Matrix
 */
// TODO: how to name ancestor nodes (sequences)
//       - if we have this we can add to results (if desired))
void SequenceGenerator::preorder_tree_traversal (Tree * tree, bool showancs, vector<double> multirates, bool mm) {

    double brlength = 0.0;
    //int count = 0; // not used for anything
    int rate_count = 0;
    string str = "";
    int check = 0;
    vector < vector <double> > QMatrix(4, vector <double>(4, 0.0));
    vector< vector <double> > rate_matrix(4, vector<double>(4, 0.33));

    //vector < vector <double> > PMatrix(4, vector <double>(4, 0.0));
    // NOTE: this uses order: A,T,C,G
    if (mm == true) {        
        rmatrix[0][2] = multirates[0]; //A->C
        rmatrix[2][0] = multirates[0]; //C->A
        rmatrix[0][3] = multirates[1]; //A->G
        rmatrix[3][0] = multirates[1]; //G->A
        rmatrix[0][1] = multirates[2]; //A->T
        rmatrix[1][0] = multirates[2]; //T->A
        rmatrix[2][3] = multirates[3]; //C->G
        rmatrix[3][2] = multirates[3]; //G->C
        rmatrix[1][2] = multirates[4]; //C->T
        rmatrix[2][1] = multirates[4]; //T->C
        rmatrix[1][3] = multirates[5]; //G->T
        rmatrix[3][1] = multirates[5]; //T->G
        rmatrix[0][0] = (multirates[0]+multirates[1]+multirates[2]) * -1;
        rmatrix[1][1] = (multirates[2]+multirates[4]+multirates[5]) * -1;
        rmatrix[2][2] = (multirates[0]+multirates[3]+multirates[4]) * -1;
        rmatrix[3][3] = (multirates[1]+multirates[3]+multirates[5]) * -1;
        for (unsigned int i = 0; i < 6; i++) {
            multirates.erase (multirates.begin() + 0); 
        }
        //QMatrix = calcQmatrix(rate_matrix);
        //QMatrix = calculate_q_matrix();
    }
    QMatrix = calculate_q_matrix();

    Node * root = tree->getRoot();
    seqs[root] = rootSequence;
    ancq[root] = QMatrix;
    
    if (showancs) {
        string tname = root->getName();
        Sequence seq(tname, rootSequence);
        res.push_back(seq);
    }
    
    // Pre-Order Traverse the tree
    for (int k = (tree->getNodeCount() - 2); k >= 0; k--) {
        brlength = tree->getNode(k)->getBL();
        /*for (unsigned int i = 0; i < QMatrix.size(); i++) {
            for (unsigned int j = 0; j < QMatrix.size(); j++) {
                cout << QMatrix[i][j] << " ";
            }
            cout << "\n";
        }
        cout << "\n";*/
        if (mm == true) {
            check = (int)round(multirates[0]);
            //cout << check << " " << rate_count << endl;
            if (tree->getNode(k)->isInternal() == true && multirates.size() != 0) {
                if (check == rate_count) {
                    rmatrix[0][2] = multirates[1];
                    rmatrix[2][0] = multirates[1];
                    rmatrix[0][3] = multirates[2];
                    rmatrix[3][0] = multirates[2];
                    rmatrix[0][1] = multirates[3];
                    rmatrix[1][0] = multirates[3];
                    rmatrix[2][3] = multirates[4];
                    rmatrix[3][2] = multirates[4];
                    rmatrix[1][2] = multirates[5];
                    rmatrix[2][1] = multirates[5];
                    rmatrix[1][3] = multirates[6];
                    rmatrix[3][1] = multirates[6];
                    rmatrix[0][0] = (multirates[1]+multirates[2]+multirates[3]) * -1;
                    rmatrix[1][1] = (multirates[3]+multirates[5]+multirates[6]) * -1;
                    rmatrix[2][2] = (multirates[1]+multirates[4]+multirates[5]) * -1;
                    rmatrix[3][3] = (multirates[2]+multirates[4]+multirates[6]) * -1;

                    for (unsigned int i = 0; i < 7; i++){
                        multirates.erase(multirates.begin() + 0);
                    }
                    //cout << "Size " << multirates.size() << endl;
                    //cout << multirates[0] << endl;
                    /*
                    for (unsigned int i = 0; i < multirates.size(); i++){
                        cout << multirates[i] << endl;
                    }*/        
                }
                rate_count++;
            }
        }
        QMatrix = calculate_q_matrix();

        //PMatrix = calculate_p_matrix(QMatrix, brlength);
        Node * dec = tree->getNode(k);
        Node * parent = tree->getNode(k)->getParent();
        //ancq[dec] = QMatrix;
        vector < vector <double> > Qparent = ancq[parent];
        string ancSeq = seqs[parent];
        string decSeq = simulate_sequence(ancSeq, Qparent, brlength);
        /*
        for (unsigned int i = 0; i < Qparent.size(); i++) {
            for (unsigned int j = 0; j < Qparent.size(); j++) {
                cout << Qparent[i][j] << " ";
            }
            cout << "\n";
        }
        cout << "\n";*/
        
        seqs[dec] = decSeq;
        ancq[dec] = QMatrix; // why store this?
        
        if (showancs == true && tree->getNode(k)->isInternal() == true) {
            string tname = tree->getNode(k)->getName();
            Sequence seq(tname, decSeq);
            res.push_back(seq);
        }
        // If its a tip print the name and the sequence
        if (tree->getNode(k)->isInternal() != true) {
            string tname = tree->getNode(k)->getName();
            Sequence seq(tname, decSeq);
            res.push_back(seq);
        }
    }
}


// this should probably be returned on its own
void SequenceGenerator::print_node_labels() {
    cout << tree->getRoot()->getNewick(true) <<";" << endl;
}


void SequenceGenerator::label_internal_nodes() {
    int count = 1;
    string str = "Node";
    string nlabel = "";
    Node * root = tree->getRoot();
    root->setName("Node_0");
    for (int k = (tree->getNodeCount() - 2); k >= 0; k--) {
        if (tree->getNode(k)->isInternal() == true) {
            //cout << k << endl;
            str = to_string(count);
            nlabel = "Node_" + str;
            tree->getNode(k)->setName(nlabel);
            count++;
        }
    }
}


// involves both gamma and pinvar
vector <float> SequenceGenerator::set_site_rates () {
    vector <float> srates(seqlen, 1.0);
    
    // invariable sites
    if (pinvar != 0.0) {
        int numsample = seqlen * pinvar + 0.5;
        // sample invariable sites
        vector <int> randsites = sample_without_replacement(seqlen, numsample);
        // must be a more elegant way of doing this
        for (int i = 0; i < numsample; i++) {
            srates[randsites[i]] = 0.0;
        }
    }
    
    // gamma-distributed rate variation. could explore other distributions...
    if (alpha != -1.0) { // default i.e. no rate variation
        for (int i = 0; i < seqlen; i++) {
            // want to skip over sites that are set to invariable
            if (srates[i] != 0.0) {
                srates[i] = get_gamma_random_deviate(alpha);
            }
        }
    }
    return srates;
} 


// initialized as string of length seqlength, all 'G'
string SequenceGenerator::generate_random_sequence () {
    
    string ancseq(seqlen, 'G');
    vector <double> cumsum(4);
    std::vector <double>::iterator low;
    // cumulative sum
    std::partial_sum(basefreqs.begin(), basefreqs.end(), cumsum.begin(), plus<double>());
    
    for (int i = 0; i < seqlen; i++) {
        float RandNumb = get_uniform_random_deviate();
        low = std::lower_bound (cumsum.begin(), cumsum.end(), RandNumb);
        ancseq[i] = nucleotides[low - cumsum.begin()];
    }
    return ancseq;
}


// rates are in order: A<->C,A<->G,A<->T,C<->G,C<->T,G<->T
vector < vector <double> > SequenceGenerator::construct_rate_matrix (vector <double> const& rates) {
    
    // initialize
    vector < vector <double> > ratemat(4, vector<double>(4, 0.33));
    
    // planning ahead here for potential non-reversible matrices
    if (rates.size() == 6) {
        ratemat[0][1] = rates[0];
        ratemat[1][0] = rates[0];
        ratemat[0][2] = rates[1];
        ratemat[2][0] = rates[1];
        ratemat[0][3] = rates[2];
        ratemat[3][0] = rates[2];
        ratemat[1][2] = rates[3];
        ratemat[2][1] = rates[3];
        ratemat[1][3] = rates[4];
        ratemat[3][1] = rates[4];
        ratemat[2][3] = rates[5];
        ratemat[3][2] = rates[5];
        
        ratemat[0][0] = (rates[0] + rates[1] + rates[2]) * -1;
        ratemat[1][1] = (rates[0] + rates[3] + rates[4]) * -1;
        ratemat[2][2] = (rates[1] + rates[3] + rates[5]) * -1;
        ratemat[3][3] = (rates[2] + rates[4] + rates[5]) * -1;
        
    } else {
        cout << "Er, we don't deal with " << rates.size() << " rates at the moment..." << endl;
        exit(0);
    }
    return ratemat;
}


// set all values
void SequenceGenerator::initialize () {
    // set the number generator being used
    if (seed != -1) { // user provided seed
        generator = mt19937(seed);
    } else {
        //random_device rand_dev;
        generator = mt19937(get_clock_seed());
    }
    if (showancs) {
        label_internal_nodes();
    }
    if (rootSequence.length() == 0) {
        rootSequence = generate_random_sequence();
    } else {
        check_valid_sequence();
        // if root sequence is provided, set length to this
        seqlen = rootSequence.size();
    }
    // set site-specific rate (pinvar and gamma)
    site_rates = set_site_rates();
}


// make sure sequence contains only valid nucleotide characters
void SequenceGenerator::check_valid_sequence () {
    // make sure uppercase
    std::transform(rootSequence.begin(), rootSequence.end(), rootSequence.begin(), ::toupper);
    std::size_t found = rootSequence.find_first_not_of(nucleotides);
    if (found != std::string::npos) {
        cout << "Error: illegal character '" << rootSequence[found] << "' at position " 
            << found+1 << " (only A,C,G,T allowed). Exiting." << endl;
        exit(0);
    }
}


SequenceGenerator::SequenceGenerator (int const &seqlength, vector <double> const& basefreq,
    vector < vector<double> >& rmatrix, Tree * tree, bool const& showancs, 
    int const& nreps, int const& seed, float const& alpha, float const& pinvar, 
    string const& ancseq, bool const& printpost, vector<double> const& multirates,
    bool const& mm):tree(tree), seqlen(seqlength), nreps(nreps), seed(seed), alpha(alpha),
    pinvar(pinvar), rootSequence(ancseq), basefreqs(basefreq), rmatrix(rmatrix), 
    multirates(multirates), showancs(showancs), printnodelabels(printpost),
    mm(mm) {
    
    initialize();
    
    // Print out the nodes names
    if (printnodelabels == true) {
        print_node_labels();
    }
    
    preorder_tree_traversal(tree, showancs, multirates, mm);
}


// not sure of a more elegant way to do this...
vector<Sequence> SequenceGenerator::get_sequences () {
    return res;
}


// call this whenever a random float is needed
float SequenceGenerator::get_uniform_random_deviate () {
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    return distribution(generator);
}


float SequenceGenerator::get_gamma_random_deviate (float alpha) {

    //default_random_engine generator;
    std::gamma_distribution<float> distribution(alpha,(1/alpha));
    return distribution(generator);
}

//SEQGEN::~SEQGEN() {
//    // TODO Auto-generated destructor stub
//}
