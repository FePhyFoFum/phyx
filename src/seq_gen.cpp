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
//using namespace arma;
using arma::randn;
using arma::mat;
using arma::expmat;

#include "seq_gen.h"
#include "utils.h"
#include "node.h"
#include "tree.h"
#include "tree_reader.h"
#include "tree_utils.h"

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

map <char, int> SequenceGenerator::nuc_map_ = {
   {'A', 0},
   {'C', 1},
   {'G', 2},
   {'T', 3}
};

string SequenceGenerator::nucleotides_ = "ACGT";
string SequenceGenerator::amino_acids_ = "ARNDCQEGHILKMFPSTWYV";

map <char, int> SequenceGenerator::aa_map_ = {

   {'A', 0},
   {'R', 1},
   {'N', 2},
   {'D', 3},
   {'C', 4},
   {'Q', 5},
   {'E', 6},
   {'G', 7},
   {'H', 8},
   {'I', 9},
   {'L', 10},
   {'K', 11},
   {'M', 12},
   {'F', 13},
   {'P', 14},
   {'S', 15},
   {'T', 16},
   {'W', 17},
   {'Y', 18},
   {'V', 19}
};

/*
map <char, int> SequenceGenerator::nucMap = {
   {'A', 0},
   {'T', 1},
   {'C', 2},
   {'G', 3}
};
string SequenceGenerator::nucleotides = "ATCG";
*/
 
SequenceGenerator::SequenceGenerator (int const &seqlength, vector <double> const& basefreq,
    vector < vector<double> >& rmatrix, Tree * tree, bool const& showancs, 
    int const& nreps, int const& seed, float const& alpha, float const& pinvar, 
    string const& ancseq, bool const& printpost, vector<double> const& multirates,
    vector <double> const& aabasefreq, bool const& is_dna):tree_(tree),
    seqlen_(seqlength), nreps_(nreps), seed_(seed), alpha_(alpha), pinvar_(pinvar),
    root_sequence_(ancseq), base_freqs_(basefreq), aa_freqs_(aabasefreq), rmatrix_(rmatrix), 
    multi_rates_(multirates), show_ancs_(showancs), print_node_labels_(printpost),
    multi_model_(false), is_dna_(is_dna)  {
    /*
     for (unsigned int i = 0; i < rmatrix.size(); i++) {
        for (unsigned int j = 0; j < rmatrix.size(); j++) {
            cout << rmatrix[i][j] << " ";
        }
        cout << "\n";
    }*/
    initialize();
    if (is_dna_) {
        nstates_ = 4;
    } else {
        nstates_ = 20;
    }
    // Print out the nodes names
    if (print_node_labels_) {
        label_internal_nodes();
        print_node_labels();
        exit(0);
    }
    preorder_tree_traversal();
}


// set all values
void SequenceGenerator::initialize () {
    // set the number generator being used
    if (seed_ != -1) { // user provided seed
        generator_ = mt19937(seed_);
    } else {
        generator_ = mt19937(get_clock_seed());
    }
    
    // construct uniform distribution from which random numbers will be generated
    // this happens to be the default distribution, but good to be explicit
    uniformDistrib_ = uniform_real_distribution<float>(0.0, 1.0);
    
    // construct gamma distribution (if necessary)
    if (alpha_ != -1.0) {
        gammaDistrib_ = gamma_distribution<float>(alpha_, (1/alpha_));
    }
    
    if (show_ancs_) {
        label_internal_nodes();
    }
    if (root_sequence_.length() == 0) {
        root_sequence_ = generate_random_sequence();
    } else {
        check_valid_sequence();
        // if root sequence is provided, set length to this
        seqlen_ = root_sequence_.size();
    }
    // set site-specific rate (pinvar and gamma)
    site_rates_ = set_site_rates();
    
    if (multi_rates_.size() != 0) {
        multi_model_ = true;
    }
}


/* Use the P matrix probabilities and randomly draw numbers to see
 * if each individual state will undergo some type of change
 */
string SequenceGenerator::simulate_sequence (string const& anc, 
    vector < vector <double> >& QMatrix, float const& brlength) {
    std::vector<double>::iterator low;
    vector < vector <double> > PMatrix(nstates_, vector <double>(nstates_, 0.0));
    //int ancChar = 0;
    string newstring = anc; // instead of building, set size and replace
    for (int i = 0; i < seqlen_; i++) {
        float RandNumb = get_uniform_random_deviate();
        int ancChar = 0;
        if (is_dna_) {
            ancChar = nuc_map_[anc[i]];
        } else {
            ancChar = aa_map_[anc[i]];
            //cout << ancChar << endl;
        }
        float brnew = brlength * site_rates_[i];
        PMatrix = calculate_p_matrix(QMatrix, brnew);
        for (int i = 0; i < nstates_; i++) {
            // this calculates a cumulative sum
            std::partial_sum(PMatrix[i].begin(), PMatrix[i].end(), PMatrix[i].begin(), plus<double>());
        }
        low = std::lower_bound (PMatrix[ancChar].begin(), PMatrix[ancChar].end(), RandNumb);
        
        if (is_dna_) {
            newstring[i] = nucleotides_[low - PMatrix[ancChar].begin()];
        } else {
            newstring[i] = amino_acids_[low - PMatrix[ancChar].begin()];    
        }
    }
    //cout << newstring << endl;
    return newstring;
}


/*
 * Calculate the Q Matrix (Substitution rate matrix)
 */
vector < vector <double> > SequenceGenerator::calculate_q_matrix () {

    vector < vector <double> > bigpi(nstates_, vector <double>(nstates_, 1.0));
    vector < vector <double> > t(nstates_, vector <double>(nstates_, 0.0));
    
    double tscale = 0.0;
    
    // doing the same looping multiple times here. simplify?
    for (unsigned int i = 0; i < rmatrix_.size(); i++) {
        for (unsigned int j = 0; j < rmatrix_.size(); j++) {
            if (i != j) {
                if (is_dna_) {
                    bigpi[i][j] *= base_freqs_[i] * base_freqs_[j] * rmatrix_[i][j];
                } else {
                    bigpi[i][j] *= aa_freqs_[i] * aa_freqs_[j] * rmatrix_[i][j];    
                }
                tscale += bigpi[i][j];
            } else {
                bigpi[i][j] = 0.0;
            }
        }
    }
    for (unsigned int i = 0; i < rmatrix_.size(); i++) {
        for (unsigned int j = 0; j < rmatrix_.size(); j++) {
            if (i != j) {
                bigpi[i][j] /= tscale;
            } else {
                // set the diagnols to zero *** are they not set to zero above?
                bigpi[i][j] = 0.0;
            }
        }
    }
    for (unsigned int i = 0; i < rmatrix_.size(); i++) {
        double diag = 0.0;
        for (unsigned int j = 0; j < rmatrix_.size(); j++) {
            if (i != j) {
                diag -= bigpi[i][j];
            }
        }
        bigpi[i][i] = diag;
    }
    //Divide and Transpose
    for (unsigned int i = 0; i < rmatrix_.size(); i++) {
        for (unsigned int j = 0; j < rmatrix_.size(); j++) {
            if (is_dna_) {
                bigpi[i][j] /= base_freqs_[i];
            } else {
                bigpi[i][j] /= aa_freqs_[i];
            }
        }
    }
    return bigpi;
}


/* Calculate the P Matrix (Probability Matrix)
 * Changes to armadillos format then back I don't like the way could be more
 * efficient but yeah...
 */
vector < vector <double> > SequenceGenerator::calculate_p_matrix (vector < vector <double> > const& QMatrix,
    float br) {

    vector < vector <double> > Pmatrix(nstates_, vector <double>(nstates_, 0.0));
    mat A = randn<mat>(nstates_, nstates_);
    mat B = randn<mat>(nstates_, nstates_); // why not just copy A?
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
void SequenceGenerator::preorder_tree_traversal () {

    double brlength = 0.0;
    int rate_count = 0;
    int check = 0;
    vector < vector <double> > QMatrix(nstates_, vector <double>(nstates_, 0.0));
    //vector < vector <double> > PMatrix(4, vector <double>(4, 0.0));
    // NOTE: this uses order: A,T,C,G
    if (multi_model_) {        
        rmatrix_[0][2] = multi_rates_[0]; // A->C
        rmatrix_[2][0] = multi_rates_[0]; // C->A
        rmatrix_[0][3] = multi_rates_[1]; // A->G
        rmatrix_[3][0] = multi_rates_[1]; // G->A
        rmatrix_[0][1] = multi_rates_[2]; // A->T
        rmatrix_[1][0] = multi_rates_[2]; // T->A
        rmatrix_[2][3] = multi_rates_[3]; // C->G
        rmatrix_[3][2] = multi_rates_[3]; // G->C
        rmatrix_[1][2] = multi_rates_[4]; // C->T
        rmatrix_[2][1] = multi_rates_[4]; // T->C
        rmatrix_[1][3] = multi_rates_[5]; // G->T
        rmatrix_[3][1] = multi_rates_[5]; // T->G
        rmatrix_[0][0] = (multi_rates_[0]+multi_rates_[1]+multi_rates_[2]) * -1;
        rmatrix_[1][1] = (multi_rates_[2]+multi_rates_[4]+multi_rates_[5]) * -1;
        rmatrix_[2][2] = (multi_rates_[0]+multi_rates_[3]+multi_rates_[4]) * -1;
        rmatrix_[3][3] = (multi_rates_[1]+multi_rates_[3]+multi_rates_[5]) * -1;
        for (unsigned int i = 0; i < 6; i++) {
            multi_rates_.erase (multi_rates_.begin() + 0); 
        }
        //QMatrix = calcQmatrix(rate_matrix);
        //QMatrix = calculate_q_matrix();
    }
    QMatrix = calculate_q_matrix();
    Node * root = tree_->getRoot();
    seqs_[root] = root_sequence_;
    ancq_[root] = QMatrix;
    
    if (show_ancs_) {
        string tname = root->getName();
        Sequence seq(tname, root_sequence_);
        res.push_back(seq);
    }
    
    // Pre-Order Traverse the tree
    for (int k = (tree_->getNodeCount() - 2); k >= 0; k--) {
        brlength = tree_->getNode(k)->getBL();
        /*
        for (unsigned int i = 0; i < QMatrix.size(); i++) {
            for (unsigned int j = 0; j < QMatrix.size(); j++) {
                cout << QMatrix[i][j] << " ";
            }
            cout << "\n";
        }
        cout << "\n";
        */
        if (multi_model_) {
            check = (int)round(multi_rates_[0]);
            //cout << check << " " << rate_count << endl;
            if (tree_->getNode(k)->isInternal() == true && multi_rates_.size() != 0) {
                if (check == rate_count) {
                    rmatrix_[0][2] = multi_rates_[1];
                    rmatrix_[2][0] = multi_rates_[1];
                    rmatrix_[0][3] = multi_rates_[2];
                    rmatrix_[3][0] = multi_rates_[2];
                    rmatrix_[0][1] = multi_rates_[3];
                    rmatrix_[1][0] = multi_rates_[3];
                    rmatrix_[2][3] = multi_rates_[4];
                    rmatrix_[3][2] = multi_rates_[4];
                    rmatrix_[1][2] = multi_rates_[5];
                    rmatrix_[2][1] = multi_rates_[5];
                    rmatrix_[1][3] = multi_rates_[6];
                    rmatrix_[3][1] = multi_rates_[6];
                    rmatrix_[0][0] = (multi_rates_[1]+multi_rates_[2]+multi_rates_[3]) * -1;
                    rmatrix_[1][1] = (multi_rates_[3]+multi_rates_[5]+multi_rates_[6]) * -1;
                    rmatrix_[2][2] = (multi_rates_[1]+multi_rates_[4]+multi_rates_[5]) * -1;
                    rmatrix_[3][3] = (multi_rates_[2]+multi_rates_[4]+multi_rates_[6]) * -1;

                    for (unsigned int i = 0; i < 7; i++) {
                        multi_rates_.erase(multi_rates_.begin() + 0);
                    }
                    //cout << "Size " << multi_rates_.size() << endl;
                    //cout << multi_rates_[0] << endl;
                    /*
                    for (unsigned int i = 0; i < multi_rates_.size(); i++) {
                        cout << multi_rates_[i] << endl;
                    }*/        
                }
                rate_count++;
            }
        }
        QMatrix = calculate_q_matrix();
        //PMatrix = calculate_p_matrix(QMatrix, brlength);
        Node * dec = tree_->getNode(k);
        Node * parent = tree_->getNode(k)->getParent();
        //ancq[dec] = QMatrix;
        vector < vector <double> > Qparent = ancq_[parent];
        string ancSeq = seqs_[parent];
        string decSeq = simulate_sequence(ancSeq, Qparent, brlength);
        /*
        for (unsigned int i = 0; i < Qparent.size(); i++) {
            for (unsigned int j = 0; j < Qparent.size(); j++) {
                cout << Qparent[i][j] << " ";
            }
            cout << "\n";
        }
        cout << "\n";*/
        
        seqs_[dec] = decSeq;
        ancq_[dec] = QMatrix; // why store this?
        if (show_ancs_ && tree_->getNode(k)->isInternal() == true) {
            string tname = tree_->getNode(k)->getName();
            Sequence seq(tname, decSeq);
            res.push_back(seq);
        }
        // If its a tip print the name and the sequence
        if (tree_->getNode(k)->isInternal() != true) {
            string tname = tree_->getNode(k)->getName();
            Sequence seq(tname, decSeq);
            res.push_back(seq);
        }
    }
}


// this should probably be returned on its own
void SequenceGenerator::print_node_labels() {
    cout << getNewickString(tree_) << endl;
    //cout << tree_->getRoot()->getNewick(true) <<";" << endl;
}


void SequenceGenerator::label_internal_nodes() {

    int count = 1;
    string str = "Node";
    string nlabel = "";
    Node * root = tree_->getRoot();
    root->setName("Node_0");
    for (int k = (tree_->getNodeCount() - 2); k >= 0; k--) {
        if (tree_->getNode(k)->isInternal() == true) {
            //cout << k << endl;
            str = to_string(count);
            nlabel = "Node_" + str;
            tree_->getNode(k)->setName(nlabel);
            count++;
        }
    }
}


// involves both gamma and pinvar
vector <float> SequenceGenerator::set_site_rates () {
    vector <float> srates(seqlen_, 1.0);
    
    // invariable sites
    if (pinvar_ != 0.0) {
        int numsample = seqlen_ * pinvar_ + 0.5;
        // sample invariable sites
        vector <int> randsites = sample_without_replacement(seqlen_, numsample);
        // must be a more elegant way of doing this
        for (int i = 0; i < numsample; i++) {
            srates[randsites[i]] = 0.0;
        }
    }
    
    // gamma-distributed rate variation. could explore other distributions...
    if (alpha_ != -1.0) { // default i.e. no rate variation
        for (int i = 0; i < seqlen_; i++) {
            // want to skip over sites that are set to invariable
            if (srates[i] != 0.0) {
                srates[i] = get_gamma_random_deviate(alpha_);
            }
        }
    }
    return srates;
} 


// initialized as string of length seqlength, all 'G'
string SequenceGenerator::generate_random_sequence () {
    
    string ancseq(seqlen_, 'G');
    if (is_dna_) {
        //string ancseq(seqlen, 'G');
        vector <double> cumsum(4);
        std::vector <double>::iterator low;
        // cumulative sum
        std::partial_sum(base_freqs_.begin(), base_freqs_.end(), cumsum.begin(), plus<double>());
    
        for (int i = 0; i < seqlen_; i++) {
            float RandNumb = get_uniform_random_deviate();
            low = std::lower_bound (cumsum.begin(), cumsum.end(), RandNumb);
            ancseq[i] = nucleotides_[low - cumsum.begin()];
        }
    } else {
        //string ancseq(seqlen, 'G');
        vector <double> cumsum(20);
        std::vector <double>::iterator low;
        // cumulative sum
        std::partial_sum(aa_freqs_.begin(), aa_freqs_.end(), cumsum.begin(), plus<double>());
    
        for (int i = 0; i < seqlen_; i++) {
            float RandNumb = get_uniform_random_deviate();
            low = std::lower_bound (cumsum.begin(), cumsum.end(), RandNumb);
            ancseq[i] = amino_acids_[low - cumsum.begin()];
        }    
    }
    //cout << ancseq << endl;
    return ancseq;
}


// rates are in order: A<->C,A<->G,A<->T,C<->G,C<->T,G<->T
vector < vector <double> > SequenceGenerator::construct_rate_matrix (vector <double> const& rates) {
    
    // initialize
    vector < vector <double> > ratemat(nstates_, vector<double>(4, 0.33));
    
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


// make sure sequence contains only valid nucleotide characters
void SequenceGenerator::check_valid_sequence () {
    // make sure uppercase
    root_sequence_ = string_to_upper(root_sequence_);
    if (is_dna_) {
        std::size_t found = root_sequence_.find_first_not_of(nucleotides_);
        if (found != std::string::npos) {
            cout << "Error: illegal character '" << root_sequence_[found] << "' at position " 
                << found+1 << " (only A,C,G,T allowed). Maybe specify AA with -c? Exiting." << endl;
            exit(0);
        }
    } else {
        std::size_t found = root_sequence_.find_first_not_of(amino_acids_);
        if (found != std::string::npos) {
            cout << "Error: illegal character '" << root_sequence_[found] << "' at position " 
                << found+1 << " (only AA chars allowed). Exiting." << endl;
            exit(0);
        }        
        
    }
}


// not sure of a more elegant way to do this...
vector<Sequence> SequenceGenerator::get_sequences () {
    return res;
}


// call this whenever a random float is needed
float SequenceGenerator::get_uniform_random_deviate () {
    return uniformDistrib_(generator_);
}


float SequenceGenerator::get_gamma_random_deviate (float alpha) {
    return gammaDistrib_(generator_);
}

//SEQGEN::~SEQGEN() {
//    // TODO Auto-generated destructor stub
//}
