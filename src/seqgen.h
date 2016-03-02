/*
 * seqgen.h
 *
 *  Created on: Jun 23, 2015
 *      Author: joe
 */

#ifndef _SEQ_GEN_H_
#define _SEQ_GEN_H_

#include <random>

#include "sequence.h"
#include "tree.h"

class SequenceGenerator {

private:
    // constant values
    int seqlen;
    bool mm;
    vector <double> basefreqs;
    vector <double> multirates;
    vector < vector <double> > rmatrix;
    Tree * tree;
    int nreps;
    int seed;
    string rootSequence;
    static map<char, int> nucMap;
    static string nucleotides;
    
    bool showancs;
    
    // set all values
    void initialize ();
    
    map <Node *, string> seqs;
    map <Node *, vector < vector <double> >> ancq;
    vector < vector <double> > Qparent;
    
    vector <Sequence> res;
    float alpha;
    float pinvar;
    
    bool printpost;
    void print_node_labels ();
    void label_internal_nodes();
    void preorder_tree_traversal (Tree * tree, bool showancs,
        vector <double> multirates, bool mm);
    vector <float> site_rates;
    vector < vector <double> > calculate_q_matrix ();
    vector <vector<double>> calcQmatrix (vector<vector <double>>);
    vector < vector <double> > calculate_p_matrix (vector < vector <double> > QMatrix,
        float br);
    string simulate_sequence (string const& anc, vector < vector <double> >& Matrix,
        float const& brlength);
    string generate_random_sequence ();
    vector < vector <double> > construct_rate_matrix (vector <double> const& rates);
    
    mt19937 generator;
    float get_uniform_random_deviate ();
    float get_gamma_random_deviate (float);
    vector <float> set_site_rates ();
    
public:
    
    SequenceGenerator (int const &seqlength, vector <double> const& basefreq,
        vector < vector<double> >& rmatrix, Tree * tree, bool const& showancs, 
        int const& nreps, int const & seed, float const& alpha, float const& pinvar,
        string const& ancseq, bool const& printpost, vector<double> const& multirates,
        bool const& mm);
    vector<Sequence> get_sequences ();
    
};

#endif /* _SEQ_GEN_H_ */
