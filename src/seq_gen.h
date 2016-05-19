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
    Tree * tree;
    
    int seqlen;
    int nreps;
    int seed;
    
    float alpha;
    float pinvar;
    
    string rootSequence;
    
    vector <double> basefreqs;
    vector < vector <double> > rmatrix;
    vector <double> multirates;
    
    bool showancs;
    bool printnodelabels;
    bool mm;
    
    // hard-coded stuff
    static map<char, int> nucMap;
    static string nucleotides;
    
    mt19937 generator;
    std::uniform_real_distribution<float> uniformDistrib;
    std::gamma_distribution<float> gammaDistrib;
    
    // set all values
    void initialize ();
    
    // intermediate results
    map <Node *, string> seqs;
    map <Node *, vector < vector <double> > > ancq;
    vector < vector <double> > Qparent;
    
    // the result to return
    vector <Sequence> res;
    
    // los funciones
    void print_node_labels ();
    void label_internal_nodes();
    void preorder_tree_traversal (Tree * tree, bool showancs,
        vector <double> multirates, bool mm);
    vector <float> site_rates;
    vector < vector <double> > calculate_q_matrix ();
    vector < vector<double> > calcQmatrix (vector<vector <double>>);
    vector < vector <double> > calculate_p_matrix (vector < vector <double> > QMatrix,
        float br);
    string simulate_sequence (string const& anc, vector < vector <double> >& Matrix,
        float const& brlength);
    string generate_random_sequence ();
    vector < vector <double> > construct_rate_matrix (vector <double> const& rates);
    void check_valid_sequence ();
    float get_uniform_random_deviate ();
    float get_gamma_random_deviate (float);
    vector <float> set_site_rates ();
    
public:
    
    SequenceGenerator (int const &seqlength, vector <double> const& basefreq,
        vector < vector<double> >& rmatrix, Tree * tree, bool const& showancs, 
        int const& nreps, int const & seed, float const& alpha, float const& pinvar,
        string const& ancseq, bool const& printpost, vector<double> const& multirates,
        bool const& mm);
    
    // return results
    vector<Sequence> get_sequences ();
    
};

#endif /* _SEQ_GEN_H_ */
