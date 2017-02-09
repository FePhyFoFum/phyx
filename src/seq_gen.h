/*
 * seqgen.h
 *
 *  Created on: Jun 23, 2015
 *      Author: joe
 */

#ifndef _SEQ_GEN_H_
#define _SEQ_GEN_H_

#include <random>

using namespace std;

#include "sequence.h"
#include "tree.h"

class SequenceGenerator {

private:
    // constant values
    Tree * tree_;
    
    int seqlen_;
    int nreps_;
    int seed_;
    int nstates_; // number of character states
    
    float alpha_;
    float pinvar_;
    
    string root_sequence_;
    
    vector <double> base_freqs_;
    vector <double> aa_freqs_;
    vector < vector <double> > rmatrix_;
    vector <double> multi_rates_;
    vector <float> site_rates_;
    
    bool show_ancs_;
    bool print_node_labels_;
    bool multi_model_;
    bool is_dna_;
    
    // hard-coded stuff
    static map<char, int> nuc_map_;
    static map<char, int> aa_map_;
    static string nucleotides_;
    static string amino_acids_;
    
    mt19937 generator_;
    std::uniform_real_distribution<float> uniformDistrib_;
    std::gamma_distribution<float> gammaDistrib_;
    
    
    // set all values
    void initialize ();
    
    // intermediate results
    map <Node *, string> seqs_;
    map <Node *, vector < vector <double> > > ancq_;
    
    // the result to return
    vector <Sequence> res;
    
    // los funciones
    void print_node_labels ();
    void label_internal_nodes();
    void preorder_tree_traversal ();
    //vector <float> site_rates;
    vector < vector <double> > calculate_q_matrix ();
    vector < vector<double> > calcQmatrix (vector<vector <double>>);
    vector < vector <double> > calculate_p_matrix (vector < vector <double> > const&QMatrix,
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
        vector <double> const& aabasefreq, bool const& is_dna);
    
    // return results
    vector<Sequence> get_sequences ();
    
};

#endif /* _SEQ_GEN_H_ */
