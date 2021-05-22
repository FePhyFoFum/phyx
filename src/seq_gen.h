#ifndef PX_SEQ_GEN_H
#define PX_SEQ_GEN_H

#include <string>
#include <vector>
#include <map>
#include <random>

#include "sequence.h"

class Tree; // forward declaration
class Node; // forward declaration

class SequenceGenerator {

private:
    // constant values
    Tree * tree_;
    
    int seqlen_;
    int nreps_;
    long int seed_;
    int nstates_; // number of character states
    
    float alpha_;
    float pinvar_;
    
    std::string root_sequence_;
    
    std::vector<double> base_freqs_;
    std::vector<double> aa_freqs_;
    std::vector< std::vector<double> > rmatrix_;
    std::vector<double> multi_rates_;
    std::vector<float> site_rates_;
    
    bool show_ancs_;
    bool print_node_labels_;
    bool multi_model_;
    bool is_dna_;
    
    // hard-coded stuff
    static std::map<char, int> nuc_map_;
    static std::map<char, int> aa_map_;
    static std::string nucleotides_;
    static std::string amino_acids_;
    
    std::mt19937 generator_;
    std::uniform_real_distribution<float> uniformDistrib_;
    std::gamma_distribution<float> gammaDistrib_;
    
    // set all values
    void initialize();
    
    // intermediate results
    std::map<Node *, std::string> seqs_;
    std::map<Node *, std::vector< std::vector<double> > > ancq_;
    
    // the result to return
    std::vector<Sequence> res;
    
    // los funciones
    void print_node_labels ();
    void label_internal_nodes ();
    void preorder_tree_traversal ();
    //std::vector<float> site_rates;
    std::vector< std::vector<double> > calculate_q_matrix ();
    std::vector< std::vector<double> > calcQmatrix (std::vector< std::vector<double> >);
    std::vector< std::vector<double> > calculate_p_matrix (const std::vector< std::vector<double> >& QMatrix,
        float br);
    std::string simulate_sequence (const std::string& anc, std::vector< std::vector<double> >& Matrix,
        const float& brlength);
    std::string generate_random_sequence ();
    std::vector< std::vector<double> > construct_rate_matrix (const std::vector<double>& rates);
    void check_valid_sequence ();
    float get_uniform_random_deviate ();
    float get_gamma_random_deviate (float);
    std::vector<float> set_site_rates ();
    
public:
    SequenceGenerator (const int&seqlength, const std::vector<double>& basefreq,
        std::vector< std::vector<double> >& rmatrix, Tree * tree, const bool& showancs, 
        const int& nreps, const long int& seed, const float& alpha, const float& pinvar,
        const std::string& ancseq, const bool& printpost, const std::vector<double>& multirates,
        const std::vector<double>& aabasefreq, const bool& is_dna);
    
    // return results
    std::vector<Sequence> get_sequences ();
};

#endif /* PX_SEQ_GEN_H */
