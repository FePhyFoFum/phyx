#include <vector>
#include <map>
#include <random>

#include "state_reconstructor_simple.h"
#include "rate_model.h"
#include "mcmc.h"
#include "tree.h"
#include "sequence.h"
#include "seq_utils.h"

#include <armadillo>
using namespace arma;

void sm0_mcmc(int reps, int sampleiter, Tree * tree, StateReconstructorSimple & sr,
    RateModel & rm, vector<Sequence> & seqs, vector<Sequence> & sr_seqs,
    map<string,vector<int> > & codon_pos, mat &bf, mat &K, mat &w, mat &inq) {
    
    int sites = (seqs[0].get_sequence().size()/3);
    double curlike = 0;
    double sw = 0.5;
    rm.selection_model = 0;
    for (int s=0; s < sites; s++) {
        create_vector_seq_codon_state_reconstructor(seqs, sr_seqs, s, codon_pos);
        sr.set_tip_conditionals(sr_seqs, s);
        curlike +=  sr.eval_likelihood(s);
        rm.set_sameQ(true);
    }
    cout << "start likelihood: " << curlike << endl;

    std::default_random_engine re;
    
    double newlike = 0;
    double startk = 1.0; //start
    double startw = 1.0; //startw
    for (int i=0; i < reps; i++) {
        std::uniform_real_distribution<double> unif1(startk - sw, startk + sw);
        std::uniform_real_distribution<double> unif2(startw - sw, startw + sw);
        double newk = fabs(unif1(re));
        double neww = fabs(unif2(re));
        newlike = 0;
        update_simple_goldman_yang_q(&inq, newk, neww, bf, K, w);
        rm.setup_Q(inq);
        sr.clear_map_ps();
        for (int s=0; s < sites; s++) {
            //create_vector_seq_codon_state_reconstructor(seqs,sr_seqs,s,codon_pos);
            //sr.set_tip_conditionals(sr_seqs,s);
            newlike +=  sr.eval_likelihood(s);
            rm.set_sameQ(true);
        }
        if (newlike < curlike) {
            curlike = newlike;
            startk = newk;
            startw = neww;
        }
        if (i%sampleiter == 0) {
            cout << "iter: "<< i << " like: " << curlike << " K: "<< startk <<" w: " << startw << endl;
        }
    }
}


/**
 * this should have 5 parameters
 * K = free, w0 = 0,w1 = 1,w2 = free ,p0 = free ,p1 = free ,p2 = free
 */
void sm2a_mcmc(int reps, int sampleiter, Tree * tree, StateReconstructorSimple & sr,
    RateModel & rm, vector<Sequence> & seqs, vector<Sequence> & sr_seqs,
    map<string,vector<int> > & codon_pos, mat &bf, mat &K, mat &w, mat & inq0,
    mat & inq1, mat &inq2) {
    
    int sites = (seqs[0].get_sequence().size()/3);
    double curlike = 0;
    double sw = 0.2;
    sr.pp0 = 0.38008; sr.pp1 = 0.28326; sr.pp2 = 0.33666;

    for (int s=0; s < sites; s++) {
        create_vector_seq_codon_state_reconstructor(seqs, sr_seqs, s, codon_pos);
        sr.set_tip_conditionals(sr_seqs, s);
        curlike +=  sr.eval_likelihood(s);
        rm.set_sameQ(true);
    }
    cout << "start likelihood: " << curlike << endl;

    std::default_random_engine re;
    
    double newlike = 0;
    double startk = 1.36714; //start
    double startw0 = 0.0; //start w0
    double startw2 = 8.0; //start w2
    double startp0 = 0.38008;
    double startp1 = 0.28326;
    double startp2 = 0.33666; //not free
    for (int i=0; i < reps; i++) {
        std::uniform_real_distribution<double> unif1(startk - sw, startk + sw);
        std::uniform_real_distribution<double> unif2(startw0 - sw, startw0 + sw);
        std::uniform_real_distribution<double> unif3(startw2 - sw, startw2 + sw);
        std::uniform_real_distribution<double> unif4(startp0 - sw, startp0 + sw);
        std::uniform_real_distribution<double> unif5(startp1 - sw, startp1 + sw);
        std::uniform_real_distribution<double> unif6(startp2 - sw, startp2 + sw);

        double newk = fabs(unif1(re));
        double neww0 = fabs(unif2(re));
        if (neww0 >= 1) {
            neww0 = 1 - (neww0-1);
        }
        double neww2 = fabs(unif3(re));
        double newp0 = fabs(unif4(re));
        double newp1 = fabs(unif5(re));
        double newp2 = fabs(unif6(re));

        double sum1 = newp0 + newp1 + newp2;
        newp0 = newp0/sum1; newp1 = newp1/sum1; newp2 = newp2/sum1;
        sr.pp0 = newp0; sr.pp1 = newp1; sr.pp2 = newp2;

        newlike = 0;
        //only do it with 
        update_simple_goldman_yang_q(&inq0, newk, neww0, bf, K, w);
        update_simple_goldman_yang_q(&inq2, newk, neww2, bf, K, w);
        rm.set_Q_which(inq0, 0);
        rm.set_Q_which(inq2, 2);
        sr.clear_map_ps();
        for (int s=0; s < sites; s++) {
            //create_vector_seq_codon_state_reconstructor(seqs,sr_seqs,s,codon_pos);
            //sr.set_tip_conditionals(sr_seqs,s);
            newlike +=  sr.eval_likelihood(s);
            rm.set_sameQ(true);
        }
        if (newlike < curlike) {
            curlike = newlike;
            startk = newk;
            startw0 = neww0;
            startw2 = neww2;
            startp0 = newp0;
            startp1 = newp1;
            startp2 = newp2;
        }
        if (i % sampleiter == 0) {
            cout << "iter: "<< i << " like: " << curlike << " K: "<< startk <<" w0: "
            << startw0 << " w2: " << startw2 << " p0: " << startp0 << " p1: " << startp1
            << " p2: "<< startp2 << endl;
        }
    }
}
