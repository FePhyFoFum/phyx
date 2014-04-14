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




void sm0_mcmc(int reps, int sampleiter, Tree * tree, StateReconstructorSimple & sr, RateModel & rm, vector<Sequence> & seqs, vector<Sequence> & sr_seqs, map<string,vector<int> > & codon_pos, mat &bf, mat &K, mat &w, mat &inq){
    int sites = (seqs[0].get_sequence().size()/3);
    double curlike = 0;
    double sw = 0.5;
    for(int s=0;s<sites;s++){
	create_vector_seq_codon_state_reconstructor(seqs,sr_seqs,s,codon_pos);
	sr.set_tip_conditionals(sr_seqs);
	curlike +=  sr.eval_likelihood();
	rm.set_sameQ(true);
    }
    cout << "start likelihood: " << curlike << endl;

    std::default_random_engine re;
    
    double newlike = 0;
    double startk = 1.0; //start
    double startw = 1.0; //startw
    for (int i=0;i<reps;i++){
	std::uniform_real_distribution<double> unif1(startk-sw,startk+sw);
	std::uniform_real_distribution<double> unif2(startw-sw,startw+sw);
	double newk = fabs(unif1(re));
	double neww = fabs(unif2(re));
	newlike = 0;
	update_simple_goldman_yang_q(&inq,newk,neww,bf,K,w);
	rm.setup_Q(inq);
       	for(int s=0;s<sites;s++){
	    create_vector_seq_codon_state_reconstructor(seqs,sr_seqs,s,codon_pos);
	    sr.set_tip_conditionals(sr_seqs);
	    newlike +=  sr.eval_likelihood();
	    rm.set_sameQ(true);
	}
	if(newlike < curlike){
	    curlike = newlike;
	    startk = newk;
	    startw = neww;
	}
	if(i%sampleiter == 0){
	    cout << "iter: "<< i << " like: " << curlike << " K: "<< startk <<" w: " << startw << endl;
	}
    }
}
