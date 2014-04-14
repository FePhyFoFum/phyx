#include <vector>
#include <map>

#include "state_reconstructor_simple.h"
#include "rate_model.h"
#include "mcmc.h"
#include "tree.h"
#include "sequence.h"
#include "seq_utils.h"


void sm0_mcmc(int reps, int sampleiter, Tree * tree, StateReconstructorSimple & sr, RateModel & rm, vector<Sequence> & seqs, vector<Sequence> & sr_seqs, map<string,vector<int> > & codon_pos){
    int sites = (seqs[0].get_sequence().size()/3);
    double curlike = 0;
    for(int s=0;s<sites;s++){
	create_vector_seq_codon_state_reconstructor(seqs,sr_seqs,s,codon_pos);
	sr.set_tip_conditionals(sr_seqs);
	curlike +=  sr.eval_likelihood();
	rm.set_sameQ(true);
    }
    cout << "start likelihood: " << curlike << endl;
    double newlike = 0;
    double K = 1.0; //start
    double w = 1.0; //startw
    for (int i=0;i<reps;i++){
       	for(int s=0;s<sites;s++){
	    create_vector_seq_codon_state_reconstructor(seqs,sr_seqs,s,codon_pos);
	    sr.set_tip_conditionals(sr_seqs);
	    newlike +=  sr.eval_likelihood();
	    rm.set_sameQ(true);
	}
	if(i%sampleiter == 0){
	    cout << "iter: "<< i << " like: " << curlike << endl;
	}
    }
}
