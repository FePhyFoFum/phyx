#ifndef _MCMC_H_
#define _MCMC_H_

#include "state_reconstructor_simple.h"
#include "rate_model.h"
#include "tree.h"
#include "sequence.h"

#include <armadillo>
using namespace arma;

void sm0_mcmc(int reps, int sampleiter, Tree * tree, StateReconstructorSimple & sr,
    RateModel & rm, vector<Sequence> & seqs, vector<Sequence> & sr_seqs,
    map<string,vector<int> > & codon_pos, mat &bf, mat &K, mat &w, mat &inq);

void sm2a_mcmc(int reps, int sampleiter, Tree * tree, StateReconstructorSimple & sr,
    RateModel & rm, vector<Sequence> & seqs, vector<Sequence> & sr_seqs,
    map<string,vector<int> > & codon_pos, mat &bf, mat &K, mat &w, mat &inq0,
    mat &inq1, mat & inq2);

void sm0_prof(int reps, int sampleiter, Tree * tree, StateReconstructorSimple & sr,
    RateModel & rm, vector<Sequence> & seqs, vector<Sequence> & sr_seqs, map<string,
    vector<int> > & codon_pos, mat &bf, mat &K, mat &w, mat &inq);

void sm0_bmcmc(int reps, int sampleiter, Tree * tree, StateReconstructorSimple & sr,
    RateModel & rm, vector<Sequence> & seqs, vector<Sequence> & sr_seqs,
    map<string, vector<int> > & codon_pos, mat &bf, mat &K, mat &w, mat &inq);

double accept_prob(double & curLike, double & propLike, double & curPrior, double & propPrior,
    int & curD, int & propD, int & N);

double likelihood_ratio(double & curLike, double & propLike);

double hastings_ratio(int & curD, int & propD, int & N);

#endif /* _MCMC_H_ */
