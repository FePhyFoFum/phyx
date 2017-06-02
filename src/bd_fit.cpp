#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <nlopt.hpp>

using namespace std;

#include "bd_fit.h"
#include "tree.h"
#include "node.h"
#include "tree_utils.h"

typedef struct {
    int N;
    vector <double> bt;
} analysis_data;

int counter = 0;

BDFit::BDFit (Tree * intree, string const& modelflavour) {
    tree_ = intree;
    model_ = modelflavour;
    fit_model();
}


void BDFit::fit_model () {
    treelength_ = get_tree_length(tree_);
    nintnodes_ = tree_->getInternalNodeCount();
    nspeciation_ = nintnodes_ - 1.0;
    ntips_ = tree_->getExternalNodeCount();
    rootheight_ = tree_->getRoot()->getHeight();
    if (model_ == "yule") {
        fit_yule();
    } else if (model_ == "bd") {
        fit_bd();
    } else if (model_ == "best") {
        get_best_model();
    } else {
        cout << "Huh?!?" << endl;
    }
}

// find the best model between yule and bd using aic
void BDFit::get_best_model () {
    treelength_ = get_tree_length(tree_);
    nintnodes_ = tree_->getInternalNodeCount();
    nspeciation_ = nintnodes_ - 1.0;
    ntips_ = tree_->getExternalNodeCount();
    rootheight_ = tree_->getRoot()->getHeight();
    fit_yule();
    fit_bd();
    if (aic_yule_ < aic_) {
        cout << "yule model fits better by " << aic_ - aic_yule_
                << " AIC units" << endl;
        model_ = "yule";
    } else {
        cout << "bd model fits better by " << aic_yule_ - aic_
                << " AIC units" << endl;
        model_ = "bd";
    }
    
}

/*
ML for lambda is simply: (nintnodes - 1) / treelength
  - start with 2 lineages i.e. do not count speciation at root
  - that is, there are (nintnodes - 1) speciation events

Likelihood = (nintnodes - 1) * log(lambda) - lambda * treelength + lfactorial(nintnodes)

for n, x = (n+1); lfactorial(x) = (x - 0.5)*log(x) - x + 0.5*log(2*PI) + 1.0/(12.0*x);

Easier: factorial(n) = gamma(n+1)
        lfactorial(n) = std::lgamma(n+1)
*/
void BDFit::fit_yule () {
    lambda_yule_ = nspeciation_ / treelength_;
    likelihood_yule_ = nspeciation_ * log(lambda_yule_) - lambda_yule_ * treelength_
        + std::lgamma(nintnodes_ + 1.0);
    get_aic (likelihood_yule_);
    aicc_yule_ = aicc_;
    aic_yule_ = aic_;
}


void BDFit::fit_bd () {
    branching_times_.resize(nintnodes_);
    for (int i=0; i < nintnodes_; i++) {
        branching_times_[i] = tree_->getInternalNode(i)->getHeight();
    }
    
    // sort in descending order
    sort(branching_times_.begin(), branching_times_.end(), std::greater<double>());
    
    analysis_data a;
    a.N = ntips_;
    a.bt = branching_times_;
    
    // explore different algorithms here
    // the following give error: nlopt roundoff-limited
    //nlopt::opt opt(nlopt::LN_COBYLA, 2);
    //nlopt::opt opt(nlopt::LN_BOBYQA,2);
    
    nlopt::opt opt(nlopt::LN_PRAXIS, 2); // example: variable iterations, most smaller
    //nlopt::opt opt(nlopt::LN_NELDERMEAD, 2); // example: 223 iterations
    //nlopt::opt opt(nlopt::LN_SBPLX, 2); // example: 316 iterations
    
    opt.set_min_objective(nlopt_bd_log_lik, &a);
    
    // parameters: r, epsilon
    std::vector<double> x(2);
    // starting values. maybe use intelligent starting values
    x[0] = 0.05; x[1] = 0.5; // lambda = 0.1, mu = 0.05
    
    // ML vals for example tree
    //x[0] = 0.7383142; x[1] = 0.3018887;
    
    // lower bounds: 0 for both r, epsilon
    opt.set_lower_bounds(0.0);
    // upper bounds: none for r, 1 for epsilon
    std::vector<double> ub(2);
    ub[0] = 50; ub[1] = 1;
    opt.set_upper_bounds(ub);
    
    double minf;
    nlopt::result result = opt.optimize(x, minf);
    
    r_ = x[0];
    epsilon_ = x[1];
    lambda_ = r_ / (1 - epsilon_);
    mu_ = lambda_ - r_;
    
    likelihood_ = -minf;
    get_aic (likelihood_);
    //cout << "found minimum after " << counter << " evaluations." << endl;
}


double nlopt_bd_log_lik (const std::vector<double> &x, std::vector<double> &grad,
    void *data) {
    // count iterations for optimization of algorithm
    counter++;
    
    analysis_data * d = (analysis_data *) data;
    
    int N = d->N;
    vector <double> bt = d->bt;
    
    double lik = std::lgamma(N) + (N - 2) * log(x[0]) + N * log(1 - x[1]);
    
    lik += (x[0] * accumulate(bt.begin()+1, bt.end(), 0.0));
    
    double tempsum = 0.0;
    for (unsigned int i = 0; i < bt.size(); i++) {
        tempsum += log(exp(bt[i] * x[0]) - x[1]);
    }
    
    lik += tempsum * (-2);
    
    return -lik;
}


// calculate model-specific raw and small-sample-corrected AIC
// 'n' here (number of data points) is taken as the number of terminals
void BDFit::get_aic (const double &lik) {
    double K = 1.0;
    double n = ntips_;
    if (model_ == "bd") {
        K = 2.0;
    }
    aic_ = (-2.0 * lik) + (2.0 * K);
    aicc_ = aic_ + (2 * K * (K + 1)) / (n - K - 1);
}


// probably want in terms of r/epsilon too
// probably want to support model comparison, too
void BDFit::get_pars (ostream* poos) {
    (*poos) << "ntips: " << ntips_ << endl;
    (*poos) << "nspeciation: " << nspeciation_ << endl;
    (*poos) << "treelength: " << treelength_ << endl;
    (*poos) << "rootheight: " << rootheight_ << endl;
    (*poos) << "model: " << model_ << endl;
    if (model_ == "yule") {
        (*poos) << "likelihood: " << likelihood_yule_ << endl;
        (*poos) << "aic: " << aic_yule_ << endl;
        (*poos) << "aicc: " << aicc_yule_ << endl;
        (*poos) << "b: " << lambda_yule_ << endl;
    } else if (model_ == "bd") {
        (*poos) << "likelihood: " << likelihood_ << endl;
        (*poos) << "aic: " << aic_ << endl;
        (*poos) << "aicc: " << aicc_ << endl;
        (*poos) << "b: " << lambda_ << endl;
        (*poos) << "d: " << mu_ << endl;
        (*poos) << "r (b-d): " << r_ << endl;
        (*poos) << "e (d/b): " << epsilon_ << endl;
    }
    
}
