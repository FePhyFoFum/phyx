#include <string>
#include <vector>
#include <iostream>
//#include <sstream>
#include <cmath>
//#include <limits>
//#include <cstdlib>
#include <algorithm>
#include <nlopt.hpp>

using namespace std;

#include "bd_fit.h"
#include "tree.h"
#include "node.h"
#include "tree_utils.h"
//#include "constants.h" // for PI. not needed anymore  since lgamma used instead

typedef struct {
    int N;
    vector <double> bt;
} analysis_data;

int counter = 0;

BDFit::BDFit (Tree * intree, string const& modelflavour) {
    tree = intree;
    model = modelflavour;
    fit_model();
}


void BDFit::fit_model () {
    treelength = get_tree_length(tree);
    nintnodes = tree->getInternalNodeCount();
    nspeciation = nintnodes - 1.0;
    ntips = tree->getExternalNodeCount();
    
    if (model == "yule") {
        fit_yule();
    } else if (model == "bd") {
        fit_bd();
    } else {
        cout << "Huh?!?" << endl;
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
    lambda = nspeciation / treelength;
    likelihood = nspeciation * log(lambda) - lambda * treelength
        + std::lgamma(nintnodes + 1.0);
    /*
    cout << "ntips: " << ntips << endl;
    cout << "nintnodes: " << nintnodes << endl;
    cout << "nspeciation: " << nspeciation << endl;
    cout << "lambda: " << lambda << endl;
    cout << "likelihood: " << likelihood << endl;
    */
}



void BDFit::fit_bd () {
    branching_times.resize(nintnodes);
    for (int i=0; i < nintnodes; i++) {
        branching_times[i] = tree->getInternalNode(i)->getHeight();
    }
    
    // sort in descending order
    sort(branching_times.begin(), branching_times.end(), std::greater<double>());
    
    analysis_data a;
    a.N = ntips;
    a.bt = branching_times;
    
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
    
    r = x[0];
    epsilon = x[1];
    lambda = r / (1 - epsilon);
    mu = lambda - r;
    
    likelihood = -minf;
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
// 'n' here (numbr of data points) is taken as the number of terminals
void BDFit::get_aic () {
    double K = 1.0;
    double n = ntips;
    if (model == "bd") {
        K = 2.0;
    }
    aic = (-2.0 * likelihood) + (2.0 * K);
    aicc = aic + (2 * K * (K + 1)) / (n - K - 1);
}


// probably want in terms of r/epsilon too
// probably want to support model comparison, too
void BDFit::get_pars (ostream* poos) {
    (*poos) << "ntips: " << ntips << endl;
    (*poos) << "nspeciation: " << nspeciation << endl;
    (*poos) << "treelength: " << treelength << endl;
    (*poos) << "model: " << model << endl;
    (*poos) << "likelihood: " << likelihood << endl;
    get_aic ();
    (*poos) << "aic: " << aic << endl;
    (*poos) << "aicc: " << aicc << endl;
    if (model == "yule") {
        (*poos) << "lambda: " << lambda << endl;
    } else {
        (*poos) << "lambda: " << lambda << endl;
        (*poos) << "mu: " << mu << endl;
        (*poos) << "r: " << r << endl;
        (*poos) << "epsilon: " << epsilon << endl;
    }
    
}
