#include <string>
#include <vector>
#include <iostream>
//#include <sstream>
#include <cmath>
//#include <limits>
//#include <cstdlib>
//#include <algorithm>

using namespace std;

#include "bd_fit.h"
#include "tree.h"
#include "node.h"
#include "tree_utils.h"
#include "constants.h" // for PI


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


// coming soon
void BDFit::fit_bd () {
    
}


void BDFit::get_pars (ostream* poos) {
    (*poos) << "ntips: " << ntips << endl;
    (*poos) << "nspeciation: " << nspeciation << endl;
    (*poos) << "treelength: " << treelength << endl;
    (*poos) << "lambda: " << lambda << endl;
    if (model == "bd") {
        (*poos) << "epsilon: " << epsilon << endl;
    }
    (*poos) << "likelihood: " << likelihood << endl;
}
