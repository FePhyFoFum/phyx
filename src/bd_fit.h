#ifndef _BD_FIT_H_
#define _BD_FIT_H_

#include <string>

using namespace std;

#include "tree.h"

class BDFit {
private:
    string model;
    Tree* tree;
    double treelength;
    double nintnodes;
    double nspeciation;
    double ntips;
    vector <double> branching_times;
    double rootheight;
    double lambda;
    double mu;
    double r;
    double epsilon;
    double likelihood;
    double aicc;
    double aic;
    
    void fit_model ();
    void fit_yule ();
    void fit_bd ();
    void get_aic ();

public:
    BDFit (Tree * intree, string const& modelflavour);
    void get_pars (ostream* poos);
};

// non-member function, as nlopt is weird with pointers...
double nlopt_bd_log_lik (const std::vector<double> &x, std::vector<double> &grad,
    void *data);

#endif /* _BD_FIT_H_ */
