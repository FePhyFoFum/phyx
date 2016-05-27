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
    double lambda;
    double mu;
    double r;
    double epsilon;
    double likelihood;
    
    void fit_model ();
    void fit_yule ();
    void fit_bd ();
    double nlopt_bd_log_lik (const std::vector<double> &x, std::vector<double> &grad,
        void *data);
    double foo (const std::vector<double> &x);

public:
    BDFit(Tree * intree, string const& modelflavour);
    void get_pars (ostream* poos);
};

#endif /* _BD_FIT_H_ */
