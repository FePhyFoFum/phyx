#ifndef _BD_FIT_H_
#define _BD_FIT_H_

using namespace std;

#include "tree.h"

class BDFit {
private:
    string model_;
    Tree* tree_;
    double treelength_;
    double nintnodes_;
    double nspeciation_;
    double ntips_;
    vector <double> branching_times_;
    double rootheight_;
    double lambda_yule_; // separate for model selection
    double lambda_;
    double mu_;
    double r_;
    double epsilon_;
    double likelihood_;
    double likelihood_yule_; // separate for model selection
    double aicc_;
    double aic_;
    double aicc_yule_; // separate for model selection
    double aic_yule_; // separate for model selection
    
    void fit_model ();
    void get_best_model ();
    void fit_yule ();
    void fit_bd ();
    void get_aic (const double &lik);

public:
    BDFit (Tree * intree, string const& modelflavour);
    void get_pars (ostream* poos);
};

// non-member function, as nlopt is weird with pointers...
double nlopt_bd_log_lik (const std::vector<double> &x, std::vector<double> &grad,
    void *data);

#endif /* _BD_FIT_H_ */
