#ifndef PX_BD_FIT_H
#define PX_BD_FIT_H

#include <vector>
#include <string>
#include <iostream>

class Tree; // forward declaration

class BDFit {
private:
    std::string model_;
    
    double lambda_bd_;
    double lambda_yule_;
    double mu_;
    double r_;
    double epsilon_;
    
    double likelihood_bd_;
    double likelihood_yule_;
    
    double aic_bd_;
    double aicc_bd_;
    double aic_yule_;
    double aicc_yule_;
    
    Tree* tree_;
    double treelength_;
    double nintnodes_;
    double nspeciation_;
    double ntips_;
    std::vector<double> branching_times_;
    double rootheight_;
    
    void fit_model();
    void get_best_model();
    void fit_yule();
    void fit_bd();
    void get_aic(const double& lik, double& aic, double& aicc);

public:
    BDFit(Tree * intree, std::string modelflavour);
    void get_pars(std::ostream* poos);
};

// non-member function, as nlopt is weird with pointers...
double nlopt_bd_log_lik (const std::vector<double>& x, std::vector<double>& grad,
    void *data);

#endif /* PX_BD_FIT_H */
