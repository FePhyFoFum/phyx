#ifndef PX_RATE_MODEL_H
#define PX_RATE_MODEL_H

#include <armadillo>
#include <string>
#include <vector>
#include <map>

using namespace arma;

class RateModel {
private:
    mat Q;
    std::vector<mat> Qs;
    std::vector<std::string> labels;
    std::vector< std::vector<double> > Q_mask;
    bool sameQ;
    cx_mat lasteigval;
    cx_mat lasteigvec;
    cx_mat eigval;
    cx_mat eigvec;
    mat lasteigval_simple;
    mat lasteigvec_simple;
    mat eigval_simple;
    mat eigvec_simple;
    bool lastImag;
    int fortran_iexph;
    double * fortran_wsp;
    int fortran_m;
    
public:
    explicit RateModel(int _nstates);
    /*
      storing once optimization has occurred
      map of bl and p matrix
      map<branch length,p matrix>
    */
    int selection_model; //0,2 = 2a
    int nstates;
    std::map<double, cx_mat> stored_p_matrices;
    bool neg_p;
    void set_n_qs (int number);
    void set_Q_which (mat& inQ, int which);
    void set_Q (mat& inQ);
    void set_Q_diag ();
    void set_Q_cell (int, int, double);
    void setup_Q ();
    void setup_Q (std::vector< std::vector<double> >& inQ);
    void setup_Q (mat& inQ);
    mat& get_Q ();
    void set_sameQ (bool);
    cx_mat setup_P (double, bool);
    void setup_P_simple (mat& p, double, bool);
    void setup_fortran_P_whichQ (int which, mat& P, double t);
    void setup_fortran_P (mat& P, double t, bool store_p_matrices);
    /*
     * get things from stmap
     */
    //this should be used for getting the eigenvectors and eigenvalues
    bool get_eigenvec_eigenval_from_Q (cx_mat * eigval, cx_mat * eigvec);
    void get_eigenvec_eigenval_from_Q_simple (mat * eigval, mat * eigvec);
    
};
void update_simple_goldman_yang_q (mat * inm, double K, double w, mat& bigpibf,
        mat& bigpiK, mat& bigpiw);
bool test_transition (char a, char b);
void generate_bigpibf_K_w (mat * bf, mat * K, mat * w, std::map<std::string,
        std::string>& codon_dict, std::map<std::string, std::vector<int> >& codon_index,
        std::vector<std::string>& codon_list);
void convert_matrix_to_single_row_for_fortran (mat& inmatrix, double t, double * H);

#endif /* PX_RATE_MODEL_H */
