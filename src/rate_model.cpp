#include <vector>
#include <string>
#include <map>

#include "rate_model.h"

#include <armadillo>
using namespace arma;

inline int signof(double d) {return d >= 0 ? 1 : -1;}
inline double roundto(double in) {return floor(in*(1000)+0.5)/(1000);}

RateModel::RateModel(int _nstates):Q(_nstates,_nstates),labels(),Q_mask(),
    lasteigval(_nstates,_nstates),lasteigvec(_nstates,_nstates),
    eigval(_nstates,_nstates),eigvec(_nstates,_nstates),
    lasteigval_simple(_nstates,_nstates),lasteigvec_simple(_nstates,_nstates),
    eigval_simple(_nstates,_nstates),eigvec_simple(_nstates,_nstates),nstates(_nstates) {
    
    setup_Q();
    sameQ = false;
    lasteigval.fill(0);
    lasteigvec.fill(0);
    eigval.fill(0);
    eigvec.fill(0);
    lasteigval_simple.fill(0);
    lasteigvec_simple.fill(0);
    eigval_simple.fill(0);
    eigvec_simple.fill(0);
}


void RateModel::set_Q_cell(int from, int to, double num) {
    Q(from,to) = num;
    sameQ = false;
}

void RateModel::set_Q_diag() {
    for (unsigned int i=0; i < Q.n_rows; i++) {
        double su = 0;
        for (unsigned int j=0; j < Q.n_cols; j++) {
            if (i != j) {
                su += Q(i, j);
            }
        }
        Q(i, i) = 0-su;
    }
    sameQ = false;
}

void RateModel::setup_Q() {
    Q.fill(0);
    for (unsigned int i=0; i < Q.n_rows; i++) {
        for (unsigned int j=0; j < Q.n_cols; j++) {
            if (i != j) {
                Q(i, j) = 1./nstates;
            } else {
                Q(i, j) = -(1./nstates * (nstates-1));
            }
        }
    }
    sameQ = false;
}

void RateModel::setup_Q(vector<vector<double> > & inQ) {
    for (unsigned int i=0; i < Q.n_rows; i++) {
        for (unsigned int j=0; j < Q.n_cols; j++) {
            Q(i, j) = inQ[i][j];
        }
    }
    for (unsigned int i=0; i < Q.n_rows; i++) {
        colvec a = (sum(Q,1));
        Q(i, i) = -(a(i)-Q(i, i));
    }
    sameQ = false;
}

void RateModel::setup_Q(mat & inQ) {
    for (unsigned int i=0; i < Q.n_rows; i++) {
        for (unsigned int j=0; j < Q.n_cols; j++) {
            Q(i, j) = inQ(i, j);
        }
    }
    set_Q_diag();
    sameQ = false;
}

void RateModel::set_n_qs(int number) {
    for (int i=0; i < number; i++) {
        mat tm(nstates,nstates);
        Qs.push_back(tm);
    }
}

void RateModel::set_Q_which(mat & inQ,int which) {
    for (unsigned int i=0; i < Qs[which].n_rows; i++) {
        for (unsigned int j=0; j < Qs[which].n_cols; j++) {
            Qs[which](i, j) = inQ(i, j);
        }
    }
    for (unsigned int i=0; i < Qs[which].n_rows; i++) {
        double su = 0;
        for (unsigned int j=0; j < Qs[which].n_cols; j++) {
            if (i != j) {
                su += Qs[which](i, j);
            }
        }
        Qs[which](i, i) = 0-su;
    }
    sameQ = false;
}

void RateModel::set_Q(mat & inQ) {
    for (unsigned int i=0; i < Q.n_rows; i++) {
        for (unsigned int j=0; j < Q.n_cols; j++) {
            Q(i, j) = inQ(i, j);
        }
    }
    sameQ = false;
}

mat & RateModel::get_Q() {
    return Q;
}

cx_mat RateModel::setup_P(double bl,bool store_p_matrices) {
    //sameQ = false;
    eigvec.fill(0);
    eigval.fill(0);
    bool isImag = get_eigenvec_eigenval_from_Q(&eigval, &eigvec); // isImag is not used?
    //cout << eigval << endl;
    //cout << eigvec << endl;
    for (int i=0; i < nstates; i++) {
        eigval(i, i) = exp(eigval(i, i) * bl);
    }
    cx_mat C_inv = inv(eigvec);
    cx_mat P = eigvec * eigval * C_inv;
    neg_p = false;
    for (unsigned int i=0; i < P.n_rows; i++) {
        for (unsigned int j=0; j < P.n_cols; j++) {
            if (real(P(i, j))<0)
            neg_p = true;
        }
    }
    if (store_p_matrices == true) {
        stored_p_matrices[bl] = P;    
    }
    return P;
}

void RateModel::setup_P_simple(mat & p,double bl,bool store_p_matrices) {
//    sameQ = false;
    eigvec_simple.fill(0);
    eigval_simple.fill(0);
    get_eigenvec_eigenval_from_Q_simple(&eigval_simple, &eigvec_simple);
    //cout << eigval << endl;
    //cout << eigvec << endl;
    for (int i=0; i < nstates; i++) {
        eigval_simple(i, i) = exp(eigval_simple(i, i) * bl);
    }
    mat C_inv = inv(eigvec_simple);
    p = eigvec_simple * eigval_simple * C_inv;
/*    if (store_p_matrices == true) {
    stored_p_matrices[bl] = P;    
    }*/
    //return P;
}

void RateModel::get_eigenvec_eigenval_from_Q_simple(mat * eigval, mat * eigvec) {
    if (sameQ == true) {
        for (unsigned int i=0; i < Q.n_rows; i++) {
            for (unsigned int j=0; j < Q.n_cols; j++) {
                (*eigval)(i, j) = lasteigval_simple(i, j);
                (*eigvec)(i, j) = lasteigvec_simple(i, j);
            }
        }
        return;
    }
    mat tQ(nstates,nstates); tQ.fill(0);
    for (unsigned int i=0; i < Q.n_rows; i++) {
        for (unsigned int j=0; j < Q.n_cols; j++) {
            tQ(i, j) = Q(i, j);
        }
    }
    colvec eigva;
    mat eigve;
    eig_sym(eigva,eigve,tQ);
    //bool isImag = false; // not used
    for (unsigned int i=0; i < Q.n_rows; i++) {
        for (unsigned int j=0; j < Q.n_cols; j++) {
            if (i == j) {
                (*eigval)(i, j) = eigva(i);
                lasteigval_simple(i, j) = eigva(i);
            } else {
                (*eigval)(i, j) = 0;
                lasteigval_simple(i, j) = 0;
            }
            (*eigvec)(i, j) = eigve(i, j);
            lasteigvec_simple(i, j) = eigve(i, j);
        }
    }
    return;
}


/*
 * this should be used to caluculate the eigenvalues and eigenvectors
 * as U * Q * U-1 -- eigen decomposition
 *
 * this should use the armadillo library
 */
bool RateModel::get_eigenvec_eigenval_from_Q(cx_mat * eigval, cx_mat * eigvec) {
    if (sameQ == true) {
        for (unsigned int i=0; i < Q.n_rows; i++) {
            for (unsigned int j=0; j < Q.n_cols; j++) {
                (*eigval)(i, j) = lasteigval(i, j);
                (*eigvec)(i, j) = lasteigvec(i, j);
            }
        }
        return lastImag;
    }
    mat tQ(nstates,nstates); tQ.fill(0);
    for (unsigned int i=0; i < Q.n_rows; i++) {
        for (unsigned int j=0; j < Q.n_cols; j++) {
            tQ(i, j) = Q(i, j);
        }
    }
    cx_colvec eigva;
    cx_mat eigve;
    eig_gen(eigva,eigve,tQ);
    bool isImag = false;
    for (unsigned int i=0; i < Q.n_rows; i++) {
        for (unsigned int j=0; j < Q.n_cols; j++) {
            if (i == j) {
                (*eigval)(i, j) = eigva(i);
                lasteigval(i, j) = eigva(i);
            } else {
                (*eigval)(i, j) = 0;
                lasteigval(i, j) = 0;
            }
            (*eigvec)(i, j) = eigve(i, j);
            lasteigvec(i, j) = eigve(i, j);
            if (imag((*eigvec)(i, j)) > 0 || imag((*eigval)(i, j))) {
                isImag = true;
            }
        }
    }
    //lasteigval = eigval;
    //lasteigvec = eigvec;
    lastImag = isImag;
    return isImag;
}

//taking out fortran
//
/*
extern"C" {
    void wrapalldmexpv_(int * n,int* m,double * t,double* v,double * w,double* tol,double* anorm,double* wsp,int * lwsp,int* iwsp,int *liwsp, int * itrace,int *iflag,int *ia, int *ja, double *a, int *nz, double * res);
    void wrapsingledmexpv_(int * n,int* m,double * t,double* v,double * w,double* tol,double* anorm,double* wsp,int * lwsp,int* iwsp,int *liwsp, int * itrace,int *iflag,int *ia, int *ja, double *a, int *nz, double * res);
    void wrapdgpadm_(int * ideg,int * m,double * t,double * H,int * ldh,double * wsp,int * lwsp,int * ipiv,int * iexph,int *ns,int *iflag );
}

void RateModel::setup_fortran_P_whichQ(int which, mat & P, double t) {
    Q = Qs[which];
    setup_fortran_P(P,t,false);
}
*/
//taking out fortran
/*
 * runs the basic padm fortran expokit full matrix exp
 *
void RateModel::setup_fortran_P(mat & P, double t, bool store_p_matrices) {
    //
    //  return P, the matrix of dist-to-dist transition probabilities,
    //  from the model's rate matrix (Q) over a time duration (t)
    //
    int ideg = 6;
    int m = Q.n_rows; // square so you only need the rows
    int ldh = m;
    double tol = 1;
    int iflag = 0;
    int lwsp = 4*m*m+6+1;
    double * wsp = new double[lwsp];
    int * ipiv = new int[m];
    int iexph = 0;
    int ns = 0;
    double * H = new double [m*m];
    convert_matrix_to_single_row_for_fortran(Q, t, H);
    wrapdgpadm_(&ideg, &m, &tol, H, &ldh, wsp, &lwsp, ipiv, &iexph, &ns, &iflag);

    for (int i=0; i < m; i++) {
        for (int j=0; j < m; j++) {
            P(i, j) = wsp[iexph + (j-1) * m + (i-1) + m];
        }
    }
    delete [] wsp;
    delete [] ipiv;
    delete [] H;
    for (int i=0; i < nstates; i++) {
        double sum = 0.0;
        for (int j=0; j < nstates; j++) {
            sum += P(i, j);
        }
        for (int j=0; j < nstates; j++) {
            P(i, j) = (P(i, j)/sum);
        }
    }
    
}
*/
void RateModel::set_sameQ(bool s) {
    sameQ = s;
}

void update_simple_goldman_yang_q(mat * inm, double K, double w, mat & bigpibf,mat &bigpiK, mat & bigpiw) {
    double s = 0;
    for (int i=0; i < 61; i++) {
        for (int j=0; j < 61; j++) {
            if (bigpibf(i, j) == 0) {
                (*inm)(i, j) = 0;
            } else {
                (*inm)(i, j) = 1/61.;
                if (bigpiK(i, j) != 0) {
                    (*inm)(i, j) *= K;
                }
                if (bigpiw(i, j) != 0) {
                    (*inm)(i, j) *= w;
                }
            } 
            if (i == j) {
                (*inm)(i, j) = 0;
            }
            s += (*inm)(i, j);
        }
    }
    (*inm) = (*inm)/s;
    for (unsigned int i=0; i < (*inm).n_rows; i++) {
        double su = 0;
        for (unsigned int j=0; j < (*inm).n_cols; j++) {
            if (i != j) {
            su += (*inm)(i, j);
            }
        }
        (*inm)(i, i) = 0-su;
    }
    (*inm) = (*inm)/(1/61.);
    //(*inm) = trans((*inm));
}

bool test_transition(char a, char b) {
    bool ret = false;
    if ((a == 'A' &&  b == 'G') || (a == 'C' &&  b == 'T') || (a == 'G' &&  b == 'A') || (a == 'T' &&  b == 'C')) {
        ret = true;
    }
    return ret;
}

void generate_bigpibf_K_w(mat * bf, mat * K, mat * w,map<string, string> & codon_dict, 
    map<string, vector<int> > & codon_index, vector<string> & codon_list) {
    for (int i=0; i < 61; i++) {
        for (int j=0; j < 61; j++) {
            int diff = 0;
            bool transit = false;
            bool nonsyn = false;
            for (int m=0; m < 3; m++) {
                if (codon_list[i][m] != codon_list[j][m]) {
                    diff += 1;
                }
                transit = test_transition(codon_list[i][m],codon_list[j][m]);
            }
            if (diff > 1) {
                (*bf)(i, j) = 0;
            } else {
                if (codon_dict[codon_list[i]] != codon_dict[codon_list[j]]) {
                    nonsyn = true;
                }
                (*bf)(i, j) = 1;
                if (transit) {
                    (*K)(i, j) = 1;
                }
                if (nonsyn) {
                    (*w)(i, j) = 1;
                }
            }
        }
    }
}

//take out fortran
/*
void convert_matrix_to_single_row_for_fortran(mat & inmatrix, double t, double * H) {
    int count = 0;
    for (unsigned int i=0; i < inmatrix.n_cols; i++) {
        for (unsigned int j=0; j < inmatrix.n_cols; j++) {
            H[i+(j*inmatrix.n_cols)] = inmatrix(i, j)*t;
            count += 1;
        }
    }
}
*/
