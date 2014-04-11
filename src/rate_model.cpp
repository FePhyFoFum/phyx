#include <vector>
#include <string>
#include <map>

#include "rate_model.h"

#include <armadillo>
using namespace arma;

inline int signof(double d){return d >= 0 ? 1 : -1;}
inline double roundto(double in){return floor(in*(1000)+0.5)/(1000);}

RateModel::RateModel(int _nstates):Q(_nstates,_nstates),labels(),Q_mask(),nstates(_nstates){
    setup_Q();
    sameQ = false;
}


void RateModel::set_Q_cell(int from, int to, double num){
    Q(from,to) = num;
    sameQ = false;
}

void RateModel::set_Q_diag(){
    for(unsigned int i=0;i<Q.n_rows;i++){
	double su = 0;
	for(unsigned int j=0;j<Q.n_cols;j++){
	    if(i!=j){
		su += Q(i,j);
	    }
	}
	Q(i,i) = 0-su;
    }
    sameQ = false;
}

void RateModel::setup_Q(){
    Q.fill(0);
    for(unsigned int i=0;i<Q.n_rows;i++){
	for(unsigned int j=0;j<Q.n_cols;j++){
	    if(i!=j){
		Q(i,j) = 1./nstates;
	    }else{
		Q(i,j) = -(1./nstates * (nstates-1));
	    }
	}
    }
    sameQ = false;
}

void RateModel::setup_Q(vector<vector<double> > & inQ){
    for(unsigned int i=0;i<Q.n_rows;i++){
	for(unsigned int j=0;j<Q.n_cols;j++){
	    Q(i,j) = inQ[i][j];
	}
    }
    for(unsigned int i=0;i<Q.n_rows;i++){
	colvec a = (sum(Q,1));
	Q(i,i) = -(a(i)-Q(i,i));
    }
    sameQ = false;
}

void RateModel::setup_Q(mat & inQ){
    for(unsigned int i=0;i<Q.n_rows;i++){
	for(unsigned int j=0;j<Q.n_cols;j++){
	    Q(i,j) = inQ(i,j);
	}
    }
    set_Q_diag();
    sameQ = false;
}

void RateModel::set_Q(mat & inQ){
    for(unsigned int i=0;i<Q.n_rows;i++){
	for(unsigned int j=0;j<Q.n_cols;j++){
	    Q(i,j) = inQ(i,j);
	}
    }
    sameQ = false;
}

mat & RateModel::get_Q(){
    return Q;
}

cx_mat RateModel::setup_P(double bl,bool store_p_matrices){
//	sameQ = false;
    cx_mat eigvec(nstates,nstates);eigvec.fill(0);
    cx_mat eigval(nstates,nstates);eigval.fill(0);
    bool isImag = get_eigenvec_eigenval_from_Q(&eigval, &eigvec);
    //cout << eigval << endl;
    //cout << eigvec << endl;
    for(int i=0;i<nstates;i++){
	eigval(i,i) = exp(eigval(i,i) * bl);
    }
    cx_mat C_inv = inv(eigvec);
    cx_mat P = eigvec * eigval * C_inv;
    neg_p = false;
    for(int i=0;i<P.n_rows;i++){
	for(int j=0;j<P.n_cols;j++){
	    if (real(P(i,j))<0)
		neg_p = true;
	}
    }
    if(store_p_matrices == true){
	stored_p_matrices[bl] = P;	
    }
    return P;
}

/*
 * this should be used to caluculate the eigenvalues and eigenvectors
 * as U * Q * U-1 -- eigen decomposition
 *
 * this should use the armadillo library
 */
bool RateModel::get_eigenvec_eigenval_from_Q(cx_mat * eigval, cx_mat * eigvec){
    if(sameQ==true){
	eigval = lasteigval;
	eigvec = lasteigvec;
	return lastImag;
    }
    mat tQ(nstates,nstates); tQ.fill(0);
    for(unsigned int i=0;i<Q.n_rows;i++){
	for(unsigned int j=0;j<Q.n_cols;j++){
	    tQ(i,j) = Q(i,j);
	}
    }
    cx_colvec eigva;
    cx_mat eigve;
    eig_gen(eigva,eigve,tQ);
    bool isImag = false;
    for(unsigned int i=0;i<Q.n_rows;i++){
	for(unsigned int j=0;j<Q.n_cols;j++){
	    if(i==j)
		(*eigval)(i,j) = eigva(i);
	    else
		(*eigval)(i,j) = 0;
	    (*eigvec)(i,j) = eigve(i,j);
	    if(imag((*eigvec)(i,j)) > 0 || imag((*eigval)(i,j)))
		isImag = true;
	}
    }
    lasteigval = eigval;
    lasteigvec = eigvec;
    lastImag = isImag;
    return isImag;
}


void update_simple_goldman_yang_q(mat * inm, double K, double w, mat & bigpibf,mat &bigpiK, mat & bigpiw){
    double s = 0;
    for (int i=0;i<61;i++){
	for (int j=0;j<61;j++){
	    if (bigpibf(i,j)==0)
		(*inm)(i,j) = 0;
	    else{
		(*inm)(i,j) = 1/61.;
		if(bigpiK(i,j)!=0)
		    (*inm)(i,j) *= K;
		if(bigpiw(i,j)!=0)
		    (*inm)(i,j) *= w;
	    } 
	    if(i == j){
		(*inm)(i,j) = 0;
	    }
	    s += (*inm)(i,j);
	}
    }
    (*inm) = (*inm)/s;
    for(unsigned int i=0;i<(*inm).n_rows;i++){
	double su = 0;
	for(unsigned int j=0;j<(*inm).n_cols;j++){
	    if(i!=j){
		su += (*inm)(i,j);
	    }
	}
	(*inm)(i,i) = 0-su;
    }
    (*inm) = (*inm)/(1/61.);
    //(*inm) = trans((*inm));
}

bool test_transition(char a, char b){
    bool ret = false;
    if ((a == 'A' &&  b == 'G') || (a == 'C' &&  b == 'T') || (a == 'G' &&  b == 'A') || (a == 'T' &&  b == 'C')){
	ret = true;
    }
    return ret;
}

void generate_bigpibf_K_w(mat * bf, mat * K, mat * w,map<string, string> & codon_dict, map<string, vector<int> > & codon_index, vector<string> & codon_list){
    for(int i=0;i<61;i++){
	for(int j=0;j<61;j++){
	    int diff = 0;
	    bool transit = false;
	    bool nonsyn = false;
	    for(int m=0;m<3;m++){
		if (codon_list[i][m] != codon_list[j][m])
		    diff += 1;
		transit = test_transition(codon_list[i][m],codon_list[j][m]);
	    }
	    if(diff > 1){
		(*bf)(i,j) = 0;
	    }else{
		if (codon_dict[codon_list[i]] != codon_dict[codon_list[j]])
		    nonsyn = true;
		(*bf)(i,j) = 1;
		if (transit)
		    (*K)(i,j) = 1;
		if(nonsyn)
		    (*w)(i,j) = 1;
	    }
	}
    }
}
