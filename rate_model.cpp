#include <vector>
#include <string>
#include <map>

#include "rate_model.h"

#include <armadillo>
using namespace arma;

inline int signof(double d){return d >= 0 ? 1 : -1;}
inline double roundto(double in){return floor(in*(1000)+0.5)/(1000);}

RateModel::RateModel(int _nstates):nstates(_nstates),Q(_nstates,_nstates),labels(),Q_mask(){
	setup_Q();
}


void RateModel::set_Q_cell(int from, int to, double num){
	Q(from,to) = num;
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
}

void RateModel::setup_Q(mat & inQ){
	for(unsigned int i=0;i<Q.n_rows;i++){
		for(unsigned int j=0;j<Q.n_cols;j++){
			Q(i,j) = inQ(i,j);
		}
	}
	set_Q_diag();
}

mat & RateModel::get_Q(){
	return Q;
}

cx_mat RateModel::setup_P(double bl,bool store_p_matrices){
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
	return isImag;
}
