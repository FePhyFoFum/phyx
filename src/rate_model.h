#ifndef RATEMODEL_H_
#define RATEMODEL_H_

#include <vector>
#include <string>
#include <map>
using namespace std;

#include <armadillo>
using namespace arma;

class RateModel{

private:

	mat Q;
	vector<string> labels;
	vector< vector<double> > Q_mask;
	bool sameQ;
	cx_mat * lasteigval;
	cx_mat * lasteigvec;
	bool lastImag;

public:
	RateModel(int);
	/*
		 storing once optimization has occured
		 map of bl and p matrix
		 map<branch length,p matrix>
	 */
	int nstates;
	map<double, cx_mat> stored_p_matrices;
	bool neg_p;

	void set_Q_diag();
	void set_Q_cell(int,int,double);
	void setup_Q();
	void setup_Q(vector<vector<double> > & inQ);
	void setup_Q(mat & inQ);
	mat & get_Q();

	cx_mat setup_P(double,bool);

	/*
	 * get things from stmap
	 */
	//this should be used for getting the eigenvectors and eigenvalues
	bool get_eigenvec_eigenval_from_Q(cx_mat * eigenvalues, cx_mat * eigenvectors);
};

#endif
