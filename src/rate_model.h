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
	cx_mat lasteigval;
	cx_mat lasteigvec;
	cx_mat eigval;
	cx_mat eigvec;
	mat lasteigval_simple;
	mat lasteigvec_simple;
	mat eigval_simple;
	mat eigvec_simple;
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
	void set_Q(mat & inQ);
	void set_Q_diag();
	void set_Q_cell(int,int,double);
	void setup_Q();
	void setup_Q(vector<vector<double> > & inQ);
	void setup_Q(mat & inQ);
	mat & get_Q();
	void set_sameQ(bool);
	cx_mat setup_P(double,bool);
	mat setup_P_simple(double,bool);

	/*
	 * get things from stmap
	 */
	//this should be used for getting the eigenvectors and eigenvalues
	bool get_eigenvec_eigenval_from_Q(cx_mat * eigenvalues, cx_mat * eigenvectors);
	void get_eigenvec_eigenval_from_Q_simple(mat * eigenvalues, mat * eigenvectors);

};
void update_simple_goldman_yang_q(mat * inm, double K, double w, mat & bigpibf,mat &bigpiK, mat & bigpiw);
bool test_transition(char a, char b);
void generate_bigpibf_K_w(mat * bf, mat * K, mat * w,map<string, string> & codon_dict, map<string, vector<int> > & codon_index, vector<string> & codon_list);

#endif
