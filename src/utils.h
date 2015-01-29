#ifndef _UTILS_H_
#define _UTILS_H_

#include <vector>
#include <string>
#include <map>
#include "superdouble.h"

using namespace std;

void tokenize(const string& str, vector<string>& tokens,const string& delimiters = " ");
vector <string> tokenize (string const& input);
void trim_spaces(string& str);
bool is_number(const string &);
double calculate_vector_double_sum(vector<double> & in);
double calculate_vector_double_mean(vector<double> & in);
int calculate_vector_int_sum(vector<int> & in);
vector<vector<double> > processRateMatrixConfigFile(string filename, int numareas);
Superdouble calculate_vector_Superdouble_sum(vector<Superdouble> & in);
int random_int_range(int min, int max);
void print_error(char * pname, char arg);

int sum_matrix_col(vector<vector<int> > & matrix,int col);
int sum_matrix_col_negs(vector<vector<int> > & matrix,int col);
bool test_logical(vector<int> & matA,vector<int> & matB);
double sum(vector<double> &inm);
double sum(vector<int> &inm);

string get_string_vector(vector<string> &sts);
string get_string_vector(vector<int> &sts);

double logn(double x,double base);

bool essentially_equal (double a, double b);
bool all_equal (vector <double> vals);

#endif /* _UTILS_H_ */
