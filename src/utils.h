#ifndef UTILS_H_
#define UTILS_H_

#include <vector>
#include <string>
#include <map>
#include "superdouble.h"

using namespace std;

void tokenize(const string& str, vector<string>& tokens,const string& delimiters = " ");
void trim_spaces( string& str);
bool is_number(const string &);
double calculate_vector_double_sum(vector<double> & in);
double calculate_vector_double_mean(vector<double> & in);
int calculate_vector_int_sum(vector<int> & in);
vector<vector<double> > processRateMatrixConfigFile(string filename, int numareas);
Superdouble calculate_vector_Superdouble_sum(vector<Superdouble> & in);
int random_int_range(int min, int max);
void print_error(char * pname, char * arg);

#endif /* UTILS_H_ */
