#ifndef UTILS_H_
#define UTILS_H_

#include <vector>
#include <string>
#include <map>
#include "superdouble.h"

using namespace std;

void Tokenize(const string& str, vector<string>& tokens,const string& delimiters = " ");
void TrimSpaces( string& str);
double calculate_vector_double_sum(vector<double> & in);
int calculate_vector_int_sum(vector<int> & in);
vector<vector<double> > processRateMatrixConfigFile(string filename, int numareas);
Superdouble calculate_vector_Superdouble_sum(vector<Superdouble> & in);
int random_int_range(int min, int max);

#endif /* UTILS_H_ */
