#ifndef UTILS_H_
#define UTILS_H_

#include <vector>
#include <string>
#include <map>

using namespace std;
void Tokenize(const string& str, vector<string>& tokens,const string& delimiters = " ");
void TrimSpaces( string& str);
double calculate_vector_double_sum(vector<double> & in);

#endif /* UTILS_H_ */
