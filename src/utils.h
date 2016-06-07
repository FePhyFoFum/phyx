#ifndef _UTILS_H_
#define _UTILS_H_

#include <vector>
#include <iterator>
#include "superdouble.h"

void check_file_exists (const string& filename);

void tokenize (const string& str, vector <string>& tokens, const string& delimiters = " ");
vector <string> tokenize (string const& input);
void trim_spaces (string& str);
bool is_number (const string &);

unsigned int get_clock_seed ();

vector <double> parse_double_comma_list (string& str);

// do stuff over vectors
double sum (vector <double> & in);
int sum (vector <int> & in);
Superdouble calculate_vector_Superdouble_sum (vector <Superdouble> & in);
double mean (vector <double> & in);
double variance (vector <double> & in);
vector <int> sum (vector <int> & vec1, vector <int> & vec2);

vector <vector <double> > processRateMatrixConfigFile (string filename, int numareas);
int random_int_range (int min, int max);

vector <int> sample_without_replacement (int const& numTotal, int const& numSample);

void print_error (char * pname, char arg);
bool test_logical (vector <int> & matA, vector <int> & matB);

int sum_matrix_col (vector <vector <int> > & matrix, int col);
int sum_matrix_col_negs (vector <vector <int> > & matrix, int col);

string get_string_vector(vector <string> &sts);
string get_string_vector(vector <int> &sts);

void replace_all (string& str, string const& origSubStr, string const& replSubStr);
void replace_each (string& str, string const& badChars, string const& replSubStr);

string get_valid_newick_label (string const& inLabel);
string get_valid_nexus_label (string const& inLabel);
string get_safe_taxon_label (string const& inLabel);
void quotify_label (string & token);

// not currently used
template<typename T> void print_vector (vector <T> & vec) {
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(std::cout, " "));
    cout << endl;
}

unsigned int calc_hamming_dist (string const& s1, string const& s2);

double logn (double x, double base);

bool essentially_equal (double a, double b);
bool all_equal (vector <double> vals);

#endif /* _UTILS_H_ */
