#ifndef _UTILS_H_
#define _UTILS_H_

#include <string>
#include <vector>
#include <iterator>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "superdouble.h"


// TODO organize this stuff


void check_file_exists (const std::string& filename);
void check_inout_streams_identical (char * in, char * out);
int string_to_int (const std::string& in, const std::string& arg);

std::string string_to_upper (const std::string& str);
float string_to_float (const std::string& in, const std::string& arg);

void tokenize (const std::string& str, std::vector<std::string>& tokens,
        const std::string& delimiters = " ");
std::vector<std::string> tokenize (const std::string& input);
void trim_spaces (std::string& str);
bool check_comment_line (const std::string& line);
bool is_number (const std::string&);

int factorial (int n);
int doublefactorial (int n);

unsigned int get_clock_seed ();

// template version. pass result vector as arg; will delete anything already in res
template<typename T> void parse_comma_list (std::string& str, std::vector<T>& res) {
    std::stringstream ss(str);
    T i;
    res.clear();
    while (ss >> i) {
        res.push_back(i);
        bool done = false;
        while (!done) { // shouldn't be any spaces, but let's be safe
            if (ss.peek() == ' ' || ss.peek() == ',') {
                ss.ignore();
            } else {
                done = true;
            }
        }
    }
}

//------------------------------------------------------------------------//
// do stuff over vectors
double sum (std::vector<double>& in);
int sum (std::vector<int>& in);

int count_zeros (std::vector<int>& in);
Superdouble calculate_vector_Superdouble_sum (std::vector<Superdouble>& in);
double mean (std::vector<double>& in);
double variance (std::vector<double>& in);
std::vector<int> sum (std::vector<int>& vec1, std::vector<int>& vec2);

template<typename T> std::vector<T> sum_vectors_elementwise (std::vector<T>& vec1, std::vector<T>& vec2) {
    // bail if sequences are of different lengths.
    if (vec1.size() != vec2.size()) {
      throw std::invalid_argument(
          "Vectors must be of equal length"
      );
    }
    std::vector<T> res = vec1;
    std::transform(res.begin(), res.end(), vec2.begin(), res.begin(), std::plus<T>());
    return res;
}

std::vector<double> average_vectors_elementwise (std::vector<double>& vec1, std::vector<double>& vec2);

//------------------------------------------------------------------------//

std::vector<std::vector<double> > processRateMatrixConfigFile (std::string filename, int numareas);
int random_int_range (int min, int max);

std::vector<int> sample_without_replacement( const int& numTotal,const int& numSample);

void print_error (char * pname, char arg);
bool test_logical (std::vector<int>& matA, std::vector<int>& matB);
bool test_logical (std::vector<int>& matA, std::vector<int>& matB, bool edgewise);

int sum_matrix_col (std::vector<std::vector<int> >& matrix, int col);
int sum_matrix_col_negs (std::vector<std::vector<int> >& matrix, int col);

std::string get_string_vector (std::vector<std::string>& sts);
std::string get_string_vector (std::vector<int>& sts);

void replace_all (std::string& str, const std::string& origSubStr, const std::string& replSubStr);
void replace_each (std::string& str, const std::string& badChars, const std::string& replSubStr);

std::string get_valid_newick_label (const std::string& inLabel);
std::string get_valid_nexus_label (const std::string& inLabel);
std::string get_safe_taxon_label (const std::string& inLabel);
void quotify_label (std::string& token);

template<typename T> void print_vector (std::vector<T>& vec) {
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(std::cout, " "));
    std::cout << std::endl;
}

unsigned int calc_hamming_dist (const std::string& s1, const std::string& s2);

double logn (double x, double base);

bool essentially_equal (double a, double b);
bool all_equal (std::vector<double> vals);

// a basic poll checker for stream inputs
bool check_for_input_to_stream ();

std::string peek_line (std::istream* pios);

#endif /* _UTILS_H_ */
