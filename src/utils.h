#ifndef _UTILS_H_
#define _UTILS_H_

#include <vector>
#include <iterator>
#include <iostream>
#include <sstream>
#include <string>


#include "superdouble.h"

void check_file_exists (const std::string& filename);
void check_inout_streams_identical (char * in, char * out);
int string_to_int (std::string const& in, std::string const& arg);

std::string string_to_upper (std::string const& str);
float string_to_float (std::string const& in, std::string const& arg);

void tokenize (const std::string& str, std::vector <std::string>& tokens, const std::string& delimiters = " ");
std::vector <std::string> tokenize (std::string const& input);
void trim_spaces (std::string & str);
bool check_comment_line (std::string const& line);
bool is_number (const std::string &);

int factorial (int n);
int doublefactorial(int n);

unsigned int get_clock_seed ();

//vector <double> parse_double_comma_list (string& str);
//vector <int> parse_int_comma_list (string& str);

// template version. pass result vector as arg; will delete anything already in res
template<typename T> void parse_comma_list (std::string& str, std::vector <T> & res) {
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

// do stuff over vectors
double sum (std::vector <double> & in);
int sum (std::vector <int> & in);
int sum_zeros (std::vector <int> & in);
Superdouble calculate_vector_Superdouble_sum (std::vector <Superdouble> & in);
double mean (std::vector <double> & in);
double variance (std::vector <double> & in);
std::vector <int> sum (std::vector <int> & vec1, std::vector <int> & vec2);

std::vector <std::vector <double> > processRateMatrixConfigFile (std::string filename, int numareas);
int random_int_range (int min, int max);

std::vector <int> sample_without_replacement (int const& numTotal, int const& numSample);

void print_error (char * pname, char arg);
bool test_logical (std::vector <int> & matA, std::vector <int> & matB);
bool test_logical (std::vector <int> & matA, std::vector <int> & matB, bool edgewise);

int sum_matrix_col (std::vector <std::vector <int> > & matrix, int col);
int sum_matrix_col_negs (std::vector <std::vector <int> > & matrix, int col);

std::string get_string_vector(std::vector <std::string> &sts);
std::string get_string_vector(std::vector <int> &sts);

void replace_all (std::string& str, std::string const& origSubStr, std::string const& replSubStr);
void replace_each (std::string& str, std::string const& badChars, std::string const& replSubStr);

std::string get_valid_newick_label (std::string const& inLabel);
std::string get_valid_nexus_label (std::string const& inLabel);
std::string get_safe_taxon_label (std::string const& inLabel);
void quotify_label (std::string & token);

// not currently used (but cool)
template<typename T> void print_vector (std::vector <T> & vec) {
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(std::cout, " "));
    std::cout << std::endl;
}

unsigned int calc_hamming_dist (std::string const& s1, std::string const& s2);

double logn (double x, double base);

bool essentially_equal (double a, double b);
bool all_equal (std::vector <double> vals);


// a basic poll checker for stream inputs
bool check_for_input_to_stream();

#endif /* _UTILS_H_ */
