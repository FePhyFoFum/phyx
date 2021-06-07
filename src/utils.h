#ifndef PX_UTILS_H
#define PX_UTILS_H

#include <string>
#include <vector>
#include <iterator>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <regex>

#include "superdouble.h"


// TODO organize this stuff!

void check_file_exists (const std::string& filename);
void check_inout_streams_identical (char * in, char * out);
int string_to_int (const std::string& in, const std::string& arg);
long int string_to_long_int (const std::string& in, const std::string& arg);

std::string string_to_upper (const std::string& instr);
std::string string_to_lower (const std::string& instr);
float string_to_float (const std::string& in, const std::string& arg);
double string_to_double (const std::string& in, const std::string& arg);

void tokenize (const std::string& str, std::vector<std::string>& tokens,
        const std::string& delimiters = " ");
std::vector<std::string> tokenize (const std::string& input);
void trim_spaces (std::string& str);
bool check_comment_line (const std::string& line);
bool is_number (const std::string&);
unsigned int get_longest_label (std::vector<std::string>& labels);

std::istream& getline_safe(std::istream& is, std::string& t);

unsigned int factorial (unsigned int n);
unsigned long int doublefactorial (unsigned int n);

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

double v_median (std::vector<double>& in);
double v_mean (std::vector<double>& in);
void v_mean_variance (std::vector<double>& in, double& mn, double& varr);
double v_variance (std::vector<double>& in);
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

// get rid of any duplicate elements
// assumes number of elements is not huge (otherwise a set might be better)
template<typename T> std::vector<T> get_unique_elements (std::vector<T> vec) {
    sort(vec.begin(), vec.end());
    vec.erase(unique(vec.begin(), vec.end()), vec.end());
    return vec;
}

// find elements present in the first vector but not the second. does not do symmetric search
template<typename T> std::vector<T> in_first_not_second (std::vector<T> firstv, std::vector<T> secondv) {
    std::vector<T> res(firstv.size());
    std::sort(firstv.begin(), firstv.end());
    std::sort(secondv.begin(), secondv.end());
    
    typename std::vector<T>::iterator it;
    it = std::set_difference(firstv.begin(), firstv.end(), secondv.begin(), secondv.end(), res.begin());
    res.resize(static_cast<unsigned long>(it-res.begin()));
    
    //std::remove_copy_if(firstv.begin(), firstv.end(), std::back_inserter(res),
    //    [&secondv](const T& arg) {return (std::find(secondv.begin(), secondv.end(), arg) != secondv.end());});
    return res;
}

std::vector<double> string_v_to_double_v (const std::vector<std::string>& in);

//------------------------------------------------------------------------//

std::vector<std::vector<double> > processRateMatrixConfigFile (const std::string& filename, int numstates);

unsigned int random_int_range (const unsigned int& min, const unsigned int& max);
std::vector<unsigned int> sample_without_replacement (const unsigned int& numTotal,
        const unsigned int& numSample);

void print_error (char * pname);
bool test_logical (std::vector<int>& matA, std::vector<int>& matB);
bool test_logical (std::vector<int>& matA, std::vector<int>& matB, bool edgewise);

int sum_matrix_col (std::vector<std::vector<int> >& matrix, int col);
int sum_matrix_col_negs (std::vector<std::vector<int> >& matrix, int col);

std::string get_string_vector (std::vector<std::string>& sts);
std::string get_string_vector (std::vector<int>& sts);

void replace_all (std::string& str, const std::string& origSubStr, const std::string& replSubStr);
void replace_each (std::string& str, const std::string& badChars, const std::string& replSubStr);

bool check_nexus_comment (std::string line);
void process_nexus_comment (std::istream& stri, std::string& tline);

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

bool essentially_equal_doubles (double a, double b);

// template version to accommodate different precision types
template<typename T> bool essentially_equal (T a, T b) {
    T EPSILON = 1e-7;
    /*
    std::cout << "fabs(a - b) = " << fabs(a - b) << std::endl;
    std::cout << "ApproximatelyEqual: "
        << (fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * EPSILON)) << std::endl;
    std::cout << "EssentiallyEqual: "
        << (fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * EPSILON)) << std::endl;
    */
    
    // not used
    /*
    if (fabs(a - b) <= std::max(EPSILON, EPSILON * std::max(std::abs(a), std::abs(b)))) {
        equal = true;
    }
    */
    bool equal = (std::abs(a - b) <= ( (std::abs(a) > std::abs(b) ? std::abs(b) : std::abs(a)) * EPSILON));
    
    return equal;
}


bool all_equal (std::vector<double> vals);
bool all_equal (std::vector<int> vals);
bool all_equal (std::vector<unsigned int> vals);


// a basic poll checker for stream inputs
bool check_for_input_to_stream ();

std::vector<std::string> get_complement (std::vector<std::string>& a, std::vector<std::string>& b);

std::string peek_line (std::istream& pios);
std::vector<std::string> peek_lines (std::istream& pios, const int& n);

std::vector<std::string> regex_search_labels (const std::vector<std::string>& names,
        const std::string& pattern);

#endif /* PX_UTILS_H */
