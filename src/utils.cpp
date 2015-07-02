
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <sstream>
//#include <iterator>

#include <cmath>
#include <limits>
#include <functional>

using namespace std;

// set threshold. should maybe use superdouble.
double EPSILON = 1e-7;

// other stuff i am playing around with. temporary. 
/*
#define MEPSILON numeric_limits<double>::min() // 2.22507e-308
#define GEPSILON numeric_limits<double>::round_error() // 0.5
#define FLEPSILON numeric_limits<float>::epsilon() // 1.19209e-07
#define MANTISSA std::numeric_limits<double>::digits10 // 15
#define FLAMANTISSA std::numeric_limits<float>::digits10 // 6
#define double THRESH = e-10;
*/

#include "utils.h"
#include "superdouble.h"

void tokenize(const string& str, vector <string>& tokens, const string& delimiters) {
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

// works a little differently than above; don't need to trim spaces
vector <string> tokenize (string const& input) {
    vector <string> tokens;
    string temp;
    istringstream str(input);
    while (str >> temp) {
        tokens.push_back(temp);
    }
    return tokens;
}

void trim_spaces(string& str) {
    // Trim Both leading and trailing spaces
    size_t startpos = str.find_first_not_of(" \t\r\n"); // Find the first character position after excluding leading blank spaces
    size_t endpos = str.find_last_not_of(" \t\r\n"); // Find the first character position from reverse af

    // if all spaces or empty return an empty string
    if ((string::npos == startpos) || ( string::npos == endpos)) {
        str = "";
    } else {
        str = str.substr(startpos, endpos - startpos + 1);
    }
    /*
    // Code for  Trim Leading Spaces only
    size_t startpos = str.find_first_not_of(\t); // Find the first character position after excluding leading blank spaces
    if (string::npos != startpos)
    str = str.substr( startpos );
    */

    /*
    // Code for Trim trailing Spaces only
    size_t endpos = str.find_last_not_of(\t); // Find the first character position from reverse af
    if (string::npos != endpos)
    str = str.substr(0, endpos + 1);
    */
}

bool is_number(const string& s) {
    string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) {
        ++it;
    }
    return !s.empty() && it == s.end();
}

vector <vector <double> > processRateMatrixConfigFile(string filename, int numstates) {
    vector <double> cols(numstates,1);
    vector <vector <double> > ratematrix = vector <vector <double> > (numstates,cols);
    //read file
    ifstream ifs(filename.c_str());
    string line;
    int fromarea = 0;
    while (getline(ifs,line)) {
        if (line.size() > 3) {
            vector <string> tokens;
            string del(" ,\t");
            tokens.clear();
            tokenize(line, tokens, del);
            for (unsigned int j=0; j < tokens.size(); j++) {
                trim_spaces(tokens[j]);
            }
            for (unsigned int j=0; j < tokens.size(); j++) {
                ratematrix[fromarea][j] = atof(tokens[j].c_str());
            }
            if (fromarea < numstates-1) {
                fromarea += 1;
            } else {
                fromarea = 0;
            }
        }
    }
    ifs.close();
    return ratematrix;
}

Superdouble calculate_vector_Superdouble_sum(vector <Superdouble> & in) {
    Superdouble sum = 0;
    for (unsigned int i=0; i < in.size(); i++) {
        sum += in[i];
        // cout << in[i] << " sum:" << sum << endl;
    }
    // cout << "endsum:" << sum << endl;
    return sum;
}

int random_int_range(int min, int max) {
    return min + (rand() % (int)(max - min + 1));
}

// arg below is always '?'. besides, getopt prints errors to sterr
void print_error(char * pname, char arg) {
    // cout << pname <<": invalid option -- '" << arg << "'" << endl;
    cout << "Try `" << pname << " --help' for more information." << endl;
}


int sum_matrix_col(vector <vector <int> > & matrix,int col) {
    int x=0;
    for (unsigned int i=0; i < matrix.size(); i++) {
        x += matrix[i][col];
    }
    return x;
}

int sum_matrix_col_negs(vector <vector <int> > & matrix, int col) {
    int x=0;
    for (unsigned int i=0; i < matrix.size(); i++) {
        if (matrix[i][col] < 0) {
            x += matrix[i][col];
        }
    }
    return x;
}

bool test_logical(vector <int> & matA, vector <int> & matB) {
    bool test = false;
    int match1 = 0;
    unsigned int numdiffs = 0;
    for (unsigned int i=0; i < matA.size(); i++) {
        //added the -1, can take out if something gets fishy
        if (((matA[i] == 1) && (matB[i] == 1))) {
            match1 += 1;
        } else {
            numdiffs += 1;
        }
    }
    //had to change because of negatives
    if ((match1 != sum(matA)) && (match1 != sum(matB))) {
        if (numdiffs != matA.size()) {
            test = true;
        }
    }
    return test;
}

// simple math on vectors
double mean (vector <double> & in) {
    return sum (in) / (double)in.size();
}

double sum (vector <double> & in) {
    return accumulate(in.begin(), in.end(), 0.0);
}

int sum (vector <int> & in) {
    return accumulate(in.begin(), in.end(), 0);
}

vector <int> sum (vector <int> & vec1, vector <int> & vec2) {
    
    // bail if sequences are of different lengths. should be caught earlier than this
    if (vec1.size() != vec2.size()) {
      throw std::invalid_argument(
          "Vectors must be of equal length"
      );
    }
    
    vector <int> res = vec1;
    std::transform(res.begin(), res.end(), vec2.begin(), res.begin(), std::plus<int>());
    return res;
}

string get_string_vector(vector <string> &sts) {
    string rets;
    for (unsigned int i=0; i < sts.size(); i++) {
        rets += sts[i]+ " ";
    }
    return rets;
}

string get_string_vector(vector <int> & sts) {
    string rets;
    for (unsigned int i=0; i < sts.size(); i++) {
        rets += to_string(sts[i]) + " ";
    }
    return rets;
}

// Hamming (edit) distance
unsigned int calc_hamming_dist (string const& s1, string const& s2) {
    
    if (s1 == s2) {
        return 0;
    }
    
    // bail if sequences are of different lengths. should be caught earlier than this
    if (s1.size() != s2.size()) {
        throw std::invalid_argument(
          "Hamming distances are only defined for strings of equal length"
      );
    }

    return std::inner_product(
        s1.begin(), s1.end(), s2.begin(), 0, std::plus<unsigned int>(),
        std::not2(std::equal_to<std::string::value_type>())
    );
}

double logn(double x, double base) {
    return log10(x)/log10(base);
}

// check if 2 doubles are equal within some tolerance.
bool essentially_equal (double a, double b) {
    bool equal = false;
    
    /*
    cout << "fabs(a - b) = " << fabs(a - b) << endl;
    cout << "ApproximatelyEqual: " <<  (fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * EPSILON)) << endl;
    cout << "EssentiallyEqual: " <<  (fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * EPSILON)) << endl;
    */
    
    if (fabs(a - b) <= max(EPSILON, EPSILON * max(abs(a), abs(b)))) {
        equal = true;
    }
    
    equal = (fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * EPSILON));
    
    return equal;
}

bool all_equal (vector <double> vals) {
    bool equal = false;
    vector <double>::iterator it = find_if_not(vals.begin()+1, vals.end(), bind(essentially_equal, placeholders::_1, vals[0]));
    if (it == end(vals)) {
        equal = true;
    }
    return equal;
}

