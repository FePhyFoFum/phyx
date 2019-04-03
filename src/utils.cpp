
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>
#include <functional>
#include <ctime>
#include <chrono>
#include <poll.h>
#include <unistd.h>

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


void check_file_exists (const string& filename) {
    std::ifstream infile(filename);
    if (!infile.good()) {
        cout << "File '" << filename << "' not found. Exiting." << endl;
        exit(0);
    }
}

void check_inout_streams_identical (char * in, char * out) {
    if ((string)in == (string)out) {
        cout << "Warning: input and output file names must differ (streams!). Exiting." << endl;
        exit(0);
    }
}

// convenience func. returns copy so original can still be used
string string_to_upper (string const& instr) {
    string outstr = instr;
    std::transform(outstr.begin(), outstr.end(), outstr.begin(), ::toupper);
    return outstr;
}

void tokenize (const string& str, vector <string>& tokens, const string& delimiters) {
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


// deprecated in favour of template version (in header file)
//vector <double> parse_double_comma_list (string& str) {
//    vector <double> res;
//    std::stringstream ss(str);
//    double i;
//    while (ss >> i) {
//        res.push_back(i);
//        bool done = false;
//        while (!done) { // shouldn't be any spaces, but let's be safe
//            if (ss.peek() == ' ' || ss.peek() == ',') {
//                ss.ignore();
//            } else {
//                done = true;
//            }
//        }
//    }
//    return res;
//}

// deprecated in favour of template version (in header file)
//vector <int> parse_int_comma_list (string& str) {
//    vector <int> res;
//    std::stringstream ss(str);
//    int i;
//    while (ss >> i) {
//        res.push_back(i);
//        bool done = false;
//        while (!done) { // shouldn't be any spaces, but let's be safe
//            if (ss.peek() == ' ' || ss.peek() == ',') {
//                ss.ignore();
//            } else {
//                done = true;
//            }
//        }
//    }
//    return res;
//}


bool is_number (const string& s) {
    string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) {
        ++it;
    }
    return !s.empty() && it == s.end();
}


// higher resolution than time( NULL );
unsigned int get_clock_seed () {
    return (std::chrono::high_resolution_clock::now().time_since_epoch().count());
}


// works a little differently than above; don't need to trim spaces
// assumes delimiter is some form of whitespace
vector <string> tokenize (string const& input) {
    vector <string> tokens;
    string temp;
    istringstream str(input);
    while (str >> temp) {
        tokens.push_back(temp);
    }
    return tokens;
}


void trim_spaces (string& str) {
    // Trim Both leading and trailing spaces
    size_t startpos = str.find_first_not_of(" \t\r\n"); // Find the first character position after excluding leading blank spaces
    size_t endpos = str.find_last_not_of(" \t\r\n"); // Find the first character position from reverse af

    // if all spaces or empty return an empty string
    if (string::npos == startpos || string::npos == endpos) {
        str = "";
    } else {
        str = str.substr(startpos, endpos - startpos + 1);
    }
    /*
    // Code for Trim Leading Spaces only
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


bool check_comment_line (string const& line) {
    bool comment = false;
    string trimmed = line;
    trim_spaces(trimmed); // get rid of leading spaces
    if (trimmed[0] == '#' || trimmed[0] == '[') {
        comment = true;
    }
    return comment;
}


vector <vector <double> > processRateMatrixConfigFile (string filename, int numstates) {
    vector <double> cols(numstates,1);
    vector <vector <double> > ratematrix = vector <vector <double> > (numstates, cols);
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


int random_int_range (int min, int max) {
    return min + (rand() % (int)(max - min + 1));
}


// given numTotal sites, sample numSample without replacement between 0 -> (numTotal-1)
// ok, this is pretty sweet, if i do say so myself
vector <int> sample_without_replacement (int const& numTotal, int const& numSample) {
    vector <int> randsites (numSample); // numchar zero-initialized elements
    vector <int> allsites (numTotal);
    iota(allsites.begin(), allsites.end(), 0); // generate sequence 0,1,2...n-1
    
    int randNum = 0;
    for (int i = 0; i < numSample; i++) {
        randNum = random_int_range(i, (numTotal - 1));
    // swap, so don't have to worry about multiple hits
        swap(allsites[i], allsites[randNum]);
        randsites[i] = allsites[i];
    }
    return randsites;
}


// arg below is always '?'. besides, getopt prints errors to sterr
void print_error (char * pname, char arg) {
    // cout << pname <<": invalid option -- '" << arg << "'" << endl;
    cout << "Try `" << pname << " --help' for more information." << endl;
}


bool test_logical (vector <int> & matA, vector <int> & matB) {
    bool test = false;
    int match1 = 0;
    unsigned int numdiffs = 0;
    for (unsigned int i=0; i < matA.size(); i++) {
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

/*
 * edgewise added because we have to check the reverse if not rooted
 *
 * TODO : should probably look at this
 */
bool test_logical (vector <int> & matA, vector <int> & matB,bool edgewise) {
    bool test = false;
    int match1 = 0;
    unsigned int numdiffs = 0;
    for (unsigned int i=0; i < matA.size(); i++) {
        if (((matA[i] == 1) && (matB[i] == 1))) {
            match1 += 1;
        } else {
            numdiffs += 1;
        }
    }
    //had to change because of negatives
    if ((match1 != sum(matA)) && (match1 != sum(matB))) {
        if (numdiffs != matA.size()) {
            if(edgewise == true){
                int match2 = 0;
                unsigned int numdiffs2 = 0;
                for (unsigned int i=0; i < matA.size(); i++) {
                    if (((matA[i] == 1) && (matB[i] == 0))) {
                        match2 += 1;
                    } else {
                        numdiffs2 += 1;
                    }
                }
                if ((match2 != sum(matA)) && (match2 != sum_zeros(matB))) {
                    if (numdiffs2 != matA.size()) {
                        test = true;
                    }
                }
            }else{
                test = true;
            }
        }
    }
    return test;
}


//------------------------------------------------------------------------//
// simple math on vectors

int sum_matrix_col (vector <vector <int> > & matrix,int col) {
    int x=0;
    for (unsigned int i=0; i < matrix.size(); i++) {
        x += matrix[i][col];
    }
    return x;
}


int sum_matrix_col_negs (vector <vector <int> > & matrix, int col) {
    int x=0;
    for (unsigned int i=0; i < matrix.size(); i++) {
        if (matrix[i][col] < 0) {
            x += matrix[i][col];
        }
    }
    return x;
}


double mean (vector <double> & in) {
    return sum (in) / (double)in.size();
}


double variance (vector <double> & in) {
    double meann = mean(in);
    
    std::vector<double> diff(in.size());
    std::transform(in.begin(), in.end(), diff.begin(),
        std::bind2nd(std::minus<double>(), meann));
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    
    double var = sq_sum / (double)in.size();
    
    return var;
}


double sum (vector <double> & in) {
    return accumulate(in.begin(), in.end(), 0.0);
}


int sum (vector <int> & in) {
    return accumulate(in.begin(), in.end(), 0);
}


int sum_zeros (vector <int> & in) {
    int x = 0;
    for(unsigned int i=0;i < in.size(); i++){
        if (in[i] == 0)
            x += 1;
    }
    return x;
}


Superdouble calculate_vector_Superdouble_sum (vector <Superdouble> & in) {
    Superdouble sum = 0;
    for (unsigned int i=0; i < in.size(); i++) {
        sum += in[i];
        // cout << in[i] << " sum:" << sum << endl;
    }
    // cout << "endsum:" << sum << endl;
    return sum;
}


// add element i in vec1 to element i in vec2
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


//------------------------------------------------------------------------//


string get_string_vector (vector <string> &sts) {
    string rets;
    for (unsigned int i=0; i < sts.size(); i++) {
        rets += sts[i]+ " ";
    }
    return rets;
}


string get_string_vector (vector <int> & sts) {
    string rets;
    for (unsigned int i=0; i < sts.size(); i++) {
        rets += to_string(sts[i]) + " ";
    }
    return rets;
}


// replace all occurrences of origSubStr to replSubStr
void replace_all (string& str, string const& origSubStr, string const& replSubStr) {
    if (origSubStr.empty()) {
        return;
    }
    size_t start_pos = 0;
    while ((start_pos = str.find(origSubStr, start_pos)) != std::string::npos) {
        str.replace(start_pos, origSubStr.length(), replSubStr);
        start_pos += replSubStr.length();
    }
}


// similar to above, but here we replace any 'bad' characters present in origSubStr
// by the replacement. each bad character will be replaced (i.e. contiguous
// characters will each be replaced). replacement string may be longer than what
// is being replaced (a single character).
// e.g. we might replace each chars in "()[]:;" by "_"
void replace_each (string& str, string const& badChars, string const& replSubStr) {
    if (badChars.empty()) {
        return;
    }
    size_t start_pos = str.find_first_of(badChars);
    while (start_pos != std::string::npos) {
        str.replace(start_pos, 1, replSubStr);
        start_pos = str.find_first_of(badChars, start_pos+replSubStr.length());
    }
}

/*------------------------------------------------------------------------/

 Phun with phormats!
 Both Nexus and newick have rules regarding punctuation and the treatment
 of whitespace characters.
 Of course they are not equivalent (that would be too easy)...
 
 Newick rules are given here: http://evolution.genetics.washington.edu/phylip/newick_doc.html
 Unquoted labels may not contain blanks, parentheses, square brackets,
 single_quotes, colons, semicolons, or commas.
 
 Nexus rules are given in Maddison et al. 1997. Syst. Biol. 46 (4): 590-621 doi: 10.1093/sysbio/46.4.590
 Nexus punctuation characters: () [] {} / \, ; := * ' "` + - <>
 
 Commonalities:
 1) if a token contains punctuation it needs to be surrounded by single quotes
 2) single quotes within a token are replaced by doubled single quotes
    e.g. the token 'John's' should be represented as 'John''s'
 3) underscores in unquoted tokens should be treated as blank spaces ('dark characters')
    e.g. B._zephyrum is equivalent to 'B. zephyrum'
 
 Where they diverge (punctuation in Nexus but not newick):  
 
/-----------------------------------------------------------------------*/
const string nexus_punct  = "()[]{}/\\,;:=*\'\"`+-<>";
const string newick_punct = "()[]\':;,";


// get a taxon label that is newick-compliant
string get_valid_newick_label (string const& inLabel) {
    string outLabel = inLabel;
    
    // if surrounded by single quotes already, assume cool
    if (outLabel[0] == '\'' && outLabel[outLabel.size() - 1] == '\'') {
        return outLabel;
    }
    // does the label require quotes?
    if (outLabel.find_first_of(newick_punct) != std::string::npos) {
        quotify_label(outLabel);
    } else if (outLabel.find_first_of(" ") != std::string::npos) {
        std::replace(outLabel.begin(), outLabel.end(), ' ', '_');
    }
    return outLabel;
}


// get a taxon label that is Nexus-compliant
string get_valid_nexus_label (string const& inLabel) {
    string outLabel = inLabel;
    
    // if surrounded by single quotes already, assume cool
    if (outLabel[0] == '\'' && outLabel[outLabel.size() - 1] == '\'') {
        return outLabel;
    }
    // does the label require quotes?
    if (outLabel.find_first_of(nexus_punct) != std::string::npos) {
        quotify_label(outLabel);
    } else if (outLabel.find_first_of(" ") != std::string::npos) {
        std::replace(outLabel.begin(), outLabel.end(), ' ', '_');
    }
    return outLabel;
}


// alters the label, but *should* open in all programs
// replace all illegal characters with an underscore
string get_safe_taxon_label (string const& inLabel) {
    string outLabel = inLabel;
    replace_each(outLabel, nexus_punct, "_");
    return outLabel;
}


// got an invalid token. replace internal quotes and underscores, surround by quotes
void quotify_label (string & token) {
    // replace internal quotes
    replace_all(token, "'", "''");
    
    // replace underscores
    std::replace(token.begin(), token.end(), '_', ' ');
    
    // add outer quotes
    token = "'" + token + "'";
}


//------------------------------------------------------------------------//

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


double logn (double x, double base) {
    return log10(x)/log10(base);
}


//------------------------------------------------------------------------//
// equality checks //

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


bool check_for_input_to_stream (){
    struct pollfd pfd = { STDIN_FILENO, POLLIN, 0 };
    int ret = 0;
    ret = poll(&pfd, 1, 500);
    if (ret == 0) {
        return false;
    } else {
        return true;
    }
}


// not using right now
// return elements in a *not* found in b
vector <string> get_complement (vector <string> & a, vector <string> & b) {
    vector <string> comp;
    for (unsigned int i=0; i < a.size(); i++) {
        if (find(b.begin(), b.end(), a[i]) == b.end()) {
            comp.push_back(a[i]);
        }
    }
    return comp;
}
