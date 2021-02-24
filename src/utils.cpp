#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <limits>
#include <functional>
#include <ctime>
#include <chrono>
#include <poll.h>
#include <unistd.h>

#include "utils.h"
#include "superdouble.h"


// TODO: use const where possible

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


void check_file_exists (const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.good()) {
        std::cerr << "Error: file '" << filename << "' not found. Exiting." << std::endl;
        exit(0);
    }
}


void check_inout_streams_identical (char * in, char * out) {
    if ((std::string)in == (std::string)out) {
        std::cerr << "Error: input and output file names must differ (streams!). Exiting." << std::endl;
        exit(0);
    }
}


// catch exception (atoi does not do this)
int string_to_int (const std::string& in, const std::string& arg) {
    int res = 0;
    try {
            res = stoi(in);
        }
        catch (const std::invalid_argument& ia) {
            std::cerr << "Error: invalid argument for " << arg << " (expecting int). Exiting." << std::endl;
            exit(0);
        }
    return res;
}


// catch exception (atof does not do this)
float string_to_float (const std::string& in, const std::string& arg) {
    float res = 0;
    try {
            res = stof(in);
        }
        catch (const std::invalid_argument& ia) {
            std::cerr << "Error: invalid argument for " << arg << " (expecting float). Exiting." << std::endl;
            exit(0);
        }
    return res;
}


// convenience func. returns copy so original can still be used
std::string string_to_upper (const std::string& instr) {
    std::string outstr = instr;
    std::transform(outstr.begin(), outstr.end(), outstr.begin(), ::toupper);
    return outstr;
}


// and the other one
std::string string_to_lower (const std::string& instr) {
    std::string outstr = instr;
    std::transform(outstr.begin(), outstr.end(), outstr.begin(), ::tolower);
    return outstr;
}


void tokenize (const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters) {
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}


// works a little differently than above; don't need to trim spaces
// assumes delimiter is some form of whitespace
std::vector<std::string> tokenize (const std::string& input) {
    std::vector<std::string> tokens;
    std::string temp;
    std::istringstream str(input);
    while (str >> temp) {
        tokens.push_back(temp);
    }
    return tokens;
}


bool is_number (const std::string& s) {
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) {
        ++it;
    }
    return !s.empty() && it == s.end();
}


// get the longest label. for printing/formatting purposes
int get_longest_label (std::vector<std::string>& labels) {
    int longest_label_ = 0;
    for (unsigned int i = 0; i < labels.size(); i++) {
        if ((int)labels[i].size() > longest_label_) {
            longest_label_ = labels[i].size();
        }
    }
    return longest_label_;
}


int factorial (int n) {
    return (n == 0 || n == 1) ? 1 : factorial(n - 1) * n;
}


// used for counting trees
int doublefactorial(int n) {
    return (n == 0 || n == 1) ? 1 : doublefactorial(n - 2) * n; 
} 


// higher resolution than time( NULL );
unsigned int get_clock_seed () {
    return (std::chrono::high_resolution_clock::now().time_since_epoch().count());
}


void trim_spaces (std::string& str) {
    // Trim Both leading and trailing spaces
    size_t startpos = str.find_first_not_of(" \t\r\n"); // Find the first character position after excluding leading blank spaces
    size_t endpos = str.find_last_not_of(" \t\r\n"); // Find the first character position from reverse af

    // if all spaces or empty return an empty string
    if (std::string::npos == startpos || std::string::npos == endpos) {
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


bool check_comment_line (const std::string& line) {
    bool comment = false;
    std::string trimmed = line;
    trim_spaces(trimmed); // get rid of leading spaces
    if (trimmed[0] == '#' || trimmed[0] == '[') {
        comment = true;
    }
    return comment;
}


std::vector<std::vector<double> > processRateMatrixConfigFile (std::string filename, int numstates) {
    std::vector<double> cols(numstates, 1);
    std::vector<std::vector<double> > ratematrix = std::vector<std::vector<double> > (numstates, cols);
    //read file
    std::ifstream ifs(filename.c_str());
    std::string line;
    int fromarea = 0;
    while (getline(ifs, line)) {
        if (line.size() > 3) {
            std::vector<std::string> tokens;
            std::string del(" ,\t");
            tokens.clear();
            tokenize(line, tokens, del);
            for (unsigned int j = 0; j < tokens.size(); j++) {
                trim_spaces(tokens[j]);
            }
            for (unsigned int j = 0; j < tokens.size(); j++) {
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
std::vector<int> sample_without_replacement (const int& numTotal, const int& numSample) {
    std::vector<int> randsites (numSample); // numchar zero-initialized elements
    std::vector<int> allsites (numTotal);
    std::iota(allsites.begin(), allsites.end(), 0); // generate sequence 0, 1, 2..., n-1
    
    int randNum = 0;
    for (int i = 0; i < numSample; i++) {
        randNum = random_int_range(i, (numTotal - 1));
    // swap, so don't have to worry about multiple hits
        std::swap(allsites[i], allsites[randNum]);
        randsites[i] = allsites[i];
    }
    return randsites;
}


// arg below is always '?'. besides, getopt prints errors to sterr
void print_error (char * pname, char arg) {
    // std::cout << pname <<": invalid option -- '" << arg << "'" << std::endl;
    std::cout << "Try `" << pname << " --help' for more information." << std::endl;
}


bool test_logical (std::vector<int>& matA, std::vector<int>& matB) {
    bool test = false;
    int match1 = 0;
    unsigned int numdiffs = 0;
    for (unsigned int i = 0; i < matA.size(); i++) {
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
bool test_logical (std::vector<int>& matA, std::vector<int>& matB, bool edgewise) {
    bool test = false;
    int match1 = 0;
    unsigned int numdiffs = 0;
    for (unsigned int i = 0; i < matA.size(); i++) {
        if (((matA[i] == 1) && (matB[i] == 1))) {
            match1 += 1;
        } else {
            numdiffs += 1;
        }
    }
    //had to change because of negatives
    if ((match1 != sum(matA)) && (match1 != sum(matB))) {
        if (numdiffs != matA.size()) {
            if (edgewise == true) {
                int match2 = 0;
                unsigned int numdiffs2 = 0;
                for (unsigned int i = 0; i < matA.size(); i++) {
                    if (((matA[i] == 1) && (matB[i] == 0))) {
                        match2 += 1;
                    } else {
                        numdiffs2 += 1;
                    }
                }
                if ((match2 != sum(matA)) && (match2 != count_zeros(matB))) {
                    if (numdiffs2 != matA.size()) {
                        test = true;
                    }
                }
            } else {
                test = true;
            }
        }
    }
    return test;
}


//------------------------------------------------------------------------//
// simple math on vectors

int sum_matrix_col (std::vector<std::vector<int> >& matrix, int col) {
    int x=0;
    for (unsigned int i = 0; i < matrix.size(); i++) {
        x += matrix[i][col];
    }
    return x;
}


int sum_matrix_col_negs (std::vector<std::vector<int> >& matrix, int col) {
    int x=0;
    for (unsigned int i = 0; i < matrix.size(); i++) {
        if (matrix[i][col] < 0) {
            x += matrix[i][col];
        }
    }
    return x;
}


// should make this a templated function
double v_mean (std::vector<double>& in) {
    return sum (in) / (double)in.size();
}


double v_variance (std::vector<double>& in) {
    double meann = v_mean(in);
    
    std::vector<double> diff(in.size());
    std::transform(in.begin(), in.end(), diff.begin(),
        std::bind2nd(std::minus<double>(), meann));
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    
    double var = sq_sum / (double)in.size();
    
    return var;
}


double sum (std::vector<double>& in) {
    return std::accumulate(in.begin(), in.end(), 0.0);
}


int sum (std::vector<int>& in) {
    return std::accumulate(in.begin(), in.end(), 0);
}


std::vector<double> string_v_to_double_v (const std::vector<std::string>& in) {
    std::vector<double> out(in.size());
    std::transform(in.begin(), in.end(), out.begin(), [](const std::string& val) {return std::stod(val);});
    return out;
}


int count_zeros (std::vector<int>& in) {
    int x = 0;
    for (unsigned int i = 0; i < in.size(); i++) {
        if (in[i] == 0) {
            x += 1;
        }
    }
    return x;
}


Superdouble calculate_vector_Superdouble_sum (std::vector<Superdouble>& in) {
    Superdouble sum = 0;
    for (unsigned int i = 0; i < in.size(); i++) {
        sum += in[i];
        // std::cout << in[i] << " sum:" << sum << std::endl;
    }
    // std::cout << "endsum:" << sum << std::endl;
    return sum;
}


// add element i in vec1 to element i in vec2
std::vector<int> sum (std::vector<int>& vec1, std::vector<int>& vec2) {
    // bail if sequences are of different lengths.
    if (vec1.size() != vec2.size()) {
      throw std::invalid_argument(
          "Vectors must be of equal length"
      );
    }
    std::vector<int> res = vec1;
    std::transform(res.begin(), res.end(), vec2.begin(), res.begin(), std::plus<int>());
    return res;
}


std::vector<double> average_vectors_elementwise (std::vector<double>& vec1, std::vector<double>& vec2) {
    // bail if sequences are of different lengths.
    if (vec1.size() != vec2.size()) {
      throw std::invalid_argument(
          "Vectors must be of equal length"
      );
    }
    std::vector<double> res = sum_vectors_elementwise(vec1, vec2);
    std::transform(res.begin(), res.end(), res.begin(), std::bind1st (std::multiplies <double> () , 0.5));
    
    return res;
}

//------------------------------------------------------------------------//


std::string get_string_vector(std::vector<std::string>& sts) {
    std::string rets;
    for (unsigned int i = 0; i < sts.size(); i++) {
        rets += sts[i]+ " ";
    }
    return rets;
}


std::string get_string_vector(std::vector<int>& sts) {
    std::string rets;
    for (unsigned int i = 0; i < sts.size(); i++) {
        rets += std::to_string(sts[i]) + " ";
    }
    return rets;
}


// replace all occurrences of origSubStr to replSubStr
void replace_all (std::string& str, const std::string& origSubStr, const std::string& replSubStr) {
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
void replace_each (std::string& str, const std::string& badChars, const std::string& replSubStr) {
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
const std::string nexus_punct  = "()[]{}/\\,;:=*\'\"`+-<>";
const std::string newick_punct = "()[]\':;,";

// given a line that begins with [, keep reading until a line where last char is ] (i.e., end of comment)
// occurs in BOTH tree and alignment files (hence, location)
// NOTE: returns the end comment line (so reader will need to read in next valid line)
void process_nexus_comment (std::istream& stri, std::string& tline) {
    bool done = false;
    std::string terp = tline;
    trim_spaces(terp);
    // check single-line comment
    if (terp.back() == ']') {
        //std::cout << "single-line comment!" << std::endl;
        return;
    }
    // if not, dealing with a multi-line comment
    while (!done) {
        getline(stri, terp);
        trim_spaces(terp);
        if (!terp.empty()) {
            if (terp.back() == ']') {
                //std::cout << "found end of multi-line comment" << std::endl;
                return;
            }
        }
    }
}


bool check_nexus_comment (std::string line) {
    bool comment = false;
    trim_spaces(line);
    if (line[0] == '[') {
        comment = true;
    }
    return comment;
}

// get a taxon label that is newick-compliant
std::string get_valid_newick_label (const std::string& inLabel) {
    std::string outLabel = inLabel;
    
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
std::string get_valid_nexus_label (const std::string& inLabel) {
    std::string outLabel = inLabel;
    
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
std::string get_safe_taxon_label (const std::string& inLabel) {
    std::string outLabel = inLabel;
    replace_each(outLabel, nexus_punct, "_");
    return outLabel;
}


// got an invalid token. replace internal quotes and underscores, surround by quotes
void quotify_label (std::string& token) {
    // replace internal quotes
    replace_all(token, "'", "''");
    
    // replace underscores
    std::replace(token.begin(), token.end(), '_', ' ');
    
    // add outer quotes
    token = "'" + token + "'";
}


//------------------------------------------------------------------------//

// Hamming (edit) distance
unsigned int calc_hamming_dist (const std::string& s1, const std::string& s2) {
    
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
    std::cout << "fabs(a - b) = " << fabs(a - b) << std::endl;
    std::cout << "ApproximatelyEqual: " <<  (fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * EPSILON)) << std::endl;
    std::cout << "EssentiallyEqual: " <<  (fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * EPSILON)) << std::endl;
    */
    
    if (fabs(a - b) <= std::max(EPSILON, EPSILON * std::max(abs(a), abs(b)))) {
        equal = true;
    }
    
    equal = (fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * EPSILON));
    
    return equal;
}


bool all_equal (std::vector<double> vals) {
    bool equal = false;
    std::vector<double>::iterator it = std::find_if_not(vals.begin()+1, vals.end(),
            std::bind(essentially_equal, std::placeholders::_1, vals[0]));
    if (it == end(vals)) {
        equal = true;
    }
    return equal;
}


bool check_for_input_to_stream () {
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
std::vector<std::string> get_complement (std::vector<std::string>& a, std::vector<std::string>& b) {
    std::vector<std::string> comp;
    for (unsigned int i = 0; i < a.size(); i++) {
        if (find(b.begin(), b.end(), a[i]) == b.end()) {
            comp.push_back(a[i]);
        }
    }
    return comp;
}


// peek at the next line in the stream
std::string peek_line (std::istream& pios) {
    std::string nextLine = "";
    // get current position
    std::streampos spt = pios.tellg();
    // read in next line
    getline(pios, nextLine);
    // return to the original position in the stream
    pios.seekg(spt, std::ios_base::beg);
    return nextLine;
}


// same as above, except read in the next n lines
std::vector<std::string> peek_lines (std::istream& pios, const int& n) {
    std::vector<std::string> peekedLines;
    std::string nextLine = "";
    // get current position
    std::streampos spt = pios.tellg();
    // read in the lines
    for (int i = 0; i < n; i++) {
        getline(pios, nextLine);
        peekedLines.push_back(nextLine);
    }
    // return to the original position in the stream
    pios.seekg(spt, std::ios_base::beg);
    return peekedLines;
}


// given a list of names and a regex pattern, return the list of names that contain the pattern
std::vector<std::string> regex_search_labels (const std::vector<std::string>& names,
        const std::string& pattern) {
    std::vector<std::string> results;
    const std::regex regexp(pattern);
    for (int i = 0; i < (int)names.size(); i++) {
        //std::cout << names[i] << ": " << std::regex_search(names[i], regexp) << std::endl;
        if (std::regex_search(names[i], regexp)) {
            results.push_back(names[i]);
        }
    }
    return results;
}
