#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <limits>
#include <functional>
#include <ctime>
#include <chrono>
#include <poll.h>
#include <unistd.h>
#include <regex>

#include "utils.h"
#include "superdouble.h"
#include "constants.h"


// TODO: use const where possible

// set threshold. should maybe use superdouble.
//extern const double EPSILON; // moved to constants.h

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
    if (in == out) {
        std::cerr << "Error: input and output streams (files) must differ. Exiting." << std::endl;
        exit(0);
    }
}


// catch exception (atoi does not do this)
// arg is passed for error messages only
int string_to_int (const std::string& in, const std::string& arg) {
    int res = 0;
    try {
            res = std::stoi(in);
    } catch (const std::invalid_argument& ia) {
        std::cerr << "Error: invalid argument for " << arg
                << "; expecting int (" << ia.what() << "). Exiting." << std::endl;
        exit(0);
    }
    return res;
}


long int string_to_long_int (const std::string& in, const std::string& arg) {
    long int res = 0;
    try {
            res = std::stol(in);
    } catch (const std::invalid_argument& ia) {
        std::cerr << "Error: invalid argument for " << arg
                << "; expecting long int (" << ia.what() << "). Exiting." << std::endl;
        exit(0);
    }
    return res;
}


// catch exception (atof does not do this)
float string_to_float (const std::string& in, const std::string& arg) {
    float res = 0;
    try {
            res = std::stof(in);
    } catch (const std::invalid_argument& ia) {
            std::cerr << "Error: invalid argument for " << arg
                    << "; expecting float (" << ia.what() << "). Exiting." << std::endl;
            exit(0);
    }
    return res;
}


double string_to_double (const std::string& in, const std::string& arg) {
    double res = 0;
    try {
            res = std::stod(in);
    } catch (const std::invalid_argument& ia) {
            std::cerr << "Error: invalid argument for " << arg
                    << "; expecting double (" << ia.what() << "). Exiting." << std::endl;
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


void tokenize (const std::string& str, std::vector<std::string>& tokens,
        const std::string& delimiters) {
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
    while (it != s.end() && (std::isdigit(*it) != 0)) {
        ++it;
    }
    return !s.empty() && it == s.end();
}


// get the longest label. for printing/formatting purposes
unsigned int get_longest_label (std::vector<std::string>& labels) {
    unsigned int longest_label_ = 0;
    for (auto & label : labels) {
        auto cur_len = static_cast<unsigned int>(label.size());
        if (cur_len > longest_label_) {
            longest_label_ = cur_len;
        }
    }
    return longest_label_;
}


// properly read in unix, old mac, and windows line endings
// pilfered from: https://stackoverflow.com/a/6089413/3389183
std::istream& getline_safe(std::istream& is, std::string& t) {
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for (;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if (sb->sgetc() == '\n') {
                sb->sbumpc();
            }
            return is;
        case std::streambuf::traits_type::eof():
            // Also handle the case when the last line has no line ending
            if (t.empty()) {
                is.setstate(std::ios::eofbit);
            }
            return is;
        default:
            t += static_cast<char>(c);
        }
    }
}


// not currently used
unsigned int factorial (unsigned int n) {
    return (n == 0 || n == 1) ? 1 : factorial(n - 1) * n;
}


// used for counting trees
unsigned long int doublefactorial (unsigned int n) {
    return (n == 0 || n == 1) ? 1 : doublefactorial(n - 2) * n; 
} 


// higher resolution than time( nullptr );
unsigned int get_clock_seed () {
    return (static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
}


void trim_spaces (std::string& str) {
    // Trim Both leading and trailing spaces
    // Find the first character position after excluding leading blank spaces
    size_t startpos = str.find_first_not_of(" \t\r\n");
    // Find the first character position from reverse af
    size_t endpos = str.find_last_not_of(" \t\r\n");

    // if all spaces or empty return an empty string
    if (std::string::npos == startpos || std::string::npos == endpos) {
        str = "";
    } else {
        str = str.substr(startpos, endpos - startpos + 1);
    }
    /*
    // Code for Trim Leading Spaces only
    // Find the first character position after excluding leading blank spaces
    size_t startpos = str.find_first_not_of(\t);
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


// used by pxstrec
std::vector<std::vector<double> > processRateMatrixConfigFile (const std::string& filename,
        int numstates) {
    auto ns = static_cast<size_t>(numstates);
    std::vector<double> cols(ns, 1);
    std::vector<std::vector<double> > ratematrix;
    ratematrix = std::vector<std::vector<double> > (ns, cols);
    //read file
    std::ifstream ifs(filename.c_str());
    std::string line;
    int fromarea = 0;
    while (getline_safe(ifs, line)) {
        if (line.size() > 3) {
            std::vector<std::string> tokens;
            std::string del(" ,\t");
            tokens.clear();
            tokenize(line, tokens, del);
            for (auto & tk : tokens) {
                trim_spaces(tk); // this will never be used, as it was split on whitespace
            }
            for (unsigned int j = 0; j < tokens.size(); j++) {
                ratematrix[static_cast<size_t>(fromarea)][j] = std::atof(tokens[j].c_str());
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


// NOTE: this assumes that srand have been seeded previously
// e.g., srand(get_clock_seed());
unsigned int random_int_range (const unsigned int& min, const unsigned int& max) {
    return min + (rand() % static_cast<unsigned int>(max - min + 1));
}


// given numTotal sites, sample numSample without replacement between 0 -> (numTotal-1)
// ok, this is pretty sweet, if i do say so myself
std::vector<unsigned int> sample_without_replacement (const unsigned int& numTotal,
        const unsigned int& numSample) {
    std::vector<unsigned int> randsites (static_cast<size_t>(numSample)); // numchar zero-initialized elements
    std::vector<unsigned int> allsites (static_cast<size_t>(numTotal));
    std::iota(allsites.begin(), allsites.end(), 0); // generate sequence 0, 1, 2..., n-1
    
    for (unsigned int i = 0; i < numSample; i++) {
        unsigned int randNum = random_int_range(i, (numTotal - 1u));
    // swap, so don't have to worry about multiple hits
        std::swap(allsites[static_cast<size_t>(i)], allsites[static_cast<size_t>(randNum)]);
        randsites[static_cast<size_t>(i)] = allsites[static_cast<size_t>(i)];
    }
    return randsites;
}


// arg below is always '?'. besides, getopt prints errors to sterr
void print_error (char * pname) {
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
            if (edgewise) {
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

int sum_matrix_col (std::vector<std::vector<int> >& matrix, unsigned int col) {
    int x = 0;
    for (auto & mati : matrix) {
        x += mati[static_cast<size_t>(col)];
    }
    return x;
}


// not used
int sum_matrix_col_negs (std::vector<std::vector<int> >& matrix, int col) {
    int x = 0;
    for (auto & mati : matrix) {
        if (mati[static_cast<size_t>(col)] < 0) {
            x += mati[static_cast<size_t>(col)];
        }
    }
    return x;
}


// NOTE: this partially rearranges elements, so make a copy if order is important
double v_median (std::vector<double>& in) {
    int n = static_cast<int>(in.size()) / 2;
    std::nth_element(in.begin(), in.begin()+n, in.end());
    double median = in[static_cast<size_t>(n)];
    
    if (n % 2 == 0) { // if the length is even
      auto max_it = std::max_element(in.begin(), in.begin()+n);
      median = (*max_it + median) / 2.0;
    }
    return median;    
}


// should make this a templated function
double v_mean (std::vector<double>& in) {
    return sum(in) / static_cast<double>(in.size());
}


// if you want one, prolly want the other too
void v_mean_variance (std::vector<double>& in, double& mn, double& varr) {
    auto n = static_cast<double>(in.size());
    mn = sum(in) / n;
    std::vector<double> diff(in.size());
    std::transform(in.begin(), in.end(), diff.begin(),
        std::bind(std::minus<double>(), std::placeholders::_1, mn));
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    varr = sq_sum / n;
}


double v_variance (std::vector<double>& in) {
    auto n = static_cast<double>(in.size());
    double meann = sum(in) / n;
    std::vector<double> diff(static_cast<size_t>(n));
    std::transform(in.begin(), in.end(), diff.begin(),
        std::bind(std::minus<double>(), std::placeholders::_1, meann));
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double var = sq_sum / n;
    return var;
}


double sum (std::vector<double>& in) {
    return std::accumulate(in.begin(), in.end(), 0.0);
}


int sum (std::vector<int>& in) {
    return std::accumulate(in.begin(), in.end(), 0);
}


// hrm look into templated function
unsigned int sum (std::vector<unsigned int>& in) {
    return static_cast<unsigned int>(std::accumulate(in.begin(), in.end(), 0));
}


std::vector<double> string_v_to_double_v (const std::vector<std::string>& in) {
    std::vector<double> out(in.size());
    std::transform(in.begin(), in.end(), out.begin(), [](const std::string& val) {return std::stod(val);});
    return out;
}


int count_zeros (std::vector<int>& in) {
    int x = 0;
    for (int i : in) {
        if (i == 0) {
            x += 1;
        }
    }
    return x;
}


Superdouble calculate_vector_Superdouble_sum (std::vector<Superdouble>& in) {
    Superdouble sum(0);
    for (const auto & i : in) {
        sum += i;
        // std::cout << i << " sum:" << sum << std::endl;
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
    std::transform(res.begin(), res.end(), vec2.begin(), res.begin(), std::plus<>());
    return res;
}


std::vector<double> average_vectors_elementwise (std::vector<double>& vec1,
        std::vector<double>& vec2) {
    // bail if sequences are of different lengths.
    if (vec1.size() != vec2.size()) {
      throw std::invalid_argument(
          "Vectors must be of equal length"
      );
    }
    std::vector<double> res = sum_vectors_elementwise(vec1, vec2);
    std::transform(res.begin(), res.end(), res.begin(), std::bind(std::multiplies<double>(), 0.5, std::placeholders::_1));
    
    return res;
}

//------------------------------------------------------------------------//


std::string get_string_vector(std::vector<std::string>& sts) {
    std::string rets;
    for (const auto & st : sts) {
        rets += st + " ";
    }
    return rets;
}


std::string get_string_vector(std::vector<int>& sts) {
    std::string rets;
    for (int st : sts) {
        rets += std::to_string(st) + " ";
    }
    return rets;
}


// replace all occurrences of origSubStr to replSubStr
void replace_all (std::string& str, const std::string& origSubStr,
        const std::string& replSubStr) {
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
void replace_each (std::string& str, const std::string& badChars,
        const std::string& replSubStr) {
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
//const std::string nexus_punct  = "()[]{}/\\,;:=*\'\"`+-<>"; // escaped string literal
const std::string nexus_punct  = R"(()[]{}/\,;:=*'"`+-<>)"; // raw string literal
const std::string newick_punct = "()[]\':;,";

// given a line that begins with [, keep reading until a line where last char is ] (i.e., end of comment)
// occurs in BOTH tree and alignment files (hence, location)
// NOTE: returns the end comment line (so reader will need to read in next valid line)
void process_nexus_comment (std::istream& stri, std::string& tline) {
    std::string terp = tline;
    trim_spaces(terp);
    // check single-line comment
    if (terp.back() == ']') {
        //std::cout << "single-line comment!" << std::endl;
        return;
    }
    // if not, dealing with a multi-line comment
    while (true) {
        getline_safe(stri, terp);
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
// this is not currently used, but probably should be
std::string get_valid_newick_label (const std::string& inLabel) {
    std::string outLabel = inLabel;
    
    // if surrounded by single quotes already, assume cool
    if (outLabel[0] == '\'' && outLabel[outLabel.size() - 1] == '\'') {
        return outLabel;
    }
    // does the label require quotes?
    if (outLabel.find_first_of(newick_punct) != std::string::npos) {
        quotify_label(outLabel);
    } else if (outLabel.find_first_of(' ') != std::string::npos) {
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
    } else if (outLabel.find_first_of(' ') != std::string::npos) {
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
          "Hamming distances are only defined for strings of equal length."
      );
    }
    
    unsigned int startv = 0;
    
    return std::inner_product(
        s1.begin(), s1.end(), s2.begin(), startv, std::plus<>(),
        std::not2(std::equal_to<std::string::value_type>())
    );
}


double logn (double x, double base) {
    return log10(x)/log10(base);
}


//------------------------------------------------------------------------//
// equality checks //

// check if 2 doubles are equal within some tolerance.
bool essentially_equal_doubles (double a, double b) {
    
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


bool all_equal (std::vector<double> vals) {
    bool equal = false;
    auto it = std::find_if_not(vals.begin()+1, vals.end(),
            std::bind(essentially_equal_doubles, std::placeholders::_1, vals[0]));
    if (it == end(vals)) {
        equal = true;
    }
    return equal;
}


bool all_equal (std::vector<int> vals) {
    bool equal = false;
    if (std::all_of(vals.begin(), vals.end(), [&] (int i) {return i == vals[0];})) {
        equal = true;
    }
    return equal;
}


bool all_equal (std::vector<unsigned int> vals) {
    bool equal = false;
    if (std::all_of(vals.begin(), vals.end(), [&] (unsigned int i) {return i == vals[0];})) {
        equal = true;
    }
    return equal;
}


bool check_for_input_to_stream () {
    bool ret = true;
    struct pollfd pfd = { STDIN_FILENO, POLLIN, 0 };
    int ret_poll = 0;
    ret_poll = poll(&pfd, 1, 500);
    if (ret_poll == 0) {
        ret = false;
    }
    return ret;
}


// do strings statr/end with suffix?!?
bool ends_with (const std::string& str, const std::string& suffix) {
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}


bool starts_with (const std::string& str, const std::string& prefix) {
    return str.size() >= prefix.size() && 0 == str.compare(0, prefix.size(), prefix);
}


void remove_last_N (std::string &str, const long unsigned int& n) {
    if (str.length() < n) {
        return;
    }
    str.erase(str.length() - n);
}

// not using right now
// return elements in a *not* found in b
std::vector<std::string> get_complement (std::vector<std::string>& a, std::vector<std::string>& b) {
    std::vector<std::string> comp;
    for (const auto & i : a) {
        if (find(b.begin(), b.end(), i) == b.end()) {
            comp.push_back(i);
        }
    }
    return comp;
}


// peek at the next line in the stream
std::string peek_line (std::istream& pios) {
    std::string nextLine;
    // get current position
    std::streampos spt = pios.tellg();
    // read in next line
    getline_safe(pios, nextLine);
    // return to the original position in the stream
    pios.seekg(spt, std::ios_base::beg);
    return nextLine;
}


// same as above, except read in the next n lines
std::vector<std::string> peek_lines (std::istream& pios, const int& n) {
    std::vector<std::string> peekedLines;
    std::string nextLine;
    // get current position
    std::streampos spt = pios.tellg();
    // read in the lines
    for (int i = 0; i < n; i++) {
        getline_safe(pios, nextLine);
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
    for (const auto & name : names) {
        //std::cout << name << ": " << std::regex_search(name, regexp) << std::endl;
        if (std::regex_search(name, regexp)) {
            results.push_back(name);
        }
    }
    return results;
}


// not currently used
std::vector<std::string> regex_replace_labels (const std::vector<std::string>& names,
        const std::string& pattern, const std::string& replacetext) {
    std::vector<std::string> results;
    const std::regex regexp(pattern);
    for (const auto & name : names) {
        std::string res = std::regex_replace(name, regexp, replacetext);
        //std::cout << name << ": " << std::regex_search(name, regexp) << std::endl;
    }
    return results;
}
