#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "superdouble.h"

using namespace std;

void tokenize(const string& str, vector<string>& tokens,
	      const string& delimiters){
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

void trim_spaces( string& str)  {
    // Trim Both leading and trailing spaces
    size_t startpos = str.find_first_not_of(" \t\r\n"); // Find the first character position after excluding leading blank spaces
    size_t endpos = str.find_last_not_of(" \t\r\n"); // Find the first character position from reverse af

    // if all spaces or empty return an empty string
    if(( string::npos == startpos ) || ( string::npos == endpos)){
	str = "";
    }
    else{
	str = str.substr( startpos, endpos-startpos+1 );
    }
    /*
    // Code for  Trim Leading Spaces only
    size_t startpos = str.find_first_not_of(\t); // Find the first character position after excluding leading blank spaces
    if( string::npos != startpos )
    str = str.substr( startpos );
    */

    /*
    // Code for Trim trailing Spaces only
    size_t endpos = str.find_last_not_of(\t); // Find the first character position from reverse af
    if( string::npos != endpos )
    str = str.substr( 0, endpos+1 );
    */
}

bool is_number(const string& s){
    string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)){ 
	++it;
    }
    return !s.empty() && it == s.end();
}

double calculate_vector_double_sum(vector<double> & in){
    double sum = 0;
    for (unsigned int i=0;i<in.size();i++){
	sum += in[i];
    }
    return sum;
}

double calculate_vector_double_mean(vector<double> & in){
    double sum = 0;
    for (unsigned int i=0;i<in.size();i++){
	sum += in[i]/in.size();
    }
    return sum;
}

int calculate_vector_int_sum(vector<int> & in){
    int sum = 0;
    for (unsigned int i=0;i<in.size();i++){
	sum += in[i];
    }
    return sum;
}

vector<vector<double> > processRateMatrixConfigFile(string filename, int numstates){
    vector<double> cols(numstates,1);
    vector<vector<double> >  ratematrix = vector<vector<double> > (numstates,cols);
    //read file
    ifstream ifs(filename.c_str());
    string line;
    int fromarea = 0;
    while(getline(ifs,line)){
	if(line.size() > 3){
	    vector<string> tokens;
	    string del(" ,\t");
	    tokens.clear();
	    tokenize(line, tokens, del);
	    for(unsigned int j=0;j<tokens.size();j++){
		trim_spaces(tokens[j]);
	    }
	    for(unsigned int j=0;j<tokens.size();j++){
		ratematrix[fromarea][j] = atof(tokens[j].c_str());
	    }
	    if(fromarea < numstates-1){
		fromarea += 1;
	    }else{
		fromarea = 0;
	    }
	}
    }
    ifs.close();
    return ratematrix;
}

Superdouble calculate_vector_Superdouble_sum(vector<Superdouble> & in){
    Superdouble sum = 0;
    for (unsigned int i=0;i<in.size();i++){
	sum += in[i];
	//cout << in[i] << " sum:" << sum << endl;
    }
    //cout << "endsum:" << sum << endl;
    return sum;
}

int random_int_range(int min, int max){
    return min + (rand() % (int)(max - min + 1));
}

void print_error(char * pname, char arg){
    cout << pname <<": invalid option -- '"<< arg << "'" << endl;
    cout << "Try `"<<pname<<" --help' for more information." << endl;
}
