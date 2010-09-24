#include <iostream>
#include <vector>
#include <map>
#include <string>

#include "utils.h"

using namespace std;

void Tokenize(const string& str, vector<string>& tokens,
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

void TrimSpaces( string& str)  {
	// Trim Both leading and trailing spaces
	size_t startpos = str.find_first_not_of(" \t\r\n"); // Find the first character position after excluding leading blank spaces
	size_t endpos = str.find_last_not_of(" \t\r\n"); // Find the first character position from reverse af

	// if all spaces or empty return an empty string
	if(( string::npos == startpos ) || ( string::npos == endpos))
	{
		str = "";
	}
	else
		str = str.substr( startpos, endpos-startpos+1 );

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

double calculate_vector_double_sum(vector<double> & in){
	double sum = 0;
	for (unsigned int i=0;i<in.size();i++){
		sum += in[i];
	}
	return sum;
}
