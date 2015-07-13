/*
 * aatocdn.h
 *
 *  Created on: Jun 15, 2015
 *      Author: joe
 */

#ifndef _AATOCDN_H_
#define _AATOCDN_H_

//#include <iostream>
//#include <string>
//#include <fstream>
//#include <vector>
//#include <sstream>
#include <iterator>
//#include <algorithm>
//#include <map>
#include <iterator>

//using namespace std;

class AAtoCDN {
private:

    //map <string, string> CodonAln; // not used
    map <string, string> CodonSequences;
    map <string, string>::iterator iter;
    string AminoAcidSequence;
    string NucleotideSequence;
    //string temp; // not used
    //string::iterator it; // not used
    //int j; // not used

public:
    AAtoCDN();
    map<string, string> FastaToOneLine (string&);
    map<string, string> ChangeToCodon (map <string, string>&, map <string,string>&);

};

#endif /* _AATOCDN_H_ */
