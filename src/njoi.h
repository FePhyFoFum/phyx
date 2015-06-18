/*
 * NJOI.h
 *
 *  Created on: Jun 12, 2015
 *      Author: joe
 */

//Compile g++ -std=c++11 NJOI.cpp main_NJOI.cpp utils.cpp superdouble.cpp -o test
#ifndef NJOI_H_
#define NJOI_H_

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <map>
#include <iterator>

using namespace std;
class NJOI {
private:

    map <string, string> sequences;
    //map <string, int> NumbKey; // not used
    map <string, string>::iterator iter; // not used
    vector<string> names;
    //string fasta;
    vector< vector<double> > Matrix;
    void set_name_key ();

    // additions:
    int ntax;
    int nchar;
    string seqfile;
    map<int, string> NameKey;
    string newickstring; // temporary

public:
    NJOI();
    NJOI (string & fasta);
    map<string, string> FastaToOneLine(string& fasta);
    vector< vector<double> > BuildMatrix(map<string, string>& sequences);
    void TREEMAKE(vector<string>&, map <int, string>&, vector< vector<double> >&);
    string get_newick ();
};

#endif /* NJOI_H_ */
