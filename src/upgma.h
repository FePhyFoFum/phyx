/*
 * upgma.h
 *
 *  Created on: Jun 10, 2015
 *      Author: joe
 */

#ifndef _UPGMA_H_
#define _UPGMA_H_

// put includes in source files
/*
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
*/

class UPGMA {
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
    UPGMA();
    UPGMA (string & fasta);
    map<string, string> FastaToOneLine(string& fasta);
    vector< vector<double> > BuildMatrix(map<string, string>& sequences);
    void TREEMAKE(vector<string>&, map <int, string>&, vector< vector<double> >&);
    string get_newick ();
    //virtual ~UPGMA();
};

#endif /* _UPGMA_H_ */
