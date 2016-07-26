/*
 * upgma.h
 *
 *  Created on: Jun 10, 2015
 *      Author: joe
 */

#ifndef _UPGMA_H_
#define _UPGMA_H_

#include <vector>
#include <map>

using namespace std;

class UPGMA {
private:

    map <string, string> sequences;
    map <string, string>::iterator iter;
    vector<string> names;
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
    UPGMA (istream* pios);
    map<string, string> FastaToOneLine(ifstream&);
    vector< vector<double> > BuildMatrix(map<string, string>& sequences);
    void TREEMAKE(vector<string>&, map <int, string>&, vector< vector<double> >&);
    string get_newick ();
};

#endif /* _UPGMA_H_ */
