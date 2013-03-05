#ifndef SEQ_MODELS_H
#define SEQ_MODELS_H

#include <vector>
#include <map>
#include <set>
#include <string>

using namespace std;

#include "sequence.h"

void read_scoring_matrix(char * filename, map<char, map<char,int> > & sc_mat);
void read_scoring_matrix_from_lines(vector<string> & lines, map<char,map<char,int> > & sc_mat);
void get_ednafull(map<char, map<char,int> > & sc_mat);
void get_blosum62(map<char, map<char,int> > & sc_mat);

#endif
