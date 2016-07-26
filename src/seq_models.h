#ifndef _SEQ_MODELS_H_
#define _SEQ_MODELS_H_

#include <map>

using namespace std;

#include "sequence.h"

void read_scoring_matrix(char * filename, map<char, map<char,int> > & sc_mat);
void read_scoring_matrix_from_lines(vector<string> & lines, map<char,map<char,int> > & sc_mat);
void get_ednafull(map<char, map<char,int> > & sc_mat);
void get_blosum62(map<char, map<char,int> > & sc_mat);

#endif /* _SEQ_MODELS_H_ */
