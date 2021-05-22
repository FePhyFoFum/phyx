#ifndef PX_SEQ_MODELS_H
#define PX_SEQ_MODELS_H

#include <string>
#include <vector>
#include <map>

#include "sequence.h"


void read_scoring_matrix (char * filename, std::map<char, std::map<char, int> >& sc_mat);
void read_scoring_matrix_from_lines (std::vector<std::string>& lines,
    std::map<char, std::map<char, int> >& sc_mat);
void get_ednafull (std::map<char, std::map<char, int> >& sc_mat);
void get_blosum62 (std::map<char, std::map<char, int> >& sc_mat);

#endif /* PX_SEQ_MODELS_H */
