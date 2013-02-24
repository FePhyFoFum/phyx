#ifndef FASTA_UTIL_H_
#define FASTA_UTIL_H_

#include <string>
#include "sequence.h"

using namespace std;

bool read_fasta_file(string,vector<Sequence>&);
bool write_fasta_file(string, vector<Sequence>&);

#endif
