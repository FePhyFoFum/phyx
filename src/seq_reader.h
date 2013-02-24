#ifndef SEQ_READER_H_
#define SEQ_READER_H_

#include <string>
#include <vector>
#include "sequence.h"

using namespace std;

bool read_fasta_file(string,vector<Sequence>&);
bool read_phylip_file(string,vector<Sequence>&);

#endif
