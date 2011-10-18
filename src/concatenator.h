#ifndef CONCAT_H_
#define CONCAT_H_

#include <string>
#include <vector>

#include "Sequence.h"

using namespace std;

class Concatenator{
private:
    vector<string> & filenames;
    vector<Sequence> concatenated_seqs;
    vector<int> gene_lengths;
    
public:
    Concatenator();
    Concatenator(const vector<string> &filenames);
    run();
}; 

#endif
