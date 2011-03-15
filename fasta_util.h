#ifndef FASTA_READER_H_
#define FASTA_READER_H_

#include <string>

using namespace std;

class FastaReader{
private:

public:
	FastaReader();
	bool readFile(string,vector<Sequence>&);
};

#endif
