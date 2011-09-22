#ifndef FASTA_UTIL_H_
#define FASTA_UTIL_H_

#include <string>

using namespace std;

class FastaUtil{
private:

public:
	FastaUtil();
	bool readFile(string,vector<Sequence>&);
	bool writeFileFromVector(string, vector<Sequence>&);
};

#endif
