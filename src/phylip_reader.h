#ifndef PHYLIP_READER_H_
#define PHYLIP_READER_H_

#include <string>
#include <vector>
using namespace std;

class PhylipReader {
private:

public:
	PhylipReader();
	bool readFile(string,vector<Sequence>&);
};

#endif
