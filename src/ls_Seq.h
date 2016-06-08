#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <map>
#include <sstream>


using namespace std;

#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"

class Stats {
private:

	string Concatenated;
	string temp_seq;
	map <char, double> Nuc_Total;
	map <char,double> AA_Total;

public:
    Stats();
    Stats (istream* pios, bool& all, bool& prot);
    void GC_Getter(string& seq);
    void AA_STAT_Getter(string& seq);
};
