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
    string Molecule;
    string type;
    bool finished;
    string Mol;
    string name;
    map <char, double> Total;
    int seqcount;


public:
    Stats();
    Stats (istream* pios, bool& all, bool& prot);
    void STAT_Getter(string& seq, bool& prot);
    void Printer(bool& prot);
};
