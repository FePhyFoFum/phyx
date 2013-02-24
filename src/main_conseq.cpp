
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

#include "seq_reader.h"
#include "sequence.h"
#include "seq_utils.h"

int main(int argc, char * argv[]){
    
    if (argc != 3){
	cout << "usage: pxconseq infile.fasta outfile.fasta" << endl;
	exit(0);
    }

    vector<Sequence> seqs;
    bool ret = read_fasta_file(argv[1],seqs);
    if (ret == false){
	cout << argv[1] << " is not a valid fasta file " << endl;
	exit(0);
    }
    cout << seqs.size() << " sequences read" << endl;
    string rets = consensus_seq(seqs,0);
    ofstream ofs;
    ofs.open(argv[2],ios::out);
    ofs << ">" << argv[1] << endl;
    ofs << rets << endl;
    ofs.close();
}
