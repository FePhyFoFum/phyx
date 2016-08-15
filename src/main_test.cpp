/*
 * main_TEST.cpp
 *
 */

/*
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
*/

#include <cstdlib>
#include <iostream>

// g++ -std=c++11 main_test.cpp -o JWB_test

//g++ -std=c++11 nj.cpp main_nj.cpp utils.cpp superdouble.cpp sequence.cpp seq_reader.cpp seq_utils.cpp log.cpp -o test


using namespace std;

#include "vcf_reader.h"

/*
#include "node.h"
#include "tree_reader.h"
#include "string_node_object.h"
#include "vector_node_object.h"
#include "tree.h"
#include "utils.h"
#include "sequence.h"
#include "fasta_util.h"
#include "log.h"
*/

int main(int argc, char * argv[]){
    
    //log_call(argc, argv);
    
    /*
    TreeReader tr;

    if (argc != 1){
        cout << "usage: phyx_test" << endl;
        exit(0);
    }
    string datafile = "../../../projects/PHLAWD_fish/12S.keep";
    vector<Sequence> seqs;
    FastaUtil pr;
    bool phyl = pr.readFile(datafile,seqs);
    cout << "sequences: " << seqs.size() << endl;
    cout << seqs[0].get_sequence() <<endl;
    seqs[0].reverse_complement();
    cout << endl;
    cout << seqs[0].get_sequence()<< endl;
    cout << "writing file" << endl;
    string outfile = "test.fasta";
    pr.writeFileFromVector(outfile,seqs);
    */
    
    cout << "werked" << endl;
    
    return EXIT_SUCCESS;
}
