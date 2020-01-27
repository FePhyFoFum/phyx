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
        std::cout << "usage: phyx_test" << std::endl;
        exit(0);
    }
    string datafile = "../../../projects/PHLAWD_fish/12S.keep";
    vector<Sequence> seqs;
    FastaUtil pr;
    bool phyl = pr.readFile(datafile,seqs);
    std::cout << "sequences: " << seqs.size() << std::endl;
    std::cout << seqs[0].get_sequence() <<std::endl;
    seqs[0].reverse_complement();
    std::cout << std::endl;
    std::cout << seqs[0].get_sequence()<< std::endl;
    std::cout << "writing file" << std::endl;
    string outfile = "test.fasta";
    pr.writeFileFromVector(outfile,seqs);
    */
    
    std::cout << "werked" << std::endl;
    
    return EXIT_SUCCESS;
}
