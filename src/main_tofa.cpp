/*
 * main_mrca.cpp
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#include "utils.h"
#include "seq_reader.h"
#include "sequence.h"

int main(int argc, char * argv[]){

    if (argc > 2){
	cout << "usage: pxtofa file" << endl;
	exit(0);
    }
    
    if (argc == 1){
	string retstring;
	int ft = test_seq_filetype_stream(cin,retstring);
	Sequence seq;
	while(read_next_seq_from_stream(std::cin,ft,retstring,seq)){
	    cout << seq.get_fasta();
	}
    }else if(argc == 2){
	string retstring;
	Sequence seq;
	ifstream fstr(argv[1]);
	int ft = test_seq_filetype_stream(fstr,retstring);
	cout << ft << endl;
	while(read_next_seq_from_stream(fstr,ft,retstring,seq)){
	    cout << seq.get_fasta();
	}
	//fasta has a trailing one
	if (ft == 2){
	    cout << seq.get_fasta();
	}
	fstr.close();
    }
    return EXIT_SUCCESS;
}
