/*
 * main_seq_test.cpp
 *
 */

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#include "utils.h"
#include "seq_reader.h"

int main(int argc, char * argv[]) {

    if (argc > 2) {
        cout << "usage: pxseqtest file" << endl;
        exit(0);
    }

    cout << test_seq_filetype(argv[1]) << endl;
    return EXIT_SUCCESS;
}
