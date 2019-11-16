#include <iostream>
#include <fstream>
#include <string>

#include "utils.h"
#include "seq_reader.h"

int main(int argc, char * argv[]) {

    if (argc > 2 || argc == 1) {
        std::cout << "usage: pxseqtest file" << std::endl;
        exit(0);
    }

    std::cout << test_seq_filetype(argv[1]) << std::endl;
    return EXIT_SUCCESS;
}
