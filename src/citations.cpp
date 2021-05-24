#include <string>
#include "citations.h"

// put citations for various methods here.

std::string get_phyx_citation () {
    std::string cite = "Brown, Joseph W., Joseph F. Walker, and Stephen A. Smith. 2017. ";
    cite += "Phyx: phylogenetic tools for unix. Bioinformatics 33(12): 1886-1888. ";
    cite += "https://doi.org/10.1093/bioinformatics/btx063";
    return cite;
}

std::string get_NW_citation () {
    std::string cite = "Needleman, Saul B., and Christian D. Wunsch. 1970. ";
    cite += "A general method applicable to the search for similarities in the amino acid sequence of two proteins. ";
    cite += "Journal of Molecular Biology. 48 (3): 443–53. ";
    cite += "https://doi.org/10.1016/0022-2836(70)90057-4";
    return cite;
}

std::string get_SW_citation () {
    std::string cite = "Smith, Temple F., and Michael S. Waterman. 1981. ";
    cite += "Identification of Common Molecular Subsequences. ";
    cite += "Journal of Molecular Biology. 147 (1): 195–197. ";
    cite += "https://doi.org/10.1016/0022-2836(81)90087-5";
    return cite;
}
