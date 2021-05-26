#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <cassert>
#include <fstream>

#include "seq_models.h"
#include "sequence.h"
#include "utils.h"


void read_scoring_matrix (char * filename, std::map<char, std::map<char, int> >& sc_mat) {
    std::ifstream fstr(filename);
    sc_mat.clear();
    std::vector<char> order;
    bool first = true;
    std::string line;
    while (getline_safe(fstr, line)) {
        if (line[0] == '#' || line.empty()) {
            continue;
        }
        std::vector<std::string> tokens;
        std::string del(" \t");
        tokenize(line, tokens, del);
        for (auto & tk : tokens) {
            trim_spaces(tk); // this will never be used, as it was split on whitespace
        }
        if (first) {
            first = false;
            for (auto & tk : tokens) {
                order.push_back(tk[0]);
            }
            for (char j : order) {
                sc_mat[j] = std::map<char, int>();
            }
            continue;
        }
        for (unsigned int j = 0; j < order.size(); j++) {
            sc_mat[tokens[0][0]][order[j]] = std::atoi(tokens[j+1].c_str()); //#changed from int to float
        }
    }
    fstr.close();
}


void read_scoring_matrix_from_lines(std::vector<std::string>& lines, std::map<char, std::map<char, int> >& sc_mat) {
    sc_mat.clear();
    std::vector<char> order;
    bool first = true;
    for (auto line : lines) {
        if (line[0] == '#') {
            continue;
        }
        std::vector<std::string> tokens;
        std::string del(" \t");
        tokenize(line, tokens, del);
        for (auto & tk : tokens) {
            trim_spaces(tk); // this will never be used, as it was split on whitespace
        }
        if (first) {
            first = false;
            for (auto & tk : tokens) {
                order.push_back(tk[0]);
            }
            for (char j : order) {
                sc_mat[j] = std::map<char, int>();
            }
            continue;
        }
        for (unsigned int j = 0; j < order.size(); j++) {
            sc_mat[tokens[0][0]][order[j]] = std::atoi(tokens[j+1].c_str()); //#changed from int to float
        }
    }
}


void get_ednafull(std::map<char, std::map<char, int> >& sc_mat) {
    std::vector<std::string> sts;
    sts.emplace_back("#");
    sts.emplace_back("# This matrix was created by Todd Lowe   12/10/92");
    sts.emplace_back("#");
    sts.emplace_back("# Uses ambiguous nucleotide codes, probabilities rounded to");
    sts.emplace_back("#  nearest integer");
    sts.emplace_back("#");
    sts.emplace_back("# Lowest score = -4, Highest score = 5");
    sts.emplace_back("#");
    sts.emplace_back("    A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N   U");
    sts.emplace_back("A   5  -4  -4  -4  -4   1   1  -4  -4   1  -4  -1  -1  -1  -2  -4");
    sts.emplace_back("T  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2   5");
    sts.emplace_back("G  -4  -4   5  -4   1  -4   1  -4   1  -4  -1  -1  -4  -1  -2  -4");
    sts.emplace_back("C  -4  -4  -4   5   1  -4  -4   1  -4   1  -1  -1  -1  -4  -2  -4");
    sts.emplace_back("S  -4  -4   1   1  -1  -4  -2  -2  -2  -2  -1  -1  -3  -3  -1  -4");
    sts.emplace_back("W   1   1  -4  -4  -4  -1  -2  -2  -2  -2  -3  -3  -1  -1  -1   1");
    sts.emplace_back("R   1  -4   1  -4  -2  -2  -1  -4  -2  -2  -3  -1  -3  -1  -1  -4");
    sts.emplace_back("Y  -4   1  -4   1  -2  -2  -4  -1  -2  -2  -1  -3  -1  -3  -1   1");
    sts.emplace_back("K  -4   1   1  -4  -2  -2  -2  -2  -1  -4  -1  -3  -3  -1  -1   1");
    sts.emplace_back("M   1  -4  -4   1  -2  -2  -2  -2  -4  -1  -3  -1  -1  -3  -1  -4");
    sts.emplace_back("B  -4  -1  -1  -1  -1  -3  -3  -1  -1  -3  -1  -2  -2  -2  -1  -1");
    sts.emplace_back("V  -1  -4  -1  -1  -1  -3  -1  -3  -3  -1  -2  -1  -2  -2  -1  -4");
    sts.emplace_back("H  -1  -1  -4  -1  -3  -1  -3  -1  -3  -1  -2  -2  -1  -2  -1  -1");
    sts.emplace_back("D  -1  -1  -1  -4  -3  -1  -1  -3  -1  -3  -2  -2  -2  -1  -1  -1");
    sts.emplace_back("N  -2  -2  -2  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2");
    sts.emplace_back("U  -4   5  -4  -4  -4   1  -4   1   1  -4  -1  -4  -1  -1  -2   5");
    read_scoring_matrix_from_lines(sts, sc_mat);
}


void get_blosum62(std::map<char, std::map<char, int> >& sc_mat) {
    std::vector<std::string> sts;
    sts.push_back("#  Matrix made by matblas from blosum62.iij");
    sts.push_back("#  * column uses minimum score");
    sts.push_back("#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units");
    sts.push_back("#  Blocks Database = /data/blocks_5.0/blocks.dat");
    sts.push_back("#  Cluster Percentage: >= 62");
    sts.push_back("#  Entropy =   0.6979, Expected =  -0.5209");
    sts.push_back("   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *");
    sts.push_back("A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 ");
    sts.push_back("R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 ");
    sts.push_back("N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 ");
    sts.push_back("D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 ");
    sts.push_back("C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 ");
    sts.push_back("Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 ");
    sts.push_back("E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 ");
    sts.push_back("G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 ");
    sts.push_back("H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 ");
    sts.push_back("I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 ");
    sts.push_back("L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 ");
    sts.push_back("K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 ");
    sts.push_back("M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 ");
    sts.push_back("F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 ");
    sts.push_back("P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 ");
    sts.push_back("S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 ");
    sts.push_back("T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 ");
    sts.push_back("W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 ");
    sts.push_back("Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 ");
    sts.push_back("V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 ");
    sts.push_back("B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 ");
    sts.push_back("Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 ");
    sts.push_back("X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 ");
    sts.push_back("* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 ");
    read_scoring_matrix_from_lines(sts, sc_mat);
}
