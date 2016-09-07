/*
 * CodonAligner.cpp
 *
 *  Created on: Jan 12, 2015
 *      Author: joe
 */

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>

using namespace std;

//Makes a vector with fasta sequences on one line
vector <string> ReadFasta(string& fasta) {
    string line, dna;
    ifstream readline;
    int count = 0;
    vector<string> parsed_fasta;
    readline.open(fasta.c_str());
    if (readline.is_open()) {
        while (getline (readline, line)) {
            if (line[0] == '>') {
                if (count != 0) {
                    parsed_fasta.push_back(dna);
                    parsed_fasta.push_back(line);
                    dna = "";
                } else {
                    parsed_fasta.push_back(line);
                }
            } else {
                count++;
                dna += line;
            }
        }
    }
    parsed_fasta.push_back(dna);
    return parsed_fasta;
}

//Change amino acid to codon
vector <string> AAtoCodon(vector<string>& parsed_nucleotide, vector<string>& parsed_aa) {
    vector <string> codon_alignment;
    vector <string> codons;
    vector <char> AA;
    string line;
    int count = 0;
    ifstream readline;
    for (vector<string>::size_type i = 0; i != parsed_nucleotide.size(); i = i + 2) {
        for (vector<string>::size_type j = 0; j != parsed_aa.size(); j =  j + 2) {

            if (parsed_aa[j] == parsed_nucleotide[i]) {
                //cout << parsed_aa[j+1] << "\n" << parsed_nucleotide[i+1] << endl;
                codons.clear();
                for (std::string::size_type k = 0; k < parsed_nucleotide[i+1].length(); k++) {
                    if (count == 2) {
                        line = line + parsed_nucleotide[j+1][k];
                        //cout << line << endl;
                        codons.push_back(line);
                        count = 0;
                        line = "";
                    } else {
                        line = line + parsed_nucleotide[j+1][k];
                        count++;
                    }
                }
                int position = 0;
                string seq = "";
                for (std::string::size_type l = 0; l < parsed_aa[j+1].length(); l++) {
                //for (std::string::size_type l = 0; l < codons.size(); l++) {
                   // cout << parsed_aa[j+1][l] << endl;
                    if (parsed_aa[j+1][l] == '-') {
                  //      cout << "---";
                        seq = seq + "---";
                    } else {
                    //    cout << codons[position];
                        seq = seq + codons[position];
                        position++;
                    }
                }
                codon_alignment.push_back(parsed_nucleotide[i]);
                codon_alignment.push_back(seq);
                //cout << endl;
            }
        }
    }
    return codon_alignment;
}

int main() {
    string unaligned_nucleotide, aligned_aa;
    vector<string> parsed_nucleotide, codon_alignment, parsed_aa;

    cout << "Please input unaligned nucleotide sequences" << endl;
    cin >> unaligned_nucleotide;
    cout << "Please input aligned Amino Acid sequences" << endl;
    cin >> aligned_aa;

    parsed_nucleotide = ReadFasta(unaligned_nucleotide);
    parsed_aa = ReadFasta(aligned_aa);
    codon_alignment = AAtoCodon(parsed_nucleotide, parsed_aa);
    for (vector<string>::size_type i = 0; i != codon_alignment.size(); i++) {
        cout << codon_alignment[i] << endl;
    }

    return EXIT_SUCCESS;
}
