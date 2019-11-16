#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "vcf_reader.h"
#include "utils.h"

VcfReader::VcfReader(std::istream* pios) {
    read_vcf(pios);
}


void VcfReader::read_vcf (std::istream* pios) {
    bool started = false;
    bool first = true;
    // these column numbers should be constant (i think?), but let's not leave anything to chance
    int refcol = 0;
    int altcol = 0;
    int taxstartcol = 0;
    int ncols = 0;
    
    std::string line;
    while (getline(*pios, line)) {
        std::vector<std::string> temp = tokenize(line);
        if (started) {
            // every line below the bottom header line should be site data
            std::string refstate = temp[refcol];
            std::string terp = temp[altcol];
            std::vector<std::string> states = get_alts(terp);
            // put all observed states in same vector
            states.insert(states.begin(), refstate);
            
            if (!first) {
                int counter = 0;
                for (int i = taxstartcol; i < ncols; i++) {
                    int idx = stoi(temp[i]);
                    seqs_[counter] += states[idx];
                    counter++;
                }
            } else {
                // construct result vector during first data row (site)
                for (int i = taxstartcol; i < ncols; i++) {
                    int idx = stoi(temp[i]);
                    seqs_.push_back(states[idx]);
                }
                first = false;
            }
        } else {
            // skip preceeding lines
            if (temp[0] == "#CHROM") {
                bool read_taxa = false;
                ncols = temp.size();
                for (unsigned int i = 1; i < temp.size(); i++) {
                    if (read_taxa) {
                        taxa_.push_back(temp[i]);    
                        //cout << " " << temp[i];
                    } else {
                        if (temp[i] =="FORMAT") {
                            // here assuming taxa cols directly follow this column
                            read_taxa = true;
                            taxstartcol = i+1;
                            //cout << std::endl << "Taxa:";
                        } else if (temp[i] =="REF") {
                            refcol = i;
                        } else if (temp[i] =="ALT") {
                            altcol = i;
                        }
                    }
                }
                //cout << std::endl;
                started = true;
            }
        }
    }
    //cout << "Read in " << taxa_.size() << " taxa, each with " << seqs_[0].size()
    //    << " characters." << std::endl;
}



// *** seems like a general function would be useful here 

// split alt states on comma
std::vector<std::string> VcfReader::get_alts (const std::string& str) {
    std::vector<std::string> res;
    std::stringstream ss(str);
    while(ss.good()) {
        std::string substr;
        getline(ss, substr, ',');
        res.push_back(substr);
    }
    return res;
}


void VcfReader::write_seqs (bool const& uppercase, std::ostream* poos) {
    for (unsigned int i = 0; i < taxa_.size(); i++) {
        (*poos) << ">" << taxa_[i] << std::endl;
        if (uppercase) {
            std::string terp = string_to_upper(seqs_[i]);
            (*poos) << terp << std::endl;
        } else {
            (*poos) << seqs_[i] << std::endl;
        }
    }
}
