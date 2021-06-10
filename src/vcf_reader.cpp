#include <string>
#include <vector>
#include <iostream>
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
    size_t refcol = 0;
    size_t altcol = 0;
    size_t taxstartcol = 0;
    size_t ncols = 0;
    
    std::string line;
    while (getline_safe(*pios, line)) {
        if (line.empty()) {
            continue;
        }
        std::vector<std::string> temp = tokenize(line);
        if (started) {
            // every line below the bottom header line should be site data
            std::string refstate = temp[refcol];
            std::string terp = temp[altcol];
            std::vector<std::string> states = get_alts(terp);
            // put all observed states in same vector
            states.insert(states.begin(), refstate);
            
            if (!first) {
                size_t counter = 0;
                for (size_t i = taxstartcol; i < ncols; i++) {
                    size_t idx = stoul(temp[i]);
                    seqs_[counter] += states[idx];
                    counter++;
                }
            } else {
                // construct result vector during first data row (site)
                for (size_t i = taxstartcol; i < ncols; i++) {
                    size_t idx = stoul(temp[i]);
                    seqs_.push_back(states[idx]);
                }
                first = false;
            }
        } else {
            // skip preceeding lines
            if (temp[0] == "#CHROM") {
                bool read_taxa = false;
                ncols = temp.size();
                for (size_t i = 1; i < ncols; i++) {
                    if (read_taxa) {
                        taxa_.push_back(temp[i]);    
                        //std::cout << " " << temp[i];
                    } else {
                        if (temp[i] =="FORMAT") {
                            // here assuming taxa cols directly follow this column
                            read_taxa = true;
                            taxstartcol = i+1;
                            //std::cout << std::endl << "Taxa:";
                        } else if (temp[i] =="REF") {
                            refcol = i;
                        } else if (temp[i] =="ALT") {
                            altcol = i;
                        }
                    }
                }
                //std::cout << std::endl;
                started = true;
            }
        }
    }
    //std::cout << "Read in " << taxa_.size() << " taxa, each with " << seqs_[0].size()
    //    << " characters." << std::endl;
}



// *** seems like a general function would be useful here 

// split alt states on comma
std::vector<std::string> VcfReader::get_alts (const std::string& str) {
    std::vector<std::string> res;
    std::stringstream ss(str);
    while (ss.good()) {
        std::string substr;
        getline(ss, substr, ',');
        res.push_back(substr);
    }
    return res;
}


void VcfReader::write_seqs (const bool& uppercase, std::ostream* poos) {
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
