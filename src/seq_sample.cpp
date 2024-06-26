#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <fstream>

#include "sequence.h"
#include "seq_reader.h"
#include "seq_sample.h"
#include "seq_utils.h"
#include "utils.h"


SequenceSampler::SequenceSampler (std::istream* pios, const long int& seed,
        const double& jackfract, std::string& partf):num_taxa_(0u), num_char_(0u),
        jkfract_(jackfract), jackknife_(false), partitioned_(false),
        num_partitioned_sites_(0u), num_partitions_(0u) {
    if (seed == -1) {
        srand(get_clock_seed());
    } else {
        srand(static_cast<unsigned int>(seed));
    }
    if (jkfract_ != 0.0) {
        jackknife_ = true;
    }
    if (!partf.empty()) {
        partitioned_ = true;
        parse_partitions(partf);
    }
    read_in_sequences(pios);
}


void SequenceSampler::read_in_sequences (std::istream* pios) {
    std::string alphaName; // not used, but required by reader
    seqs_ = ingest_alignment(pios, alphaName);
    num_taxa_ = static_cast<unsigned int>(seqs_.size());
    num_char_ = static_cast<unsigned int>(seqs_[0].get_length());
    
    // check that it is aligned (doesn't make sense otherwise)
    if (!is_aligned(seqs_)) {
        std::cerr << "Error: sequences are not aligned. Exiting." << std::endl;
        exit(0);
    }
    
    // check that partitions lengths match alignment
    if (partitioned_) {
        if (num_char_ != get_num_partitioned_sites()) {
            std::cerr << "Error: the number of characters in the sequence (" << num_char_ <<
                ") does not match that in partition file (" << get_num_partitioned_sites() <<
                "). Exiting." << std::endl;
        }
    }
    
    // good to go
    sample_sites(num_char_);
}


// not currently used
std::vector<unsigned int> SequenceSampler::get_sampled_sites () const {
    return sample_sites_;
}


std::string SequenceSampler::get_resampled_seq (const std::string& origseq) {
    std::string seq;
    for (size_t i = 0; i < sample_sites_.size(); i++) {
        if (i == 0) {
            seq = origseq[sample_sites_[i]];
        } else {
            seq += origseq[sample_sites_[i]];
        }
    }
    return seq;
}


// this is done once, to populate site sample vector
void SequenceSampler::sample_sites (const unsigned int& numchar) {
    if (partitioned_) {
        sample_sites_ = get_partitioned_bootstrap_sites();
    } else if (!jackknife_) {
        sample_sites_ = get_bootstrap_sites(numchar);
    } else {
        sample_sites_ = get_jackknife_sites(numchar);
    }
}


// sample with replacement.
std::vector<unsigned int> SequenceSampler::get_bootstrap_sites (const unsigned int& numchar) {
    // numchar zero-initialized elements
    size_t nc = static_cast<size_t>(numchar);
    std::vector<unsigned int> randsites(nc);
    
    for (size_t i = 0; i < nc; i++) {
        unsigned int randnum = random_int_range(0, numchar - 1);
        randsites[i] = static_cast<unsigned int>(randnum);
    }
    sort(randsites.begin(), randsites.end());
    
    return randsites;
}


// set up so same composition as original
std::vector<unsigned int> SequenceSampler::get_partitioned_bootstrap_sites () {
    std::vector<unsigned int> master(num_partitioned_sites_, 0);
    for (size_t i = 0; i < num_partitions_; i++) {
        size_t curNum = partitions_[i].size();
        //std::cout << "Partition #" << i << " contains " << curNum << " sites." << std::endl;
        std::vector<unsigned int> randsites = get_bootstrap_sites(static_cast<unsigned int>(curNum));
        for (size_t j = 0; j < curNum; j ++) {
        // put partitions back in same spot as original, so partition file does not need to change
            master[partitions_[i][j]] = partitions_[i][randsites[j]];
        }
    }
    return master;
}

// sample WITHOUT replacement. not with partitioned models
std::vector<unsigned int> SequenceSampler::get_jackknife_sites (const unsigned int& numchar) {
    auto numsample = static_cast<unsigned int>(std::lround(numchar * jkfract_ + 0.5));
    std::vector<unsigned int> randsites = sample_without_replacement(numchar, numsample);
    return randsites;
}


/*
grab partition information from separate file. example:
mtDNA_1st = 1-2066\3
TGFb2 = 5059-5721
TODO: update to 1) RAxML and 2) Nexus style
 - RAxML: DNA, gene0 = 1-4243
 - MrBayes: CHARSET gene0 = 1-4243;
 * So, thrown out first token, and should work for both
*/
void SequenceSampler::parse_partitions (std::string& partf) {
    std::vector<unsigned int> temp;
    std::string line;
    std::ifstream infile(partf.c_str());
    
    while (getline_safe(infile, line)) {
        if (line.empty()) {
            continue;
        }
        temp = get_partition_sites(line);
        partitions_.push_back(temp);
        temp.clear();
    }
    infile.close();
    
    num_partitions_ = static_cast<unsigned int>(partitions_.size());
    //std::cout << "Found " << num_partitions_ << " partitions." << std::endl;
    calculate_num_partitioned_sites();
    
    // do error-checking here:
    check_valid_partitions();
}


// expecting pattern: (CHARSET/DNA,) name = start-end[\3][,]
// want to be flexible with spaces/lackthereof
// guaranteed to be ordered
std::vector<unsigned int> SequenceSampler::get_partition_sites (const std::string& part) {
    std::vector<unsigned int> sites;
    std::vector<std::string> tokens;
    unsigned int start = 0;
    unsigned int stop = 0;
    unsigned int interval = 1;
    
    std::string delim(" -=;\t\\");
    tokenize(part, tokens, delim);
    get_partition_parameters(tokens, start, stop, interval);
    
    unsigned int i = start;
    while (i <= stop) {
        sites.push_back(i);
        i += interval;
    }
    return sites;
}


// opposite of get_partition_sites. want single vector listing site-specific partitions. e.g.
// 1231231231231234444444444444
// not currently used
/*
void SequenceSampler::get_site_partitions () {
    std::vector<unsigned int> sites(static_cast<size_t>(num_partitioned_sites_), 0);
    for (unsigned long i = 0; i < static_cast<size_t>(num_partitions_); i++) {
        for (unsigned long j = 0; j < partitions_[i].size(); j++) {
            sites[partitions_[i][j]] = static_cast<unsigned int>(i);
        }
    }
    site_partitions_ = sites;
}
*/


// CHARSET GADPH = 2991-3406\3
// after being tokenized, should be of length 4 or 5 (latter when interval)
// convert from 1-start to 0-start
void SequenceSampler::get_partition_parameters (std::vector<std::string>& tokens,
        unsigned int& start, unsigned int& stop, unsigned int& interval) {
    auto tsize = static_cast<unsigned int>(tokens.size());
    if (tsize < 4u || tsize > 5u) {
        std::cerr << "Error: invalid/unsupported partition specification. Exiting." << std::endl;
        exit(0);
    }
    
    // ignore first token. will be either CHARSET or DNA,
    partition_names_.push_back(tokens[1]);
    //std::cout << "Processing partition '" << tokens[1] << "': ";
    
    // why is there not std::stou?!?
    
    if (is_number(tokens[2])) {
        start = static_cast<unsigned int>(std::stoul(tokens[2])) - 1;
        //std::cout << "start = " << start;
    }
    if (is_number(tokens[3])) {
        stop =  static_cast<unsigned int>(std::stoul(tokens[3])) - 1;
        //std::cout << "; stop = " << stop;
    }
    if ((static_cast<int>(tokens.size()) == 5) && is_number(tokens[4])) {
        interval =  static_cast<unsigned int>(std::stoul(tokens[4]));
        //std::cout << "; interval = " << interval;
    }
    //std::cout << std::endl;
}


void SequenceSampler::calculate_num_partitioned_sites () {
    num_partitioned_sites_ = 0;
    for (unsigned int i = 0; i < num_partitions_; i++) {
        num_partitioned_sites_ += static_cast<unsigned int>(partitions_[static_cast<size_t>(i)].size());
//         std::cout << "Partition #" << i << " contains " << (int)partitions[i].size()
//        << " sites:" << std::endl;
//         for (unsigned int j = 0; j < partitions[i].size(); j++) {
//             std::cout << partitions[i][j] << " ";
//         }
//         std::cout << std::endl;
    }
}


unsigned int SequenceSampler::get_num_partitioned_sites () const {
    return num_partitioned_sites_;
}


// should do some error-checking e.g. for 1) missing sites, 2) overlapping partitions
void SequenceSampler::check_valid_partitions () {
    std::vector<unsigned int> allSites = partitions_[0];
    for (size_t i = 1; i < static_cast<size_t>(num_partitions_); i++) {
        allSites.insert(allSites.end(), partitions_[i].begin(), partitions_[i].end());
    }
    sort(allSites.begin(), allSites.end());
    
    unsigned int max = allSites.back();
    auto count = static_cast<unsigned int>(allSites.size());
    unsigned int diff = max - count + 1;
    
    if (diff != 0) { // sites are duplicated
        //std::cerr << "Error in partitioning: maximum site value " << max
        //  << " does not equal site count " << count << "." << std::endl;
        find_duplicates_missing(allSites);
    }
}


// find sites 1) present in more than 1 partition or 2) not present in any
[[ noreturn ]] void SequenceSampler::find_duplicates_missing (const std::vector<unsigned int>& allSites) {
    std::vector<unsigned int> unique;
    std::vector<unsigned int> duplicates;
    std::vector<unsigned int> missing;
    
    unsigned int maxVal = allSites.back();
    
    std::vector<unsigned int> counts(maxVal, 0);
    for (unsigned int allSite : allSites) {
        counts[static_cast<size_t>(allSite)]++;
    }
    
    for (unsigned int i = 0; i < maxVal; i++) {
        switch (counts[i]) {
            case 0:
                missing.push_back(i);
                break;
            case 1:
                unique.push_back(i);
                break;
            default:
                duplicates.push_back(i);
        }
    }
    
    if (!duplicates.empty()) {
        std::cerr << "Error: the following " << duplicates.size()
                << " sites are found in more than one partition: ";
        for (unsigned int duplicate : duplicates) {
            std::cerr << duplicate << " ";
        }
        std::cerr << std::endl;
    }
    if (!missing.empty()) {
        std::cerr << "Error: the following " << missing.size()
                << " sites are not found in any partition: ";
        for (unsigned int ms : missing) {
            std::cerr << ms << " ";
        }
        std::cerr << std::endl;
    }
    std::cerr << "Exiting." << std::endl;
    exit(0);
}


void SequenceSampler::write_resampled_seqs (std::ostream* poos) {
    for (size_t i = 0; i < static_cast<size_t>(num_taxa_); i++) {
        (*poos) << ">" << seqs_[i].get_id() << std::endl;
        (*poos) << get_resampled_seq(seqs_[i].get_sequence()) << std::endl;
    }
}
