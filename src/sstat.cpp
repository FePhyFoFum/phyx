#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#include "utils.h"
#include "sequence.h"
#include "sstat.h"


MultinomialSeqStat::MultinomialSeqStat (std::vector<Sequence>& seqs):num_char_(0),
        test_statistic_(0.0), seqs_(seqs) {
    num_taxa_ = static_cast<int>(seqs_.size());
    
    //std::cout << "Read in " << num_taxa_ << " sequences!" << std::endl;
    
    if (!checked_aligned()) {
        std::cerr << "Error: cannot calculate statistic as sequences are not aligned. Exiting." << std::endl;
        exit(0);
    }
    
    collect_site_patters();
    calculateTestStatistic();
}


bool MultinomialSeqStat::checked_aligned () {
    bool is_aligned_ = true;
    std::vector<int> seq_lengths(static_cast<unsigned long>(num_taxa_), 0);
    
    // gather all lengths
    for (unsigned long i = 0; i < static_cast<unsigned long>(num_taxa_); i++) {
        seq_lengths[i] = static_cast<int>(seqs_[i].get_length());
    }
    
    // check if all seqs are the same length
    if (std::adjacent_find( seq_lengths.begin(), seq_lengths.end(), std::not_equal_to<int>()) == seq_lengths.end() ) {
        is_aligned_ = true;
        num_char_ = seq_lengths[0];
    } else {
        is_aligned_ = false;
    }
    return is_aligned_;
}


void MultinomialSeqStat::collect_site_patters () {
    std::vector<std::string> inputPatterns;
    for (unsigned long i = 0; i < static_cast<unsigned long>(num_char_); i++) {
        char a = seqs_[0].get_sequence().at(i);
        std::string pat(1, a);
        for (unsigned long j = 1; j < static_cast<unsigned long>(num_taxa_); j++) {
            char residue = seqs_[j].get_sequence().at(i);
            pat += residue;
        }
        inputPatterns.push_back(pat);
    }
    // *** Sort - enables more rapid pattern count enumeration later on ***
    sort(inputPatterns.begin(), inputPatterns.end());
    
    // process raw counts
    //int numUniquePatterns = 0;
    
// Iterate from first site to last site
    for (std::vector<std::string>::const_iterator iterSites = inputPatterns.begin(); iterSites < inputPatterns.end(); ++iterSites) {
        std::string currentPattern = *iterSites;
        int currentPatternCount = 1;
        
// Check if pattern n+1 == pattern n; if so, update count, test next site; if not, log current count/pattern; MUCH FASTER!!
//     - Need to check if last site, lest risk overshooting vector boundary
        bool siteMatched = true;
        while (siteMatched) {
            if (iterSites != inputPatterns.end() - 1) { // Last site
                ++iterSites;
                if (*iterSites == currentPattern) {
                    currentPatternCount++;
                } else { // Exit loop, counted all matching patterns
                    --iterSites;
                    siteMatched = false;
                }
            } else { // At last pattern
                siteMatched = false;
            }
        }
// No need to check against stored sites, as vector has been sorted in function 'getInputPatterns'
        //patterns_and_counts_.push_back(std::make_pair(currentPattern, currentPatternCount));
        patterns_and_counts_.emplace_back(currentPattern, currentPatternCount);
        //numUniquePatterns++;
    }
    //std::cout << "Found " << numUniquePatterns << " unique site patterns." << std::endl;
}


// Calculate T(X) test statistic
void MultinomialSeqStat::calculateTestStatistic () {
/*----------------------------------------------------------------------------------------------------------
| Calculate summary statistics from Bollback (2002)
| 1) Calculate n*ln(n) for each pattern
| 2) Sum values above
| 3) Subtract N*ln(N) from sum
*/
    long double tempResult = 0.0L;
    long double finalResult = 0.0L;
    auto numPatterns = static_cast<unsigned long>(patterns_and_counts_.size());

    for (unsigned long iterStoredPat = 0; iterStoredPat < numPatterns; ++iterStoredPat) {
        if (patterns_and_counts_[iterStoredPat].second != 0) {
            tempResult = (patterns_and_counts_[iterStoredPat].second)*(log(patterns_and_counts_[iterStoredPat].second));
            finalResult += tempResult;
        }
    }
    tempResult = static_cast<long double>(num_char_*(log(num_char_)));
    /*
        std::cout << "'numIncludedChar*log(numIncludedChar)' calculated as: " << tempResult << std::endl;
        std::cout << "'final sum' is: " << finalResult << std::endl;
        std::cout << "'numIncludedChar' is: " << num_char_ << std::endl;
        std::cout << "'log(numIncludedChar)' is: " << log(num_char_) << std::endl;
    */
    finalResult = finalResult - tempResult;
    //std::cout << numPatterns << " distinct site patterns. T(X) = " << finalResult << std::endl;
    
    test_statistic_ = finalResult;
}


double MultinomialSeqStat::get_test_statistic () const {
    return test_statistic_;
}
