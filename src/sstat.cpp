#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

#include "utils.h"
#include "sequence.h"
#include "sstat.h"


MultinomialSeqStat::MultinomialSeqStat (vector<Sequence> & seqs) {
    seqs_ = seqs;
    num_taxa_ = seqs.size();
    
    //cout << "Read in " << num_taxa_ << " sequences!" << endl;
    
    if (!checked_aligned()) {
        cout << "Cannot calculate statistic as sequences are not aligned. Exiting." << endl;
        exit(0);
    }
    
    collect_site_patters();
    calculateTestStatistic();
}

bool MultinomialSeqStat::checked_aligned () {
    bool is_aligned_ = true;
    vector <int> seq_lengths (num_taxa_, 0);
    
    // gather all lengths
    for (int i = 0; i < num_taxa_; i++) {
        seq_lengths[i] = (int)seqs_[i].get_sequence().size();
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
    vector <string> inputPatterns;
    for (int i = 0; i < num_char_; i++) {
        char a = seqs_[0].get_sequence().at(i);
        string pat(1, a);
        for (int j = 1; j < num_taxa_; j++) {
            char residue = seqs_[j].get_sequence().at(i);
            pat += residue;
        }
        inputPatterns.push_back(pat);
    }
    // *** Sort - enables more rapid pattern count enumeration later on ***
    sort(inputPatterns.begin(), inputPatterns.end());
    
    // process raw counts
    int numUniquePatterns = 0;
    bool siteMatched = true;
    
// Iterate from first site to last site
    for (vector<string>::const_iterator iterSites = inputPatterns.begin(); iterSites < inputPatterns.end(); ++iterSites) {
        string currentPattern = *iterSites;
        int currentPatternCount = 1;
        
// Check if pattern n+1 == pattern n; if so, update count, test next site; if not, log current count/pattern; MUCH FASTER!!
//     - Need to check if last site, lest risk overshooting vector boundary
        siteMatched = true;
        while (siteMatched) {
            if (iterSites != inputPatterns.end() - 1) { // Last site
                iterSites++;
                if (*iterSites == currentPattern) {
                    currentPatternCount++;
                } else {    // Exit loop, counted all matching patterns
                    iterSites--;
                    siteMatched = false;
                }
            } else {    // At last pattern
                siteMatched = false;
            }
        }
// No need to check against stored sites, as vector has been sorted in function 'getInputPatterns'
        patterns_and_counts_.push_back(make_pair(currentPattern, currentPatternCount));
        numUniquePatterns++;
    }
    //cout << "Found " << numUniquePatterns << " unique site patterns." << endl;
}

// Calculate T(X) test statistic
void MultinomialSeqStat::calculateTestStatistic () {
/*----------------------------------------------------------------------------------------------------------
| Calculate summary statistics from Bollback (2002)
| 1) Calculate n*ln(n) for each pattern
| 2) Sum values above
| 3) Subtract N*ln(N) from sum
*/
    long double tempResult = 0.0;
    long double finalResult = 0.0;
    int numPatterns = patterns_and_counts_.size();

    for (int iterStoredPat = 0; iterStoredPat < numPatterns; ++iterStoredPat) {
        if (patterns_and_counts_[iterStoredPat].second != 0) {
            tempResult = (patterns_and_counts_[iterStoredPat].second)*(log(patterns_and_counts_[iterStoredPat].second));
            finalResult += tempResult;
            tempResult = 0.0;
        }
    }
    tempResult = (double)num_char_*(log(num_char_));
    /*
        cout << "'numIncludedChar*log(numIncludedChar)' calculated as: " << tempResult << endl;
        cout << "'final sum' is: " << finalResult << endl;
        cout << "'numIncludedChar' is: " << num_char_ << endl;
        cout << "'log(numIncludedChar)' is: " << log(num_char_) << endl;
    */
    finalResult = finalResult - tempResult;
    //cout << numPatterns << " distinct site patterns. T(X) = " << finalResult << endl;
    
    test_statistic_ = finalResult;
}

double MultinomialSeqStat::get_test_statistic () {
    return test_statistic_;
}
