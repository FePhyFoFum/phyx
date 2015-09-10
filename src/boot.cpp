#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <ctime>
#include <fstream>

using namespace std;

#include "boot.h"
#include "utils.h"

SequenceSampler::SequenceSampler (int const& seed, float const& jackfract, string & partf)
:jkfract(jackfract), jackknife(false), partitioned(false), numPartitionedSites(0), numPartitions(0)
{
    if (seed == -1) {
        srand((unsigned)time(NULL));
    } else {
        srand(seed);
    }
    if (jkfract != 0.0) {
        jackknife = true;
    }
    if (!partf.empty()) {
        partitioned = true;
        parsePartitions(partf);
    }
}

// not used
vector <int> SequenceSampler::get_sampled_sites () {
    return samplesites;
}

string SequenceSampler::get_resampled_seq (string const& origseq) {
    string seq;
    for (int i = 0; i < (int)samplesites.size(); i++) {
        if (i == 0) {
            seq = origseq[samplesites[i]];
        } else {
            seq += origseq[samplesites[i]];
        }
    }
    return seq;
}

// this is done once, to populate site sample vector
void SequenceSampler::sample_sites (int const& numchar) {
    if (partitioned) {
        samplesites = get_partitioned_bootstrap_sites();
    } else if (!jackknife) {
        samplesites = get_bootstrap_sites(numchar);
    } else {
        samplesites = get_jackknife_sites(numchar);
    }
}

// sample with replacement.
vector <int> SequenceSampler::get_bootstrap_sites (int const& numchar) {
    vector <int> randsites (numchar); // numchar zero-initialized elements
    int randnum = 0;
    
    for (int i = 0; i < numchar; i++) {
        randnum = random_int_range(0, (numchar - 1));
        randsites[i] = randnum;
    }
    sort(randsites.begin(), randsites.end());
    
    return randsites;
}

// set up so same composition as original
vector <int> SequenceSampler::get_partitioned_bootstrap_sites () {
    vector <int> master(numPartitionedSites, 0);
    
    for (int i = 0; i < numPartitions; i++) {
        int curNum = (int)partitions[i].size();
        //cout << "Partition #" << i << " contains " << curNum << " sites." << endl;
        vector <int> randsites = get_bootstrap_sites(curNum);
        for (int j = 0; j < curNum; j ++) {
        // put partitions back in same spot as original, so partition file does not need to change
            master[partitions[i][j]] = partitions[i][randsites[j]];
        }
    }
    
    return master;
}

// sample WITHOUT replacement. not with partitioned models
vector <int> SequenceSampler::get_jackknife_sites (int const& numchar) {
    int numsample = numchar * jkfract + 0.5;
    int randnum = 0;
    
    if (numsample == 0) {
        cout << "Jackknife fraction " << jkfract << " leaves no characters remaining." << endl;
        exit(0);
    }
    
    vector <int> randsites (numsample); // numchar zero-initialized elements
    
// ugh. must be a more succinct way to do this.
    vector <int> allsites (numchar);
    for (int i = 0; i < numchar; i++) {
        allsites[i] = i;
    }
    
    for (int i = 0; i < numsample; i++) {
        randnum = random_int_range(i, (numchar - 1));
    // swap, so don't have to worry about multiple hits
        swap(allsites[i], allsites[randnum]);
        randsites[i] = allsites[i];
    }
    
    return randsites;
}


/*
grab partition information from separate file. example:
mtDNA_1st = 1-2066\3
TGFb2 = 5059-5721
*/
void SequenceSampler::parsePartitions (string & partf) {
    vector <int> temp;
    string line;
    ifstream infile(partf.c_str());
    
    while (getline(infile,line)) {
        if (line.size() < 1) {
            continue;
        } else {
            temp = getPartitionSites(line);
            partitions.push_back(temp);
            temp.clear();
        }
    }
    infile.close();
    
    numPartitions = (int)partitions.size();
    calNumPartitionedSites();
    
    // do error-checking here:
    checkValidPartitions();
}

// expecting pattern: name = start-end[\3]
// want to be flexible with spaces/lackthereof
// guaranteed to be ordered
vector <int> SequenceSampler::getPartitionSites (string const& part) {
    vector <int> sites;
    vector <string> tokens;
    int start = 0;
    int stop = 0;
    int interval = 1;
    
    string delim(" -=;\t\\");
    tokenize(part, tokens, delim);
    getPartitionParameters (tokens, start, stop, interval);
    
    int i = start;
    
    while (i <= stop) {
        sites.push_back(i);
        i += interval;
    }
    
    return sites;
}

// opposite of getPartitionSites. want single vector listing site-specific partitions. e.g.
// 1231231231231234444444444444
// not used
void SequenceSampler::getSitePartitions () {
    vector <int> sites(numPartitionedSites, 0);
    
    for (int i = 0; i < numPartitions; i++) {
        for (int j = 0; j < (int)partitions[i].size(); j++) {
            sites[partitions[i][j]] = i;
        }
    }
    sitePartitions = sites;
}

// GADPH = 2991-3406\3
// after being tokenized, should be of length 3 or 4 (latter when interval)
// convert from 1-start to 0-start
void SequenceSampler::getPartitionParameters (vector <string> & tokens, int & start, int & stop, int & interval) {
    if ((int)tokens.size() < 3 || (int)tokens.size() > 4) {
        cout << "Error: invalid/unsupported partition specification." << endl;
        exit (0);
    }
    
    partitionNames.push_back(tokens[0]);
    //cout << "Processing partition '" << tokens[0] << "': ";
    if (is_number(tokens[1])) {
        start = atoi(tokens[1].c_str()) - 1;
        //cout << "start = " << start;
    }
    if (is_number(tokens[2])) {
        stop = atoi(tokens[2].c_str()) - 1;
        //cout << "; stop = " << stop;
    }
    if (((int)tokens.size() == 4) && is_number(tokens[3])) {
        interval = atoi(tokens[3].c_str());
        //cout << "; interval = " << interval;
    }
    //cout << endl;
}

void SequenceSampler::calNumPartitionedSites () {
    numPartitionedSites = 0;
    
    for (int i = 0; i < numPartitions; i++) {
        numPartitionedSites += (int)partitions[i].size();
//         cout << "Partition #" << i << " contains " << (int)partitions[i].size() << " sites:" << endl;
//         for (int j = 0; j < (int)partitions[i].size(); j++) {
//             cout << partitions[i][j] << " ";
//         }
//         cout << endl;
    }
}

int SequenceSampler::getNumPartitionedSites () {
    return numPartitionedSites;
}

// should do some error-checking e.g. for 1) missing sites, 2) overlapping partitions
void SequenceSampler::checkValidPartitions () {
    vector <int> allSites = partitions[0];
    for (int i = 1; i < numPartitions; i++) {
        allSites.insert(allSites.end(), partitions[i].begin(), partitions[i].end());
    }
    sort(allSites.begin(), allSites.end());
    
    int max = allSites.back();
    int count = (int)allSites.size();
    int diff = max - count + 1;
    
    if (diff != 0) { // sites are duplicated
        //cout << "Error in partitioning: maximum site value " << max << " does not equal site count " << count << "." << endl;
        findDuplicatesMissing(allSites);
    }
}

void SequenceSampler::findDuplicatesMissing (vector <int> const& allSites) {
    vector <int> unique;
    vector <int> duplicates;
    vector <int> missing;
    
    int maxVal = allSites.back();
    
    vector <int> counts(maxVal, 0);
    for (int i = 0; i < (int)allSites.size(); i++) {
        counts[allSites[i]]++;
    }
    
    for (int i = 0; i < maxVal; i++) {
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
    
    if (duplicates.size() != 0) {
        cout << "Error: the following " << duplicates.size() << " sites are found in more than one partition: ";
        for (int i = 0; i < (int)duplicates.size(); i++) {
            cout << duplicates[i] << " ";
        }
        cout << endl;
    }
    if (missing.size() != 0) {
        cout << "Error: the following " << missing.size() << " sites are not found in any partition: ";
        for (int i = 0; i < (int)missing.size(); i++) {
            cout << missing[i] << " ";
        }
        cout << endl;
    }
    cout << "Exiting." << endl;
    exit (0);
}
