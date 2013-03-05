#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <ctime>

using namespace std;

#include "seq_samp.h"
#include "utils.h"

SequenceSampler::SequenceSampler(int const& seed, float const& jackfract):jkfract(jackfract), jackknife(false) {
    if (seed == -1) {
        srand(time(NULL));
    } else {
        srand(seed);
    }
    if (jkfract != 0.0) {
    	jackknife = true;
    }
}

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

void SequenceSampler::sample_sites (int const& numchar) {
	if (!jackknife) {
		get_bootstrap_sites(numchar);
	} else {
		get_jackknife_sites(numchar);
	}
}

// sample with replacement
void SequenceSampler::get_bootstrap_sites (int const& numchar) {
	vector <int> randsites (numchar); // numchar zero-initialized elements
	int randnum = 0;
	
	for (int i = 0; i < numchar; i++) {
		randnum = random_int_range(0, (numchar - 1));
		randsites[i] = randnum;
	}
	sort(randsites.begin(), randsites.end());
	
	samplesites = randsites;
}

// sample WITHOUT replacement
void SequenceSampler::get_jackknife_sites (int const& numchar) {
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
	
	samplesites = randsites;
}
