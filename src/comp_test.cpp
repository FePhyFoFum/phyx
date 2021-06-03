#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <cmath>

#include "comp_test.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "seq_utils.h"


CompTest::CompTest (std::istream* pios, std::ostream* poos):num_taxa_(0u), seq_length_(0u),
        total_(0u), df_(0u), test_stat_(0.0), prob_(0.0), seq_chars_(""), alpha_name_(""),
        alpha_set_(false), datatype_set_(false), is_multi_(false), gap_('-'), missing_('?'),
        pios_(pios), poos_(poos), longest_tax_label_(0) {
    read_in_alignment();
    count_chars();
    calc_chi_square();
}


void CompTest::read_in_alignment () {
    seqs_ = ingest_alignment(pios_, alpha_name_);
    
    if (!is_aligned(seqs_)) {
        std::cerr << "Error: sequences must be aligned. Exiting." << std::endl;
        exit(0);
    }
    
    num_taxa_ = static_cast<unsigned int>(seqs_.size());
    seq_length_ = seqs_[0].get_length();
    set_datatype();
    
    // if datatype is multi, but alphabet not set, get from entire concatenated sequence
    if (is_multi_ && !alpha_set_) {
        // grab all unique characters from the input string
        // here, seqs from all individuals are concatenated, so represents all sampled characters
        std::string concatenated;
        for (unsigned int i = 0; i < num_taxa_; i++) {
            concatenated += seqs_[static_cast<unsigned long>(i)].get_sequence();
        }
        seq_chars_ = get_alphabet_from_sequence(concatenated);
        // remove gap and missing (if present)
        seq_chars_.erase(std::remove(seq_chars_.begin(), seq_chars_.end(), gap_), seq_chars_.end());
        seq_chars_.erase(std::remove(seq_chars_.begin(), seq_chars_.end(), missing_), seq_chars_.end());
        alpha_set_ = true;
    }
    
    // resize
    col_totals_.resize(seq_chars_.size(), 0);
    taxon_labels_ = collect_names(seqs_);
}


void CompTest::set_datatype () {
    if (alpha_name_ == "DNA") {
        seq_chars_ = "ACGT";
        alpha_set_ = true;
    } else if (alpha_name_ == "AA") {
        seq_chars_ = "ACDEFGHIKLMNPQRSTVWYX";
        alpha_set_ = true;
    } else if (alpha_name_ == "BINARY") {
        seq_chars_ = "01";
        alpha_set_ = true;
    } else if (alpha_name_ == "MULTI") {
        is_multi_ = true;
        // alphabet is either 1) supplied (nexus) or 2) comes from entire alignment
    } else {
        std::cerr << "Error: cannot determine alignment type. Exiting." << std::endl;
        exit(0);
    }
    datatype_set_ = true;
}


// get counts of all valid character states per taxon
void CompTest::count_chars () {
    std::string seq;
    for (unsigned int i = 0; i < num_taxa_; i++) {
        unsigned int sum = 0;
        seq = string_to_upper(seqs_[static_cast<unsigned long>(i)].get_sequence());
        std::vector<int> icounts(seq_chars_.length(), 0);
        
        for (unsigned int j = 0; j < seq_chars_.length(); j++) {
            auto num = static_cast<unsigned int>(std::count(seq.begin(), seq.end(), seq_chars_[j]));
            icounts[j] += num;
            sum += num;
            col_totals_[j] += num;
        }
        indiv_char_counts_.push_back(icounts);
        // row totals need to be all the same. if not, could be:
        // 1. if not aligned - this is checked
        // 2. missing chars (-, ?, ambig) which are not yet supported here
        row_totals_.push_back(sum);
        total_ += sum;
    }
    // checking row totals
    if (!all_equal(row_totals_)) {
        std::cerr << "Error: missing/ambiguous characters not currently supported. Exiting."
            << std::endl;
        exit(0);
    }
}


// get the longest label. for printing purposes
void CompTest::get_longest_taxon_label () {
    longest_tax_label_ = 0;
    for (unsigned int i = 0; i < num_taxa_; i++) {
        unsigned int cur_len = static_cast<unsigned int>(taxon_labels_[static_cast<unsigned int>(i)].size());
        if (cur_len > longest_tax_label_) {
            longest_tax_label_ = cur_len;
        }
    }
}


void CompTest::return_freq_table () {
    const int colWidth = 12;
    // need to take into account longest_tax_label_
    get_longest_taxon_label();
    std::string pad = std::string(longest_tax_label_, ' ');
    // header
    (*poos_) << "Observed character counts:" << std::endl;
    (*poos_) << pad << " ";
    for (char seq_char : seq_chars_) {
        (*poos_) << std::right << std::setw(colWidth) << seq_char << " ";
    }
    (*poos_) << std::right << std::setw(colWidth) << "Nchar" << std::endl;
    for (unsigned long i = 0; i < static_cast<unsigned long>(num_taxa_); i++) {
        unsigned int diff = longest_tax_label_ - static_cast<unsigned int>(taxon_labels_[i].size());
        (*poos_) << taxon_labels_[i];
        if (diff > 0) {
            pad = std::string(diff, ' ');
            (*poos_) << pad;
        }
        (*poos_) << " ";
        for (unsigned int j = 0; j < seq_chars_.length(); j++) {
            (*poos_) << std::right << std::setw(colWidth) << indiv_char_counts_[i][j] << " ";
        }
        (*poos_) << std::right << std::setw(colWidth) << row_totals_[i] << std::endl;
    }
    unsigned int diff = longest_tax_label_ - 5;
    pad = std::string(diff, ' ');
    (*poos_) << "Total" << pad << " ";
    for (unsigned int col_total : col_totals_) {
        (*poos_) << std::right << std::setw(colWidth) << col_total << " ";
    }
    (*poos_) << std::right << std::setw(colWidth) << total_ << std::endl;
} 


void CompTest::calc_chi_square () {
    test_stat_ = 0.0;
    df_ = static_cast<unsigned int>((num_taxa_ - 1) * (col_totals_.size() - 1));
    double observed = 0.0;
    double expected = 0.0;
    double cellv = 0.0;
    for (unsigned long i = 0; i < static_cast<unsigned long>(num_taxa_); i++) {
        for (unsigned int j = 0; j < col_totals_.size(); j++) {
            observed = static_cast<double>(indiv_char_counts_[i][j]);
            expected = static_cast<double>(col_totals_[j]) * static_cast<double>(row_totals_[i])
                    / static_cast<double>(total_);
            cellv = get_cell_value(observed, expected);
            test_stat_ += cellv;
        }
    }
    prob_ = calc_chi_square_prob(df_, test_stat_);
}

void CompTest::print_results () {
    return_freq_table();
    (*poos_) << "chi-square test stat. = " << test_stat_ << std::endl;
    (*poos_) << "df = " << df_ << std::endl;
    (*poos_) << "prob = " << prob_ << std::endl;
}

// A chi square distribution with n degrees of freedom is equal to
// a gamma distribution with a = n / 2 and b = 0.5 (or β = 2).
// s = df/2, t = chi2/2
double CompTest::calc_chi_square_prob (const double& df, const double& xstat) {
    // prob given by igf(df/2, x/2) / gamma(df/2)
    double s = df / 2;
    double t = xstat / 2;
    double prob = lower_incomplete_gamma_function(s, t);
    prob /= tgamma(s); // from cmath
    prob = 1 - prob;
    return prob;
}


double CompTest::get_cell_value (const double& observed,
        const double& expected) const {
    double res = pow((observed - expected), 2.0) / expected;
    return res;
}


// https://en.wikipedia.org/wiki/Chi-squared_distribution#Cumulative_distribution_function
// s is df/2, t is chi-sq.stat/2
double CompTest::lower_incomplete_gamma_function (double s, double t) {
    if (t < 0.0)  {
        return 0.0;
    }
    double Sc = (1.0 / s);
    Sc *= pow(t, s);
    Sc *= exp(-t);
     
    double sum = 1.0;
    double numerator = 1.0;
    double denominator = 1.0;
    double stop = 1e-30;
    
    bool done = false;
    //int counter = 0;
    while (!done) {
        //counter++;
        numerator *= t;
    	s++;
    	denominator *= s;
    	double curr = (numerator / denominator);
    	sum += curr;
        if (curr < stop) { // get a better stop criterion
            done = true;
        }
    }
    //std::cerr << "LIGF run for " << counter++ << " iterations" << std::endl;
    return sum * Sc;
}
