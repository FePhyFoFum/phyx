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


CompTest::CompTest (std::istream* pios, std::ostream* poos):num_taxa_(0), seq_length_(0),
        total_(0), df_(0), test_stat_(0.0), prob_(0.0), seq_chars_(""), alpha_name_(""),
        alpha_set_(false), datatype_set_(false), is_multi_(false), gap_('-'), missing_('?'),
        pios_(pios), poos_(poos) {
    read_in_alignment();
    count_chars();
    calc_chi_square();
}


void CompTest::read_in_alignment () {
    Sequence seq;
    std::string retstring;
    int ft = test_seq_filetype_stream(*pios_, retstring);
    int file_ntax = 0; // ntax declared in the file itself (nexus or phylip)
    std::string file_type_ = get_filetype_string(ft);
    
    // if nexus, grab metadata
    if (file_type_.compare("nexus") == 0) {
        bool is_interleaved = false;
        get_nexus_alignment_properties(*pios_, file_ntax, seq_length_,
                is_interleaved, alpha_name_, seq_chars_, gap_, missing_);
        // std::cout << "alpha_name_ = " << alpha_name_ << std::endl;
        set_datatype();
        if (is_multi_ && seq_chars_.compare("") != 0) {
            seq_chars_ += gap_;
            seq_chars_ += missing_;
            alpha_set_ = true;
        }
        retstring = ""; // have to set so seq_reader knows we are mid-file
        if (!is_interleaved) {
            while (read_next_seq_from_stream(*pios_, ft, retstring, seq)) {
                seqs_.push_back(seq);
            }
        } else {
            // need to read in everything at once
            seqs_ = read_interleaved_nexus(*pios_, file_ntax, seq_length_);
        }
    } else {
        bool complicated_phylip = false;
        // check if we are dealing with a complicated phylip format
        if (file_type_.compare("phylip") == 0) {
            get_phylip_dimensions(retstring, file_ntax, seq_length_);
            complicated_phylip = is_complicated_phylip(*pios_, seq_length_);
        }
        if (complicated_phylip) {
            seqs_ = read_phylip(*pios_, file_ntax, seq_length_);
            if (!datatype_set_) {
                alpha_name_ = seqs_[0].get_alpha_name();
                set_datatype();
            }
        } else {
            while (read_next_seq_from_stream(*pios_, ft, retstring, seq)) {
                seqs_.push_back(seq);
                if (!datatype_set_) {
                    alpha_name_ = seq.get_alpha_name();
                    //std::cout << "Setting alpha to: " << alpha_name_ << std::endl;
                    set_datatype();
                }
            }
            if (ft == 2) { // fasta has an trailing one
                seqs_.push_back(seq);
            }
        }
    }
    
    if (!is_aligned(seqs_)) {
        std::cerr << "Error: sequences must be aligned. Exiting." << std::endl;
        exit(0);
    }
    
    num_taxa_ = (int)seqs_.size();
    
    // if datatype is multi, but alphabet not set, get from entire concatenated sequence
    if (is_multi_ && !alpha_set_) {
        // grab all unique characters from the input string
        // here, seqs from all individuals are concatenated, so represents all sampled characters
        std::string concatenated = "";
        for (int i = 0; i < num_taxa_; i++) {
            concatenated += seqs_[i].get_sequence();
        }
        seq_chars_ = get_alphabet_from_sequence(concatenated);
        // remove gap and missing (if present)
        seq_chars_.erase(std::remove(seq_chars_.begin(), seq_chars_.end(), gap_), seq_chars_.end());
        seq_chars_.erase(std::remove(seq_chars_.begin(), seq_chars_.end(), missing_), seq_chars_.end());
        alpha_set_ = true;
    }
    // some error checking
    if (file_ntax != 0) {
        if (file_ntax != (int)seqs_.size()) {
            std::cerr << "Error: number of taxa declared in the file ("
                << ") does not match the number read (" << seqs_.size()
                << "). Exiting." << std::endl;
            exit(1);
        }
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
    int sum = 0;
    std::string seq = "";
    for (int i = 0; i < num_taxa_; i++) {
        sum = 0;
        seq = string_to_upper(seqs_[i].get_sequence());
        std::vector<int> icounts(seq_chars_.length(), 0);
        
        for (unsigned int i = 0; i < seq_chars_.length(); i++) {
            int num = std::count(seq.begin(), seq.end(), seq_chars_[i]);
            icounts[i] += num;
            sum += num;
            col_totals_[i] += num;
        }
        indiv_char_counts_.push_back(icounts);
        row_totals_.push_back(sum);
        total_ += sum;
    }
}


// get the longest label. for printing purposes
void CompTest::get_longest_taxon_label () {
    longest_tax_label_ = 0;
    for (int i = 0; i < num_taxa_; i++) {
        if ((int)taxon_labels_[i].size() > longest_tax_label_) {
            longest_tax_label_ = taxon_labels_[i].size();
        }
    }
}


void CompTest::return_freq_table () {
    const char separator = ' ';
    const int colWidth = 8;
    // need to take into account longest_tax_label_
    get_longest_taxon_label();
    std::string pad = std::string(longest_tax_label_, ' ');
    // header
    (*poos_) << "Observed character counts:" << std::endl;
    (*poos_) << pad << " ";
    for (unsigned int i = 0; i < seq_chars_.length(); i++) {
        (*poos_) << std::right << std::setw(colWidth) << std::setfill(separator)
            << seq_chars_[i] << " ";
    }
    (*poos_) << std::right << std::setw(colWidth) << std::setfill(separator) << "Nchar" << std::endl;
    for (int i = 0; i < num_taxa_; i++) {
        int diff = longest_tax_label_ - taxon_labels_[i].size();
        (*poos_) << taxon_labels_[i];
        if (diff > 0) {
            pad = std::string(diff, ' ');
            (*poos_) << pad;
        }
        (*poos_) << " ";
        for (unsigned int j = 0; j < seq_chars_.length(); j++) {
            (*poos_) << std::right << std::setw(colWidth) << std::setfill(separator)
                << indiv_char_counts_[i][j] << " ";
        }
        (*poos_) << std::right << std::setw(colWidth) << std::setfill(separator) << row_totals_[i] << std::endl;
    }
    int diff = longest_tax_label_ - 5;
    pad = std::string(diff, ' ');
    (*poos_) << "Total" << pad << " ";
    for (unsigned int i = 0; i < col_totals_.size(); i++) {
        (*poos_) << std::right << std::setw(colWidth) << std::setfill(separator)
            << col_totals_[i] << " ";
    }
    (*poos_) << std::right << std::setw(colWidth) << std::setfill(separator) << total_ << std::endl;
} 


void CompTest::calc_chi_square () {
    test_stat_ = 0.0;
    df_ = (num_taxa_ - 1) * (col_totals_.size() - 1);
    for (int i = 0; i < num_taxa_; i++) {
        for (unsigned int j = 0; j < col_totals_.size(); j++) {
            double observed = (double)indiv_char_counts_[i][j];
            double expected = (double)col_totals_[j] * (double)row_totals_[i]
                / (double) total_;
            double cellv = get_cell_value(observed, expected);
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
    double prob = 0.0;
    double s = df / 2;
    double t = xstat / 2;
    prob = lower_incomplete_gamma_function(s, t);
    prob /= tgamma(s); // from cmath
    prob = 1 - prob;
    return prob;
}


double CompTest::get_cell_value (const double& observed, const double& expected) {
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
    double curr = 0.0;
    double stop = 1e-30;
    
    bool done = false;
    int counter = 0;
    while (!done) {
        counter++;
        numerator *= t;
    	s++;
    	denominator *= s;
    	curr = (numerator / denominator);
    	sum += curr;
    	//std::cout << counter << ". Adding (Nom / Denom) = " << terp << "; Sum = " << Sum << std::endl;
        if (curr < stop) { // get a better stop criterion
            done = true;
        }
    }
    //std::cerr << "LIGF run for " << counter++ << " iterations" << std::endl;
    return sum * Sc;
}
