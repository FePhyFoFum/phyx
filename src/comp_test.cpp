#include <iostream>
#include <algorithm>
#include <iomanip>
#include <cmath>

using namespace std;

#include "comp_test.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"

CompTest::CompTest (istream* pios, ostream* poos, bool const& force_protein) {
    // set parameters
    is_protein_ = false;
    if (force_protein) {
        is_protein_ = true;
    }
    pios_ = pios;
    poos_ = poos;
    total_ = 0;
    read_seqs();
    return_freq_table(poos_);
    calc_chi_square();
}

void CompTest::read_seqs () {
    bool first = true;
    Sequence seq;
    string retstring;
    seqcount_ = 0;
    int ft = test_seq_filetype_stream(*pios_, retstring);
    while (read_next_seq_from_stream(*pios_, ft, retstring, seq)) {
        if (first) {
            if (!is_protein_) {
                string alpha_name = seq.get_alpha_name();
                if (alpha_name == "AA") {
                    is_protein_ = true;
                }
            }
            set_alphabet ();
            first = false;
        }
        seqcount_++;
        string temp_seq = seq.get_sequence();
        string name = seq.get_id();
        count_chars(temp_seq);
        taxon_labels_.push_back(name);
    }
    if (ft == 2) {
        seqcount_++;
        string temp_seq = seq.get_sequence();
        string name = seq.get_id();
        count_chars(temp_seq);
        taxon_labels_.push_back(name);
    }
}

// count occurrences of each valid character state in current sequence
void CompTest::count_chars (string& seq) {
    int sum = 0;
    seq = string_to_upper(seq);
    vector <int> icounts(seq_chars_.length(), 0);
        
    for (unsigned int i = 0; i < seq_chars_.length(); i++) {
        int num = count(seq.begin(), seq.end(), seq_chars_[i]);
        icounts[i] += num;
        sum += num;
        col_totals_[i] += num;
    }
    
    indiv_char_counts_.push_back(icounts);
    row_totals_.push_back(sum);
    total_ += sum;
}

// get the longest label. for printing purposes
void CompTest::get_longest_taxon_label () {
    longest_tax_label_ = 0;
    for (int i = 0; i < seqcount_; i++) {
        if ((int)taxon_labels_[i].size() > longest_tax_label_) {
            longest_tax_label_ = taxon_labels_[i].size();
        }
    }
}

// do not include gaps/ambiguous states (-, N, X)
void CompTest::set_alphabet () {
    if (is_protein_) {
        seq_chars_ = "ACDEFGHIKLMNPQRSTVWY";
    } else {
        seq_chars_ = "ACGT";
    }
    col_totals_.resize(seq_chars_.size(), 0);
}

void CompTest::return_freq_table (ostream* poos) {
    const char separator = ' ';
    const int colWidth = 10;
    // need to take into account longest_tax_label_
    get_longest_taxon_label();
    string pad = std::string(longest_tax_label_, ' ');
    // header
    (*poos) << pad << " ";
    for (unsigned int i = 0; i < seq_chars_.length(); i++) {
        (*poos) << right << setw(colWidth) << setfill(separator)
            << seq_chars_[i] << " ";
    }
    (*poos) << right << setw(colWidth) << setfill(separator) << "Nchar" << endl;
    for (int i = 0; i < seqcount_; i++) {
        int diff = longest_tax_label_ - taxon_labels_[i].size();
        (*poos_) << taxon_labels_[i];
        if (diff > 0) {
            pad = std::string(diff, ' ');
            (*poos_) << pad;
        }
        (*poos_) << " ";
        for (unsigned int j = 0; j < seq_chars_.length(); j++) {
            (*poos) << right << setw(colWidth) << setfill(separator)
                << indiv_char_counts_[i][j] << " ";
        }
        (*poos) << right << setw(colWidth) << setfill(separator) << row_totals_[i] << endl;
    }
    int diff = longest_tax_label_ - 5;
    pad = std::string(diff, ' ');
    (*poos_) << "Total" << pad << " ";
    for (unsigned int i = 0; i < col_totals_.size(); i++) {
        (*poos) << right << setw(colWidth) << setfill(separator)
            << col_totals_[i] << " ";
    }
    (*poos_) << right << setw(colWidth) << setfill(separator) << total_ << endl;
}

double CompTest::calc_chi_square () {
    test_stat_ = 0.0;
    df_ = (seqcount_ - 1) * (col_totals_.size() - 1);
    for (unsigned int i = 0; i < seqcount_; i++) {
        for (unsigned int j = 0; j < col_totals_.size(); j++) {
            double observed = (double)indiv_char_counts_[i][j];
            double expected = (double)col_totals_[j] * (double)row_totals_[i]
                / (double) total_;
            double cellv = get_cell_value(observed, expected);
            test_stat_ += cellv;
        }
    }
    cout << "Test statistic = " << test_stat_ << endl;
    cout << "DF = " << df_ << endl;
    
}

double CompTest::calc_chi_square_prob () {
    // prob given by igf(df/2, x/2) / gamma(k/2)
}

double CompTest::get_cell_value (double const& observed, double const& expected) {
    double res = pow((observed - expected), 2.0) / expected;
    return res;
}
