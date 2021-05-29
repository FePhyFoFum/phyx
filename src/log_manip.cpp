#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <numeric>

#include "log_manip.h"
#include "utils.h"


LogManipulator::LogManipulator (const std::string& logtype, const std::vector<std::string>& input_files,
    std::ostream* poos, const bool& verbose):burnin_(0), nthin_(0), nrandom_(0),
    seed_(0), count_(false), ntotal_samples_(0), num_cols_(0), num_cols_retain_(0),
    files_(input_files) {
    num_files_ = static_cast<int>(input_files.size());
    logtype_ = logtype;
    poos_ = poos;
    verbose_ = verbose;
}


void LogManipulator::sample (const int& burnin, const int& nthin, const int& nrandom,
    const long int& seed) {
    burnin_ = burnin;
    nthin_= nthin;
    nrandom_ = nrandom;
    seed_ = seed;
    if (logtype_ == "parameter") {
        sample_parameters();
    } else {
        sample_trees();
    }
}


void LogManipulator::summarize (const int& burnin, const int& nthin) {
    // need to:
    // 1. read in and store values. may involve burnin/thinning/multiple values
    // 2. for each, calculate ECC, ACT, mean, median, SD, etc.
    burnin_ = burnin;
    nthin_= nthin;
    
    collect_parameter_samples ();
    return_statistics_table ();
}


void LogManipulator::return_statistics_table () {
    // summarize shit
    std::vector<std::string> stat_names{"Mean", "Median", "Variance", "ESS", "ACT", "N"};
    double mean = 0.0;
    double variance = 0.0;
    double median = 0.0;
    double ESS = 0.0;
    double ACT = 0.0;
    // n should be constant (i.e. a rectangular log)
    auto nsamples = static_cast<int>(parm_samples_[0].size());
    // need interval between samples i.e. how thinned
    // first column in parm log is the state number
    // converted to int since all read in as doubles
    // note this assumes constant thinning. should be safe
    // could only possibly be a problem when combining logs from analyses with different sampling settings
    auto step_size = static_cast<int>(parm_samples_[0][1] - parm_samples_[0][0]);
    
    // need to find the longest parameter name to set width of first column
    const int colWidth = 12; // this value works nice with decimal and sci. notation
    std::cout.precision(6);
    //(*poos_) << std::fixed;
    int longest_label = get_longest_label(parm_names_);
    std::string pad = std::string(longest_label, ' ');
    // header
    (*poos_) << pad << " ";
    for (const auto & stat_name : stat_names) {
        (*poos_) << std::right << std::setw(colWidth) << stat_name << " ";
    }
    (*poos_) << std::endl;
    for (unsigned long i = 1; i < static_cast<unsigned long>(num_cols_); i++) {
        calculate_summary_statistics (parm_samples_[i], mean, variance, median,
            ESS, ACT, nsamples, step_size);
        int diff = longest_label - static_cast<int>(parm_names_[i].size());
        (*poos_) << parm_names_[i];
        if (diff > 0) {
            pad = std::string(diff, ' ');
            (*poos_) << pad;
        }
        (*poos_) << " ";
        (*poos_) << std::right << std::setw(colWidth) << mean << " ";
        (*poos_) << std::right << std::setw(colWidth) << median << " ";
        (*poos_) << std::right << std::setw(colWidth) << variance << " ";
        (*poos_) << std::right << std::setw(colWidth) << ESS << " ";
        (*poos_) << std::right << std::setw(colWidth) << ACT << " ";
        (*poos_) << std::right << std::setw(colWidth) << nsamples << std::endl;
    }
}


/*
want things like:
1. mean - done
1.2. geometric mean? (only positive numbers)
1.3. harmonic mean?
1.4. std. error of mean? - nah
2. std. dev. and/or variance - done
3. median - done
4. range = min, max <- trivial to add, but desired?
5. 95% HPD range
6. ESS - done
7. ACT - done
8. n - done
*/
// since each parameter is independent, could use multiple threads to speed things up...
void LogManipulator::calculate_summary_statistics (std::vector<double>& vals,
        double& mean, double& variance, double& median, double& ESS, double& ACT,
        const int& n_samples, const int& step_size) {
    // leave median for last, since it alters the ordering, which is needed for ACT
    v_mean_variance(vals, mean, variance);
    
    //-----------------------------------------------------------//
    // effective sample size (ESS) and autocorrelation time (ACT)
    // NOTE: there does not appear to be A correct way to compute ESS:
    //  https://stats.stackexchange.com/a/441628/263112
    // the (preliminary) stuff below uses the BEAST-flavour, in part because values can be validated
    // indeed, at present this is just C++-ified code from the java source:
    // https://tinyurl.com/847bpbpm
    
    // maximum lag to consider for ACT. may never be reached (i.e., smart early exit)
    // increasing MAXIMUM_LAG greatly increases computation time (see loops below)
    // it is not clear to me that 2000 is a reasonable value...
    //      - setting to 5000 gives different ESS & ACT in a predictable direction
    //          - i.e. ACT goes up, and ESS (obviously) goes down
    //          - this obviously only occurs when early exit is not triggered
    //              - although it did so with the first random log i cam across...
    // anyway, keep it like this for now to make sure identical values are computed
    
    int MAXIMUM_LAG = 2000;
    unsigned long max_lag = static_cast<unsigned long>(std::min(n_samples - 1, MAXIMUM_LAG));
    std::vector<double> gamma_stat(max_lag);
    
    double var_stat = 0.0;
    double del1 = 0.0;
    double del2 = 0.0;

    for (unsigned long i = 0; i < max_lag; i++) {
        for (unsigned long j = 0; j < n_samples - i; j++) {
            del1 = vals[j] - mean;
            del2 = vals[j + i] - mean;
            gamma_stat[i] += (del1 * del2);
        }

        gamma_stat[i] /= static_cast<double>(n_samples - i);

        if (i == 0) {
            var_stat = gamma_stat[0];
        } else if (i % 2 == 0) {
            if (gamma_stat[i - 1] + gamma_stat[i] > 0) {
                var_stat += 2.0 * (gamma_stat[i - 1] + gamma_stat[i]);
            } else {
                //std::cout << "quitting at i = " << i << std::endl;
                // quit since we've gone past where samples are autocorrelated
                max_lag = i;
            }
        }
    }

    // ACT
    if (essentially_equal(gamma_stat[0], 0.0)) {
        ACT = 0;
    } else {
        ACT = step_size * var_stat / gamma_stat[0];
    }

    // finally, ESS
    if (essentially_equal(ACT, 0.0)) {
        ESS = 1;
    } else {
        ESS = (step_size * n_samples) / ACT;
    }
    
    // now we can compute the median (partially rearranges order
    median = v_median(vals);
}


// like sample_parameters below, but store retained samples/values
void LogManipulator::collect_parameter_samples () {
    if (!files_.empty()) { // hrm this should not be necessary...
        bool first_entry = true; // use for initialization
        std::vector <double> terp; // a vector to reuse throughout
        ntotal_samples_ = 0;
        for (unsigned long i = 0; i < static_cast<unsigned long>(num_files_); i++) {
            std::string curfile = files_[i];
            infilestr_.open(curfile.c_str());
            std::string line;
            bool first_line = true;
            int par_counter = 0; // this is the raw number of parameter lines in a file
            int sample_counter = 0;
            while (getline_safe(infilestr_, line)) {
                if (line.empty() || check_comment_line(line)) {
                    continue;
                }
                if (first_line) {
                    std::vector<std::string> header = tokenize(line);
                    auto curpars = static_cast<int>(header.size());
                    if (i == 0) { // first header
                        num_cols_ = curpars;
                        parm_names_ = header;
                        if (first_entry) {
                            // let us initialize the sample collector
                            unsigned long n_expected = 10000000; // purposely an overestimate to avoid reallocation
                            for (unsigned long j = 0; j < static_cast<unsigned long>(num_cols_); j++) {
                                parm_samples_.push_back(terp);
                                parm_samples_[j].reserve(n_expected);
                            }
                            first_entry = false;
                        }
                    } else {
                        // check that we've still got the same number of parameters i.e. files match
                        if (curpars != num_cols_) {
                            std::cerr << "Error: number of parameters in file " << (i + 1)
                                << "(" << curpars << ") does not match that from first file ("
                                << num_cols_ << "). Exiting." << std::endl;
                            exit(0);
                        } else if (header != parm_names_) {
                            // check that headers are identical
                            std::cerr << "Error: header for file " << (i + 1)
                                << "does not match that from first file. Exiting." << std::endl;
                            exit(0);
                        }
                    }
                    first_line = false;
                    continue;
                }
                if ((par_counter - burnin_) > 0 && (par_counter - burnin_) < nthin_) {
                    // skip because does not match sampling parameters
                    par_counter++;
                    continue;
                }
                if ((par_counter - burnin_) == 0) {
                    // keep first post-burnin sample from a file
                    par_counter++;
                    store_sample(line);
                    sample_counter++;
                    ntotal_samples_++;
                    continue;
                }
                if ((par_counter - burnin_) > 0 && (par_counter - burnin_) % nthin_ == 0) {
                    par_counter++;
                    store_sample(line);
                    sample_counter++;
                    ntotal_samples_++;
                    continue;
                }
                // skip because have not yet exceeded burnin
                par_counter++;
            }
            indiv_raw_counts_.push_back(par_counter);
            indiv_sample_totals_.push_back(sample_counter);
            infilestr_.close();
        }
        if (verbose_) {
            for (unsigned long i = 0; i < static_cast<unsigned long>(num_files_); i++) {
                std::cout << files_[i] << ": " << indiv_sample_totals_[i]
                    << " samples retained (from original " << indiv_raw_counts_[i]
                    << " samples) for " << (num_cols_ - 1) << " variables." << std::endl;
            }
            std::cout << "Retained " << ntotal_samples_ << " total samples and " << (num_cols_ - 1)
                << " variables across " << num_files_ << " input files." << std::endl;
        }
        if (ntotal_samples_ == 0) {
            std::cerr << "Error: no samples remain to summarize. Check burnin/thinning. Exiting." << std::endl;
        }
    }
}


// convenience function to avoid code duplication
void LogManipulator::store_sample (std::string const& line) {
    // first, tokenize
    std::vector<std::string> sample = tokenize(line);
    // convert to double. note that this include the first entry (state) which is an int
    std::vector<double> converted = string_v_to_double_v(sample);
    // finally, store these puppies in the initialized container parm_samples_
    for (unsigned long i = 0; i < static_cast<unsigned long>(num_cols_); i++) {
        parm_samples_[i].push_back(converted[i]);
    }
}


void LogManipulator::count () {
    if (logtype_ == "parameter") {
        count_parameter_samples();
    } else {
        count_tree_samples();
    }
}


// TODO: should counting allow burnin/thinning? nah.
void LogManipulator::count_parameter_samples () {
    num_cols_ = 0;
    if (!files_.empty()) {
        for (int i = 0; i < num_files_; i++) {
            std::string curfile = files_[static_cast<unsigned long>(i)];
            infilestr_.open(curfile.c_str());
            std::string line;
            bool first_line = true;
            int num_samps = 0;
            while (getline_safe(infilestr_, line)) {
                if (line.empty() || check_comment_line(line)) {
                    continue;
                }
                if (first_line) {
                    std::vector<std::string> header = tokenize(line);
                    auto curpars = static_cast<int>(header.size());
                    if (i == 0) { // first header
                        num_cols_ = curpars;
                        parm_names_ = header;
                    } else {
                        // check that we've still got the same number of parameters i.e. files match
                        if (curpars != num_cols_) {
                            std::cerr << "Error: number of parameters in file " << (i + 1)
                                << "(" << curpars << ") does not match that from first file ("
                                << num_cols_ << "). Exiting." << std::endl;
                            exit(0);
                        } else if (header != parm_names_) {
                            // check that headers are identical
                            std::cerr << "Error: header for file " << (i + 1)
                                << "does not match that from first file. Exiting." << std::endl;
                            exit(0);
                        }
                    }
                    first_line = false;
                    continue;
                }
                num_samps++;
                continue;
            }
            indiv_sample_totals_.push_back(num_samps);
            infilestr_.close();
        }
        //ntotal_samples_ = accumulate(indiv_totals_.begin(), indiv_totals_.end(), 0);
        //(*poos_) << "Counted " << ntotal_samples_ << " total samples and " << (num_cols_ - 1)
        //    << " variables across " << num_files_ << " files." << std::endl;
    } else {
        
        // stream stuff will go here (maybe)
        
    }
    ntotal_samples_ = std::accumulate(indiv_sample_totals_.begin(), indiv_sample_totals_.end(), 0);
}


void LogManipulator::count_tree_samples () {
    if (!files_.empty()) {
        for (int i = 0; i < num_files_; i++) {
            std::string curfile = files_[static_cast<unsigned long>(i)];
            infilestr_.open(curfile.c_str());
            std::string line;
            int num_samps = 0;
            while (getline_safe(infilestr_, line)) {
                if (line.empty() || check_comment_line(line)) {
                    continue;
                }
                std::vector<std::string> tokens = tokenize(line);
                std::string first = tokens[0];
                std::transform(first.begin(), first.end(), first.begin(), ::tolower);
                if (first == "tree") {
                    num_samps++;
                }
            }
            indiv_sample_totals_.push_back(num_samps);
            infilestr_.close();
        }
        //ntotal_samples_ = accumulate(indiv_totals_.begin(), indiv_totals_.end(), 0);
        //(*poos_) << "Counted " << ntotal_samples_ << " total samples and " << (num_cols_ - 1)
        //    << " variables across " << num_files_ << " files." << std::endl;
    } else {
        
        // stream stuff will go here (maybe)
        
    }
    ntotal_samples_ = std::accumulate(indiv_sample_totals_.begin(), indiv_sample_totals_.end(), 0);
}


void LogManipulator::get_sample_counts () {
    if (logtype_ == "parameter") {
        for (unsigned long i = 0; i < static_cast<unsigned long>(num_files_); i++) {
            (*poos_) << files_[i] << ": " << indiv_sample_totals_[i] << " samples of "
                << (num_cols_ - 1) << " variables." << std::endl;
        }
        if (num_files_ > 1) {
            (*poos_) << "Counted " << ntotal_samples_ << " total samples of "
                << (num_cols_ - 1) << " variables across " << num_files_ << " files." << std::endl;
        }
    } else {
        for (unsigned long i = 0; i < static_cast<unsigned long>(num_files_); i++) {
            (*poos_) << files_[i] << ": " << indiv_sample_totals_[i] << " trees." << std::endl;
        }
        if (num_files_ > 1) {
            (*poos_) << "Counted " << ntotal_samples_ << " total tree samples across "
                << num_files_ << " files." << std::endl;
        }
    }
}


void LogManipulator::get_column_names () {
    if (!files_.empty()) {
        if (num_files_ == 1) {
            std::string curfile = files_[0];
            infilestr_.open(curfile.c_str());
            std::string line;
            bool first_line = true;
            while (getline_safe(infilestr_, line)) {
                if (line.empty() || check_comment_line(line)) {
                    continue;
                }
                if (first_line) {
                    std::vector<std::string> header = tokenize(line);
                    auto curpars = static_cast<int>(header.size());
                    num_cols_ = curpars;
                    parm_names_ = header;
                    first_line = false;
                    break;
                }
            }
            for (int i = 0; i < num_cols_; i++) {
                (*poos_) << i+1 << ". " << parm_names_[static_cast<unsigned long>(i)] << std::endl;
            }
        } else {
        // multiple files. make sure number/name/order of columns is identical
            
        }
    }
}


// delete variables (columns) from a large parameter log
// assumes a single log file
// assume that burnin etc. is done separately?
void LogManipulator::delete_columns (const std::vector<int>& col_ids) {
    if (!files_.empty()) {
        ntotal_samples_ = 0;
        std::string curfile = files_[0];
        infilestr_.open(curfile.c_str());
        std::string line;
        bool first_line = true;
        int sample_counter = 0;
        std::vector<int> cols_to_retain;
        
        while (getline_safe(infilestr_, line)) {
            if (line.empty() || check_comment_line(line)) {
                continue;
            }
            if (first_line) {
                std::vector<std::string> header = tokenize(line);
                auto curpars = static_cast<int>(header.size());
                num_cols_ = curpars;
                parm_names_ = header;
                
                cols_to_retain.resize(static_cast<unsigned long>(num_cols_));
                
                // check that right end of col_ids is valid (0 already checked upstream)
                // vector has been sorted
                if (col_ids.back() > num_cols_) {
                    std::cerr << "Warning: column numbers are 1-indexed. Exiting." << std::endl;
                    exit(0);
                }
                
                std::iota(cols_to_retain.begin(), cols_to_retain.end(), 0);
                
                // remove unwanted column indices (reverse order)
                for (auto i = col_ids.rbegin(); i < col_ids.rend(); i++) {
                    // subtract 1 because input is 1-indexed
                    cols_to_retain.erase(cols_to_retain.begin()+(*i)-1);
                }
                
                num_cols_retain_ = static_cast<int>(cols_to_retain.size());
                for (int i = 0; i < num_cols_retain_; i++) {
                    (*poos_) << header[cols_to_retain[static_cast<unsigned long>(i)]];
                    if (i < (num_cols_retain_ - 1)) {
                        (*poos_) << "\t";
                    }
                }
                (*poos_) << std::endl;
                first_line = false;
                continue;
            }
            std::vector<std::string> samp = tokenize(line);
            for (int i = 0; i < num_cols_retain_; i++) {
                (*poos_) << samp[cols_to_retain[static_cast<unsigned long>(i)]];
                if (i < (num_cols_retain_ - 1)) {
                    (*poos_) << "\t";
                }
            }
            (*poos_) << std::endl;
            sample_counter++;
        }
            
        if (verbose_) {
            std::cout << "Retained " << sample_counter << " samples for "
                << num_cols_retain_ << " variables." << std::endl;
        }
    }
}


// complement of above: retain columns passed in
void LogManipulator::retain_columns (const std::vector<int>& col_ids) {
    if (!files_.empty()) {
        ntotal_samples_ = 0;
        std::string curfile = files_[0];
        infilestr_.open(curfile.c_str());
        std::string line;
        bool first_line = true;
        int sample_counter = 0;
        std::vector<int> cols_to_retain;
        for (int col_id : col_ids) {
            // subtract 1 because input is 1-indexed
            cols_to_retain.push_back(col_id - 1);
        }
        num_cols_retain_ = static_cast<int>(cols_to_retain.size());
        
        while (getline_safe(infilestr_, line)) {
            if (line.empty() || check_comment_line(line)) {
                continue;
            }
            if (first_line) {
                std::vector<std::string> header = tokenize(line);
                auto curpars = static_cast<int>(header.size());
                num_cols_ = curpars;
                
                // check that right end of col_ids is valid (0 already checked upstream)
                // vector has been sorted
                if (col_ids.back() > num_cols_) {
                    std::cerr << "Warning: column numbers are 1-indexed. Exiting." << std::endl;
                    exit(0);
                }
                
                parm_names_ = header;
                for (int i = 0; i < num_cols_retain_; i++) {
                    (*poos_) << header[cols_to_retain[static_cast<unsigned long>(i)]];
                    if (i < (num_cols_retain_ - 1)) {
                        (*poos_) << "\t";
                    }
                }
                (*poos_) << std::endl;
                first_line = false;
                continue;
            }
            std::vector<std::string> samp = tokenize(line);
            for (int i = 0; i < num_cols_retain_; i++) {
                (*poos_) << samp[cols_to_retain[static_cast<unsigned long>(i)]];
                if (i < (num_cols_retain_ - 1)) {
                    (*poos_) << "\t";
                }
            }
            (*poos_) << std::endl;
            sample_counter++;
        }
            
        if (verbose_) {
            std::cout << "Retained " << sample_counter << " samples for "
                << num_cols_retain_ << " variables." << std::endl;
        }
    }
}


void LogManipulator::sample_parameters () {
    if (!files_.empty()) { // hrm this should not be necessary...
        ntotal_samples_ = 0;
        for (int i = 0; i < num_files_; i++) {
            std::string curfile = files_[static_cast<unsigned long>(i)];
            infilestr_.open(curfile.c_str());
            std::string line;
            bool first_line = true;
            int par_counter = 0; // this is the raw number of parameter lines in a file
            int sample_counter = 0;
            while (getline_safe(infilestr_, line)) {
                if (line.empty() || check_comment_line(line)) {
                    continue;
                }
                if (first_line) {
                    std::vector<std::string> header = tokenize(line);
                    auto curpars = static_cast<int>(header.size());
                    if (i == 0) { // first header
                        num_cols_ = curpars;
                        parm_names_ = header;
                        for (int j = 0; j < num_cols_; j++) {
                            (*poos_) << parm_names_[static_cast<unsigned long>(j)];
                            if (j < (num_cols_ - 1)) {
                                (*poos_) << "\t";
                            }
                        }
                        (*poos_) << std::endl;
                    } else {
                        // check that we've still got the same number of parameters i.e. files match
                        if (curpars != num_cols_) {
                            std::cerr << "Error: number of parameters in file " << (i + 1)
                                << "(" << curpars << ") does not match that from first file ("
                                << num_cols_ << "). Exiting." << std::endl;
                            exit(0);
                        } else if (header != parm_names_) {
                            // check that headers are identical
                            std::cerr << "Error: header for file " << (i + 1)
                                << "does not match that from first file. Exiting." << std::endl;
                            exit(0);
                        }
                    }
                    first_line = false;
                    continue;
                }
                if ((par_counter - burnin_) > 0 && (par_counter - burnin_) < nthin_) {
                    // skip because does not match sampling parameters
                    par_counter++;
                    continue;
                }
                if ((par_counter - burnin_) == 0) {
                    // keep first post-burnin sample from a file
                    par_counter++;
                    write_reformatted_sample(line, ntotal_samples_);
                    sample_counter++;
                    ntotal_samples_++;
                    continue;
                }
                if ((par_counter - burnin_) > 0 && (par_counter - burnin_) % nthin_ == 0) {
                    par_counter++;
                    write_reformatted_sample(line, ntotal_samples_);
                    sample_counter++;
                    ntotal_samples_++;
                    continue;
                }
                // skip because have not yet exceeded burnin
                par_counter++;
            }
            indiv_raw_counts_.push_back(par_counter);
            indiv_sample_totals_.push_back(sample_counter);
            infilestr_.close();
        }
        if (verbose_) {
            for (unsigned long i = 0; i < static_cast<unsigned long>(num_files_); i++) {
                std::cout << files_[i] << ": " << indiv_sample_totals_[i]
                    << " samples retained (from original " << indiv_raw_counts_[i]
                    << " samples) for " << (num_cols_ - 1) << " variables." << std::endl;
            }
            std::cout << "Retained " << ntotal_samples_ << " total samples and " << (num_cols_ - 1)
                << " variables across " << num_files_ << " input files." << std::endl;
        }
    }
}


void LogManipulator::sample_trees () {
    if (!files_.empty()) {
        ntotal_samples_ = 0;
        for (int i = 0; i < num_files_; i++) {
            std::string curfile = files_[static_cast<unsigned long>(i)];
            infilestr_.open(curfile.c_str());
            std::string line;
            int tree_counter = 0; // this is the raw number of tree lines in a file
            int sample_counter = 0;
            bool trees_encountered = false;
            while (getline_safe(infilestr_, line)) {
                if (line.empty() || check_comment_line(line)) {
                    if (i == 0) {
                        // keep comment information from _top_ of first file
                        if (!trees_encountered) {
                            (*poos_) << line << std::endl;
                        }
                    }
                    continue;
                }
                std::vector<std::string> tokens = tokenize(line);
                std::string first = tokens[0];
                std::transform(first.begin(), first.end(), first.begin(), ::tolower);
                if (first == "tree") {
                    if (tree_counter == 0) {
                        trees_encountered = true;
                    }
                    if (ntotal_samples_ == 0) {
                        // grab tree naming scheme
                        get_tree_name_prefix(line);
                    }
                    if ((tree_counter - burnin_) > 0 && (tree_counter - burnin_) < nthin_) {
                        // skip because does not match sampling parameters
                        tree_counter++;
                        continue;
                    }
                    if ((tree_counter - burnin_) == 0) {
                        // keep first post-burnin sample from a file
                        tree_counter++;
                        write_reformatted_sample(line, ntotal_samples_);
                        sample_counter++;
                        ntotal_samples_++;
                        continue;
                    }
                    if ((tree_counter - burnin_) > 0 && (tree_counter - burnin_) % nthin_ == 0) {
                        tree_counter++;
                        write_reformatted_sample(line, ntotal_samples_);
                        sample_counter++;
                        ntotal_samples_++;
                        continue;
                    }
                    // skip because have not yet exceeded burnin
                    tree_counter++;
                } else { // not a tree line. only care about first file here. don't want anything below trees
                    // keep header from first file
                    // likely includes translation table
                    if (i == 0 && !trees_encountered) {
                        (*poos_) << line << std::endl;
                    }
                }
            }
            indiv_raw_counts_.push_back(tree_counter);
            indiv_sample_totals_.push_back(sample_counter);
            infilestr_.close();
        }
        (*poos_) << "End;" << std::endl;
        
        if (verbose_) {
            for (unsigned long i = 0; i < static_cast<unsigned long>(num_files_); i++) {
                std::cerr << files_[i] << ": " << indiv_sample_totals_[i]
                    << " tree samples retained (from original " << indiv_raw_counts_[i]
                    << " samples)." << std::endl;
            }
            std::cerr << "Retained " << ntotal_samples_ << " total tree samples across "
                << num_files_ << " input files." << std::endl;
        }
    }
}


// gen.NNN from Mrbayes, STATE_NNN from BEAST
// try to be general here, though: take whatever precedes the sample number
void LogManipulator::get_tree_name_prefix (std::string& sample) {
    std::vector<std::string> terp = tokenize(sample);
    std::string tree_name = terp[1];
    std::size_t found = tree_name.find_first_of("0123456789");
    tree_name.replace(tree_name.begin(), tree_name.end(), tree_name.begin(), tree_name.begin()+found);
    tree_name_prefix_ = tree_name;
}


void LogManipulator::write_reformatted_sample (std::string& sample, int& sample_num) {
    std::vector<std::string> terp = tokenize(sample);
    if (logtype_ == "parameter") {
        (*poos_) << sample_num;
        for (int i = 1; i < num_cols_; i++) {
            (*poos_) << "\t" << terp[static_cast<unsigned long>(i)];
        }
        (*poos_) << std::endl;
    } else {
        // format should be: tree tree_name lots_of_other_optional_things
        (*poos_) << "tree " << tree_name_prefix_ << sample_num;
        for (unsigned int i = 2; i < terp.size(); i++) {
            (*poos_) << " " << terp[i];
        }
        (*poos_) << std::endl;
    }
}
