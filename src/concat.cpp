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
#include <cstring>

using namespace std;

#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"
#include "concat.h"

SequenceConcatenater::SequenceConcatenater (string & seqf, bool & toupcase):num_partitions_(0),
    num_char_(0), num_taxa_(0), ft_(0) {
    toupcase_ = toupcase;
    read_sequences(seqf);
}

SequenceConcatenater::SequenceConcatenater ():num_partitions_(0), num_char_(0),
    num_taxa_(0), ft_(0), interleave_(false) {
}

void SequenceConcatenater::read_sequences (string & seqf) {
    filename_ = seqf;
    string retstring;
    istream * pios = new ifstream(filename_);
    ft_ = test_seq_filetype_stream(*pios, retstring);
    Sequence seq;
    int counter = 0;
    int length = 0;
    
    // phylip (1) NEXUS (0)
    if (ft_ == 1 || ft_ == 0) {
        if (ft_ == 1) {
            vector <string> fileDim = tokenize(retstring);
            num_taxa_ = stoi(fileDim[0]);
            num_char_ = stoi(fileDim[1]);
        } else {
            get_nexus_dimensions_file(seqf, num_taxa_, num_char_, interleave_);
        }
        if (!interleave_) {
            while (read_next_seq_from_stream(*pios, ft_, retstring, seq)) {
                length = (int)seq.get_sequence().size();
                if (length != num_char_) {
                    cout << "Sequence '" << seq.get_id() << "' has " << length << " characters, but the file '"
                        << filename_ << "' specified " << num_char_ << " characters. Exiting." << endl;
                    delete pios;
                    exit(1);
                }
                if (toupcase_) {
                    seq.set_sequence(seq.seq_to_upper());
                }
                seqs_.push_back(seq);
                counter++;
            }
            if (counter != num_taxa_) {
                cout << "Read " << counter << " taxa, but the file '" << filename_ << "' specified "
                    << num_taxa_ << " taxa. Exiting." << endl;
                delete pios;
                exit(1);
            }
        } else {
            seqs_ = read_interleaved_nexus_file(seqf, num_taxa_, num_char_);
            if (toupcase_) {
                for (int i = 0; i < num_taxa_; i++) {
                    seqs_[i].set_sequence(seqs_[i].seq_to_upper());
                }
            }
        }
        
    } else if (ft_ == 2) { // fasta
        bool first = true;
        while (read_next_seq_from_stream(*pios, ft_, retstring, seq)) {
            int curr = (int)seq.get_sequence().size();
            if (!first) {
                if (curr != length) {
                    cout << "Error: current sequence has " << curr << " characters, but previous sequence had "
                        << length << " characters. Exiting." << endl;
                    delete pios;
                    exit(1);
                }
            } else {
                length = curr;
                first = false;
            }
            if (toupcase_) {
                seq.set_sequence(seq.seq_to_upper());
            }
            seqs_.push_back(seq);
            counter++;
        }
        // fasta has a trailing one
        if (toupcase_) {
            seq.set_sequence(seq.seq_to_upper());
        }
        seqs_.push_back(seq);
        counter++;
        num_taxa_ = counter;
        num_char_ = length;
    } else {
        cout << "I don't know what that alignment file format is! Exiting." << endl;
        exit(0);
    }
    num_partitions_ = 1;
    partition_sizes_.push_back(num_char_);
    delete pios;
}

// where stuff actually happens
void SequenceConcatenater::concatenate(SequenceConcatenater & newSeqs) {
    string old_filler(num_char_, '-');
    int new_seq_len = newSeqs.get_sequence_length();
    string new_filler(new_seq_len, '-');
    num_char_ += new_seq_len;
    for (int i = 0; i != num_taxa_; i++) {
        bool match_found = false;
        if (newSeqs.num_taxa_ > 0) {
            for (int j = 0; j != newSeqs.num_taxa_; j++) {
                if (seqs_[i].get_id() == newSeqs.seqs_[j].get_id()) {
                    seqs_[i].set_sequence(seqs_[i].get_sequence() + newSeqs.seqs_[j].get_sequence());
                    match_found = true;
                    // erase matched entry so it won't have to be compared against again.
                    // erase in reverse order.
                    delete_sequence(newSeqs, j);
                    break;
                }
            }
        }
        if (!match_found) { // taxon is missing from present locus.
            seqs_[i].set_sequence(seqs_[i].get_sequence() + new_filler);
        }
    }

    // now, all that should be left are the novel sequences from the new file
    if (newSeqs.num_taxa_ > 0) {
        for (int i = 0; i != newSeqs.num_taxa_; i++) {
            newSeqs.seqs_[i].set_sequence(old_filler + newSeqs.seqs_[i].get_sequence());
            seqs_.push_back(newSeqs.seqs_[i]);
            num_taxa_++;
        }
    }
    num_partitions_++;
    partition_sizes_.push_back(newSeqs.get_sequence_length());
}

int SequenceConcatenater::get_sequence_length () {
    return num_char_;
}

int SequenceConcatenater::get_num_taxa () {
    return num_taxa_;
}

void SequenceConcatenater::delete_sequence (SequenceConcatenater & newSeqs, int const& index) {
    newSeqs.seqs_.erase(newSeqs.seqs_.begin() + index);
    newSeqs.num_taxa_--;
}

Sequence SequenceConcatenater::get_sequence (int const & index) {
    return seqs_[index];
}

vector <int> SequenceConcatenater::get_partition_sizes () {
    return partition_sizes_;
}

void SequenceConcatenater::write_partition_information (vector <string> const& inputFiles,
    string & partfile) {
    ofstream outfile(partfile.c_str());
    int charIndex = 1;
    int stopIndex = 1;
    
    // need to check seq type when writing this
    // use infer_alpha / get_alpha_name
    // but: are mixed seq types allowed? prolly...
    //     - so: need to check each one
    
    for (unsigned int i = 0; i < partition_sizes_.size(); i++) {
        stopIndex = charIndex + partition_sizes_[i] - 1;
        bool going = true;
        string alpha = "";
        int j = 0;
        while (going) {
            Sequence terp = seqs_[j];
            string subseq = terp.get_sequence().substr((charIndex - 1), partition_sizes_[i]);
            // check if all are the same character (presumably all N, but useful either way)
            if (subseq.find_first_not_of(subseq.front()) != std::string::npos) {
                terp.set_sequence(subseq);
                alpha = terp.get_alpha_name();
                going = false;
            }
            j++;
        }
        outfile << alpha << ", " << inputFiles[i] << " = " << charIndex << "-" << stopIndex << endl;
        charIndex = stopIndex + 1;
    }
    outfile.close();
}
