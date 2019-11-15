#ifndef _CONCAT_H_
#define _CONCAT_H_

#include <string>
#include <vector>

//#include "sequence.h"

class SequenceConcatenater {
private:
    std::vector<Sequence> seqs_;
    int num_partitions_;
    int num_char_;
    int num_taxa_;
    int ft_;
    bool toupcase_;
    bool interleave_;
    std::string filename_;
    std::vector<int> partition_sizes_;
    
    void read_sequences(std::string& seqf);
    void delete_sequence(SequenceConcatenater& newSeqs, const int& index);

public:
    SequenceConcatenater();
    SequenceConcatenater(std::string& seqf, bool& toupcase);
    void concatenate(SequenceConcatenater& newSeqs);
    int get_sequence_length();
    int get_num_taxa();
    Sequence get_sequence(const int& index);
    std::vector<int> get_partition_sizes();
    void write_partition_information(const std::vector<std::string>& inputFiles,
        std::string& partfile);
    //~SequenceConcatenater();
};

#endif /* _CONCAT_H_ */
