#ifndef PX_CONCAT_H
#define PX_CONCAT_H

#include <string>
#include <vector>

class SequenceConcatenater {
private:
    std::vector<Sequence> seqs_;
    int num_partitions_;
    int num_char_;
    int num_taxa_;
    bool toupcase_;
    std::string filename_;
    std::vector<int> partition_sizes_;
    
    void read_sequences ();
    int get_sequence_length ()const;
    
    std::vector<int> get_partition_sizes ()const;
    void delete_sequence (SequenceConcatenater& newSeqs, const int& index);

public:
    SequenceConcatenater (const bool& toupcase);
    SequenceConcatenater (std::string& seqf, const bool& toupcase);
    void concatenate (SequenceConcatenater& newSeqs);
    
    int get_num_taxa ()const;
    Sequence get_sequence (const int& index)const;
    void write_partition_information (const std::vector<std::string>& inputFiles,
        std::string& partfile);
    //~SequenceConcatenater ();
};

#endif /* PX_CONCAT_H */
