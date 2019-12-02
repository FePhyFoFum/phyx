#ifndef _CONCAT_H_
#define _CONCAT_H_

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
    
    void read_sequences (std::string& seqf);
    int get_sequence_length ();
    
    std::vector<int> get_partition_sizes ();
    void delete_sequence (SequenceConcatenater& newSeqs, const int& index);

public:
    SequenceConcatenater (const bool& toupcase);
    SequenceConcatenater (std::string& seqf);
    void concatenate (SequenceConcatenater& newSeqs);
    
    int get_num_taxa ();
    Sequence get_sequence (const int& index);
    void write_partition_information (const std::vector<std::string>& inputFiles,
        std::string& partfile);
    //~SequenceConcatenater ();
};

#endif /* _CONCAT_H_ */
