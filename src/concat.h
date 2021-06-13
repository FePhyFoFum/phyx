#ifndef PX_CONCAT_H
#define PX_CONCAT_H

#include <string>
#include <vector>

class SequenceConcatenater {
private:
    std::vector<Sequence> seqs_;
    unsigned int num_partitions_;
    unsigned int num_char_;
    unsigned int num_taxa_;
    bool toupcase_;
    std::string filename_;
    std::vector<unsigned int> partition_sizes_;
    
    void read_sequences ();
    unsigned int get_sequence_length () const;
    
    //std::vector<unsigned int> get_partition_sizes ()const; // not currently used
    void delete_sequence (SequenceConcatenater& newSeqs, const unsigned int& index);

public:
    explicit SequenceConcatenater (const bool& toupcase);
    SequenceConcatenater (std::string& seqf, const bool& toupcase);
    void concatenate (SequenceConcatenater& newSeqs);
    
    unsigned int get_num_taxa () const;
    Sequence get_sequence (const unsigned int& index) const;
    void write_partition_information (const std::vector<std::string>& inputFiles,
        std::string& partfile);
    //~SequenceConcatenater ();
};

#endif /* PX_CONCAT_H */
