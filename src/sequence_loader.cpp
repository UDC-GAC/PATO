#include "sequence_loader.hpp"

#include <string>
#include <iostream>
#include <algorithm>

#include <seqan/seq_io.h>
#include <seqan/stream.h>

bool load_sequences(triplex_set_t& sequences,
                    name_set_t& names,
                    const char *file_name)
{
    seqan::SeqFileIn fasta_file;

    if (!seqan::open(fasta_file, file_name)) {
        std::cerr << "PATO: error opening input file '" << file_name << "'\n";
        return false;
    }
    seqan::readRecords(names, sequences, fasta_file);

    // crop sequence name
    for (auto& name : names) {
        std::string tmp_name(name.data_begin, seqan::length(name));
        std::size_t num_chars = std::min(tmp_name.find_first_of(' '),
                                         tmp_name.size());
        name = tmp_name.substr(0, num_chars);
    }

    return sequences.size() > 0;
}
