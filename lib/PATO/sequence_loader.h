// MIT License
//
// Copyright (c) 2022-onwards Iñaki Amatria-Barral
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef PATO_LIB_SEQUENCE_LOADER_H
#define PATO_LIB_SEQUENCE_LOADER_H

#include <PATO/options.h>
#include <seqan/seq_io.h>

#include <optional>

#include "types.h"

namespace pato {

class sequence_loader_t {
public:
  static std::optional<sequence_loader_t>
  create(const seqan::CharString &file_name);

  bool load_sequences(triplex_vector_t &sequences, name_vector_t &names,
                      unsigned num_sequences);

private:
  sequence_loader_t(seqan::SeqFileIn *fasta_file_) : fasta_file{fasta_file_} {}

  std::shared_ptr<seqan::SeqFileIn> fasta_file;
};

} // namespace pato

#endif // PATO_LIB_SEQUENCE_LOADER_H
