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

#include "sequence_loader.h"

std::optional<pato::sequence_loader_t>
pato::sequence_loader_t::create(const seqan::CharString &file_name) {
  seqan::SeqFileIn *fasta_file{new seqan::SeqFileIn{}};
  if (seqan::open(*fasta_file, seqan::toCString(file_name))) {
    return pato::sequence_loader_t{fasta_file};
  }
  return std::nullopt;
}

static void crop_sequence_name(seqan::CharString &name) {
  std::string tmp_name{name.data_begin, seqan::length(name)};
  std::size_t num_chars =
      std::min(tmp_name.find_first_of(' '), tmp_name.size());
  name = tmp_name.substr(0, num_chars);
}

bool pato::sequence_loader_t::load_sequences(pato::triplex_vector_t &sequences,
                                             pato::name_vector_t &names,
                                             unsigned num_sequences) {
  seqan::readRecords(names, sequences, *fasta_file, num_sequences);
  for (seqan::CharString &name : names) {
    crop_sequence_name(name);
  }
  return !sequences.empty();
}
