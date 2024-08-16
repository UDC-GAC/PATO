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

#ifndef PATO_LIB_OUTPUT_WRITER_H
#define PATO_LIB_OUTPUT_WRITER_H

#include <PATO/options.h>

#include <optional>

#include "types.h"

namespace pato {

class output_writer_t {
public:
  static std::optional<output_writer_t> create(const options_t &opts);

  void print_motifs(const motif_vector_t &motifs, const name_vector_t &names);
  void print_motifs_summary(const motif_potential_vector_t &potentials,
                            const name_vector_t &names);
  void print_triplexes(
#if !defined(_OPENMP)
      const match_vector_t &matches, const motif_vector_t &tfo_motifs,
      const name_vector_t &tfo_names, const motif_vector_t &tts_motifs,
      const name_vector_t &tts_names
#else
      const match_vector_vector_t &matches, const motif_vector_t &tfo_motifs,
      const name_vector_t &tfo_names, const motif_vector_t &tts_motifs,
      const name_vector_t &tts_names
#endif
  );
  void print_triplex_summary(const potential_map_t &potentials,
                             const name_vector_t &tfo_names,
                             const name_vector_t &tts_names);

  // FIXME: Use RAII!
  void destroy();

private:
  output_writer_t(std::FILE *output_file_, std::FILE *summary_file_,
                  const options_t &opts_)
      : output_file{output_file_}, summary_file{summary_file_}, opts{opts_} {}

  std::FILE *output_file;
  std::FILE *summary_file;

  const options_t &opts;
};

} // namespace pato

#endif // PATO_LIB_OUTPUT_WRITER_H
