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

#include <PATO/options.h>

#include "types.h"

#ifndef PATO_LIB_GUANINE_FILTER_H
#define PATO_LIB_GUANINE_FILTER_H

namespace pato {

enum class orientation_t : int {
  antiparallel = -1,
  both = 0,
  parallel = 1,
};

struct guanine_filter_args_t {
  motif_vector_t &motifs;

  char_vector_vector_t &block_runs;
  char_vector_vector_t &encoded_seq;

  bool reduce_set;

  char filter_char;
  char interrupt_char;

  orientation_t ornt;

  guanine_filter_args_t(motif_vector_t &motifs_,
                        char_vector_vector_t &block_runs_,
                        char_vector_vector_t &encoded_seq_, bool reduce_set_,
                        char filter_char_, char interrupt_char_)
      : motifs{motifs_}, block_runs{block_runs_}, encoded_seq{encoded_seq_},
        reduce_set{reduce_set_}, filter_char{filter_char_},
        interrupt_char{interrupt_char_} {}
};

using tts_t = seqan::Tag<struct _tts>;
using tfo_t = seqan::Tag<struct _tfo>;
using mixed_motif_t = seqan::Tag<struct _mixed_motif>;
using purine_motif_t = seqan::Tag<struct _purine_motif>;
using pyrimidine_motif_t = seqan::Tag<struct _pyrimidine_motif>;

template <typename tag_t>
unsigned filter_guanine_error_rate(motif_t &motif, guanine_filter_args_t &args,
                                   const tag_t &tag, const options_t &opts);

} // namespace pato

#endif // PATO_LIB_GUANINE_FILTER_H
