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

#include <PATO/tts_finder.h>

#include "finders.h"
#include "guanine_filter.h"
#include "output_writer.h"
#include "repeat_filter.h"
#include "segment_parser.h"
#include "sequence_loader.h"
#include "types.h"

namespace {

struct tts_finder_args_t {
  pato::graph_t plus_parser;
  pato::graph_t minus_parser;

  pato::repeat_vector_t repeats;

#if !defined(_OPENMP)
  pato::motif_vector_t &motifs;
  pato::motif_potential_vector_t &potentials;
#else
  pato::motif_vector_t motifs;
  pato::motif_potential_vector_t potentials;
#endif

  pato::char_vector_vector_t block_runs;
  pato::char_vector_vector_t encoded_seq;

  pato::segment_vector_t segments;

  pato::guanine_filter_args_t filter_args;

#if !defined(_OPENMP)
  explicit tts_finder_args_t(pato::motif_vector_t &motifs_,
                             pato::motif_potential_vector_t &potentials_)
      : motifs{motifs_}, potentials{potentials_},
        filter_args{motifs, block_runs, encoded_seq, true, 'G', 'Y'} {
    filter_args.ornt = pato::orientation_t::both;
  }
#else
  tts_finder_args_t()
      : filter_args{motifs, block_runs, encoded_seq, true, 'G', 'Y'} {
    filter_args.ornt = pato::orientation_t::both;
  }
#endif
};

} // namespace

static void make_tts_parsers(tts_finder_args_t &args, unsigned max_interrupts) {
  pato::triplex_t valid_chars, invalid_chars;
  valid_chars = "GAR";
  invalid_chars = "TCYN";
  pato::make_parser(args.plus_parser, valid_chars, invalid_chars,
                    max_interrupts);

  valid_chars = "TCY";
  invalid_chars = "GARN";
  pato::make_parser(args.minus_parser, valid_chars, invalid_chars,
                    max_interrupts);
}

static void find_tts_motifs(pato::triplex_t &sequence, unsigned id,
                            tts_finder_args_t &args,
                            const pato::options_t &opts) {
  unsigned matches_plus = 0;
  unsigned matches_minus = 0;

  // + motif
  pato::parse_segments(args.plus_parser, args.segments, sequence,
                       opts.max_interruptions, opts.min_length);

  for (auto &segment : args.segments) {
    pato::motif_t motif{segment, true, id, false, '+'};
    matches_plus +=
        filter_guanine_error_rate(motif, args.filter_args, pato::tts_t{}, opts);
  }
  args.segments.clear();

  // - motif
  pato::parse_segments(args.minus_parser, args.segments, sequence,
                       opts.max_interruptions, opts.min_length);

  for (auto &segment : args.segments) {
    pato::motif_t motif{segment, true, id, false, '-'};
    matches_minus +=
        filter_guanine_error_rate(motif, args.filter_args, pato::tts_t{}, opts);
  }
  args.segments.clear();

  if (opts.run_mode == pato::run_mode_t::tts_search) {
    pato::motif_potential_t potential{id};

    seqan::addCount(potential, matches_plus, '+');
    seqan::addCount(potential, matches_minus, '-');

    seqan::setNorm(potential, seqan::length(sequence), opts.max_length,
                   opts.min_length);

    args.potentials.push_back(potential);
  }
}

void pato::find_tts_motifs(pato::motif_vector_t &motifs,
                           pato::motif_potential_vector_t &potentials,
                           pato::triplex_vector_t &sequences,
                           const pato::options_t &opts) {
  pato::index_vector_t indices(sequences.size(), 0);
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(),
            [&](unsigned i, unsigned j) -> bool {
              return seqan::length(sequences[i]) > seqan::length(sequences[j]);
            });

#pragma omp parallel
  {
#if !defined(_OPENMP)
    tts_finder_args_t args{motifs, potentials};
#else
    tts_finder_args_t args;
#endif
    if (opts.run_mode == pato::run_mode_t::tts_search) {
      args.filter_args.reduce_set = opts.merge_features;
    }

    make_tts_parsers(args, opts.max_interruptions);

#pragma omp for schedule(dynamic) nowait
    for (unsigned i = 0; i < sequences.size(); ++i) {
      if (opts.filter_repeats) {
        pato::filter_repeats(args.repeats, sequences[indices[i]],
                             opts.min_repeat_length, opts.max_repeat_period);
      }
      ::find_tts_motifs(sequences[indices[i]], indices[i], args, opts);
    }

#if defined(_OPENMP)
#pragma omp critical
    {
      motifs.reserve(motifs.size() + args.motifs.size());
      std::move(args.motifs.begin(), args.motifs.end(),
                std::back_inserter(motifs));
    } // #pragma omp critical
    if (opts.run_mode == run_mode_t::tts_search) {
#pragma omp critical(potentials)
      {
        potentials.reserve(potentials.size() + args.potentials.size());
        std::move(args.potentials.begin(), args.potentials.end(),
                  std::back_inserter(potentials));
      } // #pragma omp critical
    }
#endif
  } // #pragma omp parallel
}

pato::find_tts_motifs_result
pato::find_tts_motifs(const pato::options_t &opts) {
  pato::sequence_loader_t sequence_loader;
  if (!sequence_loader.init(opts.tts_file)) {
    return pato::find_tts_motifs_result::cannot_open_tts_file;
  }

  auto output_writer = pato::output_writer_t::create(opts);
  if (!output_writer) {
    return pato::find_tts_motifs_result::cannot_create_output_file;
  }

  pato::name_vector_t tts_names;
  pato::triplex_vector_t tts_sequences;
  pato::motif_vector_t tts_motifs;
  pato::motif_potential_vector_t tts_potentials;

  while (true) {
    if (!sequence_loader.load_sequences(tts_sequences, tts_names,
                                        opts.chunk_size)) {
      break;
    }
    pato::find_tts_motifs(tts_motifs, tts_potentials, tts_sequences, opts);

#pragma omp parallel sections num_threads(2)
    {
#pragma omp section
      output_writer->print_motifs(tts_motifs, tts_names);
#pragma omp section
      output_writer->print_motifs_summary(tts_potentials, tts_names);
    } // #pragma omp parallel sections num_threads(2)

    tts_names.clear();
    tts_motifs.clear();
    tts_sequences.clear();
    tts_potentials.clear();
  }
  output_writer->destroy();

  return pato::find_tts_motifs_result::success;
}
