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

#include <PATO/tfo_finder.h>

#include "finders.h"
#include "guanine_filter.h"
#include "output_writer.h"
#include "repeat_filter.h"
#include "segment_parser.h"
#include "sequence_loader.h"
#include "types.h"

namespace {

struct tfo_finder_args_t {
  pato::graph_t tc_parser;
  pato::graph_t ga_parser;
  pato::graph_t gt_parser;

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
  explicit tfo_finder_args_t(pato::motif_vector_t &motifs_,
                             pato::motif_potential_vector_t &potentials_)
      : motifs{motifs_}, potentials{potentials_},
        filter_args{motifs, block_runs, encoded_seq, true, 'G', 'N'} {}
#else
  tfo_finder_args_t()
      : filter_args{motifs, block_runs, encoded_seq, true, 'G', 'N'} {}
#endif
};

} // namespace

static void make_tfo_parsers(tfo_finder_args_t &args, unsigned max_interrupts) {
  pato::triplex_t valid_chars, invalid_chars;
  valid_chars = "TCY";
  invalid_chars = "GARN";
  pato::make_parser(args.tc_parser, valid_chars, invalid_chars, max_interrupts);

  valid_chars = "GAR";
  invalid_chars = "TCYN";
  pato::make_parser(args.ga_parser, valid_chars, invalid_chars, max_interrupts);

  valid_chars = "GTK";
  invalid_chars = "CAMN";
  pato::make_parser(args.gt_parser, valid_chars, invalid_chars, max_interrupts);
}

static void find_tfo_motifs(pato::triplex_t &sequence, unsigned id,
                            tfo_finder_args_t &args,
                            const pato::options_t &opts) {
  unsigned matches_y = 0;
  unsigned matches_r = 0;
  unsigned matches_m = 0;

  // TC motif
  if (opts.tc_motif) {
    pato::parse_segments(args.tc_parser, args.segments, sequence,
                         opts.max_interruptions, opts.min_length);

    args.filter_args.ornt = pato::orientation_t::parallel;

    for (auto &segment : args.segments) {
      pato::motif_t motif{segment, true, id, true, 'Y'};
      matches_y += filter_guanine_error_rate(motif, args.filter_args,
                                             pato::pyrimidine_motif_t{}, opts);
    }
    args.segments.clear();
  }

  // GA motif
  if (opts.ga_motif) {
    pato::parse_segments(args.ga_parser, args.segments, sequence,
                         opts.max_interruptions, opts.min_length);

    args.filter_args.ornt = pato::orientation_t::antiparallel;

    for (auto &segment : args.segments) {
      pato::motif_t motif{segment, false, id, true, 'R'};
      matches_r += filter_guanine_error_rate(motif, args.filter_args,
                                             pato::purine_motif_t{}, opts);
    }
    args.segments.clear();
  }

  // GT motif
  if (opts.gt_a_motif || opts.gt_p_motif) {
    pato::parse_segments(args.gt_parser, args.segments, sequence,
                         opts.max_interruptions, opts.min_length);

    if (opts.gt_a_motif && opts.gt_p_motif &&
        opts.run_mode != pato::run_mode_t::tfo_search) {
      args.filter_args.ornt = pato::orientation_t::both;
    } else if (opts.gt_p_motif ||
               opts.run_mode == pato::run_mode_t::tfo_search) {
      args.filter_args.ornt = pato::orientation_t::parallel;
    } else {
      args.filter_args.ornt = pato::orientation_t::antiparallel;
    }

    for (auto &segment : args.segments) {
      if ((args.filter_args.ornt == pato::orientation_t::both ||
           args.filter_args.ornt == pato::orientation_t::parallel) &&
          opts.mixed_parallel_max_guanine >= opts.min_guanine_rate) {
        pato::motif_t motif{segment, true, id, true, 'M'};
        matches_m += filter_guanine_error_rate(motif, args.filter_args,
                                               pato::mixed_motif_t{}, opts);
      }
      if ((args.filter_args.ornt == pato::orientation_t::both ||
           args.filter_args.ornt == pato::orientation_t::antiparallel) &&
          opts.mixed_antiparallel_min_guanine <= opts.max_guanine_rate) {
        pato::motif_t motif{segment, false, id, true, 'M'};
        matches_m += filter_guanine_error_rate(motif, args.filter_args,
                                               pato::mixed_motif_t{}, opts);
      }
    }
    args.segments.clear();
  }

  if (opts.run_mode == pato::run_mode_t::tfo_search) {
    pato::motif_potential_t potential{id};

    seqan::addCount(potential, matches_y, 'Y');
    seqan::addCount(potential, matches_r, 'R');
    seqan::addCount(potential, matches_m, 'M');

    seqan::setNorm(potential, seqan::length(sequence), opts.max_length,
                   opts.min_length);

    args.potentials.push_back(potential);
  }
}

void pato::find_tfo_motifs(pato::motif_vector_t &motifs,
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
    tfo_finder_args_t args{motifs, potentials};
#else
    tfo_finder_args_t args;
#endif

    if (opts.run_mode == run_mode_t::tfo_search) {
      args.filter_args.reduce_set = opts.merge_features;
    }
    make_tfo_parsers(args, opts.max_interruptions);

#pragma omp for schedule(dynamic) nowait
    for (unsigned i = 0; i < sequences.size(); ++i) {
      if (opts.filter_repeats) {
        pato::filter_repeats(args.repeats, sequences[indices[i]],
                             opts.min_repeat_length, opts.max_repeat_period);
      }
      ::find_tfo_motifs(sequences[indices[i]], indices[i], args, opts);
    }

#if defined(_OPENMP)
#pragma omp critical(motifs)
    {
      motifs.reserve(motifs.size() + args.motifs.size());
      std::move(args.motifs.begin(), args.motifs.end(),
                std::back_inserter(motifs));
    } // #pragma omp critical
    if (opts.run_mode == run_mode_t::tfo_search) {
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

pato::find_tfo_motifs_result
pato::find_tfo_motifs(const pato::options_t &opts) {
  pato::sequence_loader_t sequence_loader;
  if (!sequence_loader.init(opts.tfo_file)) {
    return pato::find_tfo_motifs_result::cannot_open_tfo_file;
  }

  auto output_writer = pato::output_writer_t::create(opts);
  if (!output_writer) {
    return pato::find_tfo_motifs_result::cannot_create_output_file;
  }

  pato::name_vector_t tfo_names;
  pato::triplex_vector_t tfo_sequences;
  sequence_loader.load_sequences(tfo_sequences, tfo_names,
                                 std::numeric_limits<unsigned>::max());
  pato::motif_vector_t tfo_motifs;
  pato::motif_potential_vector_t tfo_potentials;
  pato::find_tfo_motifs(tfo_motifs, tfo_potentials, tfo_sequences, opts);

#pragma omp parallel sections num_threads(2)
  {
#pragma omp section
    output_writer->print_motifs(tfo_motifs, tfo_names);
#pragma omp section
    output_writer->print_motifs_summary(tfo_potentials, tfo_names);
  } // #pragma omp parallel sections num_threads(2)
  output_writer->destroy();

  return pato::find_tfo_motifs_result::success;
}
