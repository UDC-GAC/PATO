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

#include <PATO/tpx_finder.h>

#include "finders.h"
#include "guanine_filter.h"
#include "output_writer.h"
#include "segment_parser.h"
#include "sequence_loader.h"
#include "types.h"

namespace {

struct tpx_finder_args_t {
  pato::graph_t tpx_parser;

  pato::motif_vector_t tpx_motifs;

  pato::char_vector_vector_t block_runs;
  pato::char_vector_vector_t encoded_seq;

  pato::segment_vector_t segments;

  pato::match_vector_t &matches;

#if !defined(_OPENMP)
  pato::potential_map_t &potentials;
#else
  pato::potential_map_t potentials;
#endif

  pato::guanine_filter_args_t filter_args;

  int min_score;

#if !defined(_OPENMP)
  tpx_finder_args_t(pato::match_vector_t &matches_,
                    pato::potential_map_t &potentials_,
                    const pato::options_t &opts)
      : matches{matches_}, potentials{potentials_},
        filter_args{tpx_motifs, block_runs, encoded_seq, false, 'G', 'Y'} {
    filter_args.ornt = pato::orientation_t::both;
    min_score = opts.min_length -
                static_cast<int>(std::ceil(opts.error_rate * opts.min_length));
  }
#else
  tpx_finder_args_t(pato::match_vector_t &matches_, const pato::options_t &opts)
      : matches{matches_}, filter_args{tpx_motifs, block_runs, encoded_seq,
                                       false,      'G',        'Y'} {
    filter_args.ornt = pato::orientation_t::both;
    min_score = opts.min_length -
                static_cast<int>(std::ceil(opts.error_rate * opts.min_length));
  }
#endif
};

} // namespace

static void make_triplex_parser(tpx_finder_args_t &args,
                                unsigned max_interrupts) {
  pato::triplex_t valid_chars = "GAR";
  pato::triplex_t invalid_chars = "TCYN";
  pato::make_parser(args.tpx_parser, valid_chars, invalid_chars,
                    max_interrupts);
}

static unsigned find_tpx_motifs(pato::triplex_t &sequence, unsigned id,
                                tpx_finder_args_t &tpx_args,
                                const pato::options_t &opts) {
  pato::parse_segments(tpx_args.tpx_parser, tpx_args.segments, sequence,
                       opts.max_interruptions, opts.min_length);

  unsigned matches = 0;
  for (auto &segment : tpx_args.segments) {
    pato::motif_t motif{segment, true, id, false, '+'};
    matches += filter_guanine_error_rate(motif, tpx_args.filter_args,
                                         pato::tts_t{}, opts);
  }
  tpx_args.segments.clear();

  return matches;
}

static void search_triplex(pato::motif_t &tfo_motif, unsigned tfo_id,
                           pato::motif_t &tts_motif, unsigned tts_id,
                           tpx_finder_args_t &tpx_args,
                           const pato::options_t &opts) {
  auto tfo_candidate = seqan::ttsString(tfo_motif);
  auto tts_candidate = seqan::ttsString(tts_motif);

  int tfo_length = seqan::length(tfo_candidate);
  int tts_length = seqan::length(tts_candidate);

  for (int diag = -(tts_length - opts.min_length);
       diag <= tfo_length - opts.min_length; ++diag) {
    int tfo_offset = 0;
    int tts_offset = 0;
    if (diag < 0) {
      tts_offset = -diag;
    } else {
      tfo_offset = diag;
    }
    int match_length =
        std::min(tts_length - tts_offset, tfo_length - tfo_offset);

    int match_score = 0;
    for (int i = 0; i < match_length; ++i) {
      if (tfo_candidate[tfo_offset + i] == tts_candidate[tts_offset + i]) {
        ++match_score;
      }
    }
    if (match_score < tpx_args.min_score) {
      continue;
    }

    pato::triplex_t tmp_tts{
        seqan::infix(tts_candidate, tts_offset, tts_offset + match_length)};

    for (int i = 0; i < match_length; ++i) {
      if (tmp_tts[i] != tfo_candidate[tfo_offset + i]) {
        tmp_tts[i] = 'N';
      }
    }

    unsigned total = find_tpx_motifs(tmp_tts, seqan::getSequenceNo(tts_motif),
                                     tpx_args, opts);
    if (total == 0) {
      tpx_args.tpx_motifs.clear();
      continue;
    }

    char strand;
    std::size_t tfo_start, tfo_end;
    std::size_t tts_start, tts_end;
    for (auto &triplex : tpx_args.tpx_motifs) {
      unsigned score = 0;
      unsigned guanines = 0;

      for (unsigned i = seqan::beginPosition(triplex);
           i < seqan::endPosition(triplex); ++i) {
        if (tmp_tts[i] != 'N') {
          ++score;
        }
        if (tmp_tts[i] == 'G') {
          ++guanines;
        }
      }

      if (seqan::isParallel(tfo_motif)) {
        tfo_start = tfo_offset + seqan::beginPosition(tfo_motif) +
                    seqan::beginPosition(triplex);
        tfo_end = tfo_start + seqan::length(triplex);
      } else {
        tfo_end = seqan::endPosition(tfo_motif) -
                  (tfo_offset + seqan::beginPosition(triplex));
        tfo_start = tfo_end - seqan::length(triplex);
      }

      if (seqan::getMotif(tts_motif) == '+') {
        tts_start = tts_offset + seqan::beginPosition(tts_motif) +
                    seqan::beginPosition(triplex);
        tts_end = tts_start + seqan::length(triplex);
        strand = '+';
      } else {
        tts_end = seqan::endPosition(tts_motif) -
                  (tts_offset + seqan::beginPosition(triplex));
        tts_start = tts_end - seqan::length(triplex);
        strand = '-';
      }

      // FIXME: Use list initialization and make all type narrowing explicit!
      pato::match_t match(tfo_id, tfo_start, tfo_end,
                          seqan::getSequenceNo(tts_motif), tts_id, tts_start,
                          tts_end, score, seqan::isParallel(tfo_motif),
                          seqan::getMotif(tfo_motif), strand, guanines);
      tpx_args.matches.push_back(match);
    }
    tpx_args.tpx_motifs.clear();

    auto key = std::make_pair(seqan::getSequenceNo(tfo_motif),
                              seqan::getSequenceNo(tts_motif));
    auto result_ptr = tpx_args.potentials.find(key);
    if (result_ptr != tpx_args.potentials.end()) {
      seqan::addCount(result_ptr->second, total, seqan::getMotif(tfo_motif));
    } else {
      pato::potential_t potential(key);
      seqan::addCount(potential, total, seqan::getMotif(tfo_motif));
      seqan::setNorm(potential, seqan::length(seqan::host(tfo_motif)),
                     seqan::length(seqan::host(tts_motif)), opts.max_length,
                     opts.min_length);
      tpx_args.potentials.insert(
          std::make_pair(std::move(key), std::move(potential)));
    }
  }
}

static void match_tfo_tts_motifs(
#if !defined(_OPENMP)
    pato::match_vector_t &matches, pato::potential_map_t &potentials,
    pato::motif_vector_t &tfo_motifs, pato::motif_vector_t &tts_motifs,
    const pato::options_t &opts
#else
    pato::match_vector_vector_t &matches, pato::potential_map_t &potentials,
    pato::motif_vector_t &tfo_motifs, pato::motif_vector_t &tts_motifs,
    const pato::options_t &opts
#endif
) {
#if defined(_OPENMP)
  matches.resize(omp_get_max_threads());
#endif

#pragma omp parallel
  {
#if !defined(_OPENMP)
    tpx_finder_args_t tpx_args{matches, potentials, opts};
#else
    tpx_finder_args_t tpx_args{matches[omp_get_thread_num()], opts};
#endif

    make_triplex_parser(tpx_args, opts.max_interruptions);

    auto tfo_size = static_cast<uint64_t>(tfo_motifs.size());
    auto tts_size = static_cast<uint64_t>(tts_motifs.size());
#if defined(_OPENMP)
    uint64_t chunk_size = std::min(tfo_size, tts_size);
#endif

#pragma omp for schedule(dynamic, chunk_size) collapse(2) nowait
    for (uint64_t i = 0; i < tfo_size; i++) {
      for (uint64_t j = 0; j < tts_size; j++) {
        search_triplex(tfo_motifs[i], i, tts_motifs[j], j, tpx_args, opts);
      }
    }

#if defined(_OPENMP)
#pragma omp critical(potential_lock)
    {
      for (auto &potential_entry : tpx_args.potentials) {
        auto result_ptr = potentials.find(potential_entry.first);
        if (result_ptr == potentials.end()) {
          potentials.insert(std::move(potential_entry));
        } else {
          seqan::mergeCount(result_ptr->second, potential_entry.second);
        }
      }
    } // #pragma omp critical
#endif
  } // #pragma omp parallel
}

pato::find_tpx_result pato::find_tpxes(const pato::options_t &opts) {
  pato::sequence_loader_t tfo_sequence_loader;
  if (!tfo_sequence_loader.init(opts.tfo_file)) {
    return pato::find_tpx_result::cannot_open_tfo_file;
  }
  pato::sequence_loader_t tts_sequence_loader;
  if (!tts_sequence_loader.init(opts.tts_file)) {
    return pato::find_tpx_result::cannot_open_tts_file;
  }

  auto output_writer = output_writer_t::create(opts);
  if (!output_writer) {
    return pato::find_tpx_result::cannot_create_output_file;
  }

  pato::name_vector_t tfo_names;
  pato::triplex_vector_t tfo_sequences;
  tfo_sequence_loader.load_sequences(tfo_sequences, tfo_names,
                                     std::numeric_limits<unsigned>::max());
  pato::motif_vector_t tfo_motifs;
  pato::motif_potential_vector_t tfo_potentials;
  pato::find_tfo_motifs(tfo_motifs, tfo_potentials, tfo_sequences, opts);

  pato::name_vector_t tts_names;
  pato::triplex_vector_t tts_sequences;
  pato::motif_vector_t tts_motifs;
  pato::motif_potential_vector_t tts_potentials;

#if !defined(_OPENMP)
  pato::match_vector_t matches;
#else
  pato::match_vector_vector_t matches;
#endif
  pato::potential_map_t potentials;

  while (true) {
    if (!tts_sequence_loader.load_sequences(tts_sequences, tts_names,
                                            opts.chunk_size)) {
      break;
    }

    pato::find_tts_motifs(tts_motifs, tts_potentials, tts_sequences, opts);
    match_tfo_tts_motifs(matches, potentials, tfo_motifs, tts_motifs, opts);

#pragma omp parallel sections num_threads(2)
    {
#pragma omp section
      output_writer->print_triplexes(matches, tfo_motifs, tfo_names, tts_motifs,
                                     tts_names);
#pragma omp section
      output_writer->print_triplex_summary(potentials, tfo_names, tts_names);
    } // #pragma omp parallel sections num_threads(2)

    tts_names.clear();
    tts_motifs.clear();
    tts_sequences.clear();

#if !defined(_OPENMP)
    matches.clear();
#else
    for (auto &local_matches : matches) {
      local_matches.clear();
    }
#endif
    potentials.clear();
  }

  return pato::find_tpx_result::success;
}
