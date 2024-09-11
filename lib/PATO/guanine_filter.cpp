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

#include "guanine_filter.h"

#include <seqan/misc/interval_tree.h>

static void increase_right(pato::char_vector_vector_t &encoded_seq,
                           unsigned &pos, unsigned &filter_chars,
                           unsigned &interrupt_chars,
                           unsigned &non_filter_chars) {
  filter_chars += encoded_seq[0][pos];
  interrupt_chars += encoded_seq[1][pos];
  non_filter_chars += encoded_seq[2][pos];

  ++pos;
}

static void increase_left(pato::char_vector_vector_t &encoded_seq,
                          unsigned &pos, unsigned &filter_chars,
                          unsigned &interrupt_chars,
                          unsigned &non_filter_chars) {
  filter_chars -= encoded_seq[0][pos];
  interrupt_chars -= encoded_seq[1][pos];
  non_filter_chars -= encoded_seq[2][pos];

  ++pos;
}

static bool is_interrupt_char(pato::char_vector_vector_t &encoded_seq,
                              unsigned pos) {
  return encoded_seq[1][pos];
}

template <typename tag_t>
static bool
motif_specific_constraint([[maybe_unused]] double filter_rate,
                          [[maybe_unused]] pato::orientation_t ornt,
                          [[maybe_unused]] const tag_t &tag,
                          [[maybe_unused]] const pato::options_t &opts) {
  return true;
}

static bool
motif_specific_constraint(double filter_rate, pato::orientation_t ornt,
                          [[maybe_unused]] const pato::mixed_motif_t &tag,
                          const pato::options_t &opts) {
  if (ornt != pato::orientation_t::parallel &&
      filter_rate >= opts.mixed_antiparallel_min_guanine) {
    return true;
  } else if (ornt != pato::orientation_t::antiparallel &&
             filter_rate <= opts.mixed_parallel_max_guanine) {
    return true;
  }
  return false;
}

static void add_match(pato::motif_vector_t &motifs, pato::motif_t &motif,
                      unsigned start, unsigned end, unsigned errors,
                      [[maybe_unused]] const pato::tfo_t &tag) {
  auto motif_length = seqan::length(motif);
  auto st = seqan::isParallel(motif) ? start : motif_length - end;
  auto nd = seqan::isParallel(motif) ? end : motif_length - start;

  pato::motif_t tmp_motif{seqan::host(motif),
                          seqan::beginPosition(motif) + st,
                          seqan::beginPosition(motif) + nd,
                          seqan::isParallel(motif),
                          seqan::getSequenceNo(motif),
                          seqan::isTFO(motif),
                          seqan::getMotif(motif)};
  seqan::setScore(tmp_motif, end - start - errors);

  motifs.push_back(std::move(tmp_motif));
}

static void add_match(pato::motif_vector_t &motifs, pato::motif_t &motif,
                      unsigned start, unsigned end, unsigned errors,
                      [[maybe_unused]] const pato::mixed_motif_t &tag) {
  add_match(motifs, motif, start, end, errors, pato::tfo_t());
}

static void add_match(pato::motif_vector_t &motifs, pato::motif_t &motif,
                      unsigned start, unsigned end, unsigned errors,
                      [[maybe_unused]] const pato::purine_motif_t &tag) {
  add_match(motifs, motif, start, end, errors, pato::tfo_t());
}

static void add_match(pato::motif_vector_t &motifs, pato::motif_t &motif,
                      unsigned start, unsigned end, unsigned errors,
                      [[maybe_unused]] const pato::pyrimidine_motif_t &tag) {
  add_match(motifs, motif, start, end, errors, pato::tfo_t());
}

static void add_match(pato::motif_vector_t &motifs, pato::motif_t &motif,
                      unsigned start, unsigned end, unsigned errors,
                      [[maybe_unused]] const pato::tts_t &tag) {
  auto motif_length = seqan::length(motif);
  auto st = seqan::getMotif(motif) == '+' ? start : motif_length - end;
  auto nd = seqan::getMotif(motif) == '+' ? end : motif_length - start;

  pato::motif_t tmp_motif{seqan::host(motif),
                          seqan::beginPosition(motif) + st,
                          seqan::beginPosition(motif) + nd,
                          seqan::isParallel(motif),
                          seqan::getSequenceNo(motif),
                          seqan::isTFO(motif),
                          seqan::getMotif(motif)};
  seqan::setScore(tmp_motif, end - start - errors);

  motifs.push_back(std::move(tmp_motif));
}

template <typename string_t>
static void encode_sequence(string_t &motif, char filter_char,
                            char interrupt_char,
                            pato::char_vector_vector_t &block_runs,
                            pato::char_vector_vector_t &encoded_seq,
                            unsigned min_block_run) {
  if (encoded_seq.empty() || seqan::length(motif) > encoded_seq[0].size()) {
    encoded_seq.clear();
    encoded_seq.resize(3, pato::char_vector_t(seqan::length(motif) * 2, false));
  } else {
    encoded_seq[0].assign(seqan::length(motif), false);
    encoded_seq[1].assign(seqan::length(motif), false);
    encoded_seq[2].assign(seqan::length(motif), false);
  }

  unsigned counter = 0;
  unsigned run_counter = 0;

  for (auto m : motif) {
    if (m == filter_char) {
      encoded_seq[0][counter] = true;
    } else if (m == interrupt_char) {
      encoded_seq[1][counter] = true;
      encoded_seq[2][counter] = true;

      if (counter - run_counter >= min_block_run) {
        for (unsigned i = 0; i <= counter - min_block_run; ++i) {
          for (unsigned j = std::max(i, run_counter) + min_block_run;
               j <= seqan::length(motif); ++j) {
            block_runs[i][j] = true;
          }
        }
      }

      run_counter = counter + 1;
    } else {
      encoded_seq[2][counter] = true;
    }
    ++counter;
  }

  if (counter - run_counter >= min_block_run) {
    for (unsigned i = 0; i <= counter - min_block_run; ++i) {
      for (unsigned j = std::max(i, run_counter) + min_block_run;
           j <= seqan::length(motif); ++j) {
        block_runs[i][j] = true;
      }
    }
  }
}

namespace {

template <typename pair_t> struct second_t {
  typename pair_t::second_type operator()(const pair_t &p) const {
    return p.second;
  }
};

} // namespace

template <typename map_t>
static second_t<typename map_t::value_type>
second([[maybe_unused]] const map_t &m) {
  return second_t<typename map_t::value_type>{};
}

static void reduce_motif_vector(pato::motif_vector_t &output,
                                pato::motif_vector_t &input) {
  using interval_value_t = seqan::Position<pato::motif_t>::Type;
  using graph_t = seqan::Graph<seqan::Undirected<void>>;
  using vertex_descriptor_t = seqan::VertexDescriptor<graph_t>::Type;
  using interval_t =
      seqan::IntervalAndCargo<interval_value_t, vertex_descriptor_t>;
  using graph_size_t = typename seqan::Size<graph_t>::Type;
  using interval_tree_t =
      seqan::IntervalTree<interval_value_t, vertex_descriptor_t>;
  using vertex_iterator_t =
      seqan::Iterator<graph_t, seqan::VertexIterator>::Type;
  using counts_t = seqan::String<vertex_descriptor_t>;
  using component_t = seqan::String<graph_size_t>;
  using prop_map_t = std::unordered_map<vertex_descriptor_t, pato::motif_t>;
  using comp_map_t = std::unordered_map<graph_size_t, pato::motif_t>;

  if (input.empty()) {
    return;
  } else if (input.size() == 1) {
    output.push_back(std::move(input[0]));
    return;
  }

  graph_t parser;
  prop_map_t prop_map;
  std::vector<interval_t> intervals;

  prop_map.reserve(input.size());
  intervals.resize(input.size());

  unsigned count = 0;
  for (auto &motif : input) {
    vertex_descriptor_t vtx = seqan::addVertex(parser);

    intervals[count].i1 =
        static_cast<interval_value_t>(seqan::beginPosition(motif));
    intervals[count].i2 =
        static_cast<interval_value_t>(seqan::endPosition(motif));
    intervals[count].cargo = vtx;

    prop_map.insert(std::make_pair(vtx, std::move(motif)));

    ++count;
  }

  counts_t tree_results;
  interval_tree_t tree(intervals);
  vertex_iterator_t vertex_it(parser);
  while (!seqan::atEnd(vertex_it)) {
    auto &motif = prop_map.find(*vertex_it)->second;

    seqan::findIntervals(tree_results, tree, seqan::beginPosition(motif),
                         seqan::endPosition(motif));
    for (auto &result : tree_results) {
      if (*vertex_it == result) {
        continue;
      }
      seqan::addEdge(parser, *vertex_it, result);
    }

    seqan::clear(tree_results);
    ++vertex_it;
  }

  component_t components;
  graph_size_t num_components = seqan::connectedComponents(components, parser);

  comp_map_t comp_map;
  comp_map.reserve(num_components);

  seqan::goBegin(vertex_it);
  while (!seqan::atEnd(vertex_it)) {
    auto result_ptr = comp_map.find(seqan::getProperty(components, *vertex_it));

    if (result_ptr != comp_map.end()) {
      auto &motif = result_ptr->second;
      merge(motif, prop_map.find(*vertex_it)->second);
    } else {
      auto &motif = prop_map.find(*vertex_it)->second;
      comp_map.insert(std::make_pair(seqan::getProperty(components, *vertex_it),
                                     std::move(motif)));
    }

    ++vertex_it;
  }

  output.reserve(output.size() + num_components);
  std::transform(comp_map.begin(), comp_map.end(), std::back_inserter(output),
                 second(comp_map));
}

template <typename tag_t>
unsigned pato::filter_guanine_error_rate(motif_t &motif,
                                         pato::guanine_filter_args_t &args,
                                         const tag_t &tag,
                                         const pato::options_t &opts) {
  pato::motif_vector_t tmp_set;
  pato::motif_vector_t &motifs_ref = args.reduce_set ? tmp_set : args.motifs;

  // The following lines optimize memory usage and runtime performance in PATO
  // by reusing allocated memory. Since this function is called frequently
  // throughout the application, minimizing memory allocations and deallocations
  // is crucial. By reusing the capacity already allocated for storing the
  // encoded sequence, we reduce memory pressure and improve overall efficiency.
  auto motif_length = seqan::length(motif);
  if (args.block_runs.empty() ||
      motif_length - opts.min_block_run + 1 > args.block_runs.size()) {
    args.block_runs.clear();
    args.block_runs.resize((motif_length - opts.min_block_run + 1) * 2,
                           pato::char_vector_t((motif_length + 1) * 2, false));
  } else {
    for (unsigned i = 0; i < motif_length - opts.min_block_run + 1; i++) {
      args.block_runs[i].assign(motif_length + 1, false);
    }
  }

  char filter_char = args.filter_char;
  char interrupt_char = args.interrupt_char;
  if (opts.min_guanine_rate <= 0.0) {
    pato::filter_t filtered_sequence{motif};
    filter_char = filter_char == 'G' ? 'R' : 'Y';
    encode_sequence(filtered_sequence, filter_char, interrupt_char,
                    args.block_runs, args.encoded_seq, opts.min_block_run);
  } else {
    encode_sequence(motif, filter_char, interrupt_char, args.block_runs,
                    args.encoded_seq, opts.min_block_run);
  }

  double max_error = std::floor(motif_length * opts.error_rate);
  double max_tolerated =
      std::floor(motif_length * (1.0 - opts.min_guanine_rate));
  if (opts.maximal_error >= 0) {
    max_error = std::min(max_error, static_cast<double>(opts.maximal_error));
  }

  unsigned filter_chars = 0;
  unsigned interrupt_chars = 0;
  unsigned non_filter_chars = 0;

  unsigned max_length = motif_length;
  if (opts.max_length >= opts.min_length) {
    max_length = static_cast<unsigned>(opts.max_length);
  }

  unsigned tmp_start = 0;
  unsigned tmp_end = 0;
  unsigned tmp_errors = 0;
  unsigned covered_end = 0;

  unsigned left = 0;
  unsigned right = 0;

  unsigned matches = 0;
  while (args.block_runs[left][motif_length] &&
         left + opts.min_length <= motif_length) {
    while (static_cast<int>(right - left) < opts.min_length &&
           right < motif_length) {
      while (static_cast<int>(right - left) < opts.min_length &&
             right < motif_length) {
        increase_right(args.encoded_seq, right, filter_chars, interrupt_chars,
                       non_filter_chars);
      }
      while (interrupt_chars > max_error) {
        increase_left(args.encoded_seq, left, filter_chars, interrupt_chars,
                      non_filter_chars);
      }
      while (non_filter_chars > max_tolerated) {
        increase_left(args.encoded_seq, left, filter_chars, interrupt_chars,
                      non_filter_chars);
      }
      while (left < motif_length && is_interrupt_char(args.encoded_seq, left)) {
        increase_left(args.encoded_seq, left, filter_chars, interrupt_chars,
                      non_filter_chars);
      }

      if (right < left) {
        right = left;

        filter_chars = 0;
        interrupt_chars = 0;
        non_filter_chars = 0;
      }
    }

    if (static_cast<int>(right - left) < opts.min_length) {
      break;
    }

    bool is_match = false;

    while (interrupt_chars <= max_error && non_filter_chars <= max_tolerated &&
           right - left <= max_length) {
      double filter_chars_rate =
          static_cast<double>(filter_chars) / (right - left);
      double interrupt_chars_rate =
          static_cast<double>(interrupt_chars) / (right - left);

      if (args.block_runs[left][right] &&
          !is_interrupt_char(args.encoded_seq, right - 1) &&
          interrupt_chars_rate <= opts.error_rate &&
          opts.min_guanine_rate <= filter_chars_rate &&
          filter_chars_rate <= opts.max_guanine_rate &&
          motif_specific_constraint(filter_chars_rate, args.ornt, tag, opts)) {
        is_match = true;
        ++matches;

        tmp_start = left;
        tmp_end = right;
        tmp_errors = interrupt_chars;

        if (opts.all_matches) {
          covered_end = tmp_end;
          is_match = false;
          add_match(motifs_ref, motif, tmp_start, tmp_end, tmp_errors, tag);
        }
      }

      if (right < motif_length) {
        increase_right(args.encoded_seq, right, filter_chars, interrupt_chars,
                       non_filter_chars);
      } else {
        break;
      }
    }

    if (is_match && tmp_end > covered_end) {
      covered_end = tmp_end;
      add_match(motifs_ref, motif, tmp_start, tmp_end, tmp_errors, tag);
    }

    ++left;
    while (left < motif_length && is_interrupt_char(args.encoded_seq, left)) {
      ++left;
    }

    right = left;

    filter_chars = 0;
    interrupt_chars = 0;
    non_filter_chars = 0;
  }

  if (args.reduce_set) {
    reduce_motif_vector(args.motifs, motifs_ref);
  }

  return matches;
}

template unsigned pato::filter_guanine_error_rate<pato::tts_t>(
    motif_t &motif, pato::guanine_filter_args_t &args, const pato::tts_t &tag,
    const pato::options_t &opts);
template unsigned pato::filter_guanine_error_rate<pato::mixed_motif_t>(
    motif_t &motif, pato::guanine_filter_args_t &args,
    const pato::mixed_motif_t &tag, const pato::options_t &opts);
template unsigned pato::filter_guanine_error_rate<pato::purine_motif_t>(
    motif_t &motif, pato::guanine_filter_args_t &args,
    const pato::purine_motif_t &tag, const pato::options_t &opts);
template unsigned pato::filter_guanine_error_rate<pato::pyrimidine_motif_t>(
    motif_t &motif, pato::guanine_filter_args_t &args,
    const pato::pyrimidine_motif_t &tag, const pato::options_t &opts);
