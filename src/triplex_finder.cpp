#include "triplex_finder.hpp"

#if defined(_OPENMP)
#define CACHE_LINE_SIZE 64

#include <stdlib.h>
#include <string.h>
#endif

#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <utility>
#include <iostream>
#include <algorithm>

#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/parallel/parallel_macros.h>

#include "tfo_finder.hpp"
#include "tts_finder.hpp"
#include "segment_parser.hpp"
#include "guanine_filter.hpp"
#include "triplex_pattern.hpp"

struct triplex_pair
{
    unsigned int tfo_id;
    unsigned int tts_id;

    int tfo_offset;
    int tts_offset;
    int length;

    triplex_pair(unsigned int _tfo_id,
                 unsigned int _tts_id,
                 int _tfo_offset,
                 int _tts_offset,
                 int _length)
        : tfo_id(_tfo_id),
          tts_id(_tts_id),
          tfo_offset(_tfo_offset),
          tts_offset(_tts_offset),
          length(_length)
    {}
};

typedef std::vector<triplex_pair> pair_set_t;
#if defined(_OPENMP)
typedef std::vector<pair_set_t> pair_set_set_t;
#endif

struct tpx_arguments
{
    graph_t tpx_parser;

    motif_set_t tpx_motifs;

    char_set_set_t block_runs;
    char_set_set_t encoded_seq;

    segment_set_t segments;

    pair_set_t& pairs;

    match_set_t& matches;

#if !defined(_OPENMP)
    potential_set_t& potentials;
#else
    potential_set_t potentials;
#endif

    filter_arguments filter_args;

    int min_score;
#if defined(_OPENMP)
    int thread_id;
    int num_threads;
#endif

#if !defined(_OPENMP)
    tpx_arguments(pair_set_t& _pairs,
                  match_set_t& _matches,
                  potential_set_t& _potentials,
                  const options& opts)
        : pairs(_pairs),
          matches(_matches),
          potentials(_potentials),
          filter_args(tpx_motifs, block_runs, encoded_seq)
    {
        filter_args.ornt = orientation::both;
        filter_args.filter_char = 'G';
        filter_args.interrupt_char = 'Y';
        filter_args.reduce_set = false;

        min_score = opts.min_length
                    - static_cast<int>(std::ceil(opts.error_rate
                                                 * opts.min_length));
    }
#else
    tpx_arguments(pair_set_t& _pairs,
                  match_set_t& _matches,
                  const options& opts)
        : pairs(_pairs),
          matches(_matches),
          filter_args(tpx_motifs, block_runs, encoded_seq)
    {
        filter_args.ornt = orientation::both;
        filter_args.filter_char = 'G';
        filter_args.interrupt_char = 'Y';
        filter_args.reduce_set = false;

        thread_id = omp_get_thread_num();
        num_threads = omp_get_max_threads();
        min_score = opts.min_length
                    - static_cast<int>(std::ceil(opts.error_rate
                                                 * opts.min_length));
    }
#endif
};

void make_triplex_parser(tpx_arguments& args, unsigned int max_interrupts)
{
    triplex_t valid_chars = "GAR";
    triplex_t invalid_chars = "TCYN";
    make_parser(args.tpx_parser, valid_chars, invalid_chars, max_interrupts);
}

unsigned int find_tpx_motifs(triplex_t& sequence,
                             unsigned int id,
                             tpx_arguments& tpx_args,
                             const options& opts)
{
    parse_segments(tpx_args.tpx_parser,
                   tpx_args.segments,
                   sequence,
                   opts.max_interruptions,
                   opts.min_length);

    unsigned int matches = 0;
    for (auto& segment : tpx_args.segments) {
        motif_t motif(segment, true, id, false, '+');
        matches += filter_guanine_error_rate(motif,
                                             tpx_args.filter_args,
                                             tts_t(),
                                             opts);
    }
    tpx_args.segments.clear();

    return matches;
}

void search_triplex_pairs(motif_t& tfo_motif,
                          unsigned int tfo_id,
                          motif_t& tts_motif,
                          unsigned int tts_id,
                          tpx_arguments& tpx_args,
                          int min_length)
{
    auto tfo_candidate = seqan::ttsString(tfo_motif);
    auto tts_candidate = seqan::ttsString(tts_motif);

    int tfo_length = seqan::length(tfo_candidate);
    int tts_length = seqan::length(tts_candidate);

    for (int diag = -(tts_length - min_length); diag <= tfo_length - min_length; diag++) {
        int tfo_offset = 0;
        int tts_offset = 0;
        if (diag < 0) {
            tts_offset = -diag;
        } else {
            tfo_offset = diag;
        }
        int match_length = std::min(tts_length - tts_offset,
                                    tfo_length - tfo_offset);

        int score = 0;
        for (int i = 0; i < match_length; i++) {
            if (tfo_candidate[tfo_offset++] == tts_candidate[tts_offset++]) {
                score++;
            }
        }
        if (score < tpx_args.min_score) {
            continue;
        }

        triplex_pair pair(tfo_id,
                          tts_id,
                          tfo_offset - match_length,
                          tts_offset - match_length,
                          match_length);
        tpx_args.pairs.push_back(std::move(pair));
    }
}

void search_triplex(motif_t& tfo_motif,
                    motif_t& tts_motif,
                    triplex_pair& pair,
                    tpx_arguments& tpx_args,
                    const options& opts)
{
    auto tts_candidate = seqan::ttsString(tts_motif);
    auto tfo_candidate = seqan::ttsString(tfo_motif);

    triplex_t tmp_tts(seqan::infix(tts_candidate,
                                   pair.tts_offset,
                                   pair.tts_offset + pair.length));

    int offset = pair.tfo_offset;
    for (int i = 0; i < pair.length; i++) {
        if (tmp_tts[i] != tfo_candidate[offset++]) {
            tmp_tts[i] = 'N';
        }
    }

    unsigned int total = find_tpx_motifs(tmp_tts,
                                         seqan::getSequenceNo(tts_motif),
                                         tpx_args,
                                         opts);
    if (total == 0) {
        return;
    }

    char strand;
    std::size_t tfo_start, tfo_end;
    std::size_t tts_start, tts_end;
    for (auto& triplex : tpx_args.tpx_motifs) {
        int score = 0;
        int guanines = 0;

        for (unsigned int i = seqan::beginPosition(triplex); i < seqan::endPosition(triplex); i++) {
            if (tmp_tts[i] != 'N') {
                score++;
            }
            if (tmp_tts[i] == 'G') {
                guanines++;
            }
        }

        if (seqan::isParallel(tfo_motif)) {
            tfo_start = pair.tfo_offset
                        + seqan::beginPosition(tfo_motif)
                        + seqan::beginPosition(triplex);
            tfo_end = tfo_start + seqan::length(triplex);
        } else {
            tfo_end = seqan::endPosition(tfo_motif)
                        - (pair.tfo_offset + seqan::beginPosition(triplex));
            tfo_start = tfo_end - seqan::length(triplex);
        }

        if (seqan::getMotif(tts_motif) == '+') {
            tts_start = pair.tts_offset
                        + seqan::beginPosition(tts_motif)
                        + seqan::beginPosition(triplex);
            tts_end = tts_start + seqan::length(triplex);
            strand = '+';
        } else {
            tts_end = seqan::endPosition(tts_motif)
                        - (pair.tts_offset + seqan::beginPosition(triplex));
            tts_start = tts_end - seqan::length(triplex);
            strand = '-';
        }

        match_t match(pair.tfo_id,
                      tfo_start,
                      tfo_end,
                      seqan::getSequenceNo(tts_motif),
                      pair.tts_id,
                      tts_start,
                      tts_end,
                      score,
                      seqan::isParallel(tfo_motif),
                      seqan::getMotif(tfo_motif), strand,
                      guanines);
        tpx_args.matches.push_back(std::move(match));
    }
    tpx_args.tpx_motifs.clear();

    auto key = std::make_pair<unsigned int, unsigned int>(seqan::getSequenceNo(tfo_motif),
                                                          seqan::getSequenceNo(tts_motif));
    auto result_ptr = tpx_args.potentials.find(key);
    if (result_ptr != tpx_args.potentials.end()) {
        seqan::addCount(result_ptr->second, total, seqan::getMotif(tfo_motif));
    } else {
        potential_t potential(key);
        seqan::addCount(potential, total, seqan::getMotif(tfo_motif));
        seqan::setNorm(potential,
                       seqan::length(seqan::host(tfo_motif)),
                       seqan::length(seqan::host(tts_motif)),
                       opts);
        tpx_args.potentials.insert(std::make_pair<std::pair<unsigned int, unsigned int>,
                                                  potential_t>(std::move(key), std::move(potential)));
    }
}

#if !defined(_OPENMP)
void match_tfo_tts_motifs(match_set_t& matches,
                          potential_set_t& potentials,
                          motif_set_t& tfo_motifs,
                          motif_set_t& tts_motifs,
                          const options& opts)
#else
void match_tfo_tts_motifs(match_set_set_t& matches,
                          potential_set_t& potentials,
                          motif_set_t& tfo_motifs,
                          motif_set_t& tts_motifs,
                          const options& opts)
#endif
{
#if !defined(_OPENMP)
    pair_set_t pairs;
#else
    matches.resize(omp_get_max_threads());

    pair_set_set_t pairs(omp_get_max_threads());

    unsigned int *indices, *are_ready;
    posix_memalign((void **) &indices,
                   CACHE_LINE_SIZE,
                   omp_get_max_threads() * sizeof(unsigned int));
    posix_memalign((void **) &are_ready,
                   CACHE_LINE_SIZE,
                   omp_get_max_threads() * sizeof(unsigned int));
    memset(indices, 0, omp_get_max_threads() * sizeof(unsigned int));
    memset(are_ready, 0, omp_get_max_threads() * sizeof(unsigned int));
#endif

    double st = omp_get_wtime();
#pragma omp parallel
{
#if !defined(_OPENMP)
    tpx_arguments tpx_args(pairs, matches, potentials, opts);
#else
    tpx_arguments tpx_args(pairs[omp_get_thread_num()],
                           matches[omp_get_thread_num()],
                           opts);
#endif
    make_triplex_parser(tpx_args, opts.max_interruptions);

    uint64_t tfo_size = static_cast<uint64_t>(tfo_motifs.size());
    uint64_t tts_size = static_cast<uint64_t>(tts_motifs.size());
#pragma omp for schedule(static) collapse(2) nowait
    for (uint64_t i = 0; i < tfo_size; i++) {
        for (uint64_t j = 0; j < tts_size; j++) {
            search_triplex_pairs(tfo_motifs[i], i, tts_motifs[j], j, tpx_args, opts.min_length);
        }
    }

#if !defined(_OPENMP)
    for (unsigned int i = 0; i < pairs.size(); i++) {
        auto& pair = pairs[i];
        search_triplex(tfo_motifs[pair.tfo_id], tts_motifs[pair.tts_id], pair, tpx_args, opts);
    }
#else
#pragma omp atomic write
    are_ready[tpx_args.thread_id] = 1;

    int is_ready;
    int are_finished = 0;
    int steal_id = tpx_args.thread_id;

    unsigned int idx;
    while (true) {
#pragma omp atomic read
        is_ready = are_ready[steal_id];

        if (is_ready == 0) {
            are_finished = 0;
        } else if (is_ready == 2) {
            if (++are_finished == tpx_args.num_threads) {
                break;
            }
        } else {
            while (true) {
#pragma omp atomic capture
                idx = indices[steal_id]++;
                if (idx >= pairs[steal_id].size()) {
#pragma omp atomic write
                    are_ready[steal_id] = 2;
                    break;
                }

                auto& pair = pairs[steal_id][idx];
                search_triplex(tfo_motifs[pair.tfo_id], tts_motifs[pair.tts_id], pair, tpx_args, opts);
            }
        }

        steal_id = (steal_id + 1) % tpx_args.num_threads;
    }

// TODO: evaluate if using a lock over a shared data structure is faster than this copy
#pragma omp critical (potential_lock)
{
    for (auto& potential_entry : tpx_args.potentials) {
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

#if defined(_OPENMP)
    free(indices);
    free(are_ready);
#endif

    double nd = omp_get_wtime();
#if !defined(_OPENMP)
    std::cout << "TPX in: " << nd - st << " seconds (" << matches.size() << ")\n";
#else
    std::size_t total = 0;
    for (auto& local_matches : matches) {
        total += local_matches.size();
    }
    std::cout << "TPX in: " << nd - st << " seconds (" << total << ")\n";
#endif
}

seqan::CharString error_string(match_t& match,
                               motif_set_t& tfo_motifs,
                               motif_set_t& tts_motifs,
                               const options& opts)
{
    std::ostringstream errors;

    motif_t tfo_motif(seqan::host(tfo_motifs[match.tfoNo]),
                      match.oBegin,
                      match.oEnd,
                      match.parallel,
                      tfo_motifs[match.tfoNo].seqNo,
                      true,
                      match.motif);
    motif_t tts_motif(seqan::host(tts_motifs[match.ttsNo]),
                      match.dBegin,
                      match.dEnd,
                      match.parallel,
                      match.ttsSeqNo,
                      false,
                      match.strand);

    seqan::CharString tfo_ps = seqan::prettyString(tfo_motif);
    seqan::CharString tts_ps = seqan::prettyString(tts_motif);

    if (match.strand == '-') {
        seqan::reverse(tts_ps);

        if (tfo_motifs[match.tfoNo].parallel) {
            seqan::reverse(tfo_ps);
        }

        auto tts_it = seqan::begin(tts_motif);
        auto tts_end = seqan::end(tts_motif);
        auto tfo_it = seqan::begin(tfo_motif);
        auto tfo_end = seqan::end(tfo_motif);

        unsigned int i = 0;
        while (tts_end != tts_it && tfo_end != tfo_it) {
            tts_end--;
            tfo_end--;

            if (*tts_end != *tfo_end) {
                if (!isupper(seqan::value(tts_ps, i))
                    && !isupper(seqan::value(tfo_ps, i))) {
                    errors << "b" << i;
                } else if (!isupper(seqan::value(tts_ps, i))) {
                    errors << "d" << i;
                } else if (!isupper(seqan::value(tfo_ps, i))) {
                    errors << "o" << i;
                } else {
                    errors << "t" << i;
                }
            }
            i++;
        }
    } else {
        if (!tfo_motifs[match.tfoNo].parallel) {
            seqan::reverse(tfo_ps);
        }

        auto tts_it = seqan::begin(tts_motif);
        auto tts_end = seqan::end(tts_motif);
        auto tfo_it = seqan::begin(tfo_motif);
        auto tfo_end = seqan::end(tfo_motif);

        unsigned int i = 0;
        while (tts_end != tts_it && tfo_end != tfo_it) {
            if (*tts_it != *tfo_it) {
                if (!isupper(seqan::value(tts_ps, i))
                    && !isupper(seqan::value(tfo_ps, i))) {
                    errors << "b" << i;
                } else if (!isupper(seqan::value(tts_ps, i))) {
                    errors << "d" << i;
                } else if (!isupper(seqan::value(tfo_ps, i))) {
                    errors << "o" << i;
                } else {
                    errors << "t" << i;
                }
            }
            tts_it++;
            tfo_it++;
            i++;
        }
    }

    return errors.str();
}

#if !defined(_OPENMP)
void print_tfo_tts_pairs(match_set_t& matches,
                         motif_set_t& tfo_motifs,
                         name_set_t& tfo_names,
                         motif_set_t& tts_motifs,
                         name_set_t& tts_names,
                         const options& opts)
#else
void print_tfo_tts_pairs(match_set_set_t& matches,
                         motif_set_t& tfo_motifs,
                         name_set_t& tfo_names,
                         motif_set_t& tts_motifs,
                         name_set_t& tts_names,
                         const options& opts)
#endif
{
    seqan::CharString output_file_name;
    seqan::append(output_file_name, opts.output_file);
    seqan::append(output_file_name, ".out");

    std::ofstream output_file(seqan::toCString(output_file_name),
                              std::ios_base::out);
    if (!output_file) {
        std::cerr << "PATO: error opening output file '"
                  << seqan::toCString(output_file_name) << "'\n";
        return;
    }

    output_file << "# Sequence-ID\tTFO start\tTFO end\tDuplex-ID\tTTS start\tTT"
                   "S end\tScore\tError-rate\tErrors\tMotif\tStrand\tOrientatio"
                   "n\tGuanine-rate\n";

#if !defined(_OPENMP)
    for (auto& match: matches) {
#else
    for (auto& local_matches : matches) {
        for (auto& match : local_matches) {
#endif
            auto tfo_seq_id = tfo_motifs[match.tfoNo].seqNo;
            auto tts_seq_id = match.ttsSeqNo;

            output_file << tfo_names[tfo_seq_id] << "\t"
                        << match.oBegin << "\t"
                        << match.oEnd << "\t"
                        << tts_names[tts_seq_id] << "\t"
                        << match.dBegin << "\t"
                        << match.dEnd << "\t"
                        << match.mScore << "\t"
                        << std::setprecision(2) << 1.0 - match.mScore / (match.dEnd - match.dBegin) << "\t"
                        << error_string(match, tfo_motifs, tts_motifs, opts) << "\t"
                        << match.motif << "\t"
                        << match.strand << "\t"
                        << (match.parallel ? 'P' : 'A') << "\t"
                        << match.guanines / (match.dEnd - match.dBegin)
                        << "\n";
#if defined(_OPENMP)
        }
#endif
    }
}

void print_summary(potential_set_t& potentials,
                   name_set_t& tfo_names,
                   name_set_t& tts_names,
                   const options& opts)
{
    seqan::CharString output_file_name;
    seqan::append(output_file_name, opts.output_file);
    seqan::append(output_file_name, ".summary");

    std::ofstream output_file(seqan::toCString(output_file_name),
                              std::ios_base::out);
    if (!output_file) {
        std::cerr << "PATO: error opening output file '"
                  << seqan::toCString(output_file_name) << "'\n";
        return;
    }

    output_file << "# Duplex-ID\tSequence-ID\tTotal (abs)\tTotal (rel)\tGA (abs"
                   ")\tGA (rel)\tTC (abs)\tTC (rel)\tGT (abs)\tGT (rel)\n";
    for (auto& potential_entry : potentials) {
        auto& potential = potential_entry.second;
        if (seqan::hasCount(potential)) {
            output_file << tts_names[seqan::getKey(potential).second] << "\t"
                        << tfo_names[seqan::getKey(potential).first] << "\t"
                        << seqan::getCounts(potential) << "\t"
                        << std::setprecision(3) << seqan::getCounts(potential) / seqan::getNorm(potential) << "\t"
                        << seqan::getCount(potential, 'R') << "\t"
                        << std::setprecision(3) << seqan::getCount(potential, 'R') / seqan::getNorm(potential) << "\t"
                        << seqan::getCount(potential, 'Y') << "\t"
                        << std::setprecision(3) << seqan::getCount(potential, 'Y') / seqan::getNorm(potential) << "\t"
                        << seqan::getCount(potential, 'M') << "\t"
                        << std::setprecision(3) << seqan::getCount(potential, 'M') /seqan:: getNorm(potential) << "\t"
                        << "\n";
        }
    }
}

void find_triplexes(const options& opts)
{
    name_set_t tfo_names;
    motif_set_t tfo_motifs;
    triplex_set_t tfo_sequences;
    if (!find_tfo_motifs(tfo_motifs, tfo_sequences, tfo_names, opts)) {
        return;
    }

    name_set_t tts_names;
    motif_set_t tts_motifs;
    triplex_set_t tts_sequences;
    if (!find_tts_motifs(tts_motifs, tts_sequences, tts_names, opts)) {
        return;
    }

#if !defined(_OPENMP)
    match_set_t matches;
#else
    match_set_set_t matches;
#endif
    potential_set_t potentials;
    match_tfo_tts_motifs(matches, potentials, tfo_motifs, tts_motifs, opts);

    print_tfo_tts_pairs(matches, tfo_motifs, tfo_names, tts_motifs, tts_names, opts);
    print_summary(potentials, tfo_names, tts_names, opts);
}
