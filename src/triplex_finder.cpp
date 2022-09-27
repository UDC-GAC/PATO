#include "triplex_finder.hpp"

#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>

#include <seqan/sequence.h>
#include <seqan/parallel/parallel_macros.h>

#include "tfo_finder.hpp"
#include "tts_finder.hpp"
#include "output_writer.hpp"
#include "segment_parser.hpp"
#include "guanine_filter.hpp"
#include "triplex_pattern.hpp"

struct tpx_arguments
{
    graph_t tpx_parser;

    motif_set_t tpx_motifs;

    char_set_set_t block_runs;
    char_set_set_t encoded_seq;

    segment_set_t segments;

    match_set_t& matches;

#if !defined(_OPENMP)
    potential_set_t& potentials;
#else
    potential_set_t potentials;
#endif

    filter_arguments filter_args;

    int min_score;

#if !defined(_OPENMP)
    tpx_arguments(match_set_t& _matches,
                  potential_set_t& _potentials,
                  const options& opts)
        : matches(_matches),
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
    tpx_arguments(match_set_t& _matches,
                  const options& opts)
        : matches(_matches), filter_args(tpx_motifs, block_runs, encoded_seq)
    {
        filter_args.ornt = orientation::both;
        filter_args.filter_char = 'G';
        filter_args.interrupt_char = 'Y';
        filter_args.reduce_set = false;

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

void search_triplex(motif_t& tfo_motif,
                    unsigned int tfo_id,
                    motif_t& tts_motif,
                    unsigned int tts_id,
                    tpx_arguments& tpx_args,
                    const options& opts)
{
    auto tfo_candidate = seqan::ttsString(tfo_motif);
    auto tts_candidate = seqan::ttsString(tts_motif);

    int tfo_length = seqan::length(tfo_candidate);
    int tts_length = seqan::length(tts_candidate);

    for (int diag = -(tts_length - opts.min_length); diag <= tfo_length - opts.min_length; diag++) {
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
            if (tfo_candidate[tfo_offset + i] == tts_candidate[tts_offset + i]) {
                score++;
            }
        }
        if (score < tpx_args.min_score) {
            continue;
        }

        triplex_t tmp_tts(seqan::infix(tts_candidate,
                                       tts_offset,
                                       tts_offset + match_length));

        for (int i = 0; i < match_length; i++) {
            if (tmp_tts[i] != tfo_candidate[tfo_offset + i]) {
                tmp_tts[i] = 'N';
            }
        }

        unsigned int total = find_tpx_motifs(tmp_tts,
                                             seqan::getSequenceNo(tts_motif),
                                             tpx_args,
                                             opts);
        if (total == 0) {
            tpx_args.tpx_motifs.clear();
            continue;
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
                tfo_start = tfo_offset
                            + seqan::beginPosition(tfo_motif)
                            + seqan::beginPosition(triplex);
                tfo_end = tfo_start + seqan::length(triplex);
            } else {
                tfo_end = seqan::endPosition(tfo_motif)
                          - (tfo_offset + seqan::beginPosition(triplex));
                tfo_start = tfo_end - seqan::length(triplex);
            }

            if (seqan::getMotif(tts_motif) == '+') {
                tts_start = tts_offset
                            + seqan::beginPosition(tts_motif)
                            + seqan::beginPosition(triplex);
                tts_end = tts_start + seqan::length(triplex);
                strand = '+';
            } else {
                tts_end = seqan::endPosition(tts_motif)
                          - (tts_offset + seqan::beginPosition(triplex));
                tts_start = tts_end - seqan::length(triplex);
                strand = '-';
            }

            match_t match(tfo_id,
                          tfo_start,
                          tfo_end,
                          seqan::getSequenceNo(tts_motif),
                          tts_id,
                          tts_start,
                          tts_end,
                          score,
                          seqan::isParallel(tfo_motif),
                          seqan::getMotif(tfo_motif),
                          strand,
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
#if defined(_OPENMP)
    matches.resize(omp_get_max_threads());
#endif

    double st = omp_get_wtime();
#pragma omp parallel
{
#if !defined(_OPENMP)
    tpx_arguments tpx_args(matches, potentials, opts);
#else
    tpx_arguments tpx_args(matches[omp_get_thread_num()], opts);
#endif
    make_triplex_parser(tpx_args, opts.max_interruptions);

    uint64_t tfo_size = static_cast<uint64_t>(tfo_motifs.size());
    uint64_t tts_size = static_cast<uint64_t>(tts_motifs.size());
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

    print_triplex_pairs(matches, tfo_motifs, tfo_names, tts_motifs, tts_names, opts);
    print_triplex_summary(potentials, tfo_names, tts_names, opts);
}
