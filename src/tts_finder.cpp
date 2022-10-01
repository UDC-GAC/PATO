#include "tts_finder.hpp"

#include <numeric>
#include <iostream>
#include <algorithm>

#include <seqan/parallel/parallel_macros.h>

#include "triplex_enums.hpp"
#include "triplex_match.hpp"
#include "output_writer.hpp"
#include "repeat_filter.hpp"
#include "guanine_filter.hpp"
#include "segment_parser.hpp"
#include "sequence_loader.hpp"
#include "duplicate_filter.hpp"

struct tts_arguments
{
    graph_t plus_parser;
    graph_t minus_parser;

    repeat_set_t repeats;

#if !defined(_OPENMP)
    motif_set_t& motifs;
    motif_potential_set_t& potentials;
#else
    motif_set_t motifs;
    motif_potential_set_t potentials;
#endif

    char_set_set_t block_runs;
    char_set_set_t encoded_seq;

    segment_set_t segments;

    filter_arguments filter_args;

#if !defined(_OPENMP)
    explicit tts_arguments(motif_set_t& _motifs,
                           motif_potential_set_t& _potentials)
        : motifs(_motifs),
          potentials(_potentials),
          filter_args(motifs, block_runs, encoded_seq)
    {
        filter_args.ornt = orientation_t::both;
        filter_args.filter_char = 'G';
        filter_args.interrupt_char = 'Y';
        filter_args.reduce_set = true;
    }
#else
    tts_arguments() : filter_args(motifs, block_runs, encoded_seq)
    {
        filter_args.ornt = orientation_t::both;
        filter_args.filter_char = 'G';
        filter_args.interrupt_char = 'Y';
        filter_args.reduce_set = true;
    }
#endif
};

void make_tts_parsers(tts_arguments& args, unsigned int max_interrupts)
{
    triplex_t valid_chars, invalid_chars;

    valid_chars = "GAR";
    invalid_chars = "TCYN";
    make_parser(args.plus_parser, valid_chars, invalid_chars, max_interrupts);

    valid_chars = "TCY";
    invalid_chars = "GARN";
    make_parser(args.minus_parser, valid_chars, invalid_chars, max_interrupts);
}

void find_tts_motifs(triplex_t& sequence,
                     unsigned int id,
                     tts_arguments& tts_args,
                     const options& opts)
{
    unsigned int matches_plus, matches_minus;
    matches_plus = matches_minus = 0;

    // + motif
    parse_segments(tts_args.plus_parser,
                   tts_args.segments,
                   sequence,
                   opts.max_interruptions,
                   opts.min_length);

    for (auto& segment : tts_args.segments) {
        motif_t motif(segment, true, id, false, '+');
        matches_plus += filter_guanine_error_rate(motif,
                                                  tts_args.filter_args,
                                                  tts_t(),
                                                  opts);
    }
    tts_args.segments.clear();

    // - motif
    parse_segments(tts_args.minus_parser,
                   tts_args.segments,
                   sequence,
                   opts.max_interruptions,
                   opts.min_length);

    for (auto& segment : tts_args.segments) {
        motif_t motif(segment, true, id, false, '-');
        matches_minus += filter_guanine_error_rate(motif,
                                                   tts_args.filter_args,
                                                   tts_t(),
                                                   opts);
    }
    tts_args.segments.clear();

    if (opts.run_mode == run_mode_t::tts_search) {
        motif_potential_t potential(id);

        seqan::addCount(potential, matches_plus, '+');
        seqan::addCount(potential, matches_minus, '-');

        seqan::setNorm(potential, seqan::length(sequence), opts);

        tts_args.potentials.push_back(potential);
    }
}

bool find_tts_motifs(motif_set_t& motifs,
                     motif_potential_set_t& potentials,
                     triplex_set_t& sequences,
                     name_set_t& names,
                     const options& opts)
{
    if (!load_sequences(sequences, names, seqan::toCString(opts.tts_file))) {
        return false;
    }

    index_set_t indices(sequences.size(), 0);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](unsigned int i, unsigned int j) -> bool {
        return seqan::length(sequences[i]) > seqan::length(sequences[j]);
    });

    double st = omp_get_wtime();
#pragma omp parallel
{
#if !defined(_OPENMP)
    tts_arguments tts_args(motifs, potentials);
#else
    tts_arguments tts_args;
#endif
    if (opts.run_mode == run_mode_t::tts_search) {
        tts_args.filter_args.reduce_set = opts.merge_features;
    }

    make_tts_parsers(tts_args, opts.max_interruptions);

#pragma omp for schedule(dynamic) nowait
    for (unsigned int i = 0; i < sequences.size(); i++) {
        if (opts.filter_repeats) {
            filter_repeats(tts_args.repeats,
                           sequences[indices[i]],
                           opts.min_repeat_length,
                           opts.max_repeat_period);
        }
        find_tts_motifs(sequences[indices[i]], indices[i], tts_args, opts);
    }

#if defined(_OPENMP)
#pragma omp critical
{
    motifs.reserve(motifs.size() + tts_args.motifs.size());
    std::move(tts_args.motifs.begin(), tts_args.motifs.end(), std::back_inserter(motifs));
} // #pragma omp critical
    if (opts.run_mode == run_mode_t::tts_search) {
#pragma omp critical (potentials)
{
        potentials.reserve(potentials.size() + tts_args.potentials.size());
        std::move(tts_args.potentials.begin(), tts_args.potentials.end(), std::back_inserter(potentials));
} // #pragma omp critical
    }
#endif
} // #pragma omp parallel

    if (opts.run_mode == run_mode_t::tts_search
        && opts.detect_duplicates != detect_duplicates_t::off) {
        count_duplicates(motifs, opts);
        if (opts.duplicate_cutoff >= 0) {
            filter_duplicates(motifs, opts.duplicate_cutoff);
        }
    }

    double nd = omp_get_wtime();
    std::cout << "TTS in: " << nd - st << " seconds (" << motifs.size() << ")\n";

    return true;
}

void find_tts_motifs(const options& opts)
{
    name_set_t tts_names;
    motif_set_t tts_motifs;
    triplex_set_t tts_sequences;
    motif_potential_set_t tts_potentials;
    if (!find_tts_motifs(tts_motifs, tts_potentials, tts_sequences, tts_names, opts)) {
        return;
    }

    print_motifs(tts_motifs, tts_names, opts);
    print_summary(tts_potentials, tts_names, opts);
}
