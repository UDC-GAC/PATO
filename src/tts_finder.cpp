#include "tts_finder.hpp"

#include <numeric>
#include <iostream>
#include <algorithm>

#include "repeat_filter.hpp"
#include "guanine_filter.hpp"
#include "segment_parser.hpp"
#include "sequence_loader.hpp"

struct tts_arguments
{
    graph_t plus_parser;
    graph_t minus_parser;

#if !defined(_OPENMP)
    motif_set_t& motifs;
#else
    motif_set_t motifs;
#endif

    char_set_set_t block_runs;
    char_set_set_t encoded_seq;

    segment_set_t segments;

    filter_arguments filter_args;

#if !defined(_OPENMP)
    tts_arguments(motif_set_t& _motifs)
        : motifs(_motifs), filter_args(motifs, block_runs, encoded_seq)
    {
        filter_args.ornt = orientation::both;
        filter_args.filter_char = 'G';
        filter_args.interrupt_char = 'Y';
        filter_args.reduce_set = true;
    }
#else
    tts_arguments() : filter_args(motifs, block_runs, encoded_seq)
    {
        filter_args.ornt = orientation::both;
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
    // + motif
    parse_segments(tts_args.plus_parser,
                   tts_args.segments,
                   sequence,
                   opts.max_interruptions,
                   opts.min_length);

    for (auto& segment : tts_args.segments) {
        motif_t motif(segment, true, id, false, '+');
        filter_guanine_error_rate(motif,
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
        filter_guanine_error_rate(motif,
                                  tts_args.filter_args,
                                  tts_t(),
                                  opts);
    }
    tts_args.segments.clear();
}

bool find_tts_motifs(motif_set_t& motifs,
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
    if (opts.filter_repeats) {
        repeat_set_t repeats;
#pragma omp for schedule(dynamic)
        for (unsigned int i = 0; i < sequences.size(); i++) {
            filter_repeats(repeats,
                           sequences[indices[i]],
                           opts.min_repeat_length,
                           opts.max_repeat_period);
        }
    }

#if !defined(_OPENMP)
    tts_arguments tts_args(motifs);
#else
    tts_arguments tts_args;
#endif
    make_tts_parsers(tts_args, opts.max_interruptions);

#pragma omp for schedule(dynamic) nowait
    for (unsigned int i = 0; i < sequences.size(); i++) {
        find_tts_motifs(sequences[indices[i]], indices[i], tts_args, opts);
    }

#if defined(_OPENMP)
#pragma omp critical
{
    motifs.reserve(motifs.size() + tts_args.motifs.size());
    std::move(tts_args.motifs.begin(), tts_args.motifs.end(), std::back_inserter(motifs));
} // #pragma omp critical
#endif
} // #pragma omp parallel
    double nd = omp_get_wtime();
    std::cout << "TTS in: " << nd - st << " seconds (" << motifs.size() << ")\n";

    return true;
}
