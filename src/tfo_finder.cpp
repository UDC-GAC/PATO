#include "tfo_finder.hpp"

#include <numeric>
#include <iostream>
#include <algorithm>

#include "repeat_filter.hpp"
#include "guanine_filter.hpp"
#include "segment_parser.hpp"
#include "sequence_loader.hpp"

struct tfo_arguments
{
    graph_t tc_parser;
    graph_t ga_parser;
    graph_t gt_parser;

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
    tfo_arguments(motif_set_t& _motifs)
        : motifs(_motifs), filter_args(motifs, block_runs, encoded_seq)
    {
        filter_args.filter_char = 'G';
        filter_args.interrupt_char = 'N';
        filter_args.reduce_set = true;
    }
#else
    tfo_arguments() : filter_args(motifs, block_runs, encoded_seq)
    {
        filter_args.filter_char = 'G';
        filter_args.interrupt_char = 'N';
        filter_args.reduce_set = true;
    }
#endif
};

void make_tfo_parsers(tfo_arguments& args, unsigned int max_interrupts)
{
    triplex_t valid_chars, invalid_chars;

    valid_chars = "TCY";
    invalid_chars = "GARN";
    make_parser(args.tc_parser, valid_chars, invalid_chars, max_interrupts);

    valid_chars = "GAR";
    invalid_chars = "TCYN";
    make_parser(args.ga_parser, valid_chars, invalid_chars, max_interrupts);

    valid_chars = "GTK";
    invalid_chars = "CAMN";
    make_parser(args.gt_parser, valid_chars, invalid_chars, max_interrupts);
}

void find_tfo_motifs(triplex_t& sequence,
                     unsigned int id,
                     tfo_arguments& tfo_args,
                     const options& opts)
{
    // TODO: use motifs opts
    // TC motif
    parse_segments(tfo_args.tc_parser,
                   tfo_args.segments,
                   sequence,
                   opts.max_interruptions,
                   opts.min_length);

    tfo_args.filter_args.ornt = orientation::parallel;

    for (auto& segment : tfo_args.segments) {
        motif_t motif(segment, true, id, true, 'Y');
        filter_guanine_error_rate(motif,
                                  tfo_args.filter_args,
                                  pyrimidine_motif_t(),
                                  opts);
    }
    tfo_args.segments.clear();

    // GA motif
    parse_segments(tfo_args.ga_parser,
                   tfo_args.segments,
                   sequence,
                   opts.max_interruptions,
                   opts.min_length);

    tfo_args.filter_args.ornt = orientation::antiparallel;

    for (auto& segment : tfo_args.segments) {
        motif_t motif(segment, false, id, true, 'R');
        filter_guanine_error_rate(motif,
                                  tfo_args.filter_args,
                                  purine_motif_t(),
                                  opts);
    }
    tfo_args.segments.clear();

    // GT motif
    parse_segments(tfo_args.gt_parser,
                   tfo_args.segments,
                   sequence,
                   opts.max_interruptions,
                   opts.min_length);

    tfo_args.filter_args.ornt = orientation::both;

    for (auto& segment : tfo_args.segments) {
        if ((tfo_args.filter_args.ornt == orientation::both
             || tfo_args.filter_args.ornt == orientation::parallel)
            && opts.mixed_parallel_max_guanine >= opts.min_guanine_rate) {
            motif_t motif(segment, true, id, true, 'M');
            filter_guanine_error_rate(motif,
                                      tfo_args.filter_args,
                                      mixed_motif_t(),
                                      opts);
        }
        if ((tfo_args.filter_args.ornt == orientation::both
             || tfo_args.filter_args.ornt == orientation::antiparallel)
            && opts.mixed_antiparallel_min_guanine <= opts.max_guanine_rate) {
            motif_t motif(segment, false, id, true, 'M');
            filter_guanine_error_rate(motif,
                                      tfo_args.filter_args,
                                      mixed_motif_t(),
                                      opts);
        }
    }
    tfo_args.segments.clear();
}

bool find_tfo_motifs(motif_set_t& motifs,
                     triplex_set_t& sequences,
                     name_set_t& names,
                     const options& opts)
{
    if (!load_sequences(sequences, names, seqan::toCString(opts.tfo_file))) {
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
    tfo_arguments tfo_args(motifs);
#else
    tfo_arguments tfo_args;
#endif
    make_tfo_parsers(tfo_args, opts.max_interruptions);

#pragma omp for schedule(dynamic) nowait
    for (unsigned int i = 0; i < sequences.size(); i++) {
        find_tfo_motifs(sequences[indices[i]], indices[i], tfo_args, opts);
    }

#if defined(_OPENMP)
#pragma omp critical
{
    motifs.reserve(motifs.size() + tfo_args.motifs.size());
    std::move(tfo_args.motifs.begin(), tfo_args.motifs.end(), std::back_inserter(motifs));
} // #pragma omp critical
#endif
} // #pragma omp parallel
    double nd = omp_get_wtime();
    std::cout << "TFO in: " << nd - st << " seconds (" << motifs.size() << ")\n";

    return true;
}
