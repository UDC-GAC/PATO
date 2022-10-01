#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include <seqan/sequence.h>

#include "triplex_enums.hpp"

struct options
{
    seqan::CharString tfo_file;
    seqan::CharString tts_file;
    seqan::CharString output_file;

    double error_rate;
    double min_guanine_rate;
    double max_guanine_rate;
    double mixed_parallel_max_guanine;
    double mixed_antiparallel_min_guanine;

    int min_length;
    int max_length;
    int maximal_error;
    int duplicate_cutoff;

    unsigned int min_block_run;
    unsigned int min_repeat_length;
    unsigned int max_repeat_period;
    unsigned int max_interruptions;

    run_mode_t run_mode;
    output_format_t output_format;
    error_reference_t error_reference;
    detect_duplicates_t detect_duplicates;

    bool tc_motif;
    bool ga_motif;
    bool gt_p_motif;
    bool gt_a_motif;
    bool all_matches;
    bool pretty_output;
    bool filter_repeats;
    bool merge_features;
    bool same_sequence_duplicates;
    bool report_duplicate_locations;
};

#endif
