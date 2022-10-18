#include "options.hpp"

#include <iostream>

#include <seqan/sequence.h>

void print_options(const options& opts)
{
    std::cout << "\e[1mPATO: high PerformAnce TriplexatOr\e[0m -- ";
    switch (opts.run_mode) {
        case run_mode_t::tfo_search:
            std::cout << "TFO search\n\n";
            break;
        case run_mode_t::tts_search:
            std::cout << "TTS search\n\n";
            break;
        case run_mode_t::tpx_search:
            std::cout << "TPX search\n\n";
            break;
    }

    if (opts.run_mode != run_mode_t::tts_search) {
        std::cout << "  (-ss) TFO file: " << opts.tfo_file << "\n";
    }
    if (opts.run_mode != run_mode_t::tfo_search) {
        std::cout << "  (-ds) TTS file: " << opts.tts_file << "\n";
    }
    std::cout << "\n";

    std::cout << "   (-l) Minimum triplex length: " << opts.min_length << "\n";
    std::cout << "   (-L) Maximum triplex length: " << opts.max_length << "\n";
    std::cout << "   (-e) Error rate: " << opts.error_rate * 100.0 << "%\n";
    std::cout << "   (-E) Maximal error rate: " << opts.maximal_error << "\n";
    std::cout << "   (-c) Maximum consecutive errors: " << opts.max_interruptions << "\n";
    std::cout << "   (-g) Minimal guanine proportion: " << opts.min_guanine_rate * 100.0 << "%\n";
    std::cout << "   (-G) Maximal guanine proportion: " << opts.max_guanine_rate * 100.0 << "%\n";
    if (opts.run_mode != run_mode_t::tts_search) {
        std::cout << "   (-m) Triplex motifs allowed \n";
        std::cout << "        - TC: " << (opts.tc_motif ? "on" : "off") << "\n";
        std::cout << "        - GA: " << (opts.ga_motif ? "on" : "off") << "\n";
        std::cout << "        - Parallel GT: " << (opts.gt_p_motif ? "on" : "off") << "\n";
        std::cout << "        - Antiprallel GT: " << (opts.gt_a_motif ? "on" : "off") << "\n"; 
    }
    std::cout << "(-mpmg) Maximum guanine content in a parallel mixed-motif: " << opts.mixed_parallel_max_guanine * 100.0 << "%\n";
    std::cout << "(-mamg) Minimum guanine content in an anti-parallel mixed-motif: " << opts.mixed_antiparallel_min_guanine * 100.0 << "%\n";
    std::cout << "   (-b) Required number of consecutive matches: " << opts.min_block_run << "\n";
    std::cout << "   (-a) Process all sub-matches: " << (opts.all_matches ? "on" : "off") << "\n\n";

    std::cout << "  (-fr) Disregard low-complex regions: " << (opts.filter_repeats ? "on" : "off") << "\n";
    if (opts.filter_repeats) {
        std::cout << " (-mrl) Minimum length to filter a low-complex region: " << opts.min_repeat_length << "\n";
        std::cout << " (-mrp) Maximum repeat period to filter a low-complex region: " << opts.max_repeat_period << "\n";
    }
    if (opts.run_mode != run_mode_t::tpx_search) {
        std::cout << "  (-mf) Merge overlapping features: " << (opts.merge_features ? "on" : "off") << "\n";
    }
    std::cout << "\n";

    std::cout << "   (-o) Output file name: " << opts.output_file << "\n";
    std::cout << "  (-of) Output format: ";
    if (opts.output_format == output_format_t::summary) {
        std::cout << "summary only\n";
    } else if (opts.output_format == output_format_t::bed) {
        std::cout << "BED\n";
    } else if (opts.run_mode == run_mode_t::tpx_search) {
        std::cout << "Triplexator\n";
    } else {
        std::cout << "FASTA\n";
    }
    std::cout << "  (-po) Pretty output: " << (opts.pretty_output ? "on" : "off") << "\n";
    std::cout << "  (-er) Error reference: ";
    switch (opts.error_reference) {
        case error_reference_t::watson_strand:
            std::cout << "watson strand\n\n";
            break;
        case error_reference_t::purine_strand:
            std::cout << "purine strand\n\n";
            break;
        case error_reference_t::last:
        case error_reference_t::third_strand:
            std::cout << "third strand\n\n";
            break;
    }
}
