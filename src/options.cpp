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
        std::cout << "  (-ss) \e[1mTFO file:\e[0m " << opts.tfo_file << "\n";
    }
    if (opts.run_mode != run_mode_t::tfo_search) {
        std::cout << "  (-ds) \e[1mTTS file:\e[0m " << opts.tts_file << "\n";
    }
    std::cout << "\n";

    std::cout << "   (-l) \e[1mMinimum triplex length:\e[0m " << opts.min_length << "\n";
    std::cout << "   (-L) \e[1mMaximum triplex length:\e[0m " << opts.max_length << "\n";
    std::cout << "   (-e) \e[1mError rate:\e[0m " << opts.error_rate * 100.0 << "%\n";
    std::cout << "   (-E) \e[1mMaximal error rate:\e[0m " << opts.maximal_error << "\n";
    std::cout << "   (-c) \e[1mMaximum consecutive errors:\e[0m " << opts.max_interruptions << "\n";
    std::cout << "   (-g) \e[1mMinimal guanine proportion:\e[0m " << opts.min_guanine_rate * 100.0 << "%\n";
    std::cout << "   (-G) \e[1mMaximal guanine proportion:\e[0m " << opts.max_guanine_rate * 100.0 << "%\n";
    if (opts.run_mode != run_mode_t::tts_search) {
        std::cout << "   (-m) \e[1mTriplex motifs allowed\e[0m \n";
        std::cout << "        - \e[1mTC:\e[0m " << (opts.tc_motif ? "on" : "off") << "\n";
        std::cout << "        - \e[1mGA:\e[0m " << (opts.ga_motif ? "on" : "off") << "\n";
        std::cout << "        - \e[1mParallel GT:\e[0m " << (opts.gt_p_motif ? "on" : "off") << "\n";
        std::cout << "        - \e[1mAntiprallel GT:\e[0m " << (opts.gt_a_motif ? "on" : "off") << "\n"; 
    }
    std::cout << "(-mpmg) \e[1mMaximum guanine content in a parallel mixed-motif:\e[0m " << opts.mixed_parallel_max_guanine * 100.0 << "%\n";
    std::cout << "(-mamg) \e[1mMinimum guanine content in an anti-parallel mixed-motif:\e[0m " << opts.mixed_antiparallel_min_guanine * 100.0 << "%\n";
    std::cout << "   (-b) \e[1mRequired number of consecutive matches:\e[0m " << opts.min_block_run << "\n";
    std::cout << "   (-a) \e[1mProcess all sub-matches:\e[0m " << (opts.all_matches ? "on" : "off") << "\n\n";

    std::cout << "  (-fr) \e[1mDisregard low-complex regions:\e[0m " << (opts.filter_repeats ? "on" : "off") << "\n";
    if (opts.filter_repeats) {
        std::cout << " (-mrl) \e[1mMinimum length to filter a low-complex region:\e[0m " << opts.min_repeat_length << "\n";
        std::cout << " (-mrp) \e[1mMaximum repeat period to filter a low-complex region:\e[0m " << opts.max_repeat_period << "\n";
    }
    if (opts.run_mode != run_mode_t::tpx_search) {
        std::cout << "  (-mf) \e[1mMerge overlapping features:\e[0m " << (opts.merge_features ? "on" : "off") << "\n";
    }
    std::cout << "\n";

    std::cout << "   (-o) \e[1mOutput file name:\e[0m " << opts.output_file << "\n";
    std::cout << "  (-of) \e[1mOutput format:\e[0m ";
    if (opts.output_format == output_format_t::summary) {
        std::cout << "summary only\n";
    } else if (opts.output_format == output_format_t::bed) {
        std::cout << "BED\n";
    } else if (opts.run_mode == run_mode_t::tpx_search) {
        std::cout << "Triplexator\n";
    } else {
        std::cout << "FASTA\n";
    }
    std::cout << "  (-po) \e[1mPretty output:\e[0m " << (opts.pretty_output ? "on" : "off") << "\n";
    std::cout << "  (-er) \e[1mError reference:\e[0m ";
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

    if (opts.run_mode != run_mode_t::tfo_search) {
        std::cout << "  (-cs) \e[1mTTS window size\e[0m: " << opts.chunk_size << "\n\n";
    }
}
