#ifndef _COMMAND_LINE_PARSER_HPP_
#define _COMMAND_LINE_PARSER_HPP_

#include <string>

#include "options.hpp"

#include <seqan/misc/misc_cmdparser.h>

void parse_motifs(options& opts, const std::string& motifs)
{
    opts.tc_motif = false;
    opts.ga_motif = false;
    opts.gt_p_motif = false;
    opts.gt_a_motif = false;

    if (motifs.find('A') != std::string::npos) {
        opts.ga_motif = true;
        opts.gt_a_motif = true;
    }
    if (motifs.find('M') != std::string::npos) {
        opts.gt_p_motif = true;
        opts.gt_a_motif = true;
    }
    if (motifs.find('P') != std::string::npos) {
        opts.tc_motif = true;
        opts.gt_p_motif = true;
    }
    if (motifs.find('R') != std::string::npos) {
        opts.ga_motif = true;
    }
    if (motifs.find('Y') != std::string::npos) {
        opts.tc_motif = true;
    }
}

bool parse_command_line(options& opts, int argc, const char *argv[])
{
    seqan::CommandLineParser parser;

    seqan::addVersionLine(parser, "PATO v0.0.0 (Compiled on " __DATE__ ")");
    seqan::addUsageLine(parser, "[options] tfo-file tts-file output-file");

    seqan::addOption(parser, seqan::CommandLineOption("l", "lower-length-bound", "minimum triplex feature length required", seqan::OptionType::Integer, 16));
    seqan::addOption(parser, seqan::CommandLineOption("L", "upper-length-bound", "maximum triplex feature length permitted (disable with -1)", seqan::OptionType::Integer, 30));
    seqan::addOption(parser, seqan::CommandLineOption("e", "error-rate", "set the maximal error rate tolerated in %", seqan::OptionType::Double, 5.0));
    seqan::addOption(parser, seqan::CommandLineOption("E", "maximal-error", "set the maximal overall error tolerated (disable with -1)", seqan::OptionType::Integer, -1));
    seqan::addOption(parser, seqan::CommandLineOption("c", "consecutive-errors", "maximum number of consecutive errors", seqan::OptionType::Integer, 1));
    seqan::addOption(parser, seqan::CommandLineOption("g", "min-guanine", "set the minimal guanine proportion required in %", seqan::OptionType::Double, 10.0));
    seqan::addOption(parser, seqan::CommandLineOption("G", "max-guanine", "set the maximal guanine proportion allowed in %", seqan::OptionType::Double, 100.0));
    seqan::addOption(parser, seqan::CommandLineOption("m", "triplex-motifs", "triplex motifs allowed [R,Y,M,P,A]", seqan::OptionType::String, "R,Y,M,P,A"));
    seqan::addOption(parser, seqan::CommandLineOption("mpmg", "mixed-parallel-max-guanine", "maximum guanine content, in %, to consider parallel binding in a mixed-motif (GT)", seqan::OptionType::Double, 100.0));
    seqan::addOption(parser, seqan::CommandLineOption("mamg", "mixed-antiparallel-min-guanine", "minimum guanine content, in %, to consider anti-parallel binding in a mixed-motif (GT)", seqan::OptionType::Double, 0.0));
    seqan::addOption(parser, seqan::CommandLineOption("b", "minimum-block-run", "required number of consecutive matches", seqan::OptionType::Integer, 1));
    seqan::addOption(parser, seqan::CommandLineOption("a", "all-matches", "process and report all sub-matches in addition to the longest match", seqan::OptionType::Boolean, false));
    seqan::addOption(parser, seqan::CommandLineOption("dd", "detect-duplicates", "indicates wheter and how duplicates should be detected", seqan::OptionType::Integer, 0));
    seqan::addOption(parser, seqan::CommandLineOption("ssd", "same-sequence-duplicates", "wheter to count a feature copy in the same sequence as duplicate or not", seqan::OptionType::String, "on"));

    seqan::addOption(parser, seqan::CommandLineOption("fr", "filter-repeats", "disregards repeated and low-complex regions if enabled", seqan::OptionType::String, "on"));
    seqan::addOption(parser, seqan::CommandLineOption("mrl", "minimum-repeat-length", "minimum length requirement for low-complex regions to be filtered", seqan::OptionType::Integer, 10));
    seqan::addOption(parser, seqan::CommandLineOption("mrp", "maximum-repeat-period", "maximum repeat period for low-complex regions to be filtered", seqan::OptionType::Integer, 4));
    seqan::addOption(parser, seqan::CommandLineOption("dc", "duplicate-cutoff", "disregard feature if it occurs more often than this cutoff (disable with -1)", seqan::OptionType::Integer, -1));
    seqan::addOption(parser, seqan::CommandLineOption("dl", "duplicate-locations", "report the location of duplicates", seqan::OptionType::Boolean, false));

    if (!seqan::parse(parser, argc, argv)) {
        return false;
    }

    // arguments
    opts.tfo_file = seqan::getArgumentValue(parser, 0);
    opts.tts_file = seqan::getArgumentValue(parser, 1);
    opts.output_file = seqan::getArgumentValue(parser, 2);

    if (opts.tfo_file == parser._null
        || opts.tts_file == parser._null
        || opts.output_file == parser._null) {
        seqan::help(parser);
        return false;
    }

    // options
    seqan::getOptionValue(parser, "l", opts.min_length);
    seqan::getOptionValue(parser, "L", opts.max_length);
    seqan::getOptionValue(parser, "e", opts.error_rate);
    seqan::getOptionValue(parser, "E", opts.maximal_error);
    seqan::getOptionValue(parser, "c", opts.max_interruptions);
    seqan::getOptionValue(parser, "g", opts.min_guanine_rate);
    seqan::getOptionValue(parser, "G", opts.max_guanine_rate);
    seqan::getOptionValue(parser, "mpmg", opts.mixed_parallel_max_guanine);
    seqan::getOptionValue(parser, "mamg", opts.mixed_antiparallel_min_guanine);
    seqan::getOptionValue(parser, "b", opts.min_block_run);
    seqan::getOptionValue(parser, "dd", opts.detect_duplicates);

    seqan::getOptionValue(parser, "mrl", opts.min_repeat_length);
    seqan::getOptionValue(parser, "mrp", opts.max_repeat_period);
    seqan::getOptionValue(parser, "dc", opts.duplicate_cutoff);

    // flags
    opts.all_matches = seqan::isSet(parser, "a");
    opts.report_duplicate_locations = seqan::isSet(parser, "dl");
    auto is_flag_set = [&parser](const std::string& name) {
        std::string value;
        seqan::getOptionValue(parser, name, value);
        return value == "on" || value == "yes" || value == "ON" || value == "1"
               || value == "YES" || value == "true" || value == "TRUE";
    };
    opts.filter_repeats = is_flag_set("fr");
    opts.same_sequence_duplicates = is_flag_set("ssd");

    // motifs
    std::string motifs;
    seqan::getOptionValue(parser, "m", motifs);
    parse_motifs(opts, motifs);

    // check options

    // prepare options
    opts.error_rate /= 100.0;
    opts.min_guanine_rate /= 100.0;
    opts.max_guanine_rate /= 100.0;
    opts.mixed_parallel_max_guanine /= 100.0;
    opts.mixed_antiparallel_min_guanine /= 100.0;

    if (opts.error_rate == 0.0 || opts.maximal_error == 0) {
        opts.max_interruptions = 0;
    }

    if (opts.max_length >= opts.min_length) {
        if (opts.maximal_error < 0) {
            opts.maximal_error = static_cast<int>(opts.error_rate * opts.max_length);
        } else {
            opts.maximal_error = std::min(opts.maximal_error,
                                          static_cast<int>(opts.error_rate * opts.max_length));
        }
    }

    return true;
}

#endif
