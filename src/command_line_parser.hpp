#ifndef _COMMAND_LINE_PARSER_HPP_
#define _COMMAND_LINE_PARSER_HPP_

#include <string>

#include "seqan.hpp"
#include "options.hpp"

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
    seqan::ArgumentParser parser("PATO");

    seqan::setShortDescription(parser, "PArallel TriplexatOr");
    seqan::addUsageLine(parser, "[options] tfo-file tts-file output-file");
    seqan::addDescription(parser, "PATO is a high performance tool for detectin"
                                  "g nucleic acid triple helices and triplex fe"
                                  "atures in nucleotide sequences. PATO is base"
                                  "d on Triplexator and functions nearly as a d"
                                  "rop in replacement to accelerate the triplex"
                                  " analyses in multicore computing clusters an"
                                  "d supercomputing facilities.");

    seqan::setDate(parser, "June 2022");
    seqan::setVersion(parser, "v0.0.0");
    seqan::setUrl(parser, "https://github.com/amatria/pato");
    seqan::setShortCopyright(parser, "2022 Iñaki Amatria-Barral.");

    seqan::addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "tfo-file"));
    seqan::addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "tts-file"));
    seqan::addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::OUTPUT_FILE, "output-file"));

    seqan::addOption(parser, seqan::ArgParseOption("l", "lower-length-bound", "Minimum triplex feature length required.", seqan::ArgParseOption::INTEGER));
    seqan::addOption(parser, seqan::ArgParseOption("L", "upper-length-bound", "Maximum triplex feature length permitted (disable with -1).", seqan::ArgParseOption::INTEGER));
    seqan::addOption(parser, seqan::ArgParseOption("e", "error-rate", "Set the maximal error rate tolerated in %.", seqan::ArgParseOption::DOUBLE));
    seqan::addOption(parser, seqan::ArgParseOption("E", "maximal-error", "Set the maximal overall error tolerated (disable with -1).", seqan::ArgParseOption::INTEGER));
    seqan::addOption(parser, seqan::ArgParseOption("c", "consecutive-errors", "Maximum number of consecutive errors.", seqan::ArgParseOption::INTEGER));
    seqan::addOption(parser, seqan::ArgParseOption("g", "min-guanine", "Set the minimal guanine proportion required in %.", seqan::ArgParseOption::DOUBLE));
    seqan::addOption(parser, seqan::ArgParseOption("G", "max-guanine", "Set the maximal guanine proportion allowed in %.", seqan::ArgParseOption::DOUBLE));
    seqan::addOption(parser, seqan::ArgParseOption("m", "triplex-motifs", "Triplex motifs allowed [R,Y,M,P,A].", seqan::ArgParseOption::STRING));
    seqan::addOption(parser, seqan::ArgParseOption("mpmg", "mixed-parallel-max-guanine", "Maximum guanine content to consider parallel binding in a mixed-motif in %.", seqan::ArgParseOption::DOUBLE));
    seqan::addOption(parser, seqan::ArgParseOption("mamg", "mixed-antiparallel-min-guanine", "Minimum guanine content to consider anti-parallel binding in a mixed-motif in %.", seqan::ArgParseOption::DOUBLE));
    seqan::addOption(parser, seqan::ArgParseOption("b", "minimum-block-run", "Required number of consecutive matches.", seqan::ArgParseOption::INTEGER));
    seqan::addOption(parser, seqan::ArgParseOption("a", "all-matches", "Process and report all sub-matches in addition to the longest match.", seqan::ArgParseOption::BOOL));
    seqan::addOption(parser, seqan::ArgParseOption("dd", "detect-duplicates", "Indicates wheter and how duplicates should be detected.", seqan::ArgParseOption::INTEGER));
    seqan::addOption(parser, seqan::ArgParseOption("ssd", "same-sequence-duplicates", "Wheter to count a feature copy in the same sequence as duplicate or not.", seqan::ArgParseOption::STRING));
    seqan::addOption(parser, seqan::ArgParseOption("fr", "filter-repeats", "Disregards repeated and low-complex regions if enabled.", seqan::ArgParseOption::STRING));
    seqan::addOption(parser, seqan::ArgParseOption("mrl", "minimum-repeat-length", "Minimum length requirement for low-complex regions to be filtered.", seqan::ArgParseOption::INTEGER));
    seqan::addOption(parser, seqan::ArgParseOption("mrp", "maximum-repeat-period", "Maximum repeat period for low-complex regions to be filtered.", seqan::ArgParseOption::INTEGER));
    seqan::addOption(parser, seqan::ArgParseOption("dc", "duplicate-cutoff", "Disregard feature if it occurs more often than this cutoff (disable with -1).", seqan::ArgParseOption::INTEGER));
    seqan::addOption(parser, seqan::ArgParseOption("dl", "duplicate-locations", "Report the location of duplicates.", seqan::ArgParseOption::BOOL));

    seqan::setDefaultValue(parser, "l", 16);
    seqan::setDefaultValue(parser, "L", 30);
    seqan::setDefaultValue(parser, "e", 5.0);
    seqan::setDefaultValue(parser, "E", -1);
    seqan::setDefaultValue(parser, "c", 1);
    seqan::setDefaultValue(parser, "g", 10.0);
    seqan::setDefaultValue(parser, "G", 100.0);
    seqan::setDefaultValue(parser, "m", "R,Y,M,P,A");
    seqan::setDefaultValue(parser, "mpmg", 100.0);
    seqan::setDefaultValue(parser, "mamg", 0.0);
    seqan::setDefaultValue(parser, "b", 1);
    seqan::setDefaultValue(parser, "a", false);
    seqan::setDefaultValue(parser, "dd", 0);
    seqan::setDefaultValue(parser, "ssd", "on");
    seqan::setDefaultValue(parser, "fr", "on");
    seqan::setDefaultValue(parser, "mrl", 10);
    seqan::setDefaultValue(parser, "mrp", 4);
    seqan::setDefaultValue(parser, "dc", -1);
    seqan::setDefaultValue(parser, "dl", false);

    if (seqan::parse(parser, argc, argv) != seqan::ArgumentParser::PARSE_OK) {
        return false;
    }

    opts.tfo_file = seqan::getArgumentValue(seqan::getArgument(parser, 0));
    opts.tts_file = seqan::getArgumentValue(seqan::getArgument(parser, 1));
    opts.output_file = seqan::getArgumentValue(seqan::getArgument(parser, 2));

    seqan::getOptionValue(opts.min_length, parser, "l");
    seqan::getOptionValue(opts.max_length, parser, "L");
    seqan::getOptionValue(opts.error_rate, parser, "e");
    seqan::getOptionValue(opts.maximal_error, parser, "E");
    seqan::getOptionValue(opts.max_interruptions, parser, "c");
    seqan::getOptionValue(opts.min_guanine_rate, parser, "g");
    seqan::getOptionValue(opts.max_guanine_rate, parser, "G");
    seqan::getOptionValue(opts.mixed_parallel_max_guanine, parser, "mpmg");
    seqan::getOptionValue(opts.mixed_antiparallel_min_guanine, parser, "mamg");
    seqan::getOptionValue(opts.min_block_run, parser, "b");
    seqan::getOptionValue(opts.detect_duplicates, parser, "dd");
    seqan::getOptionValue(opts.min_repeat_length, parser, "mrl");
    seqan::getOptionValue(opts.max_repeat_period, parser, "mrp");
    seqan::getOptionValue(opts.duplicate_cutoff, parser, "dc");

    opts.all_matches = seqan::isSet(parser, "a");
    opts.report_duplicate_locations = seqan::isSet(parser, "dl");

    auto is_flag_set = [&parser](const std::string& name) -> bool {
        std::string value;
        seqan::getOptionValue(value, parser, name);
        return value == "on" || value == "yes" || value == "ON" || value == "1"
               || value == "YES" || value == "true" || value == "TRUE";
    };

    opts.filter_repeats = is_flag_set("fr");
    opts.same_sequence_duplicates = is_flag_set("ssd");

    // motifs
    std::string motifs;
    seqan::getOptionValue(motifs, parser, "m");
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
