// MIT License
//
// Copyright (c) 2022-onwards Iñaki Amatria-Barral
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <PATO/tfo_finder.h>
#include <PATO/tpx_finder.h>
#include <PATO/tts_finder.h>
#include <seqan/arg_parse.h>

#include <variant>

static void parse_motifs(const std::string &motifs, pato::options_t &opts) {
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

static std::variant<int, pato::options_t> parse_command_line(int argc,
                                                             char *argv[]) {
  seqan::ArgumentParser parser{"PATO"};

  seqan::setShortDescription(parser, "high PerformAnce TriplexatOr");
  seqan::addUsageLine(
      parser,
      "[options] {-ss tfo_file | -ds tts_file | -ss tfo_file -ds tts_file}");
  seqan::addDescription(
      parser,
      "PATO is a high performance tool for the fast and efficient detection of "
      "triple helices and triplex features in nucleotide sequences. PATO is a "
      "modern altenative to Triplexator and TDF and functions as a drop in "
      "replacement to accelerate the triplex analyses in multicore computers. "
      "Moreover, PATO's efficiency allows a more exhaustive search of the "
      "triplex-forming solution space, and so PATO achieves higher levels of "
      "prediction accuracy in far less time than other tools in the state of "
      "the art.");

  seqan::setDate(parser, "August 2024");
  seqan::setVersion(parser, "v1.0.0");
  seqan::setUrl(parser, "https://github.com/UDC-GAC/PATO");
  seqan::setShortCopyright(parser,
                           "2022-onwards Iñaki Amatria-Barral, Universidade da "
                           "Coruña; released under the MIT License.");
  seqan::setCitation(
      parser, "Iñaki Amatria-Barral, Jorge González-Domínguez and Juan "
              "Touriño. \"PATO: genome-wide prediction of lncRNA-DNA triple "
              "helices\". Bioinformatics, 39(3):btad134, March 2023.");

  seqan::addSection(parser, "Input options");
  seqan::addOption(
      parser,
      seqan::ArgParseOption{
          "ss", "single-strand-file",
          "File in FASTA format that is searched for TFOs (e.g. RNA or DNA).",
          seqan::ArgParseOption::STRING});
  seqan::addOption(
      parser, seqan::ArgParseOption{
                  "ds", "duplex-file",
                  "File in FASTA format that is searched for TTSs (e.g. DNA).",
                  seqan::ArgParseOption::STRING});
  seqan::addSection(parser, "Main options");
  seqan::addOption(
      parser, seqan::ArgParseOption{"l", "lower-length-bound",
                                    "Minimum triplex feature length required.",
                                    seqan::ArgParseOption::INTEGER});
  seqan::addOption(
      parser, seqan::ArgParseOption{
                  "L", "upper-length-bound",
                  "Maximum triplex feature length permitted (disable with -1).",
                  seqan::ArgParseOption::INTEGER});
  seqan::addOption(parser,
                   seqan::ArgParseOption{"e", "error-rate",
                                         "Maximal error rate tolerated in %.",
                                         seqan::ArgParseOption::DOUBLE});
  seqan::addOption(parser,
                   seqan::ArgParseOption{
                       "E", "maximal-error",
                       "Maximal overall error tolerated (disable with -1).",
                       seqan::ArgParseOption::INTEGER});
  seqan::addOption(
      parser, seqan::ArgParseOption{"c", "consecutive-errors",
                                    "Maximum number of consecutive errors.",
                                    seqan::ArgParseOption::INTEGER});
  seqan::addOption(
      parser, seqan::ArgParseOption{"g", "min-guanine",
                                    "Minimal guanine proportion required in %.",
                                    seqan::ArgParseOption::DOUBLE});
  seqan::addOption(
      parser, seqan::ArgParseOption{"G", "max-guanine",
                                    "Maximal guanine proportion allowed in %.",
                                    seqan::ArgParseOption::DOUBLE});
  seqan::addOption(parser,
                   seqan::ArgParseOption{"m", "triplex-motifs",
                                         "Triplex motifs allowed [R,Y,M,P,A].",
                                         seqan::ArgParseOption::STRING});
  seqan::addOption(
      parser, seqan::ArgParseOption{"mpmg", "mixed-parallel-max-guanine",
                                    "Maximum guanine content to consider "
                                    "parallel binding in a mixed-motif in %.",
                                    seqan::ArgParseOption::DOUBLE});
  seqan::addOption(parser, seqan::ArgParseOption{
                               "mamg", "mixed-antiparallel-min-guanine",
                               "Minimum guanine content to consider "
                               "anti-parallel binding in a mixed-motif in %.",
                               seqan::ArgParseOption::DOUBLE});
  seqan::addOption(
      parser, seqan::ArgParseOption{"b", "minimum-block-run",
                                    "Required number of consecutive matches.",
                                    seqan::ArgParseOption::INTEGER});
  seqan::addOption(parser,
                   seqan::ArgParseOption{"a", "all-matches",
                                         "Process and report all sub-matches "
                                         "in addition to the longest match.",
                                         seqan::ArgParseOption::BOOL});
  seqan::addSection(parser, "Filtering options");
  seqan::addOption(parser, seqan::ArgParseOption{
                               "fr", "filter-repeats",
                               "Disregard repeated and low-complex regions.",
                               seqan::ArgParseOption::BOOL});
  seqan::addOption(
      parser, seqan::ArgParseOption{
                  "mrl", "minimum-repeat-length",
                  "Minimum length requirement to filter low-complex regions.",
                  seqan::ArgParseOption::INTEGER});
  seqan::addOption(parser,
                   seqan::ArgParseOption{
                       "mrp", "maximum-repeat-period",
                       "Maximum repeat period to filter low-complex regions.",
                       seqan::ArgParseOption::INTEGER});
  seqan::addOption(
      parser, seqan::ArgParseOption{"mf", "merge-features",
                                    "Merge overlapping features into a cluster "
                                    "and report the spanning region.",
                                    seqan::ArgParseOption::BOOL});
  seqan::addSection(parser, "Output options");
  seqan::addOption(parser,
                   seqan::ArgParseOption{"o", "output", "Output file name",
                                         seqan::ArgParseOption::STRING});
  seqan::addOption(parser, seqan::ArgParseOption{
                               "of", "output-format", "Output format [0,1,2].",
                               seqan::ArgParseOption::INTEGER});
  seqan::addOption(
      parser,
      seqan::ArgParseOption{
          "po", "pretty-output",
          "Indicate matching/mismatching characters with upper/lower case.",
          seqan::ArgParseOption::BOOL});
  seqan::addOption(parser,
                   seqan::ArgParseOption{"er", "error-reference",
                                         "Error reference format [0,1,2].",
                                         seqan::ArgParseOption::INTEGER});
  seqan::addSection(parser, "Performance options");
  seqan::addOption(
      parser,
      seqan::ArgParseOption{
          "cs", "chunk-size",
          "Maximum number of sequences that may be processed simultaneously "
          "(less simultaneous sequences equals less memory usage).",
          seqan::ArgParseOption::INTEGER});

  // input options
  seqan::setDefaultValue(parser, "ss", "(null)");
  seqan::setDefaultValue(parser, "ds", "(null)");

  // main options
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

  // filtering options
  seqan::setDefaultValue(parser, "fr", true);
  seqan::setDefaultValue(parser, "mrl", 10);
  seqan::setDefaultValue(parser, "mrp", 4);
  seqan::setDefaultValue(parser, "mf", false);

  // output options
  seqan::setDefaultValue(parser, "o", "pato");
  seqan::setDefaultValue(parser, "of",
                         static_cast<unsigned>(pato::output_format_t::bed));
  seqan::setDefaultValue(parser, "po", false);
  seqan::setDefaultValue(
      parser, "er",
      static_cast<unsigned>(pato::error_reference_t::watson_strand));

  // performance options
  seqan::setDefaultValue(parser, "cs", 128);

  seqan::ArgumentParser::ParseResult result = seqan::parse(parser, argc, argv);
  switch (result) {
  case seqan::ArgumentParser::PARSE_OK:
    break;
  case seqan::ArgumentParser::PARSE_ERROR:
    return 1;
  case seqan::ArgumentParser::PARSE_HELP:
  case seqan::ArgumentParser::PARSE_VERSION:
  case seqan::ArgumentParser::PARSE_COPYRIGHT:
  case seqan::ArgumentParser::PARSE_WRITE_CTD:
  case seqan::ArgumentParser::PARSE_EXPORT_HELP:
    return 0;
  }

  pato::options_t opts;

  // input options
  seqan::getOptionValue(opts.tfo_file, parser, "ss");
  seqan::getOptionValue(opts.tts_file, parser, "ds");

  // main options
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
  seqan::getOptionValue(opts.all_matches, parser, "a");

  std::string motifs;
  seqan::getOptionValue(motifs, parser, "m");
  parse_motifs(motifs, opts);

  // filtering options
  seqan::getOptionValue(opts.filter_repeats, parser, "fr");
  seqan::getOptionValue(opts.min_repeat_length, parser, "mrl");
  seqan::getOptionValue(opts.max_repeat_period, parser, "mrp");
  seqan::getOptionValue(opts.merge_features, parser, "mf");

  // output options
  seqan::getOptionValue(opts.output_file, parser, "o");
  seqan::getOptionValue(opts.pretty_output, parser, "po");

  {
    unsigned tmp;
    seqan::getOptionValue(tmp, parser, "of");
    opts.output_format = pato::output_format_t(tmp);
    seqan::getOptionValue(tmp, parser, "er");
    opts.error_reference = pato::error_reference_t(tmp);
  }

  // performance options
  seqan::getOptionValue(opts.chunk_size, parser, "cs");

  // check options
  if (opts.tfo_file == "(null)" && opts.tts_file == "(null)") {
    std::cerr << "PATO: at least one type of input file must be supplied\n";
    return 1;
  }
  if (opts.min_length < 5) {
    std::cerr
        << "PATO: the minimum triplex length must be greater or equal to 5\n";
    return 1;
  }
  if (opts.max_length > 1000) {
    std::cerr << "PATO: the maximum triplex length must be smaller or equal to "
                 "1000\n";
    return 1;
  }
  if (opts.error_rate < 0.0 || opts.error_rate > 20.0) {
    std::cerr << "PATO: the error rate must be a value between 0 and 20\n";
    return 1;
  }
  if (opts.max_interruptions > 3) {
    std::cerr << "PATO: the maximum consecutive interruptions must be smaller "
                 "or equal to 3\n";
    return 1;
  }
  if (opts.min_guanine_rate < 0.0 || opts.min_guanine_rate > 100.0) {
    std::cerr << "PATO: the minimum guanine proportion in the triplex target "
                 "site must be a value between 0 and 100\n";
    return 1;
  }
  if (opts.max_guanine_rate < 0.0 || opts.max_guanine_rate > 100.0) {
    std::cerr << "PATO: the maximum guanine proportion in the triplex target "
                 "site must be a value between 0 and 100\n";
    return 1;
  }
  if (opts.min_guanine_rate > opts.max_guanine_rate) {
    std::cerr << "PATO: the maximum guanine proportion cannot be smaller than "
                 "the minimum guanine proportion\n";
    return 1;
  }
  if (opts.mixed_parallel_max_guanine < 0.0 ||
      opts.mixed_parallel_max_guanine > 100.0) {
    std::cerr << "PATO: the maximum guanine proportion in a parallel mixed "
                 "motif must be a value between 0 and 100\n";
    return 1;
  }
  if (opts.mixed_antiparallel_min_guanine < 0.0 ||
      opts.mixed_antiparallel_min_guanine > 100.0) {
    std::cerr << "PATO: the minimum guanine proportion in an anti-parallel "
                 "mixed motif must be a value between 0 and 100\n";
    return 1;
  }
  if (opts.error_reference >= pato::error_reference_t::last) {
    std::cerr << "PATO: unknown error reference format\n";
    return 1;
  }
  if (opts.output_format >= pato::output_format_t::last) {
    std::cerr << "PATO: unknown output format\n";
    return 1;
  }
  if (opts.chunk_size == 0) {
    std::cerr << "PATO: the chunk size must be a positive number\n";
    return 1;
  }

  unsigned tolerated_error =
      static_cast<unsigned>(std::floor(opts.error_rate * opts.min_length));
  if (opts.min_block_run > opts.min_length - 2 * tolerated_error) {
    std::cerr << "PATO: block match too large given minimum length constraint "
                 "and error rate\n";
    return 1;
  }

  // prepare options
  if (opts.tts_file == "(null)") {
    opts.run_mode = pato::run_mode_t::tfo_search;
  } else if (opts.tfo_file == "(null)") {
    opts.run_mode = pato::run_mode_t::tts_search;
  } else {
    opts.run_mode = pato::run_mode_t::tpx_search;
  }

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
      opts.maximal_error =
          std::min(opts.maximal_error,
                   static_cast<int>(opts.error_rate * opts.max_length));
    }
  }

  return opts;
}

static int handle_tfo_search(const pato::options_t &opts) {
  switch (pato::find_tfo_motifs(opts)) {
  case pato::find_tfo_motifs_result::success:
    return 0;
  case pato::find_tfo_motifs_result::cannot_open_tfo_file:
    std::cerr << "PATO: error opening input file '" << opts.tfo_file << "'\n";
    break;
  case pato::find_tfo_motifs_result::cannot_create_output_file:
    std::cerr << "PATO: error opening output file '" << opts.output_file
              << "'\n";
    break;
  }
  return 1;
}

static int handle_tts_search(const pato::options_t &opts) {
  switch (pato::find_tts_motifs(opts)) {
  case pato::find_tts_motifs_result::success:
    return 0;
  case pato::find_tts_motifs_result::cannot_open_tts_file:
    std::cerr << "PATO: error opening input file '" << opts.tts_file << "'\n";
    break;
  case pato::find_tts_motifs_result::cannot_create_output_file:
    std::cerr << "PATO: error opening output file '" << opts.output_file
              << "'\n";
    break;
  }
  return 1;
}

static int handle_tpx_search(const pato::options_t &opts) {
  switch (pato::find_tpxes(opts)) {
  case pato::find_tpx_result::success:
    return 0;
  case pato::find_tpx_result::cannot_open_tfo_file:
    std::cerr << "PATO: error opening input file '" << opts.tfo_file << "'\n";
    break;
  case pato::find_tpx_result::cannot_open_tts_file:
    std::cerr << "PATO: error opening input file '" << opts.tts_file << "'\n";
    break;
  case pato::find_tpx_result::cannot_create_output_file:
    std::cerr << "PATO: error opening output file '" << opts.output_file
              << "'\n";
    break;
  }
  return 1;
}

namespace {

template <typename... Tys> struct visitors_t : Tys... {
  using Tys::operator()...;
};
template <typename... Tys> visitors_t(Tys...) -> visitors_t<Tys...>;

} // namespace

int main(int argc, char *argv[]) {
  auto result = parse_command_line(argc, argv);
  return std::visit(visitors_t{[](int return_code) { return return_code; },
                               [](pato::options_t &opts) {
                                 int result;
                                 switch (opts.run_mode) {
                                 case pato::run_mode_t::tfo_search:
                                   result = handle_tfo_search(opts);
                                   break;
                                 case pato::run_mode_t::tts_search:
                                   result = handle_tts_search(opts);
                                   break;
                                 case pato::run_mode_t::tpx_search:
                                   result = handle_tpx_search(opts);
                                   break;
                                 }
                                 return result;
                               }},
                    result);
}
