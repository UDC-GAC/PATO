#include "output_writer.hpp"

#include <cctype>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iostream>

#include <seqan/modifier.h>
#include <seqan/sequence.h>

#include "triplex_enums.hpp"
#include "triplex_match.hpp"
#include "triplex_pattern.hpp"

seqan::CharString triplex_alignment_string(match_t& match,
                                           motif_set_t& tfo_motifs,
                                           motif_set_t& tts_motifs)
{
    std::ostringstream alignment;

    motif_t tfo_motif(seqan::host(tfo_motifs[match.tfoNo]),
                      match.oBegin,
                      match.oEnd,
                      match.parallel,
                      tfo_motifs[match.tfoNo].seqNo,
                      true,
                      match.motif);
    motif_t tts_motif(seqan::host(tts_motifs[match.ttsNo]),
                      match.dBegin,
                      match.dEnd,
                      match.parallel,
                      match.ttsSeqNo,
                      false,
                      match.strand);

    seqan::CharString tfo_ps = seqan::prettyString(tfo_motif);
    seqan::CharString tts_ps = seqan::prettyString(tts_motif);

    seqan::CharString opp(tts_ps);
    seqan::complement(opp);

    alignment << "\n";
    if (match.strand == '-') {
        seqan::reverse(opp);
        seqan::reverse(tts_ps);
        alignment << "     5'- " << opp << " -3'" << "\n";
        alignment << "TTS: 3'- " << tts_ps << " -5'" << "\n";
        alignment << "         ";

        auto tts_it = seqan::begin(tts_motif);
        auto tts_end = seqan::end(tts_motif);
        auto tfo_it = seqan::begin(tfo_motif);
        auto tfo_end = seqan::end(tfo_motif);
        while (tts_end != tts_it && tfo_end != tfo_it) {
            tts_end--;
            tfo_end--;
            if (*tts_end == *tfo_end) {
                alignment << "|";
            } else {
                alignment << "*";
            }
        }
        alignment << "\n";

        if (!tfo_motifs[match.tfoNo].parallel) {
            alignment << "TFO: 5'- " << tfo_ps << " -3'" << "\n";
        } else {
            seqan::reverse(tfo_ps);
            alignment << "TFO: 3'- " << tfo_ps << " -5'" << "\n";
        }
    } else {
        if (!tfo_motifs[match.tfoNo].parallel) {
            seqan::reverse(tfo_ps);
            alignment << "TFO: 3'- " << tfo_ps << " -5'" << "\n";
        } else {
            alignment << "TFO: 5'- " << tfo_ps << " -3'" << "\n";
        }
        alignment << "         ";

        auto tts_it = seqan::begin(tts_motif);
        auto tts_end = seqan::end(tts_motif);
        auto tfo_it = seqan::begin(tfo_motif);
        auto tfo_end = seqan::end(tfo_motif);
        while (tts_it != tts_end && tfo_it != tfo_end) {
            if (*tts_it == *tfo_it) {
                alignment << "|";
            } else {
                alignment << "*";
            }
            tts_it++;
            tfo_it++;
        }
        alignment << "\n";

        alignment << "TTS: 5'- " << tts_ps << " -3'" << "\n";
        alignment << "     3'- " << opp << " -5'" << "\n";
    }

    return alignment.str();
}

seqan::CharString triplex_error_string(match_t& match,
                                       motif_set_t& tfo_motifs,
                                       motif_set_t& tts_motifs,
                                       const options& opts)
{
    std::ostringstream errors;

    motif_t tfo_motif(seqan::host(tfo_motifs[match.tfoNo]),
                      match.oBegin,
                      match.oEnd,
                      match.parallel,
                      tfo_motifs[match.tfoNo].seqNo,
                      true,
                      match.motif);
    motif_t tts_motif(seqan::host(tts_motifs[match.ttsNo]),
                      match.dBegin,
                      match.dEnd,
                      match.parallel,
                      match.ttsSeqNo,
                      false,
                      match.strand);

    seqan::CharString tfo_ps = seqan::prettyString(tfo_motif);
    seqan::CharString tts_ps = seqan::prettyString(tts_motif);

    if (opts.error_reference == error_reference_t::purine_strand) {
        if (!tfo_motifs[match.tfoNo].parallel) {
            seqan::reverse(tfo_ps);
        }
    } else if (opts.error_reference == error_reference_t::third_strand) {
        if (!tfo_motifs[match.tfoNo].parallel) {
            seqan::reverse(tts_ps);
        }
    } else {
        if (match.strand == '-') {
            seqan::reverse(tts_ps);
            if (tfo_motifs[match.tfoNo].parallel) {
                seqan::reverse(tfo_ps);
            }
        } else {
            if (!tfo_motifs[match.tfoNo].parallel) {
                seqan::reverse(tfo_ps);
            }
        }
    }

    auto tts_it = seqan::begin(tts_motif);
    auto tts_end = seqan::end(tts_motif);
    auto tfo_it = seqan::begin(tfo_motif);
    auto tfo_end = seqan::end(tfo_motif);
    unsigned int i = 0;

    if (opts.error_reference == error_reference_t::purine_strand
        || (opts.error_reference == error_reference_t::watson_strand && match.strand == '+')
        || (opts.error_reference == error_reference_t::third_strand && tfo_motifs[match.tfoNo].parallel)) {
        while (tts_it != tts_end && tfo_it != tfo_end) {
            if (*tts_it != *tfo_it) {
                if (!isupper(seqan::value(tts_ps, i))
                    && !isupper(seqan::value(tfo_ps, i))) {
                    errors << "b" << i;
                } else if (!isupper(seqan::value(tts_ps, i))) {
                    errors << "d" << i;
                } else if (!isupper(seqan::value(tfo_ps, i))) {
                    errors << "o" << i;
                } else {
                    errors << "t" << i;
                }
            }
            tts_it++;
            tfo_it++;
            i++;
        }
    } else {
        while (tts_end != tts_it && tfo_end != tfo_it) {
            tts_end--;
            tfo_end--;
            if (*tts_end != *tfo_end) {
                if (!isupper(seqan::value(tts_ps, i))
                    && !isupper(seqan::value(tfo_ps, i))) {
                    errors << "b" << i;
                } else if (!isupper(seqan::value(tts_ps, i))) {
                    errors << "d" << i;
                } else if (!isupper(seqan::value(tfo_ps, i))) {
                    errors << "o" << i;
                } else {
                    errors << "t" << i;
                }
            }
            i++;
        }
    }

    return errors.str();
}

void print_motifs(motif_set_t& motifs,
                  name_set_t& names,
                  const options& opts)
{
    if (opts.output_format == output_format_t::summary) {
        return;
    }

    seqan::CharString output_file_name;
    seqan::append(output_file_name, opts.output_file);
    seqan::append(output_file_name, ".out");

    std::ofstream output_file(seqan::toCString(output_file_name),
                              std::ios_base::out);
    if (!output_file) {
        std::cerr << "PATO: error opening output file '"
                  << seqan::toCString(output_file_name) << "'\n";
        return;
    }

    if (opts.output_format == output_format_t::bed) {
        if (opts.run_mode == run_mode_t::tfo_search) {
            output_file << "# Sequence-ID\tStart\tEnd\tScore\tMotif\tError-rate"
                           "\tErrors\tGuanine-rate\tDuplicates\tTFO\tDuplicate "
                           "locations\n";
        } else {
            output_file << "# Duplex-ID\tStart\tEnd\tScore\tStrand\tError-rate"
                           "\tErrors\tGuanine-rate\tDuplicates\tTTS\tDuplicate "
                           "locations\n";
        }

    }

    unsigned int counter = 1;
    for (auto& m : motifs) {
        switch (opts.output_format) {
            case output_format_t::bed:
                output_file << names[seqan::getSequenceNo(m)] << "\t"
                            << seqan::beginPosition(m) << "\t"
                            << seqan::endPosition(m) << "\t"
                            << seqan::score(m) << "\t"
                            << m.motif << "\t"
                            << std::setprecision(2) << 1.0 - seqan::score(m) / (seqan::endPosition(m) - seqan::beginPosition(m)) << "\t"
                            << seqan::errorString(m) << "\t"
                            << seqan::guanineRate(m) << "\t"
                            << seqan::duplicates(m) << "\t"
                            << (opts.pretty_output ? seqan::prettyString(m) : seqan::outputString(m)) << "\t";
                if (!opts.report_duplicate_locations
                    || seqan::duplicates(m) < 1
                    || opts.duplicate_cutoff <= seqan::duplicates(m)) {
                    output_file << "-";
                } else {
                    for (int d = 0; d < seqan::duplicates(m); d++) {
                        output_file << names[seqan::getDuplicateAt(m, d).i1] << ":"
                                    << seqan::getDuplicateAt(m, d).i2 << "-"
                                    << seqan::getDuplicateAt(m, d).i2 + seqan::length(m) << ";";
                    }
                }
                output_file << "\n";
                break;
            case output_format_t::triplex:
                output_file << ">"
                            << names[seqan::getSequenceNo(m)] << "_"
                            << counter++ << "\t"
                            << seqan::beginPosition(m) << "-"
                            << seqan::endPosition(m) << " "
                            << m.motif << "\t"
                            << seqan::score(m) << "\t"
                            << seqan::errorString(m) << "\t"
                            << seqan::duplicates(m) << "\t"
                            << seqan::guanineRate(m) << "\t";
                if (!opts.report_duplicate_locations
                    || seqan::duplicates(m) < 1
                    || opts.duplicate_cutoff <= seqan::duplicates(m)) {
                    output_file << "-";
                } else {
                    for (int d = 0; d < seqan::duplicates(m); d++) {
                        output_file << names[seqan::getDuplicateAt(m, d).i1] << ":"
                                    << seqan::getDuplicateAt(m, d).i2 << "-"
                                    << seqan::getDuplicateAt(m, d).i2 + seqan::length(m) << ";";
                    }
                }
                output_file << "\n"
                            << (opts.pretty_output ? seqan::prettyString(m) : seqan::outputString(m))
                            << "\n";
                break;
            default:
                break;
        }
    }
}

void print_summary(motif_potential_set_t& potentials,
                   name_set_t& names,
                   const options& opts)
{
    seqan::CharString output_file_name;
    seqan::append(output_file_name, opts.output_file);
    seqan::append(output_file_name, ".summary");

    std::ofstream output_file(seqan::toCString(output_file_name),
                              std::ios_base::out);
    if (!output_file) {
        std::cerr << "PATO: error opening output file '"
                  << seqan::toCString(output_file_name) << "'\n";
        return;
    }

    if (opts.run_mode == run_mode_t::tfo_search) {
        output_file << "# Sequence-ID\tTFOs (abs)\tTFOs (rel)\tGA (abs)\tGA (rel)\t"
                    "TC (abs)\tTC (rel)\tGT (abs)\tGT (rel)\n";
    } else {
        output_file << "# Duplex-ID\tTTSs (abs)\tTTSs (rel)\n";
    }

    for (auto& potential : potentials) {
        if (seqan::hasCount(potential)) {
            output_file << names[seqan::getKey(potential)] << "\t"
                        << seqan::getCounts(potential) << "\t"
                        << std::setprecision(3) << seqan::getCounts(potential) / seqan::getNorm(potential);
            if (opts.run_mode == run_mode_t::tfo_search) {
                output_file << "\t" << seqan::getCount(potential, 'R') << "\t"
                            << std::setprecision(3) << seqan::getCount(potential, 'R') / seqan::getNorm(potential) << "\t"
                            << seqan::getCount(potential, 'Y') << "\t"
                            << std::setprecision(3) << seqan::getCount(potential, 'Y') / seqan::getNorm(potential) << "\t"
                            << seqan::getCount(potential, 'M') << "\t"
                            << std::setprecision(3) << seqan::getCount(potential, 'M') /seqan:: getNorm(potential) << "\t";
            }
            output_file << "\n";
        }
    }
}

#if !defined(_OPENMP)
void print_triplex_pairs(match_set_t& matches,
                         motif_set_t& tfo_motifs,
                         name_set_t& tfo_names,
                         motif_set_t& tts_motifs,
                         name_set_t& tts_names,
                         const options& opts)
#else
void print_triplex_pairs(match_set_set_t& matches,
                         motif_set_t& tfo_motifs,
                         name_set_t& tfo_names,
                         motif_set_t& tts_motifs,
                         name_set_t& tts_names,
                         const options& opts)
#endif
{
    if (opts.output_format == output_format_t::summary) {
        return;
    }

    seqan::CharString output_file_name;
    seqan::append(output_file_name, opts.output_file);
    seqan::append(output_file_name, ".out");

    std::ofstream output_file(seqan::toCString(output_file_name),
                              std::ios_base::out);
    if (!output_file) {
        std::cerr << "PATO: error opening output file '"
                  << seqan::toCString(output_file_name) << "'\n";
        return;
    }

    output_file << "# Sequence-ID\tTFO start\tTFO end\tDuplex-ID\tTTS start\tTT"
                   "S end\tScore\tError-rate\tErrors\tMotif\tStrand\tOrientatio"
                   "n\tGuanine-rate\n";
#if !defined(_OPENMP)
    for (auto& match: matches) {
#else
    for (auto& local_matches : matches) {
        for (auto& match : local_matches) {
#endif
            auto tfo_seq_id = tfo_motifs[match.tfoNo].seqNo;
            auto tts_seq_id = match.ttsSeqNo;

            output_file << tfo_names[tfo_seq_id] << "\t"
                        << match.oBegin << "\t"
                        << match.oEnd << "\t"
                        << tts_names[tts_seq_id] << "\t"
                        << match.dBegin << "\t"
                        << match.dEnd << "\t"
                        << match.mScore << "\t"
                        << std::setprecision(2) << 1.0 - match.mScore / (match.dEnd - match.dBegin) << "\t"
                        << triplex_error_string(match, tfo_motifs, tts_motifs, opts) << "\t"
                        << match.motif << "\t"
                        << match.strand << "\t"
                        << (match.parallel ? 'P' : 'A') << "\t"
                        << match.guanines / (match.dEnd - match.dBegin);
            if (opts.output_format == output_format_t::triplex) {
                output_file << triplex_alignment_string(match, tfo_motifs, tts_motifs);
            }
            output_file << "\n";
#if defined(_OPENMP)
        }
#endif
    }
}

void print_triplex_summary(potential_set_t& potentials,
                           name_set_t& tfo_names,
                           name_set_t& tts_names,
                           const options& opts)
{
    seqan::CharString output_file_name;
    seqan::append(output_file_name, opts.output_file);
    seqan::append(output_file_name, ".summary");

    std::ofstream output_file(seqan::toCString(output_file_name),
                              std::ios_base::out);
    if (!output_file) {
        std::cerr << "PATO: error opening output file '"
                  << seqan::toCString(output_file_name) << "'\n";
        return;
    }

    output_file << "# Duplex-ID\tSequence-ID\tTotal (abs)\tTotal (rel)\tGA (abs"
                   ")\tGA (rel)\tTC (abs)\tTC (rel)\tGT (abs)\tGT (rel)\n";
    for (auto& potential_entry : potentials) {
        auto& potential = potential_entry.second;
        if (seqan::hasCount(potential)) {
            output_file << tts_names[seqan::getKey(potential).second] << "\t"
                        << tfo_names[seqan::getKey(potential).first] << "\t"
                        << seqan::getCounts(potential) << "\t"
                        << std::setprecision(3) << seqan::getCounts(potential) / seqan::getNorm(potential) << "\t"
                        << seqan::getCount(potential, 'R') << "\t"
                        << std::setprecision(3) << seqan::getCount(potential, 'R') / seqan::getNorm(potential) << "\t"
                        << seqan::getCount(potential, 'Y') << "\t"
                        << std::setprecision(3) << seqan::getCount(potential, 'Y') / seqan::getNorm(potential) << "\t"
                        << seqan::getCount(potential, 'M') << "\t"
                        << std::setprecision(3) << seqan::getCount(potential, 'M') /seqan:: getNorm(potential) << "\t"
                        << "\n";
        }
    }
}
