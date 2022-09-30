#include "duplicate_filter.hpp"

#include <algorithm>

#include <seqan/map.h>
#include <seqan/index.h>
#include <seqan/sequence.h>

#include "triplex_pattern.hpp"

void count_duplicates(motif_set_t& motifs, const options& opts)
{
    triplex_string_set_t strings;
    if (opts.run_mode == 1 && opts.detect_duplicates == duplicate::permissive) {
        for (auto& motif : motifs) {
            if (seqan::getMotif(motif) == '+') {
                seqan::appendValue(strings, seqan::getSegment(motif));
            } else {
                seqan::appendValue(strings, mod_rev_t(seqan::getSegment(motif)));
            }
        }
    } else {
        for (auto& motif : motifs) {
            seqan::appendValue(strings, opts.detect_duplicates == duplicate::permissive
                                        ? seqan::tfoString(motif)
                                        : seqan::ttsString(motif));
        }
    }
    index_t index(strings);

    for (auto& motif : motifs) {
        duplicate_t identity(seqan::getSequenceNo(motif),
                             seqan::beginPosition(motif));
        duplicate_map_t duplicates;

        int copies = 0;
        finder_t finder(index);

        triplex_t pattern;
        if (opts.run_mode == 1 && opts.detect_duplicates == duplicate::permissive) {
            if (seqan::getMotif(motif) == '+') {
                pattern = seqan::getSegment(motif);
            } else {
                pattern = mod_rev_t(seqan::getSegment(motif));
            }
        } else {
            pattern = opts.detect_duplicates == duplicate::permissive
                      ? seqan::tfoString(motif)
                      : seqan::ttsString(motif);
        }

        while (seqan::find(finder, pattern)) {
            auto& tmp_motif = motifs[seqan::position(finder).i1];

            if (!opts.same_sequence_duplicates
                && seqan::getSequenceNo(motif) == seqan::getSequenceNo(tmp_motif)) {
                continue;
            }

            duplicate_pos_t pos_value = seqan::getMotif(tmp_motif) == '-'
                                        ? seqan::endPosition(tmp_motif) - seqan::position(finder).i2 - seqan::length(motif)
                                        : seqan::beginPosition(tmp_motif) + seqan::position(finder).i2;
            duplicate_t pos(seqan::getSequenceNo(tmp_motif), pos_value);
            if (pos != identity && !seqan::hasKey(duplicates, pos)) {
                copies++;
                seqan::add(duplicates, pos);
            }

            if (opts.duplicate_cutoff >= 0 && copies > opts.duplicate_cutoff) {
                break;
            }
        }
        seqan::duplicates(motif, copies);

        if (opts.report_duplicate_locations
            && opts.duplicate_cutoff >= 0
            && copies <= opts.duplicate_cutoff) {
            for (auto& duplicate : duplicates) {
                seqan::addDuplicate(motif,
                                    seqan::getSequenceNo(duplicate),
                                    seqan::getPosition(duplicate));
            }
        }
    }
}

void filter_duplicates(motif_set_t& motifs, int cutoff)
{
    motifs.erase(std::remove_if(motifs.begin(), motifs.end(), [&](const auto& motif) -> bool {
        return seqan::duplicates(motif) > cutoff;
    }), motifs.end());
}
