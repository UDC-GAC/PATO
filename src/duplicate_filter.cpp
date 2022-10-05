/*
 * MIT License
 *
 * Copyright (c) 2022 Iñaki Amatria-Barral
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "duplicate_filter.hpp"

#include <algorithm>

#include <seqan/map.h>
#include <seqan/index.h>
#include <seqan/sequence.h>

#include "triplex_enums.hpp"
#include "triplex_pattern.hpp"

void count_duplicates(motif_set_t& motifs, const options& opts)
{
    triplex_string_set_t strings;
    if (opts.run_mode == run_mode_t::tts_search
        && opts.detect_duplicates == detect_duplicates_t::permissive) {
        for (auto& motif : motifs) {
            if (seqan::getMotif(motif) == '+') {
                seqan::appendValue(strings, seqan::getSegment(motif));
            } else {
                seqan::appendValue(strings, mod_rev_t(seqan::getSegment(motif)));
            }
        }
    } else {
        for (auto& motif : motifs) {
            seqan::appendValue(strings, opts.detect_duplicates == detect_duplicates_t::permissive
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
        if (opts.run_mode == run_mode_t::tts_search
            && opts.detect_duplicates == detect_duplicates_t::permissive) {
            if (seqan::getMotif(motif) == '+') {
                pattern = seqan::getSegment(motif);
            } else {
                pattern = mod_rev_t(seqan::getSegment(motif));
            }
        } else {
            pattern = opts.detect_duplicates == detect_duplicates_t::permissive
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
            && opts.run_mode != run_mode_t::tpx_search
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
