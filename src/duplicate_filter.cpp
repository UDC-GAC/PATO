#include "duplicate_filter.hpp"

#include <algorithm>

#include <seqan/map.h>
#include <seqan/index.h>

void count_duplicates(motif_set_t& motifs, const options& opts)
{
    triplex_string_set_t strings;
    for (auto& motif : motifs) {
        if (opts.detect_duplicates == duplicate::permissive) {
            seqan::appendValue(strings, seqan::tfoString(motif));
        } else {
            seqan::appendValue(strings, seqan::ttsString(motif));
        }
    }
    index_t index(strings);

    for (auto& motif : motifs) {
        duplicate_t identity(seqan::getSequenceNo(motif),
                             seqan::beginPosition(motif));
        duplicate_map_t duplicates;

        int copies = 0;
        finder_t finder(index);
        while (seqan::find(finder,
                           opts.detect_duplicates == duplicate::permissive
                           ? seqan::tfoString(motif)
                           : seqan::ttsString(motif))) {
            auto& tmp_motif = motifs[seqan::position(finder).i1];

            if (!opts.same_sequence_duplicates
                && seqan::getSequenceNo(motif) == seqan::getSequenceNo(tmp_motif)) {
                continue;
            }

            duplicate_t pos(seqan::getSequenceNo(tmp_motif),
                            seqan::beginPosition(tmp_motif) + seqan::position(finder).i2);
            if (pos != identity && !seqan::hasKey(duplicates, pos)) {
                copies++;
                seqan::add(duplicates, pos);
            }

            if (opts.duplicate_cutoff >= 0 && copies > opts.duplicate_cutoff) {
                break;
            }
        }

        seqan::duplicates(motif, copies);
    }
}

void filter_duplicates(motif_set_t& motifs, int cutoff)
{
    motifs.erase(std::remove_if(motifs.begin(), motifs.end(), [&](const auto& motif) -> bool {
        return seqan::duplicates(motif) > cutoff;
    }), motifs.end());
}
