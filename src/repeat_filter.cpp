#include "repeat_filter.hpp"

#include <seqan/index.h>

void filter_repeats(repeat_set_t& repeats,
                    triplex_t& sequence,
                    unsigned int min_repeat_length,
                    unsigned int max_repeat_period)
{
    seqan::findRepeats(repeats, sequence, min_repeat_length, max_repeat_period);

    for (auto& repeat : repeats) {
        for (unsigned int i = repeat.beginPosition; i < repeat.endPosition; i++) {
            sequence[i] = 'N';
        }
    }
}
