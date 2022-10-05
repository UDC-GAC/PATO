#ifndef REPEAT_FILTER_HPP
#define REPEAT_FILTER_HPP

#include "triplex_definitions.hpp"

void filter_repeats(repeat_set_t& repeats,
                    triplex_t& sequence,
                    unsigned int min_repeat_length,
                    unsigned int max_repeat_period);

#endif
