#ifndef DUPLICATE_FILTER_HPP
#define DUPLICATE_FILTER_HPP

#include "options.hpp"
#include "triplex_definitions.hpp"

void count_duplicates(motif_set_t& motifs, const options& opts);
void filter_duplicates(motif_set_t& motifs, int cutoff);

#endif
