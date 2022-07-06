#ifndef TFO_FINDER_HPP
#define TFO_FINDER_HPP

#include "options.hpp"
#include "triplex_definitions.hpp"

bool find_tfo_motifs(motif_set_t& motifs,
                     triplex_set_t& sequences,
                     name_set_t& names,
                     const options& opts);

#endif
