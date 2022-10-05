#ifndef TFO_FINDER_HPP
#define TFO_FINDER_HPP

#include "options.hpp"
#include "triplex_definitions.hpp"

bool find_tfo_motifs(motif_set_t& motifs,
                     motif_potential_set_t& potentials,
                     triplex_set_t& sequences,
                     name_set_t& names,
                     const options& opts);
void find_tfo_motifs(const options& opts);

#endif
