#ifndef OUTPUT_WRITER_HPP
#define OUTPUT_WRITER_HPP

#include "options.hpp"
#include "triplex_definitions.hpp"

void print_tfo_motifs(motif_set_t& tfo_motifs,
                      name_set_t& tfo_names,
                      const options& opts);
void print_tfo_potentials(motif_potential_set_t& tfo_potentials,
                          name_set_t& tfo_names,
                          const options& opts);

#if !defined(_OPENMP)
void print_triplex_pairs(match_set_t& matches,
                         motif_set_t& tfo_motifs,
                         name_set_t& tfo_names,
                         motif_set_t& tts_motifs,
                         name_set_t& tts_names,
                         const options& opts);
#else
void print_triplex_pairs(match_set_set_t& matches,
                         motif_set_t& tfo_motifs,
                         name_set_t& tfo_names,
                         motif_set_t& tts_motifs,
                         name_set_t& tts_names,
                         const options& opts);
#endif
void print_triplex_summary(potential_set_t& potentials,
                           name_set_t& tfo_names,
                           name_set_t& tts_names,
                           const options& opts);

#endif
