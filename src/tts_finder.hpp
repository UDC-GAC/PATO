#ifndef TTS_FINDER_HPP
#define TTS_FINDER_HPP

#include "options.hpp"
#include "triplex_definitions.hpp"

bool find_tts_motifs(motif_set_t& motifs,
                     triplex_set_t& sequences,
                     name_set_t& names,
                     const options& opts);

#endif
