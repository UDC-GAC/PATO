#ifndef SEGMENT_PARSER_HPP
#define SEGMENT_PARSER_HPP

#include "triplex_definitions.hpp"

void make_parser(graph_t& parser,
                 triplex_t& valid_chars,
                 triplex_t& invalid_chars,
                 unsigned int max_interrupts);
void parse_segments(graph_t& parser,
                    segment_set_t& segments,
                    triplex_t& sequence,
                    unsigned int max_interrupts,
                    int min_length);

#endif
