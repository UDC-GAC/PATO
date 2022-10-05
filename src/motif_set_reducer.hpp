#ifndef MOTIF_SET_REDUCER_HPP
#define MOTIF_SET_REDUCER_HPP

#include "triplex_definitions.hpp"

template<typename pair_t>
struct second_t
{
    typename pair_t::second_type operator()(const pair_t& p) const
    {
        return p.second;
    }
};

template<typename map_t>
second_t<typename map_t::value_type> second(const map_t& m)
{
    return second_t<typename map_t::value_type>();
}

void reduce_motif_set(motif_set_t& output, motif_set_t& input);

#endif
