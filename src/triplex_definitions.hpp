#ifndef TRIPLEX_DEFINITIONS_HPP
#define TRIPLEX_DEFINITIONS_HPP

#include <vector>
#include <unordered_map>

#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/graph_types.h>

#include "triplex_match.hpp"
#include "triplex_pattern.hpp"
#include "triplex_alphabet.hpp"
#include "triplex_functors.hpp"

enum orientation
{
    antiparallel = -1,
    both = 0,
    parallel = 1
};

struct pair_hash_t
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

struct _tts;
struct _tfo;
struct _mixed_motif;
struct _purine_motif;
struct _pyrimidine_motif;

typedef seqan::Tag<_tts> tts_t;
typedef seqan::Tag<_tfo> tfo_t;
typedef seqan::Tag<_mixed_motif> mixed_motif_t;
typedef seqan::Tag<_purine_motif> purine_motif_t;
typedef seqan::Tag<_pyrimidine_motif> pyrimidine_motif_t;

typedef seqan::TriplexString triplex_t;
typedef seqan::ModStringTriplex<triplex_t, triplex_t> motif_t;
typedef seqan::Repeat<unsigned int, unsigned int> repeat_t;
typedef typename seqan::Infix<triplex_t>::Type segment_t;
typedef seqan::ModifiedString<motif_t, seqan::ModView<seqan::FunctorRYFilter>> filter_t;
typedef seqan::TriplexMatch<seqan::Difference<seqan::TriplexString>::Type, unsigned int, double> match_t;
typedef seqan::TriplexPotential<std::pair<unsigned int, unsigned int>> potential_t;

typedef std::vector<triplex_t> triplex_set_t;
typedef std::vector<seqan::CharString> name_set_t;
typedef std::vector<motif_t> motif_set_t;
typedef std::vector<repeat_t> repeat_set_t;
typedef std::vector<segment_t> segment_set_t;
typedef std::vector<char> char_set_t;
typedef std::vector<char_set_t> char_set_set_t;
typedef std::vector<match_t> match_set_t;
typedef std::vector<unsigned int> index_set_t;
typedef std::unordered_map<std::pair<unsigned int, unsigned int>, potential_t, pair_hash_t> potential_set_t;

#if defined(_OPENMP)
typedef std::vector<match_set_t> match_set_set_t;
#endif

typedef seqan::Graph<seqan::Automaton<seqan::Triplex, seqan::Triplex>> graph_t;
typedef seqan::VertexDescriptor<graph_t>::Type vertex_descriptor_t;

#endif
