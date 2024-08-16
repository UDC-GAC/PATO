// MIT License
//
// Copyright (c) 2022-onwards Iñaki Amatria-Barral
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef PATO_LIB_TYPES_H
#define PATO_LIB_TYPES_H

#include <seqan/graph_types.h>
#include <seqan/index.h>

#include <Triplexator/alphabet.h>
#include <Triplexator/match.h>
#include <Triplexator/pattern.h>

namespace pato {

using triplex_t = seqan::TriplexString;
using motif_t = seqan::ModStringTriplex<triplex_t, triplex_t>;
using motif_potential_t = seqan::TriplexPotential<unsigned>;
using repeat_t = seqan::Repeat<unsigned, unsigned>;
using segment_t = seqan::Infix<triplex_t>::Type;
using graph_t = seqan::Graph<seqan::Automaton<seqan::Triplex, seqan::Triplex>>;
using vertex_descriptor_t = seqan::VertexDescriptor<graph_t>::Type;
using filter_t =
    seqan::ModifiedString<motif_t, seqan::ModView<seqan::FunctorRYFilter>>;
using match_t =
    seqan::TriplexMatch<seqan::Difference<seqan::TriplexString>::Type, unsigned,
                        unsigned>;
using potential_t = seqan::TriplexPotential<std::pair<unsigned, unsigned>>;

using motif_vector_t = std::vector<motif_t>;
using motif_potential_vector_t = std::vector<motif_potential_t>;
using name_vector_t = std::vector<seqan::CharString>;
using triplex_vector_t = std::vector<triplex_t>;
using index_vector_t = std::vector<unsigned>;
using repeat_vector_t = std::vector<repeat_t>;
using char_vector_t = std::vector<char>;
using char_vector_vector_t = std::vector<char_vector_t>;
using segment_vector_t = std::vector<segment_t>;
using match_vector_t = std::vector<match_t>;
#if defined(_OPENMP)
using match_vector_vector_t = std::vector<match_vector_t>;
#endif

struct pair_hash_t {
  template <class Ty, class Tx>
  std::size_t operator()(const std::pair<Ty, Tx> &pair) const {
    return std::hash<Ty>()(pair.first) ^ std::hash<Tx>()(pair.second);
  }
};

using potential_map_t =
    std::unordered_map<std::pair<unsigned, unsigned>, potential_t, pair_hash_t>;

} // namespace pato

#endif // PATO_LIB_TYPES_H
