// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_AGGREGATE_CONCEPT_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_AGGREGATE_CONCEPT_H_

namespace seqan {

// TODO(holtgrew): What about empty base class optimization as does Boost's compressed_pair? Not useful?

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Concept.Aggregate
..summary:Aggregate types contain a fixed number of fixed-size values.
..remarks:Stream output operators are not shown in the function list below, but required.
..remarks:Comparison operators are not shown in the function list below, but required.

..Function.clear.concept:Concept.Aggregate
..Function.value.concept:Concept.Aggregate
..Function.assignValue.concept:Concept.Aggregate

..Metafunction.LENGTH.concept:Concept.Aggregate
..Metafunction.Value.concept:Concept.Aggregate
 */

/**
.Tag.Compressed
..cat:Aggregates
..summary:Tag to marke a "compressed" specialization.
..signature:Compressed
..include:seqan/basic.h
 */

struct Compressed_;
typedef Tag<Compressed_> Compressed;

/**
.Tag.BitCompressed
..cat:Aggregates
..summary:Tag to marke a "compressed" specialization.
..signature:Compressed<BITSIZE1, BITSIZE2>
..param.BITSIZE1:Number of bits used for first element.
...type:nolink:$unsigned$
..param.BITSIZE2:Number of bits used for second element.
...type:nolink:$unsigned$
..include:seqan/basic.h
 */

template <unsigned BITSIZE1 = 16, unsigned BITSIZE2 = 16>
struct BitCompressed;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_AGGREGATE_CONCEPT_H_
