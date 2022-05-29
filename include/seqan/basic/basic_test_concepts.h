// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Test for fulfilled basic concepts.
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS

#ifndef CORE_INCLUDE_SEQAN_BASIC_TEST_CONCEPTS_H_
#define CORE_INCLUDE_SEQAN_BASIC_TEST_CONCEPTS_H_

namespace seqan {

// ============================================================================
// Test basic concepts
// ============================================================================

inline void testAlphabetConcepts()
{
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<bool>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<char>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<unsigned>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<int>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<double>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<Pair<int, double> >));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<Dna>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<Dna5>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<AminoAcid>));
    SEQAN_CONCEPT_ASSERT((AlphabetConcept<Ascii>));
}

template <typename T>
inline void testInteger()
{
    SEQAN_CONCEPT_ASSERT((IntegerConcept<T>));
    SEQAN_STATIC_ASSERT_MSG(Is< IntegerConcept<T> >::VALUE, "Type is not marked to be an integer");
}

template <typename T>
inline void testSignedInteger()
{
    SEQAN_CONCEPT_ASSERT((SignedIntegerConcept<T>));
    SEQAN_STATIC_ASSERT_MSG(Is< SignedIntegerConcept<T> >::VALUE, "Type is not marked to be a signed integer");
    testInteger<T>();
}

template <typename T>
inline void testUnsignedInteger()
{
    SEQAN_CONCEPT_ASSERT((UnsignedIntegerConcept<T>));
    SEQAN_STATIC_ASSERT_MSG(Is< UnsignedIntegerConcept<T> >::VALUE, "Type is not marked to be an unsigned integer");
    testInteger<T>();
}

inline void testIntegers()
{
    testInteger<char>();
    testSignedInteger<signed char>();
    testSignedInteger<signed short>();
    testSignedInteger<signed int>();
    testSignedInteger<signed long>();
    testSignedInteger<__int64>();
    testUnsignedInteger<unsigned char>();
    testUnsignedInteger<unsigned short>();
    testUnsignedInteger<unsigned int>();
    testUnsignedInteger<unsigned long>();
    testUnsignedInteger<__uint64>();
//  testSignedInteger<unsigned long>();   // <== this should fail
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BASIC_TEST_CONCEPTS_H_
