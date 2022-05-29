// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Basic comparison code.
// ==========================================================================

// TODO(holtgrew): Rename to basic_comparison.h?
// TODO(holtgrew): Necessary? Documentation is somewhere else..

#ifndef SEQAN_BASIC_BASIC_COMPARE_H_
#define SEQAN_BASIC_BASIC_COMPARE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TLeft, typename TRight>
struct CompareType;

template <typename TLeft, typename TRight>
struct CompareType<TLeft const, TRight>
{
    typedef typename CompareType<TLeft, TRight>::Type const Type;
};

template <typename TLeft, typename TRight>
struct CompareType<TLeft, TRight const>
{
    typedef typename CompareType<TLeft, TRight>::Type const Type;
};

template <typename TLeft, typename TRight>
struct CompareType<TLeft const, TRight const>
{
    typedef typename CompareType<TLeft, TRight>::Type const Type;
};

// ============================================================================
// Functions
// ============================================================================

// TODO(holtgrew): What does this do / extend / complement?

template<typename T_> inline
bool lexLess(const T_& _Left, const T_& Right_)
{
    typedef typename MakeUnsigned_<T_>::Type TUnsigned;
    return (TUnsigned)_Left < (TUnsigned)Right_;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_BASIC_COMPARE_H_
