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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Various useful bit-twiddling routines, mostly taken from the website
// http://www-graphics.stanford.edu/~seander/bithacks.html
// ==========================================================================

#ifndef SEQAN_MISC_MISC_BIT_TWIDDLING_H_
#define SEQAN_MISC_MISC_BIT_TWIDDLING_H_

// TODO(holtgrew): Test this!

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Classes, Structs, Enums, Tags
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setBitTo()
// ----------------------------------------------------------------------------

/**
.Function.setBitTo
..cat:Bit Twiddling
..summary:Set the bit with the given index to the given value.
..signature:setBitTo(word, index, value)
..param.word:The number.
..param.index:The index of the bit in the word.
...type:nolink:$unsigned$
..param.value:The value to set the bit to.
...type:nolink:$bool
..returns:$void$
..include:seqan/misc/misc_bit_twiddling.h
..see:Function.setBit
..see:Function.clearBit
..see:Function.clearAllBits
..see:Function.isBitSet
 */

template <typename TWord>
inline
void
setBitTo(TWord & word, unsigned index, bool value)
{
    // See http://www-graphics.stanford.edu/~seander/bithacks.html#ConditionalSetOrClearBitsWithoutBranching
    word = (word & ~(1u << index)) | (-value & (1u << index));
}

// ----------------------------------------------------------------------------
// Function setBit()
// ----------------------------------------------------------------------------

/**
.Function.setBit
..cat:Bit Twiddling
..summary:Set the bit with the given index to 1.
..signature:setBit(word, index)
..param.word:The number.
..param.index:The index of the bit in the word.
...type:nolink:$unsigned$
..returns:$void$
..include:seqan/misc/misc_bit_twiddling.h
..see:Function.setBitTo
..see:Function.clearBit
..see:Function.clearAllBits
..see:Function.isBitSet
 */

template <typename TWord>
inline
void
setBit(TWord & word, unsigned index)
{
    word |= (1u << index);
}

// ----------------------------------------------------------------------------
// Function clearBit()
// ----------------------------------------------------------------------------

/**
.Function.clearBit
..cat:Bit Twiddling
..summary:Set the bit with the given index to 0.
..signature:clearBit(word, index)
..param.word:The number.
..param.index:The index of the bit in the word.
...type:nolink:$unsigned$
..returns:$void$
..include:seqan/misc/misc_bit_twiddling.h
..see:Function.setBitTo
..see:Function.setBit
..see:Function.clearAllBits
..see:Function.isBitSet
 */

template <typename TWord>
inline
void
clearBit(TWord & word, unsigned index)
{
    word &= ~(1u << index);
}

// ----------------------------------------------------------------------------
// Function clearAllBits()
// ----------------------------------------------------------------------------

/**
.Function.clearAllBits
..cat:Bit Twiddling
..summary:Set all bits to 0.
..signature:clearAllBits(word)
..param.word:The number.
..returns:$void$
..include:seqan/misc/misc_bit_twiddling.h
..see:Function.setBitTo
..see:Function.setBit
..see:Function.clearBit
..see:Function.isBitSet
 */

template <typename TWord>
inline
void
clearBits(TWord & word)
{
    word = 0;
}

// ----------------------------------------------------------------------------
// Function isBitSet()
// ----------------------------------------------------------------------------

/**
.Function.isBitSet
..cat:Bit Twiddling
..summary:Returns whether the bit with the given index is set to 1.
..signature:isBitSet(word, index)
..param.word:The number.
..param.index:The index.
...type:nolink:$unsigned$
..returns:$true$ if the bit was set and $false$ otherwise.
...type:nolink:$bool$
..include:seqan/misc/misc_bit_twiddling.h
..see:Function.setBitTo
..see:Function.setBit
..see:Function.clearBit
..see:Function.clearAllBits
 */

template <typename TWord>
inline
bool
isBitSet(TWord const & word, unsigned index)
{
    return (word & (1u << index)) != static_cast<TWord>(0);
}

}  // namespace seqan

#endif // #ifndef SEQAN_MISC_MISC_BIT_TWIDDLING_H_
