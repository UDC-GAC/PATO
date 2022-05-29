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
// Base class for streams.
// ==========================================================================

#ifndef SEQAN_STREAM_STREAM_BASE_H_
#define SEQAN_STREAM_STREAM_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Class.Stream
..cat:Input/Output
..signature:Stream<TSpec>
..summary:Abstract base class to fulfill the @Concept.Stream@ concept.
..concept:Concept.Stream
..include:seqan/stream.h
 */

template <typename TPointer = char *>
struct CharArray;

#if SEQAN_HAS_ZLIB  // Enable Stream<GZFile> if available.
struct GZFile_;
typedef Tag<GZFile_> GZFile;
#endif  // #if SEQAN_HAS_ZLIB

#if SEQAN_HAS_BZIP2  // Enable Stream<BZ2File> if available.
struct BZ2File_;
typedef Tag<BZ2File_> BZ2File;
#endif  // #if SEQAN_HAS_ZLIB

template <typename TSpec>
class Stream;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

///.Function.atEnd.iterator.type:Class.Stream

template <typename TSpec>
inline bool
atEnd(Stream<TSpec> & stream)
{
    return streamEof(stream);
}

template <typename TSpec>
inline bool
atEnd(Stream<TSpec> const & stream)
{
    return streamEof(stream);
}

// ----------------------------------------------------------------------------
// Function streamPut()
// ----------------------------------------------------------------------------

// this is generic for all specializations of Stream<> right now

template <typename TStream>
inline int
streamPut(Stream<TStream> & stream, char const c)
{
    return streamWriteChar(stream, c);
}

template <typename TStream>
inline int
streamPut(Stream<TStream> & stream, char const * source)
{
    return (streamWriteBlock(stream, source, strlen(source))
                == strlen(source) )  ?   0 : 1;
}

template <typename TStream, typename TSpec>
inline int
streamPut(Stream<TStream> & stream, String<char, TSpec> const & source)
{
    return (streamWriteBlock(stream, toCString(source), length(source))
                == length(source))  ?   0 : 1;
}

template <typename TStream, typename TSource>
inline int
streamPut(Stream<TStream> & stream, TSource const & source)
{
    char buffer[1024] = "";
    ::std::stringstream s;

    s << source;
    if (s.fail())
        return s.fail();

    s >> buffer;
    if (s.fail())
        return s.fail();

    buffer[1023] = 0;

//TODO(h4nn3s): we should be able to use the following and then s.str() directly
// so we wouldnt need an extra buffer at all. but it doesnt work
//     s << source << std::ends;
//     if (s.fail())
//         return s.fail();

    return (streamWriteBlock(stream, buffer, strlen(buffer))
                == strlen(buffer) )  ?   0 : 1;
}

}  // namespace seqean

#endif  // #ifndef SEQAN_STREAM_STREAM_BASE_H_
