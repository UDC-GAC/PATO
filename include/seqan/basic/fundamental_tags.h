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
// Author: Andres Gogol-Döring <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Global (future: generic) tag definitions.
// ==========================================================================

// TODO(holtgrew): This should probably be minimalized, tags should be moved to single modules whenever possible.

#ifndef SEQAN_BASIC_BASIC_TAGS_H_
#define SEQAN_BASIC_BASIC_TAGS_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Tag<T>
// ----------------------------------------------------------------------------

/**
.Tag.Tag
..cat:Fundamental
..signature:Tag<T>
..param.T:Any parameter less types.
..summary:Template for tag definition.
..remarks:
This $struct$ is defined such that parameter less tags are easier recognizeable.
This is best explained with the example below.
..example.text:Usually, tags are defined in the following way.
..example.code:
struct SomeTag_;
typedef Tag<SomeTag_> SomeTag;
..example.text:They are then used as follows.
..example.code:
template <typename T>
void f(T const & x, SomeTag const & tag)
{
    // ...
}

// Somewhere else:
f(3, SomeTag());
..example.text:
This has the advantages that (1) the type of tag parameters is printed as $Tag<SomeTag_>$ in compiler error traces.
Furthermore, (2) parameter less tags can be defined redundantly in multiple headers and we can still instantiate them anywhere where they are declared.
The latter (2) cannot be achieved with only forward declaration ($struct SomeTag;$) or full declarations ($struct SomeTag {};$) everywhere.
..include:seqan/basic.h
 */

template <typename T>
struct Tag {};

// ----------------------------------------------------------------------------
// Tag Default
// ----------------------------------------------------------------------------

/**
.Tag.Default:
..cat:Basic
..summary:Tag that specifies default behavior.
..tag.Default:Use default behavior. 
..include:seqan/basic.h
*/
struct Default_;
typedef Tag<Default_> Default;

// ----------------------------------------------------------------------------
// Tag Nothing
// ----------------------------------------------------------------------------

/**
.Tag.Nothing:
..cat:Basic
..summary:Tag that represents an absent parameter or an absent type.
..tag.Nothing:Omit parameter.
..include:seqan/basic.h
*/

struct Nothing {};

// ----------------------------------------------------------------------------
// Tag Move
// ----------------------------------------------------------------------------

// TODO(holtgrew): This should probably go into basic_transport.h.

/**
.Tag.Move Switch:
..cat:Basic
..summary:Switch to force move.
..tag.Move:Move instead of assign. 
..remarks.text:The difference between move constructor and copy constructor
is that the source object is not copied but moved into the target object.
The source object can lose its content and will be empty after
this operation in this case.
A move constructor can sigificantly faster than a copy constructor.
..example.code:String source("hello");
String target(source, Move()); // source is moved to target
std::cout << source; //nothing printed since source lost content
std::cout << target; //"hello"
..see:Function.move
.example.text:Move constructors are like copy-constructors. However, their argument is not const.
.example.code:
class Klass
{
public:
    seqan::String m;

    // Copy constructor, other is left untouched.
    Klass(Klass const & other) { ... }

    // Move constructor, leaves other and its members in an "empty" state.
    Klass(Klass & other, seqan::Move const &) { ... }
};
..include:seqan/basic.h
*/

struct Move_;
typedef Tag<Move_> Move;

// ----------------------------------------------------------------------------
// Tag MinimalCtor
// ----------------------------------------------------------------------------

// TODO(holtgrew): This should go into initialization part of alphabet header set.

//construct without initializing
struct MinimalCtor_;
typedef Tag<MinimalCtor_> MinimalCtor;

// ----------------------------------------------------------------------------
// Tag MinimalCtor
// ----------------------------------------------------------------------------

// TODO(holtgrew): This should go into initialization part of alphabet header set.

//construct with initializing
struct NonMinimalCtor_;
typedef Tag<NonMinimalCtor_> NonMinimalCtor;

// ----------------------------------------------------------------------------
// Tag MinimalCtor
// ----------------------------------------------------------------------------

// TODO(holtgrew): This should go into the iterators header set?

//Pass to c'tor of iterator to move it to the end
struct GoEnd_;
typedef Tag<GoEnd_> GoEnd;

// ----------------------------------------------------------------------------
// Tag TagList<TTag, TNext>
// ----------------------------------------------------------------------------

/**
.Tag.TagList:
..cat:Basic
..summary:A structure to represent a list of tags.
..signature:TagList<TTag1>
..signature:TagList<TTag1, TagList<TTag2> >
..signature:TagList<TTag1, TagList<TTag2, TagList<TTag3[...]> > >
..param.TTag1:The first tag of the list.
..param.TTag2:The second tag of the list.
..param.TTag3:The third tag of the list.
..include:seqan/basic.h
*/

template <typename TTag = void, typename TSubList = void>
struct TagList
{
    typedef TTag Type;
};

// ----------------------------------------------------------------------------
// Class TagSelector
// ----------------------------------------------------------------------------

/**
.Class.TagSelector:
..cat:Basic
..summary:A structure to select a tag from a @Tag.TagList@.
..signature:TagSelector<TTagList>
..param.TTagList:A tag list.
...type:Tag.TagList
.Memvar.TagSelector#tagId:
..class:Class.TagSelector
..type:nolink:int
..cat:Basic
..summary:Stores the index of a @Page.Glossary.Tag@ in the tag list.
..include:seqan/basic.h
*/

template <typename TTagList = void>
struct TagSelector
{
    int tagId;

    TagSelector()
    {
        tagId = 0;
    }

    inline bool
    operator==(TagSelector const & other) const
    {
        return other.tagId == tagId;
    }
};

template <typename TTag, typename TSubList>
struct TagSelector< TagList<TTag, TSubList> >
        : TagSelector<TSubList>
{
    typedef TTag                    Type;
    typedef TagSelector<TSubList>   Base;
};

// ----------------------------------------------------------------------------
// Tag DotDrawing
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should probably not be defined here.

/**
.Tag.DotDrawing
..cat:Input/Output
..summary:Switch to trigger drawing in dot format.
..value.DotDrawing:Graphs in dot format.
..include:seqan/basic.h
*/

struct DotDrawing_;
typedef Tag<DotDrawing_> DotDrawing;

// TODO(holtgrew): Should probably not be defined here.
// TODO(holtgrew): Are these used at all?

/**
.Tag.HammingDistance
..cat:Basic
..summary:Switch to trigger Hamming distance, which is a measure of character substitutions.
..include:seqan/basic.h
*/

/**
.Tag.LevenshteinDistance
..cat:Basic
..summary:Switch to trigger Levenshtein distance, which is a measure of edit operations (character substitutions, deletions or insertions).
..remarks:$EditDistance$ is a synonym for $LevenshteinDistance$.
..see:Spec.EditDistance
..include:seqan/basic.h
*/

struct HammingDistance_;
struct LevenshteinDistance_;

typedef Tag<HammingDistance_>       HammingDistance;
typedef Tag<LevenshteinDistance_>   LevenshteinDistance;
typedef Tag<LevenshteinDistance_>   EditDistance; 

// ----------------------------------------------------------------------------
// Tag NeedlemanWunsch
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should probably not be defined here.

//Sollte eigentlich nach align/, aber da jetzt ja so viele
//alignment algorithmen in graph/ gelandet sind...

/**
.Tag.Global Alignment Algorithms:
..cat:Alignments
..summary:Global alignment algorithm used by globalAlignment.
..see:Function.globalAlignment
..see:Tag.Local Alignment Algorithms
..include:seqan/basic.h
*/

/**
.Tag.Global Alignment Algorithms.value.NeedlemanWunsch:
    Dynamic programming algorithm for alignments by Needleman and Wunsch.
..include:seqan/basic.h
*/

struct NeedlemanWunsch_;
typedef Tag<NeedlemanWunsch_> NeedlemanWunsch;

// ----------------------------------------------------------------------------
// Tag BandedNeedlemanWunsch
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should probably not be defined here.

/**
.Tag.Global Alignment Algorithms.value.BandedNeedlemanWunsch:
    The Needleman-Wunsch alignment algorithm in a banded version.
..include:seqan/basic.h
*/
struct BandedNeedlemanWunsch_;
typedef Tag<BandedNeedlemanWunsch_> BandedNeedlemanWunsch;

// ----------------------------------------------------------------------------
// Tag Gotoh
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should probably not be defined here.

/**
.Tag.Global Alignment Algorithms.value.Gotoh:
    Gotoh's affine gap cost alignment algorithm.
..include:seqan/basic.h
*/
struct Gotoh_;
typedef Tag<Gotoh_> Gotoh;

// ----------------------------------------------------------------------------
// Tag BandedGotoh
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should probably not be defined here.

/**
.Tag.Global Alignment Algorithms.value.BandedGotoh:
    Gotoh's affine gap cost alignment algorithm in a banded version.
..include:seqan/basic.h
*/
struct BandedGotoh_;
typedef Tag<BandedGotoh_> BandedGotoh;

// ----------------------------------------------------------------------------
// Tag MyersBitVector
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should probably not be defined here.

/**
.Tag.Global Alignment Algorithms.value.MyersBitVector:
    Myers' bit vector alignment algorithm for edit distance.
    Note that this algorithm does not returns the alignment itself, but only computes the score.
..include:seqan/basic.h
*/
struct MyersBitVector_;
typedef Tag<MyersBitVector_> MyersBitVector;

// ----------------------------------------------------------------------------
// Tag MyersHirschberg
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should probably not be defined here.

/**
.Tag.Global Alignment Algorithms.value.MyersHirschberg:
    Myers' bit vector algorithm for edit distance combined with Hirschberg's linear space alignment algorithm.
..include:seqan/basic.h
*/
struct MyersHirschberg_;
typedef Tag<MyersHirschberg_> MyersHirschberg;

// TODO(holtgrew): Should probably not be defined here.

// ----------------------------------------------------------------------------
// Tag Hirschberg
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should probably not be defined here.

/**
.Tag.Global Alignment Algorithms.value.Hirschberg:
    Hirschberg's linear space global alignment algorithm.
..include:seqan/basic.h
*/
struct Hirschberg_;
typedef Tag<Hirschberg_> Hirschberg;

// ----------------------------------------------------------------------------
// Tag Lcs
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should probably not be defined here.

/**
.Tag.Global Alignment Algorithms.value.Lcs:
    Longest common subsequence algorithm.
..include:seqan/basic.h
*/
struct Lcs_;
typedef Tag<Lcs_> Lcs;

// ----------------------------------------------------------------------------
// Tag SmithWaterman
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should probably not be defined here.

/**
.Tag.Local Alignment Algorithms:
..cat:Alignments
..summary:Local alignment algorithm used by localAlignment.
..see:Function.localAlignment
..include:seqan/basic.h
*/

/**
.Tag.Local Alignment Algorithms.value.SmithWaterman:
    Triggers a Smith Waterman local alignment algorithm.
..include:seqan/basic.h
*/
struct SmithWaterman_;
typedef Tag<SmithWaterman_> SmithWaterman;

// ----------------------------------------------------------------------------
// Tag BandedSmithWaterman
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should probably not be defined here.

/**
.Tag.Local Alignment Algorithms.value.BandedSmithWaterman:
    Triggers a banded version of the Smith Waterman local alignment algorithm.
..include:seqan/basic.h
*/
struct BandedSmithWaterman_;
typedef Tag<BandedSmithWaterman_> BandedSmithWaterman;

// ----------------------------------------------------------------------------
// Tags SmithWatermanClump, WatermanEggert
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should probably not be defined here.

/**
.Tag.Local Alignment Algorithms.value.WatermanEggert:
    Local alignment algorithm by Waterman and Eggert with "declumping" (i.e. only non-overlapping local alignments are computed).
.Tag.Local Alignment Algorithms.value.SmithWatermanClump:
    Same as $WatermanEggert$.
..include:seqan/basic.h
*/
struct SmithWatermanClump_;
typedef Tag<SmithWatermanClump_> SmithWatermanClump;
typedef Tag<SmithWatermanClump_> WatermanEggert;

// ----------------------------------------------------------------------------
// Tags BandedWatermanEggert, BandedSmithWatermanClump
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should probably not be defined here.

/**
.Tag.Local Alignment Algorithms.value.BandedWatermanEggert:
    Triggers a banded version of the local alignment algorithm by Waterman and Eggert with "declumping".
.Tag.Local Alignment Algorithms.value.BandedSmithWatermanClump:
    Same as $BandedWatermanEggert$.
..include:seqan/basic.h
*/
struct BandedWatermanEggert_;
typedef Tag<BandedWatermanEggert_> BandedSmithWatermanClump;
typedef Tag<BandedWatermanEggert_> BandedWatermanEggert;

// ----------------------------------------------------------------------------
// Tag Blat
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should probably not be defined here.

struct Blat_;
typedef Tag<Blat_> Blat;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction LENGTH;  For TagList.
// ----------------------------------------------------------------------------

// TODO(holtgrew): Is this defined here or is it just a forward?

///.Metafunction.LENGTH.param.T.type:Tag.TagList

template <>
struct LENGTH<void>
{
    enum { VALUE = 0 };
};

template <typename TTag>
struct LENGTH<TagList<TTag, void> >
{
    enum { VALUE = 1 };
};

template <typename TTag, typename TSubList>
struct LENGTH<TagList<TTag, TSubList> >
{
    enum { VALUE = LENGTH<TSubList>::VALUE + 1 };
};

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_BASIC_TAGS_H_
