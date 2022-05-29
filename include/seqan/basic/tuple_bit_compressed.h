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
// Bit-compressed tuple specialization.
// ==========================================================================

// TODO(holtgrew): Should this be called Packed and the tag be BitPacked?

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_TUPLE_BIT_COMPRESSED_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_TUPLE_BIT_COMPRESSED_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Spec.Bit Compressed Tuple:
..cat:Aggregates
..general:Class.Tuple
..summary:A plain fixed-length string. Saves memory by packing bits.
..signature:Tuple<T, SIZE, Compressed>
..param.T:The value type, that is the type of characters stored in the tuple.
..param.SIZE:The size/length of the tuple.
...remarks:In contrast to @Class.String@ the length of Tuple is fixed.
..notes:The characters are stored as a bit sequence in an ordinal type (char, ..., __int64).
..remarks:Only useful for small alphabets and small tuple sizes (|Sigma|^size <= 2^64) as for @Spec.Dna@ or @Spec.AminoAcid@ m-grams)
..see:Spec.Sampler
..include:seqan/basic.h
 */

template <unsigned char _size>
struct BitVector_
{
    typedef typename BitVector_<_size + 1>::Type Type;
};

template <> struct BitVector_<8> { typedef unsigned char Type; };
template <> struct BitVector_<16> { typedef unsigned short Type; };
template <> struct BitVector_<32> { typedef unsigned int Type; };
template <> struct BitVector_<64> { typedef __uint64 Type; };
template <> struct BitVector_<255> { typedef __uint64 Type; };

// TODO(holtgrew): There is a lot of stuff defined within the class itself. A lot of it could be moved into global functions.

// template <typename T_, unsigned _size>
// const unsigned Tuple<T_, _size, Compressed>::BIT_SIZE = BitsPerValue<T_>::VALUE;
// template <typename T_, unsigned _size>
// const unsigned Tuple<T_, _size, Compressed>::BIT_MASK = (1 << Tuple<T_, _size, Compressed>::BIT_SIZE) - 1;
// template <typename T_, unsigned _size>
// const unsigned Tuple<T_, _size, Compressed>::MASK = (1 << (_size * Tuple<T_, _size, Compressed>::BIT_SIZE)) - 1;
// enum { size = _size };
// enum { BIT_SIZE = BitsPerValue<T_>::VALUE };
// enum { bitMASK = (1 << BIT_SIZE) - 1 };
// // TODO(holtgrew): The following two computations are bogus, in cases of overflow, MaxValue<CT>::VALUE should be used.
// enum { MASK = (1 << (size * BIT_SIZE)) - 1 };

// bit-compressed storage (space efficient)
#ifdef PLATFORM_WINDOWS
    #pragma pack(push,1)
#endif
template <typename T_, unsigned _size>
struct Tuple<T_, _size, Compressed>
{
    typedef T_ T;
    static const unsigned SIZE = _size;
    static const unsigned BIT_SIZE = BitsPerValue<T_>::VALUE;
    static const unsigned BIT_MASK = (1 << Tuple<T_, SIZE, Compressed>::BIT_SIZE) - 1;
    static const unsigned MASK = (1 << (SIZE * Tuple<T_, SIZE, Compressed>::BIT_SIZE)) - 1;
    // enum { size = _size };
    // enum { BIT_SIZE = BitsPerValue<T_>::VALUE };
    // enum { bitMask = (1 << BIT_SIZE) - 1 };
    // // TODO(holtgrew): The following two computations are bogus, in cases of overflow, MaxValue<CT>::VALUE should be used.
    // enum { mask = (1 << (size * BIT_SIZE)) - 1 };
    typedef typename BitVector_<BIT_SIZE * SIZE>::Type CT;

    // -----------------------------------------------------------------------
    // Members
    // -----------------------------------------------------------------------

    CT i;

    // -----------------------------------------------------------------------
    // Constructors
    // -----------------------------------------------------------------------

    // TODO(holtgrew): There is the unresolved issue whether the initialize costs critical performance. Since Tuples are PODs, it should be able to initialize Strings/arrays of them with memset().
    inline Tuple() : i(0)
    {
        SEQAN_ASSERT_LEQ(static_cast<__uint64>(BIT_SIZE * SIZE), static_cast<__uint64>(sizeof(CT) * 8));
    }

    // -----------------------------------------------------------------------
    // Subscription Operators;  Have to be declared in class.
    // -----------------------------------------------------------------------

    template <typename TPos>
    inline const T_
    operator[](TPos k) const
    {
        SEQAN_ASSERT_GEQ(static_cast<__int64>(k), 0);
        SEQAN_ASSERT_LT(static_cast<__int64>(k), static_cast<__int64>(SIZE));
        return (i >> (SIZE - 1 - k) * BIT_SIZE) & BIT_MASK;
    }

    // -----------------------------------------------------------------------
    // Assignment Operators;  Have to be declared in class.
    // -----------------------------------------------------------------------

    template <unsigned size__>
    inline Tuple operator=(Tuple<T_, size__, Compressed> const & _right)
    {
        i = _right.i;
        return *this;
    }

    // TODO(holtgrew): Move the following to global functions?

    template <typename TShiftSize>
    inline CT operator<<=(TShiftSize shift)
    {
        return i = (i << (shift * BIT_SIZE)) & MASK;
    }

    template <typename TShiftSize>
    inline CT operator<<(TShiftSize shift) const
    {
        return (i << (shift * BIT_SIZE)) & MASK;
    }

    template <typename TShiftSize>
    inline CT operator>>=(TShiftSize shift)
    {
        return i = (i >> (shift * BIT_SIZE));
    }

    template <typename TShiftSize>
    inline CT operator>>(TShiftSize shift) const
    {
        return i >> (shift * BIT_SIZE);
    }

    template <typename T>
    inline void operator|=(T const & t)
    {
        i |= t;
    }

    template <typename T, typename TSpec>
    inline void operator|=(SimpleType<T, TSpec> const & t)
    {
        i |= t.value;
    }

    inline CT* operator&()
    {
        return &i;
    }

    inline const CT* operator&() const
    {
        return &i;
    }

    // This to be inline because elements (like this tuple) of packed structs
    // can't be arguments.
    template <typename TPos, typename tmpS>
    inline tmpS const
    assignValueAt(TPos k, tmpS const source)
    {
        typedef Tuple<T_, _size, Compressed> Tup;
        typename Tup::CT MASK = Tup::BIT_MASK << ((_size - 1 - k) * BIT_SIZE);
        i = (i & ~MASK) | ((CT)ordValue(source) << ((_size - 1 - k) * BIT_SIZE));
        return source;
    }
}
#ifndef PLATFORM_WINDOWS
    __attribute__((packed))
#endif
    ;
#ifdef PLATFORM_WINDOWS
    #pragma pack(pop)
#endif

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// -----------------------------------------------------------------------
// Function assignValueAt()
// -----------------------------------------------------------------------

template <typename T_, unsigned _size, typename TPos>
T_
getValue(Tuple<T_, _size, Compressed> const & me,
         TPos k)
{
    SEQAN_ASSERT_GEQ(static_cast<__int64>(k), 0);
    SEQAN_ASSERT_LT(static_cast<__int64>(k), static_cast<__int64>(_size));
    // TODO(holtgrew): The following code is bogus, overflows can happen.
    return (me.i >> (_size - 1 - k) * BitsPerValue<T_>::VALUE) & ((1 << BitsPerValue<T_>::VALUE) - 1);
}

template <typename T_, unsigned _size, typename TPos>
T_
getValue(Tuple<T_, _size, Compressed> & me,
         TPos k)
{
    return getValue(const_cast<Tuple<T_, _size, Compressed> const &>(me), k);
}

// -----------------------------------------------------------------------
// Function assignValueAt()
// -----------------------------------------------------------------------

// TODO(holtgrew): Remove in favour of assignValue()?
// TODO(holtgrew): This is specialized for SimpleTypes, do those specializations not belong there?

template <typename T_, unsigned _size, typename tmpS, typename TPos>
inline tmpS const
assignValueAt(Tuple<T_, _size, Compressed> & me,
              TPos k,
              tmpS const source)
{
    typedef Tuple<T_, _size, Compressed> Tup;
    typename Tup::CT MASK = Tup::BIT_MASK << ((_size - 1 - k) * me.BIT_SIZE);
    me.i = (me.i & ~MASK) | source << ((_size - 1 - k) * me.BIT_SIZE);
    return source;
}

template <typename T_, typename tmpS, typename Spec_, unsigned _size, typename TPos>
inline SimpleType<tmpS, Spec_> const &
assignValueAt(Tuple<T_, _size, Compressed> & me,
              TPos k,
              SimpleType<tmpS, Spec_> const & source)
{
    typedef Tuple<T_, _size, Compressed> Tup;
    typename Tup::CT MASK = Tup::BIT_MASK << ((_size - 1 - k) * me.BIT_SIZE);
    me.i = (me.i & ~MASK) | source.value << ((_size - 1 - k) * me.BIT_SIZE);
    return source;
}

// -----------------------------------------------------------------------
// Function assignValue()
// -----------------------------------------------------------------------

template <typename T_, unsigned _size, typename tmpS, typename TPos>
inline tmpS const
assignValue(Tuple<T_, _size, Compressed> & me,
            TPos k,
            tmpS const source)
{
    return assignValueAt(me, k, source);
}

template <typename T_, typename tmpS, typename Spec_, unsigned _size, typename TPos>
inline SimpleType<tmpS, Spec_> const &
assignValue(Tuple<T_, _size, Compressed> & me,
            TPos k,
            SimpleType<tmpS, Spec_> const & source)
{
    return assignValueAt(me, k, source);
}

// -----------------------------------------------------------------------
// Function setValue()
// -----------------------------------------------------------------------

template <typename T_, unsigned _size, typename tmpS, typename TPos>
inline tmpS const
setValue(Tuple<T_, _size, Compressed> & me,
         TPos k,
         tmpS const source)
{
    return assignValue(me, k, source);
}

template <typename T_, typename tmpS, typename Spec_, unsigned _size, typename TPos>
inline SimpleType<tmpS, Spec_> const &
setValue(Tuple<T_, _size, Compressed> & me,
         TPos k,
         SimpleType<tmpS, Spec_> const & source)
{
    return assignValue(me, k, source);
}

// -----------------------------------------------------------------------
// Function moveValue()
// -----------------------------------------------------------------------

template <typename T_, unsigned _size, typename tmpS, typename TPos>
inline tmpS const
moveValue(Tuple<T_, _size, Compressed> & me,
          TPos k,
          tmpS const source)
{
    return assignValue(me, k, source);
}

template <typename T_, typename tmpS, typename Spec_, unsigned _size, typename TPos>
inline SimpleType<tmpS, Spec_> const &
moveValue(Tuple<T_, _size, Compressed> & me,
          TPos k,
          SimpleType<tmpS, Spec_> const & source)
{
    return assignValue(me, k, source);
}

// -----------------------------------------------------------------------
// Function shiftLeft()
// -----------------------------------------------------------------------

// Optimized version for compressed tuple using just one word.
template <typename T_, unsigned _size>
inline void shiftLeft(Tuple<T_, _size, Compressed> & me)
{
    me <<= 1;
}

// -----------------------------------------------------------------------
// Function shiftRight()
// -----------------------------------------------------------------------

template <typename T_, unsigned _size>
inline void shiftRight(Tuple<T_, _size, Compressed> & me)
{
    me >>= 1;
}

// -----------------------------------------------------------------------
// Function clear()
// -----------------------------------------------------------------------
 
template <typename T_, unsigned _size>
inline void clear(Tuple<T_, _size, Compressed> & me)
{
    me.i = 0; 
}

// -----------------------------------------------------------------------
// Function operator<()
// -----------------------------------------------------------------------

// Optimized version for compressed tuple using just one word.
template <typename T_, unsigned _size>
inline bool operator<(Tuple<T_, _size, Compressed> const &_left,
                      Tuple<T_, _size, Compressed> const & _right)
{
    return _left.i < _right.i;
}

template <typename T_, unsigned _size>
inline bool operator<(Tuple<T_, _size, Compressed> &_left,
                      Tuple<T_, _size, Compressed> & _right)
{
    return _left.i < _right.i;
}

// -----------------------------------------------------------------------
// Function operator>()
// -----------------------------------------------------------------------

// Optimized version for compressed tuple using just one word.
template <typename T_, unsigned _size>
inline bool operator>(Tuple<T_, _size, Compressed> const &_left,
                      Tuple<T_, _size, Compressed> const & _right)
{
    return _left.i > _right.i;
}

template <typename T_, unsigned _size>
inline bool operator>(Tuple<T_, _size, Compressed> &_left,
                      Tuple<T_, _size, Compressed> & _right)
{
    return _left.i > _right.i;
}

// -----------------------------------------------------------------------
// Function operator<=()
// -----------------------------------------------------------------------

// Optimized version for compressed tuple using just one word.
template <typename T_, unsigned _size>
inline bool operator<=(Tuple<T_, _size, Compressed> const &_left,
                       Tuple<T_, _size, Compressed> const & _right)
{
    return !operator>(_left, _right);
}

template <typename T_, unsigned _size>
inline bool operator<=(Tuple<T_, _size, Compressed> &_left,
                       Tuple<T_, _size, Compressed> & _right)
{
    return !operator>(_left, _right);
}

// -----------------------------------------------------------------------
// Function operator>()
// -----------------------------------------------------------------------

// Optimized version for compressed tuple using just one word.
template <typename T_, unsigned _size>
inline bool operator>=(Tuple<T_, _size, Compressed> const &_left,
                       Tuple<T_, _size, Compressed> const & _right)
{
    return !operator<(_left, _right);
}

template <typename T_, unsigned _size>
inline bool operator>=(Tuple<T_, _size, Compressed> &_left,
                       Tuple<T_, _size, Compressed> & _right)
{
    return !operator<(_left, _right);
}

// -----------------------------------------------------------------------
// Function operator==()
// -----------------------------------------------------------------------

// Optimized version for compressed tuple using just one word.
template <typename T_, unsigned _size>
inline bool operator==(Tuple<T_, _size, Compressed> const & _left,
                       Tuple<T_, _size, Compressed> const & _right)
{
    return _left.i == _right.i;
}

template <typename T_, unsigned _size>
inline bool operator==(Tuple<T_, _size, Compressed> & _left,
                       Tuple<T_, _size, Compressed> & _right)
{
    return _left.i == _right.i;
}

// -----------------------------------------------------------------------
// Function operator!=()
// -----------------------------------------------------------------------

// Optimized version for compressed tuple using just one word.
template <typename T_, unsigned _size>
inline bool operator!=(Tuple<T_, _size, Compressed> const & _left,
                       Tuple<T_, _size, Compressed> const & _right)
{
    return !operator==(_left, _right);
}

template <typename T_, unsigned _size>
inline bool operator!=(Tuple<T_, _size, Compressed> & _left,
                       Tuple<T_, _size, Compressed> & _right)
{
    return !operator==(_left, _right);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_TUPLE_BIT_COMPRESSED_H_
