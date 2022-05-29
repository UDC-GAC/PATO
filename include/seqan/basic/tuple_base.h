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
// Tuple base class.
// ==========================================================================

// TODO(holtgrew): What about move construction? Useful for pairs of strings and such. Tricky to implement since ints have no move constructor, for example.

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_TUPLE_BASE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_TUPLE_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Class.Tuple:
..cat:Aggregates
..concept:Concept.Aggregate
..summary:A plain fixed-length string.
..signature:Tuple<T, SIZE[, TSpec]>
..param.T:The value type, that is the type of characters stored in the tuple.
..param.SIZE:The size/length of the tuple.
...remarks:In contrast to @Class.String@ the length of Tuple is fixed.
..param.TSpec:The specializing type.
...default:$void$, no compression (faster access).
..include:seqan/basic.h
*/

template <typename T_, unsigned _size, typename TSpec = void>
struct Tuple
{
    typedef T_ T;
    static const unsigned SIZE;

    // -----------------------------------------------------------------------
    // Members
    // -----------------------------------------------------------------------

    T_ i[_size];
    
    // -----------------------------------------------------------------------
    // Subscription Operators;  Have to be declared in class.
    // -----------------------------------------------------------------------

    // TODO(holtgrew): Return Value<>::Type?

    template <typename TPos>
    inline T_ &
    operator[](TPos k)
    {
        SEQAN_ASSERT_GEQ(static_cast<__int64>(k), 0);
        SEQAN_ASSERT_LT(static_cast<__int64>(k), static_cast<__int64>(SIZE));
        return i[k];
    }

    template <typename TPos>
    inline const T_ &
    operator[](TPos k) const
    {
        SEQAN_ASSERT_GEQ(static_cast<__int64>(k), 0);
        SEQAN_ASSERT_LT(static_cast<__int64>(k), static_cast<__int64>(SIZE));
        return i[k];
        
    }

    // TODO(holtgrew): What's this?

    inline T_ *
    operator&() { return i; }

    inline const T_ *
    operator&() const { return i; }

    // This has to be inline because elements (like this tuple) of packed
    // structs can't be arguments.
    template <typename TPos, typename tmpS>
    inline tmpS const
    assignValueAt(TPos k, tmpS const source)
    {
        return i[k] = source;
    }
};

template <typename T_, unsigned _size, typename TSpec>
const unsigned Tuple<T_, _size, TSpec>::SIZE = _size;

// ============================================================================
// Metafunctions
// ============================================================================

// -----------------------------------------------------------------------
// Metafunction LENGTH
// -----------------------------------------------------------------------

///.Metafunction.LENGTH.param.T.type:Class.Tuple

template <typename T_, unsigned _size, typename TSpec>
struct LENGTH<Tuple<T_, _size, TSpec> >
{
    enum { VALUE = _size };
};

// -----------------------------------------------------------------------
// Metafunction Value
// -----------------------------------------------------------------------

///.Metafunction.Value.param.T.type:Class.Tuple

template <typename T_, unsigned _size, typename TSpec>
struct Value<Tuple<T_, _size, TSpec> >
{
    typedef T_ Type;
};

// -----------------------------------------------------------------------
// Metafunction Spec
// -----------------------------------------------------------------------

///.Metafunction.Spec.param.T.type:Class.Tuple

template <typename T_, unsigned _size, typename TSpec>
struct Spec<Tuple<T_, _size, TSpec> >
{
    typedef TSpec Type;
};

// ============================================================================
// Functions
// ============================================================================

// -----------------------------------------------------------------------
// Function operator<<();  Stream Output.
// -----------------------------------------------------------------------

template <typename T_, unsigned _size, typename TSpec>
inline std::ostream &
operator<<(std::ostream & out, Tuple<T_,_size,TSpec> const &a) {
    out << "[";
    if (a.SIZE > 0)
            out << a[0];
    for(unsigned j = 1; j < a.SIZE; ++j)
        out << " " << a[j];
    out << "]";
    return out;
}

// ----------------------------------------------------------------------------
// Function set().
// ----------------------------------------------------------------------------

template <typename TTuple1, typename TTuple2>
struct TupleMoveSetWorkerContext_
{
    TTuple1 & t1;
    TTuple2 & t2;

    TupleMoveSetWorkerContext_(TTuple1 & _t1, TTuple2 & _t2)
            : t1(_t1), t2(_t2)
    {}
};

struct TupleSetWorker_
{
    template <typename TArg>
    static inline void body(TArg & arg, unsigned I)
    {
        set(arg.t1.i[I - 1], arg.t2.i[I - 1]);
    }
};

template <typename T_, unsigned _size>
inline void
set(Tuple<T_, _size, void> & t1, Tuple<T_, _size, void> const & t2)
{
    typedef Tuple<T_, _size, void> TTuple1;
    typedef Tuple<T_, _size, void> const TTuple2;
    TupleMoveSetWorkerContext_<TTuple1, TTuple2> context(t1, t2);
    Loop<TupleSetWorker_, _size>::run(context);
}

template <typename T_, unsigned _size>
inline void
set(Tuple<T_, _size, void> & t1, Tuple<T_, _size, void> & t2)
{
    set(t1, const_cast<Tuple<T_, _size, void> const &>(t2));
}

// ----------------------------------------------------------------------------
// Function move().
// ----------------------------------------------------------------------------

struct TupleMoveWorker_
{
    template <typename TArg>
    static inline void body(TArg & arg, unsigned I)
    {
        move(arg.t1.i[I - 1], arg.t2.i[I - 1]);
    }
};

template <typename T_, unsigned _size>
inline void
move(Tuple<T_, _size, void> & t1, Tuple<T_, _size, void> & t2)
{
    typedef Tuple<T_, _size, void> TTuple1;
    typedef Tuple<T_, _size, void> TTuple2;
    TupleMoveSetWorkerContext_<TTuple1, TTuple2> context(t1, t2);
    Loop<TupleMoveWorker_, _size>::run(context);
}

// -----------------------------------------------------------------------
// Function assignValueAt()
// -----------------------------------------------------------------------

// TODO(holtgrew): Remove in favour of assignValue()! Remove function definition here, AT LEAST

template <typename TObject, typename TPos, typename TSource>
inline TSource & 
assignValueAt(TObject & me, TPos k, TSource &source)
{
    assign(value(me, k), source);
    return source;
}

template <typename TObject, typename TPos, typename TSource>
inline TSource const & 
assignValueAt(TObject & me, TPos k, TSource const & source)
{
    assign(value(me, k), source);
    return source;
}

template <typename T_, unsigned _size, typename tmpS, typename TPos>
inline tmpS const
assignValueAt(Tuple<T_, _size, void> & me, TPos k, tmpS const source)
{
    SEQAN_CHECK((k < _size), "Invalid position, k = %u, _size = %u.", unsigned(k), unsigned(_size));
    return me.i[k] = source;
}

// -----------------------------------------------------------------------
// Function getValue()
// -----------------------------------------------------------------------

template <typename T_, unsigned _size, typename TPos>
inline T_
getValue(Tuple<T_, _size, void> & me, TPos k)
{
    SEQAN_CHECK((unsigned(k) < _size), "Invalid position, k = %u, _size = %u.", unsigned(k), unsigned(_size));
    return me.i[k];
}

// -----------------------------------------------------------------------
// Function assignValue()
// -----------------------------------------------------------------------

template <typename T_, unsigned _size, typename tmpS, typename TPos>
inline void
assignValue(Tuple<T_, _size, void> & me, TPos k, tmpS const & source)
{
    SEQAN_CHECK((unsigned(k) < _size), "Invalid position, k = %u, _size = %u.", unsigned(k), unsigned(_size));
    assign(me.i[k], source);
}

// -----------------------------------------------------------------------
// Function setValue()
// -----------------------------------------------------------------------

template <typename T_, unsigned _size, typename tmpS, typename TPos>
inline void
setValue(Tuple<T_, _size, void> & me, TPos k, tmpS const & source)
{
    SEQAN_CHECK((unsigned(k) < _size), "Invalid position, k = %u, _size = %u.", unsigned(k), unsigned(_size));
    set(me.i[k], source);
}

// -----------------------------------------------------------------------
// Function moveValue()
// -----------------------------------------------------------------------

template <typename T_, unsigned _size, typename tmpS, typename TPos>
inline void
moveValue(Tuple<T_, _size, void> & me, TPos k, tmpS & source)
{
    SEQAN_CHECK((unsigned(k) < _size), "Invalid position, k = %u, _size = %u.", unsigned(k), unsigned(_size));
    move(me.i[k], source);
}

// -----------------------------------------------------------------------
// Function shiftLeft()
// -----------------------------------------------------------------------

// TODO(holtgrew): Document!

struct TupleShiftLeftWorker_
{
    template <typename TArg>
    static inline void body(TArg & arg, unsigned I)
    {
        arg[I-1] = arg[I];  // TODO(holtgrew): Do we really want assignment or movement here?
    }
};

template <typename T_, unsigned _size, typename TSpec>
inline void shiftLeft(Tuple<T_, _size, TSpec> &me)
{
    Loop<TupleShiftLeftWorker_, _size - 1>::run(me.i);
}

// -----------------------------------------------------------------------
// Function shiftRight()
// -----------------------------------------------------------------------

// TODO(holtgrew): Document!

struct TupleShiftRightWorker_
{
    template <typename TArg>
    static inline void body(TArg & arg, unsigned I)
    {
        arg[I] = arg[I - 1];  // TODO(holtgrew): Do we really want assignment or movement here?
    }
};

template <typename T_, unsigned _size, typename TSpec>
inline void shiftRight(Tuple<T_, _size, TSpec> & me)
{
    LoopReverse<TupleShiftRightWorker_, _size - 1>::run(me.i);
}

// -----------------------------------------------------------------------
// Function length()
// -----------------------------------------------------------------------

template <typename T_, unsigned _size, typename TSpec>
inline unsigned length(Tuple<T_, _size, TSpec> const &)
{
    return _size;
}

// -----------------------------------------------------------------------
// Function clear()
// -----------------------------------------------------------------------

template <typename T_, unsigned _size, typename TSpec>
inline void clear(Tuple<T_, _size, TSpec> & me)
{
   memset<sizeof(me.i), 0>(&(me.i));
}

// -----------------------------------------------------------------------
// Function operator==()
// -----------------------------------------------------------------------

template <typename TTuple>
struct ComparisonWorkerContext_
{
    int result;
    TTuple const & left;
    TTuple const & right;

    ComparisonWorkerContext_(int b, TTuple const & l, TTuple const & r)
            : result(b), left(l), right(r)
    {}
};

struct TupleComparisonWorkerEq_
{
    template <typename TArg>
    static inline void body(TArg & arg, unsigned I)
    {
        if (arg.result != 1)
            return;
        if (arg.left.i[I - 1] != arg.right.i[I - 1])
            arg.result = 0;
    }
};

template <typename T_, unsigned _size, typename TSpec>
inline bool
operator==(Tuple<T_, _size, TSpec> const & _left,
           Tuple<T_, _size, TSpec> const & _right)
{
    typedef Tuple<T_, _size, TSpec> TTuple;
    ComparisonWorkerContext_<TTuple> context(1, _left, _right);
    Loop<TupleComparisonWorkerEq_, _size>::run(context);
    return context.result == 1;
}

// -----------------------------------------------------------------------
// Function operator!=()
// -----------------------------------------------------------------------


template <typename T_, unsigned _size, typename TSpec>
inline bool
operator!=(Tuple<T_, _size, TSpec> const & _left,
           Tuple<T_, _size, TSpec> const & _right)
{
    return !operator==(_left, _right);
}

// -----------------------------------------------------------------------
// Function operator<()
// -----------------------------------------------------------------------

struct TupleComparisonWorkerLt_
{
    template <typename TArg>
    static inline void body(TArg & arg, unsigned I)
    {
        if (arg.result != -1)
            return;
        if (arg.left.i[I - 1] == arg.right.i[I - 1])
            return;
        if (arg.left.i[I - 1] < arg.right.i[I - 1])
            arg.result = 1;
        if (arg.left.i[I - 1] > arg.right.i[I - 1])
            arg.result = 0;
    }
};

template <typename T_, unsigned _size, typename TSpec>
inline bool
operator<(Tuple<T_, _size, TSpec> const & _left,
          Tuple<T_, _size, TSpec> const & _right)
{
    typedef Tuple<T_, _size, TSpec> TTuple;
    ComparisonWorkerContext_<TTuple> context(-1, _left, _right);
    Loop<TupleComparisonWorkerLt_, _size>::run(context);
    return context.result == 1;
}

// -----------------------------------------------------------------------
// Function operator>()
// -----------------------------------------------------------------------

struct TupleComparisonWorkerGt_
{
    template <typename TArg>
    static inline void body(TArg & arg, unsigned I)
    {
        if (arg.result != -1)
            return;
        if (arg.left.i[I - 1] == arg.right.i[I - 1])
            return;
        if (arg.left.i[I - 1] > arg.right.i[I - 1])
            arg.result = 1;
        if (arg.left.i[I - 1] < arg.right.i[I - 1])
            arg.result = 0;
    }
};

template <typename T_, unsigned _size, typename TSpec>
inline bool
operator>(Tuple<T_, _size, TSpec> const & _left,
          Tuple<T_, _size, TSpec> const & _right)
{
    typedef Tuple<T_, _size, TSpec> TTuple;
    ComparisonWorkerContext_<TTuple> context(-1, _left, _right);
    Loop<TupleComparisonWorkerGt_, _size>::run(context);
    return context.result == 1;
}

// -----------------------------------------------------------------------
// Function operator<=()
// -----------------------------------------------------------------------

template <typename T_, unsigned _size, typename TSpec>
inline bool
operator<=(Tuple<T_, _size, TSpec> const & _left,
           Tuple<T_, _size, TSpec> const & _right)
{
    return !operator>(_left, _right);
}

// -----------------------------------------------------------------------
// Function operator>=()
// -----------------------------------------------------------------------

template <typename T_, unsigned _size, typename TSpec>
inline bool
operator>=(Tuple<T_, _size, TSpec> const & _left,
           Tuple<T_, _size, TSpec> const & _right)
{
    return !operator<(_left, _right);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_TUPLE_BASE_H_
