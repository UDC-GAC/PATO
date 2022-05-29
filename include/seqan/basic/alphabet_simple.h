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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// The SimpleType alphabet type.
// ==========================================================================

#ifndef SEQAN_BASIC_ALPHABET_SIMPLE_H_
#define SEQAN_BASIC_ALPHABET_SIMPLE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Class.SimpleType:
..cat:Basic
..summary:Implementation for "simple" types.
..signature:SimpleType<TValue, TSpec>
..param.TValue:Type that stores the values of an instance.
...remarks:TValue must be a simple type.
...metafunction:Metafunction.Value
..param.TSpec:Specialization tag.
...metafunction:Metafunction.Spec
..remarks:
...text:A "simple type" is a C++ type that can be constructed without constructor,
destructed without destructor and copied without copy constructor or assignment operator.
All basic types (like $char$, $int$ or $float$) are simple. Pointers, references and arrays of
simple types are simple.
POD types ("plain old data types"), that are - simplified spoken - C++-types that already existed in C,
are simple too. 
...text:Arrays of simple types can be copied very fast by memory manipulation routines, 
but the default implementation of functions like @Function.arrayCopyForward@ and @Function.arrayCopy@
are not optimized for simple types this way.
But for classes derived from $SimpleType$, optimized variants of array manipulation functions are applied. 
...text:Note that simple types need not to be derived or specialized from $SimpleType$, but
it could be convenient to do so.
..implements:Concept.Simple Type
..include:seqan/basic.h
*/

template <typename TValue, typename TSpec>
class SimpleType
{
public:
    // ------------------------------------------------------------------------
    // Members;  Have to be defined in class.
    // ------------------------------------------------------------------------

    TValue value;

    // ------------------------------------------------------------------------
    // Constructors;  Have to be defined in class.
    // ------------------------------------------------------------------------

    SimpleType() 
    {
        SEQAN_CHECKPOINT;
    }

    SimpleType(SimpleType const & other)
    {
        SEQAN_CHECKPOINT;
        assign(*this, other);
    }

    template <typename T> 
    SimpleType(T const & other) 
    {
        SEQAN_CHECKPOINT;
        assign(*this, other);
    }

    // ------------------------------------------------------------------------
    // Assignment Operator;  Have to be defined in class.
    // ------------------------------------------------------------------------

    inline SimpleType &
    operator=(SimpleType const & other) 
    { 
        SEQAN_CHECKPOINT;
        assign(*this, other);
        return *this;
    }

    template <typename T>
    inline SimpleType &
    operator=(T const & other) 
    { 
        SEQAN_CHECKPOINT;
        assign(*this, other);
        return *this;
    }

    // ------------------------------------------------------------------------
    // Conversion Operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

    // Class.SimpleType specifies type conversion operators for all built-in
    // integer types since there is no way to extend the build-in types with
    // copy and assignment constructors in C++.
    //
    // This cannot be a template since it would conflict to the template
    // constructor.

    operator __int64() const
    {
        SEQAN_CHECKPOINT;
        __int64 c;
        assign(c, *this);
        return c;
    }

    operator __uint64() const
    {
        SEQAN_CHECKPOINT;
        __uint64 c;
        assign(c, *this);
        return c;
    }

    operator int() const
    {
        SEQAN_CHECKPOINT;
        int c;
        assign(c, *this);
        return c;
    }

    operator unsigned int() const
    {
        SEQAN_CHECKPOINT;
        unsigned int c;
        assign(c, *this);
        return c;
    }

    operator short() const
    {
        SEQAN_CHECKPOINT;
        short c;
        assign(c, *this);
        return c;
    }

    operator unsigned short() const
    {
        SEQAN_CHECKPOINT;
        unsigned short c;
        assign(c, *this);
        return c;
    }

    operator char() const
    {
        SEQAN_CHECKPOINT;
        char c;
        assign(c, *this);
        return c;
    }

    operator signed char() const
    {
        SEQAN_CHECKPOINT;
        signed char c;
        assign(c, *this);
        return c;
    }

    operator unsigned char() const
    {
        SEQAN_CHECKPOINT;
        unsigned char c;
        assign(c, *this);
        return c;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction IsSimple
// ----------------------------------------------------------------------------

///.Metafunction.IsSimple.param.T.type:Class.SimpleType

template <typename TValue, typename TSpec>
struct IsSimple<SimpleType<TValue, TSpec> >
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

///.Metafunction.Value.param.T.type:Class.SimpleType

template <typename TValue, typename TSpec>
struct Value<SimpleType<TValue, TSpec> >
{
    typedef TValue Type;
};

template <typename TValue, typename TSpec>
struct Value<SimpleType<TValue, TSpec> const>
{
    typedef TValue const Type;
};

// ----------------------------------------------------------------------------
// Metafunction MinValue
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct MinValue<SimpleType<TValue, TSpec> > : MinimumValueUnsigned_<TValue> {};

template <typename TValue, typename TSpec>
inline SimpleType<TValue, TSpec>
infimumValueImpl(SimpleType<TValue, TSpec> *)
{
    SEQAN_CHECKPOINT;
    return MinValue<SimpleType<TValue, TSpec> >::VALUE;
}

// ----------------------------------------------------------------------------
// Metafunction MaxValue
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct MaxValue<SimpleType<TValue, TSpec> >
{
    static const TValue VALUE;
};

template <typename TValue, typename TSpec>
const TValue MaxValue<SimpleType<TValue, TSpec> >::VALUE = ((TValue)ValueSize<SimpleType<TValue, TSpec> >::VALUE - 1);

template <typename TValue, typename TSpec>
inline SimpleType<TValue, TSpec>
supremumValueImpl(SimpleType<TValue, TSpec> *)
{
    SEQAN_CHECKPOINT;
    return MaxValue<SimpleType<TValue, TSpec> >::VALUE;
}

// ----------------------------------------------------------------------------
// Metafunction Spec
// ----------------------------------------------------------------------------

///.Metafunction.Spec.param.T.type:Class.SimpleType

template <typename TValue, typename TSpec>
struct Spec<SimpleType<TValue, TSpec> >
{
    typedef TSpec Type;
};

template <typename TValue, typename TSpec>
struct Spec<SimpleType<TValue, TSpec> const>
{
    typedef TSpec Type;
};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

// TODO(holtgrew): Remove auto-sequence!
template <typename TValue, typename TSpec>
struct Iterator<SimpleType<TValue, TSpec>, Standard>
{
    typedef SimpleType<TValue, TSpec> * Type;
//  typedef Iter<SimpleType<TValue, TSpec>, SimpleIterator> * Type;
};

template <typename TValue, typename TSpec>
struct Iterator<SimpleType<TValue, TSpec> const, Standard>
{
    typedef SimpleType<TValue, TSpec> const * Type;
//  typedef Iter<SimpleType<TValue, TSpec> const, SimpleIterator> * Type;
};

// ----------------------------------------------------------------------------
// Metafunction CompareType
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Metafunction CompareType
// ----------------------------------------------------------------------------

/**
.Metafunction.CompareType:
..cat:Basic
..summary:Type to convert other types for comparisons.
..signature:CompareType<TLeft, TRight>::Type
..param.TLeft:Type of the left operand of a comparison.
..param.TRight:Type of the right operand of a comparison.
..return.Type:The Type in which the arguments are converted in order to compare them.
..remarks:Comparisons are for example operators like $==$ or $<$.
..remarks.text:Note that there is no rule that guarantees that $CompareType<T1, T2>::Type$
is the same as $CompareType<T2, T1>::Type$. It is also possible, that only one of these
two types is defined.
..remarks.text:This metafunction is used for the implementation of
comparisons that involve @Class.SimpleType@.
..include:seqan/basic.h
*/
// TODO(doering): muss geprueft werden, ob diese Metafunktion noch ausgeweitet oder aber versteckt wird.
// TODO(holtgrew): Is this at the right place here? Do we need it for all types?
// TODO(holtgrew): CompareType is not symmetric because of class instantiation conflicts. Evaluate these problems and possibly fix them.

template <typename TLeft, typename TRight>
struct CompareType;

template <typename T>
struct CompareType<T, T>
{
    typedef T Type;
};

template <typename TValue, typename TSpec, typename TRight>
struct CompareType<SimpleType<TValue, TSpec>, TRight>
{
    typedef TRight Type;
};

// ============================================================================
// Functions
// ============================================================================

// TODO(holtgrew): Why are functions with Proxies defined here? Better in Proxy header? Or externalize?

// ----------------------------------------------------------------------------
// Function convertImpl()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Document

template <typename TTarget, typename T, typename TSourceValue, typename TSourceSpec>
inline typename RemoveConst_<TTarget>::Type
convertImpl(Convert<TTarget, T> const,
            SimpleType<TSourceValue, TSourceSpec> const & source_)
{
    SEQAN_CHECKPOINT;
    typename RemoveConst_<TTarget>::Type target_;
    assign(target_, source_);
    return target_;
}

// ----------------------------------------------------------------------------
// Function operator<<()
// ----------------------------------------------------------------------------

template <typename TStream, typename TValue, typename TSpec>
inline TStream &
operator<<(TStream & stream, 
           SimpleType<TValue, TSpec> const & data)
{
    SEQAN_CHECKPOINT;
    stream << convert<char>(data);
    return stream;
}

// ----------------------------------------------------------------------------
// Function operator>>()
// ----------------------------------------------------------------------------

template <typename TStream, typename TValue, typename TSpec>
inline TStream &
operator>>(TStream & stream, 
           SimpleType<TValue, TSpec> & data)
{
    SEQAN_CHECKPOINT;
    char c;
    stream >> c;
    assign(data, c);
    return stream;
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

///.Function.assign.param.target.type:Class.SimpleType
///.Function.assign.param.source.type:Class.SimpleType

template <typename TTargetValue, typename TTargetSpec, typename TSourceValue, typename TSourceSpec>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
       SimpleType<TSourceValue, TSourceSpec> & source)
{
    SEQAN_CHECKPOINT;
    target.value = source.value;
}

template <typename TTargetValue, typename TTargetSpec, typename TSourceValue, typename TSourceSpec>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
       SimpleType<TSourceValue, TSourceSpec> const & source)
{
    SEQAN_CHECKPOINT;
    target.value = source.value;
}

template <typename TTargetValue, typename TTargetSpec, typename TSource>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
       TSource & source)
{
    SEQAN_CHECKPOINT;
    target.value = source;
}

template <typename TTargetValue, typename TTargetSpec, typename TSource>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
       TSource const & source)
{
    SEQAN_CHECKPOINT;
    target.value = source;
}

// Assign Proxy to SimpleType 
// NOTE(doering): Diese Funktionen wurden noetig wegen eines seltsamen VC++-Verhaltens
// TODO(holtgrew): Still necessary with dropped 2003 support?

template <typename TTargetValue, typename TTargetSpec, typename TSourceSpec>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
       Proxy<TSourceSpec> & source)
{
    SEQAN_CHECKPOINT;
    target.value = getValue(source);
}

template <typename TTargetValue, typename TTargetSpec, typename TSourceSpec>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
       Proxy<TSourceSpec> const & source)
{
    SEQAN_CHECKPOINT;
    target.value = getValue(source);
}

//INTEGRAL TYPES
// NOTE(doering): It is not possible to write a single function here since "assign" must be specialized for the first argument at the first place

template <typename TValue, typename TSpec>
inline void 
assign(__int64 & c_target, 
       SimpleType<TValue, TSpec> & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

template <typename TValue, typename TSpec>
inline void 
assign(__int64 & c_target, 
       SimpleType<TValue, TSpec> const & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

template <typename TValue, typename TSpec>
inline void 
assign(__uint64 & c_target, 
       SimpleType<TValue, TSpec> & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

template <typename TValue, typename TSpec>
inline void 
assign(__uint64 & c_target, 
       SimpleType<TValue, TSpec> const & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

template <typename TValue, typename TSpec>
inline void 
assign(int & c_target, 
       SimpleType<TValue, TSpec> & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

template <typename TValue, typename TSpec>
inline void 
assign(int & c_target, 
       SimpleType<TValue, TSpec> const & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

template <typename TValue, typename TSpec>
inline void 
assign(unsigned int & c_target, 
       SimpleType<TValue, TSpec> & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

template <typename TValue, typename TSpec>
inline void 
assign(unsigned int & c_target, 
       SimpleType<TValue, TSpec> const & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

template <typename TValue, typename TSpec>
inline void 
assign(short & c_target, 
       SimpleType<TValue, TSpec> & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

template <typename TValue, typename TSpec>
inline void 
assign(short & c_target, 
       SimpleType<TValue, TSpec> const & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

template <typename TValue, typename TSpec>
inline void 
assign(unsigned short & c_target, 
       SimpleType<TValue, TSpec> & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

template <typename TValue, typename TSpec>
inline void 
assign(unsigned short & c_target, 
       SimpleType<TValue, TSpec> const & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

template <typename TValue, typename TSpec>
inline void 
assign(char & c_target, 
       SimpleType<TValue, TSpec> & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

template <typename TValue, typename TSpec>
inline void 
assign(char & c_target, 
       SimpleType<TValue, TSpec> const & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

template <typename TValue, typename TSpec>
inline void 
assign(signed char & c_target, 
       SimpleType<TValue, TSpec> & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

template <typename TValue, typename TSpec>
inline void 
assign(signed char & c_target, 
       SimpleType<TValue, TSpec> const & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

template <typename TValue, typename TSpec>
inline void 
assign(unsigned char & c_target, 
       SimpleType<TValue, TSpec> & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

template <typename TValue, typename TSpec>
inline void 
assign(unsigned char & c_target, 
       SimpleType<TValue, TSpec> const & source)
{
    SEQAN_CHECKPOINT;
    c_target = source.value;
}

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator==(SimpleType<TValue, TSpec> const & left_, 
           TRight const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TValue, TSpec> TLeft;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator==(TLeft const & left_, 
           SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TValue, TSpec> TRight;
    typedef typename CompareType<TRight, TLeft>::Type TCompareType;
    return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator==(SimpleType<TLeftValue, TLeftSpec> const & left_, 
           SimpleType<TRightValue, TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
    typedef SimpleType<TRightValue, TRightSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator==(SimpleType<TValue, TSpec> const & left_, 
           SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    return convert<TValue>(left_) == convert<TValue>(right_);
}

template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator==(Proxy<TSpec> const & left_, 
           SimpleType<TValue, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TLeft;
    typedef SimpleType<TValue, TSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator==(SimpleType<TValue, TSpec2> const & left_,
           Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TValue, TSpec> TLeft;
    typedef Proxy<TSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

// ----------------------------------------------------------------------------
// Function operator!=()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator!=(SimpleType<TValue, TSpec> const & left_, 
           TRight const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TValue, TSpec> TLeft;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator!=(TLeft const & left_, 
           SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TValue, TSpec> TRight;
    typedef typename CompareType<TRight, TLeft>::Type TCompareType;
    return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator!=(SimpleType<TLeftValue, TLeftSpec> const & left_, 
           SimpleType<TRightValue, TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
    typedef SimpleType<TRightValue, TRightSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator!=(SimpleType<TValue, TSpec> const & left_, 
           SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    return convert<TValue>(left_) != convert<TValue>(right_);
}

template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator!=(Proxy<TSpec> const & left_, 
           SimpleType<TValue, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TLeft;
    typedef SimpleType<TValue, TSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}
template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator!=(SimpleType<TValue, TSpec2> const & left_,
           Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TValue, TSpec> TLeft;
    typedef Proxy<TSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

// ----------------------------------------------------------------------------
// Function operator<()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator<(SimpleType<TValue, TSpec> const & left_, 
          TRight const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TValue, TSpec> TLeft;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator<(TLeft const & left_, 
          SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TValue, TSpec> TRight;
    typedef typename CompareType<TRight, TLeft>::Type TCompareType;
    return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator<(SimpleType<TLeftValue, TLeftSpec> const & left_, 
          SimpleType<TRightValue, TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
    typedef SimpleType<TRightValue, TRightSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator<(SimpleType<TValue, TSpec> const & left_, 
          SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    return convert<TValue>(left_) < convert<TValue>(right_);
}

template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator<(Proxy<TSpec> const & left_, 
          SimpleType<TValue, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TLeft;
    typedef SimpleType<TValue, TSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator<(SimpleType<TValue, TSpec2> const & left_,
          Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TValue, TSpec> TLeft;
    typedef Proxy<TSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

// ----------------------------------------------------------------------------
// Function operator<=()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator<=(SimpleType<TValue, TSpec> const & left_, 
           TRight const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TValue, TSpec> TLeft;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator<=(TLeft const & left_, 
           SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TValue, TSpec> TRight;
    typedef typename CompareType<TRight, TLeft>::Type TCompareType;
    return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator<=(SimpleType<TLeftValue, TLeftSpec> const & left_, 
           SimpleType<TRightValue, TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
    typedef SimpleType<TRightValue, TRightSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator<=(SimpleType<TValue, TSpec> const & left_, 
           SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    return convert<TValue>(left_) <= convert<TValue>(right_);
}

template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator<=(Proxy<TSpec> const & left_, 
           SimpleType<TValue, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TLeft;
    typedef SimpleType<TValue, TSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}
template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator<=(SimpleType<TValue, TSpec2> const & left_,
           Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TValue, TSpec> TLeft;
    typedef Proxy<TSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

// ----------------------------------------------------------------------------
// Function operator>()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator>(SimpleType<TValue, TSpec> const & left_, 
          TRight const & right_)
{
    SEQAN_CHECKPOINT
            typedef SimpleType<TValue, TSpec> TLeft;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator>(TLeft const & left_, 
          SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TValue, TSpec> TRight;
    typedef typename CompareType<TRight, TLeft>::Type TCompareType;
    return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator>(SimpleType<TLeftValue, TLeftSpec> const & left_, 
          SimpleType<TRightValue, TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
    typedef SimpleType<TRightValue, TRightSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator>(SimpleType<TValue, TSpec> const & left_, 
          SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    return convert<TValue>(left_) > convert<TValue>(right_);
}

template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator>(Proxy<TSpec> const & left_, 
          SimpleType<TValue, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TLeft;
    typedef SimpleType<TValue, TSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}
template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator>(SimpleType<TValue, TSpec2> const & left_,
          Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TValue, TSpec> TLeft;
    typedef Proxy<TSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

// ----------------------------------------------------------------------------
// Function operator>=()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRight>
inline bool
operator>=(SimpleType<TValue, TSpec> const & left_, 
           TRight const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TValue, TSpec> TLeft;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TLeft, typename TValue, typename TSpec>
inline bool
operator>=(TLeft const & left_, 
           SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TValue, TSpec> TRight;
    typedef typename CompareType<TRight, TLeft>::Type TCompareType;
    return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TLeftValue, typename TLeftSpec, typename TRightValue, typename TRightSpec>
inline bool
operator>=(SimpleType<TLeftValue, TLeftSpec> const & left_, 
           SimpleType<TRightValue, TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TLeftValue, TLeftSpec> TLeft;
    typedef SimpleType<TRightValue, TRightSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TValue, typename TSpec>
inline bool
operator>=(SimpleType<TValue, TSpec> const & left_, 
           SimpleType<TValue, TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    return convert<TValue>(left_) >= convert<TValue>(right_);
}

template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator>=(Proxy<TSpec> const & left_, 
           SimpleType<TValue, TSpec2> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TLeft;
    typedef SimpleType<TValue, TSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}
template <typename TSpec, typename TValue, typename TSpec2>
inline bool
operator>=(SimpleType<TValue, TSpec2> const & left_,
           Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TValue, TSpec> TLeft;
    typedef Proxy<TSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

// ----------------------------------------------------------------------------
// Function lexLess()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Remove auto-sequence feature.
template<typename T_, typename TSpec> 
inline bool
lexLess(SimpleType<T_, TSpec> const &_Left, SimpleType<T_, TSpec> const &Right_)
{   // return lexicographical _Left < Right_
    typedef typename MakeUnsigned_<T_>::Type TUnsigned;
    return (TUnsigned)(_Left.value) < (TUnsigned)(Right_.value);
}

// ----------------------------------------------------------------------------
// Function operator++()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline SimpleType<TValue, TSpec> &
operator++(SimpleType<TValue, TSpec> & me)
{
    ++me.value;
    return me;
}

template <typename TValue, typename TSpec>
inline SimpleType<TValue, TSpec>
operator++(SimpleType<TValue, TSpec> & me, int)
{
    SimpleType<TValue, TSpec> dummy = me;
    ++me.value;
    return dummy;
}

// ----------------------------------------------------------------------------
// Function operator--()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline SimpleType<TValue, TSpec> &
operator--(SimpleType<TValue, TSpec> & me)
{
    --me.value;
    return me;
}

template <typename TValue, typename TSpec>
inline SimpleType<TValue, TSpec>
operator--(SimpleType<TValue, TSpec> & me, int)
{
    SimpleType<TValue, TSpec> dummy = me;
    --me.value;
    return dummy;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_ALPHABET_SIMPLE_H_
