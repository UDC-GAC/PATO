#ifndef TRIPLEX_MATCH_HPP
#define TRIPLEX_MATCH_HPP

#include <seqan/basic.h>

#include "options.hpp"

namespace seqan
{

template <typename _TGPos, typename TSize, typename TScore>
struct TriplexMatch
{
    typedef typename MakeSigned_<_TGPos>::Type TGPos;

    TSize tfoNo;
    TGPos oBegin;
    TGPos oEnd;
    TSize ttsSeqNo;
    TSize ttsNo;
    TGPos dBegin;
    TGPos dEnd;
    TScore mScore;
    bool parallel;
    char motif;
    char strand;
    TScore guanines;

    TriplexMatch()
    {
        tfoNo = -1;
        oBegin = -1;
        oEnd = -1;
        ttsSeqNo = 0;
        ttsNo = -1;
        dBegin = -1;
        dEnd = -1;
        mScore = 0;
        parallel = false;
        motif = ' ';
        strand = ' ';
        guanines = 0;
    }

    TriplexMatch(TSize _tfoNo,
                 TGPos _oBegin,
                 TGPos _oEnd,
                 TSize _ttsSeqNo,
                 TSize _ttsNo,
                 TGPos _dBegin,
                 TGPos _dEnd,
                 TScore _mScore,
                 bool _parallel,
                 char _motif,
                 char _strand,
                 TScore _guanines):
        tfoNo(_tfoNo),
        oBegin(_oBegin),
        oEnd(_oEnd),
        ttsSeqNo(_ttsSeqNo),
        ttsNo(_ttsNo),
        dBegin(_dBegin),
        dEnd(_dEnd),
        mScore(_mScore),
        parallel(_parallel),
        motif(_motif),
        strand(_strand),
        guanines(_guanines)
    {}

    template <typename TSource>
    inline TriplexMatch&
    operator=(const TSource& source) {
        assign(*this, source);
        return *this;
    }
};

template <typename TId>
struct TriplexPotential
{
    typedef unsigned int TCount;
    typedef double TNorm;

    TId key;
    TCount count_R;
    TCount count_M;
    TCount count_Y;
    TNorm norm;

    TriplexPotential(TId _key):
        key(_key)
    {
        count_R = 0;
        count_M = 0;
        count_Y = 0;
        norm = 0.;
    }

    TriplexPotential()
    {
        count_R = 0;
        count_M = 0;
        count_Y = 0;
        norm = 0.;
    }

    bool operator==(const TriplexPotential<TId>& b) const;
    bool operator!=(const TriplexPotential<TId>& b) const;
    bool operator<(const TriplexPotential<TId>& b) const;
    bool operator>(const TriplexPotential<TId>& b) const;

    template <typename TSource>
    inline TriplexPotential&
    operator=(const TSource& source) {
        assign(*this, source);
        return *this;
    }
};

template <typename TId>
bool TriplexPotential<TId>::operator==(const TriplexPotential<TId>& b) const {
    if (key != b.key)
        return false;
    return true;
}

template <typename TId>
bool TriplexPotential<TId>::operator!=(const TriplexPotential<TId>& b) const {
    return !(*this == b);
}

template <typename TId>
bool TriplexPotential<TId>::operator<(const TriplexPotential<TId>& b) const {
    if (key < b.key)
        return true;
    return false;
}

template <typename TId>
bool TriplexPotential<TId>::operator>(const TriplexPotential<TId>& b) const {
    return b < *this;
}

template <typename TId>
inline unsigned getCount(TriplexPotential<TId>& me, char motif)
{
    if (motif == 'R' || motif == '+')
        return me.count_R;
    else if (motif == 'Y' || motif == '-')
        return me.count_Y;
    else if (motif == 'M')
        return me.count_M;
    else
        return me.count_R + me.count_Y + me.count_M;
}

template <typename TId>
inline unsigned getCounts(TriplexPotential<TId>& me)
{
    return me.count_R + me.count_Y + me.count_M;
}

template <typename TId>
inline double getNorm(TriplexPotential<TId>& me)
{
    return me.norm;
}

template <typename TId>
inline TId getKey(TriplexPotential<TId>& me)
{
    return me.key;
}

template <typename TId>
inline void addCount(TriplexPotential<TId>& me, unsigned int count, char motif)
{
    if (motif == 'R' || motif == '+')
        me.count_R += count;
    else if (motif == 'Y' || motif == '-')
        me.count_Y += count;
    else if (motif == 'M')
        me.count_M += count;
}

template <typename TId>
inline void mergeCount(TriplexPotential<TId>& me, TriplexPotential<TId>& other)
{
    me.count_R += other.count_R;
    me.count_Y += other.count_Y;
    me.count_M += other.count_M;
}

template <typename TId>
inline bool hasCount(TriplexPotential<TId>& me)
{
    if (me.count_R != 0 || me.count_Y != 0 || me.count_M != 0)
        return true;
    else
        return false;
}

template <typename TId, typename TSize>
inline void setNorm(TriplexPotential<TId>& me,
                    const TSize seq_length,
                    const options &opts)
{
    TSize max_len = opts.max_length;
    if (opts.max_length < opts.min_length || seq_length < max_len) {
        max_len = seq_length;
    }

    double norm = 0.0;
    for (TSize i = opts.min_length; i <= max_len; i++) {
        norm += static_cast<double>(seq_length - i + 1);
    }

    me.norm = norm;
}

template <typename TId>
inline void setNorm(TriplexPotential<TId>& me,
                    const unsigned int oligo_length,
                    const unsigned int duplex_length,
                    const options& opts)
{
    unsigned int max_len = opts.max_length;
    if (opts.max_length < opts.min_length
        || oligo_length < max_len
        || duplex_length < max_len) {
        max_len = std::min(duplex_length, oligo_length);
    }

    double norm = 0.0;
    for (unsigned int i = opts.min_length; i <= max_len; i++) {
        norm += (duplex_length - i + 1) * (oligo_length - i + 1);
    }
    norm *= 2.0;

    me.norm = norm;
}

}

#endif
