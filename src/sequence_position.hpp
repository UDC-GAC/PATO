#ifndef SEQUENCE_POSITION_HPP
#define SEQUENCE_POSITION_HPP

namespace seqan
{

template <typename TSeq, typename TPos>
struct SeqPos
{
    TSeq seqnr;
    TPos position;

    bool operator==(const SeqPos<TSeq, TPos>& b) const;
    bool operator!=(const SeqPos<TSeq, TPos>& b) const;
    bool operator<(const SeqPos<TSeq, TPos>& b) const;
    bool operator>(const SeqPos<TSeq, TPos>& b) const;

    SeqPos(TSeq _seqnr, TPos _pos) : seqnr(_seqnr), position(_pos)
    {}

    SeqPos()
    {}

    template <typename TSource>
    inline SeqPos &
    operator = (TSource const & source) {
        assign(*this, source);
        return *this;
    }
};

template <typename TSeq, typename TPos>
bool SeqPos<TSeq, TPos>::operator==(const SeqPos<TSeq, TPos>& b) const {
    if (seqnr != b.seqnr) return false;
    if (position != b.position) return false;
    return true;
}

template <typename TSeq, typename TPos>	
bool SeqPos<TSeq, TPos>::operator!=(const SeqPos<TSeq, TPos>& b) const {
    return !(*this == b);
}

template <typename TSeq, typename TPos>
bool SeqPos<TSeq, TPos>::operator<(const SeqPos<TSeq, TPos>& b) const {
    if (seqnr < b.seqnr) return true;
    if (seqnr > b.seqnr) return false;
    if (position < b.position) return true;
    return false;
}

template <typename TSeq, typename TPos>
bool SeqPos<TSeq, TPos>::operator>(const SeqPos<TSeq, TPos>& b) const {
    return b<*this;
}

template <typename TSeq, typename TPos>
inline TSeq getSequenceNo(SeqPos<TSeq, TPos> & me){
    return me.seqnr;
}

template <typename TSeq, typename TPos>
inline TSeq getSequenceNo(SeqPos<TSeq, TPos> const & me){
    return me.seqnr;
}

template <typename TSeq, typename TPos>
inline TPos getPosition(SeqPos<TSeq, TPos> & me){
    return me.position;
}

template <typename TSeq, typename TPos>
inline TPos getPosition(SeqPos<TSeq, TPos> const & me){
    return me.position;
}

}

#endif
