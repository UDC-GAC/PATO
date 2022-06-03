#ifndef _TRIPLEX_HPP_
#define _TRIPLEX_HPP_

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "seqan.hpp"
#include "options.hpp"
#include "triplex_pattern.hpp"
#include "triplex_alphabet.hpp"

namespace seqan
{

enum orientation
{
    antiparallel = -1,
    both = 0,
    parallel = 1
};

enum duplicates
{
    off = 0,
    permissive = 1,
    strict = 2
};

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
inline bool hasCount(TriplexPotential<TId>& me)
{
    if (me.count_R != 0 || me.count_Y != 0 || me.count_M != 0)
        return true;
    else
        return false;
}

template <typename TId, typename TSize>
inline void setNorm(TriplexPotential<TId>& me,
                    const TSize& seq_length,
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
                    const unsigned int& oligo_length,
                    const unsigned int& duplex_length,
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

template <typename TId>
struct Key<TriplexPotential<TId>>
{
    typedef TId Type;
};

template <typename TSeq, typename TPos>
struct SeqPos
{
    TSeq seqnr;
    TPos position;

    SeqPos(TSeq _seqnr, TPos _pos): seqnr(_seqnr), position(_pos) {}
    SeqPos() {}

    bool operator==(const SeqPos<TSeq, TPos>& b) const;
    bool operator!=(const SeqPos<TSeq, TPos>& b) const;
    bool operator<(const SeqPos<TSeq, TPos>& b) const;
    bool operator>(const SeqPos<TSeq, TPos>& b) const;

    template <typename TSource>
    inline SeqPos&
    operator=(const TSource& source)
    {
        assign(*this, source);
        return *this;
    }
};

template <typename TSeq, typename TPos>
bool SeqPos<TSeq, TPos>::operator==(const SeqPos<TSeq, TPos>& b) const
{
    if (seqnr != b.seqnr)
        return false;
    if (position != b.position)
        return false;
    return true;
}

template <typename TSeq, typename TPos>
bool SeqPos<TSeq, TPos>::operator!=(const SeqPos<TSeq, TPos>& b) const
{
    return !(*this == b);
}

template <typename TSeq, typename TPos>
bool SeqPos<TSeq, TPos>::operator<(const SeqPos<TSeq, TPos>& b) const
{
    if (seqnr < b.seqnr)
        return true;
    if (seqnr > b.seqnr)
        return false;
    if (position < b.position)
        return true;
    return false;
}

template <typename TSeq, typename TPos>
bool SeqPos<TSeq, TPos>::operator>(const SeqPos<TSeq, TPos>& b) const
{
    return b < *this;
}

template <typename TSeq, typename TPos>
inline TSeq getSequenceNo(SeqPos<TSeq, TPos>& me)
{
    return me.seqnr;
}

template <typename TSeq, typename TPos>
inline TSeq getSequenceNo(const SeqPos<TSeq, TPos>& me)
{
    return me.seqnr;
}

template <typename TSeq, typename TPos>
inline TPos getPosition(SeqPos<TSeq, TPos>& me)
{
    return me.position;
}

template <typename TSeq, typename TPos>
inline TPos getPosition(const SeqPos<TSeq, TPos>& me)
{
    return me.position;
}

struct _tts;
typedef Tag<_tts> tts;

struct _tfo;
typedef Tag<_tfo> tfo;

struct _mixed_motif;
typedef Tag<_mixed_motif> mixed_motif;

struct _purine_motif;
typedef Tag<_purine_motif> purine_motif;

struct _pyrimidine_motif;
typedef Tag<_pyrimidine_motif> pyrimidine_motif;

typedef TriplexString TTriplex;
typedef TriplexString TDuplex;
typedef StringSet<TTriplex> TTriplexSet;
typedef StringSet<TDuplex> TDuplexSet;
typedef ModStringTriplex<TTriplex, TTriplex> TOligoMotif;
typedef ModStringTriplex<TDuplex, TDuplex> TDuplexMotif;
typedef StringSet<TOligoMotif> TMotifSet;
typedef StringSet<TDuplexMotif> TDuplexMotifSet;
typedef Graph<Automaton<Triplex, Triplex>> TGraph;
typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
typedef TriplexMatch<Difference<TriplexString>::Type, unsigned int, double> TMatch;
typedef std::vector<TMatch> TMatches;
typedef Pair<unsigned, unsigned> TPotKey;
typedef TriplexPotential<TPotKey> TPotPair;
typedef Map<Pair<TPotKey, TPotPair>, Skiplist<>> TPotentials;
typedef typename Value<TPotentials>::Type TPotential;
typedef typename Cargo<TPotential>::Type TPotCargo;

template <typename TSequenceSet, typename TNameSet>
inline bool load_sequences(TSequenceSet& sequences,
                           TNameSet& names,
                           const char *file_name)
{
    SeqFileIn fasta_file;

    if (!open(fasta_file, file_name)) {
        std::cerr << "PATO: error opening input file '" << file_name << "'\n";
        return false;
    }
    readRecords(names, sequences, fasta_file);

    // crop sequence name
    for (auto& name : names) {
        std::string tmp_name(name.data_begin, length(name));
        size_t num_chars = std::min(tmp_name.find_first_of(' '),
                                    tmp_name.size());
        tmp_name = tmp_name.substr(0, num_chars);
        name = tmp_name;
    }

    return length(sequences) > 0;
}

template <typename TString>
inline void filter_repeats(TString& sequence,
                           unsigned int min_repeat_length,
                           unsigned int max_repeat_period)
{
    typedef Repeat<unsigned int, unsigned int> TRepeat;
    typedef String<TRepeat> TRepeatString;
    typedef typename Iterator<TRepeatString, Rooted>::Type TRepeatIter;

    TRepeatString repeats;
    findRepeats(repeats, sequence, min_repeat_length, max_repeat_period);

    for (TRepeatIter it = begin(repeats, Rooted());
         it != end(repeats, Rooted()); it++) {
        CharString replacement = std::string(it->endPosition
                                             - it->beginPosition, 'N');
        replace(sequence, it->beginPosition, it->endPosition, replacement);
    }
}

template <typename TGraph, typename TString>
inline void make_parser(TGraph& parser,
                        TString& valid_chars,
                        TString& invalid_chars,
                        unsigned int max_interrupts)
{
    typedef typename Iterator<TString, Rooted>::Type TIter;

    TVertexDescriptor root = addVertex(parser);
    assignRoot(parser, root);

    TVertexDescriptor valid_child = addVertex(parser);
    for (TIter valid_it = begin(valid_chars, Rooted());
         valid_it != end(valid_chars, Rooted()); valid_it++) {
        addEdge(parser, root, valid_child, *valid_it);
        addEdge(parser, valid_child, valid_child, *valid_it);
    }

    TVertexDescriptor last_child = valid_child;
    for (unsigned int i = 1; i <= max_interrupts; i++) {
        TVertexDescriptor invalid_child = addVertex(parser);

        for (TIter invalid_it = begin(invalid_chars, Rooted());
             invalid_it != end(invalid_chars, Rooted()); invalid_it++) {
            addEdge(parser, last_child, invalid_child, *invalid_it);
        }

        for (TIter valid_it = begin(valid_chars, Rooted());
             valid_it != end(valid_chars, Rooted()); valid_it++) {
            addEdge(parser, invalid_child, valid_child, *valid_it);
        }

        last_child = invalid_child;
    }
}

template<typename TGraph>
inline void make_tfo_parsers(TGraph& tc_parser,
                             TGraph& ga_parser,
                             TGraph& gt_parser,
                             unsigned int max_interrupts)
{
    TTriplex valid_chars = "TCY";
    TTriplex invalid_chars = "GARN";
    make_parser(tc_parser, valid_chars, invalid_chars, max_interrupts);

    valid_chars = "GAR";
    invalid_chars = "TCYN";
    make_parser(ga_parser, valid_chars, invalid_chars, max_interrupts);

    valid_chars = "GTK";
    invalid_chars = "CAMN";
    make_parser(gt_parser, valid_chars, invalid_chars, max_interrupts);
}

template <typename TGraph, typename TSegString, typename TString>
inline void parse_segments(TGraph& parser,
                           TSegString& segments,
                           TString& sequence,
                           unsigned int max_interrupts,
                           int min_length)
{
    typedef typename Value<TSegString>::Type TSegment;
    typedef typename Iterator<TString, Standard>::Type TIter;

    TIter it = begin(sequence, Standard());
    TIter run_it = begin(sequence, Standard());
    TIter end_it = end(sequence, Standard());

    parseString(parser, getRoot(parser), run_it, end_it);

    while (run_it != end_it) {
        unsigned int shift = std::min(max_interrupts,
                                      static_cast<unsigned int>(run_it - it));
        run_it -= shift;

        unsigned int size = std::max(run_it - run_it, run_it - it);
        if (size >= min_length) {
            TSegment segment = infix(sequence, it, run_it);
            appendValue(segments, segment, Generous());
        }

        run_it += shift + 1;
        it = run_it;

        if (run_it != end_it) {
            parseString(parser, getRoot(parser), run_it, end_it);
        }
    }

    unsigned int size = std::max(end_it - end_it, end_it - it);
    if (size >= min_length) {
        TSegment segment = infix(sequence, it, end_it);
        appendValue(segments, segment, Generous());
    }
}

template <typename TString, typename TChar>
inline std::vector<std::vector<bool>> encode_sequence(TString& motif,
                                                      TChar filter_char,
                                                      TChar interrupt_char,
                                                      std::vector<std::vector<bool>>& block_runs,
                                                      unsigned int min_block_run)
{
    typedef typename Iterator<TString, Standard>::Type TIter;

    std::vector<std::vector<bool>> encoded_seq(3, std::vector<bool>(length(motif), false));

    unsigned int counter = 0;
    unsigned int run_counter = 0;

    for (TIter it = begin(motif, Standard());
         it != end(motif, Standard()); it++, counter++) {
        if (*it == filter_char) {
            encoded_seq[0][counter] = true;
        } else if (*it == interrupt_char) {
            encoded_seq[1][counter] = true;
            encoded_seq[2][counter] = true;

            if (counter - run_counter >= min_block_run) {
                for (unsigned int i = 0; i <= counter - min_block_run; i++) {
                    for (unsigned int j = std::max(i, run_counter) + min_block_run;
                         j <= length(motif); j++) {
                        block_runs[i][j] = true;
                    }
                }
            }

            run_counter = counter + 1;
        } else {
            encoded_seq[2][counter] = true;
        }
    }

    if (counter - run_counter >= min_block_run) {
        for (unsigned int i = 0; i <= counter - min_block_run; i++) {
            for (unsigned int j = std::max(i, run_counter) + min_block_run;
                 j <= length(motif); j++) {
                block_runs[i][j] = true;
            }
        }
    }

    return encoded_seq;
}

template<typename TArray, typename TNum, typename TPos>
inline void increase_right(TArray& encoded_seq,
                           TPos& pos,
                           TNum& filter_chars,
                           TNum& interrupt_chars,
                           TNum& non_filter_chars)
{
    filter_chars += encoded_seq[0][pos];
    interrupt_chars += encoded_seq[1][pos];
    non_filter_chars += encoded_seq[2][pos];

    pos++;
}

template<typename TArray, typename TNum, typename TPos>
inline void increase_left(TArray& encoded_seq,
                          TPos& pos,
                          TNum& filter_chars,
                          TNum& interrupt_chars,
                          TNum& non_filter_chars)
{
    filter_chars -= encoded_seq[0][pos];
    interrupt_chars -= encoded_seq[1][pos];
    non_filter_chars -= encoded_seq[2][pos];

    pos++;
}

template <typename TArray, typename TPos>
inline bool is_interrupt_char(TArray &encoded_seq, TPos pos)
{
    return encoded_seq[1][pos];
}

template <typename TSize, typename TTag>
inline bool motif_specific_constraint(__attribute__((unused)) const TSize& filter_rate,
                                      __attribute__((unused)) const orientation& ornt,
                                      __attribute__((unused)) const TTag& tag,
                                      __attribute__((unused)) const options& opts)
{
    return true;
}

template <typename TSize>
inline bool motif_specific_constraint(const TSize& filter_rate,
                                      const orientation& ornt,
                                      const mixed_motif& tag,
                                      const options& opts)
{
    if (ornt != orientation::parallel
        && filter_rate >= opts.mixed_antiparallel_min_guanine) {
        return true;
    }
    if (ornt != orientation::antiparallel
        && filter_rate >= opts.mixed_parallel_max_guanine) {
        return true;
    }
    return false;
}

template <typename TMotifSet>
inline void add_match(TMotifSet& motif_set,
                      typename Value<TMotifSet>::Type& motif,
                      unsigned int start,
                      unsigned int end,
                      unsigned int errors,
                      const mixed_motif& tag)
{
    add_match(motif_set, motif, start, end, errors, tfo());
}

template <typename TMotifSet>
inline void add_match(TMotifSet& motif_set,
                      typename Value<TMotifSet>::Type& motif,
                      unsigned int start,
                      unsigned int end,
                      unsigned int errors,
                      const purine_motif& tag)
{
    add_match(motif_set, motif, start, end, errors, tfo());
}

template <typename TMotifSet>
inline void add_match(TMotifSet& motif_set,
                      typename Value<TMotifSet>::Type& motif,
                      unsigned int start,
                      unsigned int end,
                      unsigned int errors,
                      const pyrimidine_motif& tag)
{
    add_match(motif_set, motif, start, end, errors, tfo());
}

template <typename TMotifSet>
inline void add_match(TMotifSet& motif_set,
                      typename Value<TMotifSet>::Type& motif,
                      unsigned int start,
                      unsigned int end,
                      unsigned int errors,
                      const tfo& tag)
{
    typedef typename Value<TMotifSet>::Type TPattern;

    auto motif_length = length(motif);

    auto st = isParallel(motif) ? start : motif_length - end;
    auto nd = isParallel(motif) ? end : motif_length - start;

    TPattern sequence_pattern(host(motif),
                              beginPosition(motif) + st,
                              beginPosition(motif) + nd,
                              isParallel(motif),
                              getSequenceNo(motif),
                              isTFO(motif),
                              getMotif(motif));
    setScore(sequence_pattern, end - start - errors);

    appendValue(motif_set, sequence_pattern);
}

template <typename TMotifSet>
inline void add_match(TMotifSet& motif_set,
                      typename Value<TMotifSet>::Type& motif,
                      unsigned int start,
                      unsigned int end,
                      unsigned int errors,
                      const tts& tag)
{
    typedef typename Value<TMotifSet>::Type TPattern;

    auto motif_length = length(motif);

    auto st = getMotif(motif) == '+' ? start : motif_length - end;
    auto nd = getMotif(motif) == '+' ? end : motif_length - start;

    TPattern sequence_pattern(host(motif),
                              beginPosition(motif) + st,
                              beginPosition(motif) + nd,
                              isParallel(motif),
                              getSequenceNo(motif),
                              isTFO(motif),
                              getMotif(motif));
    setScore(sequence_pattern, end - start - errors);

    appendValue(motif_set, sequence_pattern);
}

template <typename TMotifSet>
inline void reduce_motif_set(TMotifSet& output, TMotifSet& input)
{
    typedef typename Value<TMotifSet>::Type TOligoMotif;
    typedef typename Position<TOligoMotif>::Type TIntervalValue;
    typedef Graph<Undirected<void>> TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef IntervalAndCargo<TIntervalValue, TVertexDescriptor>	TInterval;
    typedef Map<Pair<TVertexDescriptor, TOligoMotif *>, Skiplist<>> TPropMap;
    typedef typename Iterator<TMotifSet, Standard>::Type TMotifSetIter;
    typedef IntervalTree<TIntervalValue, TVertexDescriptor> TIntervalTree;
    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef String<TVertexDescriptor> TCounts;
    typedef typename Iterator<TCounts>::Type TCountIter;
    typedef typename Size<TGraph>::Type TSize;
    typedef String<TSize> TComponent;
    typedef Map<Pair<TSize, TOligoMotif *>, Skiplist<>> TCompMap;

    if (length(input) == 0) {
        return;
    }
    if (length(input) == 1) {
        appendValue(output, *begin(input));
        return;
    }

    String<TInterval> intervals;
    resize(intervals, length(input));

    TGraph parser;
    TPropMap prop_map;

    int count = 0;
    for (TMotifSetIter it = begin(input, Standard());
         it != end(input, Standard()); it++, count++) {
        TVertexDescriptor vtx = addVertex(parser);
        insert(prop_map, vtx, &(*it));

        intervals[count].i1 = static_cast<TIntervalValue>(beginPosition(*it));
        intervals[count].i2 = static_cast<TIntervalValue>(endPosition(*it));
        intervals[count].cargo = vtx;
    }

    TIntervalTree tree(intervals);

    TVertexIterator vertex_it(parser);
    while (!atEnd(vertex_it)) {
        TOligoMotif *motif = cargo(prop_map, getValue(vertex_it));

        TCounts tree_result;
        findIntervals(tree_result, tree, beginPosition(*motif),
                      endPosition(*motif));
        for (TCountIter it = begin(tree_result); it != end(tree_result); it++) {
            if (getValue(vertex_it) == *it) {
                continue;
            }
            addEdge(parser, getValue(vertex_it), *it);
        }

        vertex_it++;
    }

    TComponent components;
    TSize num_components = connectedComponents(components, parser);

    TCompMap comp_map;

    goBegin(vertex_it);
    while (!atEnd(vertex_it)) {
        TVertexDescriptor vtx = getValue(vertex_it);

        if (hasKey(comp_map, getProperty(components, vtx))) {
            TOligoMotif *motif(cargo(comp_map, getProperty(components, vtx)));
            merge(*motif, *cargo(prop_map, vtx));
        } else {
            TOligoMotif *motif(cargo(prop_map, vtx));
            insert(comp_map, getProperty(components, vtx), motif);
        }

        vertex_it++;
    }

    for (TSize s = 0; s < num_components; s++) {
        appendValue(output, *cargo(comp_map, s));
    }
}

template <typename TMotifSet, typename TTag>
inline unsigned int filter_with_guanine_and_error_rate(TMotifSet& motif_set,
                                                       typename Value<TMotifSet>::Type& motif,
                                                       char filter_char,
                                                       char interrupt_char,
                                                       bool reduce_set,
                                                       const orientation ornt,
                                                       const TTag& tag,
                                                       const options& opts)
{
    typedef typename Value<TMotifSet>::Type TPattern;
    typedef ModifiedString<TPattern, ModView<FunctorRYFilter>> TFilter;

    TMotifSet tmp_set;
    TMotifSet *ptr_set;
    ptr_set = reduce_set ? &tmp_set : &motif_set;

    auto motif_length = length(motif);

    std::vector<std::vector<bool>> encoded_seq;
    std::vector<std::vector<bool>> block_runs(motif_length - opts.min_block_run + 1,
                                              std::vector<bool>(motif_length + 1, false));
    if (opts.min_guanine_rate <= 0.0) {
        TFilter filter_sequence(motif);
        filter_char = filter_char == 'G' ? 'R' : 'Y';
        encoded_seq = encode_sequence(filter_sequence, filter_char,
                                      interrupt_char, block_runs,
                                      opts.min_block_run);
    } else {
        encoded_seq = encode_sequence(motif, filter_char, interrupt_char,
                                      block_runs, opts.min_block_run);
    }

    double max_error = std::floor(motif_length * opts.error_rate);
    double max_tolerated = std::floor(motif_length * (1.0 - opts.min_guanine_rate));
    if (opts.maximal_error >= 0) {
        max_error = std::min(max_error, static_cast<double>(opts.maximal_error));
    }

    unsigned int filter_chars = 0;
    unsigned int interrupt_chars = 0;
    unsigned int non_filter_chars = 0;

    bool is_match = false;

    unsigned int max_length = motif_length;
    if (opts.max_length >= opts.min_length) {
        max_length = static_cast<unsigned int>(opts.max_length);
    }

    unsigned int tmp_start = 0;
    unsigned int tmp_end = 0;
    unsigned int tmp_errors = 0;
    unsigned int covered_end = 0;

    unsigned int left = 0;
    unsigned int right = 0;

    unsigned int matches = 0;
    while (block_runs[left][motif_length]
           && left + opts.min_length <= motif_length) {
        while (right - left < opts.min_length && right < motif_length) {
            while (right - left < opts.min_length && right < motif_length) {
                increase_right(encoded_seq, right, filter_chars,
                               interrupt_chars, non_filter_chars);
            }
            while (interrupt_chars > max_error) {
                increase_left(encoded_seq, left, filter_chars, interrupt_chars,
                              non_filter_chars);
            }
            while (non_filter_chars > max_tolerated) {
                increase_left(encoded_seq, left, filter_chars, interrupt_chars,
                              non_filter_chars);
            }
            while (left < motif_length
                   && is_interrupt_char(encoded_seq, left)) {
                increase_left(encoded_seq, left, filter_chars, interrupt_chars,
                              non_filter_chars);
            }

            if (right < left) {
                right = left;

                filter_chars = 0;
                interrupt_chars = 0;
                non_filter_chars = 0;
            }
        }

        if (right - left < opts.min_length) {
            break;
        }

        is_match = false;

        while (interrupt_chars <= max_error
               && non_filter_chars <= max_tolerated
               && right - left <= max_length) {
            double filter_chars_rate = static_cast<double>(filter_chars)
                                       / (right - left);
            double interrupt_chars_rate = static_cast<double>(interrupt_chars)
                                          / (right - left);

            if (block_runs[left][right]
                && !is_interrupt_char(encoded_seq, right - 1)
                && interrupt_chars_rate <= opts.error_rate
                && opts.min_guanine_rate <= filter_chars_rate
                && filter_chars_rate <= opts.max_guanine_rate
                && motif_specific_constraint(filter_chars_rate, ornt, tag, opts)) {
                is_match = true;
                matches++;

                tmp_start = left;
                tmp_end = right;
                tmp_errors = interrupt_chars;

                if (opts.all_matches) {
                    covered_end = tmp_end;
                    is_match = false;
                    add_match(*ptr_set, motif, tmp_start, tmp_end, tmp_errors,
                              tag);
                }
            }

            if (right < motif_length) {
                increase_right(encoded_seq, right, filter_chars,
                               interrupt_chars, non_filter_chars);
            } else {
                break;
            }
        }

        if (is_match && tmp_end > covered_end) {
            covered_end = tmp_end;
            add_match(*ptr_set, motif, tmp_start, tmp_end, tmp_errors, tag);
        }

        is_match = false;

        left++;
        while (left < motif_length && is_interrupt_char(encoded_seq, left)) {
            left++;
        }

        right = left;

        filter_chars = 0;
        interrupt_chars = 0;
        non_filter_chars = 0;
    }

    if (reduce_set) {
        reduce_motif_set(motif_set, tmp_set);
    }

    return matches;
}

template <typename TOligoMotifSet, typename TString>
inline void process_tc_motif(TOligoMotifSet& motif_set,
                             TString& sequence,
                             TGraph& parser,
                             unsigned int id,
                             const options& opts)
{
    typedef typename Value<TOligoMotifSet>::Type TTfoMotif;
    typedef typename Infix<TString>::Type TSegment;
    typedef String<TSegment> TSegString;
    typedef typename Iterator<TSegString, Standard>::Type TSegStringIter;

    TSegString segments;
    parse_segments(parser, segments, sequence, opts.max_interruptions,
                   opts.min_length);

    for (TSegStringIter it = begin(segments, Standard());
         it != end(segments, Standard()); it++) {
        TTfoMotif motif(*it, true, id, true, 'Y');
        filter_with_guanine_and_error_rate(motif_set, motif, 'G', 'N', true,
                                           orientation::parallel,
                                           pyrimidine_motif(), opts);
    }
}

template <typename TOligoMotifSet, typename TString>
inline void process_ga_motif(TOligoMotifSet& motif_set,
                             TString& sequence,
                             TGraph& parser,
                             unsigned int id,
                             const options& opts)
{
    typedef typename Value<TOligoMotifSet>::Type TTfoMotif;
    typedef typename Infix<TString>::Type TSegment;
    typedef String<TSegment> TSegString;
    typedef typename Iterator<TSegString, Standard>::Type TSegStringIter;

    TSegString segments;
    parse_segments(parser, segments, sequence, opts.max_interruptions,
                   opts.min_length);

    for (TSegStringIter it = begin(segments, Standard());
         it != end(segments, Standard()); it++) {
        TTfoMotif motif(*it, false, id, true, 'R');
        filter_with_guanine_and_error_rate(motif_set, motif, 'G', 'N', true,
                                           orientation::antiparallel,
                                           purine_motif(), opts);
    }
}

template <typename TOligoMotifSet, typename TString>
inline void process_gt_motif(TOligoMotifSet& motif_set,
                             TString& sequence,
                             TGraph& parser,
                             unsigned int id,
                             const orientation ornt,
                             const options& opts)
{
    typedef typename Value<TOligoMotifSet>::Type TTfoMotif;
    typedef typename Infix<TString>::Type TSegment;
    typedef String<TSegment> TSegString;
    typedef typename Iterator<TSegString, Standard>::Type TSegStringIter;

    TSegString segments;
    parse_segments(parser, segments, sequence, opts.max_interruptions,
                   opts.min_length);

    for (TSegStringIter it = begin(segments, Standard());
         it != end(segments, Standard()); it++) {
        if ((ornt == orientation::both || ornt == orientation::parallel)
            && opts.mixed_parallel_max_guanine >= opts.min_guanine_rate) {
            TTfoMotif motif(*it, true, id, true, 'M');
            filter_with_guanine_and_error_rate(motif_set, motif, 'G', 'N', true,
                                               ornt, mixed_motif(), opts);
        }
        if ((ornt == orientation::both || ornt == orientation::antiparallel)
            && opts.mixed_antiparallel_min_guanine <= opts.max_guanine_rate) {
            TTfoMotif motif(*it, false, id, true, 'M');
            filter_with_guanine_and_error_rate(motif_set, motif, 'G', 'N', true,
                                               ornt, mixed_motif(), opts);
        }
    }
}

template <typename TSpec, typename TTag>
inline void count_duplicates_strict(StringSet<ModStringTriplex<TSpec, TSpec>>& strings,
                                    const TTag& tag,
                                    const options& opts)
{
    typedef ModStringTriplex<TSpec, TSpec> TString;
    typedef StringSet<TString> TStringSet;
    typedef StringSet<TriplexString> TText;
    typedef typename Iterator<TStringSet, Standard>::Type TIter;
    typedef Index<TText> TIndex;
    typedef typename Id<TString>::Type TId;
    typedef typename Position<TString>::Type TPos;
    typedef SeqPos<TId, TPos> TKey;
    typedef Map<TKey, Skiplist<>> TMap;
    typedef typename Iterator<TMap>::Type TMapIter;

    TText string_set;
    for (TIter it = begin(strings, Standard());
         it != end(strings, Standard()); it++) {
        appendValue(string_set, ttsString(*it));
    }

    TIndex index_esa(string_set);
    for (TIter it = begin(strings, Standard());
         it != end(strings, Standard()); it++) {
        int occ = 0;
        TMap loci_map;

        Finder<TIndex> finder_esa(index_esa);
        TKey identity(getSequenceNo(*it), beginPosition(*it));

        while (find(finder_esa, ttsString(*it))) {
            TString ts = value(strings, position(finder_esa).i1);

            if (!opts.same_sequence_duplicates
                && getSequenceNo(*it) == getSequenceNo(ts)) {
                continue;
            }

            if (isTFO(ts) || getMotif(ts) == '+') {
                TKey key(getSequenceNo(ts),
                         beginPosition(ts) + position(finder_esa).i2);

                if (key != identity && !hasKey(loci_map, key)) {
                    add(loci_map, key);
                    occ++;
                }
            } else {
                TKey key(getSequenceNo(ts),
                         endPosition(ts) - position(finder_esa).i2
                         - length(*it));

                if (key != identity && !hasKey(loci_map, key)) {
                    add(loci_map, key);
                    occ++;
                }
            }

            if (opts.duplicate_cutoff >= 0 && occ > opts.duplicate_cutoff) {
                break;
            }
        }

        duplicates(*it, occ);

        if (opts.report_duplicate_locations
            && opts.duplicate_cutoff >= 0
            && occ <= opts.duplicate_cutoff) {
                for (TMapIter map_it = begin(loci_map);
                     map_it != end(loci_map); map_it++) {
                    addDuplicate(*it, getSequenceNo(*map_it),
                                 getPosition(*map_it));
                }
        }
    }
}

template <typename TSpec>
inline void count_duplicates_permissive(StringSet<ModStringTriplex<TSpec, TSpec>>& strings,
                                        const tfo& tag,
                                        const options& opts)
{
    typedef ModStringTriplex<TSpec, TSpec> TString;
    typedef StringSet<TString> TStringSet;
    typedef StringSet<TriplexString> TText;
    typedef typename Iterator<TStringSet>::Type	TIter;
    typedef Index<TText> TIndex;
    typedef typename Id<TString>::Type TId;
    typedef typename Position<TString>::Type TPos;
    typedef SeqPos<TId, TPos> TKey;
    typedef Map<TKey, Skiplist<>> TMap;
    typedef typename Iterator<TMap>::Type TMapIter;

    TText string_set;
    for (TIter it = begin(strings); it != end(strings); it++) {
        appendValue(string_set, tfoString(*it));
    }

    TIndex index_esa(string_set);
    for (TIter it = begin(strings); it != end(strings); it++) {
        int occ = 0;
        TMap loci_map;

        Finder<TIndex> finder_esa(index_esa);
        TKey identity(getSequenceNo(*it), beginPosition(*it));

        while (find(finder_esa, tfoString(*it))) {
            TString ts = value(strings, position(finder_esa).i1);

            if (!opts.same_sequence_duplicates
                && getSequenceNo(*it) == getSequenceNo(ts)) {
                continue;
            }

            TKey key(getSequenceNo(ts),
                     beginPosition(ts) + position(finder_esa).i2);
            if (key != identity && !hasKey(loci_map, key)) {
                add(loci_map, key);
                occ++;
            }

            if (opts.duplicate_cutoff >= 0 && occ > opts.duplicate_cutoff) {
                break;
            }
        }

        duplicates(*it, occ);

        if (opts.report_duplicate_locations
            && opts.duplicate_cutoff >= 0
            && occ <= opts.duplicate_cutoff) {
                for (TMapIter map_it = begin(loci_map);
                     map_it != end(loci_map); map_it++) {
                    addDuplicate(*it, getSequenceNo(*map_it),
                                 getPosition(*map_it));
                }
        }
    }
}

template <typename TMotifSet>
inline void filter_duplicates_with_cutoff(TMotifSet& motif_set, int cutoff)
{
    typedef typename Iterator<TMotifSet, Standard>::Type TIter;

    TMotifSet tmp_set;

    for (TIter it = begin(motif_set, Standard());
         it != end(motif_set, Standard()); it++) {
        if (duplicates(*it) <= cutoff) {
            appendValue(tmp_set, *it);
        }
    }

    if (length(motif_set) - length(tmp_set) == 0) {
        return;
    }

    clear(motif_set);
    for (TIter it = begin(tmp_set, Standard());
         it != end(tmp_set, Standard()); it++) {
        appendValue(motif_set, *it);
    }
}

template <typename TMotifSet, typename TSeqs, typename TNames>
inline bool compute_tfo_motifs(TMotifSet& tfo_motif_set,
                               TSeqs& oligo_sequences,
                               TNames& oligo_names,
                               const options& opts)
{
    typedef typename Iterator<TTriplexSet, Standard>::Type TOligoIter;

    if (!load_sequences(oligo_sequences, oligo_names,
                        toCString(opts.tfo_file))) {
        return false;
    }

    TGraph tc_parser;
    TGraph ga_parser;
    TGraph gt_parser;
    make_tfo_parsers(tc_parser, ga_parser, gt_parser, opts.max_interruptions);

    unsigned int tfo_id = 0;
    for (TOligoIter it = begin(oligo_sequences, Standard());
         it != end(oligo_sequences, Standard()); it++) {
        if (opts.filter_repeats) {
            filter_repeats(*it, opts.min_repeat_length, opts.max_repeat_period);
        }

        if (opts.tc_motif) {
            process_tc_motif(tfo_motif_set, *it, tc_parser, tfo_id, opts);
        }
        if (opts.ga_motif) {
            process_ga_motif(tfo_motif_set, *it, ga_parser, tfo_id, opts);
        }
        if (opts.gt_p_motif && opts.gt_a_motif) {
            process_gt_motif(tfo_motif_set, *it, gt_parser, tfo_id,
                             orientation::both, opts);
        } else if (opts.gt_p_motif) {
            process_gt_motif(tfo_motif_set, *it, gt_parser, tfo_id,
                             orientation::parallel, opts);
        } else if (opts.gt_a_motif) {
            process_gt_motif(tfo_motif_set, *it, gt_parser, tfo_id,
                             orientation::antiparallel, opts);
        }

        tfo_id++;
    }

    if (opts.detect_duplicates != duplicates::off) {
        if (opts.detect_duplicates == duplicates::strict) {
            count_duplicates_strict(tfo_motif_set, tfo(), opts);
        } else {
            count_duplicates_permissive(tfo_motif_set, tfo(), opts);
        }

        if (opts.duplicate_cutoff >= 0) {
            filter_duplicates_with_cutoff(tfo_motif_set, opts.duplicate_cutoff);
        }
    }

    return true;
}

template<
    typename TDuplexMotifSet,
    typename TString,
    typename TGraph,
    typename TId
>
inline void process_tts_motif(TDuplexMotifSet& tts_set,
                              TString& duplex,
                              TGraph& plus_parser,
                              TGraph& minus_parser,
                              const TId& id,
                              const options& opts)
{
    typedef typename Value<TDuplexMotifSet>::Type TTtsMotif;
    typedef typename Infix<TString>::Type TSegment;
    typedef String<TSegment> TSegString;
    typedef typename Iterator<TSegString, Standard>::Type TSegStringIter;

    TSegString plus_segments;
    parse_segments(plus_parser, plus_segments, duplex, opts.max_interruptions,
                   opts.min_length);

    for (TSegStringIter it = begin(plus_segments, Standard());
         it != end(plus_segments, Standard()); it++) {
        TTtsMotif motif(*it, true, id, false, '+');
        filter_with_guanine_and_error_rate(tts_set, motif, 'G', 'Y', true,
                                           orientation::both, tts(), opts);
    }

    TSegString minus_segments;
    parse_segments(minus_parser, minus_segments, duplex, opts.max_interruptions,
                   opts.min_length);

    for (TSegStringIter it = begin(minus_segments, Standard());
         it != end(minus_segments, Standard()); it++) {
        TTtsMotif motif(*it, true, id, false, '-');
        filter_with_guanine_and_error_rate(tts_set, motif, 'G', 'Y', true,
                                           orientation::both, tts(), opts);
    }
}

template<typename TGraph>
inline void make_tts_parsers(TGraph& plus_parser,
                             TGraph& minus_parser,
                             unsigned int max_interrupts)
{
    TDuplex valid_chars = "GAR";
    TDuplex invalid_chars = "TCYN";
    make_parser(plus_parser, valid_chars, invalid_chars, max_interrupts);

    valid_chars = "TCY";
    invalid_chars = "GARN";
    make_parser(minus_parser, valid_chars, invalid_chars, max_interrupts);
}

template<typename TDuplexMotifSet, typename TSeqs, typename TNames>
inline bool compute_tts_motifs(TDuplexMotifSet& tts_motif_set,
                               TSeqs& duplex_sequences,
                               TNames& duplex_names,
                               const options& opts)
{
    typedef typename Iterator<TDuplexSet, Standard>::Type TDuplexIter;

    if (!load_sequences(duplex_sequences, duplex_names,
                        toCString(opts.tts_file))) {
        return false;
    }

    TGraph plus_parser;
    TGraph minus_parser;
    make_tts_parsers(plus_parser, minus_parser, opts.max_interruptions);

    unsigned int tts_id = 0;
    for (TDuplexIter it = begin(duplex_sequences, Standard());
         it != end(duplex_sequences, Standard()); it++) {
        if (opts.filter_repeats) {
            filter_repeats(*it, opts.min_repeat_length, opts.max_repeat_period);
        }

        process_tts_motif(tts_motif_set, *it, plus_parser, minus_parser, tts_id,
                          opts);

        tts_id++;
    }

    return true;
}

template <typename TGraph>
inline void make_triplex_parser(TGraph& triplex_parser,
                                unsigned int max_interrupts)
{
    TDuplex valid_chars = "GAR";
    TDuplex invalid_chars = "TCYN";
    make_parser(triplex_parser, valid_chars, invalid_chars, max_interrupts);
}

template <typename TGraph, typename TMotifSet, typename TTts>
inline unsigned int compute_triplexes(TGraph& parser,
                                      TMotifSet& triplex_set,
                                      TTts& tts_seq,
                                      unsigned int tts_id,
                                      const options& opts)
{
    typedef typename Infix<TTts>::Type TSegment;
    typedef String<TSegment> TSegString;
    typedef typename Iterator<TSegString>::Type TSegStringIter;

    TSegString segments;
    parse_segments(parser, segments, tts_seq, opts.max_interruptions,
                   opts.min_length);

    unsigned int matches = 0;
    for (TSegStringIter it = begin(segments, Standard());
         it != end(segments, Standard()); it++) {
        TDuplexMotif triplex_motif(*it, true, tts_id, false, '+');
        matches += filter_with_guanine_and_error_rate(triplex_set, triplex_motif,
                                                      'G', 'Y', false,
                                                      orientation::both, tts(),
                                                      opts);
    }

    return matches;
}

template <typename TGraph, typename TTfoMotif, typename TTtsMotif>
inline void match_tfo_tts_pair(TMatches& matches,
                               TPotentials& potentials,
                               TGraph& parser,
                               TTfoMotif& tfo_motif,
                               unsigned int tfo_id,
                               TTtsMotif& tts_motif,
                               unsigned int tts_id,
                               int min_score,
                               const options& opts)
{
    typedef typename Iterator<TDuplexMotifSet>::Type TTtsIter;
    typedef typename Value<TDuplexMotifSet>::Type TTtsValue;
    typedef typename Host<TTtsValue>::Type TTts;
    typedef typename Value<TMotifSet>::Type TTfoValue;
    typedef typename Host<TTfoValue>::Type TTfo;

    auto tfo_candidate = ttsString(tfo_motif);
    auto tts_candidate = ttsString(tts_motif);

    int tfo_length = length(tfo_candidate);
    int tts_length = length(tts_candidate);

    for (int diag = -(tts_length - opts.min_length);
         diag <= tfo_length - opts.min_length; diag++) {
        int tfo_offset = 0;
        int tts_offset = 0;

        if (diag < 0) {
            tts_offset = -diag;
        } else {
            tfo_offset = diag;
        }

        int match_length = std::min(tts_length - tts_offset,
                                    tfo_length - tfo_offset);
        TTts tmp_tts(infix(tts_candidate, tts_offset,
                           tts_offset + match_length));
        TTfo tmp_tfo(infix(tfo_candidate, tfo_offset,
                           tfo_offset + match_length));

        int m = 0;
        for (int i = 0; i < match_length; i++) {
            if (tmp_tts[i] != tmp_tfo[i]) {
                tmp_tts[i] = 'N';
            } else {
                m++;
            }
        }

        if (m < min_score) {
            continue;
        }

        TDuplexMotifSet triplex_set;
        unsigned int total = compute_triplexes(parser, triplex_set, tmp_tts,
                                               getSequenceNo(tts_motif), opts);

        if (total == 0) {
            continue;
        }

        size_t tfo_start;
        size_t tfo_end;
        size_t tts_start;
        size_t tts_end;

        char strand;

        for (TTtsIter it = begin(triplex_set); it != end(triplex_set); it++) {
            int score = 0;
            int guanines = 0;

            for (unsigned int i = beginPosition(*it);
                 i < endPosition(*it); i++) {
                if (tmp_tts[i] != 'N') {
                    score++;
                }
                if (tmp_tts[i] == 'G') {
                    guanines++;
                }
            }

            if (isParallel(tfo_motif)) {
                tfo_start = tfo_offset + beginPosition(tfo_motif)
                            + beginPosition(*it);
                tfo_end = tfo_start + length(*it);
            } else {
                tfo_end = endPosition(tfo_motif)
                          - (tfo_offset + beginPosition(*it));
                tfo_start = tfo_end - length(*it);
            }

            if (getMotif(tts_motif) == '+') {
                tts_start = tts_offset + beginPosition(tts_motif)
                            + beginPosition(*it);
                tts_end = tts_start + length(*it);
                strand = '+';
            } else {
                tts_end = endPosition(tts_motif)
                          - (tts_offset + beginPosition(*it));
                tts_start = tts_end - length(*it);
                strand = '-';
            }

            TMatch match(tfo_id, tfo_start, tfo_end, getSequenceNo(tts_motif),
                         tts_id, tts_start, tts_end, score,
                         isParallel(tfo_motif), getMotif(tfo_motif), strand,
                         guanines);
            matches.push_back(std::move(match));
        }

        TPotKey key(getSequenceNo(tfo_motif), getSequenceNo(tts_motif));
        if (hasKey(potentials, key)) {
            TPotCargo *potential = &cargo(potentials, key);
            addCount(*potential, total, getMotif(tfo_motif));
        } else {
            TPotCargo potential(key);
            addCount(potential, total, getMotif(tfo_motif));
            setNorm(potential, length(host(tfo_motif)), length(host(tts_motif)),
                    opts);
            insert(potentials, TPotential(key, potential));
        }
    }
}

template <typename TOligoMotifSet, typename TDuplexMotifSet, typename TNames>
inline void match_tfo_tts_motifs(TMatches &matches,
                                 TPotentials &potentials,
                                 TOligoMotifSet& tfo_motif_set,
                                 TNames& tfo_names,
                                 TDuplexMotifSet& tts_motif_set,
                                 TNames& tts_names,
                                 const options& opts)
{
    typedef typename Iterator<TOligoMotifSet>::Type TOligoMotifIter;
    typedef typename Iterator<TDuplexMotifSet>::Type TDuplexMotifIter;

    TGraph triplex_parser;
    make_triplex_parser(triplex_parser, opts.max_interruptions);

    int min_score = opts.min_length
                    - static_cast<int>(std::ceil(opts.error_rate
                                                 * opts.min_length));

    unsigned int tfo_id = 0;
    for (TOligoMotifIter tfo_it = begin(tfo_motif_set);
         tfo_it != end(tfo_motif_set); tfo_it++, tfo_id++) {
        unsigned int tts_id = 0;
        for (TDuplexMotifIter tts_it = begin(tts_motif_set);
             tts_it != end(tts_motif_set); tts_it++, tts_id++) {
            match_tfo_tts_pair(matches, potentials, triplex_parser,
                               *tfo_it, tfo_id, *tts_it, tts_id, min_score,
                               opts);
        }
    }
}

template <typename TOligoMotifSet, typename TDuplexMotifSet>
inline CharString error_string(TMatch &match,
                               TOligoMotifSet& tfo_motif_set,
                               TDuplexMotifSet& tts_motif_set,
                               const options& opts)
{
    typedef typename Value<TOligoMotifSet>::Type TTfoMotif;
    typedef typename Value<TDuplexMotifSet>::Type TTtsMotif;

    std::ostringstream errors;

    TTfoMotif tfo_motif(host(value(tfo_motif_set, match.tfoNo)), match.oBegin,
                        match.oEnd, match.parallel,
                        value(tfo_motif_set, match.tfoNo).seqNo, true,
                        match.motif);
    TTtsMotif tts_motif(host(value(tts_motif_set, match.ttsNo)), match.dBegin,
                        match.dEnd, match.parallel, match.ttsSeqNo, false,
                        match.strand);

    CharString tfo_ps = prettyString(tfo_motif);
    CharString tts_ps = prettyString(tts_motif);

    if (match.strand == '-') {
        reverse(tts_ps);

        if (value(tfo_motif_set, match.tfoNo).parallel) {
            reverse(tfo_ps);
        }

        auto tts_it = begin(tts_motif);
        auto tts_end = end(tts_motif);
        auto tfo_it = begin(tfo_motif);
        auto tfo_end = end(tfo_motif);

        unsigned int i = 0;
        while (tts_end != tts_it && tfo_end != tfo_it) {
            tts_end--;
            tfo_end--;

            if (*tts_end != *tfo_end) {
                if (!isupper(value(tts_ps, i)) && !isupper(value(tfo_ps, i))) {
                    errors << "b" << i;
                } else if (!isupper(value(tts_ps, i))) {
                    errors << "d" << i;
                } else if (!isupper(value(tfo_ps, i))) {
                    errors << "o" << i;
                } else {
                    errors << "t" << i;
                }
            }
            i++;
        }
    } else {
        if (!value(tfo_motif_set, match.tfoNo).parallel) {
            reverse(tfo_ps);
        }

        auto tts_it = begin(tts_motif);
        auto tts_end = end(tts_motif);
        auto tfo_it = begin(tfo_motif);
        auto tfo_end = end(tfo_motif);

        unsigned int i = 0;
        while (tts_end != tts_it && tfo_end != tfo_it) {
            if (*tts_it != *tfo_it) {
                if (!isupper(value(tts_ps, i)) && !isupper(value(tfo_ps, i))) {
                    errors << "b" << i;
                } else if (!isupper(value(tts_ps, i))) {
                    errors << "d" << i;
                } else if (!isupper(value(tfo_ps, i))) {
                    errors << "o" << i;
                } else {
                    errors << "t" << i;
                }
            }
            tts_it++;
            tfo_it++;
            i++;
        }
    }

    return errors.str();
}

template <typename TNames>
inline void print_summary(TPotentials& potentials,
                          TNames& oligo_names,
                          TNames& duplex_names,
                          const options& opts)
{
    typedef typename Iterator<TPotentials>::Type TPotIter;

    CharString output_file_name;
    append(output_file_name, opts.output_file);
    append(output_file_name, ".summary");

    std::ofstream output_file(toCString(output_file_name), std::ios_base::out);
    if (!output_file) {
        std::cerr << "PATO: error opening output file '"
                  << toCString(output_file_name) << "'\n";
        return;
    }

    output_file << "# Duplex-ID\tSequence-ID\tTotal (abs)\tTotal (rel)\tGA (abs"
                   ")\tGA (rel)\tTC (abs)\tTC (rel)\tGT (abs)\tGT (rel)\n";

    for (TPotIter it = begin(potentials); it != end(potentials); it++) {
        TPotCargo pot = cargo(*it);
        if (hasCount(pot)) {
            output_file << duplex_names[getKey(pot).i2] << "\t"
                        << oligo_names[getKey(pot).i1] << "\t"
                        << getCounts(pot) << "\t"
                        << std::setprecision(3) << getCounts(pot) / getNorm(pot) << "\t"
                        << getCount(pot, 'R') << "\t"
                        << std::setprecision(3) << getCount(pot, 'R') / getNorm(pot) << "\t"
                        << getCount(pot, 'Y') << "\t"
                        << std::setprecision(3) << getCount(pot, 'Y') / getNorm(pot) << "\t"
                        << getCount(pot, 'M') << "\t"
                        << std::setprecision(3) << getCount(pot, 'M') / getNorm(pot) << "\t"
                        << "\n";
        }
    }
}

template <typename TOligoMotifSet, typename TDuplexMotifSet, typename TNames>
inline void print_tfo_tts_pairs(TMatches &matches,
                                TOligoMotifSet& tfo_motif_set,
                                TNames& tfo_names,
                                TDuplexMotifSet& tts_motif_set,
                                TNames& tts_names,
                                const options& opts)
{
    CharString output_file_name;
    append(output_file_name, opts.output_file);
    append(output_file_name, ".out");

    std::ofstream output_file(toCString(output_file_name), std::ios_base::out);
    if (!output_file) {
        std::cerr << "PATO: error opening output file '"
                  << toCString(output_file_name) << "'\n";
        return;
    }

    output_file << "# Sequence-ID\tTFO start\tTFO end\tDuplex-ID\tTTS start\tTT"
                   "S end\tScore\tError-rate\tErrors\tMotif\tStrand\tOrientatio"
                   "n\tGuanine-rate\n";

    for (auto& match : matches) {
        auto tfo_seq_id = value(tfo_motif_set, match.tfoNo).seqNo;
        auto tts_seq_id = match.ttsSeqNo;

        output_file << tfo_names[tfo_seq_id] << "\t"
                    << match.oBegin << "\t"
                    << match.oEnd << "\t"
                    << tts_names[tts_seq_id] << "\t"
                    << match.dBegin << "\t"
                    << match.dEnd << "\t"
                    << match.mScore << "\t"
                    << std::setprecision(2) << 1.0 - match.mScore / (match.dEnd - match.dBegin) << "\t"
                    << error_string(match, tfo_motif_set, tts_motif_set, opts) << "\t"
                    << match.motif << "\t"
                    << match.strand << "\t"
                    << (match.parallel ? 'P' : 'A') << "\t"
                    << match.guanines / (match.dEnd - match.dBegin)
                    << "\n";
    }
}

inline void search_triplexes(const options& opts)
{
    TMotifSet tfo_motif_set;
    TTriplexSet oligo_sequences;
    StringSet<CharString> oligo_names;
    if (!compute_tfo_motifs(tfo_motif_set, oligo_sequences, oligo_names,
                            opts)) {
        return;
    }

    TDuplexMotifSet tts_motif_set;
    TDuplexSet duplex_sequences;
    StringSet<CharString> duplex_names;
    if (!compute_tts_motifs(tts_motif_set, duplex_sequences, duplex_names,
                            opts)) {
        return;
    }

    TMatches matches;
    TPotentials potentials;
    match_tfo_tts_motifs(matches, potentials, tfo_motif_set, oligo_names,
                         tts_motif_set, duplex_names, opts);

    print_tfo_tts_pairs(matches, tfo_motif_set, oligo_names, tts_motif_set,
                        duplex_names, opts);
    print_summary(potentials, oligo_names, duplex_names, opts);
}

}

#endif
