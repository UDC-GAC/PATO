#include "guanine_filter.hpp"

void increase_right(char_set_set_t& encoded_seq,
                    unsigned int& pos,
                    unsigned int& filter_chars,
                    unsigned int& interrupt_chars,
                    unsigned int& non_filter_chars)
{
    filter_chars += encoded_seq[0][pos];
    interrupt_chars += encoded_seq[1][pos];
    non_filter_chars += encoded_seq[2][pos];

    pos++;
}

void increase_left(char_set_set_t& encoded_seq,
                   unsigned int& pos,
                   unsigned int& filter_chars,
                   unsigned int& interrupt_chars,
                   unsigned int& non_filter_chars)
{
    filter_chars -= encoded_seq[0][pos];
    interrupt_chars -= encoded_seq[1][pos];
    non_filter_chars -= encoded_seq[2][pos];

    pos++;
}

bool is_interrupt_char(char_set_set_t& encoded_seq, unsigned int pos)
{
    return encoded_seq[1][pos];
}

bool motif_specific_constraint(double filter_rate,
                               orientation ornt,
                               const mixed_motif_t& tag,
                               const options& opts)
{
    if (ornt != orientation::parallel
        && filter_rate >= opts.mixed_antiparallel_min_guanine) {
        return true;
    } else if (ornt != orientation::antiparallel
               && filter_rate <= opts.mixed_parallel_max_guanine) {
        return true;
    }
    return false;
}

void add_match(motif_set_t& motifs,
               motif_t& motif,
               unsigned int start,
               unsigned int end,
               unsigned int errors,
               const tfo_t& tag)
{
    auto motif_length = seqan::length(motif);
    auto st = seqan::isParallel(motif) ? start : motif_length - end;
    auto nd = seqan::isParallel(motif) ? end : motif_length - start;

    motif_t tmp_motif(seqan::host(motif),
                      seqan::beginPosition(motif) + st,
                      seqan::beginPosition(motif) + nd,
                      seqan::isParallel(motif),
                      seqan::getSequenceNo(motif),
                      seqan::isTFO(motif),
                      seqan::getMotif(motif));
    seqan::setScore(tmp_motif, end - start - errors);

    motifs.push_back(tmp_motif);
}

void add_match(motif_set_t& motifs,
               motif_t& motif,
               unsigned int start,
               unsigned int end,
               unsigned int errors,
               const mixed_motif_t& tag)
{
    add_match(motifs, motif, start, end, errors, tfo_t());
}

void add_match(motif_set_t& motifs,
               motif_t& motif,
               unsigned int start,
               unsigned int end,
               unsigned int errors,
               const purine_motif_t& tag)
{
    add_match(motifs, motif, start, end, errors, tfo_t());
}

void add_match(motif_set_t& motifs,
               motif_t& motif,
               unsigned int start,
               unsigned int end,
               unsigned int errors,
               const pyrimidine_motif_t& tag)
{
    add_match(motifs, motif, start, end, errors, tfo_t());
}

void add_match(motif_set_t& motifs,
               motif_t& motif,
               unsigned int start,
               unsigned int end,
               unsigned int errors,
               const tts_t& tag)
{
    auto motif_length = seqan::length(motif);
    auto st = seqan::getMotif(motif) == '+' ? start : motif_length - end;
    auto nd = seqan::getMotif(motif) == '+' ? end : motif_length - start;

    motif_t tmp_motif(seqan::host(motif),
                      seqan::beginPosition(motif) + st,
                      seqan::beginPosition(motif) + nd,
                      seqan::isParallel(motif),
                      seqan::getSequenceNo(motif),
                      seqan::isTFO(motif),
                      seqan::getMotif(motif));
    seqan::setScore(tmp_motif, end - start - errors);

    motifs.push_back(tmp_motif);
}
