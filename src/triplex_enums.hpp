#ifndef TRIPLEX_ENUMS_HPP
#define TRIPLEX_ENUMS_HPP

enum class orientation_t : int
{
    antiparallel = -1,
    both = 0,
    parallel = 1
};

enum class detect_duplicates_t : unsigned int
{
    off = 0,
    permissive = 1,
    strict = 2,
    last = 3
};

enum class error_reference_t : unsigned int
{
    watson_strand = 0,
    purine_strand = 1,
    third_strand = 2,
    last = 3
};

enum class output_format_t : unsigned int
{
    bed = 0,
    triplex = 1,
    summary = 2,
    last = 3
};

enum class run_mode_t : unsigned int
{
    tfo_search = 0,
    tts_search = 1,
    tpx_search = 2
};

#endif
