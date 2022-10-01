#include "options.hpp"
#include "tfo_finder.hpp"
#include "tts_finder.hpp"
#include "triplex_enums.hpp"
#include "triplex_finder.hpp"
#include "command_line_parser.hpp"

int main(int argc, char *argv[])
{
    options opts;
    if (!parse_command_line(opts, argc, argv)) {
        return 0;
    }

    switch (opts.run_mode) {
        case run_mode_t::tfo_search:
            find_tfo_motifs(opts);
            break;
        case run_mode_t::tts_search:
            find_tts_motifs(opts);
            break;
        case run_mode_t::tpx_search:
            find_triplexes(opts);
            break;
    }

    return 0;
}
