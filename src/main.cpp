#include "options.hpp"
#include "tfo_finder.hpp"
#include "tts_finder.hpp"
#include "triplex_finder.hpp"
#include "command_line_parser.hpp"

int main(int argc, char *argv[])
{
    options opts;
    if (!parse_command_line(opts, argc, argv)) {
        return 0;
    }

    switch (opts.run_mode) {
        case 0:
            find_tfo_motifs(opts);
            break;
        case 1:
            find_tts_motifs(opts);
            break;
        case 2:
            find_triplexes(opts);
            break;
    }

    return 0;
}
