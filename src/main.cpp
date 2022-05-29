#include "options.hpp"
#include "triplex.hpp"
#include "command_line_parser.hpp"

int main(int argc, char *argv[])
{
    options opts;
    if (!parse_command_line(opts, argc, (const char **) argv)) {
        return 0;
    }

    seqan::search_triplexes(opts);

    return 0;
}
