#include "options.hpp"
#include "triplex_finder.hpp"
#include "command_line_parser.hpp"

int main(int argc, char *argv[])
{
    options opts;
    if (!parse_command_line(opts, argc, argv)) {
        return 0;
    }

    find_triplexes(opts);

    return 0;
}
