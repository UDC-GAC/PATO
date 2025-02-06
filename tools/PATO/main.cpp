// MIT License
//
// Copyright (c) 2022-onwards Iñaki Amatria-Barral
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <PATO/command_line_parser.h>
#include <PATO/tfo_finder.h>
#include <PATO/tpx_finder.h>
#include <PATO/tts_finder.h>

#include <iostream>

namespace {

template <typename... Tys> struct visitors_t : Tys... {
  using Tys::operator()...;
};
template <typename... Tys> visitors_t(Tys...) -> visitors_t<Tys...>;

} // namespace

int main(int argc, char *argv[]) {
  pato::parse_result_t result =
      pato::parse_command_line(argc, argv, std::cout, std::cerr);
  return std::visit(
      visitors_t{
          [](int return_code) { return return_code; },
          [](const pato::options_t &opts) {
            bool success = false;
            switch (opts.run_mode) {
            case pato::run_mode_t::tfo_search: {
              pato::find_tfo_motifs_result result = pato::find_tfo_motifs(opts);
              success =
                  pato::handle_find_tfo_motifs_result(result, std::cerr, opts);
              break;
            }
            case pato::run_mode_t::tts_search: {
              pato::find_tts_motifs_result result = pato::find_tts_motifs(opts);
              success =
                  pato::handle_find_tts_motifs_result(result, std::cerr, opts);
              break;
            }
            case pato::run_mode_t::tpx_search: {
              pato::find_tpx_result result = pato::find_tpxes(opts);
              success = pato::handle_find_tpx_result(result, std::cerr, opts);
              break;
            }
            }
            return success ? 0 : 1;
          }},
      result);
}
