gnu:
	$(MAKE) -f makefiles/Makefile.gnu

clang:
	$(MAKE) -f makefiles/Makefile.clang

intel:
	$(MAKE) -f makefiles/Makefile.intel

c: clean
clean:
	@sed -i `[ $(shell uname -s) = "Darwin" ] && echo "-E"` "s/seqan::setCompilationOpts(parser,.*);/seqan::setCompilationOpts(parser, \"COMPILATIONOPTS_PLACEHOLDER\");/g" src/command_line_parser.hpp
	@[ $(shell uname -s) = "Darwin" ] && rm -rf src/command_line_parser.hpp-E || :
	@rm -rf obj target
