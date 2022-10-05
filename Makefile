gnu:
	$(MAKE) -f Makefile.gnu

clang:
	$(MAKE) -f Makefile.clang

intel:
	$(MAKE) -f Makefile.intel

clean:
	rm -rf obj target

