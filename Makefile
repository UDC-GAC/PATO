gnu:
	$(MAKE) -f makefiles/Makefile.gnu

clang:
	$(MAKE) -f makefiles/Makefile.clang

intel:
	$(MAKE) -f makefiles/Makefile.intel

clean:
	rm -rf obj target

