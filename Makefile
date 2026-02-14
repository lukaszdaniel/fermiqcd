.PHONY: all clean doc

all:
	$(MAKE) -C Libraries
	$(MAKE) -C progs
	$(MAKE) -C Examples
	$(MAKE) -C Converters
	$(MAKE) -C SciDac2007
	$(MAKE) -C Vegas
	$(MAKE) -C Junk
	$(MAKE) -C tests
	$(MAKE) -C tests check

doc:
	cd Doxygen && doxygen Doxyfile

clean:
	$(MAKE) -C Libraries clean
	$(MAKE) -C progs clean
	$(MAKE) -C Examples clean
	$(MAKE) -C Converters clean
	$(MAKE) -C SciDac2007 clean
	$(MAKE) -C Vegas clean
	$(MAKE) -C Junk clean
	$(MAKE) -C tests clean
