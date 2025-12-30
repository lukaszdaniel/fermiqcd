all:
	cd Libraries; make all
	cd Examples; make all
	cd Converters; make all
	cd SciDac2007; make all
	cd Vegas; make all
	cd Junk; make all
	cd tests; make all
	cd tests; make check

doc:
	cd Doxygen; doxygen Doxyfile

clean:
	cd Libraries; make clean
	cd Examples; make clean
	cd Converters; make clean
	cd SciDac2007; make clean
	cd Vegas; make clean
	cd Junk; make clean
	cd tests; make clean
