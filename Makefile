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
	rm */*.exe */*.o */*~
