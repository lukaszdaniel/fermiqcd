all:
	cd Examples; make all
	cd Converters; make all
	cd SciDac2007; make all
	cd Vegas; make all

doc:
	cd Doxygen; doxygen Doxyfile

clean:
	rm */*.exe */*.o */*~
