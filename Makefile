all:
	cd Examples; make all
	cd Converters; make all

doc:
	cd Doxygen; doxygen Doxyfile

clean:
	rm */*.exe */*.o */*~
