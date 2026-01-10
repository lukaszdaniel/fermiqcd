include ../Makefile.common


check:
	for file in `ls *.exe` ; do \
       echo "Testing" $${file:-3} "..."; \
	   ./$${file} > $${file%.*}".out"; \
	   diff -I "^PSIM LOG" -I "^\* Process" -I "^PROCESS " -I "... Loading time" -I "... Saving time" -I "<time>" $${file%.*}".out" $${file%.*}".save"; \
	   echo "... Done"; \
	   echo; \
	done
