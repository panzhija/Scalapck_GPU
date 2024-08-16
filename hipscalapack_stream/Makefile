include SLmake_hipBLAS.inc

.PHONY: all clean cleanall

all:
	$(MAKE) -C PBLAS/SRC double complex16

clean:
	$(MAKE) -C PBLAS/SRC clean

cleanall: clean
	-rm -f $(SCALAPACKLIB)
