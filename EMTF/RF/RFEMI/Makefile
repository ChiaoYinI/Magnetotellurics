FC = f77
FFLAGS= -g
BIN_DIR=	$(HOME)/bin
#BIN_DIR=	../../test/bin
TARGETS=	rfemi install clean

rfemi:	rfemi.o inputEMI.o ../out_bin.f mksp.o rdblk.o outchk.o ../../include/datsz.inc ../../include/nchmx.inc
	$(FC) $(FFLAGS) -o rfemi rfemi.o inputEMI.o ../out_bin.f mksp.o rdblk.o outchk.o

rdblk.o:	rdblk.c
	cc -c -g rdblk.c

clean:
	rm -f *.o

install:	
		\mv rfemi $(BIN_DIR)/.
