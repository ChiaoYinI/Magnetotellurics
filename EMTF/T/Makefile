#	Makefile for TRANMT
# Did not compile succesfully. 
#
FC = gfortran 
FFLAGS= -g
#BIN_DIR= /home/server/homes/pi/eisel/bin
#BIN_DIR= ${MTbin}
BIN_DIR= ../examples/bin 

TARGETS=	tranmtLR install clean

OBJS=		bset.o chdec.o cmove.o edfwts.o matmult.o mkfdir.o \
		mkhat.o mkrec.o mkrhat.o modobs.o movec.o mult.o \
		rbstrn.o rdcndwt.o rfhead.o mcohc.o rbstSdm.o \
		savtf.o setup.o sindex.o sort.o stack.o trlsrr.o \
		unstack.o wt.o sdmsort.o cohsrt.o cohband.o \
		rtreref.o tranmt.o fop.o pwrvf.o readfvg.o wrt_z.o \
		reorder.o sdmwrite.o corchng.o rotsdm.o badToUse.o

INCLUDE=	iosize.inc

tranmtlr:		$(OBJS) $(INCLUDE)
		$(FC) $(FFLAGS) -o tranmtlr $(OBJS)

setup.o:	setup.f iosize.inc

mkfdir.o:	mkfdir.f iosize.inc

mkrec.o:	mkrec.f iosize.inc

readfvg.o:	readfvg.f iosize.inc

fop.o:		fop.f iosize.inc

badToUse.o:	badToUse.f

tranmt.o:	tranmt.f iosize.inc

clean:
	rm *.o

install:
		\mv tranmtlr $(BIN_DIR)/.
