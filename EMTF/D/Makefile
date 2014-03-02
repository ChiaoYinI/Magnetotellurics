#	Makefile for DNFF
#	edit the BIN_DIR directory to install dnff in your binary directory
#	use make all to compile, link, install dnff and remove all .o files
#	for debugging exchange the -O2 and -g options in FFLAGS
#
# Compiled with out of bounds warnings.

#FC = g95
FC=gfortran
FFLAGS= 
# -g
#BIN_DIR = ${MTbin}
BIN_DIR = ../examples/bin


TARGETS=	dnff clean

OBJS=		afcor.o autocor.o badrec.o bdrcsu.o ft_subs.o chdec.o \
		cldrft.o cmove.o dcfilt.o dcimte.o decset.o demean.o \
		dnff.o dtrnd.o fcorsu.o filtcor.o freqout.o frstdif.o \
		getsp.o ltslv.o mkset.o ptdist.o pterst.o sort.o \
		inpu_bin.o out_pack.o mk_offst.o phs_shft.o \
		inpu_asc.o

all:		$(TARGETS)

clean:		
		rm -f *.o *~ core

install:
		\mv dnff $(BIN_DIR)/.

dnff:		$(OBJS)
		$(FC) $(FFLAGS) -o $@ $(OBJS)



afcor.o: afcor.f

autocor.o: autocor.f

badrec.o: badrec.f

bdrcsu.o: bdrcsu.f iounits.inc

ft_subs.o: ft_subs.f

chdec.o: chdec.f

cldrft.o: cldrft.f

cmove.o: cmove.f

dcfilt.o: dcfilt.f

dcimte.o: dcimte.f decimate.inc params2.inc

decset.o: decset.f decimate.inc iounits.inc params2.inc

demean.o: demean.f

dnff.o: dnff.f decimate.inc input.inc iounits.inc params1.inc params2.inc

dtrnd.o: dtrnd.f

fcorsu.o: fcorsu.f decimate.inc params2.inc

filtcor.o: filtcor.f

freqout.o: freqout.f

frstdif.o: frstdif.f

getsp.o: getsp.f iounits.inc params2.inc

ltslv.o: ltslv.f

mkset.o: mkset.f decimate.inc params2.inc

ptdist.o: ptdist.f

pterst.o: pterst.f decimate.inc params2.inc

sort.o: sort.f

inpu_bin.o: inpu_bin.f ../include/datsz.inc ../include/nchmx.inc input.inc iounits.inc

out_pack.o: out_pack.f iounits.inc

mk_offst.o: mk_offst.f decimate.inc params2.inc

phs_shft.o: phs_shft.f

inpu_asc.o: inpu_asc.f iounits.inc input.inc ../include/datsz.inc ../include/nchmx.inc
