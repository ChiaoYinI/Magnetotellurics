##  Multiple station processing program

#  change BIN_DIR to desired bin directory 
BIN_DIR= ${MTbin}
#FC = g95 
FC=f77

#  options for debug compilation
FCFLAGS = -g 

#  compiler options for Sun Solaris
 FFLAGS = -O -L/usr/local/lib

#  following are options for compiling on IBM RS 6000/SP2 cluster
#FFLAGS = -O -qrealsize=4 -qfixed -qmaxmem=-1 \
#         -L/usr/local/lib -bmaxdata:512000000 \
#         -bmaxstack:512000000 -qSPILL=3270
# end of IBM options

#  Libraries to link to (assuming liblapack.a and liblas.a are in
#				/usr/local/lib  )
LIBS = -llapack -lblas

### ALTERNATIVE :: local, very partial lapack library
###    that you can make with source code include with EMTF/MMT
###    (only QR, general linear solver, BLAS)
#FFLAGS = -O -Llapack
#LIBS = -lqr -lblas

TARGETS=	multmtrn install clean

OBJS=		anfld.o mk_fd_rec.o mmt_mtrx.o mmt_subs.o multmtrn.o \
		prteig.o rbstk2.o readfvg.o rsp.o setup.o timerr.o \
		ln_rbst.o sindex.o extra.o n_rbst.o grad.o rbstreg.o \
		var_adj.o rr_rbst.o pc_out.o wrtx.o \
		filtpc.o canonic.o cn_wt.o sep_s_n.o wrt_z.o \
		refproj.o mk_list.o corchng.o minpwr.o badToUse.o

all:		$(TARGETS)

clean:
		rm -f *.o *~ core


install:
		\mv multmtrn $(BIN_DIR)/.

multmtrn:	$(OBJS) 
		$(FC) $(FFLAGS) -o $@ $(OBJS) $(LIBS)

canonic.o:	canonic.f nstamx.inc
	$(FC) $(FCFLAGS) -c canonic.f     

cn_wt.o:	cn_wt.f nstamx.inc
	$(FC) $(FCFLAGS) -c cn_wt.f      

emap.o:	emap.f nstamx.inc
	$(FC) $(FCFLAGS) -c emap.f      

extra.o:	extra.f nstamx.inc
	$(FC) $(FCFLAGS) -c extra.f      

filtpc.o:	filtpc.f nstamx.inc
	$(FC) $(FCFLAGS) -c filtpc.f      

ln_rbst.o:	ln_rbst.f nstamx.inc iosize.inc
	$(FC) $(FCFLAGS) -c ln_rbst.f      

mk_list.o:	mk_list.f iosize.inc
	$(FC) $(FCFLAGS) -c mk_list.f      

lnoise.o:	lnoise.f nstamx.inc
	$(FC) $(FCFLAGS) -c lnoise.f      

mmt_subs.o:	mmt_subs.f nstamx.inc
	$(FC) $(FCFLAGS) -c mmt_subs.f      

pc_coeff.o:	pc_coeff.f nstamx.inc
	$(FC) $(FCFLAGS) -c pc_coeff.f      

pc_out.o:	pc_out.f nstamx.inc
	$(FC) $(FCFLAGS) -c pc_out.f      

pcregtf.o:	pcregtf.f nstamx.inc
	$(FC) $(FCFLAGS) -c pcregtf.f      

rbstk2.o:	rbstk2.f nstamx.inc
	$(FC) $(FCFLAGS) -c rbstk2.f      

rr_rbst.o:	rr_rbst.f nstamx.inc
	$(FC) $(FCFLAGS) -c rr_rbst.f      

wrt_z.o:		wrt_z.f nstamx.inc
	$(FC) $(FCFLAGS) -c wrt_z.f      

mk_fd_rec.o:	mk_fd_rec.f iosize.inc nstamx.inc
	$(FC) $(FCFLAGS) -c mk_fd_rec.f      

multmtrn.o:	multmtrn.f iosize.inc nstamx.inc
	$(FC) $(FCFLAGS) -c multmtrn.f      

n_rbst.o:	n_rbst.f iosize.inc nstamx.inc
	$(FC) $(FCFLAGS) -c n_rbst.f      

readfvg.o:	readfvg.f iosize.inc nstamx.inc
	$(FC) $(FCFLAGS) -c readfvg.f      

setup.o:		setup.f iosize.inc nstamx.inc
	$(FC) $(FCFLAGS) -c setup.f      

timerr.o:	timerr.f iosize.inc nstamx.inc
	$(FC) $(FCFLAGS) -c timerr.f      

refproj.o:	refproj.f
	$(FC) $(FCFLAGS) -c refproj.f

corchng.o:	corchng.f
	$(FC) $(FCFLAGS) -c corchng.f

minpwr.o:	minpwr.f
	$(FC) $(FCFLAGS) -c minpwr.f

badToUse.o:	badToUse.f
	$(FC) $(FCFLAGS) -c badToUse.f
