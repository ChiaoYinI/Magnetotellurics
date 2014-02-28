OBJECT = cgeqrf.o ctrtrs.o cunmqr.o sgesv.o
FC = f77

libqr.a:	$(OBJECT)
	ar ruv libqr.a $(OBJECT)
	ranlib libqr.a
	rm *.o
	mv libqr.a ..
