OBJECT1 = caxpy.o ccopy.o cdotc.o cdotu.o crotg.o cscal.o csrot.o csscal.o cswap.o dasum.o daxpy.o dcabs1.o dcopy.o ddot.o dnrm2.o drot.o drotg.o drotm.o drotmg.o dscal.o dswap.o dzasum.o dznrm2.o icamax.o idamax.o isamax.o izamax.o sasum.o saxpy.o scasum.o scnrm2.o scopy.o sdot.o snrm2.o srot.o srotg.o srotm.o srotmg.o sscal.o sswap.o zaxpy.o zcopy.o zdotc.o zdotu.o zdrot.o zdscal.o zrotg.o zscal.o zswap.o

OBJECT2 = cgbmv.o cgemv.o cgerc.o cgeru.o chbmv.o chemv.o cher.o cher2.o chpmv.o chpr.o chpr2.o ctbmv.o ctbsv.o ctpmv.o ctpsv.o ctrmv.o ctrsv.o dgbmv.o dgemv.o dger.o dsbmv.o dspmv.o dspr.o dspr2.o dsymv.o dsyr.o dsyr2.o dtbmv.o dtbsv.o dtpmv.o dtpsv.o dtrmv.o dtrsv.o lsame.o sgbmv.o sgemv.o sger.o ssbmv.o sspmv.o sspr.o sspr2.o ssymv.o ssyr.o ssyr2.o stbmv.o stbsv.o stpmv.o stpsv.o strmv.o strsv.o xerbla.o zgbmv.o zgemv.o zgerc.o zgeru.o zhbmv.o zhemv.o zher.o zher2.o zhpmv.o zhpr.o zhpr2.o ztbmv.o ztbsv.o ztpmv.o ztpsv.o ztrmv.o ztrsv.o

OBJECT3 = cgemm.o chemm.o cher2k.o cherk.o csymm.o csyr2k.o csyrk.o ctrmm.o ctrsm.o dgemm.o dsymm.o dsyr2k.o dsyrk.o dtrmm.o dtrsm.o lsame.o sgemm.o ssymm.o ssyr2k.o ssyrk.o strmm.o strsm.o xerbla.o zgemm.o zhemm.o zher2k.o zherk.o zsymm.o zsyr2k.o zsyrk.o ztrmm.o ztrsm.o

FC = f77

libblas.a:	$(OBJECT1) $(OBJECT2) $(OBJECT3)
	ar ruv libblas.a $(OBJECT1) $(OBJECT2) $(OBJECT3)
	ranlib libblas.a
	rm *.o
	mv libblas.a ..
