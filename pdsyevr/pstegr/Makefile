include ../../make.inc

dobj1 = dlarra.o dlarrd.o disnan.o \
        dlarrb.o dstegr2.o dlar1v.o \
        dlarrc.o dlarrf.o dlarrf2.o dlarrv.o dlarre2.o dlarrk.o \
	dlarrd2.o dlarre2a.o dstegr2a.o dstegr2b.o \
	dlarrb2.o dlarrv2.o dlar1va.o



all: $(dobj1)

.f.o:
	$(F77) -c $(FLAGS) $*.f

clean:
	rm -f *.o 
