include ../../../../../GNUmakefile.chssi

#Compilers
FC=ftn
CC=cc
CPP=CC
SERIAL_CPP=icpc

#Compiler flags
FFLAGS=-D__fortran -D__fortran77 -DMPIF2C -Zp8 -zero -traceback -DMPI2 -shared-intel -O3 -fp-model strict
CFLAGS=-DMPIF2C -DC_SUFFIX  -DMPI2 -shared-intel -O3 -fp-model strict -DCB_SUFFIX
CPPFLAGS=-DMPIF2C -DC_SUFFIX -DCB_SUFFIX  -DMPI2 -shared-intel -O3 -fp-model strict

#ACESLIBS variable will be replaced by the autoconf script with appropriate libraries
LIBS = -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core  -Wl,--end-group
LIBS:= -lstdc++ -lsip1 -lsip2 -lsip_shared -lframelib -laces2 -lgeopt -lsymcor -laces2 -lerd -loed -ldup -lsip1 -lsip2 $(LIBS)

#ACESFLAGS variable will be replaced by the autoconf script with location of libraries
LIB_DIRS= -L.

SIAL_COMPILER_LIBS= -lsial -lsip_shared -laces2  -lifcore

INCLUDE_DIRS := -I../../../../sia/include -I../aces2/include $(INCLUDE_DIRS)
LIB_DIRS:=-L../../../lib $(LIB_DIRS)

SYMCOR_LIB=../../../../../lib/libsymcor.a

DEFINES= 

FSOURCE=$(wildcard *.F)
F_OBJS=$(FSOURCE:.F=.o)
fSOURCE=$(wildcard *.f)
f_OBJS:=$(fSOURCE:.f=.o)
cSOURCE=$(wildcard *.c)
c_OBJS:=$(cSOURCE:.c=.o)
CPP_SOURCE=$(wildcard *.cpp)
CPP_OBJS:=$(CPP_SOURCE:.cpp=.o)

all: tabula_rasa $(SYMCOR_LIB)

tabula_rasa: ;
	rm -f make.out make.err

$(SYMCOR_LIB): libsymcor.a
	cp libsymcor.a $(SYMCOR_LIB)
	rm -f make.out make.err

libsymcor.a: $(F_OBJS) $(f_OBJS) $(c_OBJS) $(CPP_OBJS)
	ar $(ARFLAGS) libsymcor.a $(F_OBJS) $(f_OBJS) $(c_OBJS) $(CPP_OBJS)

libclean: ;
	rm -f libsymcor.a

ppclean clean distclean: % : tabula_rasa libclean
	rm -f *.o

rebuild: libclean $(SYMCOR_LIB)

binclean relink archive: % : ;

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.F.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.cpp.o: 
	$(CPP) -c $(CPPFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err


