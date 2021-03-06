include ../../GNUmakefile.chssi

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
LIBS:= -lstdc++ -lsip1 -lsip2 -lsialinst -lsip_shared -lframelib -laces2 -lgeopt -lsymcor -laces2 -lprop_ints -lerd -loed -lecp -ldup -lsip1 -lsip2 -lsialinst -lsip1 -lsip2 $(LIBS)

#F12 
#LIBS:= -lstdc++ -lsip1 -lsip2 -lsialinst -lsip_shared -lframelib -laces2 -lgeopt -lsymcor -laces2 -lprop_ints -lerd -loed -lecp -ldup -lsip1 -lsip2 -lsialinst -lerd_F12 -lF12grid -loed_F12 -lsip1 -lsip2 $(LIBS)

#ACESFLAGS variable will be replaced by the autoconf script with location of libraries
LIB_DIRS= -L.

SIAL_COMPILER_LIBS= -lsial -lsip_shared -laces2  -lifcore




INCLUDE_DIRS := -I../sia/include -I../aces/aces_library/a2driver/aces2/include -I../sia/sip/include $(INCLUDE_DIRS)
LIB_DIRS:=-L../../lib $(LIB_DIRS)

FRAMELIB=../../lib/libframelib.a
ERD_LIB=../../lib/liberd.a
SIP_LIB1=../../lib/libsip1.a
SIP_LIB2=../../lib/libsip2.a
LIB_SIALINST=../../lib/libsialinst.a
SIAL_LIB=../../lib/libsial.a
SIP_SHARED_LIB=../../lib/libsip_shared.a
ACES2_LIB=../../lib/libaces2.a
GEOPT_LIB=../../lib/libgeopt.a
SYMCOR_LIB=../../lib/libsymcor.a
OED_LIB=../../lib/liboed.a
ECP_LIB=../../lib/libecp.a
#OEDF12_LIB=../../lib/liboed_F12.a
#ERDF12_LIB=../../lib/liberd_F12.a
#F12GRID_LIB=../../lib/libF12grid.a
PROP_INTS_LIB=../../lib/libprop_ints.a
DUP_LIB=../../lib/libdup.a
INTEGRAL_OBJS=$(shell ls ../aces/aces_library/special_directory/direct_integral/*.o 2>/dev/null)
WORKER_OBJS=$(shell ls ../sia/worker/*.o 2>/dev/null)
MANAGER_OBJS=$(shell ls ../sia/manager/*.o 2>/dev/null)

ACES3_EXE=../../bin/xaces3

FSOURCE=$(wildcard *.F)
F_OBJS=$(FSOURCE:.F=.o)
fSOURCE=$(wildcard *.f)
f_OBJS:=$(fSOURCE:.f=.o)
cSOURCE=$(wildcard *.c)
c_OBJS:=$(cSOURCE:.c=.o)
CPP_SOURCE=$(wildcard *.cpp)
CPP_OBJS:=$(CPP_SOURCE:.cpp=.o)

MAIN_OBJS=$(F_OBJS) $(f_OBJS) $(c_OBJS) $(CPP_OBJS)

all: tabula_rasa $(ACES3_EXE) 

tabula_rasa: ;
	rm -f make.out make.err

$(ACES3_EXE): xaces3
	cp xaces3 $(ACES3_EXE)
	rm -f make.out make.err

xaces3: $(FRAMELIB) $(ERD_LIB) $(SIP_LIB1) $(SIP_LIB2) $(LIB_SIALINST) $(SIAL_LIB) $(SIP_SHARED_LIB) \
        $(ACES2_LIB) $(GEOPT_LIB) $(SYMCOR_LIB) $(OED_LIB) $(ECP_LIB) $(DUP_LIB) $(INTEGRAL_OBJS) \
        $(OEDF12_LIB) $(ERDF12_LIB) $(F12GRID_LIB) $(MAIN_OBJS) $(WORKER_OBJS) $(MANAGER_OBJS)
	$(FC) $(FFLAGS) $(LFLAGS) -o xaces3 \
	       $(MAIN_OBJS) $(MANAGER_OBJS) $(WORKER_OBJS) $(INTEGRAL_OBJS) \
	      $(LIB_DIRS) $(LIBS) && test -x xaces3  || rm -f xaces3

binclean: ;
	rm -f xaces3

ppclean clean distclean: % : tabula_rasa binclean
	rm -f *.o

rebuild relink: % : binclean $(ACES3_EXE)

libclean archive: % : ;

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.F.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.cpp.o: 
	$(CPP) -c $(CPPFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

