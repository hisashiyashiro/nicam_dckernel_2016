################################################################################
#
# Makefile grand central
#
################################################################################

TOPDIR   = ../..
BUILDDIR = .

VPATH = $(BUILDDIR)

include $(TOPDIR)/sysdep/Makedef.$(NICAM_SYS)
include $(TOPDIR)/Mkinclude

all:
	$(MAKE) makedir
	$(MAKE) makebin

allclean: clean cleandir

PRJ = physicskernel_microphysics.exe

PRG = main

MODS = \
	dumpio.o               \
	mod_precision.o        \
	mod_debug.o            \
	mod_thrmdyn.o          \
	mod_satadjust.o        \
	mod_vadv1d.o           \
	mod_precip_transport.o \
	mod_mp_nsw6.o          \
	mod_mp_driver.o

makedir:
	mkdir -p $(BUILDDIR)
	mkdir -p $(NICAM_BINDIR)

makebin: makedir $(PRJ)



dumpio.o               : dumpio.c
mod_precision.o        : mod_precision.f90
mod_debug.o            : mod_debug.f90 problem_size.inc mod_precision.o

mod_thrmdyn.o          : mod_thrmdyn.f90          mod_precision.o mod_debug.o
mod_satadjust.o        : mod_satadjust.f90        mod_precision.o mod_debug.o mod_thrmdyn.o
mod_vadv1d.o           : mod_vadv1d.f90           mod_precision.o mod_debug.o
mod_precip_transport.o : mod_precip_transport.f90 mod_precision.o mod_debug.o mod_thrmdyn.o mod_vadv1d.o

mod_mp_driver.o        : mod_mp_driver.f90        mod_precision.o mod_debug.o mod_satadjust.o mod_mp_nsw6.o
mod_mp_nsw6.o          : mod_mp_nsw6.f90          mod_precision.o mod_debug.o mod_satadjust.o mod_thrmdyn.o mod_precip_transport.o

main.o : main.f90 mod_precision.o dumpio.o mod_debug.o mod_mp_driver.o

$(BUILDDIR)/$(PRJ) : $(PRG).o $(BUILDDIR)/$(MODS)
	$(LD) $(LFLAGS) -o $@ $^
	$(INSTALL) $@ $(NICAM_BINDIR)
	@echo; echo "##### PROGRAM $(notdir $@) is installed Successfully #####"; echo



cleandir:
	-rm $(NICAM_BINDIR)/*

clean:
	rm -f $(PRJ)
	rm -f *~ *.mod *.o *.lst *.L
	cd $(BUILDDIR); rm -f $(PRJ); rm -f *~ *.mod *.o *.lst *.L

.SUFFIXES: .o .f90
.f90.o:
	$(FC) $(FFLAGS) $(MODDIROPT) -o $(OBJDIR)/$@ -c $<

%.o: %.mod

$(BUILDDIR)/%.mod: %.f90
	$(make) $(patsubst %.f90,$(BUILDDIR)/%.o,$<)

.PHONY : clean cleandir allclean
