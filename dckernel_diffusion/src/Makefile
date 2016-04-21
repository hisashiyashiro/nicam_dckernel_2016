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

PRJ = dckernel_diffusion.exe

PRG = main

MODS = \
	mod_precision.o \
	mod_misc.o      \
	mod_oprt.o

makedir:
	mkdir -p $(BUILDDIR)
	mkdir -p $(NICAM_BINDIR)

makebin: makedir $(PRJ)



mod_precision.o : mod_precision.f90
mod_misc.o      : mod_misc.f90 problem_size.inc mod_precision.o
mod_oprt.o      : mod_oprt.f90 mod_misc.o mod_precision.o

main.o : main.f90 mod_precision.o mod_misc.o mod_oprt.o

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
