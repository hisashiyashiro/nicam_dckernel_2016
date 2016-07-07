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

PRJ = dckernel_vi_rhow_solver.exe

PRG = main

MODS = \
	dumpio.o \
	mod_precision.o \
	mod_misc.o      \
	mod_vi.o

makedir:
	mkdir -p $(BUILDDIR)
	mkdir -p $(NICAM_BINDIR)

makebin: makedir $(PRJ)



dumpio.o        : dumpio.c
mod_precision.o : mod_precision.f90
mod_misc.o      : mod_misc.f90 problem_size.inc mod_precision.o
mod_vi.o        : mod_vi.f90 mod_misc.o mod_precision.o

main.o : main.f90 mod_precision.o dumpio.o mod_misc.o mod_vi.o

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

.SUFFIXES: .o .f90 .c .mod
%.o : %.f90
	$(FC) $(FFLAGS) $(MODDIROPT) -o $(OBJDIR)/$@ -c $<

%.o : %.c
	$(CC) $(CFLAGS) -o $(OBJDIR)/$@ -c $<

%.o: %.mod

%.mod: %.f90
	$(make) $(patsubst %.f90,%.o,$<)

.PHONY : clean cleandir allclean
