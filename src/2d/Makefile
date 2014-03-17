#                                                                           
#
include make.inc
#
#==========================================#
# List object files needed for compilation #
#==========================================# 
OBJS0 = sn-2d.o sph-source.o spharm.o ext_source.o f_moments.o mfd_solution.o \
	manufactured.o error_2.o error_3.o output.o
#
#  
#==========================#
# Custum libraries to link #
#==========================#
LIBRARY  = 
#
#================================================
#
#
#=====================================#
# "DPdriver-a" is the executable name #
#=====================================#
a: $(OBJS0) $(LIBRARY)
	$(FORTRAN)  $(FFLAGS) -o driver  $(OBJS0) $(LIBRARY)  $(OTHERLIB)
#
%.o: %.f90
	$(FORTRAN) -c $(FFLAGS) $<

%.o: %.f
	$(FORTRAN) -c $(FFLAGS) $<

































