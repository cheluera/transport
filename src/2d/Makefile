#                                                                           
#
include make.inc
#
#==========================================#
# List object files needed for compilation #
#==========================================#
OBJS0   = sn-2d.o    
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
	/usr/local/lf6481/bin/lf95  $(FFLAGS) -o driver  $(OBJS0) $(LIBRARY)  $(OTHERLIB)
#
f.o: f.f90
	/usr/local/lf6481/bin/lf95 -c $(FFLAGS) $*.f.f90

































