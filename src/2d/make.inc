##############################
#  make include file         #
#  to change parameters      #
#  for different platforms   #
##############################
#
#
#=================================#
# Name of Fortran and C compilers #
#=================================#
FORTRAN = lf95
#
#
#================#
# Compiler flags #
#================#
##FFLAGS =  -KOMP -SSL2 --parallel
FFLAGS = --o2 --parallel --ocl --threads 4 
#FFLAGS = --o2 --staticlink
#FFLAGS =  --chk  -g  --info
####
####
CFLAGS     = $(FFLAGS) 
#
#====================#
# Personal libraries #
#====================#
DIRLIB  = $(HOME)/usr/lib/lf
#
#=====================#
# 3rd party libraries #
#=====================#
OTHERLIB = -lnsl -lblas -llapack
#
#



