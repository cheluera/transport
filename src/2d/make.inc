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
#FORTRAN = lf95
FORTRAN = ifort
#
#
#================#
# Compiler flags #
#================#
##FFLAGS =  -KOMP -SSL2 --parallel
#FFLAGS = --o2 --parallel --ocl --threads 4 
#FFLAGS = --o2 --staticlink
#FFLAGS =  --chk  -g  --info
#FFLAGS = -g -Wall -fbounds-check -fcheck=all
#FFLAGS = -O3 -g -traceback -check all -fp-stack-check
FFLAGS = -O3
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
#OTHERLIB = -lnsl -lblas -llapack
OTHERLIB = 
#
##



