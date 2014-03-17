# Set compiler commands
#------------------------------------
FC90 = ifort #Change if necessary
MP90 = mpif90 
FLAGS = -O3 -g -traceback -check all -fp-stack-check
#FLAGS = -03 
#------------------------------------
main_files  = func_ang.f90 func_mfd.f90 func_out.f90 sn-2d-parallel.f90

all: sn-2d_exe

sn-2d_exe: $(main_files)
	$(MP90) $(FLAGS) $(main_files) -o $@

clean: 
	rm sn-2d_exe *~ *.mod *.txt *.e* *.o*

#------------------------------------
# modules
#------------------------------------
func_ang.mod: func_ang.f90
	$(FC90) $(FLAGS) -c func_ang.f90

func_mfd.mod: func_mfd.f90
	$(FC90) $(FLAGS) -c func_mfd.f90

func_out.mod: func_out.f90
	$(FC90) $(FLAGS) -c func_out.f90
