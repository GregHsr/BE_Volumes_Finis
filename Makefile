FC = gfortran
# Avec options de debuggage activees
OPT = -g -O0 -fbounds-check 
# Sans option de debuggage
#OPT = -g

OBJ = m_type.o main.o subroutines.o VTSWriter.o
EXE = prog.exe

prog:	$(OBJ)
	$(FC) $(OPT) $(OBJ) -o $(EXE)

m_type.o:	m_type.f90
	$(FC) $(OPT) m_type.f90 -c	

main.o :	main.f90
	$(FC) $(OPT) main.f90 -c

subroutines.o :	subroutines.f90
	$(FC) $(OPT) subroutines.f90 -c

VTSWriter.o :	VTSWriter.f90
	$(FC) $(OPT) VTSWriter.f90 -c

clean :
	/bin/rm -f $(OBJ) $(EXE) $(MODS) *.mod


