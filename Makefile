all:
	ifort -g -check all -c Utility.F90
	ifort -g -check all -c ExpIntegral.F90
	ifort -g -check all -c Bickley.F90
	ifort -g -check all -c MatrixInverse.F90
	ifort -g -check all -c Eigen.F90
	ifort -g -check all -c Input.F90
	ifort -g -check all -c ColProbGeom.F90
	ifort -g -check all  Main.F90 Utility.o MatrixInverse.o Eigen.o ColProbGeom.o ExpIntegral.o Bickley.o Input.o
