all:
	gfortran -O3 -c Utility.F90
	gfortran -O3 -c ExpIntegral.F90
	gfortran -O3 -c Bickley.F90
	gfortran -O3 -c MatrixInverse.F90
	gfortran -O3 -c Eigen.F90
	gfortran -O3 -c Input.F90
	gfortran -O3 -c ColProbGeom.F90
	gfortran -O3  Main.F90 Utility.o MatrixInverse.o Eigen.o ColProbGeom.o ExpIntegral.o Bickley.o Input.o
