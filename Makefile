all:
	gfortran -O3 -c Utility.F90
	gfortran -O3 -c ExpIntegral.F90
	gfortran -O3 -c Bickley.F90
	gfortran -O3 -c MatrixInverse.F90
	gfortran -O3 -c Eigen.F90
	gfortran -O3 -c ColProbSlab.F90
	gfortran -O3  Main.F90 Utility.o MatrixInverse.o Eigen.o ColProbSlab.o ExpIntegral.o Bickley.o
