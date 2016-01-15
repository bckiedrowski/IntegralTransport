module Input
  implicit none

  integer, parameter :: ProbDim = 2
  integer, parameter :: nIter   = 1

  ! >>>>> common data
  real(8), parameter :: base_SigmaT   = 1.0d0
  real(8), parameter :: base_pScatter = 0.9d0
  real(8), parameter :: base_SigmaF   = base_SigmaT * ( 1.d0 - base_pScatter )
  real(8), parameter :: base_nubar    = 1.d0

  ! >>>>> 1d data
  real(8), parameter :: SlabWidth_1D = 20.d0
  integer, parameter :: NMesh_1D     = 1280 !1000
  integer, parameter :: Reflect_1D   = 1
  logical, parameter :: split_1D     = .false.
  logical, parameter :: coarsen_1D   = .false.

  real(8), parameter :: split_XLeft_1D  = 9.d0
  real(8), parameter :: split_XRight_1D = 11.1d0

  ! >>>>> 2d data
  real(8), parameter :: SlabWidth_X_2D = 5.d0
  real(8), parameter :: SlabWidth_Y_2D = 5.d0
  integer, parameter :: NMesh_X_2D     = 10
  integer, parameter :: NMesh_Y_2D     = 10
  integer, parameter :: NAngles_2D     = 128
  integer, parameter :: Reflect_X_2D   = 1
  integer, parameter :: Reflect_Y_2D   = 1
  real(8), parameter :: RaySpacing_2D  = 0.01d0
  logical, parameter :: nonuniform_2D  = .false.
  logical, parameter :: coarsen_2D     = .false.

end module Input
