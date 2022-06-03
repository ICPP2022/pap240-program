!****************************************************************************************
!> Module constants (list of all constants).
!****************************************************************************************
module mod_constants
    implicit none

    integer, parameter :: SPECIES = 3			            !<Number of species (must >=3), the first bit is (+)I/(-)V, the last bit is for im_SIA loop
    integer, parameter :: LEVELS = 8                        !<LEVELS=2^numSpecies
    !integer, parameter :: DIMENSIONS = 3                   !<The dimension (x/y/z)
    !integer, parameter :: SECTORS = 8                      !<Number of sectors
    double precision, parameter :: INITAU = 1d-2            !<Initial tau
    double precision, parameter :: NEVENT = 10d0            !<Number of events in a synchronous cycle

    double precision, parameter :: KB = 8.625d-5	        !<Boltzmann's constant (eV/K)
    double precision, parameter :: PI = 3.141592653589793d0	!<Pi
    double precision, parameter :: ZI = 1.2d0				!<Dislocation absorption efficiency for  interstitials
    double precision, parameter :: Zv = 1.0d0               !<Dislocation absorption efficiency for vacancies
    double precision, parameter :: ED = 40d0                !<Threshold energy of displacement (Fe:40eV, W:)

    !<input files
    integer, parameter :: PARAFILE = 10                     !<Used to read parameter.txtx file
    integer, parameter :: ATTRFILE = 11                     !<Used to read *_defects_*.txtx file
    !integer, parameter :: MESHFILE = 12                     !<Used to read Mesh_*.txt file
    integer, parameter :: CASFILE = 13                      !<Used to read cascades.txt File
    integer, parameter :: CPDFFILE = 14                     !<PKA spectrum
    integer, parameter :: IMPLFILE = 15                     !<Nonuniform DPA implantation profile
    !<output files
    integer, parameter :: TOTFILE = 81                      !<Used to write totdat.out file
    integer, parameter :: DEFFILE = 82                      !<Used to write defect.out file
    integer, parameter :: STAFILE = 83                      !<Used to write stadat.out file
    integer, parameter :: VTKFILE = 84                      !<Used to write defVTK.vtk file
contains
end module mod_constants

