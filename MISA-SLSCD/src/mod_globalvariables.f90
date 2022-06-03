!****************************************************************************************
!> Module globalVariables (list of globally shared variables and pointers). Includes:
! 1) MPI, processor and mesh information
! 2) Defect list, reaction list and cascade list
! 3) Defect attributes
! 4) Simulation parameters read in from configure.in
! 5) Other simulation parameters that are computed during simulation
!****************************************************************************************
module mod_globalvariables
    use mod_constants
    use mod_structuretype
    implicit none

    !*************************************************************
    !>MPI, processor and mesh information
    !************************************************************
    !>Used for MPI
    integer :: comm                                                 !<New communication domain: created by MPI_CART_CREAT
    integer :: dims(3)                                              !<Number of processors in x, y, z.
    logical, dimension(3) :: periods=(/.true., .true., .true./)     !<Boundary conditions in x, y, z.  Default: periods(3) = {.true., .true., .true.}
    integer :: ierr							                        !<Return error information
    double precision :: allreduceTime
    double precision :: casCommTime, numCasTime
    double precision :: commTime11, commTime12, commTime21, commTime22
    integer :: commCount_1, commCount_2
    integer :: countEvent

    !>Processor information
    integer :: numtasks                                 !<Number of processors in this simulation
    double precision :: systemVolume				    !<Volume of the whole system (nm^3)
    double precision :: systemCoord(6)                 !<Global boundaries of system (xmin, xmax, ymin, ymax, zmin, zmax)
    double precision :: localVolume			            !<Volume of single processor's domain (nm^3)
    type(processor) :: myProc					            !<Contains processor information

    !>Mesh information
    integer :: meshType(3)                              !<Boundary conditions in x, y and z directions (the value is 1(periodic) or 0(freeSurfaces))
    integer :: totalMeshes                              !<Total meshes in the sysytem
    integer :: numMeshes                                !<Number of local meshes
    double precision :: meshLength                      !<Length of a coarse mesh
    integer :: numxGlobal, numyGlobal, numzGlobal       !<Number of global meshes in x-direction, y-direction, z-direction
    integer :: numxLocal, numyLocal, numzLocal          !<Number of local meshes in x-direction, y-direction, z-direction
    integer :: maxElement                               !<Used to determine the size of myBoundary
    type(mesh), allocatable :: myMesh(:)				!<Local meshes. Array size: myMesh(numMeshes)

    type(ghostMesh) :: myGhost(6)                       !<Boundary cells. Array size -- myBound(six directions)
    type(sectorType) :: mySector(8)                     !<Sectors

    !>Defect list
    type(defectHead), allocatable :: myDRL(:,:)          !<Defect-reaction list--(LEVELS, numMeshes)

    !>Reaction list
    type(reaction), pointer :: implantation(:)
    double precision, allocatable :: totalRateVol(:)	!<Total reaction rate in each mesh
    !double precision :: totalRate						!<Total reaction rate in this processor
    !double precision :: maxRate						    !<Max reaction rate in all processors
    double precision :: sinks(2)                        !<S=Sd+Sg, sinks(1)--SIA, sinks(2)--Vac

    !>Cascade list
    type(cascade),pointer :: cascadeList			    !<List of cascades (read from file) that can be implanted
    double precision :: PKAtemperature                  !<Temperature of this cascade collision (K)
    double precision :: PKAenergy                       !<Energy of this PKA (eV)
    double precision :: numDisplacedAtoms		        !<Number of atoms displaced per cascade, read from cascade file
    integer :: numCascades							    !<number of cascades in the cascade file
    type(cascadeEvent) :: activeCascades                !<List of fine meshes that are active due to recent cascade implantation. Contains defect lists and reaction lists)

    !*************************************************************
    !>Defect attributes read in from input file
    !************************************************************
    integer :: numFormSingle                            !<Number of formation energies of point defects in input file
    integer :: numDiffSingle 	                        !<Number of diffusivities of small defects in input file
    integer :: numDiffFunc	                            !<Number of functional forms for diffusivities in input files
    integer :: numBindSingle 	                        !<Number of binding energies of small defects in input file
    integer :: numBindFunc		                        !<Number of functional forms for binding energies in input files
    type(formationSingle),allocatable :: formSingle(:)  !<Parameters for formation energies of point defects. Array size: (numForm)
    type(diffusionSingle),allocatable :: diffSingle(:)	!<Parameters for diffusivities of small defects. Array size: (numDiff)
    type(diffusionFunction), allocatable :: diffFunc(:) !<Parameters for functional forms of diffusivities for defects. Array size: (numDiffFunc)
    type(bindingSingle), allocatable :: bindSingle(:)	!<Parameters for binding energies of small defects. Array size: (numBind)
    type(bindingFunction), allocatable :: bindFunc(:)	!<Parameters for functional forms of binding energies defects. Array size: (numBindFunc)

    !integer :: numImplantReaction	                    !<Number of implantation reactions in input file
    !integer :: numDissocReaction	                    !<Number of dissociation reactions in input file
    !integer :: numSinkReaction		                    !<Number of sink reactions in input file
    !integer :: numImpurityReaction	                    !<Number of impurity reactions in input file
    !integer :: numDiffReaction		                    !<Number of diffusion reactions in input file
    integer :: numClusterReaction	                    !<Number of clustering reactions in input file
    !type(reactionParameters), allocatable :: implantReactions(:)	!<List of allowed implantation reactions. Array size: (numImplantReaction)
    !type(reactionParameters), allocatable :: dissocReactions(:)	    !<List of allowed dissociation reactions. Array size: (numDissocReaction)
    !type(reactionParameters), allocatable :: sinkReactions(:)		!<List of allowed sink reactions. Array size: (numSinkReaction)
    !type(reactionParameters), allocatable :: impurityReactions(:)	!<List of allowed impurity reactions. Array size: (numImpurityReaction)
    !type(reactionParameters), allocatable :: diffReactions(:)		!<List of allowed diffusion reactions. Array size: (numDiffReaction)
    type(reactionParameters), allocatable :: clusterReactions(:)	!<List of allowed clustering reactions. Array size: (numClusterReaction)

    !*************************************************************
    !>Simulation parameters read in from configure.in
    !************************************************************
    !>Toggles
    character(len=20) :: irradiationType                !<('Cascade' for neutron irradiation, 'FrenkelPair' for electron irradiation, 'None' for no irradiation)
    character(len=20) :: implantScheme			        !<('MonteCarlo' or 'explicit'), used to determine if cascades are implanted through Monte Carlo algorithm or explicitly
    character(len=20) :: implantType			        !<('uniform' or 'nonUniform'), used to determine if defects are implanted uniformly or if DPA rate / He implant rate are given for each mesh
    character(len=5) :: grainBoundToggle	            !<('yes' or 'no'), used to determine whether or not we are using grain boundaries to remove defects from simulation
    character(len=5) :: SIAPinToggle                    !<('yes' or 'no') Toggle whether or not allowing m_SIA loops to pin m_SIA loops
    character(len=5) :: HeSIAToggle                     !<('yes' or 'no') Toggle whether or not allow He-SIA clusters to form

    !>Simulation parameters
    double precision :: temperature			            !<Temperature read in (K) - used when temp. changes several times during a simulation
    double precision :: soluteConc                      !<The initial solute concentration (Cu or Cr)
    double precision :: MnConc                          !<The initial Mn concentration
    double precision :: NiConc                          !<The initial Ni concentration
    double precision :: SiConc                          !<The initial Si concentration
    double precision :: PConc                           !<The initial P concentration
    integer :: initialVacancies                         !<The number of vacancies put in. Used in V&V
    integer :: initialSIAs                              !<The number of vacancies put in. Used in V&V
    double precision :: dpaRate				            !<DPA rate in dpa/s
    double precision :: totalDPA				        !<Total DPA in simulation
    double precision :: HeDPAratio                      !<Helium to dpa ratio (atoms per atom)
    double precision :: agingTime                       !<For thermal aging simulation
    !*********************
    !<Annealing parameters
    double precision :: annealTime                      !<Annealing time
    double precision :: annealTemperature               !<Initial temperature of anneal stage (K)
    character(len=20) :: annealType                     !<('mult' or 'add' or 'constant') toggles additive or multiplicative or constant temperature in anneal steps
    double precision :: annealTempInc                   !<Temperature increment at each annealing step (additive or multipliciative)
    integer :: annealSteps                              !<Number of annealing steps
    integer :: annealIter                               !<Current anneal step
    !logical annealToggle                                !<(.TRUE. if in annealing phase, .FALSE. otherwise) used to determine how to reset reaction rates (should we include implantation or not)
    !*********************
    !<Materials parameters
    double precision :: lattice				            !<Lattice constant (nm)  atomVol = lattice^3/2
    double precision :: burgers				            !<Magnitude of burgers vector, equal to lattice constant
    double precision :: reactionRadius                  !<Material parameter used for reaction distances (impacts reaction rates) (nm)
    double precision :: grainSize			            !<Mean free path before a defect is absorbed by a grain boundary (equal to grain size)
    double precision :: dislocDensity		            !<Density of dislocations (sinks for point defects)
    double precision :: impurityDensity		            !<Denstiy of impurity atoms (traps for SIA loops)
    double precision :: cascadeVolume			        !<Volume of cascade (used for cascade mixing probability)
    integer :: max3D			                        !<Largest SIA size that can diffuse in 3D as spherical cluster
    integer :: numGrains			                    !<Number of grains inside polycrystal (default 1)
    integer :: numSims				                    !<Number of times to repeat simulation
    !*********************
    !Output  parameters
    character(len=20) :: outputToggle                   !<Output mark, used to mark the current moment output information is in irradiation or aging or annealing
    character(len=5) :: totdatToggle			        !<('yes' or 'no'), used to toggle whether we output the totdat.out data file
    character(len=5) :: defectToggle			        !<('yes' or 'no'), used to toggle whether we output the defect.out data file
    character(len=5) :: stadatToggle			        !<('yes' or 'no'), used to toggle whether we output the stadat.out data file
    character(len=5) :: vtkToggle			            !<('yes' or 'no'), used to toggle whether we output the .vtk file
    integer :: minVoid                                  !<Only n>minVoid nV_0_0 clusters are used for calculating the average cluster radius and number density
    integer :: minLoop                                  !<Only n>minLoop 0_0_nI clusters are used for calculating the average cluster radius and number density
    integer :: minS                                     !<Only n>minS (V)_nS_0 clusters are used for calculating the average cluster radius and number density
    integer :: minSV                                    !<Only (n+m)>minSV nV_mS_0 clusters are used for calculating the average cluster radius and number density
    integer :: minBubble

    !*************************************************************
    !<Simulation parameters, to be computed during simulation
    !*************************************************************
    !>Constants used for clustering rates
    double precision :: omega                           !<Geometric constant for 3D spherical clustering (see Dunn et al. JNM 2013)
    double precision :: omega2D				            !<Geometric constant for clustering with dislocation loops (see Dunn et al. JNM 2013)
    double precision :: omega1D                         !<Geometric constant for clustering with dislocation loops (see Dunn et al. JNM 2013)
    double precision :: omegacircle1D                   !<Geometric constant for clustering with dislocation loops (see Dunn et al. JNM 2013)
    double precision :: atomVol				            !<Atomic volume (nm^3)

    !>Simulation parameters, to be computed during simulation
    double precision :: time                            !<Elapsed time
    double precision :: stopTime                        !<End time of run
    integer :: step                                     !<Current steps
    integer :: annealStep                               !<Current anneal steps
    double precision :: DPA					            !<Current total DPA tracker (not a parameter)
    double precision :: DPA_low
    double precision :: rateTau(2)                      !<Used for collective communication. rateTau(1)=maxRate, rateTau(2)=time increment
    integer :: numImpEvents                             !<Number of Frenkel pairs or cascades implantation events (local)
    integer :: totImpEvents                             !<Total number of Frenkel pairs or cascades implantation events (global)
    integer :: numImpHe                                 !<Number of He implantation events (local)
    integer :: totImpHe                                 !<Total number of He implantation events (global)
    double precision :: numDamages                     !<Number of displaced atoms (local)
    double precision :: totDamages                     !<Total number of displaced atoms (global)
    !<Record the runing time of the program
    double precision :: runTime
    !>Cu solubility CeqCu(T) = exp(DelatS/kB)*exp(-Omega/(kB*T))  Reference: (F. Christien and A. Barbu, 2004)
    double precision :: ceqSIA                          !<Thermal equilibrium concentration of SIA
    double precision :: ceqVac                          !<Thermal equilibrium concentration of vacancy
    double precision :: concSIA                         !<SIA concentration
    double precision :: concVac                         !<Vacancy concentration
    integer :: numSevermesh                             !<Initial number of Solute atoms in one mesh (Cu or Cr)
    integer :: initNumSIA                               !<Initial number of self-interstitial atoms in the whole system
    integer :: initNumVac                               !<Initial number of vacancies in the whole system
    integer, allocatable :: initPointDefects(:,:)       !<List the globalID of the mesh where initial SIA and vacancies are located--(max(initNumSIA, initNumVac),2): 1-SIA, 2-Vac

    !<PKA spectrum
    character(len=20) :: PKAspectrum                    !<('yes' or 'no') Whether to use the PKA spectrum
    integer :: numCascadeFiles                          !<cascadeFile: materials_temperature_PKAenergy_*.txt
    type(cascadeFileList), allocatable :: cascadeLists(:)   !<List of cascades (read from file) that can be implanted
    type(cpdf_t) :: EPKAlist                            !<list of PKA spectrum (cpdf and PKA energies)
    integer :: numImpDatas                              !<Number of values of DPA rate
    double precision, allocatable :: impDPARates(:,:)   !List of dpa rates

contains
end module mod_globalvariables

