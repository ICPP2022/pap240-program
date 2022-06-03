!***************************************************************************************************
!>Module StructureType: contains all derived structures created for MISA-SCD.
!>Information read in from input file (
!										Formation energies: formation
!										Diffusivities: diffusion, diffusionFunction
!										Binding energies: binding, bindingFunction
!										Reactions: reactionParameters)
!										Cascade: cascade
!>Information about processors and meshes (processor, mesh, boundaryMesh)
!>Information about defects and reactions (defect, defectUpdate, reaction)
!>Information about cascade (cascadeDefect, cascadeEvent)
!***************************************************************************************************
module mod_structuretype
	implicit none
	!**********************************************************************
	!>Constant type: list of defect attributes read in from an input file.
	!**********************************************************************
	type :: formationSingle			!<List of formation energies for point defects read in from input file
		integer, allocatable :: defectType(:)			!<Type of point defect
		double precision :: Ef							!<Formation energy (eV)
	end type formationSingle

	type :: diffusionSingle			!<List of diffusion prefactors and migration energies for mobile defects read in from input file
		integer, allocatable :: defectType(:)			!<Type of defect that is allowed to diffuse
		double precision :: D0							!<Diffusion prefactor (nm^2/s)
		double precision :: Em							!<Migration energy (eV)
	end type diffusionSingle

	type :: diffusionFunction	!<List of functional forms for defects diffusivities read in from input file
		integer, allocatable :: defectType(:)			!<Type of defect that is allowed to diffuse using this functional form
		integer, allocatable :: min(:)					!<Minimum defect size allowed to diffuse using this functional form
		integer, allocatable :: max(:)					!<Maximum defect size allowed to diffuse using this functional form (-1 indicates infinity)
		integer :: fType								!<ID number of the function form to use to calculate diffusivity
		integer :: numParam								!<Number of parameters required for this functional form
		double precision, allocatable :: para(:)		!<Parameters required for this functional form
	end type diffusionFunction

	type :: bindingSingle				!<List of binding energies for small defects read in form input file
		integer, allocatable :: defectType(:)			!<Type of defect that can dissociate
		integer, allocatable :: product(:)				!<Type of defect dissociating away (point defect)
		double precision :: Eb							!<Binding energy (eV)
	end type bindingSingle

	type bindingFunction		!List of functional forms for defects binding energies read in from input file
		integer, allocatable :: defectType(:)			!<Type of cluster that is dissociating
		integer, allocatable :: product(:)				!<Type of defect that dissociates from cluster
		integer, allocatable :: min(:)					!<Minimum cluster size allowed to use this functional form
		integer, allocatable :: max(:)					!<Maximum cluster size allowed to use this functional form
		integer :: fType								!<ID number of the functional form to use to calculate binding energy
		integer :: numParam								!<Number of parameters required for this functional form
		double precision, allocatable :: para(:)		!<Parameters required for this functional form
	end type bindingFunction

	type :: reactionParameters	!<List of allowed reactions read in from input file
		integer :: numReactants							!<Number of reactants (-1:cascade  implantation
														!<						0:Frenkel Pair implantation
														!<						1: dissociation, diffusion, sinkRemoval, impurityTrapping
														!<						2: clustering)
		integer :: numProducts							!<Number of products
		integer, allocatable :: reactants(:,:)			!<Defect types for reactants in this reaction. Array size:(numSpecies, numReactants)
		integer, allocatable :: products(:,:)			!<Defect types for products in this reaction. Array size:(numSpecies, numProducts)
		integer, allocatable :: min(:)					!<Smallest defect size for reactants. Array size:(numProducts) for 1th reaction, (numProducts*2) for 2th reaction
		integer, allocatable :: max(:)					!<Largest defect size for reactants. Array size:(numProducts) for 1th reaction, (numProducts*2) for 2th reaction
		integer :: fType								!<ID number of function form used to calculate reaction rate,
	end type reactionParameters

	!**********************************************************************
	!>Structures related to processes and meshes
	!**********************************************************************
	type :: processor			!<Contains information about this processor
		integer :: taskid								!<ID number of this processor
		integer :: cartCoords(3)						!<Cartesian coordinates
		integer :: neighborProc(6)						!<Neighboring procesor (right, left, front, back, up, down).
		double precision :: localCoord(6)				!<Local boundaries of this processor (xmin, xmax, ymin, ymax, zmin, zmax)
	end type processor

	type :: mesh				!<Contains the mesh and connectivity of the meshes in the local processor as well as neighboring processors
		integer :: gid									!<Global ID of this mesh
		!integer :: coord(3)							!<Coordinates (integer) of center of this mesh
		!integer :: grainID								!<Grain ID number of this mesh (currently only set up for one grain)
		integer :: neighbor(6)							!<ID number of neighboring meshes, regardless of if they are in this processor or not. Six directions.
		integer :: neighborProc(6)						!<Processor ID numbers of neighboring mesh. Six directions.
	end type mesh

	type :: ghostMesh								    !<ghostMesh(6) -- 6 directions
		!integer :: proc									!<Processor ID of this mesh (typically different from the local processor)
        integer :: numCells						        !<Number of cells in this direction
		integer, allocatable :: cell(:)					!<Neighbor cells that are not in the processor -- (numCells)
		integer, allocatable :: local(:)				!<ID number of the local cell that borders the 'cell'
		!integer, allocatable :: grainID(:)				!<Grain ID of the 'cell' (currently only set up for one grain)
		type(defectHead),pointer :: myDRL_ghost(:,:)    !<DRL without reactions--(LEVELS, numCells)
	end type ghostMesh

	!**********************************************************************
	!>defect and reaction
	!**********************************************************************
	type :: defectHead									!<defectHead(LEVELS,numMeshes)
		type(defect), pointer :: defect_all				!<all defects, contain mobile and immobile defects
		type(defect), pointer :: defect_mobile			!<only mobile defects
	end type

	type :: defect				!<List of defects in a mesh--(LEVELS,numMeshes)???
		integer, allocatable :: defectType(:) 			!<Array containing the number of particles of each defect species in this defect type.
		integer :: num									!<Number of defects of this type inside this mesh
        double precision :: diff                      	!<Diffusivity
		type(reaction), pointer :: reactionList			!<first reaction related to a defect of type defectType
		type(defect), pointer :: next					!<Pointer to the next defect in the same mest
	end type defect

	type :: defectUpdate		!<Contains a list of defects that need to be updated
		integer, allocatable :: defectType(:) 			!<Type of defect that needs to be updated
		integer :: cell									!<ID number of the mesh that this defect is located in
		integer :: proc									!<Processor ID that this mesh is located inside
		integer :: dir									!<If defect diffused from a different mesh, indicates direction from which defect came
		integer :: neighbor								!<If defect diffused from a different mesh, indicates mesh number from which defect came
		integer :: cascadeID							!<If defect is inside a cascade, indicates cascade ID that defect is located inside
		type(defectUpdate), pointer :: next				!<Pointer to the next defect that needs to be updated
	end type defectUpdate

	type :: reaction			!<List of reactions in a mesh
		integer :: numReactants							!<Number of reactants in this reaction
		integer :: numProducts							!<Number of products in this reaction
		integer, allocatable :: reactants(:,:) 			!<Reactants in this reaction - array size:(numSpecies, numReactants)
		integer, allocatable :: products(:,:)			!<Products in this reaction - array size:(numSpecies, numProducts)
		integer, allocatable :: cell(:) 				!<Mesh numbers of defects involved in this reaction. Array size: (numReactants+numProducts). Important for diffusion reactions
		integer, allocatable :: taskid(:) 				!<Processor number of defects involved in this reaction. Array size: (numReactants+numProducts). May not all be the same processor number in the case of diffusion across processor boundaries.
		double precision :: rate
		type(reaction), pointer :: next					!<Pointer to the next reaction in the list in the same mesdh
	end type reaction

	!**********************************************************************
	!>defect and reaction
	!**********************************************************************
	type :: cascade				!<List of the cascades that are read in from a input file
		integer :: numDefects							!<Number of total defects in the cascade list
		integer :: numDisplacedAtoms					!<How much this cascade contributes to the DPA (how many lattice atoms are displaced)
		type(cascadeDefect), pointer :: listDefects		!<Pointer pointing to the list of defects in this cascade
		type(cascade), pointer:: next					!<Pointer pointing to the next cascade in the list
	end type cascade

	!<cascade files
	type :: cascadeFileList		!<List of cascades that are read in from a number of cascade files
		double precision :: temperature					!<temperature of the cascade collision
		double precision :: PKAenergy					!<PKA energy of the cascade collision
		double precision :: averDisAtoms				!<average number of atoms displaced per cascade in this cascade file
		integer :: numCascades							!<number of cascade in this cascade file
		type(cascade), pointer :: listCascades			!<list of cascades in this cascade file
	end type cascadeFileList

	type :: cascadeDefect			!<Used to read in defects from the cascade file. Contains a list of defects as well as their locations within a cascade
		integer, allocatable :: defectType(:) 			!<Type of defect (numSpecies)
		double precision :: coordinates(3) 				!<Location of defect relative to center of cascade (nm)
		type(cascadeDefect), pointer :: next			!<Pointer pointing to the next defect in the list
	end type cascadeDefect

    !Temporarily useless
	type :: cascadeEvent			!<Pointer list containing defects and reactions for cascades that are currently present in the system
		integer :: cell									!<Coarse mesh number that this cascade is located inside
		integer :: cascadeID							!<ID number of this cascade
		type(defect), pointer :: localDefects(:)		!<Array of defect lists (one for each element in the cascade)
		type(reaction), pointer :: reactionList(:) 		!<Array of reaction lists (one for each element in the cascade)
		double precision, allocatable :: totalRate(:)	!<Sum of the rates of all reactions in each cascade element
		type(cascadeEvent), pointer :: prev				!<Pointer pointing to the previous cascade currently present in the system
		type(cascadeEvent), pointer :: next				!<Pointer pointing to the next cascade currently present in the system
	end type cascadeEvent

	!**********************************************************************
	!>Used to output all defects in the system
	!**********************************************************************
	type :: postPtr
		type(postDefect), pointer :: head				!<Head of postdefectList
	end type

	type :: postDefect
		integer, allocatable :: defectType(:) 			!<Array containing the number of particles of each defect species in this defect type.
		integer :: num									!<Number of defects of this type inside this mesh
		type(postDefect), pointer :: next					!<Pointer to the next defect in the same mest
	end type postDefect

	!**********************************************************************
	!>Used for SLKMC
	!**********************************************************************
	type :: sectorType
		!integer :: numCells
		integer :: commDir(2,3)							!<Communication dir (send/recv, dim)
		integer :: numCell(3)							!<Number of meshes in x, y, z directions of the sector
		!double precision :: totalRate					!<Total rate of all reactions in the sector
	end type sectorType

	type :: defectCommHead
		type(defectComm), pointer :: commHead
	end type

	type :: defectComm		!<Contains a list of defects that need to be updated
		integer, allocatable :: defectType(:) 			!<Type of defect that needs to be communicated
		integer :: num									!<Number of defects of this type
		integer :: cell									!<Cell where the defect is located
		integer :: neighbor								!<If defect diffused to a different mesh, indicates the cell number that defect diffused to
		integer :: dir									!<If defect diffused to a different mesh, indicates direction of defect diffused
		type(defectComm), pointer :: next				!<Pointer to the next defect that needs to be communicated
	end type defectComm

    type :: buffer                                      !<recv(3), 3 dims
        integer :: num
        integer, allocatable :: datas(:,:)
    end type buffer

	!<PKA spectrum
	type :: cpdf_t
		integer :: size
		double precision, allocatable :: energy(:)
		double precision, allocatable :: cpdf(:)
	end type cpdf_t

contains
end module
