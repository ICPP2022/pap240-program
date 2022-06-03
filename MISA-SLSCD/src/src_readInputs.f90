! Created by ${USER_NAME} on 2019/11/6.
!***************************************************************************************************
!> Subroutine ReadInputs():  - reads in simulation parameters from configure.in
!This subroutine reads in all simulation parameters located in configure.in as well file names
!for all other input files (defect attributes, mesh, cascades, implantation, etc).
!***************************************************************************************************
subroutine ReadInputs()
    use mod_constants
    use mod_globalvariables
    implicit none

    character(len=20) :: char
    logical alive1, alive2, flag, flag1
    integer :: i, numArgs
    character(len=100) :: arg                            !<command argument
    character(len=100) :: defectFilename                 !<Filename of defect attributes file
    character(len=100) :: meshFilename                   !<Filename of mesh file
    character(len=100) :: cascadeFilename                !<Filename of cascade file
    character(len=100) :: pkaFilename                    !<Filename of PKA spectrum file
    character(len=100) :: implantFilename                !<Filename of implant file

    flag= .false.
    !<read in filename of defectFile
    numArgs=iargc()
    if(numArgs > 0) then
        call getarg(numArgs, arg)
    else
        write(*,*) 'ERROR! No input file is provided before running. To run, you should type: command inputFileName'
    end if
    inquire(file=arg, exist=alive1)
    if(.not. alive1) then
        write(*,*) arg, 'does not exist'
    else
        open(unit=PARAFILE, file=arg, status='old',action='read')
    end if

    !*******************************************************
    !<read in toogles
    !*******************************************************
    !<set default valuse for toogles
    implantScheme = 'MonteCarlo'
    implantType = 'uniform'
    grainBoundToggle = 'no'
    SIAPinToggle = 'no'
    HeSIAToggle = 'no'
    meshType = 1                !<meshType(3), periodic
    PKAspectrum = 'no'
    numCascadeFiles = 1

    !<read in filename of defect attributes file
    do while(flag .eqv. .false.)
        read(PARAFILE,*) char
        if(char=='defectFile') then
            read(PARAFILE,*) defectFilename
            flag=.true.
        end if
    end do
    flag= .false.

    !<read in filename of
!    do  while (flag .eqv. .false.)
!        read(PARAFILE,*) char
!        if(char=='meshFile') then
!            read(PARAFILE,*) meshFilename
!            flag=.true.
!        end if
!    end do
!    flag=.false.

    !read in irradiation type
    do while(flag .eqv. .false.)
        read(PARAFILE,*) char
        if(char=='irradiationType') then
            read(PARAFILE,*) irradiationType
            flag=.true.
        end if
    end do
    flag=.false.

    !read in cascade implantation scheme (Monte Carlo or explicit)
    do while(flag .eqv. .false.)
        read(PARAFILE,*) char
        if(char=='implantScheme') then
            read(PARAFILE,*) implantScheme
            flag=.true.
        end if
    end do
    flag=.false.

    !Whether to use the PKA spectrum
    do while (flag .eqv. .false.)
        read(PARAFILE,*) char
        if(char=='PKAspectrum') then
            read(PARAFILE,*) PKAspectrum
            flag=.true.
        end if
    end do
    flag=.false.

    !read in filename of PKA spectrum file
    do while (flag .eqv. .false.)
        read(PARAFILE,*) char
        if(char=='pkaFile') then
            read(PARAFILE,*) pkaFilename
            flag=.true.
        end if
    end do
    flag=.false.

    !read in the number of cascade files
    do while (flag .eqv. .false.)
        read(PARAFILE,*) char
        if(char=='numCascadeFiles') then
            read(PARAFILE,*) numCascadeFiles
            flag=.true.
        end if
    end do
    flag=.false.

    !read in filename of cascade file
    do while (flag .eqv. .false.)
        read(PARAFILE,*) char
        if(char=='cascadeFile') then
            read(PARAFILE,*) cascadeFilename
            flag=.true.
        end if
    end do
    flag=.false.

    !read in toggle for nonUnifrm distribution
    do while(flag .eqv. .false.)
        read(PARAFILE,*) char
        if(char=='implantType') then
            read(PARAFILE,*) implantType
            flag=.TRUE.
        end if
    end do
    flag=.false.

    !read in filename of implant file
    do while(flag .eqv. .false.)
        read(PARAFILE,*) char
        if(char=='implantFile') then
            read(PARAFILE,*) implantFilename
            flag=.true.
        end if
    end do
    flag=.false.

    !read in grain boundary toggle
    do while(flag .eqv. .false.)
        read(PARAFILE,*) char
        if(char=='grainBoundaries') then
            read(PARAFILE,*) grainBoundToggle
            flag=.true.
        end if
    end do
    flag=.false.

    !read in SIAPin toggle
    do while(flag .eqv. .false.)
        read(PARAFILE,*) char
        if(char=='SIAPinToggle') then
            read(PARAFILE,*) SIAPinToggle
            flag=.true.
        end if
    end do
    flag=.false.

    !read in HeSIA toggle
    do while(flag .eqv. .false.)
        read(PARAFILE,*) char
        if(char=='HeSIAToggle') then
            read(PARAFILE,*) HeSIAToggle
            flag=.true.
        end if
    end do
    flag=.false.

    !*******************************************************
    !<read in simulartion parameters
    !*******************************************************
    !<set default valuse for toogles
    temperature = 273d0
    soluteConc = 0d0
    initialVacancies = 0
    initialSIAs = 0
    dpaRate = 1d-4
    totalDPA = 0d0
    HeDPAratio = 0d0
    agingTime = 0d0
    lattice = 0.2867d0
    burgers = 0.287d0
    reactionRadius = 0.65d0
    grainSize = 3.0d4
    dislocDensity = 0d0
    impurityDensity = 0d0
    cascadeVolume = 1000d0
    max3D = 1
    numGrains = 1
    numSims = 1

    !<read in variables
    do while(flag .eqv. .false.)
        read(PARAFILE,*) char
        if(char=='SimulationStart') then
            flag=.true.
        end if
    end do
    flag=.false.

    do while(flag .eqv. .false.)
        flag1=.false.
        do while(flag1 .eqv. .false.)
            read(PARAFILE,*) char
            if(char=='SimulationEnd') then
                flag1=.true.
                flag=.true.
            else if(char=='temperature') then
                flag1=.true.
                read(PARAFILE,*) temperature
            else if(char=='soluteConc') then
                flag1=.true.
                read(PARAFILE,*)  soluteConc
            else if(char=='MnContent') then
                flag1=.true.
                read(PARAFILE,*)  MnConc
            else if(char=='NiContent') then
                flag1=.true.
                read(PARAFILE,*)  NiConc
            else if(char=='SiContent') then
                flag1=.true.
                read(PARAFILE,*)  SiConc
            else if(char=='PContent') then
                flag1=.true.
                read(PARAFILE,*)  PConc
            else if(char=='initialVacancies') then
                flag1=.true.
                read(PARAFILE,*) initialVacancies
            !else if(char=='initialSIAs') then
            !    flag1=.true.
            !    read(PARAFILE,*) initialSIAs
            else if(char=='dpaRate') then
                flag1=.true.
                read(PARAFILE,*) dpaRate
            else if(char=='totalDPA') then
                flag1=.true.
                read(PARAFILE,*) totalDPA
            else if(char=='HeDPAratio') then
                flag1=.true.
                read(PARAFILE,*) HeDPAratio
            else if(char=='agingTime') then
                flag1=.true.
                read(PARAFILE,*) agingTime
            else if(char=='lattice') then
                flag1=.true.
                read(PARAFILE,*) lattice
            else if(char=='burgers') then
                flag1=.true.
                read(PARAFILE,*) burgers
            else if(char=='reactionRadius') then
                flag1=.true.
                read(PARAFILE,*) reactionRadius
            else if(char=='grainSize') then
                flag1=.true.
                read(PARAFILE,*) grainSize
            else if(char=='dislocDensity') then
                flag1=.true.
                read(PARAFILE,*) dislocDensity
            else if(char=='impurityDensity') then
                flag1=.true.
                read(PARAFILE,*) impurityDensity
            else if(char=='cascadeVolume') then
                flag1=.true.
                read(PARAFILE,*) cascadeVolume
            else if(char=='max3D') then
                flag1=.true.
                read(PARAFILE,*) max3D
            else if(char=='numGrains') then
                flag1=.true.
                read(PARAFILE,*) numGrains
            else if(char=='numSims') then
                flag1=.true.
                read(PARAFILE,*) numSims
            else
                write(*,*) 'error read unrecognized parameter: ', char
            end if
        end do
        flag1=.false.
    end do
    flag=.false.

    !*******************************************************
    !<read in anneal parameters
    !*******************************************************
    !<Anneal parameters
    annealTime = 0d0
    annealTemperature = 273d0
    annealType = 'constant'
    annealTempInc = 0d0
    annealSteps = 1 !<when annealType = 'constant'

    do while(flag .eqv. .FALSE.)
        read(PARAFILE,*) char
        if(char=='AnnealStart') then
            flag=.TRUE.
        end if
    end do
    flag=.FALSE.

    do while(flag .eqv. .FALSE.)
        flag1=.FALSE.
        do while(flag1 .eqv. .FALSE.)
            read(PARAFILE,*) char
            if(char=='AnnealEnd') then
                flag1=.TRUE.
                flag=.TRUE.
            else if(char=='annealTime') then
                flag1=.TRUE.
                read(PARAFILE,*) annealTime
            else if(char=='annealTemperature') then
                flag1=.TRUE.
                read(PARAFILE,*) annealTemperature
            else if(char=='annealType') then
                flag1=.TRUE.
                read(PARAFILE,*) annealType
            else if(char=='annealTempInc') then
                flag1=.TRUE.
                read(PARAFILE,*) annealTempInc
            else if(char=='annealSteps') then
                flag1=.TRUE.
                read(PARAFILE,*) annealSteps
            else
                write(*,*) 'error parameter: ', char
            end if
        end do
        flag1=.FALSE.
    end do
    flag=.FALSE.


    !*******************************************************
    !<read in output parameters
    !*******************************************************
    !<set default valuse for output parameters
    totdatToggle ='no'
    defectToggle ='no'
    stadatToggle = 'no'
    vtkToggle = 'no'
    minVoid = 10
    minLoop = 10
    minS = 10
    minSV = 10
    minBubble = 2

    do while(flag .eqv. .false.)
        read(PARAFILE,*) char
        if(char=='OutputStart') then
            flag=.true.
        end if
    end do
    flag=.false.

    do while(flag .eqv. .false.)
        flag1=.false.
        do while(flag1 .eqv. .false.)
            read(PARAFILE,*) char
            if(char=='OutputEnd') then
                flag1=.true.
                flag=.true.
            else if(char=='totdatToggle') then
                flag1=.true.
                read(PARAFILE,*) totdatToggle
            else if(char=='defectToggle') then
                flag1=.true.
                read(PARAFILE,*) defectToggle
            else if(char=='stadatToggle') then
                flag1=.true.
                read(PARAFILE,*) stadatToggle
            else if(char=='vtkToggle') then
                flag1=.true.
                read(PARAFILE,*) vtkToggle
            else if(char=='minLoop') then
                flag1=.true.
                read(PARAFILE,*) minLoop
            else if(char=='minVoid') then
                flag1=.true.
                read(PARAFILE,*) minVoid
            else if(char=='minS') then
                flag1=.true.
                read(PARAFILE,*) minS
            else if(char=='minSV') then
                flag1=.true.
                read(PARAFILE,*) minSV
            else if(char=='minBubble') then
                flag1=.true.
                read(PARAFILE,*) minBubble
            else
                write(*,*) 'error readParameters() unrecognized parameter: '
            end if
        end do
        flag1=.false.
    end do
    flag=.false.

    !***********************************************************************
    !<read in mesh parameters
    !***********************************************************************
    do while(flag .eqv. .FALSE.)
        read(PARAFILE,*) char
        if(char=='MeshStart') then
            flag=.TRUE.
        endif
    end do
    flag=.FALSE.

    do while(flag .eqv. .FALSE.)
        flag1=.FALSE.
        do while(flag1 .eqv. .FALSE.)
            read(PARAFILE,*) char
            if(char=='MeshEnd') then
                flag1=.TRUE.
                flag=.TRUE.
            else if(char=='length') then
                flag1=.TRUE.
                read(PARAFILE,*) meshLength
            else if(char=='numx') then
                flag1=.TRUE.
                read(PARAFILE,*) numxGlobal
            else if(char=='numy') then
                flag1=.TRUE.
                read(PARAFILE,*) numyGlobal
            else if(char=='numz') then
                flag1=.TRUE.
                read(PARAFILE,*) numzGlobal
        !    else if(char=='fineLength') then
        !        flag1=.TRUE.
        !        read(PARAFILE,*) fineLength
        !    else if(char=='numxFine') then
        !        flag1=.TRUE.
        !        read(PARAFILE,*) numxCascade
        !    else if(char=='numyFine') then
        !        flag1=.TRUE.
        !        read(PARAFILE,*) numyCascade
        !    else if(char=='numzFine') then
        !        flag1=.TRUE.
        !        read(PARAFILE,*) numzCascade
            else
                write(*,*) 'error'
            end if
        end do
        flag1=.FALSE.
    end do
    flag=.FALSE.

    close(PARAFILE)

    !*******************************************************
    !<Read other files
    !*******************************************************
    !<read *_Defects_*.txt
    call readDefectAttributes(defectFilename)

    !<read cpdf.*
    if(PKAspectrum=='yes') then
        call readPKAspectrum(pkaFilename)
    end if

    !<read cascade.txt
    if(irradiationType=='Cascade') then
        if(numCascadeFiles > 1) then
            call readCascadeFiles(cascadeFilename)
        else
            call readCascadeList(cascadeFilename)
        end if
    else if(irradiationType=='FrenkelPair') then
        numDisplacedAtoms=1.0
    end if

    !<read implant file
    if(implantType=='nonUniform') then
        call readImplantFile(implantFilename)
    end if

    !*******************************************************
    !<Calculate some constants
    !*******************************************************
    atomVol = lattice**3d0/2d0         !Atomic volume (nm^3)
    omega = (48d0*PI**2d0/atomVol**2d0)**(1d0/3d0)
    omega2D = (4d0*PI/(atomVol*burgers))**(1d0/2d0)
    omega1D = (9d0*PI/(16d0*atomVol))**(1d0/6d0)
    omegacircle1D = (1d0/burgers)**(1d0/2d0)

end subroutine ReadInputs

!***************************************************************************************************
!>Subroutine readDefectAttributes(filename). Information read in includes:
!   1) Number of allowed defect types
!   2) Formation energies, diffusion prefactors, migration energies and binding energies
!   3) Allowed reactions
!***************************************************************************************************
subroutine readDefectAttributes(filename)
    use mod_constants
    use mod_globalvariables
    implicit none

    character(len=100), intent(in) :: filename
    character(len=20) :: char
    integer :: i, j
    logical flag

    open(ATTRFILE, file=filename, status='old', action='read')
    flag=.false.

    !*******************************************************
    !<Read in formation energies
    !*******************************************************
    do while(flag .eqv. .false.)
        read(ATTRFILE,*) char
        if(char=='formationEnergies') then
            flag=.true.
        end if
    end do
    flag=.false.

    do while(flag .eqv. .false.)
        read(ATTRFILE,*) char
        if(char=='numSingle') then
            flag=.true.
            read(ATTRFILE,*) numFormSingle
        end if
    end do
    flag=.false.

    allocate(formSingle(numFormSingle))
    do i=1,numFormSingle
        allocate(formSingle(i)%defectType(SPECIES))
        read(ATTRFILE,*) (formSingle(i)%defectType(j),j=1,SPECIES)
        read(ATTRFILE,*) char, formSingle(i)%Ef
    end do

    !*******************************************************
    !<Read in diffusion prefactors and migration energies
    !*******************************************************
    do while(flag .eqv. .false.)
        read(ATTRFILE,*) char
        if(char=='diffusivities') then
            flag=.true.
        end if
    end do
    flag=.false.

    do while(flag .eqv. .false.)
        read(ATTRFILE,*) char
        if(char=='numSingle') then
            flag=.true.
            read(ATTRFILE,*) numDiffSingle
        end if
    end do
    flag=.false.

    allocate(diffSingle(numDiffSingle))
    do i=1,numDiffSingle
        allocate(diffSingle(i)%defectType(SPECIES))
        read(ATTRFILE,*) (diffSingle(i)%defectType(j),j=1,SPECIES)
        read(ATTRFILE,*) char, diffSingle(i)%D0, char, diffSingle(i)%Em
    end do

    do while(flag .eqv. .false.)
        read(ATTRFILE,*) char
        if(char=='numFunction') then
            flag=.true.
            read(ATTRFILE,*) numDiffFunc
        end if
    end do
    flag=.false.

    allocate(diffFunc(numDiffFunc))
    do i=1,numDiffFunc
        allocate(diffFunc(i)%defectType(SPECIES))
        read(ATTRFILE,*) (diffFunc(i)%defectType(j),j=1,SPECIES)
        allocate(diffFunc(i)%min(SPECIES))
        allocate(diffFunc(i)%max(SPECIES))
        read(ATTRFILE,*) char, (diffFunc(i)%min(j),j=1,SPECIES)
        read(ATTRFILE,*) char, (diffFunc(i)%max(j),j=1,SPECIES)
        read(ATTRFILE,*) char, diffFunc(i)%fType
        read(ATTRFILE,*) char, diffFunc(i)%numParam
        allocate(diffFunc(i)%para(diffFunc(i)%numParam))
        if(diffFunc(i)%numParam /= 0) then
            read(ATTRFILE,*) (diffFunc(i)%para(j),j=1,diffFunc(i)%numParam)
        end if
    end do

    !*******************************************************
    !<Read in binding energies
    !*******************************************************
    do while(flag .eqv. .false.)
        read(ATTRFILE,*) char
        if(char=='bindingEnergies') then
            flag=.true.
        end if
    end do
    flag=.false.

    do while(flag .eqv. .false.)
        read(ATTRFILE,*) char
        if(char=='numSingle') then
            flag=.true.
            read(ATTRFILE,*) numBindSingle
        end if
    end do
    flag=.false.

    allocate(bindSingle(numBindSingle))
    do i=1,numBindSingle
        allocate(bindSingle(i)%defectType(SPECIES))
        allocate(bindSingle(i)%product(SPECIES))
        read(ATTRFILE,*) (bindSingle(i)%defectType(j),j=1,SPECIES), (bindSingle(i)%product(j),j=1,SPECIES)
        read(ATTRFILE,*) char, bindSingle(i)%Eb
    end do

    do while(flag .eqv. .false.)
        read(ATTRFILE,*) char
        if(char=='numFunction') then
            flag=.true.
            read(ATTRFILE,*) numBindFunc
        end if
    end do
    flag=.false.

    allocate(bindFunc(numBindFunc))
    do i=1,numBindFunc
        allocate(bindFunc(i)%defectType(SPECIES))
        allocate(bindFunc(i)%product(SPECIES))
        read(ATTRFILE,*) (bindFunc(i)%defectType(j),j=1,SPECIES), (bindFunc(i)%product(j),j=1,SPECIES)
        allocate(bindFunc(i)%min(SPECIES))
        allocate(bindFunc(i)%max(SPECIES))
        read(ATTRFILE,*) char, (bindFunc(i)%min(j),j=1,SPECIES)
        read(ATTRFILE,*) char, (bindFunc(i)%max(j),j=1,SPECIES)
        read(ATTRFILE,*) char, bindFunc(i)%fType
        read(ATTRFILE,*) char, bindFunc(i)%numParam
        allocate(bindFunc(i)%para(bindFunc(i)%numParam))
        if(bindFunc(i)%numParam /= 0) then
            read(ATTRFILE,*) (bindFunc(i)%para(j),j=1,bindFunc(i)%numParam)
        end if
    end do

    !*******************************************************
    !<Read in reactions
    !*******************************************************
    flag=.false.

    !<Read in implantation
    !do while(flag .eqv. .false.)
    !    read(ATTRFILE,*) char
    !    if(char=='Implantation') then
    !        flag=.true.
    !        read(ATTRFILE,*) numImplantReaction
    !    end if
    !end do
    !flag=.false.
    !allocate(implantReactions(numImplantReaction))

    !do i=1,numImplantReaction
    !    do while(flag .eqv. .false.)
    !        read(ATTRFILE,*) char
    !        if(char=='FrenkelPair') then
    !            flag=.TRUE.
    !        else if(char=='Cascade') then
    !            flag=.TRUE.
    !        else if(char=='HeImplant') then
    !            flag=.TRUE.
    !        else
    !            write(*,*) 'error numImplantReac1'
    !            flag=.TRUE.
    !        end if
    !    end do
    !    flag = .false.

    !    if(char=='FrenkelPair') then
    !        implantReactions(i)%numReactants=0
    !        implantReactions(i)%numProducts=2
    !        allocate(implantReactions(i)%products(SPECIES,2))
    !        read(ATTRFILE,*) (implantReactions(i)%products(j,1),j=1,SPECIES),&
    !                (implantReactions(i)%products(j,2),j=1,SPECIES)
    !        read(ATTRFILE,*) char,implantReactions(i)%fType
    !    else if(char=='Cascade') then
    !        implantReactions(i)%numReactants=-10
    !        implantReactions(i)%numProducts=0
    !        read(ATTRFILE,*) char,implantReactions(i)%fType
    !    else if(char=='HeImplant') then
    !        implantReactions(i)%numReactants=0
    !        implantReactions(i)%numProducts=1
    !        allocate(implantReactions(i)%products(SPECIES,1))
    !        read(ATTRFILE,*) (implantReactions(i)%products(j,1),j=1,SPECIES)
    !        read(ATTRFILE,*) char,implantReactions(i)%fType
    !    else
    !        write(*,*) 'error numImplantReac'
    !    end if
    !end do

    !<Read in dissociation
    !do while(flag .eqv. .false.)
    !    read(ATTRFILE,*) char
    !    if(char=='Dissociation') then
    !        flag=.true.
    !        read(ATTRFILE,*) numDissocReaction
    !    end if
    !end do
    !flag=.false.
    !allocate(dissocReactions(numDissocReaction))

    !do i=1,numDissocReaction
    !    dissocReactions(i)%numReactants=1
    !    dissocReactions(i)%numProducts=1
    !    allocate(dissocReactions(i)%reactants(SPECIES,1))
    !    allocate(dissocReactions(i)%products(SPECIES,1))
    !    read(ATTRFILE,*) (dissocReactions(i)%reactants(j,1),j=1,SPECIES),&
    !            (dissocReactions(i)%products(j,1),j=1,SPECIES)
    !    allocate(dissocReactions(i)%min(SPECIES))
    !    allocate(dissocReactions(i)%max(SPECIES))
    !    read(ATTRFILE,*) char, (dissocReactions(i)%min(j),j=1,SPECIES)
    !    read(ATTRFILE,*) char, (dissocReactions(i)%max(j),j=1,SPECIES)
    !    read(ATTRFILE,*) char, dissocReactions(i)%fType
    !end do

    !<Read in SinkRemoval
    !do while(flag .eqv. .false.)
    !    read(ATTRFILE,*) char
    !    if(char=='SinkRemoval') then
    !        flag=.true.
    !        read(ATTRFILE,*) numSinkReaction
    !    end if
    !end do
    !flag=.false.
    !allocate(sinkReactions(numSinkReaction))

    !do i=1,numSinkReaction
    !    sinkReactions(i)%numReactants=1
    !    sinkReactions(i)%numProducts=0
    !    allocate(sinkReactions(i)%reactants(SPECIES,1))
    !    read(ATTRFILE,*) (sinkReactions(i)%reactants(j,1),j=1,SPECIES)
    !    allocate(sinkReactions(i)%min(SPECIES))
    !    allocate(sinkReactions(i)%max(SPECIES))
    !    read(ATTRFILE,*) char, (sinkReactions(i)%min(j),j=1,SPECIES)
    !    read(ATTRFILE,*) char, (sinkReactions(i)%max(j),j=1,SPECIES)
    !    read(ATTRFILE,*) char, sinkReactions(i)%fType
    !end do

    !<Read in Diffusion
    !do while(flag .eqv. .false.)
    !    read(ATTRFILE,*) char
    !    if(char=='Diffusion') then
    !        flag=.true.
    !        read(ATTRFILE,*) numDiffReaction
    !    end if
    !end do
    !flag=.false.
    !allocate(diffReactions(numDiffReaction))

    !do i=1,numDiffReaction
    !    diffReactions(i)%numReactants=1
    !    diffReactions(i)%numProducts=1
    !    allocate(diffReactions(i)%reactants(SPECIES,1))
    !    allocate(diffReactions(i)%products(SPECIES,1))
    !    read(ATTRFILE,*) (diffReactions(i)%reactants(j,1),j=1,SPECIES),&
    !            (diffReactions(i)%products(j,1),j=1,SPECIES)
    !    allocate(diffReactions(i)%min(SPECIES))
    !    allocate(diffReactions(i)%max(SPECIES))
    !    read(ATTRFILE,*) char, (diffReactions(i)%min(j),j=1,SPECIES)
    !    read(ATTRFILE,*) char, (diffReactions(i)%max(j),j=1,SPECIES)
    !    read(ATTRFILE,*) char, diffReactions(i)%fType
    !end do

    !<Read in Clustering
    do while(flag .eqv. .false.)
        read(ATTRFILE,*) char
        if(char=='Clustering') then
            flag=.true.
            read(ATTRFILE,*) numClusterReaction
        end if
    end do
    flag=.false.
    allocate(clusterReactions(numClusterReaction))

    do i=1,numClusterReaction
        clusterReactions(i)%numReactants=2
        clusterReactions(i)%numProducts=1
        allocate(clusterReactions(i)%reactants(SPECIES,2))
        read(ATTRFILE,*) (clusterReactions(i)%reactants(j,1),j=1,SPECIES),&
                (clusterReactions(i)%reactants(j,2),j=1,SPECIES)
        allocate(clusterReactions(i)%min(SPECIES*2))
        allocate(clusterReactions(i)%max(SPECIES*2))
        read(ATTRFILE,*) char, (clusterReactions(i)%min(j),j=1,SPECIES*2)
        read(ATTRFILE,*) char, (clusterReactions(i)%max(j),j=1,SPECIES*2)
        read(ATTRFILE,*) char, clusterReactions(i)%fType
    end do

    close(ATTRFILE)

end subroutine readDefectAttributes

!***************************************************************************************************
!> Subroutine readCascadeList(filename): reads defects in a cascade file
!This subroutine reads a list of cascades from a cascade file and stored in derived type cascade.
!***************************************************************************************************
subroutine readCascadeList(filename)
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none

    character(len=100), intent(in) :: filename
    character(len=20) :: char
    integer :: i, j, k
    type(cascade), pointer :: cascadeCurrent
    type(cascadeDefect), pointer :: defectCurrent

    open(CASFILE, file=filename, status='old', action='read')

    allocate(cascadeList)
    nullify(cascadeList%next)
    nullify(cascadeCurrent)
    cascadeCurrent=>cascadeList
    nullify(cascadeCurrent%next)
    read(CASFILE,*) PKAtemperature
    read(CASFILE,*) PKAenergy
    read(CASFILE,*) numDisplacedAtoms
    read(CASFILE,*) numCascades

    do i=1,numCascades
        read(CASFILE,*)
        read(CASFILE,*) cascadeCurrent%numDefects
        read(CASFILE,*) cascadeCurrent%numDisplacedAtoms
        allocate(cascadeCurrent%listDefects)
        nullify(cascadeCurrent%listDefects%next)
        nullify(defectCurrent)
        defectCurrent=>cascadeCurrent%listDefects
        nullify(defectCurrent%next)

        do j=1,cascadeCurrent%numDefects
            allocate(defectCurrent%defectType(SPECIES))
            read(CASFILE,*) (defectCurrent%defectType(k),k=1,SPECIES)
            read(CASFILE,*) (defectCurrent%coordinates(k), k=1,3)

            if(j /= cascadeCurrent%numDefects) then
                allocate(defectCurrent%next)
                defectCurrent=>defectCurrent%next
                nullify(defectCurrent%next)
            end if
        end do

        if(i /= numCascades) then
            allocate(cascadeCurrent%next)
            cascadeCurrent=>cascadeCurrent%next
            nullify(cascadeCurrent%next)
        end if
    end do

    close(CASFILE)

end subroutine readCascadeList

!***************************************************************************************************
!> Subroutine readCascadeFiles(filename): reads cascade defects in a number of cascade files
!This subroutine reads a list of cascades from each cascade file and stored in derived type cascade.
!***************************************************************************************************
subroutine readCascadeFiles(filename)
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none

    character(len=100), intent(in) :: filename       !< filename = '../../inputs/cascades/Fe_*.txt'
    character(len=100) :: systemChars
    character(len=100), allocatable :: cascadeFilename(:)
    integer :: fileID, stat, totalDisAtoms, totalCascades, i, j, k
    type(cascade), pointer :: casCurrent
    type(cascadeDefect), pointer :: casDef

    allocate(cascadeLists(numCascadeFiles))
    allocate(cascadeFilename(numCascadeFiles))
    totalDisAtoms=0
    totalCascades=0

    !! by wdj
    !systemChars = "ls -R "//trim(filename)// "|awk -F_ '{print $0,int($NF)}'|sort -n -k2|awk '{print $1}'  > cas.dat"
    !call system(trim(systemChars))
    !call MPI_Barrier(comm, ierr)

    systemChars = 'ls -R '//trim(filename)//' > cas.dat'
    call system(trim(systemChars))
    open(20, file='cas.dat')
    do fileID=1, numCascadeFiles
        read(20, fmt='(A)', iostat=stat) cascadeFilename(fileID)
        if(stat /= 0) then
            write(*,*) 'Error: read cascade filename'
        end if
    end do
    !close(20, status='delete')      !<close  and delete cas.dat
    close(20)

    do fileID=1, numCascadeFiles
        open(CASFILE, file=cascadeFilename(fileID), status='old', action='read')
        !<read number of atoms displaced and number of cascades in this cascade file
        read(CASFILE,*) cascadeLists(fileID)%temperature   !(K)
        read(CASFILE,*) cascadeLists(fileID)%PKAenergy     !(eV)
        read(CASFILE,*) cascadeLists(fileID)%averDisAtoms
        read(CASFILE,*) cascadeLists(fileID)%numCascades
        totalDisAtoms=totalDisAtoms+cascadeLists(fileID)%averDisAtoms
        totalCascades=totalCascades+cascadeLists(fileID)%numCascades
        allocate(cascadeLists(fileID)%listCascades)

        nullify(casCurrent)
        casCurrent=>cascadeLists(fileID)%listCascades
        nullify(casCurrent%next)

        do i=1, cascadeLists(fileID)%numCascades
            read(CASFILE,*)
            read(CASFILE,*) casCurrent%numDefects
            read(CASFILE,*) casCurrent%numDisplacedAtoms
            allocate(casCurrent%listDefects)
            nullify(casDef)
            casDef=>casCurrent%listDefects
            nullify(casDef%next)

            do j=1,casCurrent%numDefects
                allocate(casDef%defectType(SPECIES))
                read(CASFILE,*) (casDef%defectType(k),k=1,SPECIES)
                read(CASFILE,*) (casDef%coordinates(k), k=1,3)

                if(j /= casCurrent%numDefects) then
                    allocate(casDef%next)
                    casDef=>casDef%next
                    nullify(casDef%next)
                end if
            end do

            if(i /= cascadeLists(fileID)%numCascades) then
                allocate(casCurrent%next)
                casCurrent=>casCurrent%next
                nullify(casCurrent%next)
            end if
        end do

        close(CASFILE)
    end do

    numDisplacedAtoms=dble(totalDisAtoms)/dble(totalCascades)
    !! by wdj
    !numDisplacedAtoms = 4.88 * 1000 * 0.8 / 2 / 128	! NRT model, if adding 100、150、180keV using 4.64757keV

    deallocate(cascadeFilename)

end subroutine

!***************************************************************************************************
!> Subroutine readEPKAs(filename): reads PKA spectrum file
!The first column is the PKA energy, and the second column is the cpdf
!***************************************************************************************************
subroutine readPKAspectrum(filename)
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none

    character(len=100), intent(in) :: filename       !< filename = '../../inputs/pkas/cpdf.w
    integer :: i

    open(CPDFFILE, file=filename, status='old', action='read')
    read(CPDFFILE, *) EPKAlist%size
    read(CPDFFILE, *)
    allocate(EPKAlist%energy(EPKAlist%size))
    allocate(EPKAlist%cpdf(EPKAlist%size))

    do i=1, EPKAlist%size
        read(CPDFFILE, *) EPKAlist%energy(i), EPKAlist%cpdf(i)
    end do
    close(CPDFFILE)

end subroutine

!***************************************************************************************************
!> Subroutine readImplantFile(filename): read nonuniform DPA implantation profile
!The first column is the z_coord, and the second column is the dpa rate
!***************************************************************************************************
subroutine readImplantFile(filename)
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none

    character(len=100), intent(in) :: filename       !< filename = '../../inputs/dpas/DPAProfile_*.txt
    integer :: i, j, numColumn

    open(IMPLFILE, file=filename, status='old', action='read')
    read(IMPLFILE, *) numColumn
    read(IMPLFILE, *) numImpDatas
    read(IMPLFILE, *)
    allocate(impDPARates(numColumn,numImpDatas))

    do i=1, numImpDatas
        read(IMPLFILE, *) (impDPARates(j,i), j=1,numColumn)
    end do
    close(IMPLFILE)

end subroutine