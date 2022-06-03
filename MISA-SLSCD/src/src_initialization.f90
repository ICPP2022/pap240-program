!*****************************************************************************************
!>Subroutine initialize random seeds - seeds random number generator differently for each processor
!Uses the computer clock to initialize the random seed of the master processor.
!Then generates several integers in the master processor and uses them to initialize the random number seed of other processors.
!*****************************************************************************************
subroutine initializeRandomSeeds()
    use mod_globalvariables
    use mod_randdp
    implicit none
    include 'mpif.h'

    integer :: randseed, i, irand
    integer :: status(MPI_STATUS_SIZE)

    if(myProc%taskid == 0) then
        call system_clock(Count=randseed)	!return randseed (integer, unit:ms)
        call sdprnd(randseed)               !randseed is a input parameter

        do i=1,numtasks-1
            randseed=irand(randseed)
            call MPI_SEND(randseed, 1, MPI_INTEGER, i, 99,comm, ierr)
        end do
    else
        call MPI_RECV(randseed, 1, MPI_INTEGER, 0, 99, comm, status, ierr)
        call sdprnd(randseed)
    end if

end subroutine initializeRandomSeeds

!***************************************************************************************
!>Subroutine initialize number of point defects. (SIAs and Vacancies)
!***************************************************************************************
subroutine initializeNumDefects()
    use mod_constants
    use mod_globalvariables
    use mod_randdp
    implicit none
    include 'mpif.h'

    integer :: cell, i, maxNumTemp
    double precision :: totalAtoms
    double precision :: rtemp, r1

    totalAtoms=2d0*systemVolume/(lattice**3d0)    !Total atoms in the whole system
    numSevermesh=anint(soluteConc*totalAtoms/dble(totalMeshes)) !Number of Cu atoms in every mesh

    do i=1, numFormSingle
        if(formSingle(i)%defectType(1)>0 ) then	!SIA
            ceqSIA=dexp(-formSingle(i)%Ef/(KB*temperature))
        else if(formSingle(i)%defectType(1)<0) then	!V
            ceqVac=dexp(-formSingle(i)%Ef/(KB*temperature))
        end if
    end do

    if(initialVacancies == 0) then
        initNumVac=nint(ceqVac*totalAtoms)
    else
        initNumVac=initialVacancies
        !firr=dble(numVac)/systemVolume*atomVol
    end if
    if(initialSIAs == 0) then
        initNumSIA=nint(ceqSIA*totalAtoms)
    else
        initNumSIA=initialSIAs
    end if

    if(initNumVac >= initNumSIA) then
        maxNumTemp=initNumVac
    else
        maxNumTemp=initNumSIA
    end if

    allocate(initPointDefects(maxNumTemp,2))
    initPointDefects=0

    !<SIA
    rtemp = 0d0
    if(myProc%taskid==0) then
        outer1: do i=1, initNumSIA
            r1 = dprand()
            inter1: do cell=1,totalMeshes
                rtemp=rtemp+dble(cell)/dble(totalMeshes)
                if(r1 <= rtemp) then
                    initPointDefects(i,1)=cell
                    rtemp = 0d0
                    exit inter1
                end if
            end do inter1
        end do outer1
    end if

    !<V
    rtemp = 0d0
    if(myProc%taskid==0) then
        outer2: do i=1, initNumVac
            r1 = dprand()
            inter2: do cell=1, totalMeshes
                rtemp=rtemp+dble(cell)/dble(totalMeshes)
                if(r1 <= rtemp) then
                    initPointDefects(i,2)=cell
                    rtemp = 0d0
                    exit inter2
                end if
            end do inter2
        end do outer2
    end if

    if(initNumSIA > 0 .OR. initNumVac > 0) then
        call MPI_BCAST(initPointDefects, maxNumTemp*2, MPI_INTEGER, 0, comm,ierr)
    end if

end subroutine initializeNumDefects

!*****************************************************************************************
!>Subroutine compute_sinks(): compute total sink strength for all defects and were stored in the array sinks()
!sinks(1) for SIAs
!sinks(2) for vacancies
!*****************************************************************************************
subroutine compute_sinks()
    use mod_constants
    use mod_globalvariables
    implicit none

    double precision :: Sd, Sgbi, Sgbv

    Sd=dislocDensity
    if(grainBoundToggle=='yes') then
        Sgbi=6*dsqrt(ZI*Sd)/grainSize
        Sgbv=6*dsqrt(ZV*Sd)/grainSize
    else
        Sgbi=0d0
        Sgbv=0d0
    end if
    sinks(1)=ZI*Sd+Sgbi
    sinks(2)=ZV*Sd+Sgbv

end subroutine compute_sinks

!***********************************************************************
!> Subroutine initializeAnneal()
!This subroutine empties implantation(numMeshes)
!***********************************************************************
subroutine initializeAnneal()
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none

    integer :: cell

    do cell=1, numMeshes
        totalRateVol(cell) = totalRateVol(cell) - implantation(cell)%rate
        implantation(cell)%rate = 0d0
        if(HeDPAratio > 0d0) then
            totalRateVol(cell) = totalRateVol(cell) - implantation(cell)%next%rate
            implantation(cell)%next%rate = 0d0
        end if
    end do

end subroutine



