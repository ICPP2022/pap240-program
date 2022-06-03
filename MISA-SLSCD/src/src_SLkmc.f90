!*****************************************************************************************
!>double precision timeStepGenerate(): chooses a timestep using random number and Monte Carlo algorithm
!(this is a global timestep)
!*****************************************************************************************
double precision function compute_thresholdTime()
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none
    include 'mpif.h'

    integer :: ntmp(8), i, j, subCell, cell, numCells, level, temp
    double precision :: nstop, ptmp(8), pmax, pmaxall, dt_kmc
    type(defect), pointer :: def
    type(reaction), pointer :: reac
    double precision :: commTime1, commTime2
    integer, external :: getLocalCell

    !nstop = 10.0  !nstop = (number of active defects in the sector) / (total defects in the sector)
    ntmp = 0      !Number of defects in the sector
    ptmp = 0d0    !The average reaction rate in a sector
    pmax = 0d0    !Maximum ptmp in the current process--max(ptmp(1)~ptmp(SECTORS))
    pmaxall = 0d0 !Maximum pmax in the whole system

    do i=1, 8
        numCells = mySector(i)%numCell(1)*mySector(i)%numCell(2)*mySector(i)%numCell(3)
        do subCell=1, numCells
            cell = getLocalCell(subCell, i)
            !implantation reaction
            nullify(reac)
            reac => implantation(cell)
            do while(associated(reac))
                if(reac%rate > 0d0) then
                    ntmp(i) = ntmp(i) + 1
                    ptmp(i) = ptmp(i) + reac%rate
                end if
                reac => reac%next
            end do

            do level=1, LEVELS
                nullify(def)
                def => myDRL(level, cell)%defect_all
                do while(associated(def))
                    nullify(reac)
                    reac => def%reactionList
                    temp = 0
                    do while(associated(reac))
                        temp = temp + 1
                        ptmp(i) = ptmp(i) + reac%rate*dble(def%num)
                        reac => reac%next
                    end do
                    ntmp(i) = ntmp(i) + def%num*temp
                    def => def%next
                end do
            end do
        end do
        ptmp(i) = ptmp(i) / dble(ntmp(i))   !当前扇区上，平均每个事件的反应速率
        pmax = max(ptmp(i), pmax)   !当前进程上的最大的平均反应速率
    end do

    commTime1 = MPI_WTIME()
    call MPI_ALLREDUCE(pmax, pmaxall, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
    commTime2 = MPI_WTIME()
    allreduceTime = allreduceTime + (commTime2-commTime1)

    if(pmaxall > 0d0) then
        !dt_kmc = nstop/pmaxall
        dt_kmc = NEVENT/pmaxall
    else
        dt_kmc = stopTime-time
    end if
    compute_thresholdTime = min(dt_kmc, stopTime-time)

end function

!*****************************************************************************************
!>Function SL_timeStepGenerate(sector)
!Chooses a timestep using a random number and Monte Carlo algorithm (this is a local timestep)
!*****************************************************************************************
double precision function SL_timeStepGenerate(sector, totRateSector)
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    use mod_randdp
    implicit none

    integer, intent(in) :: sector
    double precision, intent(inout) :: totRateSector
    integer :: subCell, cell, numCells
    double precision :: r1
    integer, external :: getLocalCell

    !<First: re-compute total rate of this sector
    !mySector(sector)%totalRate = 0d0
    totRateSector = 0d0
    numCells = mySector(sector)%numCell(1)*mySector(sector)%numCell(2)*mySector(sector)%numCell(3)
    do subCell=1, numCells
        cell = getLocalCell(subCell, sector)
        !mySector(sector)%totalRate = mySector(sector)%totalRate + totalRateVol(cell)
        totRateSector = totRateSector + totalRateVol(cell)
    end do

    r1 = dprand()
    !SL_timeStepGenerate = dlog(1d0/r1) / mySector(sector)%totalRate
    SL_timeStepGenerate = dlog(1d0/r1)/totRateSector

end function

!*****************************************************************************************
!>subroutine chooseReaction(): chooses a reaction in each processor according to the MC algorithm
!Note: choose a mesh first and then choose a reaction within that mesh.
!拟优化：为每个def增加%totalRate，先选cell，再选def，最后再选reaction
!*****************************************************************************************
subroutine SL_chooseReaction(sector, totRateSector, defCH, reacCH, cellCH)
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    use mod_randdp
    implicit none

    !<formal parameters
    integer, intent(in) :: sector
    double precision, intent(in) :: totRateSector
    type(defect), pointer, intent(inout) :: defCH
    type(reaction), pointer, intent(inout) :: reacCH
    integer, intent(inout) :: cellCH
    !<local variables
    type(defect), pointer :: def
    double precision :: r2, r2times, atemp, atemp_cell
    integer :: subCell, cell, j, level, numCells
    integer, external :: getLocalCell

    r2 = dprand()
    atemp = 0d0
    atemp_cell = 0d0
    !r2times = r2 * mySector(sector)%totalRate
    r2times = r2 * totRateSector
    cellCH = 0
    nullify(reacCH)
    nullify(defCH)
    numCells = mySector(sector)%numCell(1)*mySector(sector)%numCell(2)*mySector(sector)%numCell(3)
    outer: do subCell=1, numCells
        cell = getLocalCell(subCell, sector)
        atemp_cell = atemp_cell + totalRateVol(cell)
        if(r2times <= atemp_cell) then    !<choose a cell
            !<First: choose 0nd reaction
            reacCH => implantation(cell)
            do while(associated(reacCH))
                atemp = atemp + reacCH%rate
                if(r2times <= atemp) then
                    cellCH = cell
                    exit outer
                end if
                reacCH => reacCH%next
            end do

            !<Second: if not choose a reaction, try to choose reaction from 1st, diff, 2nd reactions
            if(r2times > atemp) then
                do level=1, LEVELS
                    nullify(defCH)
                    defCH => myDRL(level,cell)%defect_all
                    do while(associated(defCH))
                        reacCH => defCH%reactionList
                        do while(associated(reacCH))
                            atemp = atemp + reacCH%rate
                            if(r2times <= atemp) then
                                cellCH = cell
                                exit outer
                            end if
                            reacCH => reacCH%next
                        end do
                        defCH => defCH%next
                    end do
                end do
            end if
        else
            atemp = atemp + totalRateVol(cell)
        end if
    end do outer

    if(associated(reacCH)) then !choosed a reaction
        if(reacCH%numReactants==0 .AND. reacCH%numProducts==2) then
            numImpEvents = numImpEvents + 1
        else if(reacCH%numReactants==-10 .AND. reacCH%numProducts==0) then
            numImpEvents = numImpEvents + 1
        else if(reacCH%numReactants==0 .AND. reacCH%numProducts==1) then    !He implantation
            numImpHe = numImpHe + 1
        end if
    end if

    !test test test
    if(associated(reacCH)) then
        do j=1, 6
            if(myMesh(cellCH)%neighborProc(j)/=myProc%taskid .AND. myMesh(cellCH)%neighborProc(j)/=-1) then
                countEvent = countEvent + 1
                exit
            end if
        end do
    end if

end subroutine

!***********************************************************************************************************************
!> Subroutine updateDefectList(reacCH, defectUpdateCurrent):
!Function:  It updates defects in the local mesh (and ghost mesh) according to the reaction chosen, and add the defects
!that may have passed into neighboring processor as well as defects that may have changed in ghost to defCommList.
!It also creates a list of defects that have been updated, to inform the SL_updateReactionList(defUpdateList) subroutine.
!NOTE: defectUpdateCurrent%num = +/- 1, indicating adding or removing one of these defects.
!Date: 2019/12/25
!Author: CDD
!Amendant Record: original code
!***********************************************************************************************************************
subroutine SL_updateDefectList(defCH,reacCH,cellCH,numUpdates,defUpdateStep,casCellList,defCommList)
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    use mod_randdp
    implicit none

    !<formal parameters
    type(defect), pointer, intent(inout) :: defCH
    type(reaction), pointer, intent(in) :: reacCH
    integer, intent(in) :: cellCH, numUpdates
    integer, intent(inout), dimension(SPECIES+6,numUpdates) :: defUpdateStep
    logical, intent(inout), dimension(numMeshes) :: casCellList
    type(defectCommHead), intent(inout), dimension(3) :: defCommList
    !<local variables
    type(defect), pointer :: def,def_prev, def2     !point to myDRL(level,cell)%defect_all
    type(defect), pointer :: def_m,def_mPrev         !point to myDRL(level,cell)%defect_mobile
    type(cascade), pointer :: cascadeChoosed
    integer :: defectType(SPECIES), numProduct
    integer :: numSame, dir, i, i1, j, k, level, dim, index
    double precision :: diff, diffRandom
    !<Used for communication between processors
    type(defectComm), pointer :: defComm, defComm_prev
    !<Functions
    double precision, external :: compute_diffusivity
    integer, external :: find_level, find_index

    interface
        subroutine chooseCascade_multiFiles(cascadeChoosed)
            use mod_structuretype
            use mod_globalvariables
            use mod_randdp
            implicit none
            type(cascade), pointer, intent(inout) :: cascadeChoosed
        end subroutine

        subroutine chooseCascade(cascadeChoosed)
            use mod_structuretype
            use mod_globalvariables
            use mod_randdp
            implicit none
            type(cascade), pointer, intent(inout) :: cascadeChoosed
        end subroutine

        subroutine locate_defect(defectType, def, def_prev, numSame)
            use mod_constants
            use mod_structuretype
            implicit none
            integer, intent(in), dimension(SPECIES) :: defectType
            type(defect), pointer, intent(inout) :: def, def_prev
            integer, intent(inout) :: numSame
        end subroutine

        subroutine cascadeCombine_one(cascadeChoosed, cellCH)
            use mod_constants
            use mod_structuretype
            use mod_globalvariables
            use mod_randdp
            implicit none
            type(cascade), pointer, intent(in) :: cascadeChoosed    !default: chosed a cascade
            integer, intent(in) :: cellCH
        end subroutine

        subroutine cascadeCombine_two(cascadeChoosed, cellCH)
            use mod_constants
            use mod_structuretype
            use mod_globalvariables
            use mod_randdp
            implicit none
            type(cascade), pointer, intent(in) :: cascadeChoosed    !default: chosed a cascade
            integer, intent(in) :: cellCH
        end subroutine
    end interface

    defUpdateStep = 0
    if(associated(reacCH)) then    !<we have choosen a reaction
        if(reacCH%numReactants == -10) then
            !******************************************
            !<Mark the cell that was implanted a cascade and delete defects that are in cascadeCell
            !******************************************
            do dim=1, 3
                dirL1: do dir=dim*2-1, dim*2
                    if(myMesh(cellCH)%neighborProc(dir)/=myProc%taskid .AND. &
                            myMesh(cellCH)%neighborProc(dir)/=-1) then  !boundary cell
                        !<Mark the cell that was implanted a cascade
                        casCellList(cellCH)=.TRUE.
                        !<delete defects in defCommList with cell= cellCH
                        nullify(defComm)
                        nullify(defComm_prev)
                        defComm => defCommList(dim)%commHead
                        do while(associated(defComm))
                            if(defComm%cell==cellCH .AND. defComm%dir==0) then
                                if(.NOT. associated(defComm_prev)) then !defComm is the first node
                                    defCommList(dim)%commHead => defComm%next
                                    deallocate(defComm%defectType)
                                    deallocate(defComm)
                                    defComm => defCommList(dim)%commHead
                                else
                                    defComm_prev%next => defComm%next
                                    deallocate(defComm%defectType)
                                    deallocate(defComm)
                                    defComm => defComm_prev%next
                                end if
                            else
                                defComm_prev => defComm
                                defComm => defComm%next
                            end if
                        end do
                        exit dirL1
                    end if
                end do dirL1
            end do

            !******************************************
            !Choose a cascade
            !******************************************
            nullify(cascadeChoosed)
            do while(.NOT. associated(cascadeChoosed))
                if(numCascadeFiles > 1) then
                    call chooseCascade_multiFiles(cascadeChoosed)
                else
                    call chooseCascade(cascadeChoosed)
                end if
            end do

            call clearReactions(cellCH)

            !<combine
            call cascadeCombine_one(cascadeChoosed, cellCH)
            !call cascadeCombine_two(cascadeChoosed, cellCH)

            call reset_defectMobileList(cellCH)

        else    !chosed FrenkelPair/1st/diff/2nd reaction
            !***********************************************************************************************
            ! Defect update for reactions chosen in the coarse mesh.
            ! For defect changes in the ghost of other processors or in the ghost of this processor,
            ! local and ghost buffers are created for communication with other processors.
            !NOTE, reacCH%numReactants is real, but reacCH%numProducts may not be real for 1st reaction
            !***********************************************************************************************
            !<First: update reactants
            do i=1, reacCH%numReactants !%numReactants = 0(FrenkelPair) or 1(1st/diff) or 2(2nd)
                defectType = defCH%defectType
                if(i == 2) then !2nd eraction
                    defectType(1:SPECIES) = reacCH%reactants(1:SPECIES,1)
                end if

                !<Second: Add reactants to defUpdateStep(:,:)
                defUpdateStep(1:SPECIES,i) = defectType(1:SPECIES)
                defUpdateStep(SPECIES+1,i) = -1                 !defectUpdatecurrent%num
                defUpdateStep(SPECIES+2,i) = cellCH         !defectUpdateCurrent%cell
                defUpdateStep(SPECIES+3,i) = myProc%taskid      !defectUpdateCurrent%proc
                defUpdateStep(SPECIES+4,i) = 0                  !defectUpdateCurrent%dir
                defUpdateStep(SPECIES+5,i) = 0                  !defectUpdateCurrent%neighbor
                defUpdateStep(SPECIES+6,i) = 0                  !defectUpdateCurrent%cascadeID

                !<Third: update defects
                level = find_level(defectType)
                nullify(def)
                def => myDRL(level, cellCH)%defect_all
                nullify(def_prev)
                do while(associated(def))
                    numSame = 0
                    do j=1, SPECIES
                        if(def%defectType(j) == defectType(j)) then
                            numSame = numSame + 1
                        end if
                    end do
                    if(numSame == SPECIES) then !find this type of defect
                        exit
                    end if
                    def_prev => def
                    def => def%next
                end do

                if(.NOT. associated(def)) then
                    write(*,*) 'Tried to delete defect that wasnt there coarse'
                    write(*,*) 'numReactants', reacCH%numReactants, 'numProducts', reacCH%numProducts
                    write(*,*) 'cellCH', cellCH
                    write(*,*) 'defectType', defectType
                    write(*,*) 'step',step,'proc',myProc%taskid,'i',i
                    write(*,*) 'rate', reacCH%rate
                    call MPI_ABORT(comm,ierr)
                else    !associated(def), this defect is in myDRL
                    def%num = def%num - 1
                end if

                !<First: Update defCommList
                do dim=1, 3
                    dirL2: do dir=dim*2-1, dim*2
                        if(myMesh(cellCH)%neighborProc(dir)/=-1 .AND. &
                                myMesh(cellCH)%neighborProc(dir)/=myProc%taskid .AND. &
                                (casCellList(cellCH) .eqv. .false.)) then   !this defect should be sent to neighbor processor
                            if(def%diff > 0d0) then
                                nullify(defComm)
                                nullify(defComm_prev)
                                defComm => defCommList(dim)%commHead
                                inter2: do while(associated(defComm))
                                    numSame = 0
                                    if(defComm%cell==cellCH .AND. defComm%dir==0) then
                                        do j=1, SPECIES
                                            if(defComm%defectType(j) == defectType(j)) then
                                                numSame = numSame + 1
                                            end if
                                        end do
                                        if(numSame == SPECIES) then !find it
                                            defComm%num = defComm%num - 1
                                            exit inter2
                                        end if
                                    end if
                                    defComm_prev => defComm
                                    defComm => defComm%next
                                end do inter2

                                if(.NOT. associated(defComm)) then  !add this defect to defCommList(dim)
                                    allocate(defComm)
                                    allocate(defComm%defectType(SPECIES))
                                    defComm%defectType = defectType
                                    defComm%num = -1
                                    defComm%cell = cellCH
                                    defComm%neighbor = myMesh(cellCH)%neighbor(dir)   !local (ghost) id of neighbor (current) processor
                                    defComm%dir = 0   !defComm%dir=0: the defect is within the boundary of current processor
                                    nullify(defComm%next)
                                    if(.NOT. associated(defComm_prev)) then
                                        defCommList(dim)%commHead => defComm
                                    else
                                        defComm_prev%next => defComm
                                    end if
                                end if
                            end if
                            exit dirL2
                        end if
                    end do dirL2
                end do

                !update myDRL(level, cellCH)%defect_mobile
                if(def%diff > 0d0) then   !this defect is mobile
                    nullify(def_m)
                    def_m => myDRL(level, cellCH)%defect_mobile
                    nullify(def_mPrev)
                    do while(associated(def_m))
                        numSame = 0
                        do j=1, SPECIES
                            if(def_m%defectType(j) == defectType(j)) then
                                numSame = numSame + 1
                            end if
                        end do
                        if(numSame == SPECIES) then
                            exit
                        end if
                        def_mPrev => def_m
                        def_m => def_m%next
                    end do
                    if(.NOT. associated(def_m)) then
                        write(*,*) 'Tried to delete defect that was not the defect_mobile'
                        write(*,*) 'numReactants', reacCH%numReactants, 'numProducts', reacCH%numProducts
                        write(*,*) 'cellCH', cellCH
                        write(*,*) 'defectType', defectType
                        write(*,*) 'step',step,'proc',myProc%taskid,'i',i
                        write(*,*) 'rate', reacCH%rate
                        call MPI_ABORT(comm,ierr)
                    else if(def_m%num == 1) then    !delete this defect
                        if(.NOT. associated(def_mPrev)) then     !def_m is the first node in the list
                            myDRL(level, cellCH)%defect_mobile => def_m%next
                        else
                            def_mPrev%next => def_m%next !remove that defect type from the system
                        end if
                        deallocate(def_m%defectType)
                        deallocate(def_m)
                    else
                        def_m%num = def_m%num - 1
                    end if
                end if
            end do

            !***********************************************************************
            !If a defect migrates from one cell to another, the defect will be given a percent chance
            !of removal from the system based on the cell size (%chance of removal = cell size/grain size)
            !This is to replicate the OKMC practice of removing SIA clusters after they migrate 1 um.
            !***********************************************************************
            if(grainBoundToggle=='yes' .AND. reacCH%numReactants==1 .AND. reacCH%numProducts==1) then
                diffRandom = dprand()
                if(diffRandom <= meshLength/grainSize) then
                    numProduct = 0  !这步使得defUpdateStep(:,2) = 0
                else
                    numProduct = reacCH%numProducts
                end if
            else
                numProduct = reacCH%numProducts
            end if

            !<Second: update products, numUpdates=1/2/3/4, %numReactants=0/1/2
            !do i=1, reacCH%numProducts
            do i=1, numProduct
                !get defectType
                defectType = 0
                if(reacCH%numReactants==0 .AND. reacCH%numProducts==2) then !FP
                    if(i == 1) then
                        defectType(1) = 1
                    else
                        defectType(1) = -1
                    end if
                else if(reacCH%numReactants==0 .AND. reacCH%numProducts==1) then    !He implant
                    defectType(2) = 1   !He atom
                else if(reacCH%numReactants==1 .AND. reacCH%numProducts==2) then    !diss
                    defectType(1:SPECIES) = reacCH%products(1:SPECIES,1)    !pointDefect
                    if(i == 2) then
                        if(SIAPinToggle=='yes' .AND. defCH%defectType(SPECIES)>0) then  !defCH is im_SIAs
                            defectType = defCH%defectType
                            defectType(SPECIES) = defectType(SPECIES) - 1
                            if(defectType(SPECIES) <= max3D) then
                                defectType(1) = defectType(SPECIES)
                                defectType(SPECIES) = 0
                            end if
                        else
                            defectType = defCH%defectType - defectType
                        end if
                    end if
                else if(reacCH%numReactants==1 .AND. reacCH%numProducts==1) then    !diff
                    defectType = defCH%defectType
                else    !2nd
                    defectType(1:SPECIES) = reacCH%products(1:SPECIES,i)
                end if

                diff = compute_diffusivity(defectType)
                i1 = i + reacCH%numReactants

                if(reacCH%numReactants==1 .and. reacCH%numProducts==1) then !chosed diff reaction
                    if(reacCH%taskid(1) == myProc%taskid) then  !product is in the local cell
                        !<First: update defCommList
                        do dim=1, 3
                            dirL3: do dir=dim*2-1, dim*2
                                if(myMesh(reacCH%cell(1))%neighborProc(dir)/=-1 .AND. &
                                        myMesh(reacCH%cell(1))%neighborProc(dir)/=myProc%taskid .AND. &
                                        (casCellList(reacCH%cell(1)) .eqv. .FALSE.)) then
                                    if(diff > 0d0) then
                                        nullify(defComm)
                                        nullify(defComm_prev)
                                        defComm => defCommList(dim)%commHead
                                        inter3: do while(associated(defComm))
                                            numSame = 0
                                            if(defComm%cell==reacCH%cell(1) .AND. defComm%dir==0) then
                                                do j=1, SPECIES
                                                    if(defComm%defectType(j) == defectType(j)) then
                                                        numSame = numSame + 1
                                                    end if
                                                end do
                                                if(numSame == SPECIES) then   !find it
                                                    defComm%num = defComm%num + 1
                                                    exit inter3
                                                end if
                                            end if
                                            defComm_prev => defComm
                                            defComm => defComm%next
                                        end do inter3
                                        if(.NOT. associated(defComm)) then
                                            allocate(defComm)
                                            allocate(defComm%defectType(SPECIES))
                                            defComm%defectType = defectType
                                            defComm%num = 1
                                            defComm%cell = reacCH%cell(1)    !local cell
                                            defComm%neighbor = myMesh(reacCH%cell(1))%neighbor(dir)  !ghost cell of current processor
                                            defComm%dir = 0
                                            nullify(defComm%next)
                                            if(.NOT. associated(defComm_prev)) then
                                                defCommList(dim)%commHead => defComm
                                            else
                                                defComm_prev%next => defComm
                                            end if
                                        end if
                                    end if
                                    exit dirL3
                                end if
                            end do dirL3
                        end do

                        !<Second: add products to defUpdateStep(:,:)
                        defUpdateStep(1:SPECIES,i1) = defectType(1:SPECIES)
                        defUpdateStep(SPECIES+1,i1) = 1               !defectUpdatecurrent%num
                        defUpdateStep(SPECIES+2,i1) = reacCH%cell(1)      !defectUpdateCurrent%cell
                        defUpdateStep(SPECIES+3,i1) = myProc%taskid   !defectUpdateCurrent%proc
                        defUpdateStep(SPECIES+4,i1) = 0               !defectUpdateCurrent%dir
                        defUpdateStep(SPECIES+5,i1) = 0               !defectUpdateCurrent%neighbor
                        defUpdateStep(SPECIES+6,i1) = 0               !defectUpdateCurrent%cascadeID

                        !<Third: uodate products in %defec_all
                        level = find_level(defectType)
                        nullify(def)
                        def => myDRL(level, reacCH%cell(1))%defect_all
                        call locate_defect(defectType, def, def_prev, numSame)
                        if(associated(def)) then
                            if(numSame == SPECIES) then   !find it
                                def%num = def%num + 1
                            else    !insert the defect in front of def
                                nullify(def2)
                                allocate(def2)
                                allocate(def2%defectType(SPECIES))
                                def2%defectType = defectType
                                def2%num = 1
                                def2%diff = diff
                                nullify(def2%reactionList)
                                def2%next => def
                                if(.NOT. associated(def_prev)) then !def is the first node
                                    myDRL(level, reacCH%cell(1))%defect_all => def2
                                else
                                    def_prev%next => def2
                                end if
                            end if
                        else    !add the defect to the end of the list
                            nullify(def2)
                            allocate(def2)
                            allocate(def2%defectType(SPECIES))
                            def2%defectType = defectType
                            def2%num = 1
                            def2%diff = diff
                            nullify(def2%reactionList)
                            nullify(def2%next)
                            if(.NOT. associated(def_prev)) then !def2 will be the first node
                                myDRL(level, reacCH%cell(1))%defect_all => def2
                            else
                                def_prev%next => def2
                            end if
                        end if

                        !<Fourth: update %defect_mobile
                        if(diff > 0d0) then
                            nullify(def_m)
                            def_m => myDRL(level, reacCH%cell(1))%defect_mobile
                            call locate_defect(defectType, def_m, def_mPrev, numSame)
                            if(associated(def_m)) then
                                if(numSame == SPECIES) then
                                    def_m%num = def_m%num + 1
                                else    !insert the defect in front of def_m
                                    nullify(def2)
                                    allocate(def2)
                                    allocate(def2%defectType(SPECIES))
                                    def2%defectType = defectType
                                    def2%num = 1
                                    def2%diff = diff
                                    nullify(def2%reactionList)
                                    def2%next => def_m
                                    if(.NOT. associated(def_mPrev)) then !def_m is the first node
                                        myDRL(level, reacCH%cell(1))%defect_mobile => def2
                                    else
                                        def_mPrev%next => def2
                                    end if
                                end if
                            else
                                nullify(def2)
                                allocate(def2)
                                allocate(def2%defectType(SPECIES))
                                def2%defectType = defectType
                                def2%num = 1
                                def2%diff = diff
                                nullify(def2%reactionList)
                                nullify(def2%next)
                                if(.NOT. associated(def_mPrev)) then !def2 will be the first node
                                    myDRL(level, reacCH%cell(1))%defect_mobile => def2
                                else
                                    def_mPrev%next => def2
                                end if
                            end if
                        end if
                    else if(reacCH%taskid(1)/=myProc%taskid .AND. reacCH%taskid(1)/=-1) then    !!product is in the ghost cell
                        do dir=1, 6
                            if(myMesh(cellCH)%neighbor(dir)==reacCH%cell(1) .AND. &
                                    myMesh(cellCH)%neighborProc(dir)==reacCH%taskid(1)) then
                                exit
                            end if
                        end do

                        !<First: update defCommList
                        dim = dir/2 + mod(dir,2)
                        nullify(defComm)
                        nullify(defComm_prev)
                        defComm => defCommList(dim)%commHead
                        inter4: do while(associated(defComm))
                            numSame = 0
                            if(defComm%cell == reacCH%cell(1) .AND. defComm%dir == dir) then
                                do j=1, SPECIES
                                    if(defComm%defectType(j) == defectType(j)) then
                                        numSame = numSame + 1
                                    end if
                                end do
                                if(numSame == SPECIES) then
                                    defComm%num = defComm%num + 1
                                    exit inter4
                                end if
                            end if
                            defComm_prev => defComm
                            defComm => defComm%next
                        end do inter4
                        if(.NOT. associated(defComm)) then
                            allocate(defComm)
                            allocate(defComm%defectType(SPECIES))
                            defComm%defectType = defectType
                            defComm%num = 1
                            defComm%cell = reacCH%cell(1)    !ghost cell, not local cell
                            defComm%neighbor = -100   !a label, the defect is in ghost in the dir direction
                            defComm%dir = dir
                            nullify(defComm%next)
                            if(.NOT. associated(defComm_prev)) then
                                defCommList(dim)%commHead => defComm
                            else
                                defComm_prev%next => defComm
                            end if
                        end if

                        !<Second: add products to defUpdateStep(:,:)
                        defUpdateStep(1:SPECIES,i1) = defectType(1:SPECIES)
                        defUpdateStep(SPECIES+1,i1) = 1                                     !defectUpdatecurrent%num
                        defUpdateStep(SPECIES+2,i1) = reacCH%cell(1)      !defectUpdateCurrent%cell
                        defUpdateStep(SPECIES+3,i1) = reacCH%taskid(1)  !defectUpdateCurrent%proc
                        defUpdateStep(SPECIES+4,i1) = dir                                   !defectUpdateCurrent%dir
                        defUpdateStep(SPECIES+5,i1) = cellCH                            !defectUpdateCurrent%neighbor
                        defUpdateStep(SPECIES+6,i1) = 0                                     !defectUpdateCurrent%cascadeID

                        !<Third: update products in myGhost
                        level = find_level(defectType)
                        index = find_index(cellCH, dir)
                        nullify(def)
                        def => myGhost(dir)%myDRL_ghost(level,index)%defect_mobile
                        call locate_defect(defectType, def, def_prev, numSame)
                        if(associated(def)) then
                            if(numSame == SPECIES) then
                                def%num = def%num + 1
                            else
                                nullify(def2)
                                allocate(def2)
                                allocate(def2%defectType(SPECIES))
                                def2%defectType = defectType
                                def2%num = 1
                                def2%diff = diff
                                nullify(def2%reactionList)
                                def2%next => def
                                if(.NOT. associated(def_prev)) then
                                    myGhost(dir)%myDRL_ghost(level,index)%defect_mobile => def2
                                else
                                    def_prev%next => def2
                                end if
                            end if
                        else
                            nullify(def2)
                            allocate(def2)
                            allocate(def2%defectType(SPECIES))
                            def2%defectType = defectType
                            def2%num = 1
                            def2%diff = diff
                            nullify(def2%reactionList)
                            nullify(def2%next)
                            if(.NOT. associated(def_prev)) then
                                myGhost(dir)%myDRL_ghost(level,index)%defect_mobile => def2
                            else
                                def_prev%next => def2
                            end if
                        end if
                    end if
                else    !chosed other reaction
                    !<First: update defCommList
                    do dim=1, 3
                        dirL5: do dir=dim*2-1, dim*2
                            if(myMesh(cellCH)%neighborProc(dir)/=-1 .AND. &
                                    myMesh(cellCH)%neighborProc(dir)/=myProc%taskid .AND. &
                                    (casCellList(cellCH) .eqv. .FALSE.)) then
                                if(diff > 0d0) then
                                    nullify(defComm)
                                    nullify(defComm_prev)
                                    defComm => defCommList(dim)%commHead
                                    inter5: do while(associated(defComm))
                                        numSame = 0
                                        if(defComm%cell==cellCH .AND. defComm%dir==0) then
                                            do j=1, SPECIES
                                                if(defComm%defectType(j) == defectType(j)) then
                                                    numSame = numSame + 1
                                                end if
                                            end do
                                            if(numSame == SPECIES) then   !find it
                                                defComm%num = defComm%num + 1
                                                exit inter5
                                            end if
                                        end if
                                        defComm_prev => defComm
                                        defComm => defComm%next
                                    end do inter5
                                    if(.NOT. associated(defComm)) then
                                        allocate(defComm)
                                        allocate(defComm%defectType(SPECIES))
                                        defComm%defectType = defectType
                                        defComm%num = 1
                                        defComm%cell = cellCH    !local cell
                                        defComm%neighbor = myMesh(cellCH)%neighbor(dir)  !ghost cell of current processor
                                        defComm%dir = 0
                                        nullify(defComm%next)
                                        if(.NOT. associated(defComm_prev)) then
                                            defCommList(dim)%commHead => defComm
                                        else
                                            defComm_prev%next => defComm
                                        end if
                                    end if
                                end if
                                exit dirL5
                            end if
                        end do dirL5
                    end do

                    !<Second: add products to defUpdateStep(:,:)
                    defUpdateStep(1:SPECIES,i1) = defectType(1:SPECIES)
                    defUpdateStep(SPECIES+1,i1) = 1               !defectUpdatecurrent%num
                    defUpdateStep(SPECIES+2,i1) = cellCH      !defectUpdateCurrent%cell
                    defUpdateStep(SPECIES+3,i1) = myProc%taskid   !defectUpdateCurrent%proc
                    defUpdateStep(SPECIES+4,i1) = 0               !defectUpdateCurrent%dir
                    defUpdateStep(SPECIES+5,i1) = 0               !defectUpdateCurrent%neighbor
                    defUpdateStep(SPECIES+6,i1) = 0               !defectUpdateCurrent%cascadeID

                    !<Third: uodate products in %defec_all
                    level = find_level(defectType)
                    nullify(def)
                    def => myDRL(level, cellCH)%defect_all
                    call locate_defect(defectType, def, def_prev, numSame)
                    if(associated(def)) then
                        if(numSame == SPECIES) then   !find it
                            def%num = def%num + 1
                        else    !insert the defect in front of def
                            nullify(def2)
                            allocate(def2)
                            allocate(def2%defectType(SPECIES))
                            def2%defectType = defectType
                            def2%num = 1
                            def2%diff = diff
                            nullify(def2%reactionList)
                            def2%next => def
                            if(.NOT. associated(def_prev)) then !def is the first node
                                myDRL(level, cellCH)%defect_all => def2
                            else
                                def_prev%next => def2
                            end if
                        end if
                    else    !add the defect to the end of the list
                        nullify(def2)
                        allocate(def2)
                        allocate(def2%defectType(SPECIES))
                        def2%defectType = defectType
                        def2%num = 1
                        def2%diff = diff
                        nullify(def2%reactionList)
                        nullify(def2%next)
                        if(.NOT. associated(def_prev)) then !def2 will be the first node
                            myDRL(level, cellCH)%defect_all => def2
                        else
                            def_prev%next => def2
                        end if
                    end if

                    !<Fourth: update %defect_mobile
                    if(diff > 0d0) then
                        nullify(def_m)
                        def_m => myDRL(level, cellCH)%defect_mobile
                        call locate_defect(defectType, def_m, def_mPrev, numSame)
                        if(associated(def_m)) then
                            if(numSame == SPECIES) then
                                def_m%num = def_m%num + 1
                            else    !insert the defect in front of def_m
                                nullify(def2)
                                allocate(def2)
                                allocate(def2%defectType(SPECIES))
                                def2%defectType = defectType
                                def2%num = 1
                                def2%diff = diff
                                nullify(def2%reactionList)
                                def2%next => def_m
                                if(.NOT. associated(def_mPrev)) then !def_m is the first node
                                    myDRL(level, cellCH)%defect_mobile => def2
                                else
                                    def_mPrev%next => def2
                                end if
                            end if
                        else
                            nullify(def2)
                            allocate(def2)
                            allocate(def2%defectType(SPECIES))
                            def2%defectType = defectType
                            def2%num = 1
                            def2%diff = diff
                            nullify(def2%reactionList)
                            nullify(def2%next)
                            if(.NOT. associated(def_mPrev)) then !def2 will be the first node
                                myDRL(level, cellCH)%defect_mobile => def2
                            else
                                def_mPrev%next => def2
                            end if
                        end if
                    end if
                end if
            end do
        end if  !<reacCH%numReactants/=-10
    end if  !<we have choosen a reaction
end subroutine

!***************************************************************************************************
!>Subroutine updateReactionList(numUpdates,defUpdateStep): defUpdateStep is an array
!Function: updates reaction rates and/or creates/deletes reactions according to defect that have changed.
!This subroutine does the following:
!1) update 1st reactions
!2) update diffusion reactions
!3) update 2nd reactions
!EDIT: a) if a defect type is delete, then delete all reactions associated with it from the DRL.
!      b) if a new defect type is added, then add reactions associated with it to the DRL.
!      c) if the defect is mobile ,then update others defects' clustering reactions that contains the defect.
!***************************************************************************************************
subroutine SL_updateReactionList(numUpdates,defUpdateStep)
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    use mod_update_reactions
    implicit none

    integer, intent(in) :: numUpdates
    integer, intent(inout), dimension(SPECIES+6,numUpdates) :: defUpdateStep

    type(defect), pointer :: def, def_prev, def_m, def2
    integer :: defectType(SPECIES)
    integer :: cell, proc, dirTemp, neighbor, cascadeID
    integer :: level, level2, level_m, dir, negDir, i, j, numSame
    double precision :: diff
    double precision, external :: compute_diffusivity
    integer, external :: find_level

    do i=1, numUpdates
        !if(defUpdateStep(SPECIES+1,i)/=1 .AND. defUpdateStep(SPECIES+1,i)/=-1) then
        !    write(*,*) 'error defUpdateStep num not equal to +/- 1', myProc%taskid, defUpdateStep(SPECIES+1,i)
        !end if
        if(defUpdateStep(SPECIES+1,i)==1 .OR. defUpdateStep(SPECIES+1,i)==-1) then

        defectType(1:SPECIES) = defUpdateStep(1:SPECIES,i)
        cell = defUpdateStep(SPECIES+2,i)
        proc = defUpdateStep(SPECIES+3,i)
        dirTemp = defUpdateStep(SPECIES+4,i)
        neighbor = defUpdateStep(SPECIES+5,i)
        cascadeID = defUpdateStep(SPECIES+6,i)

        level = find_level(defectType)
        !*****************************************************
        !<if the defect is within the local mesh, update all relevant reactions
        !*****************************************************
        if(proc == myProc%taskid) then
            !*****************************************************
            !<defUpdateStep(SPECIES+6,i)==0 means a defect has changed in the coarse mesh
            !*****************************************************
            if(cascadeID == 0) then
                !<locate this defect
                nullify(def)
                def => myDRL(level,cell)%defect_all
                nullify(def_prev)
                do while(associated(def))
                    numSame = 0
                    do j=1, SPECIES
                        if(def%defectType(j) == defectType(j)) then
                            numSame = numSame + 1
                        end if
                    end do
                    if(numSame == SPECIES) then
                        exit
                    else
                        def_prev => def
                        def => def%next
                    end if
                end do

                !*****************************************
                !<Update reactions
                !*****************************************
                !<delete def%reactionList
                if(associated(def) .AND. def%num==0) then
                    call delete_reactions(def, def_prev, cell, level)
                !<add reactionList for def
                else if(associated(def) .AND. def%num/=0) then  !add reactionList
                    !<Update dissociation reactions and sinkRemovel reactions
                    call update_1st_Reaction(def, cell)

                    !<Update diffusion reactions, contain this cell to neighbor cell and neighbor cell to this cell
                    if(def%diff > 0d0) then
                        !call update_diff_reaction_cellTOsix(def, cell)
                        !<Update diff reactions in six neighbors
                        do dir=1,6
                            call update_diff_reaction_cellTOcell(def,cell,myMesh(cell)%neighbor(dir),&
                                    dir,myProc%taskid,myMesh(cell)%neighborProc(dir))

                            if(myMesh(cell)%neighborProc(dir) == myProc%taskid) then
                                if(mod(dir,2) == 0) then
                                    negDir = dir - 1
                                else
                                    negDir = dir + 1
                                end if

                                !find this type of defect
                                nullify(def2)
                                def2 => myDRL(level,myMesh(cell)%neighbor(dir))%defect_all
                                defLoop: do while(associated(def2))
                                    numSame = 0
                                    do j=1, SPECIES
                                        if(def2%defectType(j) == defectType(j)) then
                                            numSame = numSame +1
                                        end if
                                    end do
                                    if(numSame == SPECIES) then
                                        exit defLoop
                                    end if
                                    def2 => def2%next
                                end do defLoop

                                if(associated(def2)) then   !find it
                                    call update_diff_reaction_cellTOcell(def2,myMesh(cell)%neighbor(dir),cell,&
                                            negDir,myProc%taskid,myProc%taskid)
                                end if
                            end if
                        end do
                    end if

                    !<Update clustering reactions
                    if(def%diff <= 0d0) then
                        do level_m=1, LEVELS
                            nullify(def_m)
                            def_m => myDRL(level_m,cell)%defect_mobile
                            do while(associated(def_m))
                                !call update_2nd_Reaction(def, def_m, cell)
                                call update_2nd_reactionCompare(def, def_m, cell)
                                def_m => def_m%next
                            end do
                        end do
                    else
                        do level2=1, LEVELS
                            def2 => myDRL(level2,cell)%defect_all
                            do while(associated(def2))
                                !call update_2nd_Reaction(def2, def, cell)  !这一步，def2可能是新增的缺陷，可能此时还没轮到更新def2%reactionList
                                call update_2nd_reactionCompare(def2, def, cell)
                                def2 => def2%next
                            end do
                        end do
                    end if
                else    !i=2 & not associated(def) or i=1 & not associated(def)
                    if(i == 1) then
                        write(*,*) 'Error: updateReactionList---step', step, 'cell', cell,'proc',proc
                    end if
                end if
            else
                !cascadeID /= 0, update reactions in fine mesh
                !adding...
            end if
        else    !<defectUpdateCurrent%proc/=myProc%taskid
            !*****************************************************
            !<The defect is in the ghost mesh, this means that a defect has changed in the ghost to this mesh but not IN the mesh.
            !<Therefore we only need to update one diffusion reaction in one direction, the direction of the mesh with the changed defect.
            !"neighbor" is local cell, "cell" is ghost cell
            !*****************************************************
            !这一步似乎不需要，该缺陷在ghost中只有一种情况，即发生的是扩散反应，则更新反应物相关的反应时，已更新了neighbor到cell的扩散反应
            !diff = compute_diffusivity(defectType)
            !if(diff > 0d0) then
            !    nullify(def)
            !    def => myDRL(level,neighbor)%defect_all !"%neighbor" is local cell, "%cell" is ghost cell``````
            !    inter3: do while(associated(def))
            !        numSame = 0
            !        do j=1, SPECIES
            !            if(def%defectType(j) == defectType(j)) then
            !                numSame = numSame + 1
            !            end if
            !        end do
            !        if(numSame == SPECIES) then
            !            exit inter3
            !        end if
            !        def => def%next
            !    end do inter3

            !    if(associated(def)) then
            !        call update_diff_reaction_cellTOcell(def,def%diff,neighbor,cell,&
            !                dirTemp,myProc%taskid,proc)
            !    end if
            !end if
        end if

        end if  !if(defUpdateStep(SPECIES+1,i)==1 .OR. defUpdateStep(SPECIES+1,i)==-1) then
    end do
end subroutine

!***************************************************************************************************
!>Subroutine SL_synDefectList(sector, defCommList, defUpdate)
!Function: communication for common defects
!***************************************************************************************************
subroutine SL_synDefectList(sector, defCommList, defUpdate, casCellList, arrayTemp)
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none
    include 'mpif.h'

    integer, intent(in) :: sector
    type(defectCommHead), intent(inout), dimension(3) :: defCommList
    type(defectUpdate), pointer, intent(inout) :: defUpdate
    logical, intent(in), dimension(numMeshes) :: casCellList
    integer, intent(inout), dimension(2,3) :: arrayTemp

    type(defectCommHead), dimension(3) :: secondCommList
    type(defectComm), pointer :: defComm, defComm_prev, secondComm,  secondComm_prev
    type(defect), pointer :: def, def_prev, def2, def_m, def_mPrev

    integer :: dim, dir, recvDir, sendDir, send2Dir, count, numSame
    integer :: i, j, k, level, index, defectType(SPECIES), cell
    double precision :: diff
    integer :: numSend(3), numRecv(3), numSendFinal(3), numRecvFinal(3), recvCount
    !type(buffer), allocatable :: send(:), recv(:), sendFinal(:), recvFinal(:)
    integer, allocatable :: sendBuff(:,:), recvBuff(:,:), sendFinalBuff(:,:), recvFinalBuff(:,:)
    integer :: sendProc, recvProc
    integer :: istat        ! status, 0: success!
    integer :: status(MPI_STATUS_SIZE)
    integer :: sendStatus(MPI_STATUS_SIZE), recvStatus(MPI_STATUS_SIZE)
    integer :: sendRequest(3), recvRequest(3)
    integer :: sendStatusFinal(MPI_STATUS_SIZE), recvStatusFinal(MPI_STATUS_SIZE)
    integer :: sendRequestFinal(3), recvRequestFinal(3)
    !用于统计cas缺陷的数量
    integer :: numCellsSend(3), totDefects(3), numDefTemp, numSend_two(2,3), numRecv_two(2,3)
    !time counter
    double precision :: commTime1, commTime2
    !<Functions
    integer, external :: find_level, find_index
    double precision, external :: compute_diffusivity

    interface
        subroutine locate_defect(defectType, def, def_prev, numSame)
            use mod_constants
            use mod_structuretype
            implicit none
            integer, intent(in), dimension(SPECIES) :: defectType
            type(defect), pointer, intent(inout) :: def, def_prev
            integer, intent(inout) :: numSame
        end subroutine
    end interface

    !<Initialization
    numSend=0
    numRecv=0
    numSendFinal=0
    numRecvFinal=0

    numCellsSend = 0
    totDefects = 0
    arrayTemp = 0

    numSend_two = 0
    numRecv_two = 0

    !<initialize secondCommList
    do dim=1,3
        nullify(secondCommList(dim)%commHead)
    end do

    do dim=1, 3
        sendDir = mySector(sector)%commDir(1,dim)
        recvDir = mySector(sector)%commDir(2,dim)
        if(sendDir == 0) then
            sendProc = MPI_PROC_NULL
        else
            sendProc = myProc%neighborProc(sendDir)
        end if
        if(recvDir == 0) then
            recvProc = MPI_PROC_NULL
        else
            recvProc = myProc%neighborProc(recvDir)
        end if

        !<Count the number of data to be sent
        nullify(defComm)
        defComm => defCommList(dim)%commHead
        do while(associated(defComm))
            if(defComm%num /= 0) then
                numSend(dim) = numSend(dim) + 1
            end if
            defComm => defComm%next
        end do

        !统计cas缺陷的数量
        if(irradiationType=='Cascade') then
            do cell=1, numMeshes
                if(casCellList(cell) .eqv. .TRUE.) then
                    do dir=dim*2-1, dim*2
                        if(myMesh(cell)%neighborProc(dir)/=myProc%taskid .AND. myMesh(cell)%neighborProc(dir)/=-1) then
                            numCellsSend(dim) = numCellsSend(dim) + 1
                            numDefTemp = 0
                            do level=1, LEVELS
                                def => myDRL(level, cell)%defect_mobile
                                do while(associated(def))
                                    numDefTemp = numDefTemp + 1
                                    def => def%next
                                end do
                            end do
                            totDefects(dim) = totDefects(dim) + numDefTemp
                        end if
                    end do
                end if
            end do
            if(totDefects(dim) > 0) then
                arrayTemp(1,dim) = numCellsSend(dim) + totDefects(dim)
            else
                arrayTemp(1,dim) = 0
            end if
        end if

        allocate(sendBuff(SPECIES+3,numSend(dim)+1))    !sendBuff(1:SPECIES+3,1) used to count cascade defects
        sendBuff(1,1) = arrayTemp(1,dim)
        sendBuff(2:SPECIES+3,1) = 0

        !numSend_two(1,dim) = numSend(dim)
        !numSend_two(2,dim) = arrayTemp(1,dim)

        !commTime1 = MPI_WTIME()
        !call MPI_SENDRECV(numSend_two(1,dim),2,MPI_INTEGER,myProc%neighborProc(sendDir),99,&
        !        numRecv_two(1,dim),2,MPI_INTEGER,myProc%neighborProc(recvDir),99,comm,status,ierr)
        !commTime2 = MPI_WTIME()
        !commTime11  = commTime11 + (commTime2-commTime1)
        !commCount_1 = commCount_1 + 1

        !numRecv(dim) = numRecv_two(1,dim)
        !arrayTemp(2,dim) = numRecv_two(2,dim)

        !<Create a send package in the dim dimension
        if(numSend(dim) > 0) then
            !allocate(sendBuff(SPECIES+3,numSend(dim)))
            count = 1
            nullify(defComm)
            defComm => defCommList(dim)%commHead
            do while(associated(defComm))
                if(defComm%num /= 0) then
                    count = count + 1
                    sendBuff(1:SPECIES, count) = defComm%defectType  ! defectType
                    sendBuff(SPECIES+1, count) = defComm%num        ! number of defects with this defectType
                    sendBuff(SPECIES+2, count) = defComm%cell       ! cell id of these defects
                    sendBuff(SPECIES+3, count) = defComm%neighbor   ! tag of cell, 1.neighbor=-100: cell belongs to ghost area, 2.neighbor=neighbor cell id of this cell in sendDir direction
                end if
                defComm => defComm%next
            end do
            if(count-1 /= numSend(dim)) then
                write(*,*) 'Error number of defects that are sent. proc',myProc%taskid,'dim',dim,'sendDir',sendDir
            end if
        end if

        commTime1 = MPI_WTIME()
        call MPI_ISEND(sendBuff,(SPECIES+3)*(numSend(dim)+1),MPI_INTEGER,sendProc,100+dim,comm,sendRequest(dim),ierr)
        commTime2 = MPI_WTIME()
        commTime12  = commTime12 + (commTime2-commTime1)

        recvCount = 0
        commTime1 = MPI_WTIME()
        call MPI_PROBE(recvProc, 100+dim, comm, status,ierr)
        call MPI_GET_COUNT(status, MPI_INTEGER, recvCount, ierr)
        numRecv(dim) = recvCount/(SPECIES+3)
        commTime2 = MPI_WTIME()
        commTime12  = commTime12 + (commTime2-commTime1)

        !send defect datas
        !if(numSend(dim) > 0) then
        !    commTime1 = MPI_WTIME()
        !    call MPI_SEND(sendBuff,(SPECIES+3)*numSend(dim),MPI_INTEGER,myProc%neighborProc(sendDir),100,comm,ierr)
        !    commTime2 = MPI_WTIME()
        !    commTime12  = commTime12 + (commTime2-commTime1)
        !    deallocate(sendBuff)
        !end if

        !recv defect datas
        if(numRecv(dim) > 0) then
            allocate(recvBuff(SPECIES+3,numRecv(dim)))
            commTime1 = MPI_WTIME()
            call MPI_RECV(recvBuff,(SPECIES+3)*numRecv(dim),MPI_INTEGER,recvProc,100+dim,comm,status,ierr)
            commTime2 = MPI_WTIME()
            commTime12  = commTime12 + (commTime2-commTime1)
            arrayTemp(2,dim) = recvBuff(1,1)
        end if

        do i=2, numRecv(dim)
            defectType(1:SPECIES) = recvBuff(1:SPECIES,i)
            level = find_level(defectType)
            diff = compute_diffusivity(defectType)

            !<update local defects
            if(recvBuff(SPECIES+3,i) == -100) then !update boundary,num>0
                if(recvBuff(SPECIES+1,i) <= 0) then  !check num
                    write(*,*) 'Error number of recv defects:', recvBuff(1:SPECIES+3,i)
                    call MPI_ABORT(comm,ierr)
                end if

                !create a new element in defUpdateList and assign all variables except for num
                allocate(defUpdate%next,STAT=istat)
                defUpdate => defUpdate%next
                allocate(defUpdate%defectType(SPECIES))
                defUpdate%defectType = defectType
                defUpdate%cell = recvBuff(SPECIES+2,i)
                defUpdate%proc = myProc%taskid
                defUpdate%dir = 0
                defUpdate%neighbor = 0
                defUpdate%cascadeID = 0	            !not inside a cascade
                nullify(defUpdate%next)

                !Add the received defects to the %defect_all
                nullify(def)
                def => myDRL(level, recvBuff(SPECIES+2,i))%defect_all
                call locate_defect(defectType, def, def_prev, numSame)
                if(associated(def)) then
                    if(numSame == SPECIES) then     !this type of defect is existed
                        def%num = def%num + recvBuff(SPECIES+1,i)
                    else    !insert the defect in front of def
                        nullify(def2)
                        allocate(def2,STAT=istat)
                        allocate(def2%defectType(SPECIES))
                        def2%defectType = defectType
                        def2%num = recvBuff(SPECIES+1,i)
                        def2%diff = diff
                        if(.NOT. associated(def_prev)) then !def is the first node
                            def2%next => myDRL(level, recvBuff(SPECIES+2,i))%defect_all
                            myDRL(level, recvBuff(SPECIES+2,i))%defect_all => def2
                        else    !def2 is inserted in the middle
                            def2%next => def
                            def_prev%next => def2
                        end if
                        nullify(def2%reactionList)
                    end if
                else    !def2 is added to the end of the list
                    nullify(def2)
                    allocate(def2, STAT=istat)
                    allocate(def2%defectType(SPECIES))
                    def2%defectType = defectType
                    def2%num = recvBuff(SPECIES+1,i)
                    def2%diff = diff
                    if(.NOT. associated(def_prev)) then !def2 will be the first node
                        myDRL(level, recvBuff(SPECIES+2,i))%defect_all => def2
                        nullify(def2%next)
                    else
                        def_prev%next => def2
                        nullify(def2%next)
                    end if
                    nullify(def2%reactionList)
                end if

                !Add the received defects to the %defect_mobile
                if(diff > 0d0) then!this defect is mobile
                    nullify(def_m)
                    def_m => myDRL(level, recvBuff(SPECIES+2,i))%defect_mobile
                    call locate_defect(defectType, def_m, def_mPrev, numSame)
                    if(associated(def_m)) then
                        if(numSame == SPECIES) then
                            def_m%num = def_m%num + recvBuff(SPECIES+1,i)
                        else
                            nullify(def2)
                            allocate(def2,STAT=istat)
                            allocate(def2%defectType(SPECIES))
                            def2%defectType = defectType
                            def2%num = recvBuff(SPECIES+1,i)
                            def2%diff = diff
                            if(.NOT. associated(def_mPrev)) then !def_m is the first node
                                def2%next => myDRL(level, recvBuff(SPECIES+2,i))%defect_mobile
                                myDRL(level, recvBuff(SPECIES+2,i))%defect_mobile => def2
                            else
                                def2%next => def_m
                                def_mPrev%next => def2
                            end if
                            nullify(def2%reactionList)
                        end if
                    else
                        nullify(def2)
                        allocate(def2,STAT=istat)
                        allocate(def2%defectType(SPECIES))
                        def2%defectType = defectType
                        def2%num = recvBuff(SPECIES+1,i)
                        def2%diff = diff
                        if(.NOT. associated(def_mPrev)) then !def2 will be the first node
                            myDRL(level, recvBuff(SPECIES+2,i))%defect_mobile => def2
                            nullify(def2%next)
                        else
                            def_mPrev%next => def2
                            nullify(def2%next)
                        end if
                        nullify(def2%reactionList)
                    end if
                end if

                !*****************************************************
                !if this cell is in the ghost of other neighboring processors
                !(not myProc%neighborProc(recvDir), a second communication is need
                !*****************************************************
                do k=1, 3
                    if(k /= dim) then
                        do send2Dir=2*k-1, 2*k
                            if(myMesh(recvBuff(SPECIES+2,i))%neighborProc(send2Dir)/=myProc%taskid .AND. &
                                        myMesh(recvBuff(SPECIES+2,i))%neighborProc(send2Dir)/=-1) then
                                nullify(secondComm)
                                nullify(secondComm_prev)
                                secondComm => secondCommList(k)%commHead
                                comm2: do while(associated(secondComm))
                                    if(secondComm%cell==recvBuff(SPECIES+2,i) .AND. secondComm%dir==0) then
                                        count = 0
                                        do j=1, SPECIES
                                            if(secondComm%defectType(j) == defectType(j)) then
                                                count = count + 1
                                            end if
                                        end do
                                        if(count == SPECIES) then !this type of defect is existed in secondCommList(k)
                                            secondComm%num = secondComm%num + recvBuff(SPECIES+1,i)
                                            exit comm2
                                        end if
                                    end if
                                    secondComm_prev => secondComm
                                    secondComm => secondComm%next
                                end do comm2
                                if(.NOT. associated(secondComm)) then   !add node in secondCommList(k)
                                    numSendFinal(k) = numSendFinal(k) + 1
                                    allocate(secondComm,STAT=istat)
                                    allocate(secondComm%defectType(SPECIES))
                                    secondComm%defectType = defectType
                                    secondComm%num = recvBuff(SPECIES+1,i)
                                    secondComm%cell = recvBuff(SPECIES+2,i)
                                    secondComm%neighbor = myMesh(recvBuff(SPECIES+2,i))%neighbor(send2Dir)
                                    secondComm%dir = 0
                                    nullify(secondComm%next)
                                    if(.NOT. associated(secondComm_prev)) then
                                        secondCommList(k)%commHead => secondComm
                                    else
                                        secondComm_prev%next => secondComm
                                    end if
                                end if
                            end if
                        end do
                    end if
                end do
            !<update ghost defects
            else    !num >0 or <0
                allocate(defUpdate%next)
                defUpdate=>defUpdate%next
                allocate(defUpdate%defectType(SPECIES))
                defUpdate%defectType = defectType
                defUpdate%cell = recvBuff(SPECIES+2,i)            !ghost cell
                defUpdate%proc = myProc%neighborProc(recvDir)
                defUpdate%dir = recvDir
                defUpdate%neighbor = recvBuff(SPECIES+3,i)        !local cell
                defUpdate%cascadeID = 0                     !not inside a cascade
                nullify(defUpdate%next)

                !Add the received defects to the myGhost%defectList()
                index = find_index(recvBuff(SPECIES+3,i),recvDir)
                nullify(def)
                def => myGhost(recvDir)%myDRL_ghost(level,index)%defect_mobile
                call locate_defect(defectType, def, def_prev, numSame)
                if(associated(def)) then
                    if(numSame == SPECIES) then
                        def%num = def%num + recvBuff(SPECIES+1,i)    !num > 0 or <0
                        if(def%num < 0) then
                            write(*,*) 'Error delete more defects in myDRL_ghost'
                            write(*,*) 'recvNum',recvBuff(SPECIES+1,i),'beforeNum', def%num-recvBuff(SPECIES+1,i), &
                                    'afterNum', def%num
                            call MPI_ABORT(comm,ierr)
                        else if(def%num == 0) then  !delete this node
                            if(.NOT. associated(def_prev)) then
                                myGhost(recvDir)%myDRL_ghost(level,index)%defect_mobile => def%next
                            else
                                def_prev%next => def%next
                            end if
                            deallocate(def%defectType)
                            deallocate(def,STAT=istat)
                        end if
                    else    !insert the defect in front of def
                        if(recvBuff(SPECIES+1,i) < 0) then
                            write(*,*) '1Error in recv negative defect numbers---step',step,'proc', myProc%taskid
                            write(*,*) 'dim', dim, 'recvDir', recvDir,'i',i, 'numRecv',numRecv(dim),'index', index
                            call MPI_ABORT(comm,ierr)
                        else    !num>0, insert in front of def
                            nullify(def2)
                            allocate(def2,STAT=istat)
                            allocate(def2%defectType(SPECIES))
                            def2%defectType = defectType
                            def2%num = recvBuff(SPECIES+1,i)
                            def2%diff = diff
                            if(.NOT. associated(def_prev)) then !def_m is the first node
                                def2%next => myGhost(recvDir)%myDRL_ghost(level,index)%defect_mobile
                                myGhost(recvDir)%myDRL_ghost(level,index)%defect_mobile => def2
                            else
                                def2%next => def
                                def_prev%next => def2
                            end if
                            nullify(def2%reactionList)
                        end if
                    end if
                else    !insert the defect after def_prev
                    if(recvBuff(SPECIES+1,i) < 0) then
                        write(*,*) '2Error in recv(dim)%datas negative defect numbers---step',step
                        write(*,*) 'dim', dim, 'recvDir', recvDir,'i',i, 'numRecv(dim)',numRecv(dim)
                        write(*,*) 'def',defectType,'num',recvBuff(SPECIES+1,i),'cell',recvBuff(SPECIES+2,i)
                        write(*,*) 'index', index, 'proc', myProc%taskid
                        call MPI_ABORT(comm,ierr)
                    else
                        nullify(def2)
                        allocate(def2,STAT=istat)
                        allocate(def2%defectType(SPECIES))
                        def2%defectType = defectType
                        def2%num = recvBuff(SPECIES+1,i)
                        def2%diff = diff
                        if(.NOT. associated(def_prev)) then !def_m is the first node
                            myGhost(recvDir)%myDRL_ghost(level,index)%defect_mobile => def2
                            nullify(def2%next)
                        else
                            def_prev%next => def2
                            nullify(def2%next)
                        end if
                        nullify(def2%reactionList)
                    end if
                end if
            end if  !update defects in boundary or ghost area
        end do  !i=1, numRecv(dim)
        if(numRecv(dim) > 0) then
            deallocate(recvBuff)
        end if

        commTime1 = MPI_WTIME()
        call MPI_WAIT(sendRequest(dim), sendStatus, ierr)
        commTime2 = MPI_WTIME()
        commTime12  = commTime12 + (commTime2-commTime1)
        deallocate(sendBuff)

    end do  !dim=1, 3

    !commTime1 = MPI_WTIME()
    !call MPI_BARRIER(comm,ierr)
    !commTime2 = MPI_WTIME()
    !commTime11  = commTime11 + (commTime2-commTime1)

    !<********************************************
    !<Second communication
    !<********************************************
    do dim=1, 3
        !!!!需要二次通信的cell不在当前sector中，但二次发送方向与当前sector的发送方向一致
        sendDir = mySector(sector)%commDir(1,dim)
        recvDir = mySector(sector)%commDir(2,dim)
        if(sendDir == 0) then
            sendProc = MPI_PROC_NULL
        else
            sendProc = myProc%neighborProc(sendDir)
        end if
        if(recvDir == 0) then
            recvProc = MPI_PROC_NULL
        else
            recvProc = myProc%neighborProc(recvDir)
        end if

        !Package datas
        !allocate(sendFinalBuff(SPECIES+3,numSendFinal(dim)))
        if(numSendFinal(dim) > 0) then
            allocate(sendFinalBuff(SPECIES+3,numSendFinal(dim)))
            count = 0
            nullify(secondComm)
            secondComm => secondCommList(dim)%commHead
            do while(associated(secondComm))
                count = count + 1
                sendFinalBuff(1:SPECIES,count) = secondComm%defectType(1:SPECIES)
                sendFinalBuff(SPECIES+1,count) = secondComm%num
                sendFinalBuff(SPECIES+2,count) = secondComm%cell        !local cell
                sendFinalBuff(SPECIES+3,count) = secondComm%neighbor    !ghost cell
                secondComm=>secondComm%next
            end do
            if(count /= numSendFinal(dim)) then     !check
                write(*,*) 'error second send'
                call MPI_ABORT(comm,ierr)
            end if
        end if

        numRecvFinal(dim) = 0
        commTime1 = MPI_WTIME()
        call MPI_SENDRECV(numSendFinal(dim),1,MPI_INTEGER,sendProc,99,&
                numRecvFinal(dim),1,MPI_INTEGER,recvProc,99,comm,status,ierr)
        commTime2 = MPI_WTIME()
        commTime21 = commTime21 + (commTime2-commTime1)
        commCount_2 = commCount_2 + 1

        !send defect datas
        if(numSendFinal(dim) > 0) then
            commTime1 = MPI_WTIME()
            call MPI_SEND(sendFinalBuff,(SPECIES+3)*numSendFinal(dim),MPI_INTEGER,sendProc,100,comm,ierr)
            commTime2 = MPI_WTIME()
            commTime22 = commTime22 + (commTime2-commTime1)
            deallocate(sendFinalBuff)
        end if

        !commTime1 = MPI_WTIME()
        !call MPI_ISEND(sendFinalBuff,(SPECIES+3)*numSendFinal(dim),MPI_INTEGER,myProc%neighborProc(sendDir),&
        !        200+dim,comm,sendRequestFinal(dim),ierr)
        !commTime2 = MPI_WTIME()
        !commTime22 = commTime22 + (commTime2-commTime1)

        !recvCount = 0
        !commTime1 = MPI_WTIME()
        !call MPI_PROBE(myProc%neighborProc(recvDir), 200+dim, comm, status,ierr)
        !call MPI_GET_COUNT(status, MPI_INTEGER, recvCount, ierr)
        !numRecvFinal(dim) = recvCount/(SPECIES+3)
        !commTime2 = MPI_WTIME()
        !commTime22 = commTime22 + (commTime2-commTime1)

        !recv defect datas
        if(numRecvFinal(dim) > 0) then
            allocate(recvFinalBuff(SPECIES+3,numRecvFinal(dim)))
            commTime1 = MPI_WTIME()
            call MPI_RECV(recvFinalBuff,(SPECIES+3)*numRecvFinal(dim),MPI_INTEGER,recvProc,&
                    100, comm,status,ierr)
            commTime2 = MPI_WTIME()
            commTime22 = commTime22 + (commTime2-commTime1)
        end if

        !<update ghost area, all num>0
        do i=1, numRecvFinal(dim)   !num > 0
            if(recvFinalBuff(SPECIES+1,i) < 0) then  !check
                write(*,*) 'Second Recv Error: error number of defects that second communication'
                call MPI_ABORT(comm,ierr)
            end if

            defectType(1:SPECIES) = recvFinalBuff(1:SPECIES,i)
            diff = compute_diffusivity(defectType)
            level = find_level(defectType)

            !create a new element in defUpdateList and assign all variables
            allocate(defUpdate%next,STAT=istat)
            defUpdate => defUpdate%next
            allocate(defUpdate%defectType(SPECIES))
            defUpdate%defectType = defectType
            defUpdate%cell = recvFinalBuff(SPECIES+2,i)
            defUpdate%proc = myProc%neighborProc(recvDir)
            defUpdate%dir = recvDir
            defUpdate%neighbor = recvFinalBuff(SPECIES+3,i)   !local cell
            defUpdate%cascadeID = 0	!not inside a cascade
            nullify(defUpdate%next)

            !Add the received defects to the myGhost%defectList()
            index = find_index(recvFinalBuff(SPECIES+3,i),recvDir)
            nullify(def)
            def => myGhost(recvDir)%myDRL_ghost(level,index)%defect_mobile
            call locate_defect(defectType, def, def_prev, numSame)
            if(associated(def)) then
                if(numSame == SPECIES) then
                    def%num = def%num + recvFinalBuff(SPECIES+1,i)
                else    !insert the defect in front of def
                    nullify(def2)
                    allocate(def2,STAT=istat)
                    allocate(def2%defectType(SPECIES))
                    def2%defectType = defectType
                    def2%num = recvFinalBuff(SPECIES+1,i)
                    def2%diff = diff
                    if(.NOT. associated(def_Prev)) then !def_m is the first node
                        def2%next => myGhost(recvDir)%myDRL_ghost(level,index)%defect_mobile
                        myGhost(recvDir)%myDRL_ghost(level,index)%defect_mobile => def2
                    else
                        def2%next => def
                        def_Prev%next => def2
                    end if
                    nullify(def2%reactionList)
                end if
            else    !insert the defect at the end of the list
                nullify(def2)
                allocate(def2,STAT=istat)
                allocate(def2%defectType(SPECIES))
                def2%defectType = defectType
                def2%num = recvFinalBuff(SPECIES+1,i)
                def2%diff = diff
                if(.NOT. associated(def_Prev)) then !def_m is the first node
                    myGhost(recvDir)%myDRL_ghost(level,index)%defect_mobile => def2
                    nullify(def2%next)
                else
                    def_Prev%next => def2
                    nullify(def2%next)
                end if
                nullify(def2%reactionList)
            end if
        end do  !i=1, numRecvFinal(dim)
        if(numRecvFinal(dim) > 0) then
            deallocate(recvFinalBuff)
        end if

        !commTime1 = MPI_WTIME()
        !call MPI_WAIT(sendRequestFinal(dim), sendStatusFinal, ierr)
        !commTime2 = MPI_WTIME()
        !commTime22  = commTime22 + (commTime2-commTime1)
        !deallocate(sendFinalBuff)

    end do  !dim=1,3
    
    !delete defCommList%next
    do dim=1,3
        nullify(defComm)
        defComm => defCommList(dim)%commHead
        do while(associated(defComm))
            defCommList(dim)%commHead => defComm%next
            deallocate(defComm%defectType)
            deallocate(defComm,STAT=istat)
            defComm => defCommList(dim)%commHead
        end do
    end do

    !delete secondCommList
    do dim=1,3
        nullify(secondComm)
        secondComm => secondCommList(dim)%commHead
        do while(associated(secondComm))
            secondCommList(dim)%commHead => secondComm%next
            deallocate(secondComm%defectType)
            deallocate(secondComm,STAT=istat)
            secondComm => secondCommList(dim)%commHead
        end do
    end do

end subroutine

!***************************************************************************************************
!>Subroutine SL_synReactionList(defUpdateList)
!Updates reaction rates and/or creates/deletes reactions according to defect numbers that have changed.
!This subroutine does the following:
!1) update 1st reactions
!2) update diffusion reactions
!3) update 2nd reactions
!EDIT:
!      b) if a new defect type is added, then add reactions associated with it to the DRL.
!      c) if the defect is mobile ,then update others defects' clustering reactions that contains the defect.
!***************************************************************************************************
subroutine SL_synReactionList(defUpdateList)
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    use mod_update_reactions
    implicit none

    type(defectUpdate), pointer, intent(in) :: defUpdateList
    type(defectUpdate), pointer :: defUpdate
    type(defect), pointer :: def, def_prev, def_m, def2
    integer :: defectType(SPECIES), cell, proc, dirTemp, neighbor, cascadeID
    integer :: level, level2, level_m, dir, negDir, numSame, j
    double precision :: diff
    integer, external :: find_level
    double precision, external :: compute_diffusivity

    nullify(defUpdate)
    defUpdate => defUpdateList%next
    do while(associated(defUpdate))
        defectType = defUpdate%defectType
        cell = defUpdate%cell
        proc = defUpdate%proc
        dirTemp = defUpdate%dir
        neighbor = defUpdate%neighbor
        cascadeID = defUpdate%cascadeID

        level = find_level(defectType)
        !*****************************************************
        !<if the defect is within the local mesh
        !*****************************************************
        if(proc == myProc%taskid) then    !neighbor proc diffuse to this proc
            !*****************************************************
            !<defUpdate%cascadeNumber==0 means a defect has changed in the coarse mesh
            !*****************************************************
            if(cascadeID == 0) then
                !<locate this defect
                nullify(def)
                def => myDRL(level,cell)%defect_all
                nullify(def_prev)
                do while(associated(def))
                    numSame = 0
                    do j=1, SPECIES
                        if(def%defectType(j) == defectType(j)) then
                            numSame = numSame + 1
                        end if
                    end do
                    if(numSame == SPECIES) then
                        exit
                    end if
                    def_prev => def
                    def => def%next
                end do

                !*****************************************
                !<Update reactions
                !*****************************************
                if(associated(def) .AND. def%num==0) then   !delete reactionList
                    call delete_reactions(def, def_prev, cell, level)
                else if(associated(def) .AND. def%num/=0) then   !add reactionList
                    !<Update dissociation reactions and sinkRemovel reactions
                    call update_1st_Reaction(def, cell)
                    !<Update diffusion reactions, contain this cell to neighbor cell and neighbor cell to this cell
                    if(def%diff > 0d0) then
                        !call update_diff_reaction_cellTOsix(def, cell)
                        !<Update diff reactions in six neighbors
                        do dir=1,6
                            call update_diff_reaction_cellTOcell(def,cell,myMesh(cell)%neighbor(dir),&
                                    dir,myProc%taskid,myMesh(cell)%neighborProc(dir))

                            if(myMesh(cell)%neighborProc(dir) == myProc%taskid) then
                                if(mod(dir,2) == 0) then
                                    negDir = dir - 1
                                else
                                    negDir = dir + 1
                                end if

                                !find this type of defect
                                nullify(def2)
                                def2 => myDRL(level,myMesh(cell)%neighbor(dir))%defect_all
                                defLoop: do while(associated(def2))
                                    numSame = 0
                                    do j=1, SPECIES
                                        if(def2%defectType(j) == defectType(j)) then
                                            numSame = numSame +1
                                        end if
                                    end do
                                    if(numSame == SPECIES) then
                                        exit defLoop
                                    end if
                                    def2 => def2%next
                                end do defLoop

                                if(associated(def2)) then   !find it
                                    call update_diff_reaction_cellTOcell(def2,myMesh(cell)%neighbor(dir),cell,&
                                            negDir,myProc%taskid,myProc%taskid)
                                end if
                            end if
                        end do
                    end if

                    !<Update clustering reactions
                    if(def%diff <= 0d0) then
                        do level_m=1, LEVELS
                            nullify(def_m)
                            def_m => myDRL(level_m,cell)%defect_mobile    !2021.10.27，将level改为了level_m
                            do while(associated(def_m))
                                !call update_2nd_Reaction(def, def_m, cell)
                                call update_2nd_reactionCompare(def, def_m, cell)
                                def_m => def_m%next
                            end do
                        end do
                    else
                        do level2=1, LEVELS
                            def2 => myDRL(level2,cell)%defect_all
                            do while(associated(def2))
                                !call update_2nd_Reaction(def2, def, cell)
                                call update_2nd_reactionCompare(def2, def, cell)
                                def2 => def2%next
                            end do
                        end do
                    end if
                else
                    write(*,*) 'Error update defects in synReactionList---cell', cell, 'step',step
                end if
            else
                !cascadeID/=0
                !adding...
            end if
        else    !<defUpdate%proc/=myProc%taskid
            !*****************************************************
            !<The defect is in the boundary mesh, this means that a defect has changed in the boundary to this mesh but not IN the mesh.
            !<Therefore we only need to update one diffusion reaction in one direction, the direction of the mesh with the changed defect.
            !No defects have changed in this mesh.
            !*****************************************************
            !diff = compute_diffusivity(defectType)
            !if(diff > 0d0) then
                !<locate this defect
                nullify(def)
                def => myDRL(level,neighbor)%defect_all !"%neighbor" is local cell, "%cell" is ghost cell
                inter2: do while(associated(def))
                    numSame = 0
                    do j=1, SPECIES
                        if(def%defectType(j) == defectType(j)) then
                            numSame = numSame + 1
                        end if
                    end do
                    if(numSame == SPECIES) then
                        exit inter2
                    end if
                    def => def%next
                end do inter2

                if(associated(def)) then
                    call update_diff_reaction_cellTOcell(def,neighbor,cell,dirTemp,myProc%taskid,proc)
                end if
            !end if
        end if

        defUpdate => defUpdate%next
    end do

    !*****************************************************
    !<Deallocate defUpdateList (memory release)
    !*****************************************************
    nullify(defUpdate)
    defUpdate => defUpdateList%next
    do while(associated(defUpdate))
        defUpdateList%next => defUpdate%next
        deallocate(defUpdate%defectType)
        deallocate(defUpdate)
        nullify(defUpdate)
        defUpdate => defUpdateList%next
    end do

end subroutine

!***************************************************************************************************
!>Subroutine SL_synCascade(sector, casCellList)
!Function: update defects in ghost cells that were implanted cascades, and update associated diffusion reactions
!***************************************************************************************************
subroutine SL_synCascade(sector, casCellList, arrayTemp)
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    use mod_update_reactions
    implicit none
    include 'mpif.h'

    integer, intent(in) :: sector
    logical, intent(in), dimension(numMeshes) :: casCellList
    integer, intent(in), dimension(2,3) :: arrayTemp

    type(defect), pointer :: def, def2
    integer :: count, i, k, k1, ghostCell, index, defectType(SPECIES)
    integer :: numDefectsLevel(LEVELS)
    integer :: cell, dim, level, dir, sendDir, recvDir, cell_index, numDefects
    integer :: numCellsSend(3), numCellsRecv(3)
    integer :: casNumSend(3), casNumRecv(3)
    integer, allocatable :: casSendBuff(:,:), casRecvBuff(:,:)
    integer :: status(MPI_STATUS_SIZE)
    !integer :: casSendStatus(MPI_STATUS_SIZE), casRecvStatus(MPI_STATUS_SIZE)
    !integer :: casSendRequest(3), casRecvRequest(3)
    double precision :: commTime1, commTime2
    !<Functions
    integer, external :: find_level, find_index
    double precision, external :: compute_diffusivity

    casNumSend(1:3) = arrayTemp(1, 1:3)
    casNumRecv(1:3) = arrayTemp(2, 1:3)

    do dim=1,3
        sendDir = mySector(sector)%commDir(1,dim)   !send Dir
        recvDir = mySector(sector)%commDir(2,dim)

        if(casNumSend(dim) > 0) then  !need send
            !<Create a send package in the dim dimension for cascade
            allocate(casSendBuff(SPECIES+1,casNumSend(dim)))
            casSendBuff = 0 !initialize send buff

            numCellsSend(dim) = 0
            cell_index=0
            do cell=1, numMeshes
                if(casCellList(cell) .eqv. .TRUE.) then !this cell must be in bundary area
                    do dir=dim*2-1, dim*2
                        if(myMesh(cell)%neighborProc(dir)/=myProc%taskid .AND. myMesh(cell)%neighborProc(dir)/=-1) then
                            numCellsSend(dim) = numCellsSend(dim) + 1
                            cell_index = cell_index + 1
                            casSendBuff(1,cell_index) = cell    !local cell id
                            casSendBuff(2,cell_index) = myMesh(cell)%neighbor(sendDir)  !neighbor cell id
                            numDefects = 0
                            do level=1, LEVELS
                                def => myDRL(level, cell)%defect_mobile
                                do while(associated(def))
                                    numDefects = numDefects + 1
                                    casSendBuff(1:SPECIES,cell_index+numDefects) = def%defectType(1:SPECIES)
                                    casSendBuff(SPECIES+1,cell_index+numDefects) = def%num
                                    def => def%next
                                end do
                            end do
                            casSendBuff(3,cell_index) = numDefects  !number of defects in tis cell
                            casSendBuff(4,cell_index) = 0
                            cell_index = cell_index + numDefects
                            exit
                        end if
                    end do
                end if
            end do
            casSendBuff(4,1) = numCellsSend(dim)    !number of cells in this dim
        end if

        !<Communication: SEND/Recv datas
        if(casNumSend(dim) > 0) then    !send datas
            commTime1 = MPI_WTIME()
            call MPI_SEND(casSendBuff,(SPECIES+1)*casNumSend(dim),MPI_INTEGER,myProc%neighborProc(sendDir),102,comm,ierr)
            commTime2 = MPI_WTIME()
            casCommTime = casCommTime + (commTime2-commTime1)
            deallocate(casSendBuff)
        end if

        if(casNumRecv(dim) > 0) then    !recv datas
            allocate(casRecvBuff(SPECIES+1,casNumRecv(dim)))
            commTime1 = MPI_WTIME()
            call MPI_RECV(casRecvBuff,(SPECIES+1)*casNumRecv(dim),MPI_INTEGER,myProc%neighborProc(recvDir),102,comm,status,ierr)
            commTime2 = MPI_WTIME()
            casCommTime = casCommTime + (commTime2-commTime1)
        end if

        if(casNumRecv(dim) > 0) then
            numCellsRecv(dim) = casRecvBuff(4,1)
        else
            numCellsRecv(dim) = 0
        end if



        !<update defects in ghost cells which were implanted with cascades
        cell_index = 0
        numDefects = 0  !number of defects in each cell
        do i=1, numCellsRecv(dim)    !traversing the cell with cascade
            if(i == 1) then
                cell_index = 1
                ghostCell = casRecvBuff(1,cell_index)
                cell = casRecvBuff(2,cell_index)   !local cell
                numDefects = casRecvBuff(3,cell_index)
            else
                cell_index = cell_index + numDefects + 1
                ghostCell = casRecvBuff(1,cell_index)
                cell = casRecvBuff(2,cell_index)   !local cell
                numDefects = casRecvBuff(3,cell_index)
            end if

            index = find_index(cell, recvDir)
            if(myGhost(recvDir)%cell(index) /= ghostCell) then  !<check
                write(*,*) 'Error ghost cell---ghostCell',ghostCell,'cell',myGhost(recvDir)%cell(index),'recvDir',&
                        recvDir, 'index',index, 'cell', cell, 'ghost',myMesh(cell)%neighbor(recvDir),'sector', sector
                call MPI_ABORT(comm,ierr)
            end if

            !<First: delete exiting defects
            do level=1, LEVELS
                def => myGhost(recvDir)%myDRL_ghost(level,index)%defect_mobile
                do while(associated(def))
                    myGhost(recvDir)%myDRL_ghost(level,index)%defect_mobile => def%next
                    deallocate(def%defectType)
                    deallocate(def)
                    def => myGhost(recvDir)%myDRL_ghost(level,index)%defect_mobile
                end do
            end do

            !<Second: find number of defects in each level
            numDefectsLevel = 0
            defectType = 0
            do k=1, numDefects
                defectType(1:SPECIES) = casRecvBuff(1:SPECIES,cell_index+k)
                level = find_level(defectType)
                numDefectsLevel(level) = numDefectsLevel(level) + 1
            end do

            !<Third: add defects to myGhost(recvDir)%myDRL_ghost(level,index)%defect_mobile
            k = 0
            do level=1, LEVELS
                if(level==1) then
                    k = 0
                else
                    k = k + numDefectsLevel(level-1)
                end if

                nullify(def)
                do k1=1, numDefectsLevel(level)
                    nullify(def2)
                    allocate(def2)
                    allocate(def2%defectType(SPECIES))
                    def2%defectType(1:SPECIES) = casRecvBuff(1:SPECIES,cell_index+k+k1)
                    def2%num = casRecvBuff(SPECIES+1,cell_index+k+k1)
                    def2%diff = compute_diffusivity(def2%defectType)
                    nullify(def2%next)
                    nullify(def2%reactionList)
                    if(.NOT. associated(myGhost(recvDir)%myDRL_ghost(level,index)%defect_mobile)) then
                        myGhost(recvDir)%myDRL_ghost(level,index)%defect_mobile => def2
                        def => myGhost(recvDir)%myDRL_ghost(level,index)%defect_mobile
                    else
                        def%next => def2
                        def => def%next
                    end if
                end do
            end do

            !<Fourth: update diffusion reaction
            do level=1, LEVELS
                def => myDRL(level, cell)%defect_all
                do while(associated(def))
                    if(def%diff > 0d0) then
                        call update_diff_reaction_cellTOcell(def,cell,ghostCell,recvDir,&
                                myProc%taskid,myProc%neighborProc(recvDir))
                    end if
                    def => def%next
                end do
            end do
        end do

        if(allocated(casRecvBuff)) then
            deallocate(casRecvBuff)
        end if
    end do

end subroutine