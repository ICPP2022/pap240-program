!*******************************************************************************************************
!> Subroutine outputTotalDefects: output the total defects and post processing for the entire volume
! NOTE: in this subroutine, each process sends its defects to process 0, which outputs and sum them.
! Outputs a list of all defects in system (defect type, num).
! Output the concentration of solute, vacancies, and SIAs of size large.
!*******************************************************************************************************
subroutine outputTotalDefects()
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none
    include 'mpif.h'

    type(defect), pointer :: def, def_prev
    type(postPtr) :: postList(LEVELS)
    type(postDefect), pointer :: postDef, postDefPrev, postDefNew
    integer :: cell, level, numSame, i, task
    !<Used for communication between processors
    integer :: numSend(LEVELS), totalSend, totalRecv, count, defectType(SPECIES)
    integer, allocatable :: numRecv(:), sendBuff(:,:), recvBuff(:,:)
    integer :: status(MPI_STATUS_SIZE), sendRequest, sendStatus(MPI_STATUS_SIZE)
    !<number of point defects (self-interstitial atom, vacancy, solute atom, gas atom)
    integer :: numSIA, numVac, numCuAtom, numHeAtom
    !<total retained vacancies/self-interstitial atoms in the whole system
    integer :: totalSIA, totalVac
    !<total number of monomers in the clusters
    integer :: monI, monV, monS, monSV, monCu, monHe
    !<total number of various clusters in the whole system
    integer :: totalLoop, totalVoid, totalS, totalSV, totalCu, totalHe, totalHeV, totalHeVClu
    !<Number density of clusters
    double precision :: ndenLoop, ndenVoid, ndenS, ndenSV, ndenCu, ndenHeV, ndenHeVClu
    !<Average radius of clusters
    double precision :: radiusLoop, radiusVoid, radiusS, radiusSV, radiusCu, radiusHe
    !<Average size of clusters
    double precision :: sizeLoop, sizeVoid, sizeS, sizeSV, sizeCu, sizeHe
    !<Other variables
    double precision :: VRetained, VAnnihilated

    !<Test code
    type(postDefect), pointer :: postTest
    type(defect), pointer :: defTest

    interface
        subroutine findDefectInPostList(defectType, postDef, postDefPrev, numSame)
            use mod_constants
            use mod_structuretype
            implicit none

            integer, intent(in), dimension(SPECIES) :: defectType
            type(postDefect), pointer, intent(inout) :: postDef, postDefPrev
            integer, intent(inout) :: numSame
        end subroutine findDefectInPostList
    end interface

    defectType=0
    !<Add defects (cell 1) to postList
    numSend=0   !<number of types per level
    do level=1, LEVELS
        nullify(postList(level)%head)

        !<Add defects (cell 1) to postList
        !nullify(postDef)
        postDef => postList(level)%head
        nullify(def)
        def => myDRl(level,1)%defect_all
        do while(associated(def))
            numSend(level) = numSend(level) + 1
            nullify(postDefNew)
            allocate(postDefNew)
            allocate(postDefNew%defectType(SPECIES))
            postDefNew%defectType = def%defectType
            postDefNew%num = def%num
            nullify(postDefNew%next)
            if(.NOT. associated(postList(level)%head)) then
                postList(level)%head => postDefNew
                postDef => postList(level)%head
            else
                postDef%next => postDefNew
                postDef => postDef%next
            end if
            def=>def%next
        end do
    end do

    !<Add defects (cell 2~numMeshes) to postList
    do cell=2,  numMeshes
        do level=1, LEVELS
            !nullify(postDef)
            postDef => postList(level)%head
            nullify(postDefPrev)

            nullify(def)
            def=>myDRL(level,cell)%defect_all
            do while(associated(def))
                numSame = 0
                call findDefectInPostList(def%defectType, postDef, postDefPrev, numSame)
                if(associated(postDef)) then
                    if(numSame == SPECIES) then
                        postDef%num = postDef%num + def%num
                    else    !<insert the defect in front of postDef
                        numSend(level) = numSend(level) + 1

                        nullify(postDefNew)
                        allocate(postDefNew)
                        allocate(postDefNew%defectType(SPECIES))
                        postDefNew%defectType = def%defectType
                        postDefNew%num = def%num
                        postDefNew%next => postDef
                        if(.NOT. associated(postDefPrev)) then
                            postList(level)%head => postDefNew
                            postDefPrev => postList(level)%head
                        else
                            postDefPrev%next => postDefNew
                            postDefPrev => postDefPrev%next
                        end if
                    end if
                else    !<insert the defect at the end of the list
                    numSend(level) = numSend(level) + 1
                    nullify(postDefNew)
                    allocate(postDefNew)
                    allocate(postDefNew%defectType(SPECIES))
                    postDefNew%defectType = def%defectType
                    postDefNew%num = def%num
                    nullify(postDefNew%next)
                    if(.NOT. associated(postDefPrev)) then
                        postList(level)%head => postDefNew
                        postDefPrev => postList(level)%head
                        postDef => postDefPrev%next
                    else
                        postDefPrev%next => postDefNew
                        postDefPrev => postDefPrev%next
                        postDef => postDefPrev%next
                    end if
                end if
                def=>def%next
            end do
        end do
    end do

    totalSend=0
    do level=1, LEVELS
        totalSend=totalSend+numSend(level)  !<total number of types in this process
    end do
    allocate(numRecv(numtasks*LEVELS))
    call MPI_GATHER(numSend,LEVELS,MPI_INTEGER,numRecv,LEVELS,MPI_INTEGER,0,comm,ierr)

    !<Update send buffer
    if(myProc%taskid /= 0) then
        allocate(sendBuff(SPECIES+1,totalSend))
        i=0
        do level=1, LEVELS
            !nullify(postDef)
            postDef => postList(level)%head
            do while(associated(postDef))
                i=i+1
                sendBuff(1:SPECIES,i) = postDef%defectType(:)
                sendBuff(SPECIES+1,i) = postDef%num
                postDef => postDef%next
            end do
        end do
        !call MPI_ISEND(sendBuff, (SPECIES+1)*totalSend, MPI_INTEGER, 0, 400, comm, sendRequest, ierr)
        call MPI_SEND(sendBuff, (SPECIES+1)*totalSend, MPI_INTEGER, 0, 1000, comm, ierr)
        deallocate(sendBuff)
    else    !<MASTER processer
        do task=1, numtasks-1
            totalRecv=0
            do level=1, LEVELS
                totalRecv=totalRecv+numRecv(LEVELS*task+level)  !<total number of types in this process
            end do
            allocate(recvBuff(SPECIES+1,totalRecv))
            call MPI_RECV(recvBuff,(SPECIES+1)*totalRecv,MPI_INTEGER,task,1000,comm,status,ierr)

            count=0
            do level=1, LEVELS
                !<Determine the starting position of each level in recvBuff
                if(level==1) then
                    count=0
                else
                    count=count+numRecv(LEVELS*task+level-1)
                end if

                !nullify(postDef)
                postDef=>postList(level)%head
                nullify(postDefPrev)
                do i=1, numRecv(LEVELS*task+level)
                    defectType(:)=recvBuff(1:SPECIES,count+i)
                    numSame=0
                    call findDefectInPostList(defectType, postDef, postDefPrev, numSame)
                    if(associated(postDef)) then
                        if(numSame==SPECIES) then
                            postDef%num=postDef%num+recvBuff(SPECIES+1,count+i)
                        else
                            nullify(postDefNew)
                            allocate(postDefNew)
                            allocate(postDefNew%defectType(SPECIES))
                            postDefNew%defectType = defectType
                            postDefNew%num = recvBuff(SPECIES+1,count+i)
                            postDefNew%next => postDef
                            if(.NOT. associated(postDefPrev)) then
                                postList(level)%head => postDefNew
                                postDefPrev => postList(level)%head
                            else
                                postDefPrev%next => postDefNew
                                postDefPrev => postDefPrev%next
                            end if
                        end if
                    else
                        nullify(postDefNew)
                        allocate(postDefNew)
                        allocate(postDefNew%defectType(SPECIES))
                        postDefNew%defectType=defectType
                        postDefNew%num=recvBuff(SPECIES+1,count+i)
                        nullify(postDefNew%next)
                        if(.NOT. associated(postDefPrev)) then
                            postList(level)%head => postDefNew
                            postDefPrev => postList(level)%head
                            postDef => postDefPrev%next
                        else
                            postDefPrev%next => postDefNew
                            postDefPrev => postDefPrev%next
                            postDef => postDefPrev%next
                        end if
                    end if
                end do
            end do
            deallocate(recvBuff)
        end do
    end if

    !<Output totdat.out
    if(myProc%taskid==0) then
        !<number of point defects
        numSIA = 0
        numVac = 0
        numCuAtom = 0
        numHeAtom = 0
        !<total number of monomers in clusters
        monI=0
        monV=0
        monS=0
        monSV=0
        monCu=0
        monHe = 0
        !<total number of clusters in the whole system
        totalLoop=0
        totalVoid=0
        totalS=0
        totalSV=0
        totalCu=0

        totalHe = 0
        totalHeV = 0
        totalHeVClu = 0
        !<number density of clusters
        ndenLoop=0d0
        ndenVoid=0d0
        ndenS=0d0
        ndenSV=0d0
        ndenCu=0d0

        ndenHeV = 0d0
        ndenHeVClu = 0d0

        !<average radius of clusters
        radiusLoop=0d0
        radiusVoid=0d0
        radiusS=0d0
        radiusSV=0d0
        radiusCu=0d0
        !<average size of clusters
        sizeLoop=0d0
        sizeVoid=0d0
        sizeS=0d0
        sizeSV=0d0
        sizeCu=0d0

        !<Output time information
        if(totdatToggle=='yes') then
            write(TOTFILE,*)
            write(TOTFILE,*)
            write(TOTFILE,*) '***********************************************************************'
            if(irradiationType=='FrenkelPair') then
                write(TOTFILE,*) 'Time',time,'DPA',DPA,'FrenkelPair',totImpEvents,'HeImpEvents',totImpHe
            else if(irradiationType=='Cascade') then
                write(TOTFILE,*) 'Time',time,'DPA',DPA,'DPA_low',DPA_low,'Cascade',totImpEvents,'HeImpEvents',totImpHe
            end if
            write(TOTFILE,*) 'Defects (_ _ _ [defectType] _ [number]):'
        end if
        if(defectToggle=='yes') then
            write(DEFFILE,*)
            write(DEFFILE,*)
            write(DEFFILE,*) '***********************************************************************'
            if(irradiationType=='FrenkelPair') then
                write(DEFFILE,*) 'Time', time, 'DPA', DPA, 'FrenkelPair', totImpEvents, 'HeImpEvents',totImpHe
            else if(irradiationType=='Cascade') then
                write(DEFFILE,*) 'Time',time,'DPA',DPA,'DPA_low',DPA_low,'Cascade',totImpEvents,'HeImpEvents',totImpHe
            end if
            write(DEFFILE,*) 'Defects (_ _ _ [defectType] _ [number]):'
        end if
        if(stadatToggle=='yes') then
            write(STAFILE,*)
            write(STAFILE,*)
            write(STAFILE,*) '***********************************************************************'
            if(irradiationType=='FrenkelPair') then
                write(STAFILE,*) 'Time', time, 'DPA', DPA, 'FrenkelPair', totImpEvents,'HeImpEvents',totImpHe
            else if(irradiationType=='Cascade') then
                write(STAFILE,*) 'Time',time,'DPA',DPA,'DPA_low',DPA_low,'Cascade',totImpEvents,'HeImpEvents',totImpHe
            end if
        end if

        do level=1, LEVELS
            nullify(postDef)
            postDef => postList(level)%head
            do while(associated(postDef))
                !<Count SIA Loops
                if(postDef%defectType(1)>0 .OR. postDef%defectType(SPECIES)>0) then !SIA clusters
                    if(postDef%defectType(1)==1) then   !SIA
                        numSIA = postDef%num
                    end if
                    if(postDef%defectType(1) > minLoop) then
                        monI = monI + postDef%defectType(1) * postDef%num
                        totalLoop = totalLoop + postDef%num
                    else if(postDef%defectType(SPECIES) > minLoop) then
                        monI = monI + postDef%defectType(SPECIES) * postDef%num
                        totalLoop = totalLoop + postDef%num
                    end if
                end if

                !<Count Voids
                if(postDef%defectType(1) < 0) then !Void
                    if(postDef%defectType(1)==-1 .AND. postDef%defectType(2)==0) then !V_0_0
                        numVac = postDef%num
                    end if
                    if(iabs(postDef%defectType(1)) > minVoid) then
                        monV = monV + iabs(postDef%defectType(1)) * postDef%num
                        totalVoid = totalVoid + postDef%num
                    end if
                end if

                !<Count Solute
                if(soluteConc>0d0) then
                    if(postDef%defectType(2) > minS) then
                        monCu = monCu + postDef%defectType(2) * postDef%num
                        totalCu = totalCu + postDef%num
                    end if
                    if(postDef%defectType(1)==0 .AND. postDef%defectType(2)==1) then
                        numCuAtom = postDef%num
                    else if(postDef%defectType(1)==0 .AND. postDef%defectType(2)>minS) then
                        monS = monS + postDef%defectType(2) * postDef%num
                        totalS = totalS + postDef%num
                    else if(postDef%defectType(1)<0 .AND. iabs(postDef%defectType(1))+postDef%defectType(2)>minSV) then
                        monSV = monSV + (iabs(postDef%defectType(1))+postDef%defectType(2))*postDef%num
                        totalSV = totalSV + postDef%num
                    end if
                end if

                !<Count He
                if(HeDPAratio>0d0) then
                    if(postDef%defectType(2) /= 0) then
                        monHe = monHe + postDef%defectType(2) * postDef%num
                        totalHe = totalHe + postDef%num
                    end if
                    if(postDef%defectType(1)==0 .AND. postDef%defectType(2)==1) then
                        numHeAtom = postDef%num
                    else if(postDef%defectType(1)==-1 .AND. postDef%defectType(2)==1) then
                        totalHeV = totalHeV + postDef%num
                    else if(postDef%defectType(1)<0 .AND. postDef%defectType(2)>0) then
                        totalHeVClu = totalHeVClu + postDef%num
                    end if
                end if

                !<output defects
                if(totdatToggle=='yes') then
                    write(TOTFILE,*) postDef%defectType, postDef%num
                end if
                if(defectToggle=='yes') then
                    write(DEFFILE,*) postDef%defectType, postDef%num
                end if
                postDef=>postDef%next
            end do
        end do

        !<self-defects
        ndenLoop = dble(totalLoop)/(systemVolume*1d-27)     !m-3
        ndenVoid = dble(totalVoid)/(systemVolume*1d-27)     !m-3
        radiusLoop = ((dble(monI)/dble(totalLoop))*atomVol/(PI*burgers))**(1d0/2d0)
        radiusVoid = (3d0*(dble(monV)/dble(totalVoid))*atomVol/(4d0*PI))**(1d0/3d0)
        sizeLoop = dble(monI)/dble(totalLoop)
        sizeVoid = dble(monV)/dble(totalVoid)

        if(soluteConc > 0d0) then
            ndenS = dble(totalS)/(systemVolume*1d-27) !m-3
            ndenSV = dble(totalSV)/(systemVolume*1d-27)         !m-3
            ndenCu = dble(totalCu)/(systemVolume*1d-27)         !m-3

            radiusS = (3d0*(dble(monS)/dble(totalS))*atomVol/(4d0*PI))**(1d0/3d0)
            radiusSV = (3d0*(dble(monSV)/dble(totalSV))*atomVol/(4d0*PI))**(1d0/3d0)
            radiusCu = (3d0*(dble(monCu)/dble(totalCu))*atomVol/(4d0*PI))**(1d0/3d0)

            sizeS = dble(monS)/dble(totalS)
            sizeSV = dble(monSV)/dble(totalSV)
            sizeCu = dble(monCu)/dble(totalCu)
        end if

        if(HeDPAratio > 0d0) then
            ndenHeV = dble(totalHeV)/(systemVolume*1d-27) !m-3
            ndenHeVClu = dble(totalHeVClu)/(systemVolume*1d-27) !m-3
        end if

        if(totdatToggle=='yes') then
            write(TOTFILE,*) 'SIAs    ',numSIA,' numDensity(m-3)',dble(numSIA)/(systemVolume*1d-27)
            write(TOTFILE,*) 'Vac     ',numVac,' numDensity(m-3)',dble(numVac)/(systemVolume*1d-27)
            write(TOTFILE,*) 'Loops   ', totalLoop, ' numDensity(m-3)', ndenLoop, ' aveRadius(nm)', radiusLoop
            write(TOTFILE,*) 'Voids   ', totalVoid, ' numDensity(m-3)', ndenVoid, ' aveRadius(nm)', radiusVoid
            if(soluteConc > 0d0) then
                write(TOTFILE,*) 'CuClus  ',totalS,' numDensity(m-3)',ndenS,' aveRadius(nm)',radiusS
                write(TOTFILE,*) 'CuVClus ',totalSV,' numDensity(m-3)',ndenSV,' aveRadius(nm)',radiusSV
                write(TOTFILE,*) 'Precip  ', totalS+totalSV,' numDensity(m-3)',&
                        dble(totalS+totalSV)/systemVolume*1d27,' aveRadius(nm)',&
                        (3*(dble(monS+monSV)/dble(totalS+totalSV))*atomVol/(4d0*PI))**(1d0/3d0),&
                        '(V in CuV is taken into account)'
                write(TOTFILE,*) 'Precip  ',totalCu,' numDensity(m-3)',ndenCu,' aveRadius(nm)',radiusCu,&
                        '(V in CuV is not taken into account)'
            end if
            if(HeDPAratio > 0d0) then
                write(TOTFILE,*) 'HeAtom  ',numHeAtom,' numDensity(m-3)',dble(numHeAtom)/(systemVolume*1d-27)
                write(TOTFILE,*) 'He      ',totalHe,' numDensity(m-3)',dble(totalHe)/(systemVolume*1d-27)
                write(TOTFILE,*) 'HeV     ', totalHeV,' numDensity(m-3)', ndenHeV
                write(TOTFILE,*) 'HeVClus ',totalHeVClu,' numDensity(m-3)',ndenHeVClu
            end if
        end if

        if(stadatToggle=='yes') then
            write(STAFILE,*) 'SIAs    ',numSIA,' numDensity(m-3)',dble(numSIA)/(systemVolume*1d-27)
            write(STAFILE,*) 'Vac     ',numVac,' numDensity(m-3)',dble(numVac)/(systemVolume*1d-27)
            write(STAFILE,*) 'Loops   ', totalLoop, ' numDensity(m-3)', ndenLoop, ' aveRadius(nm)', radiusLoop
            write(STAFILE,*) 'Voids   ', totalVoid, ' numDensity(m-3)', ndenVoid, ' aveRadius(nm)', radiusVoid
            if(soluteConc > 0d0) then
                write(STAFILE,*) 'CuClus  ',totalS,' numDensity(m-3)',ndenS,' aveRadius(nm)',radiusS
                write(STAFILE,*) 'CuVClus ',totalSV,' numDensity(m-3)',ndenSV,' aveRadius(nm)',radiusSV
                write(STAFILE,*) 'Precip  ', totalS+totalSV,' numDensity(m-3)',&
                        dble(totalS+totalSV)/systemVolume*1d27,' aveRadius(nm)',&
                        (3*(dble(monS+monSV)/dble(totalS+totalSV))*atomVol/(4d0*PI))**(1d0/3d0),&
                        '(V in CuV is taken into account)'
                write(STAFILE,*) 'Precip  ',totalCu,' numDensity(m-3)',ndenCu,' aveRadius(nm)',radiusCu,&
                        '(V in CuV is not taken into account)'
            end if
            if(HeDPAratio > 0d0) then
                write(STAFILE,*) 'HeAtom  ',numHeAtom,' numDensity(m-3)',dble(numHeAtom)/(systemVolume*1d-27)
                write(STAFILE,*) 'He      ',totalHe,' numDensity(m-3)',dble(totalHe)/(systemVolume*1d-27)
                write(STAFILE,*) 'HeV     ', totalHeV,' numDensity(m-3)', ndenHeV
                write(STAFILE,*) 'HeVClus ',totalHeVClu,' numDensity(m-3)',ndenHeVClu
            end if
        end if
    end if

    !<Deallocate postList
    do level=1, LEVELS
        !nullify(postDef)
        postDef => postList(level)%head
        do while(associated(postDef))
            postList(level)%head => postDef%next
            deallocate(postDef%defectType)
            deallocate(postDef)
            postDef=>postList(level)%head
        end do
    end do
    deallocate(numRecv)

    !if(myProc%taskid /= 0) then
    !    call MPI_WAIT(sendRequest, sendStatus, ierr)
    !    deallocate(sendBuff)
    !end if

end subroutine

!*******************************************************************************************************
!> Subroutine outputTotalDefects_alternative: output the total defects and post processing for the entire volume
! NOTE: in this subroutine, each process counts its own Void, Loop, Precipitates, etc.,
! and then sends them to process 0, which finally sum and outputs.
!*******************************************************************************************************
subroutine outputTotalDefects_alternative()
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none
    include 'mpif.h'

    type(defect), pointer :: def, def_prev
    type(postDefect), pointer :: postList(:)
    type(postDefect), pointer :: postDef, postDefPrev, postDefNew
    integer :: cell, level, numSame, i, task
    !<Used for communication between processors
    integer :: numSend(LEVELS), totalSend, count, defectType(SPECIES), maxRecv
    integer, allocatable :: totalRecv(:)
    integer, allocatable :: numRecv(:), sendBuff(:,:), recvBuff(:,:,:)
    integer :: status(MPI_STATUS_SIZE), sendRequest, sendStatus(MPI_STATUS_SIZE), recvStatus(MPI_STATUS_SIZE)
    integer, allocatable :: recvRequest(:)
    !<I, V, S, SV
    integer :: reduce(SPECIES+9), totalReduce(SPECIES+9)
    !<number of point defects (self-interstitial atom, vacancy, solute atom, gas atom)
    integer :: totalPointDefect(SPECIES+1)
    !<total retained vacancies/self-interstitial atoms in the whole system
    integer :: totalSIA, totalVac
    !<total number of monomers in the clusters
    integer :: numI, numV, numS, numSV
    integer :: level_I, level_V, level_S, level_SV  !<Levels for statistics
    !<total number of various clusters in the whole system
    integer :: totalLoop, totalVoid, totalSolute, totalSV
    !<Number density of clusters
    double precision :: ndenLoop, ndenVoid, ndenSolute, ndenSV
    !<Average radius of clusters
    double precision :: radiusLoop, radiusVoid, radiusSolute, radiusSV
    !<Average size of clusters
    double precision :: sizeLoop, sizeVoid, sizeSolute, sizeSV
    !<Other variables
    double precision :: VRetained, VAnnihilated, conPointVac, conPointSIA
    !<Function
    integer, external :: find_level

    !<Test code
    type(postDefect), pointer :: postTest
    type(defect), pointer :: defTest

    interface
        subroutine findDefectInPostList(defectType, postDef, postDefPrev, numSame)
            use mod_constants
            use mod_structuretype
            implicit none

            integer, intent(in), dimension(SPECIES) :: defectType
            type(postDefect), pointer, intent(inout) :: postDef, postDefPrev
            integer, intent(inout) :: numSame
        end subroutine findDefectInPostList
    end interface

    defectType=0
    !<Initialize postList
    allocate(postList(LEVELS))
    !<Add defects (cell 1) to postList
    numSend=0   !<number of types per level
    do level=1, LEVELS
        !<initialize the first node
        postList(level)%num=0
        nullify(postList(level)%next)

        !<Add defects (cell 1) to postList
        nullify(postDefPrev)
        postDefPrev=>postList(level)
        nullify(def)
        def=>myDRl(level,1)%defect_all
        do while(associated(def))
            nullify(postDefNew)
            allocate(postDefNew)
            allocate(postDefNew%defectType(SPECIES))
            postDefNew%defectType=def%defectType
            postDefNew%num=def%num
            nullify(postDefNew%next)
            postDefPrev%next=>postDefNew
            postDefPrev=>postDefPrev%next
            numSend(level)=numSend(level)+1
            def=>def%next
        end do
    end do

    !<Add defects (cell 2~numMeshes) to postList
    do cell=2,  numMeshes
        do level=1, LEVELS
            nullify(postDef)
            postDef=>postList(level)
            nullify(postDefPrev)
            postDefPrev=>postDef
            postDef=>postDef%next

            nullify(def)
            def=>myDRL(level,cell)%defect_all
            do while(associated(def))
                numSame=0
                call findDefectInPostList(def%defectType, postDef, postDefPrev, numSame)
                if(associated(postDef)) then
                    if(numSame==SPECIES) then
                        postDef%num=postDef%num+def%num
                    else    !<insert the defect in front of postDef
                        nullify(postDefNew)
                        allocate(postDefNew)
                        allocate(postDefNew%defectType(SPECIES))
                        postDefNew%defectType=def%defectType
                        postDefNew%num=def%num
                        postDefNew%next=>postDef
                        postDefPrev%next=>postDefNew
                        postDefPrev=>postDefPrev%next
                        numSend(level)=numSend(level)+1
                    end if
                else    !<insert the defect at the end of the list
                    nullify(postDefNew)
                    allocate(postDefNew)
                    allocate(postDefNew%defectType(SPECIES))
                    postDefNew%defectType=def%defectType
                    postDefNew%num=def%num
                    nullify(postDefNew%next)
                    postDefPrev%next=>postDefNew
                    postDefPrev=>postDefPrev%next
                    numSend(level)=numSend(level)+1
                end if
                def=>def%next
            end do
        end do
    end do

    totalSend=0
    do level=1, LEVELS
        totalSend=totalSend+numSend(level)  !<total number of types in this process
    end do
    allocate(numRecv(numtasks*LEVELS))
    call MPI_GATHER(numSend,LEVELS,MPI_INTEGER,numRecv,LEVELS,MPI_INTEGER,0,comm,ierr)

    if(myProc%taskid /= 0) then
        allocate(sendBuff(SPECIES+1,totalSend))
        i=0
        do level=1, LEVELS
            nullify(postDef)
            postDef=>postList(level)%next
            do while(associated(postDef))
                i=i+1
                sendBuff(1:SPECIES,i)=postDef%defectType(:)
                sendBuff(SPECIES+1,i)=postDef%num
                postDef=>postDef%next
            end do
        end do
        call MPI_ISEND(sendBuff, (SPECIES+1)*totalSend, MPI_INTEGER, 0, 400+myProc%taskid, comm, sendRequest, ierr)
    else    !<process 0
        allocate(totalRecv(numtasks-1))
        totalRecv=0
        maxRecv=0
        do task=1, numtasks-1
            do level=1, LEVELS
                totalRecv(task)=totalRecv(task)+numRecv(LEVELS*task+level)  !<total number of types in this process
            end do
            if(totalRecv(task) >= maxRecv) then
                maxRecv=totalRecv(task)
            end if
        end do

        allocate(recvRequest(numtasks-1))
        allocate(recvBuff(SPECIES+1, maxRecv, numtasks-1))
        do task=1, numtasks-1
            call MPI_IRECV(recvBuff(1,1,task),(SPECIES+1)*totalRecv(task),MPI_INTEGER,task, &
                    400+task,comm,recvRequest(task),ierr)
        end do
    end if

    !<Statistical defects
    !<number of point defects
    totalPointDefect=0
    !<total number of monomers in clusters
    numI=0
    numV=0
    numS=0
    numSV=0
    !<total number of clusters in the whole system
    totalLoop=0
    totalVoid=0
    totalSolute=0
    totalSV=0
    !<number density of clusters
    ndenLoop=0d0
    ndenVoid=0d0
    ndenSolute=0d0
    ndenSV=0d0
    !<average radius of clusters
    radiusLoop=0d0
    radiusVoid=0d0
    radiusSolute=0d0
    radiusSV=0d0
    !<average size of clusters
    sizeLoop=0d0
    sizeVoid=0d0
    sizeSolute=0d0
    sizeSV=0d0

    conPointVac=0d0
    conPointSIA=0d0

    !<determine the level
    defectType=(/1,0,0/)
    level_I=find_level(defectType)
    defectType=(/-1,0,0/)
    level_V=find_level(defectType)
    defectType=(/0,1,0/)
    level_S=find_level(defectType)
    defectType=(/-1,1,0/)
    level_SV=find_level(defectType)

    do level=1, LEVELS
        nullify(postDef)
        postDef=>postList(level)%next
        do while(associated(postDef))
            if(level==level_I) then   !<Loop: I_0_0
                if(postDef%defectType(1)==1) then !I_0_0
                    totalPointDefect(1)=postDef%num
                end if
                if(postDef%defectType(1)>minLoop) then
                    numI=numI+postDef%defectType(1)*postDef%num
                    totalLoop=totalLoop+postDef%num
                end if
            else if(level==level_V) then !<Void: V_0_0
                if(postDef%defectType(1)==-1) then !V_0_0
                    totalPointDefect(2)=postDef%num
                end if
                if(iabs(postDef%defectType(1))>minVoid) then
                    numV=numV+iabs(postDef%defectType(1))*postDef%num
                    totalVoid=totalVoid+postDef%num
                end if
            else if(level==level_S) then !<S clusters: 0_S_0
                if(postDef%defectType(2)==1) then !0_S_0
                    totalPointDefect(3)=postDef%num
                end if
                if(postDef%defectType(2)>minS) then
                    numS=numS+postDef%defectType(2)*postDef%num
                    totalSolute=totalSolute+postDef%num
                end if
            else if(level==level_SV) then    !<V_S_0 clusters
                if(iabs(postDef%defectType(1))+postDef%defectType(2)>minSV) then
                    !numSV=numSV+max(iabs(postDef%defectType(1)),postDef%defectType(2))%postDef%num
                    numSV=numSV+(iabs(postDef%defectType(1))+postDef%defectType(2))*postDef%num
                    totalSV=totalSV+postDef%num
                end if
            end if
            postDef=>postDef%next
        end do
    end do

    !<Update send buffer
    if(myProc%taskid == 0) then
        do task=1, numtasks-1
            call MPI_WAIT(recvRequest(task), recvStatus, ierr)
            count=0
            do level=1, LEVELS
                !<Determine the starting position of each level in recvBuff
                if(level==1) then
                    count=0
                else
                    count=count+numRecv(LEVELS*task+level-1)
                end if

                nullify(postDef)
                postDef=>postList(level)
                nullify(postDefPrev)
                postDefPrev=>postDef
                postDef=>postDef%next

                do i=1, numRecv(LEVELS*task+level)
                    defectType(:)=recvBuff(1:SPECIES,count+i,task)
                    numSame=0
                    call findDefectInPostList(defectType, postDef, postDefPrev, numSame)
                    if(associated(postDef)) then
                        if(numSame==SPECIES) then
                            postDef%num=postDef%num+recvBuff(SPECIES+1,count+i,task)
                        else
                            nullify(postDefNew)
                            allocate(postDefNew)
                            allocate(postDefNew%defectType(SPECIES))
                            postDefNew%defectType=defectType
                            postDefNew%num=recvBuff(SPECIES+1,count+i,task)
                            postDefNew%next=>postDef
                            postDefPrev%next=>postDefNew
                            postDefPrev=>postDefPrev%next
                        end if
                    else
                        nullify(postDefNew)
                        allocate(postDefNew)
                        allocate(postDefNew%defectType(SPECIES))
                        postDefNew%defectType=defectType
                        postDefNew%num=recvBuff(SPECIES+1,count+i,task)
                        nullify(postDefNew%next)
                        postDefPrev%next=>postDefNew
                        postDefPrev=>postDefPrev%next
                    end if
                end do
            end do
        end do
        deallocate(recvBuff)
    end if

    !<Output totdat.out
    if(myProc%taskid==0) then
        !<Output time information
        if(totdatToggle=='yes') then
            write(TOTFILE,*)
            write(TOTFILE,*)
            write(TOTFILE,*) '***********************************************************************'
            write(TOTFILE,*) 'Time', time, 'DPA', DPA, outputToggle, totImpEvents
            write(TOTFILE,*) 'Defects (_ _ _ [defectType] _ [number]):'
        end if
        if(defectToggle=='yes') then
            write(DEFFILE,*)
            write(DEFFILE,*)
            write(DEFFILE,*) '***********************************************************************'
            write(DEFFILE,*) 'Time', time, 'DPA', DPA, outputToggle, totImpEvents
            write(DEFFILE,*) 'Defects (_ _ _ [defectType] _ [number]):'
        end if
        if(stadatToggle=='yes') then
            write(STAFILE,*)
            write(STAFILE,*)
            write(STAFILE,*) '***********************************************************************'
            write(STAFILE,*) 'Time', time, 'DPA', DPA, outputToggle, totImpEvents
        end if

        do level=1, LEVELS
            nullify(postDef)
            postDef=>postList(level)%next
            do while(associated(postDef))
                !<output defects
                if(totdatToggle=='yes') then
                    write(TOTFILE,*) postDef%defectType, postDef%num
                end if
                if(defectToggle=='yes') then
                    write(DEFFILE,*) postDef%defectType, postDef%num
                end if
                postDef=>postDef%next
            end do
        end do
    end if

    !<Output statistics
    reduce=0
    reduce(1:SPECIES+1)=totalPointDefect(:)
    reduce(SPECIES+2)=numI
    reduce(SPECIES+3)=numV
    reduce(SPECIES+4)=numS
    reduce(SPECIES+5)=numSV
    reduce(SPECIES+6)=totalLoop
    reduce(SPECIES+7)=totalVoid
    reduce(SPECIES+8)=totalSolute
    reduce(SPECIES+9)=totalSV
    totalReduce=0
    call MPI_REDUCE(reduce, totalReduce, SPECIES+9, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

    if(myProc%taskid==0) then
        !<number density of clusters
        ndenLoop=dble(totalReduce(SPECIES+6))/systemVolume
        ndenVoid=dble(totalReduce(SPECIES+7))/systemVolume
        ndenSolute=dble(totalReduce(SPECIES+8))/systemVolume
        ndenSV=dble(totalReduce(SPECIES+9))/systemVolume
        !<average radius of clusters
        radiusLoop=((dble(totalReduce(SPECIES+2))/dble(totalReduce(SPECIES+6)))*atomVol/(PI*burgers))**(1d0/2d0)
        radiusVoid=(3d0*(dble(totalReduce(SPECIES+3))/dble(totalReduce(SPECIES+7)))*atomVol/(4d0*PI))**(1d0/3d0)
        radiusSolute=(3*(dble(totalReduce(SPECIES+4))/dble(totalReduce(SPECIES+8)))*atomVol/(4*PI))**(1d0/3d0)
        radiusSV=(3*(dble(totalReduce(SPECIES+5))/dble(totalReduce(SPECIES+9)))*atomVol/(4*PI))**(1d0/3d0)
        !<average size of clusters
        sizeLoop=dble(totalReduce(SPECIES+2))/dble(totalReduce(SPECIES+6))
        sizeVoid=dble(totalReduce(SPECIES+3))/dble(totalReduce(SPECIES+7))
        sizeSolute=dble(totalReduce(SPECIES+4))/dble(totalReduce(SPECIES+8))
        sizeSV=dble(totalReduce(SPECIES+5))/dble(totalReduce(SPECIES+9))

        conPointSIA=dble(totalReduce(1))/systemVolume*atomVol
        conPointVac=dble(totalReduce(2))/systemVolume*atomVol

        if(totdatToggle=='yes') then
            write(TOTFILE,*) 'numClusters (Loop, Void, Precipitates, Solutes, V_S clusters):'
            write(TOTFILE,*) totalReduce(SPECIES+6), totalReduce(SPECIES+7), &
                    (totalReduce(SPECIES+8)+totalReduce(SPECIES+9)), totalReduce(SPECIES+8), totalReduce(SPECIES+9)
            write(TOTFILE,*) 'NumberDensity(m-3) (Loop, Void, Precipitates, Solutes, V_S clusters):'
            write(TOTFILE,*) ndenLoop*1d27, ndenVoid*1d27, &
                    dble(totalReduce(SPECIES+8)+totalReduce(SPECIES+9))/systemVolume*1d27, ndenSolute*1d27, ndenSV*1d27
            write(TOTFILE,*) 'Concentration(atom-1) (Loop, Void, Precipitates, Solutes, V_S clusters):'
            write(TOTFILE,*) ndenLoop*atomVol, ndenVoid*atomVol, &
                    dble(totalReduce(SPECIES+8)+totalReduce(SPECIES+9))/systemVolume*atomVol, &
                    ndenSolute*atomVol, ndenSV*atomVol
            write(TOTFILE,*) 'AverageRadius(nm) (Loop, Void, Precipitates, Solutes, V_S clusters):'
            write(TOTFILE,*) radiusLoop, radiusVoid, (3*(dble(totalReduce(SPECIES+4)+totalReduce(SPECIES+5))/&
                            dble(totalReduce(SPECIES+8)+totalReduce(SPECIES+9)))*atomVol/(4*PI))**(1d0/3d0),&
                    radiusSolute,radiusSV
            write(TOTFILE,*) 'AverageSize (Loop, Void, Precipitates, Solutes, V_S clusters):'
            write(TOTFILE,*) sizeLoop, sizeVoid, dble(totalReduce(SPECIES+4)+totalReduce(SPECIES+5))/&
                    dble(totalReduce(SPECIES+8)+totalReduce(SPECIES+9)), sizeSolute, sizeSV
            write(TOTFILE,*) 'ConcenPointDefects (SIA, V):'
            write(TOTFILE,*) conPointSIA, conPointVac
        end if
        if(stadatToggle=='yes') then
            write(STAFILE,*) 'numClusters (Loop, Void, Precipitates, Solutes, V_S clusters):'
            write(STAFILE,*) totalReduce(SPECIES+6), totalReduce(SPECIES+7), &
                    (totalReduce(SPECIES+8)+totalReduce(SPECIES+9)), totalReduce(SPECIES+8), totalReduce(SPECIES+9)
            write(STAFILE,*) 'NumberDensity(m-3) (Loop, Void, Precipitates, Solutes, V_S clusters):'
            write(STAFILE,*) ndenLoop*1d27, ndenVoid*1d27, &
                    dble(totalReduce(SPECIES+8)+totalReduce(SPECIES+9))/systemVolume*1d27, ndenSolute*1d27, ndenSV*1d27
            write(STAFILE,*) 'Concentration(atom-1) (Loop, Void, Precipitates, Solutes, V_S clusters):'
            write(STAFILE,*) ndenLoop*atomVol, ndenVoid*atomVol, &
                    dble(totalReduce(SPECIES+8)+totalReduce(SPECIES+9))/systemVolume*atomVol, &
                    ndenSolute*atomVol, ndenSV*atomVol
            write(STAFILE,*) 'AverageRadius(nm) (Loop, Void, Precipitates, Solutes, V_S clusters):'
            write(STAFILE,*) radiusLoop, radiusVoid, (3*(dble(totalReduce(SPECIES+4)+totalReduce(SPECIES+5))/&
                    dble(totalReduce(SPECIES+8)+totalReduce(SPECIES+9)))*atomVol/(4*PI))**(1d0/3d0),&
                    radiusSolute,radiusSV
            write(STAFILE,*) 'AverageSize (Loop, Void, Precipitates, Solutes, V_S clusters):'
            write(STAFILE,*) sizeLoop, sizeVoid, dble(totalReduce(SPECIES+4)+totalReduce(SPECIES+5))/&
                    dble(totalReduce(SPECIES+8)+totalReduce(SPECIES+9)), sizeSolute, sizeSV
            write(STAFILE,*) 'ConcenPointDefects (SIA, V):'
            write(STAFILE,*) conPointSIA, conPointVac
        end if

    end if

    !<Deallocate postList
    do level=1, LEVELS
        postDef=>postList(level)
        nullify(postDefPrev)
        do while(associated(postDef))
            postDefPrev=>postDef
            postDef=>postDef%next
            if(allocated(postDefPrev%defectType)) deallocate(postDefPrev%defectType)
            nullify(postDefPrev)
        end do
    end do
    nullify(postDef)
    deallocate(postList)
    deallocate(totalRecv)
    deallocate(numRecv)
    deallocate(recvRequest)

    if(myProc%taskid /= 0) then
        call MPI_WAIT(sendRequest, sendStatus, ierr)
        deallocate(sendBuff)
    end if

end subroutine

!*******************************************************************************************************
!> Subroutine find defect in post defect list:
!*******************************************************************************************************
subroutine outputVTK(fileNumber)
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none
    include 'mpif.h'

    integer, intent(in) :: fileNumber
    integer :: i, j, k, xGindex, yGindex, zGindex, cell, level, temp, task
    !total number of SIAs, Vs, Cu atoms in the system
    integer, allocatable :: totDefects(:,:,:,:), numDefects(:)
    character(len=14) :: fileName
    type(defect), pointer :: def, def_prev
    integer, allocatable :: numCells(:), sendBuff(:,:), recvBuff(:,:)
    integer :: status(MPI_STATUS_SIZE)

    if(myProc%taskid == 0) then
        if(soluteConc>0d0 .OR. HeDPAratio>0d0) then
            allocate(totDefects(4,numxGlobal,numyGlobal,numzGlobal))    !I, V, Cu/He, CuClu/HeClu
        else
            allocate(totDefects(2,numxGlobal,numyGlobal,numzGlobal))    !I, V
        end if
    else
        if(soluteConc>0d0 .OR. HeDPAratio>0d0) then
            allocate(sendBuff(numMeshes,7))    !x,y,z index of the cell, and number of I,V,Cu atoms
        else
            allocate(sendBuff(numMeshes,5))    !x,y,z index of the cell, and number of I,V
        end if
    end if

    if(soluteConc>0d0 .OR. HeDPAratio>0d0) then
        allocate(numDefects(4))    !I, V, Cu/He, CuClu/HeClu
    else
        allocate(numDefects(2))    !I, V
    end if

    do cell=1, numMeshes
        !indetify the global x, y, z index of this cell
        if(mod(myMesh(cell)%gid, numyGlobal*numzGlobal) /= 0) then
            xGindex = myMesh(cell)%gid / (numyGlobal*numzGlobal) + 1
        else
            xGindex = myMesh(cell)%gid / (numyGlobal*numzGlobal)
        end if
        temp = myMesh(cell)%gid - numyGlobal*numzGlobal*(xGindex-1)
        if(mod(temp, numzGlobal) /= 0) then
            yGindex = temp/numzGlobal + 1
        else
            yGindex = temp/numzGlobal
        end if
        zGindex = temp - numzGlobal*(yGindex-1)

        !< count number of SIAs, Vs, Cu atoms in this cell
        numDefects = 0

        do level=1, LEVELS
            nullify(def)
            def => myDRL(level, cell)%defect_all
            do while(associated(def))
                if(def%defectType(1) > 0) then
                    numDefects(1) = numDefects(1) + def%defectType(1) * def%num
                else if(def%defectType(SPECIES) > 0) then
                    numDefects(1) = numDefects(1) + def%defectType(SPECIES) * def%num
                end if

                if(def%defectType(1) < 0) then
                    numDefects(2) = numDefects(2) + iabs(def%defectType(1)) * def%num
                end if

                if(soluteConc > 0d0) then
                    if(def%defectType(2) > 0) then
                        numDefects(3) = numDefects(3) + def%defectType(2) * def%num
                        if(def%defectType(2) > minS) then
                            numDefects(4) = numDefects(4) + def%num
                        end if
                    end if
                end if

                if(HeDPAratio > 0d0) then
                    if(def%defectType(2) > 0) then
                        numDefects(3) = numDefects(3) + def%defectType(2) * def%num
                        numDefects(4) = numDefects(4) + def%num
                    end if
                end if

                def => def%next
            end do
        end do

        if(myProc%taskid == 0) then
            totDefects(:,xGindex,yGindex,zGindex) = numDefects(:)
        else
            sendBuff(cell,1) = xGindex
            sendBuff(cell,2) = yGindex
            sendBuff(cell,3) = zGindex
            sendBuff(cell,4:5) = numDefects(1:2)
            if(soluteConc>0d0 .OR. HeDPAratio>0d0) then
                sendBuff(cell,6:7) = numDefects(3:4)
            end if
        end if
    end do

    allocate(numCells(numtasks))
    call MPI_GATHER(numMeshes,1,MPI_INTEGER,numCells,1,MPI_INTEGER,0,comm,ierr)

    if(myProc%taskid /= 0) then !send
        if(soluteConc>0d0 .OR. HeDPAratio>0d0) then
            call MPI_SEND(sendBuff, numMeshes*7, MPI_INTEGER, 0, 99, comm, ierr)
        else
            call MPI_SEND(sendBuff, numMeshes*5, MPI_INTEGER, 0, 99, comm, ierr)
        end if
        deallocate(sendBuff)
    else    !0 proc, recv
        do task=1, numtasks-1
            if(soluteConc>0d0 .OR. HeDPAratio>0d0) then
                allocate(recvBuff(numCells(task+1),7))
                call MPI_RECV(recvBuff, numCells(task+1)*7, MPI_INTEGER, task, 99, comm,status, ierr)
            else
                allocate(recvBuff(numCells(task+1),5))
                call MPI_RECV(recvBuff, numCells(task+1)*5, MPI_INTEGER, task, 99, comm,status, ierr)
            end if

            do cell=1, numCells(task+1)
                if(soluteConc>0d0 .OR. HeDPAratio>0d0) then
                    totDefects(1:4,recvBuff(cell,1),recvBuff(cell,2),recvBuff(cell,3)) = recvBuff(cell,4:7)
                else
                    totDefects(1:2,recvBuff(cell,1),recvBuff(cell,2),recvBuff(cell,3)) = recvBuff(cell,4:5)
                end if
            end do
            deallocate(recvBuff)
        end do

        !< write VTK file
        fileName(1:7) = 'defVTK_'
        write(unit=fileName(8:10), fmt='(I3.3)') fileNumber
        fileName(11:14) = '.vtk'
        open(VTKFILE, file=fileName, action='write', status='Unknown')

        write(VTKFILE,'(a)') '# vtk DataFile Version 1.0'
        write(VTKFILE,'(a)') '3D Structured Grid of Linear Cubes'
        write(VTKFILE,'(a)') 'ASCII'
        write(VTKFILE,*)
        write(VTKFILE,'(a)') 'DATASET STRUCTURED_POINTS'
        write(VTKFILE,'(a)', advance='no') 'DIMENSIONS '
        write(VTKFILE,'(I4,I4,I4)') numxGlobal,numyGlobal,numzGlobal
        write(VTKFILE,'(a)', advance='no') 'ORIGIN '
        write(VTKFILE,'(I4,I4,I4)') 0, 0, 0
        write(VTKFILE,'(a)', advance='no') 'SPACING '
        write(VTKFILE,*) meshLength, meshLength, meshLength
        write(VTKFILE,'(a)', advance='no') 'POINT_DATA '
        write(VTKFILE,*) totalMeshes

        write(VTKFILE,'(a,a,a,I4)') 'SCALARS ', 'SIAs', ' float', 1
        write(VTKFILE,'(a,a)') 'LOOKUP_TABLE ', 'DEFAULT'
        do k=1, numzGlobal
            do j=1, numyGlobal
                do i=1, numxGlobal
                    write(VTKFILE,*) totDefects(1,i,j,k)
                end do
            end do
        end do

        write(VTKFILE,'(a,a,a,I4)') 'SCALARS ', 'Vacancies', ' float', 1
        write(VTKFILE,'(a,a)') 'LOOKUP_TABLE ', 'DEFAULT'
        do k=1, numzGlobal
            do j=1, numyGlobal
                do i=1, numxGlobal
                    write(VTKFILE,*) totDefects(2,i,j,k)
                end do
            end do
        end do

        if(soluteConc > 0d0) then
            write(VTKFILE,'(a,a,a,I4)') 'SCALARS ', 'Cu', ' float', 1
        else if(HeDPAratio > 0d0) then
            write(VTKFILE,'(a,a,a,I4)') 'SCALARS ', 'He', ' float', 1
        end if
        if(soluteConc>0d0 .OR. HeDPAratio>0d0) then
            write(VTKFILE,'(a,a)') 'LOOKUP_TABLE ', 'DEFAULT'
            do k=1, numzGlobal
                do j=1, numyGlobal
                    do i=1, numxGlobal
                        write(VTKFILE,*) totDefects(3,i,j,k)
                    end do
                end do
            end do
        end if

        if(soluteConc > 0d0) then
            write(VTKFILE,'(a,a,a,I4)') 'SCALARS ', 'CuClusters', ' float', 1
        else if(HeDPAratio > 0d0) then
            write(VTKFILE,'(a,a,a,I4)') 'SCALARS ', 'HeClusters', ' float', 1
        end if
        if(soluteConc>0d0 .OR. HeDPAratio>0d0) then
            write(VTKFILE,'(a,a)') 'LOOKUP_TABLE ', 'DEFAULT'
            do k=1, numzGlobal
                do j=1, numyGlobal
                    do i=1, numxGlobal
                        write(VTKFILE,*) totDefects(4,i,j,k)
                    end do
                end do
            end do
        end if

        write(VTKFILE,'(a,a,a,I4)') 'SCALARS ', 'Cube', ' float', 1
        write(VTKFILE,'(a,a)') 'LOOKUP_TABLE ', 'DEFAULT'
        do k=1, numzGlobal
            do j=1, numyGlobal
                do i=1, numxGlobal
                    write(VTKFILE,*) 1
                end do
            end do
        end do

        close(VTKFILE)
        deallocate(totDefects)
    end if

    deallocate(numDefects)
    deallocate(numCells)

end subroutine

!*******************************************************************************************************
!> Subroutine find defect in post defect list:
!*******************************************************************************************************
subroutine findDefectInPostList(defectType, postDef, postDefPrev, same)
    use mod_constants
    use mod_structuretype
    implicit none

    integer, intent(in), dimension(SPECIES) :: defectType
    type(postDefect), pointer, intent(inout) :: postDef, postDefPrev
    integer, intent(inout) :: same
    integer :: j

    !nullify(postDefPrev)
    outer: do while(associated(postDef))
        same=0
        inter: do j=1, SPECIES
            if(postDef%defectType(j)==defectType(j)) then
                same=same+1
            else if(iabs(postDef%defectType(j))>iabs(defectType(j))) then
                exit inter
            end if
        end do inter
        if(same==SPECIES) then
            exit outer
        else if(same==j-1 .AND. iabs(postDef%defectType(j))>iabs(defectType(j))) then
            exit outer
        else
            postDefPrev=>postDef
            postDef=>postDef%next
        end if
    end do outer

end subroutine