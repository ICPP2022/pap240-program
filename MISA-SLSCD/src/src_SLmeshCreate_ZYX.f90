!***************************************************************************************************
!>Subroutine SL_initializeMesh():
!The subroutine creates meshes and sectors, and establishes the mapping relationship among meshes, sectors and processes.
!Each process is divided into eight sectors.
!NOTE: Each sector should have at least two meshes per dimension (x, y, z dimension).
!***************************************************************************************************
subroutine SL_initializeMesh()
    use mod_constants
    use mod_globalvariables
    implicit none

    integer :: i, j, k, dir, localElem, dim
    integer :: numxTemp, numyTemp, numzTemp
    integer :: numXmin,numXmax,numYmin,numYmax,numZmin,numZmax	!<used to determine numxLocal, numyLocal, numzLocal
    integer :: subNumXmin,subNumXmax,subNumYmin,subNumYmax,subNumZmin,subNumZmax
    double precision :: lCoordTemp(6)

    totalMeshes = numxGlobal*numyGlobal*numzGlobal
    systemVolume = dble(totalMeshes)*(meshLength**3d0)
    systemCoord=0d0                          !xmin,ymin,zmin
    systemCoord(1)=meshLength*dble(numxGlobal)      !<xmax
    systemCoord(3)=meshLength*dble(numyGlobal)      !<ymax
    systemCoord(5)=meshLength*dble(numzGlobal)      !<zmax

    !*******************************************************
    !Boundaries of processor, temporary
    !*******************************************************
    lCoordTemp(2) = dble(myProc%cartCoords(1))*systemCoord(1)/dble(dims(1))     !<xmin
    lCoordTemp(1) = lCoordTemp(2)+systemCoord(1)/dble(dims(1))                  !<xmax
    lCoordTemp(4) = dble(myProc%cartCoords(2))*systemCoord(3)/dble(dims(2))     !<ymin
    lCoordTemp(3) = lCoordTemp(4)+systemCoord(3)/dble(dims(2))                  !<ymax
    lCoordTemp(6) = dble(myProc%cartCoords(3))*systemCoord(5)/dble(dims(3))     !<zmin
    lCoordTemp(5) = lCoordTemp(6)+systemCoord(5)/dble(dims(3))                  !<zmax

    !*******************************************************
    !Find how many meshes are in local processor and allocate myMesh accordingly,
    !numxLocal, numyLocal and numzLocal are used to determine the size of the mesh inside the local processor.
    !*******************************************************
    numxLocal=0
    numyLocal=0
    numzLocal=0

    !get numxLocal
    if(dims(1) == 1) then  !only 1 proc in x-direction
        numXmin = 0
        numXmax = numxGlobal
        numxLocal = numxGlobal
        myProc%localCoord(1:2) = systemCoord(1:2)
    else if(myProc%cartCoords(1) == 0) then	!at xmin
        numXmin = 0
        numXmax = lCoordTemp(1)/(meshLength/2d0)    !double赋值给integer
        if((dble(numXmax)*(meshLength/2d0))<=lCoordTemp(1) .AND. (dble(numXmax+1)*(meshLength/2d0))>lCoordTemp(1)) then
            if(mod(numXmax,2)==0) then
                numXmax = numXmax/2
            else
                numXmax = (numXmax+1)/2
            end if
        end if
        numxLocal = numXmax
        myProc%localCoord(1) = dble(numXmax)*meshLength
        myProc%localCoord(2) = 0d0
    else if(myProc%cartCoords(1) == (dims(1)-1)) then	!at xmax
        numXmin = lCoordTemp(2)/(meshLength/2d0)
        numXmax = numxGlobal
        if((dble(numXmin)*(meshLength/2d0))<=lCoordTemp(2) .AND. (dble(numXmin+1)*(meshLength/2d0))>lCoordTemp(2)) then
            if(mod(numXmin,2)==0) then
                numXmin = numXmin/2
            else
                numXmin = (numXmin+1)/2
            end if
        end if
        numxLocal = numXmax - numXmin
        myProc%localCoord(1) = systemCoord(1)
        myProc%localCoord(2) = dble(numXmin)*meshLength
    else	!in the middle
        numXmin = lCoordTemp(2)/(meshLength/2d0)
        numXmax = lCoordTemp(1)/(meshLength/2d0)
        if((dble(numXmin)*(meshLength/2d0))<=lCoordTemp(2) .AND. (dble(numXmin+1)*(meshLength/2d0))>lCoordTemp(2)) then
            if(mod(numXmin,2)==0) then
                numXmin = numXmin/2
            else
                numXmin = (numXmin+1)/2
            end if
        end if
        if((dble(numXmax)*(meshLength/2d0))<=lCoordTemp(1) .AND. (dble(numXmax+1)*(meshLength/2d0))>lCoordTemp(1)) then
            if(mod(numXmax,2)==0) then
                numXmax = numXmax/2
            else
                numXmax = (numXmax+1)/2
            end if
        end if
        numxLocal = numXmax-numXmin
        myProc%localCoord(1) = dble(numXmax)*meshLength
        myProc%localCoord(2) = dble(numXmin)*meshLength
    end if

    !get numyLocal
    if(dims(2) == 1) then    !only 1 proc in y-direction
        numYmin = 0
        numYmax = numyGlobal
        numyLocal = numyGlobal
        myProc%localCoord(3:4) = systemCoord(3:4)
    else if(myProc%cartCoords(2) == 0) then	!at ymin
        numYmin = 0
        numYmax = lCoordTemp(3)/(meshLength/2d0)
        if((dble(numYmax)*(meshLength/2d0))<=lCoordTemp(3) .AND. (dble(numYmax+1)*(meshLength/2d0))>lCoordTemp(3)) then
            if(mod(numYmax,2)==0) then
                numYmax = numYmax/2
            else
                numYmax = (numYmax+1)/2
            end if
        end if
        numyLocal = numYmax
        myProc%localCoord(3) = dble(numYmax)*meshLength
        myProc%localCoord(4) = 0d0
    else if(myProc%cartCoords(2) == (dims(2)-1)) then	!at ymax
        numYmin = lCoordTemp(4)/(meshLength/2d0)
        numYmax = numyGlobal
        if((dble(numYmin)*(meshLength/2d0))<=lCoordTemp(4) .AND. (dble(numYmin+1)*(meshLength/2d0))>lCoordTemp(4)) then
            if(mod(numYmin,2)==0) then
                numYmin = numYmin/2
            else
                numYmin = (numYmin+1)/2
            end if
        end if
        numyLocal = numYmax - numYmin
        myProc%localCoord(3) = systemCoord(3)
        myProc%localCoord(4) = dble(numYmin)*meshLength
    else	!in the middle
        numYmin = lCoordTemp(4)/(meshLength/2d0)
        numYmax = lCoordTemp(3)/(meshLength/2d0)
        if((dble(numYmin)*(meshLength/2d0))<=lCoordTemp(4) .AND. (dble(numYmin+1)*(meshLength/2d0))>lCoordTemp(4)) then
            if(mod(numYmin,2)==0) then
                numYmin = numYmin/2
            else
                numYmin = (numYmin+1)/2
            end if
        end if
        if((dble(numYmax)*(meshLength/2d0))<=lCoordTemp(3) .AND. (dble(numYmax+1)*(meshLength/2d0))>lCoordTemp(3)) then
            if(mod(numYmax,2)==0) then
                numYmax = numYmax/2
            else
                numYmax = (numYmax+1)/2
            end if
        end if
        numyLocal = numYmax - numYmin
        myProc%localCoord(3) = dble(numYmax)*meshLength
        myProc%localCoord(4) = dble(numYmin)*meshLength
    end if

    !get numzLocal
    if(dims(3) == 1) then   !only 1 proc in z-direction
        numZmin = 0
        numZmax = numzGlobal
        numzLocal = numzGlobal
        myProc%localCoord(5:6) = systemCoord(5:6)
    else if(myProc%cartCoords(3) == 0) then	!at zmin
        numZmin = 0
        numZmax = lCoordTemp(5)/(meshLength/2d0)
        if((dble(numZmax)*(meshLength/2d0))<=lCoordTemp(5) .AND. (dble(numZmax+1)*(meshLength/2d0))>lCoordTemp(5)) then
            if(mod(numZmax,2)==0) then
                numZmax = numZmax/2
            else
                numZmax = (numZmax+1)/2
            end if
        end if
        numzLocal = numZmax
        myProc%localCoord(5) = dble(numZmax)*meshLength
        myProc%localCoord(6) = 0d0
    else if(myProc%cartCoords(3) == (dims(3)-1)) then	!at zmax
        numZmin = lCoordTemp(6)/(meshLength/2d0)
        numZmax = numzGlobal
        if((dble(numZmin)*(meshLength/2d0))<=lCoordTemp(6) .AND. (dble(numZmin+1)*(meshLength/2d0))>lCoordTemp(6)) then
            if(mod(numZmin,2)==0) then
                numZmin = numZmin/2
            else
                numZmin = (numZmin+1)/2
            end if
        end if
        numzLocal = numZmax - numZmin
        myProc%localCoord(5) = systemCoord(5)
        myProc%localCoord(6) = dble(numZmin)*meshLength
    else	!in the middle
        numZmin = lCoordTemp(6)/(meshLength/2d0)
        numZmax = lCoordTemp(5)/(meshLength/2d0)
        if((dble(numZmin)*(meshLength/2d0))<=lCoordTemp(6) .AND. (dble(numZmin+1)*(meshLength/2d0))>lCoordTemp(6)) then
            if(mod(numZmin,2)==0) then
                numZmin = numZmin/2
            else
                numZmin = (numZmin+1)/2
            end if
        end if
        if((dble(numZmax)*(meshLength/2d0))<=lCoordTemp(5) .AND. (dble(numZmax+1)*(meshLength/2d0))>lCoordTemp(5)) then
            if(mod(numZmax,2)==0) then
                numZmax = numZmax/2
            else
                numZmax = (numZmax+1)/2
            end if
        end if
        numzLocal = numZmax - numZmin
        myProc%localCoord(5) = dble(numZmax)*meshLength
        myProc%localCoord(6) = dble(numZmin)*meshLength
    end if

    !<total meshes of this processor
    numMeshes = numxLocal*numyLocal*numzLocal
    localVolume = dble(numMeshes)*meshLength**3d0
    !allocate(totalRateVol(numMeshes))		!<Create array of total rates in each mesh
    if(numMeshes==0) then
        write(*,*) 'Error: there is no mesh in this processor'
    end if

    !*******************************************************
    !Initialize the meshes in the processor
    !*******************************************************
    allocate(myMesh(numMeshes))
    localElem=0
    do i=1, numxLocal
        do j=1, numyLocal
            do k=1, numzLocal
                localElem=localElem+1

                numxTemp = numXmin+i
                numyTemp = numYmin+j
                numzTemp = numZmin+k

                myMesh(localElem)%gid=(numxTemp-1)*numzGlobal*numyGlobal+(numyTemp-1)*numzGlobal+numzTemp
            end do
        end do
    end do

    !*******************************************************
    !<Assign neighbors and processor numbers for neighbors (periodic or free surfaces in x/y/z)
    !*******************************************************
    call createMeshConnect()
    !*******************************************************
    !<If SLKMC is used, this subroutine is needed (create connectivity for sectors).
    !*******************************************************
    call createSectorConnect()

end subroutine SL_initializeMesh

!**************************************************************************************************
!>Subroutine createConnectLocal(): create LOCAL connectivity.
!It identifies the mesh and processor of neighboring meshes for each mesh.
!The connectivity scheme is the same as in the global case, but neighboring processor numbers are used here.
!**************************************************************************************************
subroutine createMeshConnect()
    use mod_globalvariables
    implicit none
    include 'mpif.h'

    integer :: cell, localCell, globalCell, neighborGID
    integer :: i, dir, recvDir, count
    integer :: x, y, z, temp, index, max
    integer, external :: findgNeighborGID
    !<Used for communication
    integer :: numSend(6), numRecv(6)
    integer, allocatable :: sendBuff(:,:), recvBuff(:)
    integer :: status(MPI_STATUS_SIZE)

    max=0
    if(numzLocal*numyLocal>=max) then
        max=numzLocal*numyLocal
    end if
    if(numyLocal*numxLocal>=max) then
        max=numyLocal*numxLocal
    end if
    if(numzLocal*numxLocal>=max) then
        max=numzLocal*numxLocal
    end if
    allocate(sendBuff(max,6))
    sendBuff = 0
    numSend = 0
    numRecv = 0


    do cell=1,numMeshes
        !<Right (+x)
        if(mod(cell,numMeshes) > numzLocal*numyLocal*(numxLocal-1) .OR. mod(cell,numMeshes)==0) then  !<Boundary of +x
            myMesh(cell)%neighborProc(1) = myProc%neighborProc(1)
            if(myProc%neighborProc(1) == -1) then     !+x free surface
                myMesh(cell)%neighbor(1) = 0
            else if(myProc%neighborProc(1) == myProc%taskid) then     !one processor in x-direction
                myMesh(cell)%neighbor(1) = cell-numzLocal*numyLocal*(numxLocal-1)
            else
                numSend(1) = numSend(1) + 1
                sendBuff(numSend(1),1) = cell
            end if
        else    !+x neighbor is in current proc
            myMesh(cell)%neighbor(1) = cell+numzLocal*numyLocal
            myMesh(cell)%neighborProc(1) = myProc%taskid
        end if

        !<Left (-x)
        if(mod(cell,numMeshes)<=numzLocal*numyLocal .AND. (mod(cell,numMeshes)/=0 .OR. numxLocal==1)) then    !<Boundary of -x
            myMesh(cell)%neighborProc(2) = myProc%neighborProc(2)
            if(myProc%neighborProc(2) == -1) then     !-x free surface
                myMesh(cell)%neighbor(2) = 0
            else if(myProc%neighborProc(2) == myProc%taskid) then     !one processor in x-direction
                myMesh(cell)%neighbor(2) = cell+(numzLocal*numyLocal*(numxLocal-1))
            else
                numSend(2) = numSend(2) + 1
                sendBuff(numSend(2),2) = cell
            end if
        else
            myMesh(cell)%neighbor(2) = cell-numzLocal*numyLocal
            myMesh(cell)%neighborProc(2) = myProc%taskid
        end if

        !<Front (+y)
        if(mod(cell,numzLocal*numyLocal)>numzLocal*(numyLocal-1) .OR. mod(cell,numzLocal*numyLocal)==0) then  !<Boundary of +y
            myMesh(cell)%neighborProc(3) = myProc%neighborProc(3)
            if(myProc%neighborProc(3) == -1) then     !+y free surface
                myMesh(cell)%neighbor(3) = 0
            else if(myProc%neighborProc(3) == myProc%taskid) then     !one processor in y-direction
                myMesh(cell)%neighbor(3) = cell-numzLocal*(numyLocal-1)
            else
                numSend(3) = numSend(3) + 1
                sendBuff(numSend(3),3) = cell
            end if
        else
            myMesh(cell)%neighbor(3) = cell+numzLocal
            myMesh(cell)%neighborProc(3) = myProc%taskid
        end if

        !<Back (-y)
        if(mod(cell,numzLocal*numyLocal)<=numzLocal.AND.(mod(cell,numzLocal*numyLocal)/=0 .OR. numyLocal==1)) then   !<Boundary of -y
            myMesh(cell)%neighborProc(4) = myProc%neighborProc(4)
            if(myProc%neighborProc(4) == -1) then     !-y free surface
                myMesh(cell)%neighbor(4) = 0
            else if(myProc%neighborProc(4) == myProc%taskid) then     !one processor in y-direction
                myMesh(cell)%neighbor(4)=cell+numzLocal*(numyLocal-1)
            else
                numSend(4) = numSend(4) + 1
                sendBuff(numSend(4),4) = cell
            end if
        else
            myMesh(cell)%neighbor(4) = cell-numzLocal
            myMesh(cell)%neighborProc(4) = myProc%taskid
        end if

        !<Up (+z)
        if(mod(cell, numzLocal)==0) then    !<boundary of +z
            myMesh(cell)%neighborProc(5) = myProc%neighborProc(5)
            if(myProc%neighborProc(5) == -1) then     !+z free surface
                myMesh(cell)%neighbor(5) = 0
            else if(myProc%neighborProc(5) == myProc%taskid) then     !one processor in z-direction
                myMesh(cell)%neighbor(5) = cell-numzLocal+1
            else
                numSend(5) = numSend(5) + 1
                sendBuff(numSend(5),5) = cell
            end if
        else
            myMesh(cell)%neighbor(5) = cell + 1
            myMesh(cell)%neighborProc(5) = myProc%taskid
        end if

        !<Down (-z)
        if(mod(cell+numzLocal-1,numzLocal)==0) then !<boundary of -z
            myMesh(cell)%neighborProc(6) = myProc%neighborProc(6)
            if(myProc%neighborProc(6) == -1) then     !-z free surface
                myMesh(cell)%neighbor(6) = 0
            else if(myProc%neighborProc(6) == myProc%taskid) then     !one processor in z-direction
                myMesh(cell)%neighbor(6) = cell+numzLocal-1
            else
                numSend(6) = numSend(6) + 1
                sendBuff(numSend(6),6) = cell
            end if
        else
            myMesh(cell)%neighbor(6) = cell-1
            myMesh(cell)%neighborProc(6) = myProc%taskid
        end if
    end do

    do dir=1,6
        if(myProc%neighborProc(dir)/=myProc%taskid .AND. myProc%neighborProc(dir)/=-1) then
            if(dir==1 .OR. dir==2) then         !<x-direction
                myGhost(dir)%numCells = numyLocal*numzLocal
            else if(dir==3 .OR. dir==4) then    !<y-direction
                myGhost(dir)%numCells = numxLocal*numzLocal
            else if(dir==5 .OR. dir==6) then    !<z-direction
                myGhost(dir)%numCells = numxLocal*numyLocal
            end if
        else
            myGhost(dir)%numCells = 0
        end if
        allocate(myGhost(dir)%cell(myGhost(dir)%numCells))  !ghost cell id
        allocate(myGhost(dir)%local(myGhost(dir)%numCells)) !local cell id
    end do

    !Send / Recv
    do dir=1,6
        if(mod(dir,2)==0) then
            recvDir = dir-1
        else
            recvDir = dir+1
        end if

        if(myProc%neighborProc(dir)/=myProc%taskid .AND. myProc%neighborProc(dir)/=-1) then
            call MPI_SEND(sendBuff(1,dir),numSend(dir),MPI_INTEGER,myProc%neighborProc(dir),99,comm,ierr)
        end if

        if(myProc%neighborProc(recvDir)/=myProc%taskid .AND. myProc%neighborProc(recvDir)/=-1) then

            if(recvDir==1 .OR. recvDir==2) then !x-dir
                numRecv(recvDir) = numyLocal*numzLocal
            else if(recvDir==3 .OR. recvDir==4) then    !y-dir
                numRecv(recvDir) = numxLocal*numzLocal
            else if(recvDir==5 .OR. recvDir==6) then    !z-dir
                numRecv(recvDir) = numxLocal*numyLocal
            end if
            allocate(recvBuff(numRecv(recvDir)))
            call MPI_RECV(recvBuff,numRecv(recvDir),MPI_INTEGER,myProc%neighborProc(recvDir),99,comm,status,ierr)


            do i=1, numRecv(recvDir)
                localCell = sendBuff(i,recvDir)
                myMesh(localCell)%neighbor(recvDir) = recvBuff(i)
                myGhost(recvDir)%cell(i) = myMesh(localCell)%neighbor(recvDir)
                myGhost(recvDir)%local(i) = localCell
            end do
            deallocate(recvBuff)
        end if
    end do

    deallocate(sendBuff)

end subroutine createMeshConnect

!***************************************************************************************
!>Subroutine createSectorConnect(): create connectivity for sectors
!If SLKMC is used, this subroutine is needed.
!***************************************************************************************
subroutine createSectorConnect()
    use mod_constants
    use mod_globalvariables
    implicit none

    integer :: i

    do i=1,8        !z->y->x
        !<Right (+x)
        if(i > 4) then              !<Boundary of +x, sector 5,6,7,8
            if(myProc%neighborProc(1) == -1) then     !+x free surface
                mySector(i)%commDir(1,1) = 0  !(sendDir, dim1), recvDir is opposite to the sendDir
            else if(myProc%neighborProc(1) == myProc%taskid) then
                mySector(i)%commDir(1,1) = 0  !(sendDir, dim1)
            else
                mySector(i)%commDir(1,1)=1  !(sendDir, dim1)
            end if
            if(myProc%neighborProc(2)/=myProc%taskid .AND. myProc%neighborProc(2)/=-1) then
                mySector(i)%commDir(2,1) = 2  !(recvDir, dim1)
            else
                mySector(i)%commDir(2,1) = 0  !(recvDir, dim1)
            end if
            mySector(i)%numCell(1) = numxLocal/2
        else    !sector 1,2,3,4
            if(mod(numxLocal,2) == 0) then
                mySector(i)%numCell(1) = numxLocal/2
            else
                mySector(i)%numCell(1) = numxLocal/2+1
            end if
        end if

        !<Left (-x)
        if(i <= 4) then     !<Boundary of -x, sector 1,2,3,4
            if(myProc%neighborProc(2) == -1) then
                mySector(i)%commDir(1,1) = 0  !(sendDir, dim1)
            else if(myProc%neighborProc(2) == myProc%taskid) then
                mySector(i)%commDir(1,1) = 0  !(sendDir, dim1)
            else
                mySector(i)%commDir(1,1) = 2  !(sendDir, dim1)
            end if
            if(myProc%neighborProc(1)/=myProc%taskid .AND. myProc%neighborProc(1)/=-1) then
                mySector(i)%commDir(2,1) = 1  !(recvDir, dim1)
            else
                mySector(i)%commDir(2,1) = 0  !(recvDir, dim1)
            end if
        end if

        !<Front (+y)
        if(mod(i,4)>2 .OR. mod(i,4)==0) then                 !<Boundary of +y, sector 3,4,7,8
            if(myProc%neighborProc(3) == -1) then
                mySector(i)%commDir(1,2) = 0  !(secdDir, dim2)
            else if(myProc%neighborProc(3) == myProc%taskid) then
                mySector(i)%commDir(1,2) = 0  !(sendDir, dim2)
            else
                mySector(i)%commDir(1,2) = 3  !(sendDir, dim2)
            end if
            if(myProc%neighborProc(4)/=myProc%taskid .AND. myProc%neighborProc(4)/=-1) then
                mySector(i)%commDir(2,2) = 4  !(recvDir, dim2)
            else
                mySector(i)%commDir(2,2) = 0  !(recvDir, dim2)
            end if
            mySector(i)%numCell(2) = numyLocal/2
        else    !sector 1,2,5,6
            if(mod(numyLocal,2) == 0) then
                mySector(i)%numCell(2) = numyLocal/2
            else
                mySector(i)%numCell(2) = numyLocal/2+1
            end if
        end if

        !<Back (-y)
        if(mod(i,4)<=2 .AND. mod(i,4)/=0) then                !<Boundary of -y, sector 1,2,5,6
            if(myProc%neighborProc(4) == -1) then
                mySector(i)%commDir(1,2) = 0  !(secdDir, dim2)
            else if(myProc%neighborProc(4) == myProc%taskid) then
                mySector(i)%commDir(1,2) = 0  !(secdDir, dim2)
            else
                mySector(i)%commDir(1,2) = 4  !(secdDir, dim2)
            end if
            if(myProc%neighborProc(3)/=myProc%taskid .AND. myProc%neighborProc(3)/=-1) then
                mySector(i)%commDir(2,2) = 3  !(recvDir, dim2)
            else
                mySector(i)%commDir(2,2) = 0  !(recvDir, dim2)
            end if
        end if

        !<Up (+z)
        if(mod(i,2)==0) then                                   !<Boundary of +z, sector 2,4,6,8
            if(myProc%neighborProc(5) == -1) then
                mySector(i)%commDir(1,3) = 0  !(secdDir, dim3)
            else if(myProc%neighborProc(5) == myProc%taskid) then
                mySector(i)%commDir(1,3) = 0  !(secdDir, dim3)
            else
                mySector(i)%commDir(1,3) = 5  !(secdDir, dim3)
            end if
            if(myProc%neighborProc(6)/=myProc%taskid .AND. myProc%neighborProc(6)/=-1) then
                mySector(i)%commDir(2,3) = 6  !(recvDir, dim3)
            else
                mySector(i)%commDir(2,3) = 0  !(recvDir, dim3)
            end if
            mySector(i)%numCell(3) = numzLocal/2
        else    !sector 1,3,5,7
            if(mod(numzLocal,2) == 0) then
                mySector(i)%numCell(3) = numzLocal/2
            else
                mySector(i)%numCell(3) = numzLocal/2+1
            end if
        end if

        !<Down (-z)
        if(mod(i+1,2)==0) then                                  !<Boundary of -z, sector 1,3,5,7
            if(myProc%neighborProc(6) == -1) then
                mySector(i)%commDir(1,3) = 0  !(secdDir, dim3)
            else if(myProc%neighborProc(6) == myProc%taskid) then
                mySector(i)%commDir(1,3) = 0  !(secdDir, dim3)
            else
                mySector(i)%commDir(1,3) = 6  !(secdDir, dim3)
            end if
            if(myProc%neighborProc(5)/=myProc%taskid .AND. myProc%neighborProc(5)/=-1) then
                mySector(i)%commDir(2,3) = 5  !(recvDir, dim3)
            else
                mySector(i)%commDir(2,3) = 0  !(recvDir, dim3)
            end if
        end if

        !mySector(i)%numCells = mySector(i)%numCell(1) * mySector(i)%numCell(2) * mySector(i)%numCell(3)
        !if(mySector(i)%numCells == 0) then
        !    write(*,*) 'Error: there is no mesh in sector', i, 'of proc', myProc%taskid
        !end if
    end do

end subroutine createSectorConnect

!***************************************************************************************
!>Function findgNeighborGID(globalID, dir)
!Inputs: globalID, direction
!Output: globalNeighborGID of the globalID in the direaction
!***************************************************************************************
integer function findNeighborGID(globalID, dir)
    use mod_globalvariables
    implicit none

    integer, intent(in) :: globalID, dir
    integer :: neighborGID

    if(dir==1) then         !<+x
        if(mod(globalID,numzGlobal*numyGlobal*numxGlobal)>numzGlobal*numyGlobal*(numxGlobal-1).OR.&
                mod(globalID,numzGlobal*numyGlobal*numxGlobal)==0) then             !<boundary of +x
            if(meshType(1)==0) then       !Free surface in x-direction
                neighborGID=0
            else
                neighborGID=globalID-(numzGlobal*numyGlobal*(numxGlobal-1))
            end if
        else
            neighborGID=globalID+numzGlobal*numyGlobal
        end if
    else if(dir==2) then    !<-x
        if(mod(globalID,numzGlobal*numyGlobal*numxGlobal)<=numzGlobal*numyGlobal.AND.&
                (mod(globalID,numzGlobal*numyGlobal*numxGlobal)/= 0.OR.numxGlobal==1)) then !<boundary of -x
            if(meshType(1)==0) then       !Free surface in x-direction
                neighborGID=0
            else
                neighborGID=globalID+(numzGlobal*numyGlobal*(numxGlobal-1))
            end if
        else
            neighborGID=globalID-numzGlobal*numyGlobal
        end if
    else if(dir==3) then    !<+y
        if((mod(globalID,numzGlobal*numyGlobal) > numzGlobal*(numyGlobal-1)).OR.&
                (mod(globalID,numzGlobal*numyGlobal)==0)) then                      !<boundary of +y
            if(meshType(2)==0) then     !Free surface in y-direction
                neighborGID=0
            else
                neighborGID=globalID-(numzGlobal*(numyGlobal-1))
            end if
        else
            neighborGID=globalID+numzGlobal
        end if
    else if(dir==4) then    !<-y
        if((mod(globalID,numzGlobal*numyGlobal) <= numzGlobal).AND.&
                (mod(globalID,numzGlobal*numyGlobal)/= 0.OR. numyGlobal==1)) then    !<boundary of -y
            if(meshType(2)==0) then     !Free surface in y-direction
                neighborGID=0
            else
                neighborGID=globalID+(numzGlobal*(numyGlobal-1))
            end if
        else
            neighborGID=globalID-numzGlobal
        end if
    else if(dir==5) then    !<+z
        if(mod(globalID,numzGlobal)==0) then                                        !<boundary of +z
            if(meshType(3)==0) then !Free surface in z-direction
                neighborGID=0
            else
                neighborGID=globalID-numzGlobal+1
            end if
        else
            neighborGID=globalID+1
        end if
    else if(dir==6) then    !<-z
        if(mod(globalID+numzGlobal-1,numzGlobal)==0) then                           !<boundary of -z
            if(meshType(3)==0) then !Free surface in z-direction
                neighborGID=0
            else
                neighborGID=globalID+numzGlobal-1
            end if
        else
            neighborGID=globalID-1
        end if
    end if
    findNeighborGID=neighborGID

end function findNeighborGID

!***************************************************************************************
!>Function getLocalCell(subCell, sector)
!Inputs: the number of a cell in a sector and the number of that sector
!Output: local number of a cell in the processor
!***************************************************************************************
integer function getLocalCell(subCell, sector)
    use mod_globalvariables
    implicit none

    integer, intent(in) :: subCell, sector
    integer :: xIndex, yIndex, zIndex, tempID

    tempID=subCell-1
    zIndex=mod(tempID,mySector(sector)%numCell(3))+1
    tempID=tempID/mySector(sector)%numCell(3)
    yIndex=mod(tempID, mySector(sector)%numCell(2))+1
    tempID=tempID/mySector(sector)%numCell(2)
    xIndex=mod(tempID, mySector(sector)%numCell(1))+1

    if(sector>4) then    !<Boundary of +x
        xIndex=xIndex+(numxLocal-mySector(sector)%numCell(1))
    end if

    if(mod(sector,4)>2 .OR. mod(sector,4)==0) then      !<Boundary of +y
        yIndex=yIndex+(numyLocal-mySector(sector)%numCell(2))
    end if

    if(mod(sector,2)==0) then                             !<Boundary of +z
        zIndex=zIndex+(numzLocal-mySector(sector)%numCell(3))
    end if
    getLocalCell=(xIndex-1)*numzLocal*numyLocal+(yIndex-1)*numzLocal+zIndex
end function getLocalCell