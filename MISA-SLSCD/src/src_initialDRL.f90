!*****************************************************************************************
!>Subroutine initializeDRL_Defect(): initializes defects in DRL for each mesh
!The first node of each level in DRL is initialized to 0.
!*****************************************************************************************
subroutine initializeDRL_Defect()
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none

    integer :: cell, i, level, numI, numV
    type(defect), pointer :: def, def_m, def2, def_m2
    double precision, external :: compute_diffusivity             !<Function declaration
    integer, external :: find_level

    allocate(myDRL(LEVELS,numMeshes))      !contain mobile and immobile defects

    do cell=1,numMeshes
        do level=1, LEVELS
            nullify(myDRL(level,cell)%defect_all)       !Initializes the header node
            nullify(myDRL(level,cell)%defect_mobile)
        end do

        numI=0
        do i=1, initNumSIA
            if(initPointDefects(i,1)==myMesh(cell)%gid) then
                numI=numI+1
            end if
        end do
        numV=0
        do i=1, initNumVac
            if(initPointDefects(i,2)==myMesh(cell)%gid) then
                numV=numV+1
            end if
        end do

        !<SIA_0_0, initialize the first node
        if(numI>0) then
            nullify(def)
            allocate(def)
            allocate(def%defectType(SPECIES))
            def%defectType(1:SPECIES)=(/1,0,0/)
            def%num=numI
            def%diff=compute_diffusivity(def%defectType)
            nullify(def%next)
            nullify(def%reactionList)
            level=find_level(def%defectType)
            myDRL(level,cell)%defect_all=>def

            nullify(def_m)
            allocate(def_m)
            allocate(def_m%defectType(SPECIES))
            def_m%defectType=def%defectType
            def_m%num=def%num
            def_m%diff=def%diff
            nullify(def_m%next)
            nullify(def_m%reactionList)
            myDRL(level,cell)%defect_mobile=>def_m
        end if

        !<V_0_0, initialize the first node
        if(numV>0) then
            nullify(def)
            allocate(def)
            allocate(def%defectType(SPECIES))
            def%defectType(1:SPECIES)=(/-1,0,0/)
            def%num=numV
            def%diff=compute_diffusivity(def%defectType)
            nullify(def%next)
            nullify(def%reactionList)
            level=find_level(def%defectType)
            myDRL(level,cell)%defect_all=>def

            nullify(def_m)
            allocate(def_m)
            allocate(def_m%defectType(SPECIES))
            def_m%defectType=def%defectType
            def_m%num=def%num
            def_m%diff=def%diff
            nullify(def_m%next)
            nullify(def_m%reactionList)
            myDRL(level,cell)%defect_mobile=>def_m
        end if

        !<0_S_0 clusters, initialize the first node
        if(soluteConc > 0d0) then
            nullify(def)
            allocate(def)
            allocate(def%defectType(SPECIES))
            def%defectType(1:SPECIES)=(/0,1,0/)
            def%num=numSevermesh
            def%diff=compute_diffusivity(def%defectType)
            nullify(def%next)
            nullify(def%reactionList)
            level=find_level(def%defectType)
            myDRL(level,cell)%defect_all=>def

            nullify(def_m)
            allocate(def_m)
            allocate(def_m%defectType(SPECIES))
            def_m%defectType = def%defectType
            def_m%num = def%num
            def_m%diff = def%diff
            nullify(def_m%next)
            nullify(def_m%reactionList)
            myDRL(level,cell)%defect_mobile=>def_m
        end if
    end do

end subroutine initializeDRL_Defect

!*****************************************************************************************
!>Subroutine initializeDRL_Ghost():
!Initializes defects for each mesh in ghost.
!Begins with a defect with type 0 0 0 and num 0.
!*****************************************************************************************
subroutine initializeDRL_Ghost()
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none

    integer :: dir, index, i, globalCell, neighborGID, level, numI, numV
    type(defect), pointer :: def, def_m
    double precision, external :: compute_diffusivity             !<Function declaration
    integer, external :: findNeighborGID
    integer, external :: find_level
    !test
    type(defect), pointer :: def2, def_m2

    do dir=1,6
        if(myProc%neighborProc(dir)/=myProc%taskid .AND. myProc%neighborProc(dir)/=-1) then
            allocate(myGhost(dir)%myDRL_ghost(LEVELS,myGhost(dir)%numCells))

            do index=1, myGhost(dir)%numCells
                do level=1, LEVELS
                    nullify(myGhost(dir)%myDRL_ghost(level,index)%defect_mobile)
                end do

                globalCell = myMesh(myGhost(dir)%local(index))%gid
                neighborGID = findNeighborGID(globalCell, dir)
                numI = 0
                do i=1, initNumSIA
                    if(initPointDefects(i,1) == neighborGID) then
                        numI = numI + 1
                    end if
                end do
                numV = 0
                do i=1, initNumVac
                    if(initPointDefects(i,2) == neighborGID) then
                        numV = numV + 1
                    end if
                end do

                !<SIA_0_0
                if(numI > 0) then
                    nullify(def)
                    allocate(def)
                    allocate(def%defectType(SPECIES))
                    def%defectType = (/1,0,0/)
                    def%num = numI
                    def%diff = compute_diffusivity(def%defectType)
                    nullify(def%next)
                    nullify(def%reactionList)
                    level = find_level(def%defectType)
                    myGhost(dir)%myDRL_ghost(level,index)%defect_mobile => def
                end if

                !<V_0_0
                if(numV > 0) then
                    nullify(def)
                    allocate(def)
                    allocate(def%defectType(SPECIES))
                    def%defectType = (/-1,0,0/)
                    def%num = numV
                    def%diff = compute_diffusivity(def%defectType)
                    nullify(def%next)
                    nullify(def%reactionList)
                    level = find_level(def%defectType)
                    myGhost(dir)%myDRL_ghost(level,index)%defect_mobile => def
                end if

                !<0_S_0
                if(soluteConc > 0d0) then
                    nullify(def)
                    allocate(def)
                    allocate(def%defectType(SPECIES))
                    def%defectType = (/0,1,0/)
                    def%num = numSevermesh
                    def%diff = compute_diffusivity(def%defectType)
                    nullify(def%next)
                    nullify(def%reactionList)
                    level = find_level(def%defectType)
                    myGhost(dir)%myDRL_ghost(level,index)%defect_mobile => def
                end if
            end do
        end if
    end do

end subroutine initializeDRL_Ghost

!*****************************************************************************************
!>Subroutine initializeDRL_Reaction(): creates a new reaction list for each mesh and
!initializes reactions based on the type of defects that already exits.
!The first reaction is implantation reaction or 0.
!*****************************************************************************************
subroutine initializeDRL_Reaction()
    use mod_constants
    use mod_structuretype
    use mod_globalVariables
    use mod_update_reactions
    implicit none

    integer :: cell, level, level2, level_m, dir
    type(defect), pointer :: def, def_m, def2
    type(reaction), pointer :: reac
    double precision :: diff, diff_m, rate

    allocate(implantation(numMeshes))
    totalRateVol = 0d0
    !totalRate = 0d0

    !<First: initialize 0th reaction
    do cell=1, numMeshes
        if(dpaRate > 0d0 .AND. totalDPA > 0d0) then !FrenkelPair or Cascade
            if(irradiationType=='FrenkelPair') then
                implantation(cell)%numReactants = 0
                implantation(cell)%numProducts = 2
            else if(irradiationType=='Cascade') then
                implantation(cell)%numReactants = -10
                implantation(cell)%numProducts = 0
            end if
            implantation(cell)%rate = compute_0nd_rate(implantation(cell)%numReactants,implantation(cell)%numProducts,cell)
            nullify(implantation(cell)%next)
        else
            implantation(cell)%numReactants = 0
            implantation(cell)%numProducts = 0
            implantation(cell)%rate = 0d0
            nullify(implantation(cell)%next)
        end if
        totalRateVol(cell) = totalRateVol(cell) + implantation(cell)%rate

        if(HeDPAratio > 0d0) then
            nullify(reac)
            allocate(reac)
            reac%numReactants = 0
            reac%numProducts = 1
            reac%rate = compute_0nd_rate(reac%numReactants, reac%numProducts, cell)
            nullify(reac%next)
            implantation(cell)%next => reac
            totalRateVol(cell) = totalRateVol(cell) + reac%rate
        end if
    end do

    !<Second: initialize 1st, diff, 2nd reactions for initial defects
    do cell=1, numMeshes
        do level=1, LEVELS
            nullify(def)
            def => myDRL(level,cell)%defect_all  !pont to the first defect node
            do while(associated(def))
                call update_1st_Reaction(def, cell)
                if(def%diff > 0d0) then
                    !call update_diff_reaction_cellTOsix(def, cell)
                    do dir=1, 6
                        call update_diff_reaction_cellTOcell(def,cell,myMesh(cell)%neighbor(dir),dir,&
                                myProc%taskid, myMesh(cell)%neighborProc(dir))
                    end do
                end if

                if(def%diff <= 0d0) then
                    do level_m=1,LEVELS
                        nullify(def_m)
                        def_m => myDRL(level_m,cell)%defect_mobile
                        do while(associated(def_m))
                            !call update_2nd_Reaction(def, def_m, cell)
                            call update_2nd_reactionCompare(def, def_m, cell)
                            def_m=>def_m%next
                        end do
                    end do
                else
                    do level2=1, LEVELS
                        nullify(def2)
                        def2 => myDRL(level2,cell)%defect_all
                        do while(associated(def2))
                            !call update_2nd_Reaction(def2, def, cell)
                            call update_2nd_reactionCompare(def2, def, cell)
                            def2=>def2%next
                        end do
                    end do
                end if

                def=>def%next
            end do
        end do
    end do

end subroutine