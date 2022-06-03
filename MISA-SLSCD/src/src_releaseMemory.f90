!***********************************************************************
!> Subroutine deallocateDefectLists: deallocate defects in the local meshes
!***********************************************************************
subroutine deallocateMyDRL()
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none

    type(defect), pointer :: def, def_m
    type(reaction), pointer :: reac
    integer :: cell, level

    do cell=1, numMeshes
        !deallocate implantation
        nullify(reac)
        reac => implantation(cell)%next
        do while(associated(reac))
            implantation(cell)%next => reac%next
            deallocate(reac)
            reac => implantation(cell)%next
        end do

        !deallocate myDRL
        do level=1, LEVELS
            !deallocate %defect_all
            nullify(def)
            def => myDRL(level,cell)%defect_all
            do while(associated(def))
                nullify(reac)
                reac => def%reactionList
                do while(associated(reac))
                    def%reactionList => reac%next
                    if(allocated(reac%reactants)) then
                        deallocate(reac%reactants)
                    end if
                    if(allocated(reac%products)) then
                        deallocate(reac%products)
                    end if
                    if(allocated(reac%cell)) then
                        deallocate(reac%cell)
                    end if
                    if(allocated(reac%taskid)) then
                        deallocate(reac%taskid)
                    end if
                    deallocate(reac)
                    reac => def%reactionList
                end do
                myDRL(level,cell)%defect_all => def%next
                deallocate(def%defectType)
                deallocate(def)
                def => myDRL(level,cell)%defect_all
            end do

            !deallocate %defect_mobile
            nullify(def_m)
            def_m => myDRL(level,cell)%defect_mobile
            do while(associated(def_m))
                myDRL(level,cell)%defect_mobile => def_m%next
                deallocate(def_m%defectType)
                deallocate(def_m)
                def_m => myDRL(level,cell)%defect_mobile
            end do
        end do
    end do
    deallocate(myDRL)
    deallocate(implantation)
    !<clear totalRate(numMeshes)
    totalRateVol = 0d0

    if(allocated(initPointDefects)) then
        deallocate(initPointDefects)
    end if

end subroutine deallocateMyDRL

!***********************************************************************
!> Subroutine deallocateDefectLists: deallocate defects in the boundary meshes
!***********************************************************************
subroutine deallocateMyGhostDRL()
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none

    type(defect), pointer :: def, defectPrev, def_m
    integer :: cell, level, dir, index

    do dir=1,6
        do index=1, myGhost(dir)%numCells
            !deallocate %defect_mobile
            do level=1, LEVELS
                nullify(def)
                def => myGhost(dir)%myDRL_ghost(level,index)%defect_mobile
                do while(associated(def))
                    myGhost(dir)%myDRL_ghost(level,index)%defect_mobile => def%next
                    deallocate(def%defectType)
                    deallocate(def)
                    def => myGhost(dir)%myDRL_ghost(level,index)%defect_mobile
                end do
            end do
        end do

        if(myProc%neighborProc(dir)/=myProc%taskid .AND. myProc%neighborProc(dir)/=-1) then
            deallocate(myGhost(dir)%myDRL_ghost)
        end if
    end do

end subroutine deallocateMyGhostDRL

!***********************************************************************
!> Subroutine deallocateInputs: deallocate input data (cascades, binding and migration energies, etc)
! Including:
! CascadeList
! formSingle(:), diffSingle(:), diffFunc(:), bindSingle(:), bindFunc(:),
! implantReactions(:), dissocReactions(:), sinkReactions(:), diffReactions(:), clusterReactions(:)
!***********************************************************************
subroutine deallocateInputs()
    use mod_structuretype
    use mod_globalvariables
    implicit none

    type(cascade), pointer :: cas
    type(cascadeDefect), pointer :: def
    integer :: i, fileID, test

    !<Deallocate cascadeList
    if(irradiationType=='Cascade') then
        if(numCascadeFiles > 1) then    !<delete cascadeLists(numCascadeFiles)
            do fileID=1, numCascadeFiles
                !<delete %listDefects
                test=1
                nullify(cas)
                cas=>cascadeLists(fileID)%listCascades%next
                do while(associated(cas))
                    test=test+1
                    nullify(def)
                    def=>cas%listDefects%next
                    do while(associated(def))
                        cas%listDefects%next=>def%next
                        deallocate(def%defectType)
                        deallocate(def)
                        def=>cas%listDefects%next
                    end do
                    deallocate(cas%listDefects%defectType)
                    deallocate(cas%listDefects)

                    cascadeLists(fileID)%listCascades%next=>cas%next
                    deallocate(cas)
                    cas=>cascadeLists(fileID)%listCascades%next
                    !cascadeLists(fileID)%listCascades=>cas%next
                    !cas=>cascadeLists(fileID)%listCascades
                end do
                !<delete the first cascade
                nullify(cas)
                cas=>cascadeLists(fileID)%listCascades
                do while(associated(cas))
                    nullify(def)
                    def=>cas%listDefects%next
                    do while(associated(def))
                        cas%listDefects%next=>def%next
                        deallocate(def%defectType)
                        deallocate(def)
                        def=>cas%listDefects%next
                    end do
                    deallocate(cas%listDefects%defectType)
                    deallocate(cas%listDefects)

                    deallocate(cas)
                    nullify(cas)
                end do
            end do
            deallocate(cascadeLists)
        else        !<delete cascadeList
            nullify(cas)
            cas=>cascadeList%next
            !cas=>cascadeList
            do while(associated(cas))
                nullify(def)
                def=>cas%listDefects%next
                !def=>cas%listDefects
                do while(associated(def))
                    cas%listDefects%next=>def%next
                    !cas%listDefects=>def%next
                    deallocate(def%defectType)
                    deallocate(def)
                    def=>cas%listDefects%next
                    !def=>cas%listDefects
                end do
                deallocate(cas%listDefects%defectType)
                deallocate(cas%listDefects)

                cascadeList%next=>cas%next
                !cascadeList=>cas%next
                deallocate(cas)
                cas=>cascadeList%next
                !cas=>cascadeList
            end do
            !<delete the first cascade
            nullify(cas)
            cas=>cascadeList
            do while(associated(cas))
                nullify(def)
                def=>cas%listDefects%next
                !def=>cas%listDefects
                do while(associated(def))
                    cas%listDefects%next=>def%next
                    !cas%listDefects=>def%next
                    deallocate(def%defectType)
                    deallocate(def)
                    def=>cas%listDefects%next
                    !def=>cas%listDefects
                end do
                deallocate(cas%listDefects%defectType)
                deallocate(cas%listDefects)

                deallocate(cas)
                nullify(cas)
            end do
        end if
    end if

    if(PKAspectrum=='yes') then
        deallocate(EPKAlist%energy)
        deallocate(EPKAlist%cpdf)
    end if

    !<Deallocate defect attributes read in from input file
    if(allocated(formSingle)) then
        do i=1, numFormSingle
            deallocate(formSingle(i)%defectType)
        end do
        deallocate(formSingle)
    end if

    if(allocated(diffSingle)) then
        do i=1,numDiffSingle
            deallocate(diffSingle(i)%defectType)
        end do
        deallocate(diffSingle)
    end if

    if(allocated(diffFunc)) then
        do i=1,numDiffFunc
            deallocate(diffFunc(i)%defectType)
            deallocate(diffFunc(i)%min)
            deallocate(diffFunc(i)%max)
            if(allocated(diffFunc(i)%para)) then
                deallocate(diffFunc(i)%para)
            end if
        end do
        deallocate(diffFunc)
    end if

    if(allocated(bindSingle)) then
        do i=1,numBindSingle
            deallocate(bindSingle(i)%defectType)
            deallocate(bindSingle(i)%product)
        end do
        deallocate(bindSingle)
    end if

    if(allocated(bindFunc)) then
        do i=1,numBindFunc
            deallocate(bindFunc(i)%defectType)
            deallocate(bindFunc(i)%product)
            deallocate(bindFunc(i)%min)
            deallocate(bindFunc(i)%max)
            if(allocated(bindFunc(i)%para)) then
                deallocate(bindFunc(i)%para)
            end if
        end do
        deallocate(bindFunc)
    end if

    !<Deallocate reaction parameters
    !if(allocated(implantReactions)) then
    !    do i=1,numImplantReaction
    !        if(allocated(implantReactions(i)%reactants)) then
    !            deallocate(implantReactions(i)%reactants)
    !        end if
    !        if(allocated(implantReactions(i)%products)) then
    !            deallocate(implantReactions(i)%products)
    !        end if
    !    end do
    !    deallocate(implantReactions)
    !end if

    !if(allocated(dissocReactions)) then
    !    do i=1,numDissocReaction
    !        deallocate(dissocReactions(i)%reactants)
    !        deallocate(dissocReactions(i)%products)
    !        deallocate(dissocReactions(i)%min)
    !        deallocate(dissocReactions(i)%max)
    !    end do
    !    deallocate(dissocReactions)
    !end if

    !if(allocated(sinkReactions)) then
    !    do i=1,numSinkReaction
    !        deallocate(sinkReactions(i)%reactants)
    !        deallocate(sinkReactions(i)%min)
    !        deallocate(sinkReactions(i)%max)
    !    end do
    !    deallocate(sinkReactions)
    !end if

    !if(allocated(diffReactions)) then
    !    do i=1,numDiffReaction
    !        deallocate(diffReactions(i)%reactants)
    !        deallocate(diffReactions(i)%products)
    !       deallocate(diffReactions(i)%min)
    !        deallocate(diffReactions(i)%max)
    !    end do
    !    deallocate(diffReactions)
    !end if

    if(allocated(clusterReactions)) then
        do i=1,numClusterReaction
            deallocate(clusterReactions(i)%reactants)
            deallocate(clusterReactions(i)%min)
            deallocate(clusterReactions(i)%max)
        end do
        deallocate(clusterReactions)
    end if

end subroutine deallocateInputs