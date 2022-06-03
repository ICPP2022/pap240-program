!***************************************************************************************************
!> Module update reactions
!***************************************************************************************************
module mod_update_reactions
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none

contains

!***************************************************************************************************
!> subroutine delete_reactions(def, def_prev, diff, cell, level): first delete all reactions associated def, then delete def
!***************************************************************************************************
subroutine delete_reactions(def, def_prev, cell, level)

    type(defect), pointer, intent(inout) :: def, def_prev
    integer, intent(in) :: cell, level

    type(defect), pointer :: def2
    type(reaction), pointer :: reac, reac_prev
    integer :: level2, i, j, count, dir, negDir

    !<First: clear def%reactionList
    nullify(reac)
    reac => def%reactionList
    do while(associated(reac))
        def%reactionList => reac%next
        totalRateVol(cell) = totalRateVol(cell) - reac%rate
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

    !<Second: clear associated 2nd reactions in other %reactionList
    !def is mobile, firstly, check the whole DRL and delete 2nd reactions associated with def,
    !secondly, update diff reaction in neighbors with the same defectTpye of def
    if(def%diff > 0d0) then
        !<delete associated 2nd reactions
        do level2=1, LEVELS
            nullify(def2)
            def2 => myDRL(level2,cell)%defect_all
            do while(associated(def2))
                nullify(reac)
                reac => def2%reactionList
                nullify(reac_prev)
                reacLoop: do while(associated(reac))
                    if(reac%numReactants == 2) then !find in 2nd reactions
                        count = 0
                        do j=1,SPECIES
                            if(reac%reactants(j,1) == def%defectType(j)) then
                                count = count + 1
                            end if
                        end do
                        if(count == SPECIES) then !find this clustering reaction
                            totalRateVol(cell) = totalRateVol(cell) - reac%rate
                            if(.NOT. associated(reac_prev)) then    !reac is the first node
                                def2%reactionList => reac%next
                            else
                                reac_prev%next => reac%next
                            end if
                            deallocate(reac%reactants)
                            if(allocated(reac%products)) then
                                deallocate(reac%products)
                            end if
                            deallocate(reac)
                            exit reacLoop
                        end if
                    end if
                    reac_prev => reac
                    reac => reac%next
                end do  reacLoop
                def2 => def2%next
            end do
        end do

        !<Third: update diff reactions in neighbors
        do dir=1, 6
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
                    count = 0
                    do j=1, SPECIES
                        if(def2%defectType(j) == def%defectType(j)) then
                            count = count +1
                        end if
                    end do
                    if(count == SPECIES) then
                        exit defLoop
                    end if
                    def2 => def2%next
                end do defLoop

                if(associated(def2)) then   !find it
                    call update_diff_reaction_cellTOcell(def2,myMesh(cell)%neighbor(dir), cell, negDir, &
                            myProc%taskid, myProc%taskid)
                end if
            end if
        end do
    end if

    !<Fourth: delete def from %defect_all
    if(.NOT. associated(def_prev)) then     !def is the first node
        myDRL(level,cell)%defect_all => def%next
    else
        def_prev%next => def%next
    end if
    deallocate(def%defectType)
    deallocate(def)

end subroutine

!*****************************************************************************************
!>Subroutine update_1st_Reaction: add or update dissociation and  sinkRemoval to the reactionList.
!Input: defectType, cell1
!*****************************************************************************************
subroutine update_1st_Reaction(def, cell)
    !<formal parameters
    type(defect), pointer, intent(inout) :: def
    integer, intent(in) :: cell
    !<local variables
    type(reaction), pointer :: reac, reac_prev, reacTemp
    integer :: i, j, count, count2, pointDefect(SPECIES)
    double precision :: sinkRate, dissRate

    !************************************
    !<Update sinkRemoval reaction
    !************************************
    count = 0
    do j=2, SPECIES
        if(def%defectType(j) == 0) then
            count = count + 1
        end if
    end do

    if((count==SPECIES-1) .AND. def%diff>0d0) then  !可动且是自缺陷
        !compute reaction rate
        sinkRate = compute_sink_rate(def%defectType, def%num, def%diff)

        !find sinkRemoval reaction
        nullify(reac)
        nullify(reac_prev)
        reac => def%reactionList
        findSink: do while(associated(reac))
            if(reac%numReactants==1 .AND. reac%numProducts==0) then
                exit findSink
            end if
            reac_prev => reac
            reac => reac%next
        end do findSink

        !update sinkRemoval reaction
        if(associated(reac) .AND. sinkRate<=0d0) then   !delete
            totalRateVol(cell) = totalRateVol(cell) - reac%rate
            if(associated(reac_prev)) then
                reac_prev%next => reac%next
            else
                def%reactionList => reac%next
            end if
            deallocate(reac)
        else if(associated(reac) .AND. sinkRate>0d0) then   !update
            totalRateVol(cell) = totalRateVol(cell) - reac%rate + sinkRate
            reac%rate = sinkRate
        else if(.NOT. associated(reac) .AND. sinkRate>0d0) then !add
            totalRateVol(cell) = totalRateVol(cell) + sinkRate
            nullify(reacTemp)
            allocate(reacTemp)
            reacTemp%numReactants = 1
            reacTemp%numProducts = 0
            reacTemp%rate = sinkRate
            nullify(reacTemp%next)
            if(associated(reac_prev)) then
                reac_prev%next => reacTemp
            else
                def%reactionList => reacTemp
            end if
        end if
    end if

    !************************************
    !<Update dissociation reaction
    !************************************
    !<compute dissociation rate
    count = 0
    do j=1, SPECIES
        count = count + iabs(def%defectType(j))
    end do
    if(count > 1) then !def is a cluster
        !compute reaction rate
        do j=1, SPECIES
            pointDefect = 0
            dissRate = 0d0
            if(def%defectType(j) /= 0) then !dissociation
                pointDefect(j) = def%defectType(j)/iabs(def%defectType(j))
                if(j == SPECIES) then
                    if(SIAPinToggle == 'yes') then  !allow m_SIAs + m_SIAs -> im_SIAs
                        pointDefect(1) = pointDefect(j)
                        pointDefect(2:SPECIES) = 0
                    end if
                end if
                dissRate = compute_diss_rate(def%defectType, def%num, pointDefect)

                !find diss reaction
                nullify(reac)
                nullify(reac_prev)
                reac => def%reactionList
                findDiss: do while(associated(reac))
                    if(reac%numReactants==1 .AND. reac%numProducts==2) then
                        count2 = 0
                        do i=1, SPECIES
                            if(reac%products(i,1) == pointDefect(i)) then
                                count2 = count2 + 1
                            end if
                        end do
                        if(count2 == SPECIES) then
                            exit findDiss
                        end if
                    end if
                    reac_prev => reac
                    reac => reac%next
                end do findDiss

                !update diss reaction
                if(associated(reac) .AND. dissRate<= 0d0) then !delete
                    totalRateVol(cell) = totalRateVol(cell) - reac%rate
                    if(associated(reac_prev)) then
                        reac_prev%next => reac%next
                    else
                        def%reactionList => reac%next
                    end if
                    deallocate(reac%products)
                    deallocate(reac)
                else if(associated(reac) .AND. dissRate>0d0) then    !update
                    totalRateVol(cell) = totalRateVol(cell) - reac%rate + dissRate
                    reac%rate = dissRate
                else if(.NOT. associated(reac) .AND. dissRate>0d0) then  !add
                    nullify(reacTemp)
                    allocate(reacTemp)
                    reacTemp%numReactants = 1
                    reacTemp%numProducts = 2
                    allocate(reacTemp%products(SPECIES,1))
                    reacTemp%products(1:SPECIES,1) = pointDefect(1:SPECIES)
                    totalRateVol(cell) = totalRateVol(cell) + dissRate
                    reacTemp%rate = dissRate
                    nullify(reacTemp%next)
                    if(.NOT. associated(reac_prev)) then
                        def%reactionList => reacTemp
                    else
                        reac_prev%next => reacTemp
                    end if
                end if
            end if

        end do
    end if

end subroutine

!*****************************************************************************************
!>Subroutine update_2nd_Reaction: add 2nd-reaction to the reactionList.
!Input: defectType, cell1
!*****************************************************************************************
subroutine update_2nd_Reaction(def1, def2, cell)

    !<formal parameters
    type(defect), pointer, intent(inout) :: def1, def2
    integer, intent(in) :: cell

    !<local variables
    type(reaction), pointer :: reac, reac_prev, reacTemp
    integer :: count, numProducts, j, countM
    integer :: productTem(SPECIES)
    integer, allocatable :: products(:,:)   !1 or 2
    logical :: tag
    double precision :: clusRate

    clusRate = 0d0

    productTem = def1%defectType + def2%defectType
    countM = 0
    do j=2, SPECIES-1
        if(productTem(j) /= 0) then
            countM = countM + 1
        end if
    end do

    if(countM/=0 .AND. iabs(productTem(1))<productTem(SPECIES)) then  !V_(X)_I or 0_(X)_I
        productTem(SPECIES) = productTem(1) + productTem(SPECIES)
        productTem(1) = 0
        if(productTem(SPECIES) <= max3D) then
            productTem(1) = productTem(SPECIES)
            productTem(SPECIES) = 0
        end if
    else if(countM/=0 .AND. iabs(productTem(1))>=productTem(SPECIES)) then  !I_(X)_0 or V_(X)_I or V_(X)_0
        productTem(1) = productTem(1) + productTem(SPECIES)
        productTem(SPECIES) = 0
    end if

    if(countM/=0 .AND. (productTem(1)>0 .OR. productTem(SPECIES)>0)) then !I_(X)_0 or 0_(X)_I
        numProducts = 2
        allocate(products(SPECIES,2))
        products = 0
        products(2:SPECIES-1,2) = productTem(2:SPECIES-1)   !0_(X)_0
        products(1:SPECIES,1) = productTem(1:SPECIES) - products(1:SPECIES,2)   !I_0_0 or 0_0_I
    else    !1 or 0 product
        if(SIAPinToggle=='yes' .AND. def1%defectType(1)>max3D .AND. def2%defectType(1)>max3D) then
            productTem(SPECIES) = def1%defectType(1) + def2%defectType(1)
            productTem(1:SPECIES-1) = 0
        else if(productTem(SPECIES) > 0) then   !0/I/V_0_I
            if(productTem(1) > 0) then
                productTem(SPECIES) = productTem(1) + productTem(SPECIES)
                productTem(1) = 0
            else if(productTem(1) < 0) then !V_0_I
                if(iabs(productTem(1)) >= productTem(SPECIES)) then
                    productTem(1) = productTem(1) + productTem(SPECIES)
                    productTem(SPECIES) = 0
                else
                    productTem(SPECIES) = productTem(1) + productTem(SPECIES)
                    productTem(1) = 0
                    if(productTem(SPECIES) <= max3D) then
                        productTem(1) = productTem(SPECIES)
                        productTem(SPECIES) = 0
                    end if
                end if
            end if
        end if
        count = 0
        do j=1, SPECIES
            if(productTem(j) == 0) then
                count = count + 1
            end if
        end do
        if(count == SPECIES) then
            numProducts = 0
        else
            numProducts = 1
            allocate(products(SPECIES,1))
            products(1:SPECIES,1) = productTem
        end if
    end if

    !tag = .true.
    !if(def1%defectType(1)>0 .AND. (def2%defectType(1)==0 .AND. def2%defectType(2)>0)) then   !def is SIAs, def2 is Cu
    !    tag = .false.
    !else if(def2%defectType(1)>0 .AND. (def1%defectType(1)==0 .AND. def1%defectType(2)>0)) then  !def2 is SIA, def is Cus
    !    tag = .false.
    !end if

    !if(tag .eqv. .true.) then
        !productTem = def1%defectType + def2%defectType
        !if(productTem(1)>0 .AND. productTem(2) >0) then !I_Cu cluster, separate
        !    numProducts = 2
        !    allocate(products(SPECIES,2))
        !    products(1:SPECIES,1) = productTem
        !    products(1:SPECIES,2) = productTem
        !    products(2,1) = 0     !<product1 SIA_0_0
        !    products(1,2) = 0     !<product2 0_S_0
        !else    !1 or 0 product
            !<check annihilation
        !    count = 0
        !    do j=1, SPECIES
        !        if(productTem(j) == 0) then
        !            count = count + 1
        !        end if
        !    end do
        !    if(count == SPECIES) then   !annihilation
        !        numProducts = 0
        !    else    !1 product
        !        numProducts = 1
        !        allocate(products(SPECIES,1))
        !        products(1:SPECIES,1) = productTem(1:SPECIES)
        !    end if
        !end if

        !clusRate = compute_clus_rate(def1%defectType,def1%num,def1%diff,def2%defectType,def2%num,def2%diff)
        clusRate = compute_clus_rate(def1, def2)

        if(def1%diff <= 0d0) then    !def is immobile
            nullify(reac)
            reac => def1%reactionList
            call find_2nd_reaction(reac,reac_prev,numProducts,def2%defectType)
            !<if the reaction already exists, but the reaction rate is 0, then delete this reaction.
            if(associated(reac) .AND. clusRate<=0d0) then   !find it, delete
                totalRateVol(cell) = totalRateVol(cell) - reac%rate
                if(associated(reac_prev)) then  !not the first 2nd reaction
                    reac_prev%next => reac%next
                else    !reac is the first node of reactinList
                    def1%reactionList => reac%next
                end if
                deallocate(reac%reactants)
                if(allocated(reac%products)) then
                    deallocate(reac%products)
                end if
                deallocate(reac)
                !<if the reaction already exists, and the reaction rate isn't 0, then update the reaction rate of this reaction.
            else if(associated(reac) .AND. clusRate>0d0) then   !find it, update
                totalRateVol(cell) = totalRateVol(cell) - reac%rate + clusRate
                reac%rate = clusRate
                !<if the reaction doesn't exist, and the reaction rate isn't 0, then add this reaction to the reaction list.
            else if(.NOT. associated(reac) .AND. clusRate>0d0) then !add reacTemp to the end of reactionList
                totalRateVol(cell) = totalRateVol(cell) + clusRate
                nullify(reacTemp)
                allocate(reacTemp)
                reacTemp%numReactants = 2
                reacTemp%numProducts = numProducts
                allocate(reacTemp%reactants(SPECIES,1)) !store reactant2
                reacTemp%reactants(1:SPECIES,1) = def2%defectType(1:SPECIES)
                if(numProducts > 0) then
                    allocate(reacTemp%products(SPECIES,numProducts))
                    reacTemp%products = products
                end if
                !reacTemp%diffDir = 0
                reacTemp%rate = clusRate
                nullify(reacTemp%next)
                if(associated(reac_prev)) then  !reac_prev指向reac的前一个2nd reaction
                    reac_prev%next => reacTemp
                else   !def1%reactionList => null
                    def1%reactionList => reacTemp
                end if
            end if
        else    !def1 is mobile defect
            nullify(reac)
            reac => def1%reactionList
            call find_2nd_reaction(reac,reac_prev,numProducts,def2%defectType)
            if(associated(reac)) then
                if(clusRate <= 0d0) then    !delete
                    totalRateVol(cell) = totalRateVol(cell) - reac%rate
                    if(associated(reac_prev)) then  !not the first 2nd reaction
                        reac_prev%next => reac%next
                    else    !reac is the first node of reactinList
                        def1%reactionList => reac%next
                    end if
                    deallocate(reac%reactants)
                    if(allocated(reac%products)) then
                        deallocate(reac%products)
                    end if
                    deallocate(reac)
                else    !update %rate
                    totalRateVol(cell) = totalRateVol(cell) - reac%rate + clusRate
                    reac%rate = clusRate
                end if
            else
                !find this type of reaction in def2%reactionList
                nullify(reac)
                reac => def2%reactionList
                call find_2nd_reaction(reac,reac_prev,numProducts,def1%defectType)

                if(associated(reac) .AND. clusRate<=0d0) then   !find it, delete
                    totalRateVol(cell) = totalRateVol(cell) - reac%rate
                    if(associated(reac_prev)) then  !not the first 2nd reaction
                        reac_prev%next => reac%next
                    else    !reac is the first node of reactinList
                        def2%reactionList => reac%next
                    end if
                    deallocate(reac%reactants)
                    if(allocated(reac%products)) then
                        deallocate(reac%products)
                    end if
                    deallocate(reac)
                else if(associated(reac) .AND. clusRate>0d0) then   !find it, update
                    totalRateVol(cell) = totalRateVol(cell) - reac%rate + clusRate
                    reac%rate = clusRate
                    !<if the reaction doesn't exist, and the reaction rate isn't 0, then add this reaction to the reaction list.
                else if(.NOT. associated(reac) .AND. clusRate>0d0) then !add reacTemp to the end of reactionList
                    totalRateVol(cell) = totalRateVol(cell) + clusRate
                    nullify(reacTemp)
                    allocate(reacTemp)
                    reacTemp%numReactants = 2
                    reacTemp%numProducts = numProducts
                    allocate(reacTemp%reactants(SPECIES,1)) !store reactant2
                    reacTemp%reactants(1:SPECIES,1) = def1%defectType(1:SPECIES)
                    if(numProducts > 0) then
                        allocate(reacTemp%products(SPECIES,numProducts))
                        reacTemp%products = products
                    end if
                    !reacTemp%diffDir = 0
                    reacTemp%rate = clusRate
                    nullify(reacTemp%next)
                    if(associated(reac_prev)) then  !reac_prev指向reac的前一个2nd reaction
                        reac_prev%next => reacTemp
                    else   !def1%reactionList => null
                        def2%reactionList => reacTemp
                    end if
                end if
            end if
        end if

        if(allocated(products)) then
            deallocate(products)
        end if
    !end if

end subroutine

    !*****************************************************************************************
    !>Subroutine update_2nd_Reaction_compare
    !Input: defectType, cell1
    !*****************************************************************************************
    subroutine update_2nd_reactionCompare(def1, def2, cell)
        !<formal parameters
        type(defect), pointer, intent(inout) :: def1, def2
        integer, intent(in) :: cell

        !<local variables
        type(defect), pointer :: defTemp
        type(reaction), pointer :: reac, reac_prev, reacTemp
        integer :: numSame1, numSame2, count, countM, numProducts
        integer :: i, j, productTem(SPECIES)
        integer, allocatable :: products(:,:)   !1 or 2
        double precision :: clusRate
        logical flag

        do i=1, numClusterReaction
            flag = .FALSE.
            !*******************************************************
            !def1%defectType = clusterReactions(i)%reactants(:,1)
            !def2%defectType = clusterReactions(i)%reactants(:,2)
            numSame1 = 0
            if(def1%defectType(1)==0 .AND. clusterReactions(i)%reactants(1,1)==0) then    !0_()_()
                numSame1 = numSame1 + 1
            else if(def1%defectType(1)>0 .AND. clusterReactions(i)%reactants(1,1)==1) then    !SIA_()_()
                if(def1%defectType(1)>=clusterReactions(i)%min(1) .AND. &
                        (def1%defectType(1)<=clusterReactions(i)%max(1) .OR. clusterReactions(i)%max(1)==-1)) then
                    numSame1 = numSame1 + 1
                end if
            else if(def1%defectType(1)<0 .AND. clusterReactions(i)%reactants(1,1)==-1) then   !V_()_()
                if(iabs(def1%defectType(1))>=clusterReactions(i)%min(1) .AND. &
                        (iabs(def1%defectType(1))<=clusterReactions(i)%max(1) .OR. clusterReactions(i)%max(1)==-1)) then
                    numSame1 = numSame1 + 1
                end if
            end if
            do j=2, SPECIES
                if(def1%defectType(j)==0 .AND. clusterReactions(i)%reactants(j,1)==0) then
                    numSame1 = numSame1 + 1
                else if(def1%defectType(j)/=0 .AND. clusterReactions(i)%reactants(j,1)==1) then
                    if(def1%defectType(j)>=clusterReactions(i)%min(j) .AND. &
                            (def1%defectType(j)<=clusterReactions(i)%max(j) .OR. clusterReactions(i)%max(j)==-1)) then
                        numSame1 = numSame1 + 1
                    end if
                end if
            end do

            numSame2 = 0
            if(def2%defectType(1)==0 .AND. clusterReactions(i)%reactants(1,2)==0) then    !0_()_()
                numSame2 = numSame2 + 1
            else if(def2%defectType(1)>0 .AND. clusterReactions(i)%reactants(1,2)==1) then    !SIA_()_()
                if(def2%defectType(1)>=clusterReactions(i)%min(1+SPECIES) .AND. &
                        (def2%defectType(1)<=clusterReactions(i)%max(1+SPECIES) .OR. &
                                clusterReactions(i)%max(1+SPECIES)==-1)) then
                    numSame2 = numSame2 + 1
                end if
            else if(def2%defectType(1)<0 .AND. clusterReactions(i)%reactants(1,2)==-1) then   !V_()_()
                if(iabs(def2%defectType(1))>=clusterReactions(i)%min(1+SPECIES) .AND. &
                        (iabs(def2%defectType(1))<=clusterReactions(i)%max(1+SPECIES) .OR. &
                                clusterReactions(i)%max(1+SPECIES)==-1)) then
                    numSame2 = numSame2 + 1
                end if
            end if
            do j=2, SPECIES
                if(def2%defectType(j)==0 .AND. clusterReactions(i)%reactants(j,2)==0) then
                    numSame2 = numSame2 + 1
                else if(def2%defectType(j)/=0 .AND. clusterReactions(i)%reactants(j,2)==1) then
                    if(def2%defectType(j)>=clusterReactions(i)%min(j+SPECIES) .AND. &
                            (def2%defectType(j)<=clusterReactions(i)%max(j+SPECIES) .OR. &
                                    clusterReactions(i)%max(j+SPECIES)==-1)) then
                        numSame2 = numSame2 + 1
                    end if
                end if
            end do

            if(numSame1==SPECIES .AND. numSame2==SPECIES) then
                flag = .TRUE.
            else
                !*******************************************************
                !def1%defectType = clusterReactions(i)%reactants(:,2)
                !def2%defectType = clusterReactions(i)%reactants(:,1)
                numSame1 = 0
                if(def2%defectType(1)==0 .AND. clusterReactions(i)%reactants(1,1)==0) then    !0_()_()
                    numSame1 = numSame1 + 1
                else if(def2%defectType(1)>0 .AND. clusterReactions(i)%reactants(1,1)==1) then    !SIA_()_()
                    if(def2%defectType(1)>=clusterReactions(i)%min(1) .AND. &
                            (def2%defectType(1)<=clusterReactions(i)%max(1) .OR. clusterReactions(i)%max(1)==-1)) then
                        numSame1 = numSame1 + 1
                    end if
                else if(def2%defectType(1)<0 .AND. clusterReactions(i)%reactants(1,1)==-1) then   !V_()_()
                    if(iabs(def2%defectType(1))>=clusterReactions(i)%min(1) .AND. &
                            (iabs(def2%defectType(1))<=clusterReactions(i)%max(1) .OR. clusterReactions(i)%max(1)==-1)) then
                        numSame1 = numSame1 + 1
                    end if
                end if
                do j=2, SPECIES
                    if(def2%defectType(j)==0 .AND. clusterReactions(i)%reactants(j,1)==0) then
                        numSame1 = numSame1 + 1
                    else if(def2%defectType(j)/=0 .AND. clusterReactions(i)%reactants(j,1)==1) then
                        if(def2%defectType(j)>=clusterReactions(i)%min(j) .AND. &
                                (def2%defectType(j)<=clusterReactions(i)%max(j) .OR. clusterReactions(i)%max(j)==-1)) then
                            numSame1 = numSame1 + 1
                        end if
                    end if
                end do

                numSame2 = 0
                if(def1%defectType(1)==0 .AND. clusterReactions(i)%reactants(1,2)==0) then    !0_()_()
                    numSame2 = numSame2 + 1
                else if(def1%defectType(1)>0 .AND. clusterReactions(i)%reactants(1,2)==1) then    !SIA_()_()
                    if(def1%defectType(1)>=clusterReactions(i)%min(1+SPECIES) .AND. &
                            (def1%defectType(1)<=clusterReactions(i)%max(1+SPECIES) .OR. &
                                    clusterReactions(i)%max(1+SPECIES)==-1)) then
                        numSame2 = numSame2 + 1
                    end if
                else if(def1%defectType(1)<0 .AND. clusterReactions(i)%reactants(1,2)==-1) then   !V_()_()
                    if(iabs(def1%defectType(1))>=clusterReactions(i)%min(1+SPECIES) .AND. &
                            (iabs(def1%defectType(1))<=clusterReactions(i)%max(1+SPECIES) .OR. &
                                    clusterReactions(i)%max(1+SPECIES)==-1)) then
                        numSame2 = numSame2 + 1
                    end if
                end if
                do j=2, SPECIES
                    if(def1%defectType(j)==0 .AND. clusterReactions(i)%reactants(j,2)==0) then
                        numSame2 = numSame2 + 1
                    else if(def1%defectType(j)/=0 .AND. clusterReactions(i)%reactants(j,2)==1) then
                        if(def1%defectType(j)>=clusterReactions(i)%min(j+SPECIES) .AND. &
                                (def1%defectType(j)<=clusterReactions(i)%max(j+SPECIES) .OR. &
                                        clusterReactions(i)%max(j+SPECIES)==-1)) then
                            numSame2 = numSame2 + 1
                        end if
                    end if
                end do

                if(numSame1==SPECIES .AND. numSame2==SPECIES) then
                    flag = .TRUE.
                end if
            end if

            if(flag .EQV. .TRUE.) then
                productTem = def1%defectType + def2%defectType
                countM = 0
                do j=2, SPECIES-1
                    if(productTem(j) /=0 ) then
                        countM = countM + 1
                    end if
                end do

                if(countM/=0 .AND. iabs(productTem(1))<productTem(SPECIES)) then  !V_(X)_I
                    productTem(SPECIES) = productTem(1) + productTem(SPECIES)
                    productTem(1) = 0
                    if(productTem(SPECIES) <= max3D) then
                        productTem(1) = productTem(SPECIES)
                        productTem(SPECIES) = 0
                    end if
                else if(countM/=0 .AND. iabs(productTem(1))>=productTem(SPECIES)) then
                    productTem(1) = productTem(1) + productTem(SPECIES)
                    productTem(SPECIES) = 0
                end if

                if(countM/=0 .AND. (productTem(1)>0 .OR. productTem(SPECIES)>0)) then
                    numProducts = 2
                    allocate(products(SPECIES,2))   !I_(X)_0 or 0_(x)_I
                    products = 0
                    products(2:SPECIES-1,2) = productTem(2:SPECIES-1)   !0_(X)_0
                    products(1:SPECIES,1) = productTem(1:SPECIES) - products(1:SPECIES,2)   !I_0_0 or 0_0_I
                else    !1 or 0 product
                    if(SIAPinToggle=='yes' .AND. (def1%defectType(1)>max3D .AND. def2%defectType(1)>max3D)) then
                        productTem(SPECIES) = def1%defectType(1) + def2%defectType(1)
                        productTem(1:SPECIES-1) = 0
                    else if(productTem(SPECIES) > 0) then
                        if(productTem(1) > 0) then
                            productTem(SPECIES) = productTem(1) + productTem(SPECIES)
                            productTem(1) = 0
                        else if(productTem(1) < 0) then   !V(x)
                            if(iabs(productTem(1)) >= productTem(SPECIES)) then
                                productTem(1) = productTem(1) + productTem(SPECIES)
                                productTem(SPECIES) = 0
                            else
                                productTem(SPECIES) = productTem(1) + productTem(SPECIES)
                                productTem(1) = 0
                                if(productTem(SPECIES) <= max3D) then
                                    productTem(1) = productTem(SPECIES)
                                    productTem(SPECIES) = 0
                                end if
                            end if
                        end if
                    end if
                    !<check annihilation
                    count = 0
                    do j=1, SPECIES
                        if(productTem(j) == 0) then
                            count = count + 1
                        end if
                    end do
                    if(count == SPECIES) then
                        numProducts = 0
                    else
                        numProducts = 1
                        allocate(products(SPECIES,1))
                        products(1:SPECIES,1) = productTem(1:SPECIES)
                    end if
                end if

                !productTem = def1%defectType + def2%defectType
                !<Vn_S_0 + SIAm_0_0 -> SIA(m-n)_S_0
                !if(productTem(1)>0 .AND. productTem(2) >0) then !2 products
                !    numProducts = 2
                !    allocate(products(SPECIES,numProducts))
                !    products(1:SPECIES,1) = productTem
                !    products(2,1) = 0     !<SIA_0_0
                !    products(1:SPECIES,2) = productTem
                !    products(1,2) = 0     !<0_S_0
                !else    !1 or 0 product
                    !<check annihilation
                !    count = 0
                !    do j=1, SPECIES
                !        if(productTem(j) == 0) then
                !            count = count + 1
                !        end if
                !    end do
                !    if(count == SPECIES) then
                !        numProducts = 0
                !    else
                !        numProducts = 1
                !        allocate(products(SPECIES,1))
                !        products(1:SPECIES,1) = productTem
                !    end if
                !end if

                !clusRate = rate_clustering(def1%defectType,def1%num,def1%diff,def2%defectType,def2%num,def2%diff,&
                !        cell,clusterReactions(i)%fType)

                clusRate = rate_clustering(def1, def2, clusterReactions(i)%fType)

                if(def1%diff <= 0d0) then    !def is immobile
                    nullify(reac)
                    reac => def1%reactionList
                    call find_2nd_reaction(reac,reac_prev,numProducts,def2%defectType)
                    !<if the reaction already exists, but the reaction rate is 0, then delete this reaction.
                    if(associated(reac) .AND. clusRate<=0d0) then   !find it, delete
                        totalRateVol(cell) = totalRateVol(cell) - reac%rate
                        if(associated(reac_prev)) then  !not the first 2nd reaction
                            reac_prev%next => reac%next
                        else    !reac is the first node of reactinList
                            def1%reactionList => reac%next
                        end if
                        deallocate(reac%reactants)
                        if(allocated(reac%products)) then
                            deallocate(reac%products)
                        end if
                        deallocate(reac)
                        !<if the reaction already exists, and the reaction rate isn't 0, then update the reaction rate of this reaction.
                    else if(associated(reac) .AND. clusRate>0d0) then   !find it, update
                        totalRateVol(cell) = totalRateVol(cell) - reac%rate + clusRate
                        reac%rate = clusRate
                        !<if the reaction doesn't exist, and the reaction rate isn't 0, then add this reaction to the reaction list.
                    else if(.NOT. associated(reac) .AND. clusRate>0d0) then !add reacTemp to the end of reactionList
                        totalRateVol(cell) = totalRateVol(cell) + clusRate
                        nullify(reacTemp)
                        allocate(reacTemp)
                        reacTemp%numReactants = 2
                        reacTemp%numProducts = numProducts
                        allocate(reacTemp%reactants(SPECIES,1)) !store reactant2
                        reacTemp%reactants(1:SPECIES,1) = def2%defectType(1:SPECIES)
                        if(numProducts > 0) then
                            allocate(reacTemp%products(SPECIES,numProducts))
                            reacTemp%products = products
                        end if
                        reacTemp%rate = clusRate
                        nullify(reacTemp%next)
                        if(associated(reac_prev)) then  !reac_prev指向reac的前一个2nd reaction
                            reac_prev%next => reacTemp
                        else   !def1%reactionList => null
                            def1%reactionList => reacTemp
                        end if
                    end if
                else    !def1 is mobile defect
                    nullify(reac)
                    reac => def1%reactionList
                    call find_2nd_reaction(reac,reac_prev,numProducts,def2%defectType)
                    if(associated(reac)) then
                        if(clusRate <= 0d0) then    !delete
                            totalRateVol(cell) = totalRateVol(cell) - reac%rate
                            if(associated(reac_prev)) then  !not the first 2nd reaction
                                reac_prev%next => reac%next
                            else    !reac is the first node of reactinList
                                def1%reactionList => reac%next
                            end if
                            deallocate(reac%reactants)
                            if(allocated(reac%products)) then
                                deallocate(reac%products)
                            end if
                            deallocate(reac)
                        else    !update %rate
                            totalRateVol(cell) = totalRateVol(cell) - reac%rate + clusRate
                            reac%rate = clusRate
                        end if
                    else
                        !find this type of reaction in def2%reactionList
                        nullify(reac)
                        reac => def2%reactionList
                        call find_2nd_reaction(reac,reac_prev,numProducts,def1%defectType)

                        if(associated(reac) .AND. clusRate<=0d0) then   !find it, delete
                            totalRateVol(cell) = totalRateVol(cell) - reac%rate
                            if(associated(reac_prev)) then  !not the first 2nd reaction
                                reac_prev%next => reac%next
                            else    !reac is the first node of reactinList
                                def2%reactionList => reac%next
                            end if
                            deallocate(reac%reactants)
                            if(allocated(reac%products)) then
                                deallocate(reac%products)
                            end if
                            deallocate(reac)
                        else if(associated(reac) .AND. clusRate>0d0) then   !find it, update
                            totalRateVol(cell) = totalRateVol(cell) - reac%rate + clusRate
                            reac%rate = clusRate
                            !<if the reaction doesn't exist, and the reaction rate isn't 0, then add this reaction to the reaction list.
                        else if(.NOT. associated(reac) .AND. clusRate>0d0) then !add reacTemp to the end of reactionList
                            totalRateVol(cell) = totalRateVol(cell) + clusRate
                            nullify(reacTemp)
                            allocate(reacTemp)
                            reacTemp%numReactants = 2
                            reacTemp%numProducts = numProducts
                            allocate(reacTemp%reactants(SPECIES,1)) !store reactant2
                            reacTemp%reactants(1:SPECIES,1) = def1%defectType(1:SPECIES)
                            if(numProducts > 0) then
                                allocate(reacTemp%products(SPECIES,numProducts))
                                reacTemp%products = products
                            end if
                            !reacTemp%diffDir = 0
                            reacTemp%rate = clusRate
                            nullify(reacTemp%next)
                            if(associated(reac_prev)) then  !reac_prev指向reac的前一个2nd reaction
                                reac_prev%next => reacTemp
                            else   !def1%reactionList => null
                                def2%reactionList => reacTemp
                            end if
                        end if
                    end if
                end if
                exit
            end if
        end do

        if(allocated(products)) then
            deallocate(products)
        end if

    end subroutine

    !*****************************************************************************************
    !>Subroutine update_diff_Reaction: add or update diff reaction (diffuse from this cell to six neighbors).
    !*****************************************************************************************
    subroutine update_diff_reaction_cellTOsix(def, cell)

        !<formal parameters
        type(defect), pointer, intent(inout) :: def
        integer, intent(in) :: cell
        !<local variables
        type(reaction), pointer :: reac, reac_prev, reacTemp
        integer :: dir
        double precision :: rate

        do dir=1, 6
            rate = compute_diff_rate(def%defectType, def%diff, def%num,cell,myMesh(cell)%neighbor(dir),dir,&
                    myProc%taskid,myMesh(cell)%neighborProc(dir))

            !<find diff reaction in def%reactionList
            nullify(reac)
            nullify(reac_prev)
            reac => def%reactionList
            findDiff: do while(associated(reac))
                if(reac%numReactants==1 .AND. reac%numProducts==1) then
                    if(reac%cell(1)==myMesh(cell)%neighbor(dir) .AND. reac%taskid(1)==myMesh(cell)%neighborProc(dir)) then
                        exit findDiff
                    end if
                end if
                reac_prev => reac
                reac => reac%next
            end do findDiff

            !存在一种情况，该cell中有该缺陷，但其数量比周围的都小，那么该缺陷不往周围扩散
            if(associated(reac) .AND. rate<=0d0) then    !有diff，但每个方向的扩散速率都为0，则删除
                totalRateVol(cell) = totalRateVol(cell) - reac%rate
                if(associated(reac_prev)) then  !not the first 2nd reaction
                    reac_prev%next => reac%next
                else    !reac is the first node of reactinList
                    def%reactionList => reac%next
                end if
                deallocate(reac%cell)
                deallocate(reac%taskid)
                deallocate(reac)
            else if(associated(reac) .AND. rate>0d0) then   !有diff，同时6个扩散速率不全为0， 则更新
                totalRateVol(cell) = totalRateVol(cell) - reac%rate + rate
                reac%rate = rate
            else if((.NOT. associated(reac)) .AND. rate>0d0) then   !没有diff，同时6个扩散速率不全为0，则添加
                nullify(reacTemp)
                allocate(reacTemp)
                reacTemp%numReactants = 1
                reacTemp%numProducts = 1
                totalRateVol(cell) = totalRateVol(cell) + rate
                allocate(reac%cell(1))
                reac%cell(1) = myMesh(cell)%neighbor(dir)
                allocate(reac%taskid(1))
                reac%taskid(1) = myMesh(cell)%neighborProc(dir)
                reacTemp%rate = rate
                nullify(reacTemp%next)
                if(.NOT. associated(reac_prev)) then
                    def%reactionList => reacTemp
                else
                    reac_prev%next => reacTemp
                end if
            end if
        end do

    end subroutine

    !*****************************************************************************************
    !>Subroutine update_diff_Reaction: add or update diff reaction (diffuse from this cell to six neighbors).
    !*****************************************************************************************
    subroutine update_diff_reaction_cellTOcell(def, cell1, cell2, dir, proc1, proc2)

        !<formal parameters
        type(defect), pointer, intent(inout) :: def
        integer, intent(in) :: cell1, cell2, dir, proc1, proc2
        !<local variables
        type(reaction), pointer :: reac, reac_prev, reacTemp
        integer :: count, i
        double precision :: diffRate

        !<compute rate (cell1 to cell2)
        diffRate = compute_diff_rate(def%defectType,def%diff,def%num,cell1,cell2,dir,proc1,proc2)

        !<find diff reaction
        nullify(reac)
        nullify(reac_prev)
        reac => def%reactionList
        findDiff: do while(associated(reac))
            if(reac%numReactants==1 .AND. reac%numProducts==1) then
                if(reac%cell(1)==cell2 .AND. reac%taskid(1)==proc2) then
                    exit findDiff
                end if
            end if
            reac_prev => reac
            reac => reac%next
        end do findDiff

        !<update diff reaction
        if(associated(reac) .AND. diffRate<=0d0) then   !find
            totalRateVol(cell1) = totalRateVol(cell1) - reac%rate
            if(associated(reac_prev)) then
                reac_prev%next => reac%next
            else
                def%reactionList => reac%next
            end if
            deallocate(reac%cell)
            deallocate(reac%taskid)
            deallocate(reac)
        else if(associated(reac) .AND. diffRate>0d0) then   !update
            totalRateVol(cell1) = totalRateVol(cell1) - reac%rate + diffRate
            reac%rate = diffRate
        else if(.NOT. associated(reac) .AND. diffRate>0d0) then  !add
            totalRateVol(cell1) = totalRateVol(cell1) + diffRate
            nullify(reacTemp)
            allocate(reacTemp)
            reacTemp%numReactants = 1
            reacTemp%numProducts = 1
            allocate(reacTemp%cell(1))
            reacTemp%cell(1) = cell2
            allocate(reacTemp%taskid(1))
            reacTemp%taskid(1) = proc2
            reacTemp%rate = diffRate
            nullify(reacTemp%next)
            if(.NOT. associated(reac_prev)) then
                def%reactionList => reacTemp
            else
                reac_prev%next => reacTemp
            end if
        end if

    end subroutine

    !*****************************************************************************************
    !>Subroutine find_2ndReaction_inDefCurrent: find 2st-reaction in def%reactionList
    !*****************************************************************************************
    subroutine find_2nd_reaction(reac,reac_prev,numProducts,defectType)

        !<formal parameters
        type(reaction), pointer, intent(inout) :: reac, reac_prev
        integer, intent(in) :: numProducts
        integer, intent(in), dimension(SPECIES) :: defectType
        !<local variables
        integer :: numSame, j

        nullify(reac_prev)
        outer: do while(associated(reac))
            if(reac%numReactants==2 .AND. reac%numProducts == numProducts) then
                numSame = 0
                do j=1, SPECIES !对比反应物2的类型
                    if(reac%reactants(j,1) == defectType(j)) then    !reac%reactants(j,2)数组越界
                        numSame = numSame + 1
                    end if
                end do
                if(numSame == SPECIES) then !find
                    exit outer
                end if
            end if
            reac_prev => reac
            reac => reac%next
        end do outer

    end subroutine

    !*****************************************************************************************
    !>Module: compute reaction rates, contains
    !       0-th reaction: Implantation
    !       1-st reaction: Dissociation, SinkRemoval
    !       2-nd reaction: Clustering
    !       Diffusion
    !*****************************************************************************************

    !*****************************************************************************************
    !>Function compute_0nd_rate(): calculates reaction rates of 0-order reactions.
    !*****************************************************************************************
    !double precision function compute_0nd_rate()

    !    double precision :: rate

    !    rate = 0d0
    !    if(irradiationType=='FrenkelPair') then
    !        if(implantType == 'uniform') then
    !            rate = (meshLength**(3d0))*dpaRate/atomVol
    !        else if(implantType == 'nonUniform') then
                !adding...
    !        else
    !            write(*,*) 'Error implant type not recognized'
    !        end if
    !    else if(irradiationType=='Cascade') then
    !        if(implantType == 'uniform') then
    !            if(implantScheme == 'MonteCarlo') then
    !                rate = (meshLength**(3d0))*dpaRate/(numDisplacedAtoms*atomVol)
    !            else if(implantScheme == 'explicit') then
    !                rate = 0d0
    !            end if
    !        else if(implantType == 'nonUniform') then
                !adding...
    !        else
    !            write(*,*) 'Error implant type not recognized'
    !        end if
    !    end if
    !    compute_0nd_rate = rate

    !end function

    !*****************************************************************************************
    !>Function compute_0nd_rate: calculates reaction rates of 0-order reactions.
    !*****************************************************************************************
    double precision function compute_0nd_rate(numReactants, numProducts, cell)

        integer, intent(in) :: numReactants, numProducts, cell
        double precision :: rate, zCoord, dpaRateLocal, HedpaRateLocal

        rate = 0d0
        if(numReactants==0 .AND. numProducts==2) then   !FrenkelPair
            if(implantType == 'uniform') then
                rate = (meshLength**(3d0))*dpaRate/atomVol
            else if(implantType == 'nonUniform') then
                zCoord = find_zCoord(cell)
                dpaRateLocal = find_dpaRate(zCoord)
                rate = (meshLength**(3d0))*dpaRateLocal/atomVol
            else
                write(*,*) 'Error implant type not recognized'
            end if
         else if(numReactants==-10 .AND. numProducts==0) then   !Cascade
            if(implantType == 'uniform') then
                if(implantScheme == 'MonteCarlo') then
                    rate = (meshLength**(3d0))*dpaRate/(numDisplacedAtoms*atomVol)
                else if(implantScheme == 'explicit') then
                    rate = 0d0
                end if
            else if(implantType == 'nonUniform') then
                if(implantScheme == 'MonteCarlo') then
                    zCoord = find_zCoord(cell)
                    dpaRateLocal = find_dpaRate(zCoord)
                    rate = (meshLength**(3d0))*dpaRateLocal/(numDisplacedAtoms*atomVol)
                else if(implantScheme == 'explicit') then
                    rate = 0d0
                end if
            else
                write(*,*) 'Error implant type not recognized'
            end if
         else if(numReactants==0 .AND. numProducts==1) then !He
            if(implantType == 'uniform') then
                rate = (meshLength**(3d0))*HeDPAratio*dpaRate/atomVol
            else if(implantType == 'nonUniform') then
                zCoord = find_zCoord(cell)
                HedpaRateLocal = find_HedpaRate(zCoord)
                rate = (meshLength**(3d0))*HedpaRateLocal/atomVol
            else
                write(*,*) 'Error implant type not recognized'
            end if
        end if
        compute_0nd_rate = rate

    end function


    !*****************************************************************************************
    !>Function compute_diss_rate: calculates reaction rates of dissociation reactions
    !Output: reaction rate of dissociation
    !*****************************************************************************************
    double precision function compute_diss_rate(defectType, num, pointDefect)

        integer, intent(in), dimension(SPECIES) :: defectType, pointDefect
        integer, intent(in) :: num

        integer :: size
        double precision :: diff, Eb, rate
        double precision, external :: compute_diffusivity

        rate = 0d0
        size = find_size(defectType)
        diff = compute_diffusivity(pointDefect)
        Eb = compute_binding(defectType, pointDefect)
        if(defectType(1)>max3D .OR. defectType(SPECIES)>0) then            !<SIA Loop
            rate = omega2D*dble(size)**(1d0/2d0)*diff*dexp(-Eb/(KB*temperature))*dble(num)
        else    !< 3D, or im_SIAs
            rate = omega*dble(size)**(1d0/3d0)*diff*dexp(-Eb/(KB*temperature))*dble(num)
        end if
        if(Eb == 0d0) then
            rate = 0d0
        end if
        !if(defectType(1)==-1 .AND. defectType(2)==1 .AND. Eb==0d0) then
        !    rate = 0d0
        !end if
        compute_diss_rate = rate
    end function

    !*****************************************************************************************
    !>Function compute_sink_rate: calculates reaction rates of sinkRemovel reactions
    !Return reaction rate of sinkRemoval
    !*****************************************************************************************
    double precision function compute_sink_rate(defectType, num, diff)

        integer, intent(in), dimension(SPECIES) :: defectType
        integer, intent(in) :: num
        double precision, intent(in) :: diff
        double precision :: rate

        rate = 0d0
        if(defectType(1) > 0) then            !<SIAs
            rate = sinks(1)*diff*dble(num)
        else if(defectType(1) < 0) then       !<V
            rate = sinks(2)*diff*dble(num)
        else
            rate = 0d0
        end if
        compute_sink_rate = rate
    end function

    !*****************************************************************************************
    !>Function compute_clus_rate: calculates reaction rates of 2-order reactions.
    !Output: reaction rate of clustering
    !*****************************************************************************************
    double precision function compute_clus_rate(def1, def2)

        type(defect), pointer, intent(in) :: def1, def2
        !integer, intent(in), dimension(SPECIES) :: defectType1, defectType2
        !integer, intent(in) :: num1, num2
        !double precision, intent(in) :: diff1, diff2

        integer :: count, j, defectType1(SPECIES), defectType2(SPECIES)
        double precision :: size1, size2, num1, num2, diff1, diff2, Ztemp, vol, rate

        rate = 0d0
        defectType1 = def1%defectType
        defectType2 = def2%defectType
        num1 = dble(def1%num)
        num2 = dble(def2%num)
        diff1 = def1%diff
        diff2 = def2%diff
        size1 = dble(find_size(defectType1))
        size2 = dble(find_size(defectType2))
        vol = meshLength**3d0
        !num11 = dble(num1)
        !num22 = dble(num2)

        if((defectType1(1)>0 .OR. defectType1(SPECIES)>0) .AND. (defectType2(1)>0 .OR. defectType2(SPECIES)>0)) then
            Ztemp = 1.2d0
        else
            Ztemp = 1.0d0
        end if

        count = 0
        do j=1, SPECIES
            if(defectType1(j) == defectType2(j)) then
                count = count + 1
            end if
        end do
        if(count == SPECIES) then   !Defects with the same type
            num2 = (num2-1.0)/2.0
            !num22 = (num22-1)/2.0
        end if

        rate = 0d0
        if((defectType1(1)>max3D .OR. defectType1(SPECIES)>0) .AND. &
                (defectType2(1)>max3D .OR. defectType2(SPECIES)>0)) then    !1D + 1D
            rate = (ZI*omegacircle1D*(size1**(1d0/2d0)+size2**(1d0/2d0)))**4d0*&
                    (diff1*num2+diff2*num1)*num1*num2*(atomVol/vol)**2d0
        else if((defectType1(1)>max3D .OR. defectType1(SPECIES)>0) .OR. &
                (defectType2(1)>max3D .OR. defectType2(SPECIES)>0)) then  !3D + 1D
            if(defectType1(1)>max3D) then   !defectType1 is 1D
                size1 = dble(find_size(defectType2))  !3D
                size2 = dble(find_size(defectType1))  !1D
                num1 = dble(def2%num)
                num2 = dble(def1%num)
                diff1 = def2%diff
                diff2 = def1%diff
            end if
            rate = (omega*size1**(1d0/3d0)+omega2D*size2**(1d0/2d0))*diff1*num1*num2*atomVol/vol+&
                    (Ztemp*(omegacircle1D*size2**(1d0/2d0)+omega1D*size1**(1d0/3d0)))**4d0*diff2*num2*&
                            num1**(2d0)*(atomVol/vol)**(2d0)
        else    !3D + 3D
            rate = Ztemp*omega*(size1**(1d0/3d0)+size2**(1d0/3d0))*(diff1+diff2)*num1*num2*atomVol/vol
        end if

        !if(defectType1(1)>max3D .AND. defectType2(1)>max3D) then        !1D + 1D
        !    rate = (ZI*omegacircle1D*(size1**(1d0/2d0)+size2**(1d0/2d0)))**4d0*&
        !            (diff1*num22+diff2*num11)*num11*num22*(atomVol/vol)**2d0
        !else if((defectType1(1)>0 .AND. defectType2(1)>max3D) .OR. (defectType1(1)>max3D .AND. defectType2(1)>0)) then  !3D(SIA) + 1D(SIA)
        !    if(defectType1(1)>max3D) then
        !        rate = (omega*size2**(1d0/3d0)+omega2D*size1**(1d0/2d0))*diff2*num22*num11*atomVol/vol+&
        !                (ZI*(omegacircle1D*size1**(1d0/2d0)+omega1D*size2**(1d0/3d0)))**4d0*diff1*num11*&
        !                        num22**(2d0)*(atomVol/vol)**(2d0)
        !    else if(defectType2(1)>max3D) then
        !        rate = (omega*size1**(1d0/3d0)+omega2D*size2**(1d0/2d0))*diff1*num11*num22*atomVol/vol+&
        !                (ZI*(omegacircle1D*size2**(1d0/2d0)+omega1D*size1**(1d0/3d0)))**4d0*diff2*num22*&
        !                        num11**(2d0)*(atomVol/vol)**(2d0)
        !    end if
        !else if((defectType1(1)<0 .AND. defectType2(1)>max3D) .OR. (defectType1(1)>max3D .AND. defectType2(1)<0)) then  !3D(V/VS) + 1D(SIA)
        !    if(defectType1(1)>max3D) then
        !        rate = (omega*size2**(1d0/3d0)+omega2D*size1**(1d0/2d0))*diff2*num22*num11*atomVol/vol+&
        !                (omegacircle1D*size1**(1d0/2d0)+omega1D*size2**(1d0/3d0))**4d0*Diff1*num11*&
        !                        num22**(2d0)*(atomVol/vol)**(2d0)
        !    else if(defectType2(1)>max3D) then
        !        rate = (omega*size1**(1d0/3d0)+omega2D*size2**(1d0/2d0))*diff1*num11*num22*atomVol/vol+&
        !                (omegacircle1D*size2**(1d0/2d0)+omega1D*size1**(1d0/3d0))**4d0*Diff2*num22*&
        !                        num1**(2d0)*(atomVol/vol)**(2d0)
        !    end if
        !else    !3D + 3D
        !    if((defectType1(1)>0 .AND. defectType1(1)<=max3D) .AND. (defectType2(1)>0 .AND. defectType2(1)<=max3D)) then
        !        rate = ZI*omega*(size1**(1d0/3d0)+size2**(1d0/3d0))*(diff1+diff2)*num11*num22*atomVol/vol
        !    else
        !        rate = omega*(size1**(1d0/3d0)+size2**(1d0/3d0))*(diff1+diff2)*num11*num22*atomVol/vol
        !    end if
        !end if

        compute_clus_rate = rate
    end function

    !*****************************************************************************************
    !>Function rate_clustering: calculates reaction rates of 2-order reactions.
    !Input: defectType1, defectType2, cell
    !Output: reaction rate of clustering
    !*****************************************************************************************
    double precision function rate_clustering(def1, def2, fType)

        type(defect), pointer, intent(in) :: def1, def2
        integer, intent(in) :: fType

        integer :: count, j, defectType1(SPECIES), defectType2(SPECIES)
        double precision :: num1, num2, diff1, diff2, size1, size2, vol, Ztemp, rate

        rate = 0d0
        defectType1 = def1%defectType
        defectType2 = def2%defectType
        num1 = dble(def1%num)
        num2 = dble(def2%num)
        diff1 = def1%diff
        diff2 = def2%diff
        size1 = dble(find_size(defectType1))
        size2 = dble(find_size(defectType2))
        vol = meshLength**3d0

        if((defectType1(1)>0 .OR. defectType1(SPECIES)>0) .AND. (defectType2(1)>0 .OR. defectType2(SPECIES)>0)) then
            Ztemp = 1.2d0
        else
            Ztemp = 1.0d0
        end if

        count = 0
        do j=1, SPECIES
            if(defectType1(j) == defectType2(j)) then
                count = count + 1
            end if
        end do
        if(count == SPECIES) then
            num2 = (num2-1.0)/2.0
        end if

        if(fType == 21) then        !<3D+3D
            rate = Ztemp*omega*(size1**(1d0/3d0)+size2**(1d0/3d0))*(diff1+diff2)*num1*num2*atomVol/vol
            !if(defectType1(1)>0  .AND. defectType2(1)>0 ) then  !<3D+3D: 3D(SIA) + 3D(SIA)
            !    rate = ZI*omega*(size1**(1d0/3d0)+size2**(1d0/3d0))*(diff1+diff2)*num1*num2*atomVol/vol
            !else
            !    rate = omega*(size1**(1d0/3d0)+size2**(1d0/3d0))*(diff1+diff2)*num1*num2*atomVol/vol
            !end if
        else if(fType == 22) then  !default: 1 is 3D, 2 is 1D
            if(defectType1(1)>max3D .OR. defectType1(SPECIES)>0) then   !defectType1 is 1D
                size1 = dble(find_size(defectType2))  !3D
                size2 = dble(find_size(defectType1))  !1D
                num1 = dble(def2%num)
                num2 = dble(def1%num)
                diff1 = def2%diff
                diff2 = def1%diff
            end if
            rate = (omega*size1**(1d0/3d0)+omega2D*size2**(1d0/2d0))*diff1*num1*num2*atomVol/vol+&
                    (Ztemp*(omegacircle1D*size2**(1d0/2d0)+omega1D*size1**(1d0/3d0)))**4d0*diff2*num2*num1**2d0&
                            *(atomVol/vol)**2d0
        !else if(fType == 22) then   !<3D(V/S) + 1D(SIA)
        !    if(defectType1(1) > max3D) then         !<defectType1 is 1D, defectType2 is 3D
        !        rate = (omega*size2**(1d0/3d0)+omega2D*size1**5d-1)*diff2*num22*num11*atomVol/vol+&
        !                (omegacircle1D*size1**5d-1+omega1D*size2**(1d0/3d0))**4d0*diff1*num11*num22**2d0&
        !                        *(atomVol/vol)**2d0
        !    else if(defectType2(1) > max3D) then    !<defectType1 is 3D, defectType2 is 1D
        !        rate = (omega*size1**(1d0/3d0)+omega2D*size2**5d-1)*diff1*num11*num22*atomVol/vol+&
        !                (omegacircle1D*size2**5d-1+omega1D*size1**(1d0/3d0))**4d0*diff2*num22*num11**2d0&
        !                        *(atomVol/vol)**2d0
        !    end if
        !else if(fType == 23) then    !<3D(SIA) + 1D(SIA)
        !    if(defectType1(1) > max3D) then   !<defectType1 is 1D, defectType2 is 3D
        !        rate = (omega*size2**(1d0/3d0)+omega2D*size1**5d-1)*diff2*num22*num11*atomVol/vol+&
        !                (ZI*(omegacircle1D*size1**5d-1+omega1D*size2**(1d0/3d0)))**4d0*diff1*num11*num22**2d0&
        !                        *(atomVol/vol)**2d0
        !    else if(defectType2(1) > max3D) then    !<defectType1 is 3D, defectType2 is 1D
        !        rate = (omega*size1**(1d0/3d0)+omega2D*size2**5d-1)*diff1*num11*num22*atomVol/vol+&
        !                (ZI*(omegacircle1D*size2**5d-1+omega1D*size1**(1d0/3d0)))**4d0*diff2*num22*num11**2d0&
        !                        *(atomVol/vol)**2d0
        !    end if
        else if(fType == 23) then    !<1D+1D
            rate = (ZI*omegacircle1D*(size1**(1d0/2d0)+size2**(1d0/2d0)))**4d0&
                    *(diff1*num2+diff2*num1)*num1*num2*(atomVol/vol)**2d0
        end if

        rate_clustering = rate
    end function

    !*****************************************************************************************
    !>Function compute_diff_rate: calculates reaction rates of diffusion.
    !Output: reaction rate of diffusion
    !*****************************************************************************************
    double precision function compute_diff_rate(defectType, diff, num,cell1,cell2,dir,proc1,proc2)

        integer, intent(in), dimension(SPECIES) :: defectType
        double precision, intent(in) :: diff
        integer, intent(in) ::  num, cell1, cell2, dir, proc1, proc2
        integer :: num2
        double precision :: rate, area, vol

        rate = 0d0
        area = meshLength**(2d0)
        vol = meshLength**(3d0)
        if(proc2 == -1) then   !<free surface
            num2 = 0
        else if(proc2 == proc1) then   !<in the same processor
            num2 = find_num_myDRL(defectType, cell2)
        else    !diffuse to neighbor processor
            num2 = find_num_myGhost(defectType, cell1, dir) !cell1 diffuse to cell2, should use cell1 to get the index
        end if
        rate = diff*area*(dble(num)/vol-dble(num2)/vol)/meshLength  !maybe rate<0d0
        if(rate < 0d0) then
            rate = 0d0
        end if
        compute_diff_rate = rate

    end function

    !*****************************************************************************************
    !>Function compute_binding: compute binding energy of defects
    !*****************************************************************************************
    double precision function compute_binding(defectType, product)

        integer, intent(in), dimension(SPECIES) :: defectType, product
        integer :: numSame, numSameProduct, i, j, k, size, numS, numV, numHe
        double precision :: Eb

        !<Initialize binding
        Eb = 0d0
        !<find defect in single defect list
        do i=1, numBindSingle
            numSame = 0
            numSameProduct = 0
            do j=1, SPECIES
                if(defectType(j) == bindSingle(i)%defectType(j)) then
                    numSame = numSame + 1
                end if
                if(product(j) == bindSingle(i)%product(j)) then
                    numSameProduct = numSameProduct + 1
                end if
            end do

            if(numSame == SPECIES .AND. numSameProduct == SPECIES) then
                Eb = bindSingle(i)%Eb
                exit
            end if
        end do

        !<find defect in function list
        k = 0
        if(i == numBindSingle+1) then
            do k=1, numBindFunc
                numSame = 0
                if(defectType(1)==0 .AND. bindFunc(k)%defectType(1)==0) then
                    numSame = numSame + 1
                else if(defectType(1)>0 .AND. bindFunc(k)%defectType(1)==1) then
                    if(defectType(1)>=bindFunc(k)%min(1) .AND. &
                            (defectType(1)<=bindFunc(k)%max(1) .OR. bindFunc(k)%max(1)==-1)) then
                        numSame = numSame + 1
                    end if
                else if(defectType(1)<0 .AND. bindFunc(k)%defectType(1)==-1) then
                    if(iabs(defectType(1))>=bindFunc(k)%min(1) .AND. &
                            (iabs(defectType(1))<=bindFunc(k)%max(1) .OR. bindFunc(k)%max(1)==-1) ) then
                        numSame = numSame + 1
                    end if
                end if
                do j=2,SPECIES
                    if(defectType(j)==0 .AND. bindFunc(k)%defectType(j)==0) then
                        numSame = numSame + 1
                    else if(defectType(j)/=0 .AND. bindFunc(k)%defectType(j)==1) then
                        if(defectType(j)>=bindFunc(k)%min(j) .AND. &
                                (defectType(j)<=bindFunc(k)%max(j) .OR. bindFunc(k)%max(j)==-1) ) then
                            numSame = numSame + 1
                        end if
                    end if
                end do

                numSameProduct = 0
                do j=1, SPECIES
                    if(product(j) == bindFunc(k)%product(j)) then
                        numSameProduct = numSameProduct + 1
                    end if
                end do

                if(numSame==SPECIES .AND. numSameProduct==SPECIES) then
                    if(bindFunc(k)%fType == 11) then           !<constant
                        Eb = bindFunc(k)%para(1)
                    else if(bindFunc(k)%fType == 12) then      !<Eb of V/SIA clusters (2 parameters)
                        size = iabs(defectType(1)) + sum(defectType(2:SPECIES))
                        Eb = bindFunc(k)%para(1)+(bindFunc(k)%para(2)-bindFunc(k)%para(1))*&
                                (dble(size)**(2d0/3d0)-dble(size-1)**(2d0/3d0))/(2d0**(2d0/3d0)-1d0)
                    else if(bindFunc(k)%fType == 13) then      !<Eb of Cun->Cu+Cu(n-1) (3 parameters)
                        numS = defectType(2)
                        Eb = bindFunc(k)%para(1)*KB-bindFunc(k)%para(2)*KB*temperature-(36d0*PI)**(1d0/3d0)*&
                                atomVol**(2d0/3d0)*bindFunc(k)%para(3)*(dble(numS)**(2d0/3d0)-dble(numS-1)**(2d0/3d0))
                    else if(bindFunc(k)%fType == 14) then      !<Eb of VmCun->Cu1+VmCu(n-1)) (3 parameters)
                        numS = defectType(2)
                        numV = iabs(defectType(1))
                        Eb = bindFunc(k)%para(1)+bindFunc(k)%para(2)*(dble(numS)**(0.85d0)-dble(numS+1)**(0.85d0))&
                                -bindFunc(k)%para(3)*(dble(numV)**(1d0/3d0)-dble(numV)**(2d0/3d0))
                    else if(bindFunc(k)%fType == 15) then      !<Eb of VmCun->V1+V(m-1)Cun (4 parameters)
                        numS = defectType(2)
                        numV = iabs(defectType(1))
                        Eb = bindFunc(k)%para(1)-bindFunc(k)%para(2)*(dble(numV)**(1d0/3d0)-dble(numV+1)**(1d0/3d0))&
                                +bindFunc(k)%para(3)*(dble(numV)**(2d0/3d0)-dble(numV+1)**(2d0/3d0))&
                                -bindFunc(k)%para(4)*dble(numS)*(dble(numV)**(1d0/3d0)-dble(numV+1)**(1d0/3d0)+&
                                        dble(numV)**(2d0/3d0)-dble(numV+1)**(2d0/3d0))
                    else if(bindFunc(k)%fType == 16) then   !<Eb of VnHem->He+VnHe(m-1)
                        numV = iabs(defectType(1))
                        numHe = defectType(2)
                        if(dble(numHe)/dble(numV) <= 0.5d0) then    !dose not dissociated
                            Eb = 0d0 !?????  !helium cannot dissociate from HeV clusters with mostly V
                        else
                            Eb = bindFunc(k)%para(1) - bindFunc(k)%para(2)*dlog(dble(numHe)/dble(numV))/dlog(10d0) - &
                                    bindFunc(k)%para(3)*(dlog(dble(numHe)/dble(numV))/dlog(10d0))**2
                        end if
                    else if(bindFunc(k)%fType == 17) then   !Eb of VnHem->V+V(n-1)Hem
                        numV = iabs(defectType(1))
                        numHe = defectType(2)
                        if(dble(numHe)/dble(numV) <= 0.5d0) then    !use vacancy cluster binding energy
                            Eb = bindFunc(k)%para(1)+(bindFunc(k)%para(2)-bindFunc(k)%para(1))*&
                                    (dble(numV)**(2d0/3d0)-dble(numV-1)**(2d0/3d0))/(2d0**(2d0/3d0)-1d0)
                        else
                            Eb = bindFunc(k)%para(3) + bindFunc(k)%para(4)*dlog(dble(numHe)/dble(numV))/dlog(10d0) + &
                                    bindFunc(k)%para(5)*(dlog(dble(numHe)/dble(numV))/dlog(10d0))**2
                        end if

                    end if
                    exit
                end if
            end do
        end if

        compute_binding = Eb
    end function compute_binding

    !*****************************************************************************************
    !>Function find_num_myDRL(defectType, cell): find number of defects
    !Output: number of defects with type of defectType
    !*****************************************************************************************
    integer function find_num_myDRL(defectType, cell)

        integer, intent(in), dimension(SPECIES) :: defectType
        integer, intent(in) :: cell
        type(defect), pointer :: def
        integer :: j, count, num, level
        integer, external :: find_level

        num = 0
        level = find_level(defectType)
        nullify(def)
        def => myDRL(level,cell)%defect_all
        do while(associated(def))
            count = 0
            do j=1, SPECIES
                if(defectType(j) == def%defectType(j)) then
                    count = count + 1
                end if
            end do
            if(count == SPECIES) then
                num = def%num
                exit
            end if
            def => def%next
        end do
        find_num_myDRL = num
    end function

    !*****************************************************************************************
    !>Function find_num_myGhost(defectType, cell): find number of defects
    !Output: number of defects with type of defectType
    !*****************************************************************************************
    integer function find_num_myGhost(defectType, cell, dir)

        integer, intent(in), dimension(SPECIES) :: defectType
        integer, intent(in) :: cell, dir

        type(defect), pointer :: def
        integer :: index, level
        integer :: j, count, num
        integer, external :: find_index, find_level

        num = 0
        index = find_index(cell, dir)
        level = find_level(defectType)
        nullify(def)
        def => myGhost(dir)%myDRL_ghost(level,index)%defect_mobile
        do while(associated(def))
            count = 0
            do j=1, SPECIES
                if(defectType(j) == def%defectType(j)) then
                    count = count + 1
                end if
            end do
            if(count == SPECIES) then
                num = def%num
                exit
            end if
            def => def%next
        end do
        find_num_myGhost = num
    end function

    !*****************************************************************************************
    !>Function find_size(defectType): find size of defectType
    !Output: size of defectType
    !*****************************************************************************************
    integer function find_size(defectType)

        integer, intent(in), dimension(SPECIES) :: defectType
        integer :: j, size, max

        max = maxval(defectType(2:SPECIES))
        if(iabs(defectType(1)) > max) then
            find_size = iabs(defectType(1))
        else
            find_size = max
        end if

    end function

    !*****************************************************************************************
    !>Function find_zCrood(cell): find the global zCoord of this cell
    !Output: global zCoord
    !*****************************************************************************************
    double precision function find_zCoord(cell)
        integer, intent(in) :: cell
        integer :: xGindex, yGindex, zGindex, temp

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

        find_zCoord = meshLength*dble(zGindex-1) + meshLength/2d0
    end function

    !*****************************************************************************************
    !>Function find_dpaRate(zCoord): find the local dpa rate
    !Output: global zCoord
    !*****************************************************************************************
    double precision function find_dpaRate(zCoord)
        double precision, intent(in) :: zCoord
        integer :: i
        double precision :: xi

        do i=1, numImpDatas
            if(impDPARates(1,i) == zCoord) then
                find_dpaRate = impDPARates(2,i)
                exit
            else if(zCoord<impDPARates(1,i) .AND. i/=1) then
                xi = (zCoord-impDPARates(1,i-1))/(impDPARates(1,i)-impDPARates(1,i-1))
                find_dpaRate = (impDPARates(2,i) - impDPARates(2,i-1)) * xi + impDPARates(2,i-1)
                exit
            else if(zCoord<impDPARates(1,i) .AND. i==1) then
                write(*,*) 'error DPA rate file starts after mesh in z-direction'
                exit
            end if
        end do

    end function

    !*****************************************************************************************
    !>Function find_HedpaRate(zCoord): find the local dpa rate
    !Output: global zCoord
    !*****************************************************************************************
    double precision function find_HedpaRate(zCoord)
        double precision, intent(in) :: zCoord
        integer :: i
        double precision :: xi

        do i=1, numImpDatas
            if(impDPARates(1,i) == zCoord) then
                find_HedpaRate = impDPARates(3,i)
                exit
            else if(zCoord<impDPARates(1,i) .AND. i/=1) then
                xi = (zCoord-impDPARates(1,i-1))/(impDPARates(1,i)-impDPARates(1,i-1))
                find_HedpaRate = (impDPARates(3,i) - impDPARates(3,i-1)) * xi + impDPARates(3,i-1)
                exit
            else if(zCoord<impDPARates(1,i) .AND. i==1) then
                write(*,*) 'error He implant rate rate file starts after mesh in z-direction'
                exit
            end if
        end do

    end function

end module



