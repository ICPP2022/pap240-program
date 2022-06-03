!****************************************************************************************
!>Subroutine: Choose Cascade: chooses one cascade randomly
!Inputs: CascadeList (global variable)
!Output: CascadeTemp (pointing at the cascade we want)
!*****************************************************************************************
subroutine chooseCascade(cascadeChoosed)
    use mod_structuretype
    use mod_globalvariables
    use mod_randdp
    implicit none

    type(cascade), pointer, intent(inout) :: cascadeChoosed
    double precision :: r, atemp
    double precision, external :: damageFun

    r=dprand()
    atemp=0d0
    cascadeChoosed=>cascadeList
    do while(associated(cascadeChoosed))
        atemp=atemp+1d0/dble(numCascades)
        if(r<=atemp) then
            numDamages = numDamages + damageFun(PKAenergy)
            exit
        else
            cascadeChoosed=>cascadeChoosed%next
        end if
    end do
end subroutine chooseCascade

!****************************************************************************************
!>Subroutine: chooseCascade_multiFiles(cascadeChoosed): chooses one cascade randomly from a number of cascade files
!Inputs: CascadeList (global variable)
!Output: CascadeTemp (pointing at the cascade we want)
!*****************************************************************************************
subroutine chooseCascade_multiFiles(cascadeChoosed)
    use mod_structuretype
    use mod_globalvariables
    use mod_randdp
    implicit none

    type(cascade), pointer, intent(inout) :: cascadeChoosed
    double precision :: energy, minEnergy
    double precision :: r, r1, atemp, atemp1, ksi1, cut_p
    integer :: i, min_i
    double precision, external :: sample_PKA_energy, damageFun

    if(PKAspectrum == 'yes') then
        energy = sample_PKA_energy()      !<energy > 0d0

        if(energy < cascadeLists(1)%PKAenergy) then !cascade file 按PKA能量大小排序
            min_i = 1
        else if(energy > cascadeLists(numCascadeFiles)%PKAenergy) then
            min_i = numCascadeFiles
        else
            do i=1, numCascadeFiles-1
                if(energy>=cascadeLists(i)%PKAenergy .AND. energy<cascadeLists(i+1)%PKAenergy) then
                    ksi1 = dprand()
                    cut_p = (energy-cascadeLists(i)%PKAenergy) / (cascadeLists(i+1)%PKAenergy-energy)
                    if(ksi1 < cut_p) then
                        min_i = i + 1
                    else
                        min_i = i
                    end if
                end if
            end do
        end if

        r = dprand()
        atemp = 0d0
        cascadeChoosed => cascadeLists(min_i)%listCascades
        do while(associated(cascadeChoosed))
            atemp = atemp + 1d0/dble(cascadeLists(min_i)%numCascades)
            if(r <= atemp) then
                numDamages = numDamages + damageFun(cascadeLists(min_i)%PKAenergy)
                exit
            else
                cascadeChoosed => cascadeChoosed%next
            end if
        end do

        !minEnergy = 1.0d8
        !outer1: do i=1, numCascadeFiles
        !    if(abs(energy-cascadeLists(i)%PKAenergy) < minEnergy) then    !<find the cascade file
        !        r = dprand()
        !        atemp = 0d0
        !        cascadeChoosed => cascadeLists(i)%listCascades
        !        do while(associated(cascadeChoosed))
        !            atemp = atemp + 1d0/dble(cascadeLists(i)%numCascades)
        !            if(r <= atemp) then
        !                numDamages = numDamages + damageFun(cascadeLists(i)%PKAenergy)
        !                exit outer1
        !            else
        !                cascadeChoosed => cascadeChoosed%next
        !            end if
        !        end do
        !    end if
        !end do outer1
    else    !not use PKA spectrum
        r1 = dprand()   !<used to choose a cascade file
        atemp1 = 0d0
        outer2: do i=1, numCascadeFiles
            atemp1 = atemp1 + 1d0/dble(numCascadeFiles)
            if(r1 <= atemp1) then   !<chosed a cascade file
                r = dprand()
                atemp = 0d0
                cascadeChoosed => cascadeLists(i)%listCascades
                do while(associated(cascadeChoosed))
                    atemp = atemp + 1d0/dble(cascadeLists(i)%numCascades)
                    if(r <= atemp) then
                        numDamages = numDamages + damageFun(cascadeLists(i)%PKAenergy)
                        exit outer2
                    else
                        cascadeChoosed => cascadeChoosed%next
                    end if
                end do
            end if
        end do outer2
    end if

end subroutine

!***************************************************************************************************
!> function sample_PKA_energy(): The PKA energy is sampled by linear interpolation
!***************************************************************************************************
double precision function sample_PKA_energy()
    use mod_structuretype
    use mod_globalvariables
    use mod_randdp
    implicit none

    double precision :: cpdf, r
    integer :: i, index
    double precision :: slope

    index=0
    r=dprand()
    do i=2, EPKAlist%size
        if(r > EPKAlist%cpdf(i-1) .AND. r < EPKAlist%cpdf(i)) then
            index=i
            exit
        end if
    end do

    slope = (EPKAlist%energy(index) - EPKAlist%energy(index-1))/(EPKAlist%cpdf(index) - EPKAlist%cpdf(index-1))
    sample_PKA_energy = slope*(r - EPKAlist%cpdf(index-1)) + EPKAlist%energy(index-1)

end function

!***************************************************************************************************
!> function damageFun(Epka): Calculate the number of displaced atoms according to the PKA energy
!***************************************************************************************************
double precision function damageFun(Epka)
    use mod_constants
    implicit none

    double precision, intent(in) :: Epka
    double precision :: vNRT

    vNRT = 0d0
    if(Epka <= ED) then
        vNRT = 0d0
    else if(Epka <= (2d0*ED/0.8d0)) then
        vNRT = 1d0
    else
        vNRT = 0.8d0*Epka/(2d0*ED)
    end if
    damageFun = vNRT

end function

!***************************************************************************************************
!> subroutine clearReactions(cell): delete all reactions in the cell
!***************************************************************************************************
subroutine clearReactions(cell)
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none

    integer, intent(in) :: cell
    type(defect), pointer :: def
    type(reaction), pointer :: reac
    integer :: level, j

    !*****************************************************
    !<Clear the reactions in myDRL.
    !*****************************************************
    do level=1, LEVELS
        nullify(def)
        def => myDRL(level,cell)%defect_all
        do while(associated(def))
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
            nullify(def%reactionList)
            def => def%next
        end do
    end do

end subroutine clearReactions

!***************************************************************************************************
!> subroutine reset_defectMobileList(cell): reset myDRL(:,cell)%defect_mobile
!***************************************************************************************************
subroutine reset_defectMobileList(cell)
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none

    integer, intent(in) :: cell
    type(defect), pointer :: def, def_m
    integer :: level

    !*****************************************************
    !<First: clear defect_mobile list in myDRL.
    !*****************************************************
    do level=1, LEVELS
        nullify(def_m)
        def_m => myDRL(level,cell)%defect_mobile
        do while(associated(def_m))
            myDRL(level,cell)%defect_mobile => def_m%next
            deallocate(def_m%defectType)
            deallocate(def_m)
            def_m => myDRL(level,cell)%defect_mobile
        end do
    end do

    !*****************************************************
    !<Second: create defect_mobile list in myDRL.
    !*****************************************************
    do level=1,LEVELS
        nullify(def)
        def => myDRL(level,cell)%defect_all
        nullify(def_m)
        do while(associated(def))
            if(def%diff > 0d0) then   !mobile defects
                if(.NOT. associated(myDRL(level,cell)%defect_mobile)) then    !add first node
                    allocate(def_m)
                    allocate(def_m%defectType(SPECIES))
                    def_m%defectType = def%defectType
                    def_m%num = def%num
                    def_m%diff = def%diff
                    nullify(def_m%next)
                    nullify(def_m%reactionList)
                    myDRL(level,cell)%defect_mobile => def_m
                else
                    allocate(def_m%next)
                    def_m => def_m%next
                    allocate(def_m%defectType(SPECIES))
                    def_m%defectType = def%defectType
                    def_m%num = def%num
                    def_m%diff = def%diff
                    nullify(def_m%next)
                    nullify(def_m%reactionList)
                end if
            end if
            def => def%next
        end do
    end do
end subroutine reset_defectMobileList

!***************************************************************************************************
!> subroutine resetReactionList: reset an entire reaction list in a single mesh
! Resets the reaction list for all reactions within a single cell (used when a cascade is created or
! deleted within that cell, all reaction rates change because volume changes) and diffReactions in neighbor cells
! Input: cell number
! Outputs: reaction rates for all reactions in a cell
!***************************************************************************************************
subroutine resetReactionList(cell)
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    use mod_update_reactions
    implicit none

    integer, intent(in) :: cell
    type(defect), pointer :: def, def2, def_m
    integer :: level, level2, level_m, dir, negDir, numSame, j

    !*****************************************************
    !<Add reactions for this cell
    !*****************************************************
    do level=1, LEVELS
        nullify(def)
        def => myDRL(level,cell)%defect_all  !pont to the first defect node
        do while(associated(def))
            !<Update dissociation reactions and sinkRemovel reactions
            call update_1st_Reaction(def, cell)
            !<Update diffusion reactions, contain this cell to neighbor cell and neighbor cell to this cell
            if(def%diff > 0d0) then
                !call update_diff_reaction_cellTOsix(def, cell)
                !update diff reaction in six neighbors
                do dir=1,6
                    call update_diff_reaction_cellTOcell(def,cell,myMesh(cell)%neighbor(dir),dir, &
                            myProc%taskid, myMesh(cell)%neighborProc(dir))

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
                                if(def2%defectType(j) == def%defectType(j)) then
                                    numSame = numSame +1
                                end if
                            end do
                            if(numSame == SPECIES) then
                                exit defLoop
                            end if
                            def2 => def2%next
                        end do defLoop

                        if(associated(def2)) then   !find it
                            call update_diff_reaction_cellTOcell(def2,myMesh(cell)%neighbor(dir),cell,negDir, &
                                    myProc%taskid, myProc%taskid)
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
                        !call update_2nd_Reaction(def2, def, cell)
                        call update_2nd_reactionCompare(def2, def, cell)
                        def2=>def2%next
                    end do
                end do
            end if

            def => def%next
        end do
    end do

end subroutine

!***************************************************************************************************
!> subroutine cascadeCombine_one: cascade defects are combined with the existing defects, and product2 is temporarily
!stored in extraList. Finally, the combined products and defects in extraList are all updated into che cellCH.
!***************************************************************************************************
subroutine cascadeCombine_one(cascadeChoosed, cellCH)
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    use mod_randdp
    implicit none

    type(cascade), pointer, intent(in) :: cascadeChoosed    !default: chosed a cascade
    integer, intent(in) :: cellCH

    type(cascadeDefect), pointer :: casDef
    type(defect), pointer :: def, def_prev, defTemp, extraList, extraDef, extraDef_prev, extraDefTemp
    integer, allocatable :: casDefArray(:,:)
    integer :: level, i, j, k, casDefType(SPECIES), extraProduct(SPECIES), count, countP, numSame
    double precision :: r1, atemp
    logical :: numTag, combineTag
    integer, external :: find_level
    double precision, external :: compute_diffusivity

    interface
        subroutine combineRules_one(defectType, casDefType, extraProduct, combineTag)
            use mod_constants
            use mod_globalvariables
            implicit none
            integer, intent(in), dimension(SPECIES) :: defectType
            integer, intent(inout) :: casDefType(SPECIES), extraProduct(SPECIES)
            logical, intent(inout) :: combineTag
        end subroutine

        subroutine locate_defect(defectType, def, def_prev, numSame)
            use mod_constants
            use mod_structuretype
            implicit none
            integer, intent(in), dimension(SPECIES) :: defectType
            type(defect), pointer, intent(inout) :: def, def_prev
            integer, intent(inout) :: numSame
        end subroutine
    end interface

    !***************************************************
    !<First: add cascade defects into casDefArray
    !***************************************************
    allocate(casDefArray(SPECIES,cascadeChoosed%numDefects))
    nullify(casDef)
    casDef => cascadeChoosed%listDefects
    i = 0
    do while(associated(casDef))
        i = i + 1
        casDefArray(1:SPECIES,i) = casDef%defectType(1:SPECIES)
        casDef => casDef%next
    end do

    !***************************************************
    !<Second: for each exited defect in cellCH, check whether it combines with defects in casDefArray
    !***************************************************
    nullify(extraList)
    do level=1, LEVELS
        nullify(def)
        def => myDRL(level, cellCH)%defect_all
        nullify(def_prev)
        do while(associated(def))
            numTag = .FALSE.

            combine: do k=1, def%num
                !<2.1: choose a casDefect
                r1 = dprand()
                atemp = 0d0
                do i=1, cascadeChoosed%numDefects
                    atemp = atemp + 1d0/dble(cascadeChoosed%numDefects)
                    if(r1 < atemp) then
                        exit
                    end if
                end do

                !2.2: combine
                if(i <= cascadeChoosed%numDefects) then
                    casDefType(1:SPECIES) = casDefArray(1:SPECIES,i)
                    call combineRules_one(def%defectType, casDefType, extraProduct, combineTag) !casDefType is changed
                    casDefArray(1:SPECIES,i) = casDefType(1:SPECIES)

                    if(combineTag .eqv. .true.) then    !combined
                        !<if extraProduct /= 0, add it into the temporary extraList
                        countP = 0
                        do j=1, SPECIES
                            if(extraProduct(j) == 0) then
                                countP = countP + 1
                            end if
                        end do
                        if(countP /= SPECIES) then !Insert extraProduct into extraList
                            nullify(extraDef)
                            extraDef => extraList
                            call locate_defect(extraProduct, extraDef, extraDef_prev, numSame)
                            if(associated(extraDef)) then
                                if(numSame == SPECIES) then !find it
                                    extraDef%num = extraDef%num + 1
                                else    !insert the extraProduct in front of extraDef
                                    nullify(extraDefTemp)
                                    allocate(extraDefTemp)
                                    allocate(extraDefTemp%defectType(SPECIES))
                                    extraDefTemp%defectType = extraProduct
                                    extraDefTemp%num = 1
                                    extraDefTemp%diff = 0d0
                                    nullify(extraDefTemp%reactionList)
                                    extraDefTemp%next => extraDef
                                    if(.NOT. associated(extraDef_prev)) then
                                        extraList => extraDefTemp
                                    else
                                        extraDef_prev%next => extraDefTemp
                                    end if
                                end if
                            else
                                nullify(extraDefTemp)
                                allocate(extraDefTemp)
                                allocate(extraDefTemp%defectType(SPECIES))
                                extraDefTemp%defectType = extraProduct
                                extraDefTemp%num = 1
                                extraDefTemp%diff = 0d0
                                nullify(extraDefTemp%reactionList)
                                nullify(extraDefTemp%next)
                                if(.NOT. associated(extraDef_prev)) then
                                    extraList => extraDefTemp
                                else
                                    extraDef_prev%next => extraDefTemp
                                end if
                            end if
                        end if

                        !<2.3: update def%num def可删除，因为前面已经将cellCH中的反应清除了
                        if(def%num <= 0) then
                            write(*,*) 'Error defect num zero combining with cascade defect'
                        else if(def%num == 1) then   !delete it
                            if(.NOT. associated(def_prev)) then !def is the first node in the list
                                myDRL(level, cellCH)%defect_all => def%next
                                deallocate(def%defectType)
                                deallocate(def)
                                def => myDRL(level, cellCH)%defect_all
                                numTag = .true.
                                exit combine
                            else
                                def_prev%next => def%next
                                deallocate(def%defectType)
                                deallocate(def)
                                def => def_prev%next
                                numTag = .true.
                                exit combine
                            end if
                        else    !def%num > 1
                            def%num = def%num - 1
                        end if
                    end if  !if(combineTag .eqv. .true.) then    !combined
                end if
            end do combine
            if(numTag .eqv. .false.) then
                def_prev => def
                def => def%next
            end if
        end do  !do while(associated(def))
    end do  !do level=1, LEVELS

    !***************************************************
    !<Third: add defects of casDefArray to the defectList
    !***************************************************
    do i=1, cascadeChoosed%numDefects
        !<check 0_0_0
        count = 0
        do j=1, SPECIES
            if(casDefArray(j,i) == 0) then
                count = count + 1
            end if
        end do
        if(count /= SPECIES) then   !add
            casDefType(1:SPECIES) = casDefArray(1:SPECIES,i)
            level = find_level(casDefType)
            nullify(def)
            def => myDRL(level, cellCH)%defect_all
            call locate_defect(casDefType, def, def_prev, numSame)
            if(associated(def)) then   !add
                if(numSame == SPECIES) then
                    def%num = def%num + 1
                else
                    nullify(defTemp)
                    allocate(defTemp)
                    allocate(defTemp%defectType(SPECIES))
                    defTemp%defectType = casDefType
                    defTemp%num = 1
                    defTemp%diff = compute_diffusivity(defTemp%defectType)
                    nullify(defTemp%reactionList)
                    defTemp%next => def
                    if(associated(def_prev)) then
                        def_prev%next => defTemp
                    else
                        myDRL(level, cellCH)%defect_all => defTemp
                    end if
                end if
            else
                nullify(defTemp)
                allocate(defTemp)
                allocate(defTemp%defectType(SPECIES))
                defTemp%defectType = casDefType
                defTemp%num = 1
                defTemp%diff = compute_diffusivity(defTemp%defectType)
                nullify(defTemp%reactionList)
                nullify(defTemp%next)
                if(associated(def_prev)) then   !def_prev is at the end of the list, def=>NULL
                    def_prev%next => defTemp
                else    !defCas will be the first node in the list
                    myDRL(level, cellCH)%defect_all => defTemp
                end if
            end if
        end if
    end do
    deallocate(casDefArray)

    !***************************************************
    !<Fourth: add defects in ectraList to the defectList
    !***************************************************
    nullify(extraDef)
    extraDef => extraList
    do while(associated(extraDef))
        level = find_level(extraDef%defectType)
        nullify(def)
        def => myDRL(level, cellCH)%defect_all
        call locate_defect(extraDef%defectType, def, def_prev, numSame)
        if(associated(def)) then   !add
            if(numSame == SPECIES) then
                def%num = def%num + extraDef%num
            else
                nullify(defTemp)
                allocate(defTemp)
                allocate(defTemp%defectType(SPECIES))
                defTemp%defectType = extraDef%defectType
                defTemp%num = extraDef%num
                defTemp%diff = compute_diffusivity(defTemp%defectType)
                nullify(defTemp%reactionList)
                defTemp%next => def
                if(associated(def_prev)) then
                    def_prev%next => defTemp
                else
                    myDRL(level, cellCH)%defect_all => defTemp
                end if
            end if
        else
            nullify(defTemp)
            allocate(defTemp)
            allocate(defTemp%defectType(SPECIES))
            defTemp%defectType = extraDef%defectType
            defTemp%num = extraDef%num
            defTemp%diff = compute_diffusivity(defTemp%defectType)
            nullify(defTemp%reactionList)
            nullify(defTemp%next)
            if(associated(def_prev)) then   !def_prev is at the end of the list, def=>NULL
                def_prev%next => defTemp
            else    !defCas will be the first node in the list
                myDRL(level, cellCH)%defect_all => defTemp
            end if
        end if
        extraDef => extraDef%next
    end do
    !<delete extraList
    extraDef => extraList
    do while(associated(extraDef))
        extraList => extraDef%next
        deallocate(extraDef%defectType)
        deallocate(extraDef)
        extraDef => extraList
    end do

end subroutine

!***************************************************************************************************
!> subroutine cascadeCombine_two: cascade defects are combined with the existing defects, and temporily allowing form
!SIACu clusters. Finally, SIA and Cu are split when the combined products are updated to cellCH.
!***************************************************************************************************
subroutine cascadeCombine_two(cascadeChoosed, cellCH)
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    use mod_randdp
    implicit none

    type(cascade), pointer, intent(in) :: cascadeChoosed    !default: chosed a cascade
    integer, intent(in) :: cellCH

    type(cascadeDefect), pointer :: casDef
    type(defect), pointer :: def, def_prev, defTemp
    integer, allocatable :: casDefArray(:,:)
    integer :: level, i, j, k, casDefType(SPECIES), extraProduct(SPECIES), count, countM, numSame
    double precision :: r1, atemp
    logical :: numTag
    integer, external :: find_level
    double precision, external :: compute_diffusivity

    interface
        subroutine combineRules_two(defectType, casDefType)
            use mod_constants
            use mod_globalvariables
            implicit none
            integer, intent(in), dimension(SPECIES) :: defectType
            integer, intent(inout) :: casDefType(SPECIES)
        end subroutine

        subroutine locate_defect(defectType, def, def_prev, numSame)
            use mod_constants
            use mod_structuretype
            implicit none
            integer, intent(in), dimension(SPECIES) :: defectType
            type(defect), pointer, intent(inout) :: def, def_prev
            integer, intent(inout) :: numSame
        end subroutine
    end interface

    !***************************************************
    !<First: add cascade defects into casDefArray
    !***************************************************
    allocate(casDefArray(SPECIES,cascadeChoosed%numDefects))
    nullify(casDef)
    casDef => cascadeChoosed%listDefects
    i = 0
    do while(associated(casDef))
        i = i + 1
        casDefArray(1:SPECIES,i) = casDef%defectType(1:SPECIES)
        casDef => casDef%next
    end do

    !***************************************************
    !<Second: for each exited defect in cellCH, check whether it combines with defects in casDefArray
    !***************************************************
    do level=1, LEVELS
        nullify(def)
        def => myDRL(level, cellCH)%defect_all
        nullify(def_prev)
        do while(associated(def))
            numTag = .FALSE.

            combine: do k=1, def%num
                !<2.1: choose a casDefect
                r1 = dprand()
                atemp = 0d0
                do i=1, cascadeChoosed%numDefects
                    atemp = atemp + 1d0/dble(cascadeChoosed%numDefects)
                    if(r1 < atemp) then
                        exit
                    end if
                end do

                !2.2: combine
                if(i <= cascadeChoosed%numDefects) then
                    casDefType(1:SPECIES) = casDefArray(1:SPECIES,i)
                    call combineRules_two(def%defectType, casDefType)   !casDefType is changed
                    casDefArray(1:SPECIES,i) = casDefType(1:SPECIES)

                    !<2.3: update def%num def可删除，因为前面已经将cellCH中的反应清除了
                    if(def%num <= 0) then
                        write(*,*) 'Error defect num zero combining with cascade defect'
                    else if(def%num == 1) then   !delete it
                        if(.NOT. associated(def_prev)) then !def is the first node in the list
                            myDRL(level, cellCH)%defect_all => def%next
                            deallocate(def%defectType)
                            deallocate(def)
                            def => myDRL(level, cellCH)%defect_all
                            numTag = .true.
                            exit combine
                        else
                            def_prev%next => def%next
                            deallocate(def%defectType)
                            deallocate(def)
                            def => def_prev%next
                            numTag = .true.
                            exit combine
                        end if
                    else    !def%num > 1
                        def%num = def%num - 1
                    end if
                end if
            end do combine
            if(numTag .eqv. .false.) then
                def_prev => def
                def => def%next
            end if
        end do  !do while(associated(def))
    end do  !do level=1, LEVELS

    !***************************************************
    !<Third: add defects of casDefArray to the defectList
    !***************************************************
    do i=1, cascadeChoosed%numDefects
        !<check 0_0_0
        count = 0
        do j=1, SPECIES
            if(casDefArray(j,i) == 0) then
                count = count + 1
            end if
        end do
        if(count /= SPECIES) then   !add
            !check casDefArray(1:SPECIES,i) is or not SIA(X) clusters
            countM = 0
            do j=2, SPECIES-1
                if(casDefArray(j,i) /= 0) then
                    countM = countM + 1
                end if
            end do
            if(countM /= 0) then    !casDefArray(1:SPECIES,i) is SIA(X) clusters
                extraProduct = 0
                extraProduct(2:SPECIES-1) = casDefArray(2:SPECIES-1,i)
                casDefArray(1:SPECIES,i) = casDefArray(1:SPECIES,i) - extraProduct(1:SPECIES)
            end if

            casDefType(1:SPECIES) = casDefArray(1:SPECIES,i)
            level = find_level(casDefType)
            nullify(def)
            def => myDRL(level, cellCH)%defect_all
            call locate_defect(casDefType, def, def_prev, numSame)
            if(associated(def)) then   !add
                if(numSame == SPECIES) then
                    def%num = def%num + 1
                else
                    nullify(defTemp)
                    allocate(defTemp)
                    allocate(defTemp%defectType(SPECIES))
                    defTemp%defectType = casDefType
                    defTemp%num = 1
                    defTemp%diff = compute_diffusivity(defTemp%defectType)
                    nullify(defTemp%reactionList)
                    defTemp%next => def
                    if(associated(def_prev)) then
                        def_prev%next => defTemp
                    else
                        myDRL(level, cellCH)%defect_all => defTemp
                    end if
                end if
            else
                nullify(defTemp)
                allocate(defTemp)
                allocate(defTemp%defectType(SPECIES))
                defTemp%defectType = casDefType
                defTemp%num = 1
                defTemp%diff = compute_diffusivity(defTemp%defectType)
                nullify(defTemp%reactionList)
                nullify(defTemp%next)
                if(associated(def_prev)) then   !def_prev is at the end of the list, def=>NULL
                    def_prev%next => defTemp
                else    !defCas will be the first node in the list
                    myDRL(level, cellCH)%defect_all => defTemp
                end if
            end if

            if(countM /= 0) then    !add extraProduct to cellCH
                level = find_level(extraProduct)
                nullify(def)
                def => myDRL(level, cellCH)%defect_all
                call locate_defect(extraProduct, def, def_prev, numSame)
                if(associated(def)) then   !add
                    if(numSame == SPECIES) then
                        def%num = def%num + 1
                    else
                        nullify(defTemp)
                        allocate(defTemp)
                        allocate(defTemp%defectType(SPECIES))
                        defTemp%defectType = extraProduct
                        defTemp%num = 1
                        defTemp%diff = compute_diffusivity(defTemp%defectType)
                        nullify(defTemp%reactionList)
                        defTemp%next => def
                        if(associated(def_prev)) then
                            def_prev%next => defTemp
                        else
                            myDRL(level, cellCH)%defect_all => defTemp
                        end if
                    end if
                else
                    nullify(defTemp)
                    allocate(defTemp)
                    allocate(defTemp%defectType(SPECIES))
                    defTemp%defectType = extraProduct
                    defTemp%num = 1
                    defTemp%diff = compute_diffusivity(defTemp%defectType)
                    nullify(defTemp%reactionList)
                    nullify(defTemp%next)
                    if(associated(def_prev)) then   !def_prev is at the end of the list, def=>NULL
                        def_prev%next => defTemp
                    else    !defCas will be the first node in the list
                        myDRL(level, cellCH)%defect_all => defTemp
                    end if
                end if

            end if
        end if
    end do
    deallocate(casDefArray)

end subroutine

!***************************************************************************************************
!>Subroutine combineRules_one: combine cascade defects with the exiting defects
!Inputs: defectType(SPECIES), casDefType(SPECIES)
!Outputs: casDefType(SPECIES), extraProduct(SPECIES)
!***************************************************************************************************
subroutine combineRules_one(defectType, casDefType, extraProduct, combineTag)
    use mod_constants
    use mod_globalvariables
    implicit none

    integer, intent(in), dimension(SPECIES) :: defectType
    integer, intent(inout) :: casDefType(SPECIES), extraProduct(SPECIES)
    logical, intent(inout) :: combineTag
    integer :: j, countM, product(SPECIES), count1, count2

    extraProduct=0
    combineTag = .true.
    if(defectType(1)==0 .AND. defectType(SPECIES)==0) then   !defectType is 0_(X)_0 and casDefType is I_0_0 or 0_0_I
        if(casDefType(1)>0 .OR. casDefType(SPECIES)>0) then
            combineTag = .false.
        end if
    else if(casDefType(1)==0 .AND. casDefType(SPECIES)==0) then  !casDefType is 0_(X)_0 and defectType is I_0_0 or 0_0_I
        if(defectType(1)>0 .OR. defectType(SPECIES)>0) then
            combineTag = .false.
        end if
    end if

    if(combineTag .eqv. .true.) then
        product = casDefType + defectType
        countM = 0
        do j=2, SPECIES-1
            if(product(j) /= 0) then
                countM = countM + 1
            end if
        end do

        if(countM/=0 .AND. product(1)<0 .AND. product(SPECIES)>0) then !product is V_(X)_0 or V_(X)_I or 0_(X)_0 or 0_(X)_I or I_(X)_0
            if(iabs(product(1))<product(SPECIES)) then
                product(SPECIES) = product(1) + product(SPECIES)
                product(1) = 0
                if(product(SPECIES) <= max3D) then
                    product(1) = product(SPECIES)
                    product(SPECIES) = 0
                end if
            else
                product(1) = product(1) + product(SPECIES)
                product(SPECIES) = 0
            end if
        end if
        if(countM/=0 .AND. (product(1)>0 .OR. product(SPECIES)>0)) then !product is V_(X)_0  or 0_(X)_0 or 0_(X)_I or I_(X)_0
            extraProduct(2:SPECIES-1) = product(2:SPECIES-1)    !extraProduct is 0_(X)_0
            casDefType(1:SPECIES) = product(1:SPECIES) - extraProduct(1:SPECIES) !defCas%defectType become I_0_0 or 0_0_I
        else    !1 product, 0_(X)_0 or V_(X)_0, V_0_0 or V_0_I or I_0_0, I_0_I, 0_0_I
            if(SIAPinToggle=='yes' .AND. casDefType(1)>max3D .AND. defectType(1)>max3D) then
                product(SPECIES) = product(1)
                product(1) = 0
            else if(product(SPECIES) > 0) then
                if(product(1) > 0) then !I_0_I
                    product(SPECIES) = product(1) + product(SPECIES)
                    product(1) = 0
                else if(product(1) <= 0) then
                    if(iabs(product(1)) < product(SPECIES)) then
                        product(SPECIES) = product(1) + product(SPECIES)
                        product(1) = 0
                        if(product(SPECIES) <= max3D) then
                            product(1) = product(SPECIES)
                            product(SPECIES) = 0
                        end if
                    else
                        product(1) = product(1) + product(SPECIES)
                        product(SPECIES) = 0
                    end if
                end if
            end if
            casDefType(1:SPECIES) = product(1:SPECIES)
        end if
    end if

end subroutine

!***************************************************************************************************
!>Subroutine combineRules_two: combine cascade defects with the exiting defects
!Inputs: defectType(SPECIES), casDefType(SPECIES)
!Outputs: casDefType(SPECIES)
!***************************************************************************************************
subroutine combineRules_two(defectType, casDefType)
    use mod_constants
    use mod_globalvariables
    implicit none

    integer, intent(in), dimension(SPECIES) :: defectType
    integer, intent(inout) :: casDefType(SPECIES)

    if(SIAPinToggle=='yes' .AND. casDefType(1)>max3D .AND. defectType(1)>max3D) then
        casDefType(SPECIES) = casDefType(1)
        casDefType(1) = 0
    end if
    casDefType = casDefType + defectType
    if(casDefType(SPECIES) > 0) then
        if(casDefType(1) > 0) then !I_(X))_I
            casDefType(SPECIES) = casDefType(1) + casDefType(SPECIES)
            casDefType(1) = 0
        else if(casDefType(1) <= 0) then
            if(iabs(casDefType(1)) < casDefType(SPECIES)) then
                casDefType(SPECIES) = casDefType(1) + casDefType(SPECIES)
                casDefType(1) = 0
                if(casDefType(SPECIES) <= max3D) then
                    casDefType(1) = casDefType(SPECIES)
                    casDefType(SPECIES) = 0
                end if
            else
                casDefType(1) = casDefType(1) + casDefType(SPECIES)
                casDefType(SPECIES) = 0
            end if
        end if
    end if

end subroutine


