!*****************************************************************************************
!>Subroutine compute_diffusivity: compute diffusivities of defects
!This function has several functional forms for diffusivity, including immobile defects, constant functions, and mobile SIA loops.
!Additional functional forms can be added as needed.
!*****************************************************************************************
double precision function compute_diffusivity(defectType)
    use mod_constants
    use mod_globalvariables
    implicit none

    integer, intent(in), dimension(SPECIES) :: defectType
    integer :: numSame, i, j, k
    double precision :: D0, Em, diff

    !<initialize diff
    diff = 0d0
    !<find defect in single defect list
    do i=1, numDiffSingle
        numSame = 0
        do j=1, SPECIES
            if(defectType(j) == diffSingle(i)%defectType(j)) then
                numSame = numSame + 1
            end if
        end do

        if(numSame == SPECIES) then
            diff = diffSingle(i)%D0*dexp(-diffSingle(i)%Em/(KB*temperature))
            exit
        end if
    end do

    !<find defect in function list
    k=0
    if(i==numDiffSingle+1) then
        do k=1, numDiffFunc
            numSame = 0
            if(defectType(1)==0 .AND. diffFunc(k)%defectType(1)==0) then
                numSame = numSame + 1
            else if(defectType(1)>0 .AND. diffFunc(k)%defectType(1)==1) then
                if(defectType(1)>=diffFunc(k)%min(1) .AND. &
                        (defectType(1)<=diffFunc(k)%max(1) .OR. diffFunc(k)%max(1)==-1)) then
                    numSame = numSame + 1
                end if
            else if(defectType(1)<0 .AND. diffFunc(k)%defectType(1)==-1) then
                if(iabs(defectType(1))>=diffFunc(k)%min(1) .AND. &
                        (iabs(defectType(1))<=diffFunc(k)%max(1) .OR. diffFunc(k)%max(1)==-1)) then
                    numSame = numSame + 1
                end if
            end if
            do j=2, SPECIES
                if(defectType(j)==0 .AND. diffFunc(k)%defectType(j)==0) then
                    numSame = numSame + 1
                else if(defectType(j)/=0 .AND. diffFunc(k)%defectType(j)==1) then
                    if(defectType(j)>=diffFunc(k)%min(j) .AND. &
                            (defectType(j)<=diffFunc(k)%max(j) .OR. diffFunc(k)%max(j)==-1)) then
                        numSame = numSame + 1
                    end if
                end if
            end do

            if(numSame == SPECIES) then     !find it
                if(diffFunc(k)%fType==1) then    !<immobile defects
                    diff = 0d0
                else if(diffFunc(k)%fType==2) then   !<constant
                    diff = diffFunc(k)%para(1)
                else if(diffFunc(k)%fType==3) then   !<mobile SIA loops
                    D0 = diffFunc(k)%para(1)+diffFunc(k)%para(2)/dble(iabs(defectType(1)))**(diffFunc(k)%para(3))
                    Em = diffFunc(k)%para(4)+diffFunc(k)%para(5)/dble(iabs(defectType(1)))**(diffFunc(k)%para(6))
                    diff = D0*dexp(-Em/(KB*temperature))
                else if(diffFunc(k)%fType==4) then   !<mobile Cu clusters (Cu). D=D(1)/n
                    diff = diffFunc(k)%para(1)*dexp(-diffFunc(k)%para(2)/(KB*temperature))/dble(defectType(2))
                end if
                exit
            end if
        end do
    end if
    compute_diffusivity = diff

end function compute_diffusivity

!***************************************************************************************************
!>Subroutine locate_defect: points def at the appropriate defect in a linked list
!This subroutine places def on the defect being added to the system if it is already in the system
!If not, def points to the place after the inserted defect and def_prev points to the place before the inserted defect
!Ordering: first by SIA/V content, then by Scontent, then by He content
!Inputs: defectType
!Outputs: def and def_prev pointed at location in list, numSame
!***************************************************************************************************
subroutine locate_defect(defectType, def, def_prev, numSame)
    use mod_constants
    use mod_structuretype
    implicit none

    integer, intent(in), dimension(SPECIES) :: defectType
    type(defect), pointer, intent(inout) :: def, def_prev
    integer, intent(inout) :: numSame
    integer :: j

    nullify(def_prev)
    numSame = 0
    outer: do while(associated(def))
        numSame=0
        inter2: do j=1, SPECIES
            if(def%defectType(j) == defectType(j)) then
                numSame = numSame + 1
            else if(iabs(def%defectType(j)) > iabs(defectType(j))) then
                exit inter2
            end if
        end do inter2

        if(numSame == SPECIES) then
            exit outer
        else if(numSame==j-1 .AND. iabs(def%defectType(j))>iabs(defectType(j))) then
            exit outer
        else
            def_prev => def
            def => def%next
        end if
    end do outer

end subroutine locate_defect

!***************************************************************************************************
!>Function find_level(defectType)
!***************************************************************************************************
integer function find_level(defectType)
    use mod_constants
    implicit none

    integer, intent(in), dimension(SPECIES) :: defectType
    integer, dimension(SPECIES) :: key
    integer :: hash_value, i, count

    key = 0
    hash_value = 0
    do i=1, SPECIES
        if(defectType(i) /= 0) then
            key(i) = 1
        else
            key(i) = 0
        end if
    end do

    do i=1, SPECIES
        hash_value = hash_value+key(i)*2**(SPECIES-i)
    end do
    hash_value = hash_value + 1

    !conflict: (V_0_0) and (SIA_0_0)
    count = 0
    if(defectType(1) > 0) then    !defectType is SIA_0_0
        do i=2, SPECIES
            if(defectType(i) == 0) then
                count = count + 1
            end if
        end do
        if(count == (SPECIES-1)) then
            hash_value = 1
        end if
    end if

    find_level = hash_value

end function find_level

!***************************************************************************************************
!>Function find_index(cell, dir)    !cell can be locaCellID or ghostCellID
!***************************************************************************************************
integer function find_index(cell, dir)  !下述cell的坐标有问题
    use mod_constants
    use mod_structuretype
    use mod_globalvariables
    implicit none

    integer, intent(in) :: cell, dir
    integer :: x, y, z, temp, index

    index = 0
    if(mod(cell,numyLocal*numzLocal) /= 0) then
        x = cell/(numyLocal*numzLocal) + 1
    else
        x = cell/(numyLocal*numzLocal)
    end if
    temp = cell - numyLocal*numzLocal*(x-1)
    if(mod(temp, numzLocal) /= 0) then
        y = temp/numzLocal + 1
    else
        y = temp/numzLocal
    end if
    z = temp - numzLocal*(y-1)
    if(dir==1 .OR. dir==2) then         !<x-direction
        index = (y-1)*numzLocal+z
    else if(dir==3 .OR. dir==4) then    !<y-direction
        index = (x-1)*numzLocal+z
    else if(dir==5 .OR. dir==6) then    !<z-direction
        index = (x-1)*numyLocal+y
    end if
    find_index = index

end function find_index



