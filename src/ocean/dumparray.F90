module m_dumparray
contains

subroutine dumparray(array1, prob_desc, arrayname)
    use, intrinsic :: iso_c_binding
    use m_par, only: rundir
    implicit none
    real, intent(in) :: array1(:)
    character(kind=c_char) :: prob_desc(*)
    character(len = *), intent(in) :: arrayname
    character(:), allocatable :: filename
    real, dimension(:), allocatable :: array2
    integer :: i, strlen
    logical :: file_exists

    strlen=0
    do
        if (prob_desc(strlen+1) == C_NULL_CHAR) exit
        strlen = strlen + 1
    end do

    allocate(character(len = strlen) :: filename)
    filename = rundir//transfer(prob_desc(1:strlen), filename(:strlen))//&
        "-"//arrayname//".dat"

    inquire(file=filename, exist=file_exists)
    open(37, file = filename)
    deallocate(filename)
    if(file_exists) then
        allocate(array2(size(array1)))
        read(37,*) array2

        if(any(array1.NE.array2)) then
            write(*,*) "Mismatch!"
            error stop
        endif
        deallocate(array2)
    else
        write(37,*) (array1(i), i=lbound(array1,1), ubound(array1,1))
    endif
    close(37)
end subroutine dumparray
end module m_dumparray
