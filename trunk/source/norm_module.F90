module norm_module
    implicit none
    contains

    subroutine L2_norm(nf,ne,A,norm)
        use my_kinddefs
        implicit none

        integer(i4),intent(in) :: nf,ne
        real(dp),   intent(in) :: A(nf*ne)
        real(dp),   intent(out):: norm

        integer(i4) :: i

        norm = 0.0

        do i = 1,nf*ne
            norm = norm + A(i)*A(i)
        end do
        norm = sqrt(norm/ne)

        if (isnan(norm)) then
            print*, 'There are NANs! Stopping.'
            stop
        end if
    end subroutine
end module