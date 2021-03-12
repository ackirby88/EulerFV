module limiter_module
    implicit none
    contains

    subroutine limiter(gradL,gradR,gradLim)
        use my_kinddefs
        use inputs_module, only: limiterScheme
        use globals_module, only: numFields
        implicit none

        real(dp),intent(in)  :: gradL(numFields,2),gradR(numFields,2)
        real(dp),intent(out) :: gradLim

        real(dp)    :: leftMax(numFields),rightMax(numFields)
        real(dp)    :: beta,r
        integer(i4) :: i

        beta = 1.5

        do i = 1,numFields
            leftMax(i)  = max(abs(gradL(i,1)),abs(gradL(i,2)))
            rightMax(i) = max(abs(gradR(i,1)),abs(gradR(i,2)))
        end do

       !Calculate r
        r = 0.0
        do i = 1,numFields
            if (rightMax(i) /= 0.0) then
                r = max(r,abs(leftMax(i) / rightMax(i)))
            else
                r = max(r,abs((1.0 + leftMax(i)) / (1.0 + rightMax(i))))
            end if
        end do

        select case(limiterScheme)
            case('vanLeer')
                call vanLeer(r,gradLim)
            case('vanAlbada1')
                call vanAlbada1(r,gradLim)
            case('vanAlbada2')
                call vanAlbada2(r,gradLim)
            case('superbee')
                call superbee(r,gradLim)
            case('minmod')
                call minmod(r,gradLim)
            case('Koren')
                call Koren(r,gradLim)
            case('Osher')
                call Osher(r,gradLim,beta)    !pick a different beta
            case('UMIST')
                call UMIST(r,gradLim)
            case('ospre')
                call ospre(r,gradLim)
            case default
                print*,'No limiter was selected. Quitting.'
                stop
        end select
    end subroutine

    Subroutine vanLeer(r,lim)
        use my_kinddefs
        implicit none

        real(dp),intent(in) :: r
        real(dp),intent(out):: lim

        lim = (r + abs(r))/(1.0 + abs(r))
    end subroutine

    Subroutine vanAlbada1(r,lim)
        use my_kinddefs
        implicit none

        real(dp),intent(in) :: r
        real(dp),intent(out):: lim

        lim = (r*r + r)/(r*r + 1.0)
    end subroutine

    Subroutine vanAlbada2(r,lim)
        use my_kinddefs
        implicit none

        real(dp),intent(in) :: r
        real(dp),intent(out):: lim

        lim = (2.0*r)/(r*r + 1.0)
    end subroutine

    Subroutine superbee(r,lim)
        use my_kinddefs
        implicit none

        real(dp),intent(in) :: r
        real(dp),intent(out):: lim

        real(dp) :: val1,val2

        val1 = min(2.0*r,1.0)
        val2 = min(r,2.0)

        lim = max(0.,max(val1,val2))
    end subroutine

    Subroutine minmod(r,lim)
        use my_kinddefs
        implicit none

        real(dp),intent(in) :: r
        real(dp),intent(out):: lim

        lim = max(0.0,min(1.0,r))
    end subroutine

    Subroutine Koren(r,lim)
        use my_kinddefs
        implicit none

        real(dp),intent(in) :: r
        real(dp),intent(out):: lim

        lim = max(0.0,min(2.0*r,min(onethird*(1.0 + 2.0*r),2.0)))
    end subroutine

    Subroutine Osher(r,lim,beta)
        use my_kinddefs
        implicit none

        real(dp),intent(in) :: r
        real(dp),intent(out):: lim
        real(dp),intent(in) :: beta

        lim = max(0.0,min(r,beta))
    end subroutine

    Subroutine UMIST(r,lim)
        use my_kinddefs
        implicit none

        real(dp),intent(in) :: r
        real(dp),intent(out):: lim

        lim = max(0.0,min(min(2.0*r,0.25 + 0.75*r),min(0.75 + 0.25*r,2.0)))
    end subroutine

    Subroutine ospre(r,lim)
        use my_kinddefs
        implicit none

        real(dp),intent(in) :: r
        real(dp),intent(out):: lim

        lim = (1.5 * (r*r + r)) / (r*r + r + 1.0)
    end subroutine
end module limiter_module